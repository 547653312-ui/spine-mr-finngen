"""
reverse_mr.py — 反向 MR：脊柱感染 / 骨髓炎  ->  免疫性状

难点与解法
----------
反向 MR 要把 FinnGen 表型当作**暴露**，因此需要它的全基因组显著位点。
FinnGen 元数据显示这些表型的显著位点极少（R11）：
    M13_OSTEOMYELVERTEB   0
    M13_SPONDYLOINFECTION 0
    M13_DISCITIS          1
    M13_DISCINFECTION     1
    M13_OSTEOMYELITIS     2
且公开接口不提供 top-hits 端点（/api/top_hits* 均 404），只能扫描汇总统计文件本身。

本脚本以**流式 BGZF 解压**扫描远端 765–808MB 文件（不落盘、不全量缓存），
按阈值抽取候选 IV；因 5e-8 下 IV 不足 3 个，主用 5e-6 次级阈值（并如实标注）。
结局端（免疫性状）通过 EBI 汇总统计 API 按 rsID 反查。

用法：
    python reverse_mr.py                 # 扫描 OM + DISCITIS
    python reverse_mr.py --pheno M13_OSTEOMYELITIS --pthresh 5e-6
"""

from __future__ import annotations

import argparse
import json
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Dict, List

import numpy as np
import pandas as pd

import gwas_io as gio
import mr_methods as mrm
import mr_pipeline as mp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS = os.path.join(ROOT, "results")
DATA = os.path.join(ROOT, "data")


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def scan_significant(phenocode: str, pthresh: float,
                     release: str = "R11", chunk: int = 8 << 20) -> pd.DataFrame:
    """
    流式扫描 FinnGen 远端 bgzip 文件，抽取 P < pthresh 的变异。
    逐块 HTTP Range 读取 -> BGZF 解块 -> 行解析，内存占用恒定。
    """
    url = f"{gio.FINNGEN_BUCKET}/finngen_{release}_{phenocode}.gz"
    cache = os.path.join(DATA, "cache", f"sig_{release}_{phenocode}_{pthresh:g}.csv")
    if os.path.exists(cache):
        log(f"  使用缓存 {os.path.basename(cache)}")
        return pd.read_csv(cache)

    bg = gio.BGZFRemote(url)
    total = bg.size
    log(f"  扫描 {phenocode} ({total/1e6:.0f} MB, P<{pthresh:g}) ...")
    rows: List[dict] = []
    pos = 0
    carry = b""          # 未消费完的压缩字节（不完整 BGZF 块）
    tail = ""            # 未消费完的文本行
    t0 = time.time()
    nline = 0
    while pos < total:
        raw = bg.read_range(pos, chunk)
        if not raw:
            break
        buf = carry + raw
        text, used = _inflate_prefix(buf)
        carry = buf[used:]
        pos += len(raw)
        if not text:
            if len(carry) > 4 << 20:      # 防御：异常时避免无限增长
                carry = b""
            continue
        text = tail + text
        lines = text.split("\n")
        tail = lines.pop()
        for ln in lines:
            if not ln or ln[0] == "#":
                continue
            f = ln.split("\t")
            if len(f) < 10:
                continue
            try:
                if float(f[6]) >= pthresh:
                    continue
            except ValueError:
                continue
            nline += 1
            rows.append({"chr": f[0], "pos": int(f[1]), "ref": f[2], "alt": f[3],
                         "rsids": f[4], "gene": f[5], "pval": float(f[6]),
                         "beta": float(f[8]), "sebeta": float(f[9]),
                         "af_alt": float(f[10]) if len(f) > 10 else np.nan})
        if pos % (80 << 20) < chunk:
            log(f"    {pos/total*100:5.1f}%  命中 {nline}  "
                f"({pos/1e6/max(time.time()-t0,1e-9):.1f} MB/s)")
    df = pd.DataFrame(rows)
    log(f"  完成：{len(df)} 个 P<{pthresh:g} 的变异，用时 {time.time()-t0:.0f}s")
    os.makedirs(os.path.dirname(cache), exist_ok=True)
    df.to_csv(cache, index=False)
    return df


def _inflate_prefix(buf: bytes) -> tuple[str, int]:
    """解压 buf 中所有完整的 BGZF 块，返回 (文本, 已消费字节数)。"""
    import struct
    import zlib
    out = bytearray()
    i, n = 0, len(buf)
    while i + 18 <= n:
        if buf[i:i + 2] != b"\x1f\x8b":
            i += 1
            continue
        xlen = struct.unpack_from("<H", buf, i + 10)[0]
        bsize = None
        j, endx = i + 12, i + 12 + xlen
        while j + 4 <= endx and j + 4 <= n:
            si1, si2, slen = buf[j], buf[j + 1], struct.unpack_from("<H", buf, j + 2)[0]
            if si1 == 66 and si2 == 67 and slen == 2:
                bsize = struct.unpack_from("<H", buf, j + 4)[0] + 1
                break
            j += 4 + slen
        if bsize is None or i + bsize > n:
            break
        try:
            out.extend(zlib.decompress(buf[i:i + bsize], 31))
        except zlib.error:
            break
        i += bsize
    return out.decode("utf-8", errors="replace"), i


_STUDY_CACHE: Dict[str, List[str]] = {}


def _forward_study(abbr: str) -> str | None:
    """正向分析实际使用的 study accession（保证正反两臂同源）。"""
    f = os.path.join(DATA, f"ivs_{abbr}.csv")
    if not os.path.exists(f):
        return None
    try:
        d = pd.read_csv(f)
        if "study" in d.columns and len(d):
            return str(d["study"].mode().iloc[0])
    except Exception:
        pass
    return None


def _studies_for(efo: str) -> List[str]:
    """EFO -> 该性状下所有 EBI study accession（带进程内缓存）。"""
    if efo not in _STUDY_CACHE:
        try:
            _STUDY_CACHE[efo] = gio.ebi_trait_studies(efo) or []
        except Exception:
            _STUDY_CACHE[efo] = []
    return _STUDY_CACHE[efo]


class RateLimited(RuntimeError):
    """EBI 连续 429：已进入限流惩罚期，应中止而不是继续打接口。"""


def _get_pos_all_studies(chrom, pos: int, pause: float = 2.5) -> List[dict]:
    """
    取 chr:pos 上**所有 study** 的关联（一次请求覆盖全部免疫性状）。

    设计要点
    --------
    * EBI summary-statistics API 会**忽略** ``variant_id`` 参数（实测返回 study
      的前 N 条而非目标变异），只有 ``/chromosomes/{chr}/associations``
      配 ``bp_lower``/``bp_upper`` 才是可靠的定点查询；坐标为 GRCh38，
      与 FinnGen R11 一致。
    * 不加 ``study_accession`` 过滤，一次拿回该位点上的全部 study，
      于是 11 个免疫性状共用同一次请求：请求数从 42×11 降到 42。
    * EBI 有较严格的速率限制，串行 + 固定间隔 + 指数退避；连续 429 直接抛
      ``RateLimited`` 中止，避免把封禁期越拖越长。
    """
    url = f"{gio.EBI_SS}/chromosomes/{chrom}/associations"
    params = {"bp_lower": int(pos), "bp_upper": int(pos), "size": 1000}
    delay = pause
    for attempt in range(5):
        try:
            r = gio._SESSION.get(url, params=params, headers=gio.HEADERS, timeout=120)
        except Exception:
            time.sleep(delay)
            delay *= 2
            continue
        if r.status_code == 200:
            a = r.json().get("_embedded", {}).get("associations", {})
            recs = list(a.values()) if isinstance(a, dict) else list(a)
            time.sleep(pause)
            return recs
        if r.status_code == 404:
            time.sleep(pause)
            return []
        if r.status_code == 429:
            time.sleep(delay)
            delay *= 2
            continue
        time.sleep(pause)
        return []
    raise RateLimited(f"EBI 连续 429（chr{chrom}:{pos}）")


def collect_positions(iv: pd.DataFrame, cache_key: str,
                      pause: float = 2.5) -> Dict[str, Dict[str, dict]]:
    """
    对 IV 列表逐位点抓取“全 study 关联”，落盘缓存。
    返回 {rsid: {study_accession: 规范化关联}}。
    """
    cpath = os.path.join(DATA, "cache", f"revpos_{cache_key}.json")
    store: Dict[str, Dict[str, dict]] = {}
    if os.path.exists(cpath):
        with open(cpath, "r", encoding="utf-8") as fh:
            store = json.load(fh)
    todo = [(r["rsid"], r["chr"], int(r["pos"]))
            for _, r in iv.iterrows() if r["rsid"] not in store]
    if todo:
        log(f"    位点级抓取：{len(todo)} 个待查（已缓存 {len(store)}）")
    for i, (rs, ch, ps) in enumerate(todo, 1):
        try:
            recs = _get_pos_all_studies(ch, ps, pause)
        except RateLimited as ex:
            log(f"    !! {ex}；已抓 {i - 1}/{len(todo)}，中止并保留缓存")
            break
        bystudy: Dict[str, dict] = {}
        for rec in recs:
            nrm = gio._norm_assoc(rec)
            if nrm and int(nrm["pos"]) == ps and rec.get("study_accession"):
                bystudy[rec["study_accession"]] = nrm
        store[rs] = bystudy
        if i % 10 == 0 or i == len(todo):
            log(f"      {i}/{len(todo)} 位点")
            with open(cpath, "w", encoding="utf-8") as fh:
                json.dump(store, fh)
    with open(cpath, "w", encoding="utf-8") as fh:
        json.dump(store, fh)
    return store


def outcome_lookup_ebi(iv: pd.DataFrame, efo: str, workers: int = 4,
                       cache_key: str | None = None,
                       store: Dict[str, Dict[str, dict]] | None = None,
                       abbr: str | None = None
                       ) -> tuple[Dict[str, dict], str | None]:
    """
    从位点级缓存里取某免疫性状的效应量。

    study 选择策略（依次）：
      1) 正向分析同一性状实际使用的 study —— 保证正反两臂数据同源；
      2) 该 EFO 名下命中位点数最多的 study —— 覆盖度优先。
    返回 ({rsid: 关联}, 实际使用的 study)。
    """
    if store is None:
        return {}, None
    cands = set(_studies_for(efo))
    fwd = _forward_study(abbr) if abbr else None
    if fwd:
        cands.add(fwd)
    if not cands:
        return {}, None
    counts: Dict[str, int] = {}
    for v in store.values():
        for acc in v:
            if acc in cands:
                counts[acc] = counts.get(acc, 0) + 1
    if not counts:
        return {}, None
    if fwd and counts.get(fwd, 0) > 0:
        acc = fwd
    else:
        acc = max(counts, key=counts.get)
    return {rs: v[acc] for rs, v in store.items() if acc in v}, acc


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pheno", nargs="*",
                    default=["M13_OSTEOMYELITIS", "M13_DISCITIS"])
    ap.add_argument("--pthresh", type=float, default=5e-6)
    ap.add_argument("--pause", type=float, default=2.5,
                    help="EBI 请求间隔秒数（防 429）")
    args = ap.parse_args()

    os.makedirs(os.path.join(RESULTS, "loo"), exist_ok=True)
    outmap = {o["pheno"]: o for o in mp.OUTCOMES}
    allrows: List[dict] = []
    ivlog: List[dict] = []

    for ph in args.pheno:
        meta = outmap.get(ph, {"abbr": ph, "label": ph, "ncase": None})
        sig = scan_significant(ph, args.pthresh)
        if sig.empty:
            log(f"  {ph}: 无 P<{args.pthresh:g} 变异，跳过")
            ivlog.append({"pheno": ph, "pthresh": args.pthresh, "n_raw": 0,
                          "n_iv": 0, "note": "no variant below threshold"})
            continue
        sig = sig.rename(columns={"sebeta": "se", "rsids": "rsid"})
        sig["rsid"] = sig["rsid"].astype(str).str.split(",").str[0]
        sig = sig[sig["rsid"].str.startswith("rs").fillna(False).astype(bool)]
        sig["effect_allele"] = sig["alt"].str.upper()
        sig["other_allele"] = sig["ref"].str.upper()
        sig["eaf"] = sig["af_alt"]
        iv = mp.distance_prune(sig)
        iv["F"] = (iv["beta"] / iv["se"]) ** 2
        iv = iv[iv["F"] > mp.F_MIN]
        log(f"  {ph}: {len(sig)} 个显著变异 -> {len(iv)} 个独立 IV (F>10)")
        iv.to_csv(os.path.join(DATA, f"ivs_reverse_{meta['abbr']}.csv"), index=False)
        ivlog.append({"pheno": ph, "abbr": meta["abbr"], "pthresh": args.pthresh,
                      "n_raw": len(sig), "n_iv": len(iv),
                      "F_mean": float(iv["F"].mean()) if len(iv) else np.nan})
        if iv.empty:
            continue

        # 一次抓取覆盖全部免疫性状：请求数 = IV 数，而非 IV 数 × 性状数
        store = collect_positions(iv, f"{meta['abbr']}_{args.pthresh:g}",
                                  pause=args.pause)
        ncov = sum(1 for v in store.values() if v)
        log(f"  {meta['abbr']}: 位点级缓存 {len(store)} 个，其中 {ncov} 个在 EBI 有数据")

        for e in mp.EXPOSURES:      # 反向：免疫性状作结局
            hits, acc = outcome_lookup_ebi(iv, e["efo"], store=store,
                                           abbr=e["abbr"])
            log(f"    {meta['abbr']} -> {e['abbr']}: {len(hits)}/{len(iv)} IV "
                f"在结局中命中 (study={acc})")
            if len(hits) < 1:
                continue
            rows = []
            for _, r in iv.iterrows():
                h = hits.get(r["rsid"])
                if not h:
                    continue
                ea, oa = r["effect_allele"], r["other_allele"]
                by, byse = h["beta"], h["se"]
                if h["effect_allele"] == ea and h["other_allele"] == oa:
                    action = "same"
                elif h["effect_allele"] == oa and h["other_allele"] == ea:
                    by, action = -by, "flip"
                else:
                    continue
                rows.append({"SNP": r["rsid"], "beta_exposure": r["beta"],
                             "se_exposure": r["se"], "pval_exposure": r["pval"],
                             "F": r["F"], "beta_outcome": by, "se_outcome": byse,
                             "pval_outcome": h["pval"], "harmonise_action": action,
                             "gene": r.get("gene")})
            h_df = pd.DataFrame(rows)
            if h_df.empty:
                continue
            h_df.to_csv(os.path.join(RESULTS, "harmonised",
                                     f"reverse_{meta['abbr']}__{e['abbr']}_harmonised.csv"),
                        index=False)
            res = mrm.run_all(h_df["beta_exposure"].values, h_df["se_exposure"].values,
                              h_df["beta_outcome"].values, h_df["se_outcome"].values,
                              h_df["SNP"].tolist())
            log(f"    {meta['abbr']} -> {e['abbr']}: {len(h_df)} IV")
            for rr in res["results"]:
                d = rr.as_dict(binary_outcome=False)
                d.update({"direction": "reverse",
                          "exposure": meta["label"], "exposure_abbr": meta["abbr"],
                          "outcome": e["name"], "outcome_abbr": e["abbr"],
                          "F_mean": float(np.mean(res["F"])),
                          "F_min": float(np.min(res["F"])),
                          "Q": res["heterogeneity"]["Q"],
                          "Q_pval": res["heterogeneity"]["Q_pval"],
                          "I2": res["heterogeneity"]["I2"],
                          "egger_intercept": res["pleiotropy"]["intercept"],
                          "egger_intercept_pval": res["pleiotropy"]["intercept_pval"],
                          "iv_pthreshold": args.pthresh,
                          "outcome_study": acc})
                allrows.append(d)
            # 增量落盘：即使后续被限流/中断，已完成部分也不丢
            pd.DataFrame(allrows).to_csv(
                os.path.join(RESULTS, "MR_results_reverse.csv"), index=False)
            if len(h_df) >= 3:
                try:
                    tag = f"{meta['abbr']}__{e['abbr']}"
                    mp.plot_scatter(h_df, res, f"{meta['abbr']} -> {e['abbr']} ({len(h_df)} IVs)",
                                    os.path.join(mp.FIGDIR, f"reverse_{tag}_scatter.png"))
                    mp.plot_loo(res["loo"], f"Leave-one-out: {meta['abbr']} -> {e['abbr']}",
                                os.path.join(mp.FIGDIR, f"reverse_{tag}_loo.png"))
                except Exception as ex:
                    log(f"    !! 作图失败 {ex}")

    pd.DataFrame(ivlog).to_csv(os.path.join(RESULTS, "reverse_MR_IV_discovery.csv"),
                               index=False)
    if allrows:
        df = pd.DataFrame(allrows)
        df.to_csv(os.path.join(RESULTS, "MR_results_reverse.csv"), index=False)
        log(f"反向结果已写出：{len(df)} 行")
        prim = df[df["method"].str.startswith("IVW (random")]
        for _, r in prim[prim["pval"] < 0.05].iterrows():
            log(f"   * {r['exposure_abbr']} -> {r['outcome_abbr']}: "
                f"b={r['b']:.3f} P={r['pval']:.3g} nSNP={r['nsnp']}")
    else:
        log("反向 MR 未产出结果（IV 不足）")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
