"""
gwas_io.py — 开放数据源取数层（无需任何 token）

两个数据源：
1) EBI GWAS Catalog Summary Statistics API  -> 暴露端（免疫/细胞因子性状 IV）
2) FinnGen R11 公开 bgzip+tabix 汇总统计    -> 结局端（脊柱感染表型），
   通过 HTTP Range + 纯 Python tabix 解析按位点精准取数，避免下载 ~800MB 全文件。

背景：OpenGWAS(IEU) 自 2024-05-01 起对 /api/* 强制 JWT，本机无 token（401），
因此改用上述两个完全开放的等价数据源。详见 results/API_NOTES.md。
"""

from __future__ import annotations

import gzip
import io
import json
import os
import struct
import threading
import time
import zlib
from typing import Dict, Iterable, List, Optional, Tuple

import requests

# --------------------------------------------------------------------------
# 常量
# --------------------------------------------------------------------------
EBI_SS = "https://www.ebi.ac.uk/gwas/summary-statistics/api"
FINNGEN_BUCKET = "https://storage.googleapis.com/finngen-public-data-r11/summary_stats"
FINNGEN_API = "https://r11.finngen.fi/api"
OPENGWAS = "https://api.opengwas.io/api"

HEADERS = {"User-Agent": "spine-MR-pipeline/1.0 (research; python-requests)"}
CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "data", "cache")

# 全局 Session：复用 TLS 连接，是远端 tabix 随机访问提速的关键
_SESSION = requests.Session()
_SESSION.mount("https://", requests.adapters.HTTPAdapter(
    pool_connections=32, pool_maxsize=32, max_retries=0))


def _cache_path(name: str) -> str:
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, name)


def http_get(url: str, params: dict | None = None, tries: int = 3,
             timeout: int = 90, headers: dict | None = None,
             stream: bool = False) -> requests.Response:
    """带重试的 GET（复用全局 Session）。"""
    h = dict(HEADERS)
    if headers:
        h.update(headers)
    last = None
    for i in range(tries):
        try:
            r = _SESSION.get(url, params=params, headers=h, timeout=timeout, stream=stream)
            if r.status_code in (429, 500, 502, 503, 504):
                last = RuntimeError(f"HTTP {r.status_code} for {r.url}")
                time.sleep(2 * (i + 1))
                continue
            return r
        except requests.RequestException as e:  # 网络层错误
            last = e
            time.sleep(2 * (i + 1))
    raise RuntimeError(f"GET failed after {tries} tries: {url} :: {last}")


# ==========================================================================
# 1. OpenGWAS 连通性探测（用于 API_NOTES 记录，不参与主分析）
# ==========================================================================
def probe_opengwas(ids: Iterable[str]) -> List[dict]:
    """探测 IEU/OpenGWAS 可达性，记录精确错误。"""
    out = []
    # 旧端点
    for url in [
        "https://gwas.mrcieu.ac.uk/api/gwasinfo/ieu-b-4975",
        f"{OPENGWAS}/status",
    ]:
        try:
            r = requests.get(url, headers=HEADERS, timeout=30, allow_redirects=False)
            out.append({"endpoint": url, "status": r.status_code,
                        "body": r.text[:300], "note": "legacy/status"})
        except Exception as e:
            out.append({"endpoint": url, "status": None, "body": f"EXC {e}", "note": "legacy/status"})

    for gid in ids:
        url = f"{OPENGWAS}/gwasinfo/{gid}"
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
            out.append({"endpoint": url, "status": r.status_code, "body": r.text[:300],
                        "note": "no-token"})
        except Exception as e:
            out.append({"endpoint": url, "status": None, "body": f"EXC {e}", "note": "no-token"})
    # 需要 token 的列表端点
    try:
        r = requests.get(f"{OPENGWAS}/gwasinfo", headers=HEADERS, timeout=30)
        out.append({"endpoint": f"{OPENGWAS}/gwasinfo", "status": r.status_code,
                    "body": r.text[:400], "note": "token-required check"})
    except Exception as e:
        out.append({"endpoint": f"{OPENGWAS}/gwasinfo", "status": None,
                    "body": f"EXC {e}", "note": "token-required check"})
    return out


# ==========================================================================
# 2. EBI GWAS Catalog Summary Statistics —— 暴露端
# ==========================================================================
def ebi_trait_studies(efo: str) -> List[str]:
    """列出某 EFO 性状下有汇总统计的研究 accession。"""
    r = http_get(f"{EBI_SS}/traits/{efo}/studies", params={"size": 200})
    if r.status_code != 200:
        return []
    d = r.json()
    studies = d.get("_embedded", {}).get("studies", [])
    if isinstance(studies, dict):
        studies = list(studies.values())
    return [s.get("study_accession") for s in studies if s.get("study_accession")]


def _norm_assoc(rec: dict) -> Optional[dict]:
    """统一 EBI association 记录字段。"""
    try:
        beta = rec.get("beta")
        orv = rec.get("odds_ratio")
        if beta is None and orv:
            import math
            beta = math.log(float(orv))
        if beta is None:
            return None
        p = rec.get("p_value")
        if p is None:
            return None
        p = float(p)
        beta = float(beta)
        # EBI 汇总统计 API 不直接给 SE，用 p 与 beta 反推（双侧正态）
        se = rec.get("standard_error")
        if se is None:
            se = se_from_beta_p(beta, p)
        if se is None or se <= 0:
            return None
        return {
            "rsid": rec.get("variant_id"),
            "chr": str(rec.get("chromosome")),
            "pos": int(rec.get("base_pair_location")) if rec.get("base_pair_location") else None,
            "effect_allele": (rec.get("effect_allele") or "").upper(),
            "other_allele": (rec.get("other_allele") or "").upper(),
            "eaf": rec.get("effect_allele_frequency"),
            "beta": beta,
            "se": float(se),
            "pval": p,
            "study": rec.get("study_accession"),
        }
    except Exception:
        return None


def se_from_beta_p(beta: float, p: float) -> Optional[float]:
    """由 beta 与双侧 P 反推 SE：|beta| / z，z = Phi^-1(1 - p/2)。"""
    from scipy.stats import norm
    if p <= 0:
        p = 1e-300
    if p >= 1:
        return None
    z = norm.isf(p / 2.0)
    if z <= 0:
        return None
    return abs(beta) / z


def ebi_variant_associations(rsid: str, efo: str | None = None,
                             tries: int = 4) -> List[dict]:
    """
    取某 SNP (rsid) 在 GWAS Catalog 的关联记录（用于 MVMR 变异性水平取数）。
    可选按 EFO 过滤。返回与 _norm_assoc 相同结构的列表。
    """
    url = f"{EBI_SS}/variants/{rsid}/associations"
    out: List[dict] = []
    try:
        r = http_get(url, params={"size": 500}, tries=tries)
    except Exception:
        return out
    if r.status_code != 200:
        return out
    try:
        d = r.json()
    except Exception:
        return out
    assoc = d.get("_embedded", {}).get("associations", {})
    recs = list(assoc.values()) if isinstance(assoc, dict) else list(assoc)
    for rec in recs:
        if efo and rec.get("studyAccession") is None:
            # 无法按 EFO 可靠过滤时，用 trait 名称近似匹配
            pass
        n = _norm_assoc(rec)
        if n and n["rsid"] and n["pos"]:
            out.append(n)
    return out


def ebi_trait_associations(efo: str, p_upper: float = 5e-8,
                           max_records: int = 4000,
                           study: str | None = None) -> List[dict]:
    """
    拉取某性状（或指定研究）中 P < p_upper 的关联。
    EBI API 按页返回，size 上限 500。
    """
    base = f"{EBI_SS}/traits/{efo}/associations"
    if study:
        base = f"{EBI_SS}/studies/{study}/associations"
    out: List[dict] = []
    start = 0
    size = 500
    while len(out) < max_records:
        r = http_get(base, params={"p_upper": p_upper, "size": size, "start": start})
        if r.status_code != 200:
            break
        d = r.json()
        assoc = d.get("_embedded", {}).get("associations", {})
        if isinstance(assoc, dict):
            recs = list(assoc.values())
        else:
            recs = list(assoc)
        if not recs:
            break
        for rec in recs:
            n = _norm_assoc(rec)
            if n and n["rsid"] and n["pos"]:
                out.append(n)
        if len(recs) < size:
            break
        start += size
    return out


# ==========================================================================
# 3. FinnGen R11 —— 结局端（BGZF + tabix over HTTP Range）
# ==========================================================================
class BGZFRemote:
    """对远端 bgzip 文件做 HTTP Range 读取 + BGZF 解块。"""

    def __init__(self, url: str):
        self.url = url
        last = None
        for i in range(4):
            try:
                r = _SESSION.head(url, headers=HEADERS, timeout=60)
                r.raise_for_status()
                self.size = int(r.headers["Content-Length"])
                return
            except Exception as e:      # 含 ConnectionReset 等瞬时错误
                last = e
                time.sleep(2 * (i + 1))
        raise RuntimeError(f"HEAD failed for {url}: {last}")

    def read_range(self, start: int, length: int) -> bytes:
        end = min(start + length - 1, self.size - 1)
        h = dict(HEADERS)
        h["Range"] = f"bytes={start}-{end}"
        r = http_get(self.url, headers=h, timeout=120, tries=4)
        if r.status_code not in (200, 206):
            raise RuntimeError(f"Range request failed HTTP {r.status_code}")
        return r.content

    @staticmethod
    def inflate_blocks(buf: bytes) -> bytes:
        """逐个 BGZF 块解压；末尾不完整的块直接丢弃。"""
        out = bytearray()
        i = 0
        n = len(buf)
        while i + 18 <= n:
            if buf[i:i + 2] != b"\x1f\x8b":
                break
            xlen = struct.unpack_from("<H", buf, i + 10)[0]
            # 在 extra field 中找 BC 子字段拿 BSIZE
            bsize = None
            j = i + 12
            endx = j + xlen
            while j + 4 <= endx:
                si1, si2, slen = buf[j], buf[j + 1], struct.unpack_from("<H", buf, j + 2)[0]
                if si1 == 66 and si2 == 67 and slen == 2:
                    bsize = struct.unpack_from("<H", buf, j + 4)[0] + 1
                    break
                j += 4 + slen
            if bsize is None:
                break
            if i + bsize > n:
                break  # 截断块
            block = buf[i:i + bsize]
            try:
                out.extend(zlib.decompress(block, 31))
            except zlib.error:
                break
            i += bsize
        return bytes(out)


class TabixIndex:
    """纯 Python tabix (.tbi) 索引解析器。"""

    def __init__(self, raw: bytes):
        data = gzip.decompress(raw)
        if data[:4] != b"TBI\x01":
            raise ValueError("not a tabix index")
        off = 4
        (self.n_ref, self.fmt, self.col_seq, self.col_beg, self.col_end,
         self.meta, self.skip, l_nm) = struct.unpack_from("<8i", data, off)
        off += 32
        names_blob = data[off:off + l_nm]
        off += l_nm
        self.names = [n.decode() for n in names_blob.split(b"\x00") if n]
        self.name2idx = {n: i for i, n in enumerate(self.names)}
        self.bins: List[Dict[int, List[Tuple[int, int]]]] = []
        self.intv: List[List[int]] = []
        for _ in range(self.n_ref):
            (n_bin,) = struct.unpack_from("<i", data, off)
            off += 4
            bd: Dict[int, List[Tuple[int, int]]] = {}
            for _ in range(n_bin):
                bin_id, n_chunk = struct.unpack_from("<Ii", data, off)
                off += 8
                chunks = []
                for _ in range(n_chunk):
                    cb, ce = struct.unpack_from("<QQ", data, off)
                    off += 16
                    chunks.append((cb, ce))
                bd[bin_id] = chunks
            (n_intv,) = struct.unpack_from("<i", data, off)
            off += 4
            iv = list(struct.unpack_from(f"<{n_intv}Q", data, off)) if n_intv else []
            off += 8 * n_intv
            self.bins.append(bd)
            self.intv.append(iv)

    def linear_offset(self, chrom: str, pos: int) -> Optional[int]:
        """返回覆盖该位置的最小 coffset（线性索引）。pos 为 1-based。"""
        idx = self.name2idx.get(chrom)
        if idx is None:
            # FinnGen 用 '1'..'23'，容错 'chr1'
            idx = self.name2idx.get(chrom.replace("chr", ""))
            if idx is None:
                idx = self.name2idx.get("chr" + chrom)
        if idx is None:
            return None
        iv = self.intv[idx]
        if not iv:
            return None
        k = max(0, (pos - 1) >> 14)
        if k >= len(iv):
            k = len(iv) - 1
        vo = iv[k]
        return vo >> 16  # 压缩文件偏移


class FinnGenSumstats:
    """FinnGen R11 单个表型的远端随机访问。"""

    COLS = ["chrom", "pos", "ref", "alt", "rsids", "nearest_genes",
            "pval", "mlogp", "beta", "sebeta", "af_alt",
            "af_alt_cases", "af_alt_controls"]

    def __init__(self, phenocode: str, release: str = "R11"):
        self.pheno = phenocode
        self.url = f"{FINNGEN_BUCKET}/finngen_{release}_{phenocode}.gz"
        self.tbi_url = self.url + ".tbi"
        self.bgzf = BGZFRemote(self.url)
        tbi_cache = _cache_path(f"finngen_{release}_{phenocode}.gz.tbi")
        if os.path.exists(tbi_cache) and os.path.getsize(tbi_cache) > 1000:
            raw = open(tbi_cache, "rb").read()
        else:
            r = http_get(self.tbi_url, timeout=180)
            r.raise_for_status()
            raw = r.content
            with open(tbi_cache, "wb") as f:
                f.write(raw)
        self.index = TabixIndex(raw)
        self._blockcache: Dict[tuple, str] = {}
        self._lock = threading.Lock()

    def _fetch_text(self, coffset: int, nbytes: int = 131072) -> str:
        key = (coffset, nbytes)
        with self._lock:
            hit = self._blockcache.get(key)
        if hit is not None:
            return hit
        raw = self.bgzf.read_range(coffset, nbytes)
        txt = BGZFRemote.inflate_blocks(raw).decode("utf-8", errors="replace")
        with self._lock:
            if len(self._blockcache) > 400:
                self._blockcache.clear()
            self._blockcache[key] = txt
        return txt

    def query(self, chrom: str, pos: int, ref: str | None = None,
              alt: str | None = None, window: int = 0) -> List[dict]:
        """按 chrom:pos 精确查找（GRCh38）。返回匹配行。"""
        co = self.index.linear_offset(str(chrom), pos)
        if co is None:
            return []
        hits: List[dict] = []
        nbytes = 131072
        for attempt in range(3):
            txt = self._fetch_text(co, nbytes)
            lines = txt.split("\n")
            passed = False
            for ln in lines[1:-1] if attempt == 0 else lines[1:-1]:
                if not ln or ln.startswith("#"):
                    continue
                f = ln.split("\t")
                if len(f) < 10:
                    continue
                try:
                    p = int(f[1])
                except ValueError:
                    continue
                if f[0] != str(chrom) and f[0] != str(chrom).replace("chr", ""):
                    continue
                if abs(p - pos) <= window:
                    rec = dict(zip(self.COLS, f))
                    hits.append(rec)
                if p > pos + max(window, 0) + 1000:
                    passed = True
                    break
            if hits or passed:
                break
            nbytes *= 4  # 该窗口内变异过多，扩大读取
        # 若给定等位基因，优先精确匹配
        if hits and ref and alt:
            exact = [h for h in hits
                     if {h["ref"].upper(), h["alt"].upper()} == {ref.upper(), alt.upper()}]
            if exact:
                return exact
        return hits


def finngen_pheno_meta(phenocode: str) -> dict:
    r = http_get(f"{FINNGEN_API}/pheno/{phenocode}", timeout=60)
    if r.status_code == 200:
        try:
            return r.json()
        except json.JSONDecodeError:
            return {}
    return {}
