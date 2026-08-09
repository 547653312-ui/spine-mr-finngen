"""
wait_and_run_reverse.py — 等 EBI 限流解除后自动跑反向 MR。

背景：反向 MR 的结局端（免疫性状）只能走 EBI summary-statistics API，
而该 API 有速率限制；密集查询后会进入一段时间的 429 惩罚期。
本脚本每 3 分钟做一次**轻量探测**（单个请求，不构成滥用），
探测通过后立即调用 reverse_mr.main()，失败则最多等 40 分钟后放弃并如实记录。
"""

from __future__ import annotations

import subprocess
import sys
import time

import gwas_io as gio

PROBE = f"{gio.EBI_SS}/chromosomes/1/associations"
PARAMS = {"bp_lower": 159712443, "bp_upper": 159712443, "size": 5}
MAX_WAIT_MIN = 40
PROBE_EVERY_S = 180


def log(m: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def probe() -> int:
    try:
        r = gio._SESSION.get(PROBE, params=PARAMS, headers=gio.HEADERS, timeout=60)
        return r.status_code
    except Exception as e:
        log(f"  探测异常：{e}")
        return -1


def main() -> int:
    t0 = time.time()
    while (time.time() - t0) / 60 < MAX_WAIT_MIN:
        code = probe()
        log(f"探测 EBI -> HTTP {code}")
        if code == 200:
            log("限流已解除，开始跑反向 MR")
            cmd = [sys.executable, "-u", "reverse_mr.py",
                   "--pheno", "M13_OSTEOMYELITIS", "M13_DISCITIS",
                   "--pthresh", "5e-6", "--pause", "3.0"]
            return subprocess.call(cmd)
        time.sleep(PROBE_EVERY_S)
    log(f"等待 {MAX_WAIT_MIN} 分钟仍被限流，放弃本次反向 MR（已在 LIMITATIONS.md 记录）")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
