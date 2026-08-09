# -*- coding: utf-8 -*-
"""
Formal LD clumping sensitivity analysis via LDlink LDmatrix API
(1000 Genomes Phase 3 v5a, EUR, r2 threshold 0.001).

For each of the six key traits (IL6R, CD40, CRP, WBC, NEUT, CCL2) the
genome-wide/cis instruments are clumped using pairwise LD from the LDlink
LDmatrix endpoint, then compared with the distance-pruned instrument set.

Usage:
    python ld_clump_ldlink.py --token YOUR_LDlink_TOKEN [--build grch38]
"""
import argparse, csv, json, os, sys, time
import urllib.request, urllib.error

LDMATRIX_URL = "https://ldlink.nih.gov/LDlinkRest/ldmatrix"
R2_THRESHOLD = 0.001
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data')
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results')
TRAITS = ['IL6R', 'CD40', 'CRP', 'WBC', 'NEUT', 'CCL2']


def load_ivs(trait):
    rows = []
    with open(os.path.join(DATA, f'ivs_{trait}.csv'), encoding='utf-8') as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def ldmatrix(token, snps, build='grch38'):
    """Query LDlink LDmatrix (GET) for pairwise r2 among snps (EUR).
    Returns a dict matrix keyed by rsID -> {rsID: r2}."""
    url = (f"{LDMATRIX_URL}?snps={'%0A'.join(snps)}&pop=EUR&r2_d=r2"
           f"&genome_build={build}&token={token}")
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=120) as r:
                text = r.read().decode('utf-8', errors='replace')
            break
        except urllib.error.HTTPError as e:
            if e.code == 429:
                wait = 15 * (attempt + 1)
                print(f'  [429] rate limited, waiting {wait}s', file=sys.stderr)
                time.sleep(wait)
                continue
            raise
        except urllib.error.URLError as e:
            wait = 10 * (attempt + 1)
            print(f'  [net] {e}, waiting {wait}s', file=sys.stderr)
            time.sleep(wait)
    else:
        raise RuntimeError('LDmatrix failed after retries')
    # LDlink returns TSV: header row RS_number\tsnp1\tsnp2..., then rows snp\tvals...
    lines = [ln for ln in text.strip().split('\n') if ln.strip()]
    header = lines[0].split('\t')
    col_snps = [h for h in header[1:]]
    mat = {}
    for ln in lines[1:]:
        parts = ln.split('\t')
        snp = parts[0]
        mat[snp] = {}
        for k, v in zip(col_snps, parts[1:]):
            mat[snp][k] = v
    # rebuild square matrix in requested order; NA -> None
    M = [[None] * len(snps) for _ in range(len(snps))]
    for i, s1 in enumerate(snps):
        for j, s2 in enumerate(snps):
            if s1 in mat and s2 in mat[s1]:
                v = mat[s1][s2]
                if v not in ('NA', 'nan', ''):
                    try:
                        M[i][j] = float(v)
                    except ValueError:
                        M[i][j] = None
    return M


def clump_from_matrix(snps, matrix, pvals, r2_thresh=R2_THRESHOLD):
    """Standard stepwise clumping: iterate by ascending P, drop r2>threshold."""
    order = sorted(range(len(snps)), key=lambda i: pvals[i])
    selected = []
    removed = set()
    for i in order:
        if i in removed:
            continue
        selected.append(i)
        for j in range(len(snps)):
            if j == i or j in removed:
                continue
            if matrix[i][j] is not None and matrix[i][j] > r2_thresh:
                removed.add(j)
    return selected, sorted(removed)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--token', required=True)
    ap.add_argument('--build', default='grch38', choices=['grch37', 'grch38'])
    args = ap.parse_args()

    os.makedirs(OUT, exist_ok=True)
    summary_rows = []
    for trait in TRAITS:
        ivs = load_ivs(trait)
        snps = [r['rsid'] for r in ivs]
        pvals = [float(r['pval']) for r in ivs]
        print(f'--- {trait}: {len(snps)} IVs ---')
        try:
            mat = ldmatrix(args.token, snps, args.build)
        except Exception as e:
            print(f'  !! LDmatrix failed for {trait}: {e}')
            continue
        selected, removed = clump_from_matrix(snps, mat, pvals)
        n_before = len(snps)
        print(f'  before={n_before} after={len(selected)} removed={len(removed)}')
        for i in removed:
            print(f'    removed: {snps[i]} (p={pvals[i]:.2e})')
        kept = [{'rsid': snps[i], 'chr': ivs[i]['chr'], 'pos': ivs[i]['pos'],
                 'pval': pvals[i]} for i in selected]
        dropped = [{'rsid': snps[i], 'chr': ivs[i]['chr'], 'pos': ivs[i]['pos'],
                    'pval': pvals[i]} for i in removed]
        with open(os.path.join(OUT, f'clump_{trait}_EUR_r2_0.001.json'), 'w',
                  encoding='utf-8') as f:
            json.dump({'trait': trait, 'build': args.build, 'r2_threshold': R2_THRESHOLD,
                       'n_before': n_before, 'n_after': len(selected),
                       'kept': kept, 'dropped': dropped}, f,
                      ensure_ascii=False, indent=1)
        summary_rows.append({'trait': trait, 'n_before': n_before,
                             'n_after': len(selected), 'n_removed': len(removed)})
        time.sleep(2)   # be polite to the API
    print('\n=== SUMMARY ===')
    for r in summary_rows:
        print(r)


if __name__ == '__main__':
    main()
