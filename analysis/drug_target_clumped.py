# -*- coding: utf-8 -*-
"""Re-run drug-target (cis-pQTL) MR with formally LD-clumped IV sets
(1000 Genomes Phase 3 EUR, r2<0.001 via LDlink LDmatrix).
Compares against the distance-pruned results in results/drug_target.tsv.
"""
import sys, os, json, csv
import pandas as pd
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)
import importlib.util
spec = importlib.util.spec_from_file_location("mr_pipeline", os.path.join(BASE, "mr_pipeline.py"))
mp = importlib.util.module_from_spec(spec); spec.loader.exec_module(mp)
sys.path.insert(0, BASE)
import mr_methods as mrm

DATA = os.path.join(BASE, '..', 'data')
RESULTS = os.path.join(BASE, '..', 'results')

def clumped_ivs(trait):
    """Load LD-clumped kept cis IVs as a DataFrame."""
    jp = os.path.join(RESULTS, f'clump_{trait}_cis_EUR_r2_0.001.json')
    with open(jp, encoding='utf-8') as f:
        d = json.load(f)
    kept = {k['rsid'] for k in d['kept']}
    rows = list(csv.DictReader(open(os.path.join(DATA, f'ivs_cis_{trait}.csv'))))
    df = pd.DataFrame([r for r in rows if r['rsid'] in kept])
    for col in ('pos', 'beta', 'se', 'pval', 'eaf', 'F'):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def main():
    import time
    outcomes = {o['abbr']: o for o in mp.OUTCOMES if o['abbr'] in ('OM', 'DISCITIS')}
    out_rows = []
    for trait in ('IL6R', 'CD40'):
        iv = clumped_ivs(trait)
        if iv.empty:
            print(f'{trait}: no clumped IVs')
            continue
        print(f'--- {trait}: {len(iv)} clumped cis IVs: {iv["rsid"].tolist()}')
        for oabbr, o in outcomes.items():
            # FinnGen remote read with retries (network resets observed)
            h = None
            for attempt in range(4):
                try:
                    fs = mp.gio.FinnGenSumstats(o['pheno'])
                    h_all = mp.harmonise(iv, fs, o['label'])
                    h = h_all[h_all['status'] == 'kept']
                    break
                except Exception as e:
                    wait = 5 * (attempt + 1)
                    print(f'  [retry {attempt+1}] {o["pheno"]}: {str(e)[:80]} waiting {wait}s', flush=True)
                    time.sleep(wait)
            if h is None or len(h) < 1:
                print(f'  {trait}->{oabbr}: {0 if h is None else len(h)} kept after retries'); continue
            bx = h['beta_exposure'].values.astype(float)
            bxse = h['se_exposure'].values.astype(float)
            by = h['beta_outcome'].values.astype(float)
            byse = h['se_outcome'].values.astype(float)
            snps = h['SNP'].tolist()
            if len(snps) == 1:
                # Wald ratio
                b = by[0] / bx[0]
                se = abs(byse[0] / bx[0])
                res = mrm._finish('Wald ratio', 1, b, se, 'single-SNP Wald ratio')
                F = (bx[0] / bxse[0]) ** 2
                d = res.as_dict(binary_outcome=True)
                d.update({'target': trait, 'outcome': oabbr, 'nsnp_cis': 1, 'F_min': F,
                          'clumped': True})
                out_rows.append(d)
                print(f'  {trait}->{oabbr}: Wald nSNP=1 OR={np.exp(b):.3f} P={res.pval:.3g} '
                      f'(95% CI {np.exp(res.lo95):.3f}-{np.exp(res.hi95):.3f})')
            else:
                res = mrm.run_all(bx, bxse, by, byse, snps)
                for r in res['results']:
                    d = r.as_dict(binary_outcome=True)
                    d.update({'target': trait, 'outcome': oabbr, 'nsnp_cis': len(snps),
                              'F_min': float(np.min(res['F'])), 'clumped': True})
                    out_rows.append(d)
                ivw = res['results'][0]
                print(f'  {trait}->{oabbr}: nSNP={len(snps)} IVW OR={np.exp(ivw.b):.3f} '
                      f'P={ivw.pval:.3g} (95% CI {np.exp(ivw.lo95):.3f}-{np.exp(ivw.hi95):.3f})')
    if out_rows:
        out = os.path.join(RESULTS, 'drug_target_clumped_LDlink.tsv')
        pd.DataFrame(out_rows).to_csv(out, sep='\t', index=False)
        print(f'[written] {out} ({len(out_rows)} rows)')

if __name__ == '__main__':
    main()
