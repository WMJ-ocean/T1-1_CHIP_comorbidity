# -*- coding: utf-8 -*-
"""73_supp_tables.py —— 生成 Table S7（两步蛋白 MR + 中介）与 S8（单核细胞 DE/GSVA），
并修正正文里对补充表的编号（S6 已被 VEP 注释占用，蛋白 MR 应为 S7，单核细胞为 S8）。"""
import io, csv

BASE = 'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity'
MS = BASE + '/10_manuscript/Manuscript_Leukemia.md'
SUP = BASE + '/10_manuscript/Supplementary_Tables.md'
MR = BASE + '/07_mr_protein'
P2 = BASE + '/07_p2_transcriptome'


def rd(p):
    return list(csv.DictReader(io.open(p, encoding='utf-8'), delimiter='\t'))


# ---------- Table S7 ----------
s1 = {r['symbol']: r for r in rd(MR + '/step1_chip_to_protein.tsv')}
s2 = {(r['symbol'], r['outcome']): r for r in rd(MR + '/step2_protein_to_outcome.tsv')}
med = {(r['symbol'], r['outcome']): r for r in rd(MR + '/mediation.tsv')}
rows = []
for sym, r in s1.items():
    a, sa, pa = float(r['beta_SD_per_logOR']), float(r['se']), float(r['p'])
    k1 = r['n_iv']
    for phe in ('MPN', 'PV'):
        m = s2.get((sym, phe)); md = med.get((sym, phe))
        b = sb = pb = k2 = prop = '—'
        if m:
            b = f"{float(m['beta_logOR_per_SD']):+.3f}"
            sb = f"{float(m['se']):.3f}"
            pb = f"{float(m['p']):.3g}"
            k2 = m['k']
        if md:
            prop = f"{float(md['prop_mediated'])*100:.1f}%"
        rows.append(f'| {sym} | {k1} | {a:+.3f} | {sa:.3f} | {pa:.3g} | {phe} | {k2} | {b} | {sb} | {pb} | {prop} |')

S7 = ['', '## Table S7. Two-step Mendelian randomisation of circulating inflammatory proteins',
      '',
      '| Protein | nSNP (step 1) | Step 1 beta | SE | P | Outcome | nSNP (step 2) | Step 2 beta | SE | P | Mediated |',
      '|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|'] + rows + [
      '',
      'Step 1: CHIP germline liability (20 instruments) to protein, inverse-variance weighted on Wald ratios. '
      'Step 2: cis protein instruments (P < 5 x 10^-8, within +/-1 Mb of the gene, pruned at 250 kb) to FinnGen '
      'outcome. Step 1 beta is in protein standard deviations per log-odds unit of CHIP liability; step 2 beta is '
      'in log-odds of the outcome per protein standard deviation. Mediated proportion = (a x b)/c with c the total '
      'CHIP-to-outcome IVW estimate (MPN beta = 1.075, PV beta = 1.248); reported as a percentage. '
      'No protein with a genome-wide significant cis instrument showed an effect on either outcome that survived '
      'correction, and the largest mediated proportion was 1.9%. Dash, not estimable (no cis instrument).',
      '',
      'Step 1 note: TNF, IFN-gamma, IL-17A and LIF have no genome-wide significant pQTL in this dataset (smallest '
      'P = 1.2 x 10^-7, 6.0 x 10^-8, 7.2 x 10^-8 and 2.9 x 10^-7 respectively), so they contribute to step 1 only.']

# ---------- Table S8 ----------
de = []
for k, lab in [('MF_NBM', 'Myelofibrosis vs normal'), ('PV_NBM', 'Polycythaemia vera vs normal'),
               ('ET_NBM', 'Essential thrombocythaemia vs normal'), ('MF_PV', 'Myelofibrosis vs polycythaemia vera')]:
    t = rd(P2 + f'/DE_{k}.tsv')
    up = sum(1 for r in t if float(r['adj.P.Val']) < 0.05 and float(r['logFC']) > 0)
    dn = sum(1 for r in t if float(r['adj.P.Val']) < 0.05 and float(r['logFC']) < 0)
    de.append(f'| {lab} | {len(t):,} | {up+dn:,} | {up:,} | {dn:,} |')

sig = rd(P2 + '/signature_score_by_group.tsv')
sigrows = [f"| {r['group']} | {r['n']} | {float(r['median']):+.3f} | {float(r['mean']):+.3f} | {float(r['sd']):.3f} |" for r in sig]

gs = rd(P2 + '/gsva_hallmark_MF_vs_NBM.tsv')
gs.sort(key=lambda r: float(r['P.Value']))
top = [f"| {r['pathway'].replace('HALLMARK_','')} | {float(r['logFC']):+.3f} | {float(r['P.Value']):.2e} | {float(r['adj.P.Val']):.2e} |"
       for r in gs[:10]]

S8 = ['', '## Table S8. CD14+ monocyte transcriptome across the MPN spectrum (GSE283711)',
      '',
      '**(a) Differential expression (limma-voom, TMM, filterByExpr; 24,698 genes retained).**',
      '',
      '| Contrast | Genes tested | FDR < 0.05 | Up | Down |',
      '|---|---:|---:|---:|---:|'] + de + [
      '',
      '**(b) CHIP-direction expression score by group** (52 of 69 candidate genes; weight = sign of the CHIP-to-expression effect). '
      'Kruskal-Wallis P = 0.86; Spearman trend across normal, essential thrombocythaemia, polycythaemia vera, myelofibrosis '
      'rho = -0.01, P = 0.93; per-gene directional concordance with myelofibrosis vs normal 48.1%.',
      '',
      '| Group | n | Median | Mean | SD |',
      '|---|---:|---:|---:|---:|'] + sigrows + [
      '',
      '**(c) GSVA, MSigDB hallmark (top 10 of 50 by P; 18 sets at FDR < 0.05 for myelofibrosis vs normal).**',
      '',
      '| Pathway | logFC | P | FDR |',
      '|---|---:|---:|---:|'] + top + [
      '',
      'Scores are GSVA enrichment scores on log-CPM with Gaussian kernel. Positive logFC denotes higher activity in '
      'myelofibrosis. Full gene-level results and all 50 pathways are provided as machine-readable tables '
      '(DE_*.tsv, gsva_hallmark_MF_vs_NBM.tsv).']

# 写入补充表
sup = io.open(SUP, encoding='utf-8').read().rstrip() + '\n'
if 'Table S7.' not in sup:
    sup += '\n' + '\n'.join(S7) + '\n'
if 'Table S8.' not in sup:
    sup += '\n' + '\n'.join(S8) + '\n'
io.open(SUP, 'w', encoding='utf-8').write(sup)

# 修正正文编号与清单
ms = io.open(MS, encoding='utf-8').read()
ms = ms.replace('**Table S6.** Two-step Mendelian randomisation of 20 circulating inflammatory proteins',
                '**Table S7.** Two-step Mendelian randomisation of 20 circulating inflammatory proteins')
ms = ms.replace('**Table S7.** Differential expression (limma-voom) and GSVA of CD14+ monocytes',
                '**Table S8.** Differential expression (limma-voom) and GSVA of CD14+ monocytes')
ms = ms.replace('mediated proportions were \u22641.9% (Table S6)', 'mediated proportions were \u22641.9% (Table S7)')
# 清单里补上 S6（VEP）
if '**Table S6.**' not in ms:
    ms = ms.replace('**Table S7.** Two-step Mendelian',
                    '**Table S6.** Functional annotation of the 21 germline instruments (Ensembl VEP).\n'
                    '**Table S7.** Two-step Mendelian')
io.open(MS, 'w', encoding='utf-8').write(ms)
print('S7/S8 已生成；正文编号已对齐（S6=VEP, S7=蛋白MR, S8=单核细胞）')
