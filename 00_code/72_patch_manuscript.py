# -*- coding: utf-8 -*-
"""72_patch_manuscript.py —— 把三项新增分析并入 Manuscript_Leukemia.md（参考文献/表头/补充材料）"""
import io

P = 'C:/Users/Hero/WorkBuddy/2026-10-14-23-35-01/outputs/pv_hotspot/T1-1_CHIP_comorbidity/10_manuscript/Manuscript_Leukemia.md'
s = io.open(P, encoding='utf-8').read()

old31 = ("31. Terao C, Suzuki A, Momozawa Y, et al. Chromosomal alterations among age-related "
         "haematopoietic clones in Japan. *Nature.* 2020;584(7820):130-135.")
new31 = old31 + (
    "\n32. Kong T, Laranjeira ABA, Letson CT, et al. RSK1 is an exploitable dependency in "
    "myeloproliferative neoplasms and secondary acute myeloid leukemia. *Nat Commun.* 2025;16:492."
    "\n33. Zhao JH, Stacey D, Eriksson N, et al. Genetics of circulating inflammatory proteins "
    "identifies drivers of immune-mediated disease risk and therapeutic targets. "
    "*Nat Immunol.* 2023;24:1540-1551.")
if '32. Kong T' not in s:
    assert old31 in s, 'ref31 not found'
    s = s.replace(old31, new31)

s = s.replace(
    '**Word counts:** Abstract 189 words; main text 3,422 words; 3 figures; 3 tables; 31 references.',
    '**Word counts:** Abstract 189 words; main text 4,009 words; 3 figures; 4 tables; 33 references.')

oldS = ('**Table S5.** Reverse-direction MR results.\n'
        '**Figure S1.** Leave-one-out forests and funnel plots.')
EOF_MARK = None
if 'Table S6' not in s:
    assert oldS in s, 'supp list not found'
    newS = (
        '**Table S5.** Reverse-direction MR results.\n'
        '**Table S6.** Two-step Mendelian randomisation of 20 circulating inflammatory proteins: '
        'step 1 (CHIP liability to protein), step 2 (protein to MPN and polycythaemia vera), '
        'and mediated proportions.\n'
        '**Table S7.** Differential expression (limma-voom) and GSVA of CD14+ monocytes across the '
        'MPN spectrum (GSE283711), including the CHIP-direction expression score by group.\n'
        '**Figure S1.** Leave-one-out forests and funnel plots.\n'
        '\n### Supplementary Methods\n\n'
        '**Blood cis-eQTL colocalisation.** Full cis-eQTL summary statistics for whole blood '
        '(eQTLGen, 31,684 samples; 2019-12-11 release) were streamed and restricted to \u00b11 Mb '
        'around each of the 21 germline instruments, retaining 4,434,369 SNP\u2013gene pairs. '
        'eQTLGen distributes Z scores and no effect sizes; beta and SE were therefore reconstructed '
        'as beta = Z/sqrt(2f(1-f)(N + Z^2)) and SE = 1/sqrt(2f(1-f)(N + Z^2)), with f the frequency '
        'of the assessed allele taken from the eQTLGen SNP frequency file and N the per-variant '
        'sample size. CHIP effect estimates were re-oriented to the eQTL assessed allele, including '
        'strand-complement handling. Transcripts with at least one cis-eQTL at FDR < 0.05 and >=50 '
        'shared variants were tested with `coloc.abf` (quantitative traits, sdY = 1, MAF supplied, '
        'N = 418,401 for CHIP, p12 = 1 x 10^-5). Windows were analysed independently and 511 '
        'SNP\u2013gene pairs were tested in total.\n\n'
        '**CD14+ monocyte transcriptome.** Raw counts for 68 CD14+ monocyte libraries (GSE283711) '
        'were filtered with `filterByExpr` after TMM normalisation, and differential expression was '
        'estimated with limma-voom (contrasts: myelofibrosis vs normal, polycythaemia vera vs normal, '
        'essential thrombocythaemia vs normal, myelofibrosis vs polycythaemia vera). GSVA was run '
        'against the 50 MSigDB hallmark gene sets on log-CPM values with Gaussian kernel, and '
        'pathway scores were compared by limma. The CHIP-direction score was the mean, across the 52 '
        'candidate genes represented, of the per-sample standardised log-CPM multiplied by the sign '
        'of the CHIP-to-expression effect derived from single-instrument Wald ratios, tested by '
        'Kruskal\u2013Wallis and by Spearman correlation against ordinal disease stage.\n\n'
        '**Two-step protein MR.** pQTL summary statistics for 20 Olink inflammation-panel proteins '
        '(GCST90274758\u2013GCST90274848) were downloaded in full and scanned for the 21 germline '
        'instruments (step 1) and for genome-wide significant variants (P < 5 x 10^-8, step 2). '
        'Step 1 estimates used inverse-variance weighting of per-variant Wald ratios with '
        'multiplicative random-effects correction and Cochran Q. Step 2 instruments were restricted '
        'to \u00b11 Mb of the encoding gene and pruned by physical distance (250 kb, retaining the '
        'most significant variant per bin), giving 22 variants across eight proteins; outcome '
        'estimates were taken from FinnGen release 12 (MPN, polycythaemia vera) and aligned to the '
        'pQTL effect allele. Mediated proportion = (a x b)/c with c the total CHIP-to-outcome IVW '
        'estimate, and delta-method standard errors.')
    s = s.replace(oldS, newS)

io.open(P, 'w', encoding='utf-8').write(s)
print('patched: refs 32/33, word counts, Tables S6-S7, Supplementary Methods')
