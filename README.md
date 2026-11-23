# Genomic Variant Association Analysis

Statistical analysis of genomic variants in *Arabidopsis thaliana*, quantifying association between variants and testing its significance with several complementary methods.

## Overview

This project studies association in Arabidopsis genomic data in two parts:

**Part 1 — Association between two genomic variants.** A binary presence/absence matrix for two genomic variants across 199 individual plants (`data/p1a.csv`). Association is quantified with three measures:

| Measure | Value | Interpretation |
|---|---|---|
| Mutual information | 0.033 | Minimal shared information |
| Jaccard index | 0.0 | No overlap — the variants never co-occur in this sample |
| Chi-squared χ² | 7.878 (p ≈ 0.096) | Not significant at the 0.05 level |

All three methods agree: there is **no significant association** between the two variants, and they appear mutually exclusive in this sample.

**Part 2 — Correlation between gene expression profiles.** Pearson correlation between paired gene measurements in three datasets:

| Dataset | n | Pearson r | p-value | Interpretation |
|---|---|---|---|---|
| `p2a.csv` | 2,400 | 0.381 | 1.04 × 10⁻⁸³ | Moderate correlation, highly significant |
| `p2b.csv` | 110 | 0.931 | 3.74 × 10⁻⁴⁹ | Strong correlation, highly significant |
| `p2c.csv` | 2,100 | 0.041 | 0.059 | Weak correlation, not significant |

Ranking by strength of association is consistent across correlation, p-value, and scatter-plot inspection: **p2b > p2a > p2c**. The p2c scatter plot shows a curvilinear pattern, which the linear Pearson coefficient does not fully capture — a reminder that correlation alone can miss non-linear structure.

## Data

- `data/p1a.csv` — Binary matrix (199 rows × 2 columns): presence (1) or absence (0) of two genomic variants per plant. No header row.
- `data/p2a.csv` — 2,400 paired gene-expression measurements.
- `data/p2b.csv` — 110 paired gene-expression measurements.
- `data/p2c.csv` — 2,100 paired gene-expression measurements.

## Methods

- **Mutual information** — information-theoretic dependence between the two binary variants.
- **Jaccard index** — overlap of variant presence across individuals.
- **Chi-squared test** — independence test on the 2 × 2 contingency table.
- **Pearson correlation** — linear correlation between paired continuous measurements, with significance testing.
- **Scatter plots** — visual comparison of the three Part 2 datasets (see the report).

## Files

```
├── README.md
├── data/
│   ├── p1a.csv
│   ├── p2a.csv
│   ├── p2b.csv
│   └── p2c.csv
└── report/
    ├── Genomic_Variant_Report.pdf
    └── Genomic_Variant_Report.docx
```

The full write-up, including the analysis, figures, and discussion, is in `report/Genomic_Variant_Report.pdf` (also available as `.docx`). The analysis code is documented inside the report.
