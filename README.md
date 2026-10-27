# ARIEM: Biobank-Scale Estimation of Recessive Metabolic Disease Burden in South Asia

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![R 4.3+](https://img.shields.io/badge/R-4.3+-blue.svg)](https://www.r-project.org/)

This repository contains the complete analytical codebase, variant curation tables, spatial population genetics pipelines, and display figure generation scripts supporting the study:

> **Biobank-Scale Estimation of Recessive Metabolic Disease Burden in South Asia: Implications for Genetic Screening**  
> *Atanu Kumar Dutta¹\*, Soumyajit Mondal¹*  
> ¹ Department of Biochemistry, All India Institute of Medical Sciences (AIIMS), Kalyani, West Bengal, India.  
> \* Correspondence: `atanu.biochem@aiimskalyani.edu.in`  
> *Submitted to: Human Genetics (Springer Nature) - Special Issue: "Genetic epidemiology in the era of biobanks and multi-omics: study designs, methods, and translation"*

---

## 📌 Study Overview

South Asian populations encompass extensive ethnolinguistic diversity and pervasive founder events resulting from millenia of endogamous social structure. Traditional newborn screening panels designed for European or East Asian ancestries fail to capture the true burden of autosomal-recessive inborn errors of metabolism (IEMs) across the subcontinent.

In this study, we evaluate **696 autosomal-recessive IEM genes** across:
1. **gnomAD v4 South Asian genomes** ($n = 14,606$ alleles / 7,303 individuals)
2. **GenomeIndia whole-genome sequencing dataset** ($n = 9,768$ individuals across 83 distinct endogamous groups)

### Key Findings
- **Random-Mating Burden**: Collective broad-set disease prevalence under panmixia is **255.1 per 100,000 live births** (1 in 392).
- **Consanguinity & Endogamy Expansion**: Under an average inbreeding coefficient of $F = 0.02$, cumulative burden expands **8.1-fold to 2,072.1 per 100,000** (1 in 48), with **87.7% of disease incidence attributable to endogamous structure**.
- **12 Keystone Genes**: Just 12 genes (*MAN2B2, DPYD, PYGL, ATP7B, PAH, MMUT, PYGM, DHCR24, SPG7, SACS, GLDC, MLYCD*) account for $249.7$ per 100,000 cases under $F = 0.02$ (12.1% of all cases).
- **Purifying Selection Signatures**: ClinVar Pathogenic/Likely Pathogenic (P/LP) alleles exhibit significant depletion of extreme homozygote excess compared to predicted missense variants ($3.29\\%$ vs $8.15\\%$, Fisher's exact $P = 1.6 \\times 10^{-5}$, OR = 0.38), confirming strong purifying selection against morbid homozygosity.
- **Isolation-by-Distance (IBD)**: Genomic differentiation across 83 groups ($F_{\\text{ST}} = 0.0086$, mean geographic distance 1,114.0 km) demonstrates weak spatial isolation (Mantel $r = 0.0589, P = 0.096$), indicating that social endogamy and deep founder events predominate over geographic separation.
- **Pilot Screening Efficiency**: Detecting 5 cases at $80\\%$ statistical power requires $\\sim 15,700$ unscreened births under broad curation vs $\\sim 26,000$ under strict ClinVar curation.

---

## 📂 Repository Structure

```text
ARIEM_GitHub_Repository/
├── README.md                           # Master repository documentation (this file)
├── LICENSE                             # MIT open-source license
├── requirements.txt                    # Python package requirements
├── environment.yml                     # Conda virtual environment specification
├── .gitignore                          # Standard git ignore definitions
│
├── docs/                               # Manuscript and provenance documentation
│   ├── manuscript_human_genetics.md    # Main text (Numbered reference format)
│   ├── manuscript_human_genetics_author_year.md # Main text (Author-Year format)
│   ├── supplement_human_genetics.md    # Supplementary material text and tables
│   ├── source_provenance.md            # Raw data accessions and IEC approvals
│   └── author_confirmation_checklist.md# Springer Nature author checklist
│
├── scripts/                            # Reproducible analytical pipelines
│   ├── actionability_tiers.R           # ACMG/ClinGen tier categorization
│   ├── ariem_scikit_pipeline.py        # Empirical homozygosity and machine learning
│   ├── homozygote_selection.R          # Non-equilibrium selection analysis
│   ├── pairwise_decomposition.R        # Gene-level population stratification
│   ├── panelapp_696_clinical_tiers.R   # PanelApp curated IEM gene mapping
│   ├── plots_tables_S12.R              # Main Figure 4 and Supplementary Fig. S12
│   ├── rebuild_figure_4.R              # 600 DPI publication figure builder
│   ├── run_three_popgen_enhancements.R # Genovis, PopGen, and spatial IBD enhancements
│   ├── run_revision_S12_validation.R   # Statistical validation and bootstrap CIs
│   └── verify_all_final.py             # Integrity check script
│
├── data/                               # Curated tables and variant catalogs
│   ├── tables/                         # Formatted Excel tables (Table 1, Table 2, S1-S14)
│   │   ├── Table_1.xlsx                # Top 12 metabolic genes under F=0.02
│   │   ├── Table_2.xlsx                # GenomeIndia cohort burden comparison
│   │   ├── Table_S1.xlsx ... S14.xlsx  # Supplementary tables S1 through S14
│   │   └── Table_IBD_Mantel_Summary.xlsx # Spatial Mantel matrix summary
│   ├── csv/                            # CSV/TSV data tables
│   │   ├── MAF_83populations.tsv       # Minor allele frequencies across 83 groups
│   │   ├── panelapp_696_clinical_actionability_tiers.csv # Clinical tier annotations
│   │   └── popgen_advanced_summary.csv # PopGen metrics
│   └── source_data/                    # Raw evidence matrices and json summaries
│
└── figures/                            # Publication-grade figures (PNG, PDF, TIFF)
    ├── Figure_1.png / .pdf             # Population architecture and study overview
    ├── Figure_2.png / .pdf             # Gene-level disease incidence and ranking
    ├── Figure_3.png / .pdf             # Consanguinity scaling and founder expansions
    ├── Figure_4.png / .pdf / .tiff     # Non-equilibrium selection and screening power
    └── Supplementary_Figure_S1 ... S15 # Complete supplementary figures S1-S15
```

---

## 🚀 Quickstart & Reproduction

### 1. Environment Setup

#### Using Conda:
```bash
git clone https://github.com/aiims-kalyani-biochem/ARIEM-South-Asia-Genomics.git
cd ARIEM-South-Asia-Genomics
conda env create -f environment.yml
conda activate ariem-popgen
```

#### Using Pip & Native R:
```bash
pip install -r requirements.txt
# Ensure R (>= 4.3.0) is installed with:
# data.table, ggplot2, ggprism, patchwork, scales, vegan
```

### 2. Reproduce Core Visualizations

To re-render **Figure 4** (600 DPI publication graphic matching empirical $\\hat{F}$ bootstrap CIs):
```bash
Rscript scripts/rebuild_figure_4.R
```

To run the complete **GENOVIS and spatial isolation-by-distance** pipeline:
```bash
Rscript scripts/run_three_popgen_enhancements.R
```

To verify cross-references, line numbers, and citations across all manuscript files:
```bash
python scripts/verify_all_final.py
```

---

## 📜 Ethical Approvals & Provenance

The GenomeIndia source study reported institutional ethics committee (IEC) approvals and written informed consent for all participants across the 83 endogamous ethnolinguistic groups (Research Square preprint doi: `10.21203/rs.3.rs-11445997/v1`). gnomAD v4 South Asian genome reference data were acquired from open-access public release buckets. Full provenance and audit details are documented in [`docs/source_provenance.md`](docs/source_provenance.md).

---

## ✍️ Citation

If you use this code, data, or methodology in your research, please cite:

```bibtex
@article{dutta2026ariem,
  title={Biobank-Scale Estimation of Recessive Metabolic Disease Burden in South Asia: Implications for Genetic Screening},
  author={Dutta, Atanu Kumar and Mondal, Soumyajit},
  journal={Human Genetics},
  year={2026},
  note={Special Issue: Genetic epidemiology in the era of biobanks and multi-omics}
}
```

---

## 📄 License

This repository is distributed under the terms of the [MIT License](LICENSE).
