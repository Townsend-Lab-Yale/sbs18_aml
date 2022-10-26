# sbs18_aml
Analysis of contribution of SBS18 to acute myeloid leukemia, including in cases with RUNX1-RUNX1T1 fusions.

Summary of analysis scripts
---------------------------
- prep_mafs.R: Read in MAF data from data sources and produce analysis-ready MAF files.
- attribution_analysis.R: Perform the main analysis.
- population_scaled_effect_per_tumor.R: Code called in main analysis to attribute cancer effects to mutational processes, as described in [Cannataro et al.](https://doi.org/10.1093/molbev/msac084) This script is equivalent to
to the file with the same name in our [cancer_causes_and_effects repository](https://github.com/Townsend-Lab-Yale/cancer_causes_and_effects), commit 79d9aec.
- make_figure_1.R: Read in analysis outputs and create Figure 1.
- make_figure_2.R: Read in analysis outputs and create Figure 2.
- signature_stability_dotplot.R: Code called by make_figure_1.R to produce 1B dotplot.

Prerequisites for reproducibility
---------------------------------
- TCGA AML MAF data from [cBioPortal](https://www.cbioportal.org) under identifier _laml_tcga_pan_can_atlas_2018_.
- Pediatric AML MAF from [Gunnarsson et al.](https://www.nature.com/articles/s41375-021-01242-0), available on request [here](https://figshare.com/s/5a1ca3f39611c39bfaae).
- cancereffectsizeR package [(installation instructions)](https://townsend-lab-yale.github.io/cancereffectsizeR/), and all other R packages used in scripts (available from CRAN or Bioconductor).
