# Run from project analysis directory
# setwd(...)

library(data.table)
library(cancereffectsizeR) # 2.5.0
library(ces.refset.hg19) # 1.1.1
library(MutationalPatterns) # 3.4.1

# Load the function that does effect share calculations (from cancer_causes_and_effects repository; see file)
source("population_scaled_effect_per_tumor.R")

sample_data_cols = c("fusion", "cohort", "cohort_mp_group") # sample-level data to retain from MAF files
cesa = CESAnalysis('ces.refset.hg19')
cesa = load_maf(cesa, maf = 'prepped_mafs/tcga_prepped.maf', sample_data_cols = sample_data_cols)
cesa = load_maf(cesa, maf = "prepped_mafs/gunnarsson_prepped.maf", coverage = 'genome', sample_data_cols = sample_data_cols)

# Estimate neutral gene mutation rates
cesa = gene_mutation_rates(cesa, covariates = 'default')

## Signature analysis
# For all analysis, we'll restrict to signatures that have been established as appearing in AML.
# For Gunnarsson (WGS), run each sample normally through MutationalPatterns.
cosmic_signature_defs = ces.refset.hg19$signatures$COSMIC_v3.2$signatures
signature_exclusions = suggest_cosmic_signature_exclusions(cancer_type = 'LAML', quiet = TRUE)
included_signatures = setdiff(rownames(cosmic_signature_defs), signature_exclusions)
signatures_for_mp = t(cosmic_signature_defs[included_signatures, ]) # MP wants transpose

# Get cohort SNV counts for TCGA
tcga_cohort_counts = sapply(c('tcga.rr', 'tcga.other', 'tcga.negative'), 
  function(x) {
    curr_samples = cesa$samples[x, Unique_Patient_Identifier, on = "cohort_mp_group"]
    maf_subset = cesa$maf[curr_samples, on = "Unique_Patient_Identifier"]
    sample_counts = trinuc_snv_counts(maf = maf_subset, genome = cesa$reference_data$genome)
    return(rowSums(sample_counts))
  })

gunnarsson_samples = cesa$samples[cohort == 'gunnarsson', Unique_Patient_Identifier]
gunnarsson_counts = trinuc_snv_counts(maf = cesa$maf[gunnarsson_samples, on = "Unique_Patient_Identifier"],
                                      genome = cesa$reference_data$genome)

counts_for_mp = cbind(gunnarsson_counts, tcga_cohort_counts)

set.seed(41822) # for reproducibility (chosen from date)
mp_output = fit_to_signatures_bootstrapped(mut_matrix = counts_for_mp, signatures = signatures_for_mp,
                                           n_boots = 100)
# save full bootstrapped output
write.table(mp_output, 'results/raw_mp_out.txt', sep = "\t")

# First, make a table of (sample-averaged) Gunnarsson attributions for use in boxplot.
mp_attribution_table = as.data.table(mp_output, keep.rownames = 'sample_boot')
mp_attribution_table[, Unique_Patient_Identifier := gsub('_\\d+$', '', sample_boot)]
gunnarsson_attribution_table = mp_attribution_table[gunnarsson_samples, on = 'Unique_Patient_Identifier']
signature_names = intersect(colnames(signatures_for_mp), colnames(mp_output))
gn_attributions_sample_averaged = gunnarsson_attribution_table[, lapply(.SD, mean), by = "Unique_Patient_Identifier", .SDcols = signature_names]
gn_attributions_sample_averaged[cesa$samples, fusion := fusion, on = "Unique_Patient_Identifier"]
setcolorder(gn_attributions_sample_averaged, c("Unique_Patient_Identifier", "fusion"))
fwrite(gn_attributions_sample_averaged, 'results/gunnarsson_MP_attribution_counts.txt', sep = "\t")
# see boxplot script for significance testing on gunnarsson_MP_attribution_counts.txt

# From raw MP output, produce a matrix of sample-averaged signature weights to load into CESAnalysis.
# Mutation counts are normalized into relative weights that sum to 1. For TCGA cohort, we turn the cohort/fusion-group
# weights into sample weights (by repeating values appropriately).
weight_mat = t(apply(mp_output, 1, function(x) x = x/sum(x)))
weight_table = as.data.table(weight_mat, keep.rownames = "sample_boot")
weight_table[, cohort_mp_group := gsub('_\\d+$', '', sample_boot)]
weight_table_averaged = weight_table[, lapply(.SD, mean), by = "cohort_mp_group", .SDcols = signature_names]

weight_table_by_sample = weight_table_averaged[cesa$samples[, .(Unique_Patient_Identifier, cohort_mp_group)], on = 'cohort_mp_group']
weight_table_by_sample[, cohort_mp_group := NULL]
weight_table_by_sample[, (signature_exclusions) := 0] # put in 0 values for excluded signatures
cesa = set_signature_weights(cesa, "COSMIC_v3.2", weight_table_by_sample)

# Estimate cancer effect size using default model of selection.
cesa = ces_variant(cesa, run_name = 'cohort', variants = cesa$variants)

# Get NRSI table
pse_cohort = population_scaled_effect_per_tumor(ces_output = cesa, run_name = 'cohort', min_variant_freq = 1)

# Save CESAnalysis and NRSI table
fwrite(pse_cohort, 'results/population_scaled_effects.txt', sep = "\t")
save_cesa(cesa, 'results/cesa.rds')

# Produce a table of effects by sample
#pse_cohort_gn = pse_cohort[gunnarsson_samples, on = "Unique_Patient_Identifier"]
nrsi_cols = names(pse_cohort[, .SD, .SDcols = patterns('SBS.*nrsi')])
effects_by_sig = data.table()
for (nrsi_col in nrsi_cols) {
  sig_name = gsub('_nrsi', '', nrsi_col)
  tmp = copy(pse_cohort)
  setnames(tmp, nrsi_col, "tmp_nrsi")
  effects = tmp[, .(signature = sig_name, effect_share = sum(tmp_nrsi)/sum(.SD)), .SDcols = patterns("_nrsi"), by = "Unique_Patient_Identifier"]
  effects_by_sig = rbind(effects_by_sig, effects)
}

melted_sig_weights = melt(cesa$mutational_signatures$biological_weights, id.vars = 'Unique_Patient_Identifier', 
                          measure.vars = signature_names, variable.name = 'signature')
effects_by_sig[melted_sig_weights, sig_weight := value, on = c("Unique_Patient_Identifier", "signature")]
effects_by_sig[cesa$samples, c("fusion", "cohort") := list(fusion, cohort), on = "Unique_Patient_Identifier"]

fwrite(effects_by_sig, 'results/effect_share_by_signature.txt', sep = "\t")

# Produce fusion group averages for TCGA.
# Since different samples have different numbers of mutations, we go back to the full pse output.
tcga_samples = cesa$samples[cohort == 'tcga', Unique_Patient_Identifier]
pse_tcga = pse_cohort[tcga_samples, on = "Unique_Patient_Identifier"]
pse_tcga[cesa$samples, fusion := fusion, on = "Unique_Patient_Identifier"]
tcga_fusion_effects_by_sig = data.table()
for (nrsi_col in nrsi_cols) {
  sig_name = gsub('_nrsi', '', nrsi_col)
  tmp = copy(pse_tcga)
  setnames(tmp, nrsi_col, "tmp_nrsi")
  effects = tmp[, .(signature = sig_name, effect_share = sum(tmp_nrsi)/sum(.SD)), .SDcols = patterns("_nrsi"), by = "fusion"]
  tcga_fusion_effects_by_sig= rbind(tcga_fusion_effects_by_sig, effects)
}
melted_sig_weights[cesa$samples, fusion := fusion, on = 'Unique_Patient_Identifier']
for_tcga = melted_sig_weights[tcga_samples, .SD[1], by = c('fusion', 'signature'), on = "Unique_Patient_Identifier"][, -"Unique_Patient_Identifier"]
tcga_fusion_effects_by_sig[for_tcga, sig_weight := value, on = c('fusion', 'signature')]
fwrite(tcga_fusion_effects_by_sig, 'results/tcga_fusion_group_effect_shares.txt', sep = "\t")

# mutation weights / effect share significance testing
gn_effects_by_sig = effects_by_sig[cohort == 'gunnarsson']

gn_effects_by_sig[signature == 'SBS1', wilcox.test(sig_weight, effect_share, paired = T)]
gn_effects_by_sig[signature == 'SBS5', wilcox.test(sig_weight, effect_share, paired = T)]
gn_effects_by_sig[signature == 'SBS18', wilcox.test(sig_weight, effect_share, paired = T)]



