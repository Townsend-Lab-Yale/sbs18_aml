# Prepare MAF data from Gunnarsson et al. and TCGA. Add a sample-level fusion annotation and write out new MAF files that are ready for analysis.

# TCGA data were acquired from CBioPortal via https://cbioportal-datahub.s3.amazonaws.com/laml_tcga_pan_can_atlas_2018.tar.gz.
# The Gunnarsson MAF files were acquired directly from the authors.

# Run from project directory
# setwd(...)

library(data.table)
library(cancereffectsizeR) # 2.5.0
library(ces.refset.hg19) # 1.1.1

if (! dir.exists('prepped_mafs')) {
  dir.create('prepped_mafs')
} else {
  stop('This script appears to have already been run.')
}

## Gunnarsson data
gunnarsson_files = list.files('Gunnarsson_MAF/', full.names = T)
gunnarsson_maf = rbindlist(lapply(gunnarsson_files, fread))

# Fusion status from Gunnarsson Table S1
gunnarsson_maf[c("DX002", "DX016", "DX018", "DX022"), 
               fusion := 'rr', on = 'Tumor_Sample_Barcode']
gunnarsson_maf[c("DX007", "DX008", "DX009", "DX011", "DX012", "DX013", "DX019"), 
               fusion := 'negative', on = 'Tumor_Sample_Barcode']
gunnarsson_maf[c("DX001", "DX003", "DX004", "DX010", "DX014", "DX015", "DX021", "DX023", "DX024"),
               fusion := 'other', on = 'Tumor_Sample_Barcode']
stopifnot(all(! is.na(gunnarsson_maf$fusion))) # verify that's all the samples

## TCGA data ##
tcga_maf = fread('laml_tcga_pan_can_atlas_2018/data_mutations.txt')
tcga_clin = fread('laml_tcga_pan_can_atlas_2018/data_clinical_sample.txt', skip = 4)
tcga_clin = tcga_clin[unique(tcga_maf$Tumor_Sample_Barcode), on = 'SAMPLE_ID']
stopifnot(tcga_clin[, .(samples_per_patient = uniqueN(SAMPLE_ID)), by = 'PATIENT_ID'][, all(samples_per_patient == 1)])

tcga_fusion = fread('laml_tcga_pan_can_atlas_2018/data_fusions.txt')
tcga_rr_fusion_samples = tcga_fusion[Fusion == 'RUNX1-RUNX1T1', unique(Tumor_Sample_Barcode)]
tcga_other_fusion_samples = setdiff(tcga_fusion$Tumor_Sample_Barcode, tcga_rr_fusion_samples)

tcga_maf[tcga_rr_fusion_samples, fusion := 'rr', on = 'Tumor_Sample_Barcode']
tcga_maf[tcga_other_fusion_samples, fusion := 'other', on = 'Tumor_Sample_Barcode']
tcga_maf[is.na(fusion), fusion := 'negative']

## Data checks and filtering
# Combine MAFs into one table, run checks and apply filters, then print study-specific analysis-ready MAFs.
# For some reason,TCGA has some duplicate mutation records, which will get filtered out.
combined_maf = rbindlist(list(tcga = tcga_maf, gunnarsson = gunnarsson_maf),
                         idcol = 'cohort', fill = T)
checked_maf = preload_maf(maf = combined_maf, refset = ces.refset.hg19, keep_extra_columns = c('cohort', 'fusion'))
checked_maf = checked_maf[is.na(problem)]

# Some samples (13/742) have no remaining SNVs (only indels). Drop these since the attibutions method uses SNVs.
all_samples = unique(checked_maf$Unique_Patient_Identifier)
snv_samples = checked_maf[variant_type == 'snv', unique(Unique_Patient_Identifier)]
checked_maf = checked_maf[snv_samples, on = 'Unique_Patient_Identifier']

checked_maf[, output_name := paste0('prepped_mafs/', cohort, '_', 'prepped.maf')]
checked_maf = checked_maf[, .(Chromosome, Start_Position, Reference_Allele,
                              Tumor_Allele, Unique_Patient_Identifier, fusion, cohort, output_name)]

# For simplicity later, define sample groups that will be used for signature extraction (MutationalPatterns)
checked_maf[cohort == 'gunnarsson', cohort_mp_group := Unique_Patient_Identifier]
checked_maf[is.na(cohort_mp_group), cohort_mp_group := paste(cohort, fusion, sep = '.')]

## Save analysis-ready MAFs.
output_mafs = split(checked_maf, by = 'output_name', keep.by = F)
output_names = names(output_mafs)
for (i in 1:length(output_mafs)) {
  fwrite(output_mafs[[i]], output_names[i], sep = "\t")
}
  