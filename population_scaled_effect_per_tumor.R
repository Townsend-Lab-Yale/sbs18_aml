# From https://github.com/Townsend-Lab-Yale/cancer_causes_and_effects/blob/master/R/population_scaled_effect_per_tumor.R,
# as in commit 79d9aec. The population_scaled_effect_per_tumor function defined in this code will likely be incorporated 
# into a published package soon.

#'   Get signature contributions to variants
#'
#'   Builds a data table containing the contributions of each SNV signature to each
#'   variant in CES selection output for all tumors with each variant.
#'
#'   For each selected variant and each tumor with the variant, set the contribution of
#'   each SNV signature to be the weight of the signature in the tumor multipled by the
#'   relative rate of the trinuc context within the signature, divided by the relative
#'   rate of the trinuc SNV in the tumor. All required information (SNV signature
#'   definitions, tumor signature weights, trinuc-specific mutation rates, variant
#'   annotations) is taken from the input CESAnalysis.
population_scaled_effect_per_tumor <- function(ces_output, 
                                               run_name = NULL, 
                                               min_variant_freq = 2, 
                                               remove_nonsplice_silent = F,
                                               cores = 1) {
  if (! require("cancereffectsizeR")) {
    stop("Could not load cancereffectsizeR pacakge; is it installed?")
  }
  if (! require("data.table")) {
    stop("Could not load data.table package; is it installed?")
  }
  if(! is(ces_output, "CESAnalysis")) {
    stop("ces_output should be a CESAnalysis object", call. = F)
  }
  
  # Get signature definition matrix (rows = signature names, columns = relative rates of trinuc-context-specific SNVs)
  # Pre-2.3, just one set of signatures could be present
  running_old_version = FALSE
  if (is.data.frame(ces_output$reference_data$snv_signatures$signatures)) {
    running_old_version = TRUE
    signature_defs = ces_output$reference_data$snv_signatures$signatures
  } else {
    num_signature_sets = length(ces_output$reference_data$snv_signatures)
    if(num_signature_sets == 0) {
      stop("Input CESAnalysis contains no associated signature definitions.")
    } 
    if(num_signature_sets > 1) {
      msg = paste0("Unusual situation: Input CESAnalysis has more than one set of associated SNV mutational signature definitions. ",
                   "It should be possible to create equivalent analysis using just one signature set.")
      stop(pretty_message(msg, emit = F))
    }
    signature_defs = ces_output$reference_data$snv_signatures[[1]]$signatures
  }
  signature_names = rownames(signature_defs)
  
  # Get data.table with signature weights (one row per tumor)
  signature_weights = get_signature_weights(ces_output)
  if (is.null(signature_weights)) {
    stop("Can't run because no signature weights are associated with this analysis.\n",
    "If you have them, you can add them to the CESAnalysis with set_signature_weights().")
  }
  if(nrow(signature_weights) != nrow(ces_output$samples)) {
    stop("Not all samples in CESAnalysis have assigned signature weights, so can't run.")
  }
  
  # For efficiency, drop signatures where weights are always 0
  not_zero = signature_weights[, (sapply(.SD, function(x) any(x != 0))), .SDcols = signature_names]
  not_zero = names(which(not_zero == T))
  signature_weights = signature_weights[, c("Unique_Patient_Identifier", ..not_zero)]
  setkey(signature_weights, "Unique_Patient_Identifier")
  
  signature_names = signature_names[signature_names %in% not_zero]
  signature_defs = signature_defs[not_zero, ]
  
  flux_colnames = paste(signature_names, "flux", sep = "_")
  nrsi_colnames = paste(signature_names, "nrsi", sep = "_")
  
  
  # Get selection intensity data
  all_si_output = ces_output$selection
  if (length(all_si_output) == 0) {
    stop("No variant effect data associated with the CESAnalysis, so can't run")
  } else {
    if (is.null(run_name)) {
      if(length(all_si_output) == 1) {
        selection_data = all_si_output[[1]]
      } else {
        stop("ces_variant() has been run on this CESAnalysis multiple times. Specify which run to use with run_name.")
      }
    } else {
      if (! is.character(run_name) || length(run_name) != 1)  {
        stop("run_name should be 1-length character")
      }
      if (! run_name %in% names(all_si_output)) {
        names_used = names(all_si_output)
        stop("Chosen run_name not found. Found ces_variant runs named ", paste(names_used, collapse = ", "), '.')
      } else {
        selection_data= all_si_output[[run_name]]
      }
    }
  }
  
  # make sure there is just one fitted selection intensity, and give it the standard name
  si_col = attr(selection_data, "si_cols", exact = T)
  if (length(si_col) > 1) {
    stop("There are multiple SI columns in variant effect data, which isn't yet supported.")
  }
  setnames(selection_data, si_col, "selection_intensity")
  
  # Handle pre-2.3 prevalence column name
  setnames(selection_data, 'total_maf_freq', 'maf_prevalence', skip_absent = T)
  setnames(selection_data, 'maf_frequency', 'maf_prevalence', skip_absent = T)
  setnames(selection_data, 'maf_prevalence', 'included_with_variant', skip_absent = T)                                       
  selection_data = selection_data[included_with_variant >= min_variant_freq]
  
  message("Using variant effect data from ", selection_data[, .N], " variants...")
  
  # if requested, remove silent variants unless they are near splice sites
  if (remove_nonsplice_silent) {
    nonsplice_silent_ind = selection_data[variant_type == "aac" & essential_splice == FALSE & aa_ref == aa_alt, which = T]
    num_nonsplice_silent = length(nonsplice_silent_ind)
    if (num_nonsplice_silent > 0) {
      selection_data = selection_data[! nonsplice_silent_ind]
      message("Leaving out ", num_nonsplice_silent, " synonymous coding variants that are not near splice sites...")
    }
  }
  
  trinuc_rates = ces_output$trinuc_rates
  setkey(trinuc_rates, "Unique_Patient_Identifier")
  
  # Ensure that we're only handling AACs and SNVs
  other_variant_type_index = selection_data[! variant_type %in% c("snv", "aac"), which = T]
  if (length(other_variant_type_index) > 0) {
    selection_data = selection_data[! other_variant_type_index, ]
    warning("Some variants are neither coding mutations nor noncoding SNVs; these were skipped.")
  }
  
  # Make sure filtering didn't remove all variants
  if(selection_data[, .N] == 0) {
    stop("No variants remain.")
  }
  
  # Collect SNV IDs: just the variant ID for SNVs; for AACs, pull information from mutation annotations
  selection_data[variant_type == "snv", snv_ids := list(list(variant_id)), by = "variant_id"]
  selection_data[variant_type == "aac", snv_ids := ces_output@mutations$amino_acid_change[variant_id, constituent_snvs, on = 'aac_id']]
  missing_aac_index = selection_data[, bad := is.null(snv_ids[[1]]), by = "variant_id"][bad == T, which = T]
  if (length(missing_aac_index) > 0) {
    missing = selection_data[missing_aac_index, variant]
    selection_data = selection_data[! missing_aac_index]
    for (variant in missing) {
      warning(sprintf("Variant %s has incomplete annotation, so it was skipped.", variant))
    }
  }
  
  variant_snv_pairings = selection_data[, .(snv = unlist(snv_ids), selection_intensity), by = "variant_id"]
  
  # get trinucleotide contexts of all SNVs
  if (running_old_version) {
    reselected_snv = select_variants(ces_output, variant_passlist = variant_snv_pairings$snv)
  } else {
    reselected_snv = select_variants(ces_output, variant_ids = variant_snv_pairings$snv)
  }
  reselected_snv = reselected_snv[variant_snv_pairings$snv, on = "variant_id", nomatch = NULL]
  
  # this shouldn't happen
  if (reselected_snv[, .N] != variant_snv_pairings[, .N]) {
    stop("Failed to select all SNVs from CESAnalysis (please report the bug).")
  }
  
  variant_snv_pairings[, trinuc_context := reselected_snv$trinuc_mut]
  
  # Ensure that all SNVs have trinuc context annotation
  bad_snv_index = variant_snv_pairings[is.na(trinuc_context), which = T]
  if (length(bad_snv_index) > 0) {
    missing = variant_snv_pairings[bad_snv_index, variant]
    for (variant in missing) {
      warning(sprintf("Variant %s has incomplete SNV annotations, so it was skipped.", variant))
    }
    variant_snv_pairings = variant_snv_pairings[! variant %in% missing]
  }
  
  # Subset MAF to just the SNVs we need for quicker access
  maf = ces_output$maf[variant_snv_pairings$snv, on = "variant_id", nomatch = NULL]
  setkey(maf, "variant_id")
  
  # For each SNV that constitutes a coding mutation (or just the one mutation for nonconding SNVs), 
  # get the contribution of each signature to each tumor that has the mutation.
  # Note that for AACs, cancereffectsizeR's exclusion of di- and tri-nucleotide substitutions (as well as indels)
  # makes it so that each tumor with the AAC has exactly one of the constituent SNVs.
  process_snv = function(variant_id, snv, trinuc_context, selection_intensity) {
    
    tumors = maf[snv, Unique_Patient_Identifier, nomatch = NULL]
    if (length(tumors) == 0) {
      return(NULL)
    }
    # relative rates of the given trinuc SNV for each signature
    sig_trinuc_rate = as.numeric(signature_defs[, trinuc_context])
    snv_output = data.table(variant = variant_id, Unique_Patient_Identifier = tumors)
    
    tumor_trinuc_rates = trinuc_rates[tumors, ..trinuc_context][[1]] # get vector, not 1-column DT
    weights = signature_weights[tumors, ..signature_names]
    
    # For each row (tumor) in weights table, do (weights * sig_trinuc_rate) / tumor_trinuc_rate
    # Opaque syntax for efficiency
    contrib = setDT(Map(`*`, weights, sig_trinuc_rate))
    if(contrib[, .N] == 1) {
      contrib = as.data.table(t(apply(contrib, 2, function(x) x/tumor_trinuc_rates)))
    } else {
      contrib = as.data.table(apply(contrib, 2, function(x) x / tumor_trinuc_rates))
    }
    
    # Multiply contributions by the scalar selection intensity
    # Absurdly, the obscure syntax runs 5x faster than just using "nrsi = contrib * selection_intensity"
    nrsi = setDT(Map(`*`, contrib, selection_intensity))
    snv_output[, (flux_colnames) := contrib]
    snv_output[, (nrsi_colnames) := nrsi]
    return(snv_output)
  }
  
  # Get results from the process_snv function, over cores
  par_results <- parallel::mcmapply(FUN = process_snv, variant_snv_pairings$variant, variant_snv_pairings$snv, 
                                    variant_snv_pairings$trinuc_context, variant_snv_pairings$selection_intensity, SIMPLIFY = F,
                                    mc.cores = cores)
  
  return(rbindlist(par_results))
  
}
