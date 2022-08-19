#' Check ondisc inputs
#'
#' Checks inputs (that come in ondisc format)
#'
#' @inheritParams run_sceptre_low_moi
#' @param low_moi does this dataset correspond to low MOI (TRUE) or high MOI (FALSE)?
#'
#' @return a thinned and verified multimodal ondisc matrix
check_ondisc_inputs <- function(mm_odm, gene_grna_group_pairs, form, gene_modality_name, grna_modality_name, grna_group_column_name, low_moi) {
  modality_names <- names(mm_odm@modalities)
  # 1. gene_modality_name must be in modality_names; set gene modality name to `gene`
  if (!(gene_modality_name %in% modality_names)) {
    stop(paste0("`", gene_modality_name, "` must be a modality name in the multimodal ondisc matrix."))
  } else {
    names(mm_odm@modalities)[names(mm_odm@modalities) == gene_modality_name] <- "gene"
  }

  # 2. "grna" must be in modality_names; set grna modality name to `grna`
  if (!(grna_modality_name %in% modality_names)) {
    stop(paste0("`", grna_modality_name, "` must be a modality name in the multimodal ondisc matrix."))
  } else {
    names(mm_odm@modalities)[names(mm_odm@modalities) == grna_modality_name] <- "grna"
  }

  # 3. "grna_group" must be a column of the grna feature covariate matrix
  grna_feature_covariates <- mm_odm |> ondisc::get_modality("grna") |> ondisc::get_feature_covariates()
  if (!(grna_group_column_name %in% colnames(grna_feature_covariates))) {
    stop(paste0("The feature covariates matrix of the `grna` ondisc matrix must have a column called `", grna_group_column_name, "`."))
  } else {
    grna_feature_covariates <- grna_feature_covariates |>
      dplyr::rename("grna_group" = !!grna_group_column_name)
  }
  # if low_moi is TRUE, check for the presence of "non-targeting" in the grna_group column
  if (low_moi) {
    if (!("non-targeting" %in% grna_feature_covariates $grna_group)) {
      stop("You are carrying out an analysis of low MOI data. The entry `non-targeting` must be present in the `", grna_group_column_name, "` column of the feature covariate matrix of the gRNA ondisc matrix.")
    }
  }
  mm_odm@modalities[["grna"]]@feature_covariates <- grna_feature_covariates

  # 4. "gene_id" must be a column of gene_grna_group_pairs
  if (!("gene_id" %in% colnames(gene_grna_group_pairs))) {
    stop("The gene_grna_group_pairs data frame must contain a column called `gene_id`.")
  }

  # 5. "grna_group" must be a column of gene_grna_group_pairs
  if (!("grna_group" %in% colnames(gene_grna_group_pairs))) {
    stop("The gene_grna_group_pairs data frame must contain a column called `grna_group`.")
  }

  # 6. check that the gene IDs in the gene_grna_group_pairs data frame are a subset of the gene IDs of the gene ODM
  odm_gene_ids <- mm_odm |> ondisc::get_modality("gene") |> ondisc::get_feature_ids()
  gene_grna_group_pairs_gene_ids <- as.character(gene_grna_group_pairs$gene_id)
  if (!all(gene_grna_group_pairs_gene_ids %in% odm_gene_ids)) {
    stop("The gene IDs in the gene_grna_group_pairs data frame are not a subset of the gene IDs in the gene ondisc matrix.")
  }

  # 7. check that the grna groups in the gene_grna_group_pairs data frame are a subset of the grna groups in the grna ODM
  odm_grna_groups <- grna_feature_covariates$grna_group
  gene_grna_group_pairs_grna_groups <- as.character(gene_grna_group_pairs$grna_group)
  if (!all(gene_grna_group_pairs_grna_groups %in% odm_grna_groups)) {
    stop("The grna groups in the gene_grna_group_pairs data frame are not a subset of the grna groups in the grna ondisc matrix.")
  }


  ###############
  # UPDATE MM ODM
  ###############
  # 1. apply a formula object to the global cell covariates
  if (!is.na(form) && form != "NA") {
    if (grepl("offset", form)) stop("Offsets are not currently supported in formulas within sceptre.")
    global_cell_covariates <- mm_odm |> ondisc::get_cell_covariates()
    global_cell_covariates_new <- stats::model.matrix(object = stats::as.formula(form),
                                               data = global_cell_covariates) |> as.data.frame()
    mm_odm@global_cell_covariates <- global_cell_covariates_new
  }

  # 2. Check for weird numbers in the global cell covariate matrix
  global_cell_covariates <- mm_odm |> ondisc::get_cell_covariates()
  if (ncol(global_cell_covariates) == 0) {
    stop("The global cell covariate matrix of the multimodal ondisc matrix should contain at least one column.")
  }
  for (col_name in colnames(global_cell_covariates)) {
    vect <- global_cell_covariates[[col_name]]
    if (any(vect == -Inf) || any(vect == Inf) || any(is.na(vect))) {
      stop(paste0("The column `", col_name, "` of the global cell covariate matrix contains entries that are either NA, -Inf, or Inf. Remove these entries (by, for example, removing the corresponding cells from the multimodal ondisc matrix)."))
    }
  }
  mm_odm <- ondisc:::thin_multimodal_odm(mm_odm, "grna", "gene", "grna_group")
  return(mm_odm)
}
