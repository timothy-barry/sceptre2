#' Run sceptre (low MOI)
#'
#' Runs sceptre on a low MOI dataset.
#'
#' @param mm_odm a multimodal ondisc matrix containing gene and grna (either expression or assignment) modalities
#' @param gene_grna_group_pairs a data frame with columns `gene_id` and `grna_group` giving the gene ID / grna group pairs to analyze.
#' @param form a string that describes how to adjust for the covariates
#' @param threshold threshold to use for assigning gRNAs to cells
#' @param B number of resamples to draw for permutation/conditional randomization test
#' @param gene_modality_name name of the gene modality (e.g., "gene_expression")
#' @param grna_modality_name name of the grna modality (e.g., "grna_expression")
#' @param grna_group_column_name name of the column within the feature covariate matrix of the gRNA ODM containing the gRNA groups
#' @param n_pairs_to_sample number of pairs to sample from `response_grna_group_pairs` for debugging purposes
#' @param full_output return the full output (TRUE) or a simplified, reduced output (FALSE)?
#'
#' @return a data frame with columns `response_id`, `grna_group`, and `p_value` (as well as other entries).
#' @export
#'
#' @examples
#' \dontrun{
#' library(ondisc)
#' # a small example
#' tap_seq_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"), "data/schraivogel/ground_truth_tapseq/")
#' mm_fp <- paste0(tap_seq_dir, "multimodal_metadata.rds")
#' odm_fps <- paste0(tap_seq_dir, c("gene/matrix.odm", "grna_expression/matrix.odm"))
#' mm_odm <- ondisc::read_multimodal_odm(odm_fps = odm_fps, multimodal_metadata_fp = mm_fp)
#' gene_grna_group_pairs <- expand.grid(gene_id = mm_odm |> get_modality("gene") |> get_feature_ids(),
#' grna_group = mm_odm |> get_modality("grna_expression") |> get_feature_covariates() |> dplyr::pull(target) |> unique())
#' gene_grna_group_pairs <- gene_grna_group_pairs |> dplyr::filter(grna_group != "non-targeting")
#' form <- "~ log(gene_n_nonzero) + log(gene_n_umis) + batch"
#' B <- 2500; gene_modality_name <- "gene"; grna_modality_name <- "grna_assignment";
#' grna_group_column_name <- "target"; n_pairs_to_sample <- 50; full_output <- FALSE; side <- "both"
#' result <- run_sceptre_low_moi(mm_odm, gene_grna_group_pairs, form, B, gene_modality_name,
#' grna_modality_name, grna_group_column_name, n_pairs_to_sample)
#' }
run_sceptre_low_moi <- function(mm_odm,
                                gene_grna_group_pairs,
                                form,
                                B = 2500,
                                gene_modality_name = "gene",
                                grna_modality_name = "grna_expression",
                                grna_group_column_name = "grna_group",
                                n_pairs_to_sample = 50,
                                full_output = FALSE) {
  # DELETE AFTER REWRITING ASSIGN GRNA FUNCT
  grna_odm <- mm_odm |> get_modality("grna_expression")

  # step 1: check inputs; get the unique genes
  cat("Checking inputs. ")
  mm_odm <- check_ondisc_inputs(mm_odm = mm_odm,
                                gene_grna_group_pairs = gene_grna_group_pairs,
                                form = form,
                                gene_modality_name = gene_modality_name,
                                grna_modality_name = grna_modality_name,
                                grna_group_column_name = grna_group_column_name,
                                low_moi = TRUE)
  cat(crayon::green(' \u2713\n'))

  # step 2: assign gRNAs to cells (REWRITE TO ELIMINATE DEP ON LOWMOI AND IMPROVE MEM EFFICIENCY)
  cat("Obtaining the cell-to-gRNA assignments.")
  grna_group_info <- lowmoi::get_target_assignments_via_max_op(grna_odm) |> get_grna_group_info()
  rm(grna_odm)
  cat(crayon::green(' \u2713\n'))

  # step 3: perform the pairwise association test
  results <- perform_association_test_lowmoi_odm(mm_odm = mm_odm,
                                                 grna_group_info = grna_group_info,
                                                 gene_grna_group_pairs = gene_grna_group_pairs,
                                                 B = B,
                                                 full_output = full_output,
                                                 side = side)
  return(results)
}
