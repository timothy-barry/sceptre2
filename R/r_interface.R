#' Run sceptre (low MOI)
#'
#' Runs sceptre on a low MOI dataset.
#'
#' @param mm_odm a multimodal ondisc matrix containing response and grna (either expression or assignment) modalities
#' @param response_grna_group_pairs a data frame with columns `response_id` and `grna_group` giving the response ID / grna group pairs to analyze.
#' @param form a formula object (or string) specifying how to adjust for the covariates
#' @param B number of resamples to draw for permutation/conditional randomization test
#' @param response_modality_name name of the response modality (e.g., "gene") within the multimodal ODM
#' @param grna_modality_name name of the grna modality (e.g., "grna_expression") within the multimodal ODM
#' @param grna_group_column_name name of the column within the feature covariate matrix of the gRNA ODM that contains the gRNA groups
#' @param side sidedness of the test (one of "left", "both", and "right")
#' @param full_output return the full output (TRUE) or a simplified, reduced output (FALSE)?
#'
#' @return a data frame with columns `response_id`, `grna_group`, `p_value`, and `log_fold_change`.
#' @export
#'
#' @examples
#' \dontrun{
#' # load ondisc
#' library(ondisc)
#'
#' # set the tap seq dir
#' tap_seq_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
#' "data/schraivogel/ground_truth_tapseq/")
#'
#' # obtain the multimodal odm
#' mm_fp <- paste0(tap_seq_dir, "multimodal_metadata.rds")
#' odm_fps <- paste0(tap_seq_dir, c("gene/matrix.odm", "grna_expression/matrix.odm"))
#' mm_odm <- ondisc::read_multimodal_odm(odm_fps = odm_fps, multimodal_metadata_fp = mm_fp)
#'
#' # set the response grna group pairs
#' response_grna_group_pairs <- expand.grid(response_id = mm_odm |> get_modality("gene") |>
#' get_feature_ids() |> sample(3),
#' grna_group = mm_odm |> get_modality("grna_expression") |>
#' get_feature_covariates() |> dplyr::pull(target) |> unique()) |>
#' dplyr::filter(grna_group != "non-targeting")
#'
#' # set the arguments to the function
#' form <- formula(~ log(gene_n_nonzero) + log(gene_n_umis) + batch)
#' response_modality_name <- "gene"
#' grna_modality_name <- "grna_expression"
#' grna_group_column_name <- "target"
#' B <- 2500
#' side <- "both"
#' full_output <- FALSE
#'
#' # call the function
#' result <- run_sceptre_low_moi(mm_odm,
#' response_grna_group_pairs,
#' form,
#' response_modality_name,
#' grna_modality_name,
#' grna_group_column_name,
#' B,
#' side,
#' full_output)
#' }
run_sceptre_low_moi <- function(mm_odm,
                                response_grna_group_pairs,
                                form,
                                response_modality_name = "gene",
                                grna_modality_name = "grna",
                                grna_group_column_name = "grna_group",
                                B = 2500,
                                side = "both",
                                full_output = FALSE) {
  # DELETE AFTER REWRITING ASSIGN GRNA FUNCT
  grna_odm <- mm_odm |> ondisc::get_modality(grna_modality_name)

  # step 1: check inputs; get the unique responses
  cat("Checking inputs. ")
  mm_odm <- check_ondisc_inputs(mm_odm = mm_odm,
                                response_grna_group_pairs = response_grna_group_pairs,
                                form = form,
                                response_modality_name = response_modality_name,
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
                                                 response_grna_group_pairs = response_grna_group_pairs,
                                                 B = B,
                                                 full_output = full_output,
                                                 side = side)
  return(results)
}
