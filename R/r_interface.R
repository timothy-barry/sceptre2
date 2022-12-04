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
#' @param output_amount return the full output (2) or a simplified, reduced output (1)?
#'
#' @return a data frame with columns `response_id`, `grna_group`, `p_value`, and `log_fold_change`.
#' @export
#'
#' @examples
#' \dontrun{
#' # load ondisc
#' library(ondisc)
#'
#' ##########################
#' # EXAMPLE ON FRANGIEH DATA
#' ##########################
#' frangieh_dir <- paste0(.get_config_path("LOCAL_SCEPTRE2_DATA_DIR"),
#' "data/frangieh/co_culture/")
#'
#' # obtain the multimodal odm
#' mm_fp <- paste0(frangieh_dir, "multimodal_metadata.rds")
#' odm_fps <- paste0(frangieh_dir, c("gene/matrix.odm", "grna_assignment/matrix.odm"))
#' mm_odm <- ondisc::read_multimodal_odm(odm_fps = odm_fps, multimodal_metadata_fp = mm_fp)
#'
#' response_grna_group_pairs <- expand.grid(response_id = mm_odm |>
#' get_modality("gene") |>
#' get_feature_ids() |>
#' sample(5),
#' grna_group = c("CDKN1A", "ISYNA1"))
#'
#' # form <- formula(~ log(gene_n_nonzero) + log(gene_n_umis) + batch)
#' form <- formula(~ log(gene_n_nonzero) + log(gene_n_umis))
#' response_modality_name <- "gene"
#' grna_modality_name <- "grna_assignment"
#' grna_group_column_name <- "target"
#' B <- 250000
#' side <- "both"
#' output_amount <- 1
#' max_b_per_batch <- 250000
#' in_memory <- TRUE
#' statistic <- "full"
#'
#' # call the function
#' result <- run_sceptre_low_moi(mm_odm,
#'                               response_grna_group_pairs,
#'                               form,
#'                               response_modality_name,
#'                               grna_modality_name,
#'                               grna_group_column_name,
#'                               B,
#'                               side,
#'                               output_amount,
#'                               max_b_per_batch,
#'                               in_memory,
#'                               statistic)
#' }
run_sceptre_low_moi <- function(mm_odm,
                                response_grna_group_pairs,
                                form,
                                response_modality_name = "gene",
                                grna_modality_name = "grna",
                                grna_group_column_name = "grna_group",
                                B = 2500,
                                side = "both",
                                output_amount = 1,
                                max_b_per_batch = 250000,
                                in_memory = TRUE,
                                statistic = "full") {
  # DELETE AFTER REWRITING ASSIGN GRNA FUNCT
  grna_odm <- mm_odm |> ondisc::get_modality(grna_modality_name)

  # step 1: check inputs; get the unique responses
  cat("Checking inputs. ")
  if (!output_amount %in% c(1, 2, 3)) stop("output_amount not recognized.")
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
  input_grna_groups <- as.character(unique(response_grna_group_pairs$grna_group))
  grna_group_assignments <- lowmoi::get_target_assignments_via_max_op(grna_odm)
  grna_group_info <- get_grna_group_info(grna_group_assignments, input_grna_groups)
  rm(grna_odm)
  cat(crayon::green(' \u2713\n'))

  # step 3: perform the pairwise association test
  results <- perform_association_test_lowmoi_odm_v3(mm_odm = mm_odm,
                                                    grna_group_info = grna_group_info,
                                                    response_grna_group_pairs = response_grna_group_pairs,
                                                    B = B,
                                                    output_amount = output_amount,
                                                    side = side,
                                                    max_b_per_batch = max_b_per_batch,
                                                    in_memory = in_memory,
                                                    statistic = statistic)
  return(results)
}
