differential_expression_filter_fn <- function(
  train_fold_dat,
  miRNA_cols,
  filter_stat_cols,
  top_n = 50,
  group_term = "CC",
  adjustment_covariates = c("Sex", "Race", "Age"),
  coef_name = "CCPDAC case"
) {
  if (!length(miRNA_cols)) {
    stop("No miRNA columns were provided to the blocked limma filter.", call. = TRUE)
  }

  if (!length(filter_stat_cols)) {
    stop("filter_stat_cols must be provided to the blocked limma filter.", call. = TRUE)
  }

  case_expr <- train_fold_dat[train_fold_dat$CC == "PDAC case", miRNA_cols, drop = FALSE]
  control_expr <- train_fold_dat[train_fold_dat$CC == "Controls", miRNA_cols, drop = FALSE]

  case_stats <- data.frame(
    miRNA = miRNA_cols,
    mean_case = colMeans(case_expr, na.rm = TRUE),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  control_stats <- data.frame(
    miRNA = miRNA_cols,
    mean_control = colMeans(control_expr, na.rm = TRUE),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  limma_runner <- match.fun("run_limma_blocked")

  limma_out <- limma_runner(
    train_dat = train_fold_dat,
    selected_miRNA = miRNA_cols,
    group_term = group_term,
    adjustment_covariates = adjustment_covariates,
    coef_name = coef_name
  )

  limma_results <- limma_out$results

  class_stats <- data.frame(
    miRNA = limma_results$miRNA,
    limma_log_fc = limma_results$logFC,
    limma_avg_expr = limma_results$AveExpr,
    limma_t = limma_results$t,
    limma_p_value = limma_results$P.Value,
    limma_adj_p_value = limma_results$adj.P.Val,
    limma_b = limma_results$B,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  class_stats <- dplyr::left_join(class_stats, control_stats, by = "miRNA")
  class_stats <- dplyr::left_join(class_stats, case_stats, by = "miRNA")
  class_stats$mean_diff <- class_stats$mean_case - class_stats$mean_control

  keep_rows <-
    is.finite(class_stats$limma_t) &
    is.finite(class_stats$limma_p_value) &
    is.finite(class_stats$mean_control) &
    is.finite(class_stats$mean_case)

  class_stats <- class_stats[keep_rows, , drop = FALSE]

  missing_filter_cols <- setdiff(filter_stat_cols, names(class_stats))
  if (length(missing_filter_cols)) {
    stop(
      sprintf(
        "filter_stat_cols included columns not produced by the blocked limma filter: %s",
        paste(missing_filter_cols, collapse = ", ")
      ),
      call. = TRUE
    )
  }

  class_stats <- class_stats[, filter_stat_cols, drop = FALSE]
  class_stats <- class_stats[order(class_stats$limma_p_value, -abs(class_stats$limma_t)), , drop = FALSE]

  if (!nrow(class_stats)) {
    stop("Blocked limma filter returned zero usable miRNAs.", call. = TRUE)
  }

  selected <- class_stats |>
    dplyr::slice_head(n = top_n) |>
    dplyr::pull("miRNA")

  if (!length(selected)) {
    stop("Blocked limma ranking yielded zero selected miRNAs.", call. = TRUE)
  }

  list(
    stats = class_stats,
    selected = selected
  )
}