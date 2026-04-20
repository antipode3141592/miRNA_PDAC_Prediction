filter_fold_paired_snr <- function(train_fold_dat, miRNA_cols, top_n = 50) {
  pair_data <- train_fold_dat[, c("matched_set", "CC", miRNA_cols), drop = FALSE]
  pair_data$group <- ifelse(pair_data$CC == "PDAC case", "case", "control")

  pair_counts <- stats::aggregate(
    cbind(
      n_rows = 1,
      n_control = as.integer(pair_data$group == "control"),
      n_case = as.integer(pair_data$group == "case")
    ),
    by = list(matched_set = pair_data$matched_set),
    FUN = sum
  )

  valid_sets <- pair_counts$matched_set[
    pair_counts$n_rows == 2 &
      pair_counts$n_control == 1 &
      pair_counts$n_case == 1
  ]

  pair_data <- pair_data[
    pair_data$matched_set %in% valid_sets,
    c("matched_set", "group", miRNA_cols),
    drop = FALSE
  ]

  pair_data <- tidyr::pivot_wider(
    pair_data,
    names_from = tidyselect::all_of("group"),
    values_from = tidyselect::all_of(miRNA_cols),
    names_sep = "__"
  )

  if (nrow(pair_data) == 0) {
    return(list(stats = data.frame(), selected = character(0)))
  }

  diff_mat <- vapply(
    miRNA_cols,
    function(miRNA) {
      pair_data[[paste0(miRNA, "__case")]] - pair_data[[paste0(miRNA, "__control")]]
    },
    numeric(nrow(pair_data))
  )

  diff_mat <- as.data.frame(diff_mat)
  names(diff_mat) <- miRNA_cols
  n_pairs <- nrow(diff_mat)

  paired_stats <- data.frame(
    miRNA = miRNA_cols,
    mean_diff = colMeans(diff_mat, na.rm = TRUE),
    sd_diff = apply(diff_mat, 2, sd, na.rm = TRUE),
    row.names = NULL
  )

  paired_stats$paired_snr <- abs(paired_stats$mean_diff) / paired_stats$sd_diff
  paired_stats$t_stat <- paired_stats$mean_diff / (paired_stats$sd_diff / sqrt(n_pairs))
  paired_stats$p_value <- 2 * pt(-abs(paired_stats$t_stat), df = n_pairs - 1)

  paired_stats <- paired_stats[
    stats::complete.cases(paired_stats[, c("paired_snr", "p_value")]),
    ,
    drop = FALSE
  ]
  paired_stats <- paired_stats[order(-paired_stats$paired_snr, paired_stats$p_value), , drop = FALSE]

  selected <- head(paired_stats$miRNA, top_n)

  list(
    stats = paired_stats,
    selected = selected
  )
}