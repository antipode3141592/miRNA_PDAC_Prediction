select_top_limma_miRNA <- function(limma_stats, top_n = NULL, p_value_threshold = NULL) {
  ranked_stats <- limma_stats[
    order(limma_stats$limma_p_value, -abs(limma_stats$limma_t), limma_stats$miRNA),
    ,
    drop = FALSE
  ]

  if (!nrow(ranked_stats)) {
    return(ranked_stats)
  }

  if (!is.null(p_value_threshold)) {
    if (!is.numeric(p_value_threshold) || length(p_value_threshold) != 1L ||
        !is.finite(p_value_threshold) || p_value_threshold <= 0 || p_value_threshold > 1) {
      stop("p_value_threshold must be a single numeric value in the interval (0, 1].", call. = TRUE)
    }

    return(ranked_stats[ranked_stats$limma_p_value < p_value_threshold, , drop = FALSE])
  }

  if (is.null(top_n) || !is.numeric(top_n) || length(top_n) != 1L ||
      !is.finite(top_n) || top_n < 1L) {
    stop("top_n must be a single positive integer when p_value_threshold is not provided.", call. = TRUE)
  }

  cutoff_index <- min(top_n, nrow(ranked_stats))
  cutoff_p_value <- ranked_stats$limma_p_value[[cutoff_index]]

  ranked_stats[ranked_stats$limma_p_value <= cutoff_p_value, , drop = FALSE]
}

rank_blocked_limma_miRNA <- function(limma_stats) {
  select_top_limma_miRNA(limma_stats, top_n = nrow(limma_stats))$miRNA
}