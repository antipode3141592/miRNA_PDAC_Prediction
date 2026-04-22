run_fold_rf <- function(
  fold_inputs,
  ntree = 500,
  mtry_grid = c(5L, 10L, 20L, 30L),
  top_m_features = NULL,
  importance_metric = c("MeanDecreaseAccuracy", "MeanDecreaseGini")
) {
  importance_metric <- match.arg(importance_metric)

  candidate_mtry <- unique(as.integer(mtry_grid))
  candidate_mtry <- candidate_mtry[
    is.finite(candidate_mtry) &
      candidate_mtry >= 1L &
      candidate_mtry <= length(fold_inputs$design_feature_cols)
  ]

  if (!length(candidate_mtry)) {
    candidate_mtry <- min(
      max(1L, floor(sqrt(length(fold_inputs$design_feature_cols)))),
      length(fold_inputs$design_feature_cols)
    )
  }

  rf_fits <- lapply(candidate_mtry, function(mtry_value) {
    randomForest::randomForest(
      fold_inputs$class_formula,
      data = fold_inputs$train_design_dat,
      ntree = ntree,
      mtry = mtry_value,
      importance = TRUE
    )
  })

  oob_error <- vapply(
    rf_fits,
    function(rf_fit) as.numeric(utils::tail(rf_fit$err.rate[, "OOB"], n = 1)),
    numeric(1)
  )

  best_fit_index <- which.min(oob_error)
  rf_search_fit <- rf_fits[[best_fit_index]]
  best_mtry <- candidate_mtry[[best_fit_index]]
  oob_error_curve <- data.frame(
    n_trees = seq_len(nrow(rf_search_fit$err.rate)),
    oob_error = as.numeric(rf_search_fit$err.rate[, "OOB"]),
    mtry = best_mtry,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  optimal_ntree <- oob_error_curve$n_trees[[which.min(oob_error_curve$oob_error)]]

  rf_fit <- if (optimal_ntree < ntree) {
    randomForest::randomForest(
      fold_inputs$class_formula,
      data = fold_inputs$train_design_dat,
      ntree = optimal_ntree,
      mtry = best_mtry,
      importance = TRUE
    )
  } else {
    rf_search_fit
  }

  final_oob_error <- as.numeric(utils::tail(rf_fit$err.rate[, "OOB"], n = 1))

  importance_matrix <- randomForest::importance(rf_fit, scale = FALSE)
  feature_importance <- if (is.null(dim(importance_matrix))) {
    data.frame(
      variable = names(importance_matrix),
      raw_importance = as.numeric(importance_matrix),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    data.frame(
      variable = rownames(importance_matrix),
      importance_matrix,
      row.names = NULL,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  if (!(importance_metric %in% names(feature_importance))) {
    fallback_metric <- c("MeanDecreaseAccuracy", "MeanDecreaseGini")
    fallback_metric <- fallback_metric[fallback_metric %in% names(feature_importance)]

    if (!length(fallback_metric)) {
      stop(
        "Random forest importance output did not include a supported ranking metric.",
        call. = TRUE
      )
    }

    importance_metric <- fallback_metric[[1]]
  }

  if (!("MeanDecreaseAccuracy" %in% names(feature_importance))) {
    feature_importance$MeanDecreaseAccuracy <- NA_real_
  }

  if (!("MeanDecreaseGini" %in% names(feature_importance))) {
    feature_importance$MeanDecreaseGini <- NA_real_
  }

  feature_importance$variable_type <- ifelse(
    feature_importance$variable %in% fold_inputs$selected_miRNA,
    "miRNA",
    "demographic"
  )
  feature_importance$importance_metric <- importance_metric
  feature_importance$importance_value <- as.numeric(feature_importance[[importance_metric]])
  feature_importance <- feature_importance[
    order(-feature_importance$importance_value, feature_importance$variable),
    c(
      "variable",
      "variable_type",
      "importance_metric",
      "importance_value",
      "MeanDecreaseAccuracy",
      "MeanDecreaseGini"
    ),
    drop = FALSE
  ]
  rownames(feature_importance) <- NULL
  feature_importance$importance_rank <- seq_len(nrow(feature_importance))

  miRNA_feature_importance <- feature_importance[
    feature_importance$variable_type == "miRNA",
    ,
    drop = FALSE
  ]

  if (!nrow(miRNA_feature_importance)) {
    stop("Random forest importance ranking produced zero miRNA features.", call. = TRUE)
  }

  if (is.null(top_m_features)) {
    n_selected_miRNA <- nrow(miRNA_feature_importance)
  } else {
    top_m_features <- as.integer(top_m_features[[1]])

    if (!is.finite(top_m_features) || top_m_features < 1L) {
      stop("top_m_features must be NULL or a single positive integer.", call. = TRUE)
    }

    n_selected_miRNA <- min(top_m_features, nrow(miRNA_feature_importance))
  }

  selected_miRNA <- utils::head(miRNA_feature_importance$variable, n_selected_miRNA)
  selected_variables <- unique(c(fold_inputs$base_feature_cols, selected_miRNA))
  feature_importance$selected_for_final <- feature_importance$variable %in% selected_variables
  feature_importance <- feature_importance[
    ,
    c(
      "importance_rank",
      "variable",
      "variable_type",
      "importance_metric",
      "importance_value",
      "MeanDecreaseAccuracy",
      "MeanDecreaseGini",
      "selected_for_final"
    ),
    drop = FALSE
  ]

  prob <- as.numeric(
    predict(rf_fit, newdata = fold_inputs$valid_design_dat, type = "prob")[, "PDAC case"]
  )

  list(
    prob = prob,
    selected_variables = selected_variables,
    selected_miRNA = selected_miRNA,
    oob_error_curve = oob_error_curve,
    optimal_ntree = optimal_ntree,
    feature_importance = feature_importance,
    tuning = data.frame(
      model = "random_forest",
      n_selected_miRNA = length(selected_miRNA),
      ntree = optimal_ntree,
      ntree_search_max = ntree,
      mtry = best_mtry,
      optimal_ntree = optimal_ntree,
      oob_error = final_oob_error,
      oob_error_search_min = min(oob_error_curve$oob_error, na.rm = TRUE),
      importance_metric = importance_metric,
      top_m_features = length(selected_miRNA),
      stringsAsFactors = FALSE
    )
  )
}