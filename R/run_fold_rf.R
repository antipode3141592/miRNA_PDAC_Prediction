run_fold_rf <- function(fold_inputs, ntree = 500, mtry_grid = c(5L, 10L, 20L, 30L)) {
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
      mtry = mtry_value
    )
  })

  oob_error <- vapply(
    rf_fits,
    function(rf_fit) as.numeric(utils::tail(rf_fit$err.rate[, "OOB"], n = 1)),
    numeric(1)
  )

  best_fit_index <- which.min(oob_error)
  rf_fit <- rf_fits[[best_fit_index]]
  best_mtry <- candidate_mtry[[best_fit_index]]

  prob <- as.numeric(
    predict(rf_fit, newdata = fold_inputs$valid_design_dat, type = "prob")[, "PDAC case"]
  )

  selected_variables <- fold_inputs$design_feature_cols
  selected_miRNA <- selected_variables[selected_variables %in% fold_inputs$selected_miRNA]

  list(
    prob = prob,
    selected_variables = selected_variables,
    selected_miRNA = selected_miRNA,
    tuning = data.frame(
      model = "random_forest",
      n_selected_miRNA = length(selected_miRNA),
      ntree = ntree,
      mtry = best_mtry,
      oob_error = oob_error[[best_fit_index]],
      stringsAsFactors = FALSE
    )
  )
}