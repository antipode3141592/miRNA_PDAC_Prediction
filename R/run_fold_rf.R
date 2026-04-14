run_fold_rf <- function(fold_inputs, ntree = 500) {
  rf_fit <- randomForest::randomForest(
    fold_inputs$class_formula,
    data = fold_inputs$train_design_dat,
    ntree = ntree
  )

  prob <- as.numeric(
    predict(rf_fit, newdata = fold_inputs$valid_design_dat, type = "prob")[, "PDAC case"]
  )

  list(
    prob = prob,
    selected_features = fold_inputs$selected_miRNA,
    tuning = data.frame(
      model = "random_forest",
      n_selected_features = length(fold_inputs$selected_miRNA),
      ntree = ntree,
      stringsAsFactors = FALSE
    )
  )
}