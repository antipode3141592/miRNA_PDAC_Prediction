run_fold_rf <- function(fold_inputs, ntree = 500) {
  rf_fit <- randomForest::randomForest(
    fold_inputs$class_formula,
    data = fold_inputs$train_design_dat,
    ntree = ntree
  )

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
      stringsAsFactors = FALSE
    )
  )
}