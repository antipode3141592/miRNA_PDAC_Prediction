run_fold_dt <- function(fold_inputs, cp_dt = 0.01, maxdepth = 10) {
  dt_fit <- rpart::rpart(
    fold_inputs$class_formula,
    data = fold_inputs$train_model_dat,
    method = "class",
    control = rpart::rpart.control(cp = cp_dt, maxdepth = maxdepth)
  )

  prob <- as.numeric(
    predict(dt_fit, newdata = fold_inputs$valid_model_dat, type = "prob")[, "PDAC case"]
  )

  selected_features <- attr(terms(dt_fit), "term.labels")

  list(
    prob = prob,
    selected_features = selected_features,
    tuning = data.frame(
      model = "decision_tree",
      n_selected_features = length(selected_features),
      cp = cp_dt,
      maxdepth = maxdepth,
      n_terms = length(selected_features),
      stringsAsFactors = FALSE
    )
  )
}