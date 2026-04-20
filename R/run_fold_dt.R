run_fold_dt <- function(fold_inputs, cp_dt = 0.001, maxdepth = 20) {
  dt_fit <- rpart::rpart(
    fold_inputs$class_formula,
    data = fold_inputs$train_design_dat,
    method = "class",
    control = rpart::rpart.control(cp = cp_dt, maxdepth = maxdepth)
  )

  prob <- as.numeric(
    predict(dt_fit, newdata = fold_inputs$valid_design_dat, type = "prob")[, "PDAC case"]
  )

  selected_variables <- attr(terms(dt_fit), "term.labels")
  selected_miRNA <- selected_variables[selected_variables %in% fold_inputs$selected_miRNA]

  list(
    prob = prob,
    selected_variables = selected_variables,
    selected_miRNA = selected_miRNA,
    tuning = data.frame(
      model = "decision_tree",
      n_selected_miRNA = length(selected_miRNA),
      cp = cp_dt,
      maxdepth = maxdepth,
      n_terms = length(selected_variables),
      stringsAsFactors = FALSE
    )
  )
}