run_fold_forward_step <- function(fold_inputs, max_steps = 10) {
  scaled_design_data <- scale_fold_predictor_data(
    fold_inputs$train_design_dat,
    fold_inputs$valid_design_dat,
    fold_inputs$selected_miRNA
  )
  train_design_dat <- scaled_design_data$train
  valid_design_dat <- scaled_design_data$valid

  null_glm <- glm(
    fold_inputs$base_model_formula,
    data = train_design_dat,
    family = binomial(),
    na.action = stats::na.fail
  )
  full_glm <- glm(
    fold_inputs$model_formula,
    data = train_design_dat,
    family = binomial(),
    na.action = stats::na.fail
  )

  forward_glm <- step(
    object = null_glm,
    scope = list(lower = formula(null_glm), upper = formula(full_glm)),
    direction = "forward",
    trace = 0,
    steps = min(max_steps, length(fold_inputs$selected_miRNA))
  )

  prob <- as.numeric(
    predict(forward_glm, newdata = valid_design_dat, type = "response")
  )

  selected_variables <- attr(terms(forward_glm), "term.labels")
  selected_miRNA <- selected_variables[selected_variables %in% fold_inputs$selected_miRNA]

  list(
    prob = prob,
    selected_variables = selected_variables,
    selected_miRNA = selected_miRNA,
    tuning = data.frame(
      model = "forward_step",
      n_selected_miRNA = length(selected_miRNA),
      n_terms = length(selected_variables),
      scaled_predictors = TRUE,
      stringsAsFactors = FALSE
    )
  )
}