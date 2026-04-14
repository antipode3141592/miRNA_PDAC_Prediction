run_fold_forward_step <- function(fold_inputs, max_steps = 10) {
  null_glm <- glm(
    fold_inputs$base_model_formula,
    data = fold_inputs$train_design_dat,
    family = binomial(),
    na.action = stats::na.fail
  )
  full_glm <- glm(
    fold_inputs$model_formula,
    data = fold_inputs$train_design_dat,
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
    predict(forward_glm, newdata = fold_inputs$valid_design_dat, type = "response")
  )

  selected_features <- attr(terms(forward_glm), "term.labels")
  selected_features <- selected_features[selected_features %in% fold_inputs$selected_miRNA]

  list(
    prob = prob,
    selected_features = selected_features,
    tuning = data.frame(
      model = "forward_step",
      n_selected_features = length(selected_features),
      n_terms = length(selected_features),
      stringsAsFactors = FALSE
    )
  )
}