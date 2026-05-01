select_lasso_hyperparameters <- function(fold_inputs, alpha_lasso = 1) {
  alpha_grid <- sort(unique(as.numeric(alpha_lasso)))
  alpha_grid <- alpha_grid[is.finite(alpha_grid)]

  if (!length(alpha_grid)) {
    stop("alpha_lasso must contain at least one finite numeric value.")
  }

  if (any(alpha_grid <= 0 | alpha_grid > 1)) {
    stop("alpha_lasso values must be in (0, 1].")
  }

  cv_fits <- lapply(alpha_grid, function(alpha_value) {
    glmnet::cv.glmnet(
      fold_inputs$x_train,
      fold_inputs$y_train,
      alpha = alpha_value,
      family = "binomial",
      penalty.factor = fold_inputs$penalty_factor
    )
  })
  cv_scores <- vapply(cv_fits, function(cv_fit) min(cv_fit$cvm), numeric(1))
  best_index <- which.min(cv_scores)

  data.frame(
    model = "lasso",
    alpha = alpha_grid[[best_index]],
    lambda = cv_fits[[best_index]]$lambda.min,
    cv_score = cv_scores[[best_index]],
    alpha_grid = paste(alpha_grid, collapse = ","),
    tuning_source = "fixed_from_full_training_set",
    row.names = NULL,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

run_fold_lasso <- function(
  fold_inputs,
  alpha_lasso = 1,
  lambda_lasso,
  tuning_cv_score = NA_real_,
  alpha_grid = NULL,
  tuning_source = "fixed_from_full_training_set"
) {
  if (missing(lambda_lasso) || is.null(lambda_lasso)) {
    selected_tuning <- select_lasso_hyperparameters(
      fold_inputs,
      alpha_lasso = alpha_lasso
    )

    alpha_lasso <- selected_tuning$alpha[[1]]
    lambda_lasso <- selected_tuning$lambda[[1]]
    tuning_cv_score <- selected_tuning$cv_score[[1]]
    alpha_grid <- selected_tuning$alpha_grid[[1]]
    tuning_source <- "nested_cv_within_fold"
  }

  alpha_value <- as.numeric(alpha_lasso[[1]])

  if (!is.finite(alpha_value) || alpha_value <= 0 || alpha_value > 1) {
    stop("alpha_lasso must be a single numeric value in (0, 1].", call. = TRUE)
  }

  lambda_value <- as.numeric(lambda_lasso[[1]])

  if (!is.finite(lambda_value) || lambda_value <= 0) {
    stop("lambda_lasso must be a single positive numeric value.", call. = TRUE)
  }

  lasso_fit <- glmnet::glmnet(
    fold_inputs$x_train,
    fold_inputs$y_train,
    alpha = alpha_value,
    family = "binomial",
    penalty.factor = fold_inputs$penalty_factor,
    lambda = lambda_value
  )

  prob <- as.numeric(
    predict(lasso_fit, s = lambda_value, newx = fold_inputs$x_valid, type = "response")
  )

  coef_mat_lasso <- as.matrix(
    predict(lasso_fit, type = "coefficients", s = lambda_value)
  )

  selected_variables <- rownames(coef_mat_lasso)[coef_mat_lasso[, 1] != 0]
  selected_variables <- setdiff(selected_variables, "(Intercept)")
  selected_miRNA <- selected_variables[selected_variables %in% fold_inputs$selected_miRNA]

  list(
    prob = prob,
    selected_variables = selected_variables,
    selected_miRNA = selected_miRNA,
    final_alpha = alpha_value,
    final_lambda = lambda_value,
    tuning = data.frame(
      model = "lasso",
      n_selected_miRNA = length(selected_miRNA),
      alpha = alpha_value,
      lambda = lambda_value,
      cv_score = as.numeric(tuning_cv_score[[1]]),
      alpha_grid = if (is.null(alpha_grid)) NA_character_ else as.character(alpha_grid[[1]]),
      tuning_source = tuning_source,
      stringsAsFactors = FALSE
    )
  )
}