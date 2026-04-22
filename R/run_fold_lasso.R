run_fold_lasso <- function(fold_inputs, alpha_lasso = 1) {
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

  cv.lasso <- cv_fits[[best_index]]
  best_alpha_lasso <- alpha_grid[[best_index]]
  bestlam_lasso <- cv.lasso$lambda.min
  prob <- as.numeric(
    predict(cv.lasso, s = bestlam_lasso, newx = fold_inputs$x_valid, type = "response")
  )

  coef_mat_lasso <- as.matrix(
    predict(cv.lasso, type = "coefficients", s = bestlam_lasso)
  )

  selected_variables <- rownames(coef_mat_lasso)[coef_mat_lasso[, 1] != 0]
  selected_variables <- setdiff(selected_variables, "(Intercept)")
  selected_miRNA <- selected_variables[selected_variables %in% fold_inputs$selected_miRNA]

  list(
    prob = prob,
    selected_variables = selected_variables,
    selected_miRNA = selected_miRNA,
    final_alpha = best_alpha_lasso,
    tuning = data.frame(
      model = "lasso",
      n_selected_miRNA = length(selected_miRNA),
      alpha = best_alpha_lasso,
      lambda = bestlam_lasso,
      cv_score = cv_scores[[best_index]],
      alpha_grid = paste(alpha_grid, collapse = ","),
      stringsAsFactors = FALSE
    )
  )
}