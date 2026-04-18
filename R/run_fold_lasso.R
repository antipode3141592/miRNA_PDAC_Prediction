run_fold_lasso <- function(fold_inputs, alpha_lasso = 1) {
  cv.lasso <- glmnet::cv.glmnet(
    fold_inputs$x_train,
    fold_inputs$y_train,
    alpha = alpha_lasso,
    family = "binomial",
    penalty.factor = fold_inputs$penalty_factor
  )

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
    tuning = data.frame(
      model = "lasso",
      n_selected_miRNA = length(selected_miRNA),
      alpha = alpha_lasso,
      lambda = bestlam_lasso,
      stringsAsFactors = FALSE
    )
  )
}