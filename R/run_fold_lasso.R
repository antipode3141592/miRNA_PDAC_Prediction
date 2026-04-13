run_fold_lasso <- function(fold_inputs, alpha_lasso = 1) {
  cv.lasso <- glmnet::cv.glmnet(
    fold_inputs$x_train,
    fold_inputs$y_train,
    alpha = alpha_lasso,
    family = "binomial"
  )

  bestlam_lasso <- cv.lasso$lambda.min
  prob <- as.numeric(
    predict(cv.lasso, s = bestlam_lasso, newx = fold_inputs$x_valid, type = "response")
  )

  coef_mat_lasso <- as.matrix(
    predict(cv.lasso, type = "coefficients", s = bestlam_lasso)
  )

  selected_features <- rownames(coef_mat_lasso)[coef_mat_lasso[, 1] != 0]
  selected_features <- setdiff(selected_features, "(Intercept)")

  list(
    prob = prob,
    selected_features = selected_features,
    tuning = data.frame(
      model = "lasso",
      n_selected_features = length(selected_features),
      alpha = alpha_lasso,
      lambda = bestlam_lasso,
      stringsAsFactors = FALSE
    )
  )
}