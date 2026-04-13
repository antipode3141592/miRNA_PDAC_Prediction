run_fold_enet <- function(
  fold_inputs,
  alpha_enet = c(0.1, 0.3, 0.5, 0.7, 0.9)
) {
  alpha_grid <- sort(unique(alpha_enet))
  cv_fits <- lapply(alpha_grid, function(alpha_value) {
    glmnet::cv.glmnet(
      fold_inputs$x_train,
      fold_inputs$y_train,
      alpha = alpha_value,
      family = "binomial"
    )
  })
  cv_scores <- vapply(cv_fits, function(cv_fit) min(cv_fit$cvm), numeric(1))
  best_index <- which.min(cv_scores)

  cv.enet <- cv_fits[[best_index]]
  best_alpha_enet <- alpha_grid[[best_index]]
  bestlam_enet <- cv.enet$lambda.min
  prob <- as.numeric(
    predict(cv.enet, s = bestlam_enet, newx = fold_inputs$x_valid, type = "response")
  )

  coef_mat_enet <- as.matrix(
    predict(cv.enet, type = "coefficients", s = bestlam_enet)
  )

  selected_features <- rownames(coef_mat_enet)[coef_mat_enet[, 1] != 0]
  selected_features <- setdiff(selected_features, "(Intercept)")

  list(
    prob = prob,
    selected_features = selected_features,
    tuning = data.frame(
      model = "enet",
      n_selected_features = length(selected_features),
      alpha = best_alpha_enet,
      lambda = bestlam_enet,
      cv_score = cv_scores[[best_index]],
      alpha_grid = paste(alpha_grid, collapse = ","),
      stringsAsFactors = FALSE
    )
  )
}