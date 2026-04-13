run_fold_lda <- function(
  fold_inputs,
  max_features = 20,
  prior = c(0.5, 0.5)
) {
  selected_miRNA <- head(fold_inputs$selected_miRNA, max_features)
  class_formula <- stats::reformulate(selected_miRNA, response = "CC")
  train_model_dat <- fold_inputs$train_model_dat[, c("CC", "CC_bin", selected_miRNA), drop = FALSE]
  valid_model_dat <- fold_inputs$valid_model_dat[, c("CC", "CC_bin", selected_miRNA), drop = FALSE]

  lda_fit <- MASS::lda(
    class_formula,
    data = train_model_dat,
    prior = prior
  )
  lda_pred <- predict(lda_fit, newdata = valid_model_dat)
  prob <- as.numeric(lda_pred$posterior[, "PDAC case"])

  list(
    prob = prob,
    selected_features = selected_miRNA,
    tuning = data.frame(
      model = "lda",
      n_selected_features = length(selected_miRNA),
      n_features = length(selected_miRNA),
      prior = paste(prior, collapse = ","),
      stringsAsFactors = FALSE
    )
  )
}