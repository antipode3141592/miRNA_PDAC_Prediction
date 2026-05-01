scale_fold_predictor_data <- function(train_dat, valid_dat, miRNA_cols) {
  if (!length(miRNA_cols)) {
    return(list(train = train_dat, valid = valid_dat))
  }

  train_scaled <- train_dat
  valid_scaled <- valid_dat

  train_matrix_raw <- as.matrix(train_dat[, miRNA_cols, drop = FALSE])
  train_center <- colMeans(train_matrix_raw)
  train_scale <- apply(train_matrix_raw, 2, stats::sd)
  train_scale[!is.finite(train_scale) | train_scale == 0] <- 1

  train_matrix <- scale(
    train_matrix_raw,
    center = train_center,
    scale = train_scale
  )
  valid_matrix <- scale(
    as.matrix(valid_dat[, miRNA_cols, drop = FALSE]),
    center = train_center,
    scale = train_scale
  )

  train_scaled[, miRNA_cols] <- train_matrix
  valid_scaled[, miRNA_cols] <- valid_matrix

  list(train = train_scaled, valid = valid_scaled)
}