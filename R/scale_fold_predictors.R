scale_fold_predictor_data <- function(train_dat, valid_dat, miRNA_cols) {
  if (!length(miRNA_cols)) {
    return(list(train = train_dat, valid = valid_dat))
  }

  train_scaled <- train_dat
  valid_scaled <- valid_dat

  train_matrix <- scale(as.matrix(train_dat[, miRNA_cols, drop = FALSE]))
  valid_matrix <- scale(
    as.matrix(valid_dat[, miRNA_cols, drop = FALSE]),
    center = attr(train_matrix, "scaled:center"),
    scale = attr(train_matrix, "scaled:scale")
  )

  train_scaled[, miRNA_cols] <- train_matrix
  valid_scaled[, miRNA_cols] <- valid_matrix

  list(train = train_scaled, valid = valid_scaled)
}