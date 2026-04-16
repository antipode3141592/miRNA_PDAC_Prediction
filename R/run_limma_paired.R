run_limma_blocked <- function(
  train_dat,
  selected_miRNA,
  group_term = "CC",
  adjustment_covariates = c("Sex", "Race", "Age"),
  coef_name = "CCPDAC case"
) {
  if (length(selected_miRNA) == 0) {
    return(
      list(
        results = data.frame(),
        metadata = data.frame(),
        expr = matrix(nrow = 0, ncol = 0),
        design = NULL,
        correlation = NULL
      )
    )
  }

  metadata_train <- data.frame(
    ID = as.character(train_dat$ID),
    matched_set = factor(train_dat$matched_set),
    CC = factor(train_dat$CC, levels = c("Controls", "PDAC case")),
    Sex = factor(train_dat$Sex, levels = c("Female", "Male")),
    Race = factor(train_dat$Race, levels = c("Asian", "Black", "White", "Others")),
    Age = as.numeric(train_dat$Age),
    row.names = NULL
  )

  expr_mat <- t(as.matrix(train_dat[, selected_miRNA, drop = FALSE]))
  colnames(expr_mat) <- train_dat$ID

  metadata_train <- metadata_train[match(colnames(expr_mat), metadata_train$ID), , drop = FALSE]
  metadata_train <- droplevels(metadata_train)

  if (nrow(metadata_train) == 0 || ncol(expr_mat) == 0) {
    stop("Blocked limma requires non-empty training data and miRNA features.", call. = TRUE)
  }

  if (nlevels(metadata_train$CC) < 2) {
    stop("Blocked limma requires both outcome classes.", call. = TRUE)
  }

  design_formula <- stats::reformulate(c(group_term, adjustment_covariates))
  design_limma <- model.matrix(design_formula, data = metadata_train)

  if (!(coef_name %in% colnames(design_limma))) {
    stop(sprintf("Coefficient %s is not present in the limma design matrix.", coef_name), call. = TRUE)
  }

  corfit <- limma::duplicateCorrelation(
    expr_mat,
    design = design_limma,
    block = metadata_train$matched_set
  )

  fit_limma <- limma::lmFit(
    expr_mat,
    design = design_limma,
    block = metadata_train$matched_set,
    correlation = corfit$consensus
  )
  fit_limma <- limma::eBayes(fit_limma)

  limma_results <- limma::topTable(
    fit_limma,
    coef = coef_name,
    number = Inf,
    sort.by = "P"
  )
  limma_results$miRNA <- rownames(limma_results)

  list(
    results = limma_results,
    metadata = metadata_train,
    expr = expr_mat,
    design = design_limma,
    correlation = corfit$consensus
  )
}