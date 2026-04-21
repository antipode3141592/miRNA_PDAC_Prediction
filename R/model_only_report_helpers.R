select_top_limma_miRNA <- function(limma_stats, top_n) {
  ranked_stats <- limma_stats[
    order(limma_stats$limma_p_value, -abs(limma_stats$limma_t), limma_stats$miRNA),
    ,
    drop = FALSE
  ]

  if (!nrow(ranked_stats)) {
    return(ranked_stats)
  }

  cutoff_index <- min(top_n, nrow(ranked_stats))
  cutoff_p_value <- ranked_stats$limma_p_value[[cutoff_index]]

  ranked_stats[ranked_stats$limma_p_value <= cutoff_p_value, , drop = FALSE]
}

rank_blocked_limma_miRNA <- function(limma_stats) {
  select_top_limma_miRNA(limma_stats, top_n = nrow(limma_stats))$miRNA
}

build_binary_roc <- function(observed, predicted_prob) {
  if (!length(observed) || !length(predicted_prob)) {
    stop("ROC calculation requires non-empty observed outcomes and predicted probabilities.", call. = TRUE)
  }

  if (length(observed) != length(predicted_prob)) {
    stop("Observed outcomes and predicted probabilities must have the same length.", call. = TRUE)
  }

  if (length(unique(stats::na.omit(observed))) < 2L) {
    stop("ROC calculation requires both outcome classes to be present.", call. = TRUE)
  }

  pROC::roc(
    response = observed,
    predictor = predicted_prob,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
}

plot_model_roc_curve <- function(prediction_df, model_label, dataset_label = "Holdout") {
  required_cols <- c("observed", "prob")

  if (!all(required_cols %in% names(prediction_df))) {
    stop(
      sprintf(
        "ROC plotting requires columns named %s.",
        paste(required_cols, collapse = " and ")
      ),
      call. = TRUE
    )
  }

  roc_obj <- build_binary_roc(prediction_df$observed, prediction_df$prob)
  auc_value <- as.numeric(pROC::auc(roc_obj))

  pROC::plot.roc(
    roc_obj,
    legacy.axes = TRUE,
    xlab = "False Positive Rate",
    ylab = "True Positive Rate",
    main = paste(model_label, "ROC Curve"),
    col = "#2C7FB8",
    lwd = 3,
    print.auc = FALSE,
    grid = TRUE,
    asp = 1
  )

  graphics::abline(a = 0, b = 1, lty = 2, lwd = 2, col = "#9A9A9A")
  graphics::legend(
    "bottomright",
    legend = sprintf("%s AUC = %.3f", dataset_label, auc_value),
    col = "#2C7FB8",
    lwd = 3,
    bty = "n"
  )

  invisible(roc_obj)
}

plot_feature_importance_ranking <- function(
  feature_importance_df,
  model_label,
  top_n = 20L,
  variable_type = "miRNA",
  selected_label = "Top M retained",
  unselected_label = "Ranked but not retained"
) {
  required_cols <- c(
    "variable",
    "variable_type",
    "importance_metric",
    "importance_value",
    "selected_for_final"
  )

  if (!all(required_cols %in% names(feature_importance_df))) {
    stop(
      sprintf(
        "Feature-importance plotting requires columns named %s.",
        paste(required_cols, collapse = ", ")
      ),
      call. = TRUE
    )
  }

  plot_df <- feature_importance_df[
    feature_importance_df$variable_type == variable_type,
    ,
    drop = FALSE
  ]

  if (!nrow(plot_df)) {
    graphics::plot.new()
    graphics::title(main = paste(model_label, "Feature Importance"))
    graphics::text(
      0.5,
      0.5,
      labels = paste("No", variable_type, "feature-importance values are available.")
    )

    return(invisible(plot_df))
  }

  top_n <- as.integer(top_n[[1]])

  if (!is.finite(top_n) || top_n < 1L) {
    stop("top_n must be a single positive integer.", call. = TRUE)
  }

  top_n <- min(top_n, nrow(plot_df))
  plot_df <- plot_df[seq_len(top_n), , drop = FALSE]
  metric_label <- unique(stats::na.omit(plot_df$importance_metric))

  if (!length(metric_label)) {
    metric_label <- "Importance"
  }

  bar_colors <- ifelse(plot_df$selected_for_final, "#2C7FB8", "#B7D4E8")
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  longest_label <- max(nchar(plot_df$variable), 0L)
  graphics::par(mar = c(5.1, max(8.1, 0.6 * longest_label), 4.1, 2.1))

  graphics::barplot(
    height = rev(plot_df$importance_value),
    names.arg = rev(plot_df$variable),
    horiz = TRUE,
    las = 1,
    col = rev(bar_colors),
    border = NA,
    main = paste(model_label, "Feature Importance"),
    xlab = metric_label[[1]]
  )

  if (any(plot_df$selected_for_final) && any(!plot_df$selected_for_final)) {
    graphics::legend(
      "bottomright",
      legend = c(selected_label, unselected_label),
      fill = c("#2C7FB8", "#B7D4E8"),
      border = NA,
      bty = "n"
    )
  }

  invisible(plot_df)
}

load_model_only_pipeline_data <- function(
  normalized_data_file,
  metadata_file,
  log2_pseudocount,
  train_fraction,
  run_limma_analysis,
  limma_rds_file,
  filter_stat_cols,
  fold_filter_top_n,
  limma_group_term,
  limma_adjustment_covariates
) {
  normalized_counts <- readxl::read_xlsx(normalized_data_file)
  normalized_counts$ID <- as.character(normalized_counts$ID)

  miRNA_value_columns <- grep("^MIMAT", names(normalized_counts), value = TRUE)
  normalized_counts[miRNA_value_columns] <- lapply(
    normalized_counts[miRNA_value_columns],
    function(column_values) log2(column_values + log2_pseudocount)
  )

  metadata <- readxl::read_xlsx(metadata_file)
  metadata$ID <- as.character(metadata$ID)
  metadata$CC <- factor(metadata$CC, levels = c("Controls", "PDAC case"))
  metadata$Sex <- factor(metadata$Sex, levels = c("Female", "Male"))
  metadata$Race <- factor(metadata$Race, levels = c("Asian", "Black", "White", "Others"))
  metadata$matched_set <- factor(metadata$matched_set)
  metadata$Age <- as.numeric(metadata$age_sam)
  metadata$CC_bin <- ifelse(metadata$CC == "PDAC case", 1, 0)

  full_dat <- dplyr::left_join(metadata, normalized_counts, by = "ID")
  miRNA_all_columns <- grep("^MIMAT", names(full_dat), value = TRUE)

  all_sets <- unique(full_dat$matched_set)
  n_train_sets <- floor(train_fraction * length(all_sets))
  train_set_ids <- sample(all_sets, n_train_sets)

  train_set <- full_dat[full_dat$matched_set %in% train_set_ids, , drop = FALSE]
  holdout <- full_dat[!(full_dat$matched_set %in% train_set_ids), , drop = FALSE]

  if (run_limma_analysis) {
    dir.create(dirname(limma_rds_file), showWarnings = FALSE, recursive = TRUE)

    train_filter_fn <- get("differential_expression_filter_fn", mode = "function")

    train_filter_results <- train_filter_fn(
      train_fold_dat = train_set,
      miRNA_cols = miRNA_all_columns,
      filter_stat_cols = filter_stat_cols,
      top_n = fold_filter_top_n,
      group_term = limma_group_term,
      adjustment_covariates = limma_adjustment_covariates
    )

    saveRDS(train_filter_results, limma_rds_file)
  } else {
    if (!file.exists(limma_rds_file)) {
      stop(
        paste(
          "The blocked limma cache was not found at",
          limma_rds_file,
          ". Set run_limma_analysis <- TRUE to rebuild it."
        ),
        call. = TRUE
      )
    }

    train_filter_results <- readRDS(limma_rds_file)
  }

  list(
    full_dat = full_dat,
    miRNA_all_columns = miRNA_all_columns,
    train_set = train_set,
    holdout = holdout,
    train_filter_results = train_filter_results
  )
}

map_to_training_levels <- function(values, training_levels) {
  if (!length(training_levels)) {
    stop("No training levels were available for factor alignment.", call. = TRUE)
  }

  mapped_values <- as.character(values)
  mapped_values[!(mapped_values %in% training_levels)] <- training_levels[[1]]
  factor(mapped_values, levels = training_levels)
}

build_penalty_factor <- function(column_names, protected_terms) {
  penalty_factor <- rep(1, length(column_names))

  for (term_name in protected_terms) {
    penalty_factor[startsWith(column_names, term_name)] <- 0
  }

  penalty_factor
}

prepare_fold_inputs <- function(
  train_fold_dat,
  valid_fold_dat,
  selected_miRNA,
  model_covariates,
  factor_model_covariates
) {
  if (!is.character(selected_miRNA)) {
    stop("selected_miRNA must be a character vector of feature names.", call. = TRUE)
  }

  selected_miRNA <- unique(selected_miRNA[!is.na(selected_miRNA) & nzchar(selected_miRNA)])

  train_fold_dat$CC_bin <- ifelse(train_fold_dat$CC == "PDAC case", 1, 0)
  valid_fold_dat$CC_bin <- ifelse(valid_fold_dat$CC == "PDAC case", 1, 0)

  if (!length(selected_miRNA)) {
    stop("Feature filtering returned zero miRNAs.", call. = TRUE)
  }

  keep_cols <- c("CC", "CC_bin", model_covariates, selected_miRNA)
  train_model_dat <- train_fold_dat |>
    dplyr::select(tidyselect::all_of(keep_cols)) |>
    droplevels()
  valid_model_dat <- valid_fold_dat |>
    dplyr::select(tidyselect::all_of(keep_cols))

  for (covariate_name in intersect(factor_model_covariates, model_covariates)) {
    valid_model_dat[[covariate_name]] <- map_to_training_levels(
      valid_model_dat[[covariate_name]],
      levels(train_model_dat[[covariate_name]])
    )
  }

  raw_model_formula <- stats::reformulate(c(model_covariates, selected_miRNA), response = "CC_bin")
  train_design_matrix <- stats::model.matrix(raw_model_formula, data = train_model_dat)
  valid_design_matrix <- stats::model.matrix(raw_model_formula, data = valid_model_dat)

  x_train <- train_design_matrix[, -1, drop = FALSE]
  x_valid <- valid_design_matrix[, -1, drop = FALSE]
  design_feature_cols <- colnames(x_train)
  base_feature_cols <- design_feature_cols[!(design_feature_cols %in% selected_miRNA)]
  model_formula <- stats::reformulate(design_feature_cols, response = "CC_bin")
  base_model_formula <- if (length(base_feature_cols)) {
    stats::reformulate(base_feature_cols, response = "CC_bin")
  } else {
    stats::as.formula("CC_bin ~ 1")
  }
  class_formula <- stats::reformulate(design_feature_cols, response = "CC")

  train_design_dat <- data.frame(
    CC = train_model_dat$CC,
    CC_bin = train_model_dat$CC_bin,
    as.data.frame(x_train, check.names = FALSE),
    check.names = FALSE
  )

  valid_design_dat <- data.frame(
    CC = valid_model_dat$CC,
    CC_bin = valid_model_dat$CC_bin,
    as.data.frame(x_valid, check.names = FALSE),
    check.names = FALSE
  )

  list(
    selected_miRNA = selected_miRNA,
    model_covariates = model_covariates,
    design_feature_cols = design_feature_cols,
    base_feature_cols = base_feature_cols,
    base_model_formula = base_model_formula,
    model_formula = model_formula,
    class_formula = class_formula,
    x_train = x_train,
    x_valid = x_valid,
    penalty_factor = build_penalty_factor(design_feature_cols, model_covariates),
    y_train = train_model_dat$CC_bin,
    train_model_dat = train_model_dat,
    valid_model_dat = valid_model_dat,
    train_design_dat = train_design_dat,
    valid_design_dat = valid_design_dat
  )
}

run_single_model_cv <- function(
  train_set,
  k,
  selected_miRNA,
  model_name,
  model_runner,
  model_covariates,
  factor_model_covariates
) {
  cv_set_ids <- unique(train_set$matched_set)
  fold_id <- sample(rep(seq_len(k), length.out = length(cv_set_ids)))

  fold_map <- data.frame(
    matched_set = cv_set_ids,
    fold = fold_id,
    row.names = NULL
  )

  cv_predictions <- vector("list", k)
  cv_tuning <- vector("list", k)
  cv_feature_counts <- integer(k)
  cv_preselected_feature_counts <- integer(k)
  cv_selected_variables <- vector("list", k)
  cv_selected_miRNA <- vector("list", k)
  cv_feature_importance <- vector("list", k)

  for (fold in seq_len(k)) {
    test_sets <- fold_map$matched_set[fold_map$fold == fold]

    cv_train <- train_set[!(train_set$matched_set %in% test_sets), , drop = FALSE]
    cv_test <- train_set[train_set$matched_set %in% test_sets, , drop = FALSE]

    fold_inputs <- prepare_fold_inputs(
      train_fold_dat = cv_train,
      valid_fold_dat = cv_test,
      selected_miRNA = selected_miRNA,
      model_covariates = model_covariates,
      factor_model_covariates = factor_model_covariates
    )

    model_out <- model_runner(fold_inputs)

    cv_predictions[[fold]] <- data.frame(
      ID = cv_test$ID,
      matched_set = cv_test$matched_set,
      observed = cv_test$CC_bin,
      prob = as.numeric(model_out$prob),
      fold = fold,
      row.names = NULL,
      check.names = FALSE
    )

    cv_tuning[[fold]] <- data.frame(
      fold = fold,
      model_out$tuning,
      row.names = NULL,
      check.names = FALSE
    )

    cv_feature_counts[[fold]] <- length(model_out$selected_miRNA)
    cv_preselected_feature_counts[[fold]] <- length(fold_inputs$selected_miRNA)
    cv_selected_variables[[fold]] <- model_out$selected_variables
    cv_selected_miRNA[[fold]] <- model_out$selected_miRNA
    cv_feature_importance[[fold]] <- model_out$feature_importance
  }

  cv_predictions_df <- dplyr::bind_rows(cv_predictions)
  cv_tuning_df <- dplyr::bind_rows(cv_tuning)
  roc_obj_cv <- build_binary_roc(cv_predictions_df$observed, cv_predictions_df$prob)

  cv_summary <- data.frame(
    model = model_name,
    n_predictions = nrow(cv_predictions_df),
    cv_auc = as.numeric(pROC::auc(roc_obj_cv)),
    cv_accuracy = mean(ifelse(cv_predictions_df$prob > 0.5, 1, 0) == cv_predictions_df$observed),
    row.names = NULL,
    check.names = FALSE
  )

  list(
    fold_map = fold_map,
    cv_predictions = cv_predictions_df,
    cv_tuning = cv_tuning_df,
    cv_feature_counts = cv_feature_counts,
    cv_preselected_feature_counts = cv_preselected_feature_counts,
    cv_selected_variables = cv_selected_variables,
    cv_selected_miRNA = cv_selected_miRNA,
    cv_feature_importance = cv_feature_importance,
    cv_summary = cv_summary
  )
}

run_single_model_holdout <- function(
  train_set,
  holdout,
  selected_miRNA,
  model_name,
  model_runner,
  model_covariates,
  factor_model_covariates
) {
  final_fold_inputs <- prepare_fold_inputs(
    train_fold_dat = train_set,
    valid_fold_dat = holdout,
    selected_miRNA = selected_miRNA,
    model_covariates = model_covariates,
    factor_model_covariates = factor_model_covariates
  )

  final_model_output <- model_runner(final_fold_inputs)

  holdout_predictions <- data.frame(
    model = model_name,
    observed = final_fold_inputs$valid_design_dat$CC_bin,
    prob = as.numeric(final_model_output$prob),
    row.names = NULL,
    check.names = FALSE
  )
  holdout_predictions$predicted_class <- ifelse(holdout_predictions$prob > 0.5, 1, 0)

  holdout_roc <- build_binary_roc(holdout_predictions$observed, holdout_predictions$prob)
  holdout_summary <- data.frame(
    model = model_name,
    test_accuracy = mean(holdout_predictions$predicted_class == holdout_predictions$observed),
    test_auc = as.numeric(pROC::auc(holdout_roc)),
    row.names = NULL,
    check.names = FALSE
  )

  list(
    final_fold_inputs = final_fold_inputs,
    final_model_output = final_model_output,
    holdout_roc = holdout_roc,
    holdout_predictions = holdout_predictions,
    holdout_summary = holdout_summary
  )
}