run_fold_dt <- function(fold_inputs, cp_dt = 0.001, maxdepth = 20) {
  dt_fit <- rpart::rpart(
    fold_inputs$class_formula,
    data = fold_inputs$train_design_dat,
    method = "class",
    control = rpart::rpart.control(cp = cp_dt, maxdepth = maxdepth)
  )

  root_class_frequency <- table(fold_inputs$train_design_dat$CC)
  root_error_rate <- 1 - max(root_class_frequency) / sum(root_class_frequency)
  cv_classification_error <- data.frame(
    tree_size = as.integer(dt_fit$cptable[, "nsplit"] + 1L),
    nsplit = as.integer(dt_fit$cptable[, "nsplit"]),
    cp = as.numeric(dt_fit$cptable[, "CP"]),
    cv_relative_error = as.numeric(dt_fit$cptable[, "xerror"]),
    cv_classification_error_rate = as.numeric(dt_fit$cptable[, "xerror"]) * root_error_rate,
    cv_error_std = as.numeric(dt_fit$cptable[, "xstd"]),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  prob <- as.numeric(
    predict(dt_fit, newdata = fold_inputs$valid_design_dat, type = "prob")[, "PDAC case"]
  )

  split_variables <- dt_fit$frame$var
  split_variables <- split_variables[!is.na(split_variables) & split_variables != "<leaf>"]
  selected_variables <- unique(split_variables)
  selected_miRNA <- selected_variables[selected_variables %in% fold_inputs$selected_miRNA]

  split_variable_counts <- sort(table(split_variables), decreasing = TRUE)
  raw_importance <- dt_fit$variable.importance
  importance_variables <- unique(c(names(raw_importance), selected_variables))

  if (!length(importance_variables)) {
    feature_importance <- data.frame(
      importance_rank = integer(0),
      variable = character(0),
      variable_type = character(0),
      importance_metric = character(0),
      importance_value = numeric(0),
      primary_split_count = integer(0),
      selected_for_final = logical(0),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    importance_values <- raw_importance[importance_variables]
    importance_values[is.na(importance_values)] <- 0

    primary_split_count <- as.integer(split_variable_counts[importance_variables])
    primary_split_count[is.na(primary_split_count)] <- 0L

    feature_importance <- data.frame(
      variable = importance_variables,
      variable_type = ifelse(
        importance_variables %in% fold_inputs$selected_miRNA,
        "miRNA",
        "demographic"
      ),
      importance_metric = "TotalLossReduction",
      importance_value = as.numeric(importance_values),
      primary_split_count = primary_split_count,
      selected_for_final = importance_variables %in% selected_variables,
      row.names = NULL,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )

    feature_importance <- feature_importance[
      order(
        -feature_importance$importance_value,
        -feature_importance$primary_split_count,
        feature_importance$variable
      ),
      ,
      drop = FALSE
    ]
    rownames(feature_importance) <- NULL
    feature_importance$importance_rank <- seq_len(nrow(feature_importance))
    feature_importance <- feature_importance[
      ,
      c(
        "importance_rank",
        "variable",
        "variable_type",
        "importance_metric",
        "importance_value",
        "primary_split_count",
        "selected_for_final"
      ),
      drop = FALSE
    ]
  }

  list(
    prob = prob,
    selected_variables = selected_variables,
    selected_miRNA = selected_miRNA,
    cv_classification_error = cv_classification_error,
    feature_importance = feature_importance,
    tuning = data.frame(
      model = "decision_tree",
      n_selected_miRNA = length(selected_miRNA),
      cp = cp_dt,
      maxdepth = maxdepth,
      n_terms = length(selected_variables),
      n_nonzero_importance = nrow(feature_importance),
      stringsAsFactors = FALSE
    )
  )
}