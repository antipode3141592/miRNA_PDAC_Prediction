write_comma_separated_list <- function(items) {
  if (length(items) == 0) {
    return("")
  } else if (length(items) == 1) {
    return(items)
  } else {
    return(paste(
      paste(items[-length(items)], collapse = ", "),
      items[length(items)],
      sep = ifelse(length(items) > 2, ", and ", " and ")
    ))
  }
}

first_non_missing <- function(values) {
  values <- values[!is.na(values) & values != ""]

  if (length(values)) {
    values[[1]]
  } else {
    NA_character_
  }
}

render_basic_table_styles <- function() {
  cat(
    paste(
      "<style>",
      "table {",
      "  border-collapse: collapse;",
      "}",
      "table th,",
      "table td {",
      "  border: 1px solid #cfcfcf;",
      "  padding: 4px 8px;",
      "}",
      "</style>",
      sep = "\n"
    )
  )
}

format_elapsed_time <- function(elapsed_seconds) {
  total_seconds <- max(0L, as.integer(round(elapsed_seconds)))
  hours <- total_seconds %/% 3600L
  minutes <- (total_seconds %% 3600L) %/% 60L
  seconds <- total_seconds %% 60L

  sprintf("%02d:%02d:%02d", hours, minutes, seconds)
}

render_report_runtime <- function(report_start_time) {
  if (missing(report_start_time) || is.null(report_start_time) || !inherits(report_start_time, "POSIXt")) {
    cat("<p>Elapsed runtime was not available because the start timestamp was not initialized.</p>\n")
    return(invisible(NULL))
  }

  report_end_time <- Sys.time()
  elapsed_seconds <- as.numeric(difftime(report_end_time, report_start_time, units = "secs"))

  cat(
    paste0(
      "<p>Report runtime: ",
      format_elapsed_time(elapsed_seconds),
      " (started ",
      format(report_start_time, "%Y-%m-%d %H:%M:%S"),
      "; finished ",
      format(report_end_time, "%Y-%m-%d %H:%M:%S"),
      ").</p>\n"
    )
  )

  invisible(NULL)
}

build_binary_confusion_summary <- function(
  observed,
  predicted_prob,
  threshold = 0.5
) {
  if (length(observed) != length(predicted_prob)) {
    stop("Observed outcomes and predicted probabilities must have the same length.", call. = TRUE)
  }

  predicted_class <- ifelse(predicted_prob > threshold, 1, 0)
  class_levels <- c("Controls", "PDAC case")

  confusion_table <- table(
    Prediction = factor(predicted_class, levels = c(0, 1), labels = class_levels),
    Reference = factor(observed, levels = c(0, 1), labels = class_levels)
  )

  true_negative <- sum(predicted_class == 0 & observed == 0, na.rm = TRUE)
  true_positive <- sum(predicted_class == 1 & observed == 1, na.rm = TRUE)
  false_positive <- sum(predicted_class == 1 & observed == 0, na.rm = TRUE)
  false_negative <- sum(predicted_class == 0 & observed == 1, na.rm = TRUE)
  total_cases <- sum(!is.na(observed) & !is.na(predicted_class))

  precision <- if ((true_positive + false_positive) > 0) {
    true_positive / (true_positive + false_positive)
  } else {
    NA_real_
  }

  recall <- if ((true_positive + false_negative) > 0) {
    true_positive / (true_positive + false_negative)
  } else {
    NA_real_
  }

  specificity <- if ((true_negative + false_positive) > 0) {
    true_negative / (true_negative + false_positive)
  } else {
    NA_real_
  }

  accuracy <- if (total_cases > 0) {
    (true_positive + true_negative) / total_cases
  } else {
    NA_real_
  }

  f1_score <- if (is.finite(precision) && is.finite(recall) && (precision + recall) > 0) {
    2 * precision * recall / (precision + recall)
  } else {
    NA_real_
  }

  list(
    table = confusion_table,
    byClass = c(
      Precision = precision,
      Recall = recall,
      Sensitivity = recall,
      Specificity = specificity,
      Accuracy = accuracy,
      F1 = f1_score
    ),
    threshold = threshold
  )
}