run_multimir_query <- function(mirna_id, table_name, summary = FALSE) {
  output_connection <- file(nullfile(), open = "wt")
  message_connection <- file(nullfile(), open = "wt")

  on.exit(
    {
      sink(type = "message")
      sink()
      close(message_connection)
      close(output_connection)
    },
    add = TRUE
  )

  sink(output_connection)
  sink(message_connection, type = "message")

  multiMiR::get_multimir(
    org = "hsa",
    mirna = mirna_id,
    table = table_name,
    summary = summary
  )
}

add_query_metadata <- function(query_row, data_frame) {
  if (!nrow(data_frame)) {
    return(data_frame)
  }

  data.frame(
    miRNA = query_row$miRNA[[1]],
    Annotation = query_row$Annotation[[1]],
    multiMiR_query_id = query_row$multiMiR_query_id[[1]],
    n_models_selected = query_row$n_models_selected[[1]],
    data_frame,
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

summary_column <- function(data_frame, column_name) {
  if (column_name %in% names(data_frame)) {
    data_frame[[column_name]]
  } else {
    rep(NA, nrow(data_frame))
  }
}

summarize_validated_targets <- function(query_row, validated_multimir_summary) {
  if (!nrow(validated_multimir_summary)) {
    return(data.frame(
      miRNA = character(0),
      Annotation = character(0),
      multiMiR_query_id = character(0),
      n_models_selected = numeric(0),
      target_symbol = character(0),
      target_entrez = character(0),
      target_ensembl = character(0),
      validated_support = numeric(0),
      total_support = numeric(0),
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    ))
  }

  validated_target_summary <- data.frame(
    miRNA = query_row$miRNA[[1]],
    Annotation = query_row$Annotation[[1]],
    multiMiR_query_id = query_row$multiMiR_query_id[[1]],
    n_models_selected = query_row$n_models_selected[[1]],
    target_symbol = summary_column(validated_multimir_summary, "target_symbol"),
    target_entrez = summary_column(validated_multimir_summary, "target_entrez"),
    target_ensembl = summary_column(validated_multimir_summary, "target_ensembl"),
    validated_support = summary_column(validated_multimir_summary, "validated.sum"),
    total_support = summary_column(validated_multimir_summary, "all.sum"),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  validated_target_summary <- validated_target_summary[
    !is.na(validated_target_summary$target_symbol) &
      validated_target_summary$target_symbol != "",
    ,
    drop = FALSE
  ]

  unique(validated_target_summary)
}

summarize_disease_drug <- function(query_row, disease_multimir_data) {
  if (!nrow(disease_multimir_data)) {
    return(data.frame())
  }

  disease_groups <- split(disease_multimir_data, disease_multimir_data$disease_drug)

  data.frame(
    miRNA = query_row$miRNA[[1]],
    Annotation = query_row$Annotation[[1]],
    multiMiR_query_id = query_row$multiMiR_query_id[[1]],
    n_models_selected = query_row$n_models_selected[[1]],
    disease_drug = names(disease_groups),
    n_records = vapply(disease_groups, nrow, integer(1)),
    n_databases = vapply(
      disease_groups,
      function(group_df) dplyr::n_distinct(group_df$database),
      integer(1)
    ),
    n_pubmed = vapply(
      disease_groups,
      function(group_df) dplyr::n_distinct(group_df$paper_pubmedID[group_df$paper_pubmedID != ""], na.rm = TRUE),
      integer(1)
    ),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

has_rows <- function(data_frame) {
  nrow(data_frame) > 0
}

normalize_disease_drug_data <- function(disease_multimir_data) {
  required_columns <- c("disease_drug", "database", "paper_pubmedID")

  for (column_name in required_columns) {
    if (!(column_name %in% names(disease_multimir_data))) {
      disease_multimir_data[[column_name]] <- rep(NA_character_, nrow(disease_multimir_data))
    }
  }

  disease_multimir_data
}

empty_multi_mir_query_details <- function() {
  data.frame(
    miRNA = character(0),
    Annotation = character(0),
    multiMiR_query_id = character(0),
    n_models_selected = numeric(0),
    row.names = NULL,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

normalize_multi_mir_query_id <- function(miRNA_id, org = "hsa") {
  query_id <- trimws(miRNA_id)

  if (is.na(query_id) || !nzchar(query_id)) {
    return(NA_character_)
  }

  if (grepl("^[[:alpha:]]{3}-(miR|mir|let)-", query_id)) {
    return(query_id)
  }

  if (grepl("^(miR|mir|let)-", query_id)) {
    return(paste0(org, "-", query_id))
  }

  query_id
}

build_multi_mir_query_details <- function(miRNA_input, n_models_selected = 1L, org = "hsa") {
  if (!length(miRNA_input)) {
    return(empty_multi_mir_query_details())
  }

  n_models_selected <- rep_len(n_models_selected, length(miRNA_input))

  expanded_inputs <- lapply(seq_along(miRNA_input), function(input_index) {
    input_value <- miRNA_input[[input_index]]

    if (is.na(input_value) || !nzchar(trimws(input_value))) {
      return(NULL)
    }

    input_tokens <- trimws(strsplit(input_value, ",", fixed = TRUE)[[1]])
    input_tokens <- input_tokens[nzchar(input_tokens)]

    if (!length(input_tokens)) {
      return(NULL)
    }

    query_ids <- vapply(
      input_tokens,
      normalize_multi_mir_query_id,
      FUN.VALUE = character(1),
      org = org
    )

    data.frame(
      miRNA = query_ids,
      Annotation = input_tokens,
      multiMiR_query_id = query_ids,
      n_models_selected = rep(n_models_selected[[input_index]], length(input_tokens)),
      row.names = NULL,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  })

  expanded_inputs <- Filter(Negate(is.null), expanded_inputs)

  if (!length(expanded_inputs)) {
    return(empty_multi_mir_query_details())
  }

  query_details <- dplyr::bind_rows(expanded_inputs)
  query_details <- query_details[
    !is.na(query_details$multiMiR_query_id) &
      nzchar(query_details$multiMiR_query_id),
    ,
    drop = FALSE
  ]
  query_details <- query_details[
    !duplicated(query_details$multiMiR_query_id),
    ,
    drop = FALSE
  ]
  row.names(query_details) <- NULL

  query_details
}

initialize_multi_mir_results <- function(consensus_miRNA_details, queried_miRNA) {
  list(
    consensus_miRNA = consensus_miRNA_details,
    queried_miRNA = queried_miRNA,
    skipped_query_count = nrow(consensus_miRNA_details) - nrow(queried_miRNA),
    validated = list(
      per_query = list(),
      raw_data = data.frame(),
      summary = data.frame(),
      total_record_count = 0L
    ),
    disease_drug = list(
      per_query = list(),
      raw_data = data.frame(),
      summary = data.frame(),
      total_record_count = 0L
    ),
    query_failures = data.frame(
      miRNA = character(0),
      multiMiR_query_id = character(0),
      table_name = character(0),
      error_message = character(0),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  )
}

run_multi_mir_queries <- function(queried_miRNA, consensus_miRNA_details = queried_miRNA) {
  consensus_miRNA_details <- consensus_miRNA_details[
    !duplicated(consensus_miRNA_details$miRNA),
    ,
    drop = FALSE
  ]

  queried_miRNA <- queried_miRNA[
    !duplicated(queried_miRNA$miRNA) &
      !is.na(queried_miRNA$multiMiR_query_id) &
      nzchar(queried_miRNA$multiMiR_query_id),
    ,
    drop = FALSE
  ]

  multimir_results <- initialize_multi_mir_results(
    consensus_miRNA_details = consensus_miRNA_details,
    queried_miRNA = queried_miRNA
  )

  if (nrow(queried_miRNA) == 0) {
    return(multimir_results)
  }

  validated_per_query <- vector("list", nrow(queried_miRNA))
  disease_per_query <- vector("list", nrow(queried_miRNA))
  validated_raw_data_list <- vector("list", nrow(queried_miRNA))
  validated_summary_list <- vector("list", nrow(queried_miRNA))
  disease_raw_data_list <- vector("list", nrow(queried_miRNA))
  disease_summary_list <- vector("list", nrow(queried_miRNA))
  query_failures <- list()

  names(validated_per_query) <- queried_miRNA$miRNA
  names(disease_per_query) <- queried_miRNA$miRNA

  for (query_index in seq_len(nrow(queried_miRNA))) {
    query_row <- queried_miRNA[query_index, , drop = FALSE]
    query_id <- query_row$multiMiR_query_id[[1]]

    validated_multimir_data <- data.frame()
    validated_multimir_summary <- data.frame()

    validated_multimir <- tryCatch(
      run_multimir_query(
        mirna_id = query_id,
        table_name = "validated",
        summary = TRUE
      ),
      error = function(err) {
        query_failures[[length(query_failures) + 1L]] <<- data.frame(
          miRNA = query_row$miRNA[[1]],
          multiMiR_query_id = query_id,
          table_name = "validated",
          error_message = conditionMessage(err),
          row.names = NULL,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

        NULL
      }
    )

    if (!is.null(validated_multimir)) {
      validated_multimir_data <- as.data.frame(validated_multimir@data)
      validated_multimir_summary <- as.data.frame(validated_multimir@summary)
    }

    disease_multimir_data <- data.frame()

    disease_multimir <- tryCatch(
      run_multimir_query(
        mirna_id = query_id,
        table_name = "disease.drug",
        summary = FALSE
      ),
      error = function(err) {
        query_failures[[length(query_failures) + 1L]] <<- data.frame(
          miRNA = query_row$miRNA[[1]],
          multiMiR_query_id = query_id,
          table_name = "disease.drug",
          error_message = conditionMessage(err),
          row.names = NULL,
          stringsAsFactors = FALSE,
          check.names = FALSE
        )

        NULL
      }
    )

    if (!is.null(disease_multimir)) {
      disease_multimir_data <- as.data.frame(disease_multimir@data)
      disease_multimir_data <- normalize_disease_drug_data(disease_multimir_data)
      disease_multimir_data <- disease_multimir_data[
        !is.na(disease_multimir_data$disease_drug) &
          disease_multimir_data$disease_drug != "",
        ,
        drop = FALSE
      ]
    }

    validated_per_query[[query_index]] <- list(
      query = query_row,
      data = validated_multimir_data,
      summary = validated_multimir_summary
    )
    disease_per_query[[query_index]] <- list(
      query = query_row,
      data = disease_multimir_data
    )

    validated_raw_data_list[[query_index]] <- add_query_metadata(query_row, validated_multimir_data)
    validated_summary_list[[query_index]] <- summarize_validated_targets(query_row, validated_multimir_summary)
    disease_raw_data_list[[query_index]] <- add_query_metadata(query_row, disease_multimir_data)
    disease_summary_list[[query_index]] <- summarize_disease_drug(query_row, disease_multimir_data)
  }

  validated_raw_data_list <- Filter(has_rows, validated_raw_data_list)
  validated_summary_list <- Filter(has_rows, validated_summary_list)
  disease_raw_data_list <- Filter(has_rows, disease_raw_data_list)
  disease_summary_list <- Filter(has_rows, disease_summary_list)

  multimir_results$validated$per_query <- validated_per_query
  multimir_results$disease_drug$per_query <- disease_per_query

  if (length(query_failures) > 0) {
    multimir_results$query_failures <- dplyr::bind_rows(query_failures)
  }

  if (length(validated_raw_data_list) > 0) {
    multimir_results$validated$raw_data <- dplyr::bind_rows(validated_raw_data_list)
    multimir_results$validated$total_record_count <- nrow(multimir_results$validated$raw_data)
  }

  if (length(validated_summary_list) > 0) {
    multimir_results$validated$summary <- dplyr::bind_rows(validated_summary_list)
    multimir_results$validated$summary <- multimir_results$validated$summary[
      order(
        -multimir_results$validated$summary$n_models_selected,
        -multimir_results$validated$summary$validated_support,
        -multimir_results$validated$summary$total_support,
        multimir_results$validated$summary$miRNA,
        multimir_results$validated$summary$target_symbol
      ),
      ,
      drop = FALSE
    ]
  }

  if (length(disease_raw_data_list) > 0) {
    multimir_results$disease_drug$raw_data <- dplyr::bind_rows(disease_raw_data_list)
    multimir_results$disease_drug$total_record_count <- nrow(multimir_results$disease_drug$raw_data)
  }

  if (length(disease_summary_list) > 0) {
    multimir_results$disease_drug$summary <- dplyr::bind_rows(disease_summary_list)
    multimir_results$disease_drug$summary <- multimir_results$disease_drug$summary[
      order(
        -multimir_results$disease_drug$summary$n_models_selected,
        -multimir_results$disease_drug$summary$n_records,
        -multimir_results$disease_drug$summary$n_databases,
        multimir_results$disease_drug$summary$miRNA,
        multimir_results$disease_drug$summary$disease_drug
      ),
      ,
      drop = FALSE
    ]
  }

  multimir_results
}

run_multi_mir_analysis <- function(miRNA_details) {
  run_multi_mir_queries(
    queried_miRNA = miRNA_details,
    consensus_miRNA_details = miRNA_details
  )
}