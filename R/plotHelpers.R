render_html_table <- function(data, caption = NULL, digits = NULL) {
  kable_args <- list(
    x = data,
    format = "html",
    row.names = FALSE,
    caption = caption
  )

  if (!is.null(digits)) {
    kable_args$digits <- digits
  }

  cat(
    do.call(knitr::kable, kable_args),
    sep = "\n"
  )
}

model_display_name <- function(model_name) {
  tools::toTitleCase(gsub("_", " ", model_name))
}

save_plot <- function (filename, plot, width = 7, height = 5, units = "in", dpi = 72) {
  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    width = width,
    height = height,
    units = units,
    dpi = dpi
  )
}