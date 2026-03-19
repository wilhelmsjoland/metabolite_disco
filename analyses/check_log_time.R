library(tidyverse)

pipeline_logs <- list.files(
  path = "V:aglycone_release_100um_24h/output",
  recursive = TRUE,
  pattern = "pipeline.log",
  full.names = TRUE
)

log_times <- purrr::map(
  .x = setNames(pipeline_logs, basename(pipeline_logs)),
  .f = ~ {
    lines <- readLines(.x)
    vec <- lines[
      grepl(
        "'\\w+' \\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}",
        lines,
        perl = TRUE
      )
    ]
    timestamps <- stringr::str_extract(
      string = vec,
      pattern = "\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}"
    )
    times <- as.POSIXct(timestamps, format = "%Y-%m-%d %H:%M:%S")
    names <- stringr::str_extract(vec, "'([^']+)'", group = 1)

    tibble::tibble(
      script = names,
      start = times,
      duration = round(difftime(dplyr::lead(times), times, units = "mins"), 2)
    ) %>%
      dplyr::add_row(
        script = "total",
        duration = sum(.$duration, na.rm = TRUE)
      )
  }
)

print(log_times)

dplyr::bind_rows(log_times, .id = "log_file") %>%
  dplyr::filter(grepl("annotation", script)) %>%
  dplyr::arrange(desc(duration))

# Annotation chromatogram times
for (i in pipeline_logs) {
  lines <- readLines(i)
  start <- grep("Generating chromatograms for annotated database", lines)
  end <- grep("'mz_predictions'", lines)
  range <- lines[start:end]
  anno_time <- grep("(\\d+)/\\1 \\(100%\\)", range, perl = TRUE, value = TRUE)
  print(anno_time)
}
