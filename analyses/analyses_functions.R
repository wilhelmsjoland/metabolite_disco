read_snakes <- function(
  path = NULL
) {
  list2env(
    readRDS(
      file.path(path)
    ),
    envir = .GlobalEnv
  )
}