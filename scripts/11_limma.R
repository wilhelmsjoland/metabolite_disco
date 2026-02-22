cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Linear models using limma ----------------------------------------------------
# ==============================================================================
cli::cli_h3("Running linear models using limma")

group_used <- factor(meta$group)
design <- model.matrix(~ 0 + group_used)
colnames(design) <- levels(group_used)

comparisons <- combn(
  x = levels(group_used),
  m = 2,
  simplify = TRUE) %>%
  t(.) %>%
  as.data.frame(.) %>%
  tibble::as_tibble(
    x = .,
    rownames = "rownumber",
    .name_repair = "universal_quiet"
  ) %>%
  dplyr::mutate(comp = paste0(V1, "-", V2)) %>%
  dplyr::pull(comp)

contrasts_mat <- limma::makeContrasts(
  contrasts = comparisons,
  levels = design
)

fit <- limma::lmFit(
  log2(SummarizedExperiment::assay(res, "norm")), # don't use imputed
  design = design
)
fit <- limma::contrasts.fit(fit, contrasts_mat)
fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

limma_res <- list()
for (i in comparisons) {
  tmp <- limma::topTable(
    fit = fit,
    coef = i,
    number = Inf,
    adjust.method = "BH",
    sort.by = "none"
  ) %>%
    tibble::as_tibble(., rownames = "feature") %>%
    dplyr::mutate(contrast = i) %>%
    dplyr::left_join(
      x = .,
      y = SummarizedExperiment::rowData(res) %>%
        tibble::as_tibble(., rownames = "feature"),
      by = "feature"
    )
  limma_res[[i]] <- tmp
}

# ==============================================================================
# Saving limma results to tables -----------------------------------------------
# ==============================================================================
full_limma <- tibble::tibble()
for (i in names(limma_res)) {
  tmp_tib <- limma_res[[i]]
  full_limma <- dplyr::bind_rows(full_limma, tmp_tib)
}

message("Saving intensity information to tables...")
# Creating tables of all output data and saving to tables
# TODO Fix this dumb logic here or use as witch statement?
assay_names <- names(SummarizedExperiment::assays(res))
for (i in assay_names) {
  if (i %in% c("norm", "norm_filled")) {
    full_data <- dplyr::left_join(
      x = SummarizedExperiment::rowData(res) %>%
        tibble::as_tibble(., rownames = "feature"),
      y = SummarizedExperiment::assay(res, i) %>%
        log2() %>%
        t() %>%
        scale(., center = TRUE, scale = TRUE) %>%
        t() %>%
        tibble::as_tibble(., rownames = "feature"),
      by = "feature"
    )
  } else if (i %in% c("raw", "raw_filled")) {
    full_data <- dplyr::left_join(
      x = SummarizedExperiment::rowData(res) %>%
        tibble::as_tibble(., rownames = "feature"),
      y = SummarizedExperiment::assay(res, i) %>%
        tibble::as_tibble(., rownames = "feature"),
      by = "feature"
    )
  }

  assign(
    x = paste0("full_", i),
    value = full_data,
    envir = .GlobalEnv
  )

  readr::write_csv(
    x = full_data,
    file = paste0(opt$output, "/tables/full_", i, ".csv"),
    na = "NA",
    col_names = TRUE,
    append = FALSE
  )
}