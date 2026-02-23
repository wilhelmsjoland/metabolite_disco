cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# For each tibble create 1. untransformed, 2. log2, 3. log2 & scaled -----------
# ==============================================================================
all_names <- setNames(assay_names, assay_names)
intensities_mat <- purrr::map(
  .x = all_names,
  .f = ~ {
    mat <- SummarizedExperiment::assay(res, .x)

    list(
      untransformed = mat,
      log2 = log2(mat),
      log2_scale = mat %>%
        log2() %>%
        t() %>%
        scale(center = TRUE, scale = TRUE) %>%
        t()
    )
  }
)
cli::cli_alert_success(
  paste0( # just took a random one for 'intensities_mat$norm' (shouldn't matter)
    "Generated dataframes of {.val {names(intensities_mat$norm)}} for: ",
    "{.val {names(intensities_mat)}}"
  )
)

# ==============================================================================
# Run fits on all normalized assays --------------------------------------------
# ==============================================================================
cli::cli_h3("Running linear models using limma")
group_used <- factor(meta$group)
design <- model.matrix(~ 0 + group_used)
colnames(design) <- levels(group_used)

comparisons <- combn(
  x = levels(group_used),
  m = 2,
  simplify = TRUE
) %>%
  t() %>%
  as.data.frame() %>%
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

norm_names <- assay_names[grepl("norm", assay_names)]
norm_names <- setNames(norm_names, norm_names)

# Do not test raws since they are not median-scaled
# which is used for removing technical variation
limma_fits <- purrr::map(
  .x = norm_names,
  .f = ~ {
    fit <- limma::lmFit(
      object = intensities_mat[[.x]][["log2"]],
      design = design
    )
    fit <- limma::contrasts.fit(fit, contrasts_mat)
    fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
  }
)
cli::cli_alert_success(
  paste0(
    "Ran linear models with limma on: ",
    "{.val {norm_names}}"
  )
)

# ==============================================================================
# Extract all fits and comparisons ---------------------------------------------
# ==============================================================================
full_limmas <- list()
for (i in norm_names) {
  tmp_all_comp <- tibble::tibble()
  for (j in comparisons) {
    tmp <- limma::topTable(
      fit = limma_fits[[i]],
      coef = j,
      number = Inf,
      adjust.method = "BH",
      sort.by = "none"
    ) %>%
      tibble::as_tibble(., rownames = "feature") %>%
      dplyr::mutate(contrast = j)

    tmp_all_comp <- dplyr::bind_rows(tmp_all_comp, tmp)
  }

  if (nrow(tmp_all_comp) != (nrow(tmp) * length(comparisons))) {
    cli::cli_abort("The lengths of limma tables are mismatched")
  } else {
    full_limmas[[i]] <- tmp_all_comp
  }
}
cli::cli_bullets(
  c(
    "v" = paste0(
      "Extracted all fits for: {.val {norm_names}} for comparisons: "
    ),
    setNames(comparisons, rep("i", length(comparisons)))
  )
)

# ==============================================================================
# Map all feature definitions to intensity tibbles -----------------------------
# ==============================================================================
# Don't map the limmas to the definitions as it isn't needed
res_defs <- tibble::as_tibble(
  x = SummarizedExperiment::rowData(res),
  rownames = "feature"
)

intensities <- purrr::modify_depth(
  .x = intensities_mat,
  .depth = 2,
  .f = ~ {
    res_defs %>%
      dplyr::left_join(
        y = tibble::as_tibble(.x, rownames = "feature"),
        by = "feature"
      )
  }
)

# ==============================================================================
# Saving intensity information to .csv tables ----------------------------------
# ==============================================================================
for (i in names(intensities)) {
  for (j in names(intensities[[i]])) {
    tmp_file_path <- paste0(opt$output, "/tables/", i, "_", j, ".csv")
    if (file.exists(tmp_file_path)) {
      cli::cli_alert_danger(
        paste0(
          "{.path {tmp_file_path} already exists. Not overwriting.}"
        )
      )
    } else {
      readr::write_csv(
        x = intensities[[i]][[j]],
        file = tmp_file_path,
        na = "NA",
        col_names = TRUE,
        append = FALSE
      )
      cli::cli_alert_success(
        paste0(
          "Saved {.val {intensities[[i]][[j]]}} to: ",
          "{.path {tmp_file_path}}"
        )
      )
    }
  }
}

cli::cli_alert_success(
  paste0( # just took a random one for 'intensities_mat$norm' (shouldn't matter)
    "Dataframes of {.val {names(intensities_mat$norm)}} for: ",
    "{.val {names(intensities_mat)}} saved in ",
    "{.path {file.path(opt$output, 'tables')}}"
  )
)

# ==============================================================================
# Saving linear model information to .csv tables -------------------------------
# ==============================================================================
for (i in names(full_limmas)) {
  tmp_file_path <- paste0(opt$output, "/tables/limma_", i, ".csv")
  if (file.exists(tmp_file_path)) {
    cli::cli_alert_danger(
      paste0(
        "{.path {tmp_file_path} already exists. Not overwriting.}"
      )
    )
  } else {
    readr::write_csv(
      x = full_limmas[[i]],
      file = paste0(opt$output, "/tables/limma_", i, ".csv"),
      na = "NA",
      col_names = TRUE,
      append = FALSE
    )

    cli::cli_alert_success(
      paste0(
        "Saved {.val {full_limmas[[i]]}} to: ",
        "{.path {tmp_file_path}}"
      )
    )
  }
}

cli::cli_alert_success(
  paste0(
    "Dataframes of {.val {names(full_limmas)}} saved in: ",
    "{.path {file.path(opt$output, 'tables')}}"
  )
)
