# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
suppressWarnings(
  suppressPackageStartupMessages({
    library(tidyverse)
    library(xcms)
    library(mzR)
    library(Spectra)
    library(MsExperiment)
    library(optparse)
    library(BiocParallel)
    library(MetaboCoreUtils)
    library(Rdisop)
    library(cli)
  })
)
source("analyses/standard_disco_args.R")
source("scripts/functions.R")
bp <- BiocParallel::MulticoreParam(workers = opt$cores, fallback = FALSE)
BiocParallel::bpstart(bp)

dir.create(file.path(opt$output, "objects"), FALSE, TRUE)
dir.create(file.path(opt$output, "tables"), FALSE, TRUE)
dir.create(file.path(opt$output, "peaks", "facet"), FALSE, TRUE)
dir.create(file.path(opt$output, "peaks", "single"), FALSE, TRUE)
dir.create(file.path(opt$output, "internal_standard"), FALSE, TRUE)

cli::cli_progress_step("Import metadata")
stds <- import_mzml(data_path = opt$data_path, meta_file = opt$meta_file)
cli::cli_progress_done()

cli::cli_progress_step("Import experiment")
stds_exp_path <- file.path(opt$output, "objects", "stds_exp.rds")
if (file.exists(stds_exp_path)) {
  stds_exp <- readRDS(stds_exp_path)
} else {
  stds_exp <- MsExperiment::readMsExperiment(
    spectraFiles = stds$path,
    sampleData = stds,
    backend = Spectra::MsBackendMemory(),
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1L
  )
  saveRDS(stds_exp, stds_exp_path)
}
cli::cli_progress_done()

# ==============================================================================
# Create wide internal standard chromatograms ------------------------------–---
# ==============================================================================
std_mz_theory <- get_theory_mz(
  chem_form = opt$is_formula,
  adduct = opt$is_adduct
)
std_mz_range <- get_short_mz_range(std_mz_theory, mz_window = 0.02)

std_chr_path <- file.path(opt$output, "objects", "std_chr.rds")
cli::cli_progress_step("Generate wide internal standard chromatogram")
if (file.exists(std_chr_path)) {
  std_chr <- readRDS(std_chr_path)
} else {
  std_chr <- xcms::chromatogram(
    BPPARAM = bp,
    chunkSize = opt$cores,
    object = stds_exp,
    mz = std_mz_range,
    aggregationFun = "sum"
  )
  saveRDS(std_chr, std_chr_path)
}
cli::cli_progress_done()


# ==============================================================================
# Plot all wide internal standard chromatograms --------------------------------
# ==============================================================================
std_chr_vals <- purrr::map(
  .x = seq_len(ncol(std_chr)),
  .f = ~ {
    tibble::tibble(
      sample = colnames(std_chr)[.x],
      rt = std_chr[1, .x]@rtime,
      intensity = std_chr[1, .x]@intensity,
      mzmin = std_chr[1, .x]@mz[1],
      mzmax = std_chr[1, .x]@mz[2]
    )
  }
) %>%
  purrr::list_rbind() %>%
  dplyr::mutate(intensity = tidyr::replace_na(intensity, replace = 0))

cli::cli_progress_step("Plot internal standard chromatograms")
purrr::walk(
  .x = unique(std_chr_vals$sample),
  .f = ~ {
    tmp_p <- std_chr_vals %>%
      dplyr::filter(sample == .x) %>%
      ggplot(
        ggplot2::aes(
          x = rt,
          y = intensity,
          group = sample
        )
      ) +
      ggplot2::geom_line() +
      ggplot2::theme_classic() +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0, 0.05))) +
      ggplot2::labs(title = .x)

    ggplot2::ggsave(
      filename = file.path(
        opt$output,
        "internal_standard",
        paste0("is_", .x, ".pdf")
      ),
      plot = tmp_p,
      device = "pdf",
      height = 6,
      width = 6,
      units = "in"
    )
  }
)
cli::cli_progress_done()

# ==============================================================================
# Define expected ions from each standard formula ------------------------------
# ==============================================================================
cli::cli_progress_step("Calculate expected ions")
sample_map <- tibble::as_tibble(stds, rownames = "sample") %>%
  dplyr::select(-path) %>%
  tidyr::separate_rows(formula, group, sep = ", ") %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    mass = Rdisop::getMass(Rdisop::getMolecule(formula))
  ) %>%
  # dplyr::select(-c("unit", "bacteria")) %>%
  dplyr::distinct(group, .keep_all = TRUE)

expected_ions <- MetaboCoreUtils::mass2mz(
  x = setNames(sample_map$mass, sample_map$group),
  adduct = MetaboCoreUtils::adducts(polarity = opt$polarity)
) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("group") %>%
  tidyr::pivot_longer(
    cols = -group,
    names_to = "adduct",
    values_to = "theory_mz"
  ) %>%
  tibble::as_tibble()
cli::cli_progress_done()

# ==============================================================================
# Extract XICS -----------------------------------------------------------------
# ==============================================================================
cli::cli_progress_step("Extract XIC")
peak_chrs_vals <- NULL
peak_chrs_list <- list()
for (g in unique(expected_ions$group)) {
  keep <- stringr::str_detect(
    MsExperiment::sampleData(stds_exp)$group,
    stringr::fixed(g)
  )

  stds_exp_sub <- stds_exp[keep]

  exp_ions_sub <- expected_ions %>%
    dplyr::filter(group == g)

  peak_chrs <- NULL

  for (i in seq_len(nrow(exp_ions_sub))) {
    peak_chr <- tryCatch(
      xcms::chromatogram(
        BPPARAM = bp,
        chunkSize = opt$cores,
        object = stds_exp_sub,
        mz = exp_ions_sub$theory_mz[i] + c(-0.05, 0.05),
        aggregationFun = "sum",
        msLevel = 1L
      ),
      error = function(e) {
        cli::cli_warn(
          paste0(
            "Skipping adduct {exp_ions_sub$adduct[i]} for group ",
            "{g}: {conditionMessage(e)}"
          )
        )
        NULL
      }
    )

    if (is.null(peak_chr)) next

    rownames(peak_chr) <- exp_ions_sub$adduct[i]

    peak_chrs <- if (is.null(peak_chrs)) peak_chr else c(peak_chrs, peak_chr)
  }

  peak_chrs_val <- purrr::map(
    seq_len(nrow(peak_chrs)),
    function(i) {
      purrr::map(
        seq_len(ncol(peak_chrs)),
        function(j) {
          tibble::tibble(
            group = g,
            adduct = rownames(peak_chrs)[i],
            sample = colnames(peak_chrs)[j],
            rt = peak_chrs[i, j]@rtime,
            intensity = tidyr::replace_na(peak_chrs[i, j]@intensity, 0),
            mzmin = peak_chrs[i, j]@mz[1],
            mzmax = peak_chrs[i, j]@mz[2]
          )
        }
      ) %>%
        purrr::list_rbind()
    }
  ) %>%
    purrr::list_rbind()

  peak_chrs_vals <- dplyr::bind_rows(peak_chrs_vals, peak_chrs_val)
  peak_chrs_list[[g]] <- peak_chrs
}
cli::cli_progress_done()

readr::write_csv(
  x = peak_chrs_vals,
  file = file.path(opt$output, "tables", "chromatogram_values.csv")
)

# ==============================================================================
# Call peaks on extracted XICS -------------------------------------------------
# ==============================================================================
cli::cli_progress_step("Call peaks")
pk_tbl <- purrr::map(
  .x = names(peak_chrs_list),
  .f = function(g) {
    tot <- peak_chrs_list[[g]]

    called <- xcms::findChromPeaks(
      object = tot,
      BPPARAM = bp,
      chunkSize = opt$cores,
      param = xcms::CentWaveParam(
        ppm = opt$ppm_global,
        peakwidth = c(1, 100), # add
        snthresh = 10, # add
        prefilter = c(4, 5000), # add
        integrate = 2, # add
        noise = 1000 # add
      )
    )

    xcms::chromPeaks(called) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(
        group  = g,
        adduct = rownames(tot)[row],
        sample = colnames(tot)[column],
        .before = 1
      )
  }
) %>%
  purrr::list_rbind()
cli::cli_progress_done()

readr::write_csv(
  x = pk_tbl,
  file = file.path(opt$output, "tables", "full_peak_table.csv")
)

# ==============================================================================
# Plot faceted adducts chromatograms per standard ------------------------------
# ==============================================================================
cli::cli_progress_step("Plot faceted chromatograms")
purrr::walk(
  .x = unique(peak_chrs_vals$group),
  .f = ~ {
    tmp_p <- peak_chrs_vals %>%
      dplyr::filter(group == .x) %>%
      ggplot(
        ggplot2::aes(
          x = rt,
          y = intensity,
          group = sample
        )
      ) +
      ggplot2::geom_line() +
      ggplot2::theme_classic() +
      ggplot2::theme(
        strip.background = ggplot2::element_blank()
      ) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0, 0.05))) +
      ggplot2::facet_wrap(~adduct, scales = "free_x") +
      ggplot2::labs(
        title = .x,
        x = "Retention time (s)",
        y = "Intensity"
      )

    ggplot2::ggsave(
      filename = file.path(
        opt$output,
        "peaks",
        "facet",
        paste0(.x, ".pdf")
      ),
      plot = tmp_p,
      device = "pdf",
      height = 6,
      width = 8,
      units = "in"
    )
  }
)
cli::cli_progress_done()

# ==============================================================================
# Plot single adducts chromatograms per standard -------------------------------
# ==============================================================================
cli::cli_progress_step("Plot single chromatograms")
purrr::walk(
  .x = unique(peak_chrs_vals$group),
  .f = ~ {
    g <- .x
    purrr::walk(
      .x = unique(peak_chrs_vals$adduct[peak_chrs_vals$group == g]),
      .f = ~ {
        tmp_p <- peak_chrs_vals %>%
          dplyr::filter(group == g, adduct == .x) %>%
          ggplot(
            ggplot2::aes(
              x = rt,
              y = intensity,
              group = sample
            )
          ) +
          ggplot2::geom_line() +
          ggplot2::theme_classic() +
          ggplot2::theme(
            strip.background = ggplot2::element_blank()
          ) +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0, 0.05))) +
          ggplot2::labs(
            title = paste(g, .x),
            x = "Retention time (s)",
            y = "Intensity"
          )

        ggplot2::ggsave(
          filename = file.path(
            opt$output,
            "peaks",
            "single",
            paste0(g, "_", .x, ".pdf")
          ),
          plot = tmp_p,
          device = "pdf",
          height = 6,
          width = 8,
          units = "in"
        )
      }
    )
  }
)
cli::cli_progress_done()

# ==============================================================================
# Write peak summary tables ----------------------------------------------------
# ==============================================================================
cli::cli_progress_step("Write output tables")
max_standard_adducts <- pk_tbl %>%
  dplyr::filter(maxo > 10000) %>%
  dplyr::group_by(group, sample) %>%
  dplyr::slice_max(into, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, dplyr::desc(into))

all_standard_adducts <- pk_tbl %>%
  dplyr::filter(maxo > 10000) %>%
  dplyr::group_by(group, adduct, sample) %>%
  dplyr::slice_max(into, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(group, dplyr::desc(into))

pk_summary <- pk_tbl %>%
  dplyr::filter(maxo > 10000) %>%
  dplyr::group_by(group, sample) %>%
  dplyr::slice_max(into, n = 1) %>%
  dplyr::ungroup() %>%
  dplyr::group_by(group, adduct) %>%
  dplyr::summarize(
    median_into = median(into),
    median_maxo = median(maxo),
    rt = median(rt),
    mz = median(mz),
    rtmin = median(rtmin),
    rtmax = median(rtmax),
    mzmin = median(mzmin),
    mzmax = median(mzmax),
    .groups = "drop"
  ) %>%
  dplyr::group_by(group) %>%
  dplyr::slice_max(median_into, n = 1) %>%
  dplyr::ungroup()

readr::write_csv(
  x = max_standard_adducts,
  file = file.path(opt$output, "tables", "max_standard_adducts.csv")
)

readr::write_csv(
  x = all_standard_adducts,
  file = file.path(opt$output, "tables", "all_standard_adducts.csv")
)

readr::write_csv(
  x = pk_summary,
  file = file.path(opt$output, "tables", "top_standard_summary.csv")
)
cli::cli_progress_done()

cli::cli_alert_success("Pipeline finished.\n")
