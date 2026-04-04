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
source("scripts/functions.R")
source("analyses/standard_disco_functions.R")
source("analyses/standard_disco_args2.R")

bp_param <- if (opt$cores > 1) {
  BiocParallel::MulticoreParam(opt$cores)
} else {
  BiocParallel::SerialParam()
}

dir.create(file.path(opt$output, "objects"), FALSE, TRUE)
dir.create(file.path(opt$output, "tables"), FALSE, TRUE)
dir.create(file.path(opt$output, "peaks"), FALSE, TRUE)
dir.create(file.path(opt$output, "internal_standard_wide"), FALSE, TRUE)
dir.create(file.path(opt$output, "internal_standard"), FALSE, TRUE)

cli::cli_progress_step("Import metadata")
stds <- import_mzml(data_path = opt$data_path, meta_file = opt$meta_file)
cli::cli_progress_done()

cli::cli_progress_step("Import experiment")
stds_exp_path <- file.path(opt$output, "objects", "stds_exp.rds")
if (file.exists(stds_exp_path)) {
  stds_exp <- readRDS(stds_exp_path)
} else {
  stds_exp <- import_msexperiment_memory(
    stds,
    failure_log_path = file.path(
      opt$output,
      "tables",
      "import_skipped_nonempty_scans.csv"
    )
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

std_chr_path <- file.path(opt$output, "objects", "std_chr_wide.rds")
cli::cli_progress_step("Generate wide internal standard chromatogram")
if (file.exists(std_chr_path)) {
  std_chr <- readRDS(std_chr_path)
} else {
  std_chr <- xcms::chromatogram(
    object = stds_exp,
    BPPARAM = bp_param,
    chunkSize = opt$cores,
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
        "internal_standard_wide",
        paste0("wide_is_", .x, ".pdf")
      ),
      plot = tmp_p,
      device = "pdf",
      height = 6,
      width = 6,
      units = "in"
    )
  }
)

# ==============================================================================
# Create narrow internal standard chromatograms ------------------------------–-
# ==============================================================================
std_ranges <- get_rt_mz_range(chromatogram = std_chr, rt_window = 0.05)
std_eic_narrow_path <- file.path(opt$output, "objects", "std_eic_narrow.rds")
cli::cli_progress_step("Generate narrow internal standard chromatograms")
if (file.exists(std_eic_narrow_path)) {
  std_eic_narrow <- readRDS(std_eic_narrow_path)
} else {
  std_eic_narrow <- xcms::chromatogram(
    object = stds_exp,
    BPPARAM = bp_param,
    chunkSize = opt$cores,
    mz = std_ranges$mz_range,
    rt = std_ranges$rt_range,
    aggregationFun = "sum"
  )
  saveRDS(
    std_eic_narrow,
    std_eic_narrow_path
  )
}
cli::cli_progress_done()

std_eic_narrow_vals <- purrr::map(
  .x = seq_len(ncol(std_eic_narrow)),
  .f = ~ {
    tibble::tibble(
      sample = colnames(std_eic_narrow)[.x],
      rt = std_eic_narrow[1, .x]@rtime,
      intensity = std_eic_narrow[1, .x]@intensity,
      mzmin = std_eic_narrow[1, .x]@mz[1],
      mzmax = std_eic_narrow[1, .x]@mz[2]
    )
  }
) %>%
  purrr::list_rbind() %>%
  dplyr::mutate(intensity = tidyr::replace_na(intensity, replace = 0))

purrr::walk(
  .x = unique(std_eic_narrow_vals$sample),
  .f = ~ {
    tmp_p <- std_eic_narrow_vals %>%
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

# ==============================================================================
# Targeted standard search -----------------------------------------------------
# ==============================================================================
sample_targets <- tibble::as_tibble(stds, rownames = "filename") %>%
  tidyr::separate_rows(formula, group, sep = ", ") %>%
  dplyr::select(filename, group, formula)

targeted_features_path <- file.path(
  opt$output,
  "objects",
  "standard_feature_candidates.rds"
)

cli::cli_progress_step("Find targeted standard features")
if (file.exists(targeted_features_path)) {
  standards_all_features <- readRDS(targeted_features_path)
} else {
  standards_all_features <- build_targeted_feature_table(
    experiment = stds_exp,
    sample_targets = sample_targets,
    ppm = opt$ppm_global,
    BPPARAM = bp_param,
    chunkSize = opt$cores
  )
  saveRDS(standards_all_features, targeted_features_path)
}
cli::cli_progress_done()

standards_best_features <- standards_all_features %>%
  dplyr::group_by(group) %>%
  dplyr::arrange(
    dplyr::desc(feature_found),
    dplyr::desc(n_pos_signal_apex),
    adduct_rank,
    dplyr::desc(n_pos_signal_area),
    dplyr::desc(n_pos_samples_detected),
    dplyr::desc(enrichment_apex),
    dplyr::desc(enrichment_area),
    dplyr::desc(median_peak_area),
    dplyr::desc(median_peak_apex),
    peak_rt_mad,
    .by_group = TRUE
  ) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

cli::cli_progress_step("Plot standard chromatograms")
purrr::walk(
  .x = seq_len(nrow(standards_best_features)),
  .f = ~ {
    tmp_group <- standards_best_features$group[.x]
    tmp_adduct <- standards_best_features$adduct[.x]
    tmp_mz_theory <- standards_best_features$mz_theory[.x]
    tmp_feature_found <- standards_best_features$feature_found[.x]
    tmp_rt <- standards_best_features$rtmed[.x]
    tmp_rtmin <- standards_best_features$rtmin[.x]
    tmp_rtmax <- standards_best_features$rtmax[.x]
    tmp_mzmin <- standards_best_features$mzmin[.x]
    tmp_mzmax <- standards_best_features$mzmax[.x]
    tmp_samples <- sample_targets %>%
      dplyr::filter(group == tmp_group) %>%
      dplyr::distinct(filename) %>%
      dplyr::pull(filename)

    if (isTRUE(tmp_feature_found) && is.finite(tmp_mzmin) && is.finite(tmp_mzmax)) {
      plot_mz <- c(tmp_mzmin, tmp_mzmax)
    } else {
      tmp_mz_window <- ppm_to_num(tmp_mz_theory, opt$ppm_global)
      plot_mz <- c(tmp_mz_theory - tmp_mz_window, tmp_mz_theory + tmp_mz_window)
    }

    tmp_chr <- xcms::chromatogram(
      object = stds_exp,
      BPPARAM = bp_param,
      chunkSize = opt$cores,
      mz = plot_mz,
      aggregationFun = "sum"
    )

    tmp_tib <- purrr::map(
      .x = seq_len(ncol(tmp_chr)),
      .f = ~ {
        tibble::tibble(
          sample = colnames(tmp_chr)[.x],
          rt = tmp_chr[1, .x]@rtime,
          intensity = tidyr::replace_na(tmp_chr[1, .x]@intensity, replace = 0)
        )
      }
    ) %>%
      purrr::list_rbind() %>%
      dplyr::filter(sample %in% tmp_samples)

    plot_window <- if (!isTRUE(tmp_feature_found) || any(!is.finite(c(tmp_rtmin, tmp_rtmax)))) {
      NULL
    } else {
      c(
        max(
          min(tmp_tib$rt, na.rm = TRUE),
          tmp_rtmin - max(30, 5 * (tmp_rtmax - tmp_rtmin))
        ),
        min(
          max(tmp_tib$rt, na.rm = TRUE),
          tmp_rtmax + max(30, 5 * (tmp_rtmax - tmp_rtmin))
        )
      )
    }

    title_label <- paste0(
      tmp_group, ", ", tmp_adduct, ", ",
      "m/z: ", round(plot_mz[1], 4), " - ", round(plot_mz[2], 4), ", ",
      if (isTRUE(tmp_feature_found)) {
        paste0("peak rt: ", round(tmp_rt, 2))
      } else {
        "no detected feature"
      }
    )
    title_label_ascii <- iconv(
      title_label,
      from = "",
      to = "ASCII//TRANSLIT",
      sub = ""
    )
    if (!is.na(title_label_ascii)) {
      title_label <- gsub("[^ -~]", "", title_label_ascii)
    }

    tmp_p <- tmp_tib %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = rt,
          y = intensity,
          group = sample
        )
      ) +
      ggplot2::geom_line() +
      ggplot2::facet_wrap(~sample, scales = "fixed") +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0, 0.05))) +
      ggplot2::theme_bw() +
      ggplot2::theme(strip.background = ggplot2::element_blank()) +
      ggplot2::labs(title = title_label)

    if (!is.null(plot_window)) {
      tmp_p <- tmp_p + ggplot2::coord_cartesian(xlim = plot_window)
    }

    if (isTRUE(tmp_feature_found) && all(is.finite(c(tmp_rtmin, tmp_rtmax)))) {
      tmp_p <- tmp_p +
        ggplot2::geom_rect(
          inherit.aes = FALSE,
          xmin = tmp_rtmin,
          xmax = tmp_rtmax,
          ymin = -Inf,
          ymax = Inf,
          alpha = 0.1,
          fill = "red"
        )
    }

    if (isTRUE(tmp_feature_found) && is.finite(tmp_rt)) {
      tmp_p <- tmp_p +
        ggplot2::geom_vline(
          xintercept = tmp_rt,
          linetype = 2,
          color = "red"
        )
    }

    ggplot2::ggsave(
      filename = file.path(opt$output, "peaks", paste0(tmp_group, ".pdf")),
      plot = tmp_p,
      device = "pdf",
      height = 6,
      width = 6,
      units = "in"
    )
  }
)
cli::cli_progress_done()

cli::cli_progress_step("Save results to tables")
readr::write_csv(
  x = standards_all_features,
  file = file.path(opt$output, "tables", "standards_all_features.csv")
)
readr::write_csv(
  x = standards_best_features,
  file = file.path(opt$output, "tables", "standards_best_features.csv")
)
cli::cli_progress_done()
