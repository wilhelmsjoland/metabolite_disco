import_mzml <- function(
  data_path = NULL,
  meta_file = NULL
) {

  mzml_files <- list.files(
    data_path,
    # pattern = ".test",
    pattern = "\\.mzML$|\\.mzml$",
    full.names = TRUE,
    recursive = TRUE,
    ignore.case = TRUE
  )

  if (length(mzml_files) == 0) {
    stop("No .mzML files found in path")
  }

  mzml_tib <- data.frame(
    sample = basename(mzml_files),
    path = mzml_files
  )

  meta <- readr::read_csv(
    file = meta_file,
    show_col_types = FALSE,
    progress = FALSE
  )

  meta_matched <- dplyr::left_join(
    x = meta,
    y = mzml_tib,
    na_matches = "never",
    relationship = "one-to-one", 
    by = dplyr::join_by(sample)
  )

  samp_lacking_file <- meta_matched$sample[
    !meta_matched$sample %in% basename(mzml_files)
  ]
  if (length(samp_lacking_file) > 0) {
    cli::cli_abort(
      c(
        "i" = "Failed to load data",
        "x" = "File(s) not found: {.path {samp_lacking_file}}"
      )
    )
  }

  # non_used_files <- mzml_files[!mzml_files %in% meta_matched$path]
  # if (length(non_used_files) > 0) {
  #   cli::cli_alert_warning("Files in path not in metadata:")
  #   cli::cli_ul(non_used_files)
  # }

  meta_matched_filt <- meta_matched %>%
    dplyr::filter(!sample %in% samp_lacking_file) %>%
    tibble::column_to_rownames(var = "sample")

  return(meta_matched_filt)
}

multiply_chem <- function(chem_formula = NULL, multiply_by = 1) {
  split_form <- stringr::str_extract_all(
    chem_formula,
    "[A-Z][a-z]?[0-9]*"
  )[[1]]
  split_form_append <- stringr::str_replace_all(
    split_form,
    "([A-Za-z])$",
    "\\11"
  )
  multiply_by <- multiply_by
  for (i in seq_along(split_form_append)) {
    chem_letter <- stringr::str_extract_all(
      split_form_append[i],
      "[A-Za-z]+"
    )[[1]]
    chem_number <- as.integer(stringr::str_extract_all(
      split_form_append[i],
      "[0-9]+$"
    )[[1]])
    split_form_append[i] <- paste0(chem_letter, chem_number * multiply_by)
  }
  split_form_append <- stringr::str_flatten(split_form_append)
  clean_form <- Rdisop::getFormula(Rdisop::getMolecule(split_form_append))
  return(clean_form)
}

find_intersect_feat <- function(data = NULL, set = NULL, full_set = NULL) {
  filt_data <- data %>%
    dplyr::filter(
      dplyr::if_all(
        .cols = tidyselect::all_of(full_set[full_set %in% set]),
        .fns = ~ .x == TRUE
      ),
      dplyr::if_all(
        .cols = tidyselect::all_of(full_set[!full_set %in% set]),
        .fns = ~ .x == FALSE
      )
    )
  return(filt_data)
}

find_y_position <- function(
  test_df,
  df,
  formula,
  fun_data,
  grouping = NULL
) {

  if (!fun_data %in% c(
    "max",
    "mean",
    "min",
    "mean_sem",
    "max_sem",
    "lower_bound",
    "upper_bound"
  )) {
    stop("Incorrect specification of fun_data")
  }

  num_col <- stringr::str_trim(stringr::str_split_i(formula, "~", 1))
  group_col <- stringr::str_trim(stringr::str_split_i(formula, "~", 2))

  sum_df <- df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::summarize(
      max = max(.data[[num_col]], na.rm = TRUE),
      mean = mean(.data[[num_col]], na.rm = TRUE),
      min = min(.data[[num_col]], na.rm = TRUE),
      mean_sem = mean + sd(.data[[num_col]]) / sqrt(length(.data[[num_col]])),
      max_sem = max + sd(.data[[num_col]]) / sqrt(length(.data[[num_col]])),
      Q1 = quantile(.data[[num_col]], 0.25, na.rm = TRUE),
      Q3 = quantile(.data[[num_col]], 0.75, na.rm = TRUE),
      IQR.val = IQR(.data[[num_col]], na.rm = TRUE),
      lower_bound = Q1 - 1.5 * IQR.val,
      upper_bound = Q3 + 1.5 * IQR.val
    )

  sum_df <- sum_df %>%
    dplyr::select(.data[[group_col]], tidyselect::all_of(fun_data)) %>%
    dplyr::rename("y.pos" = .data[[fun_data]])

  group1_names <- setNames(group_col, "group1")
  group2_names <- setNames(group_col, "group2")

  final_df <- test_df %>%
    dplyr::left_join(
      x = .,
      y = sum_df,
      by = group1_names
    ) %>%
    dplyr::left_join(
      x = .,
      y = sum_df,
      by = group2_names
    ) %>%
    {
      if (!is.null(grouping)) {
        dplyr::group_by(., dplyr::across(tidyselect::all_of(grouping)))
      } else {
        .
      }
    } %>%
    # dplyr::group_by(comparisons) %>%
    dplyr::mutate(y.pos = base::max(y.pos.x, y.pos.y, na.rm = TRUE)) %>%
    dplyr::select(-c("y.pos.x", "y.pos.y"))

  return(final_df)
}

import_biotransform_meta <- function(file = NULL) {
  biotransf_meta <- readr::read_csv(
    file = file,
    comment = "#",
    show_col_types = FALSE,
    progress = FALSE
  ) %>%
    tidyr::uncount(
      .,
      weights = allowed_n,
      .id = "multiplier",
      .remove = FALSE
    ) %>%
    dplyr::mutate(name = paste0(multiplier, " x ", name)) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(delta_formula = multiply_chem(delta_formula, multiplier)) %>%
    dplyr::mutate(delta_mass = get_theory_deltamass(delta_formula)) %>%
    # dplyr::mutate(delta_formula = paste0("± ", delta_formula)) %>%
    dplyr::ungroup()

  return(biotransf_meta)
}

ppm_to_num <- function(mz, ppm) {
  abs(mz * ppm * 1e-06)
}

num_to_ppm <- function(mz, diff) {
  abs(diff / mz * 1e06)
}

get_theory_mz <- function(
  chem_form = NULL,
  adduct = NULL
) {
  chem_mass <- MetaboCoreUtils::calculateMass(chem_form)
  chem_theory_mz <- MetaboCoreUtils::mass2mz(chem_mass, adduct)
  chem_theory_mz <- chem_theory_mz[1, 1]

  return(chem_theory_mz)
}

get_theory_deltamass <- function(chem_form = NULL) {
  chem_mass <- MetaboCoreUtils::calculateMass(chem_form)
  chem_theory_mass <- unname(chem_mass)

  return(chem_theory_mass)
}

get_short_mz_range <- function(
  chem_theory_mz = NULL,
  mz_window = 0.01
) {
  chem_mz_range <- c(chem_theory_mz - mz_window, chem_theory_mz + mz_window)
  return(chem_mz_range)
}

get_rt_mz_range <- function(
  chromatogram = NULL,
  rt_window = 0.02 # percentage increase of window
) {
  ref_peak_intensity <- purrr::map_dbl(
    .x = seq_len(ncol(chromatogram)),
    .f = ~ max(intensity(chromatogram[1, .x]), na.rm = TRUE)
  ) %>% median(na.rm = TRUE)

  ref_peak_rt <- purrr::map_dbl(
    .x = seq_len(ncol(chromatogram)),
    .f = ~ {
      samp_max <- max(intensity(chromatogram[1, .x]), na.rm = TRUE)
      idx <- which(intensity(chromatogram[1, .x]) == samp_max)
      Spectra::rtime(chromatogram[1, .x])[[idx]]
    }
  ) %>% median(na.rm = TRUE)

  rt_window_min <- (1 - rt_window) * ref_peak_rt
  rt_window_max <- (1 + rt_window) * ref_peak_rt

  ref_tib <- tibble::tibble()
  for (samp_n in 1:length(chromatogram)) {
    #ref_peak_intensity <- max(intensity(chromatogram[1, samp_n]), na.rm = TRUE)
    samp_peak_intensity <- max(intensity(chromatogram[1, samp_n]), na.rm = TRUE)
    max_intensity_idx <- which(
      intensity(chromatogram[1, samp_n]) == samp_peak_intensity
    )
    samp_peak_rt <- Spectra::rtime(chromatogram[1, samp_n])[[max_intensity_idx]]
    mzmin <- min(Spectra::mz(chromatogram[1, samp_n]), na.rm = TRUE)
    mzmax <- max(Spectra::mz(chromatogram[1, samp_n]), na.rm = TRUE)
    file <- colnames(chromatogram)[samp_n]

    tmp_tib <- tibble::tibble(
      file,
      ref_peak_intensity,
      ref_peak_rt,
      samp_peak_intensity,
      samp_peak_rt,
      rt_window_min,
      rt_window_max,
      max_intensity_idx,
      mzmin,
      mzmax,
      samp_n
    )

    ref_tib <- dplyr::bind_rows(ref_tib, tmp_tib)
  }

  rt_range <- c(min(ref_tib$rt_window_min), max(ref_tib$rt_window_max))
  mz_range <- c(min(ref_tib$mzmin), max(ref_tib$mzmax))

  ranges <- list(
    "mz_range" = mz_range,
    "rt_range" = rt_range,
    "data" = ref_tib
  )

  return(ranges)
}

inspect_peak <- function(
  chrom_obj = NULL,
  peak = TRUE,
  save_loc = NULL
) {

  peak_tib <- tibble::as_tibble(
    xcms::chromPeaks(chrom_obj),
    rownames = "peak"
  )
  peak_row  <- peak_tib[peak_tib$peak == peak, ]
  filt_chr <- chrom_obj[peak_row$sample] %>%
    Spectra::filterMzRange(
      mz = c(
        peak_row$mzmin,
        peak_row$mzmax
      )
    ) %>%
    Spectra::filterRt(
      rt = c(
        peak_row$rtmin,
        peak_row$rtmax
      )
    ) %>%
    xcms::chromatogram(
      BPPARAM = BiocParallel::SerialParam(),
      chunkSize = 1L,
      aggregationFun = "sum"
    )

  # Different ways of doing the same thing
  # basename(fileNames(xchr)[peak_row$sample])
  # rownames(meta)[peak_row$sample]
  samp_name <- rownames(MsExperiment::sampleData(chrom_obj)[peak_row$sample, ])
  samp_group <- meta$group[rownames(meta) == samp_name]
  samp_color <- group_colors[names(group_colors) == samp_group]

  if (!is.null(save_loc)) {
    pdf(
      file.path(
        save_loc,
        paste0(
          peak,
          ".pdf"
        )
      )
    )
    plot(
      x = filt_chr,
      main = paste0(
        "Peak: ", peak_row$peak,
        "\n m/z: ", round(peak_row$mzmin, 3),
        " - ", round(peak_row$mzmax, 3),
        "\n RT: ", round(peak_row$rtmin, 3),
        " - ", round(peak_row$rtmax, 3)
      ),
      sub = paste0("Sample: ", samp_name),
      peakType = "polygon",
      peakCol = samp_color,
      peakBg = NA,
      lwd = 3
    )
    legend(
      x = "topleft",
      legend = samp_group,
      col = samp_color,
      lty = 1,
      lwd = 3
    )
  } else {
    plot(
      x = filt_chr,
      main = paste0(
        "Peak: ", peak_row$peak,
        "\n m/z: ", round(peak_row$mzmin, 3),
        " - ", round(peak_row$mzmax, 3),
        "\n RT: ", round(peak_row$rtmin, 3),
        " - ", round(peak_row$rtmax, 3)
      ),
      sub = paste0("Sample: ", samp_name),
      peakType = "polygon",
      peakCol = samp_color,
      peakBg = NA,
      lwd = 3
    )

    legend(
      x = "topleft",
      legend = samp_group,
      col = samp_color,
      lty = 1,
      lwd = 3
    )
  }

  dev.off()
  return(filt_chr)
}

plot_experiment_intensities <- function(
  chrom_obj = NULL,
  value = "into", # into, intb, maxo
  save_loc = NULL,
  device = "pdf"
) {
  # Extract a list of per-sample peak intensities (in log2 scale)
  # TODO Do the samples match the colors?
  xchr_peaks <- tibble::as_tibble(
    xcms::chromPeaks(chrom_obj),
    rownames = "peak"
  )
  xchr_meta <- tibble::as_tibble(
    MsExperiment::sampleData(chrom_obj),
    rownames = "sample_id"
  ) %>%
    dplyr::mutate(sample = dplyr::row_number())

  xchr_data_comb <- dplyr::left_join(
    x = xchr_peaks,
    y = xchr_meta,
    by = "sample"
  ) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = c("into", "intb"),
        .fns = ~ log2(.)
      )
    )

  # lower number of detected peaks = smaller width of the boxes
  xchr_data_p <- xchr_data_comb %>%
    dplyr::mutate(comb = paste0(sample, "_", sample_id)) %>%
    dplyr::mutate(comb = forcats::fct_reorder(
      .f = comb,
      .x = sample
    )) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = comb,
        y = {{ value }},
        fill = group
      )
    ) +
    ggplot2::geom_boxplot(varwidth = TRUE) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      y = "Log2 intensity",
      subtitle = paste0(
        deparse(substitute(value)),
        ", box width: number of detected peaks"
      )
    )

  file_nm <- paste0(save_loc, "/chrom_detected_peaks.", device)

  if (!is.null(save_loc)) {
    ggplot2::ggsave(
      filename = file_nm,
      plot = xchr_data_p,
      device = device,
      height = 6,
      width = 6,
      units = "in"
    )
  }

  return(xchr_data_p)
}

plot_twenty_feats <- function(
  chromatogram = NULL,
  save_loc = "output/apiin_bu_25_ppm/graphs/features/"
) {
  chrom_feats <- xcms::featureDefinitions(chromatogram) %>%
    tibble::as_tibble(
      x = .,
      rownames = "feature"
    ) %>%
    dplyr::arrange(mzmed)

  splitted <- split(
    chrom_feats$row,
    ceiling(seq_along(chrom_feats$row) / 20)
  )

  for (i in seq_along(splitted)) {
    idx_feats <- dplyr::filter(chrom_feats, row %in% splitted[[i]])
    min_mz <- min(idx_feats$mzmed, na.rm = TRUE)
    max_mz <- max(idx_feats$mzmed, na.rm = TRUE)
    file_name <- paste0(round(min_mz, 3), "-", round(max_mz, 3))
    plot_name <- paste0(save_loc, file_name, ".pdf")
    if (interactive() && file.exists(plot_name)) {
      cli::cli_alert_info(
        "{.path {plot_name}} already exists, skipping"
      )
      next
    } else {
      sliced_chr <- chromatogram[splitted[[i]], ]
      peaks <- xcms::chromPeaks(sliced_chr)
      peak_colors <- group_colors[
        Biobase::pData(chromatogram)$group[peaks[, "column"]]
      ]

      pdf(
        file = plot_name,
        width = 12,
        height = 10,
      )

      plot(
        sliced_chr,
        peakType = "polygon",
        col = group_colors[Biobase::pData(chromatogram)$group],
        peakCol = peak_colors,
        peakBg = NA,
        lwd = 3
      )
      invisible(dev.off())
      cli::cli_alert_success(
        "Saved {.path {plot_name}}"
      )
    }
  }
}

plot_feat_chrom_int <- function(
  feature_chrom = NULL,
  feature = NULL,
  method = NULL,
  value = NULL,
  filled = FALSE,
  missing = NULL,
  ms_level = 1,
  save_loc = NULL,
  device = "pdf",
  feat_pairs = FALSE,
  overwrite = FALSE
) {

  file_nm <- paste0(save_loc, feature, ".", device)
  if (
    interactive() && file.exists(file_nm) &&
      !overwrite &&
      !is.null(save_loc)
  ) {
    cli::cli_alert_info(
      "{.path {file_nm}} already exists, skipping"
    )
    return(invisible(NULL))
  }

  feat_tib <- tibble::as_tibble(
    x = xcms::featureDefinitions(feature_chrom),
    rownames = "feature"
  )
  feat_idx <- which(feat_tib$feature == feature)
  lone_feat <- feature_chrom[feat_idx, ]
  lone_feat_def <- tibble::as_tibble(
    x = xcms::featureDefinitions(lone_feat),
    rownames = "feature"
  )

  # Only take the ones with peaks or the coloring in the base plot is wrong
  keep_peaks <- xcms::hasChromPeaks(lone_feat)
  lone_feat_peaks <- xcms::hasChromPeaks(lone_feat)[, keep_peaks]
  base_group_colors <- meta$group[match(names(lone_feat_peaks), meta$sample)]

  if (isTRUE(feat_pairs)) {
    p1 <- function() {
      par(
        mar = c(bottom = 2.3, left = 0, top = 0, right = 0),
        mgp = c(1.5, 0.5, 0),
        tck = -0.02
      )
      plot(
        lone_feat,
        peakType = "polygon",
        peakCol = group_colors[base_group_colors],
        peakBg = NA,
        lwd = 3,
        main = "",
        ylab = "",
        xlab = "",
        xaxt = "n",
        yaxt = "n"
      )
      axis(1, las = 0, cex.axis = 0.8)
      axis(2, las = 2, cex.axis = 0.8)
      mtext("Intensity", side = 2, line = 2.8)
      mtext("Retention time", side = 1, line = 1.5)
    }
  } else {
    p1 <- function() {
      par(mar = c(5.1, 4.1, 4.1, 8.6))
      plot(
        lone_feat,
        peakType = "polygon",
        peakCol = group_colors[base_group_colors],
        peakBg = NA,
        lwd = 3,
        main = paste0(
          lone_feat_def$feature, ", ",
          "m/z: ",
          round(lone_feat_def$mzmin, 3),
          " - ",
          round(lone_feat_def$mzmax, 3),
          ", RT: ",
          round(lone_feat_def$rtmin, 2),
          " - ",
          round(lone_feat_def$rtmax, 2)
        ),
        ylab = "",
        xlab = "",
        xaxt = "n",
        yaxt = "n"
      )
      axis(1, las = 0, cex.axis = 0.8)
      axis(2, las = 2, cex.axis = 0.8)
      mtext("Intensity", side = 2, line = 3.2)
      mtext("Retention time (s)", side = 1, line = 2.5)
    }
  }

  p2_data <- xcms::featureValues(
    lone_feat,
    method = method,
    value = value,
    filled = filled,
    missing = missing,
    ms_level = ms_level
  ) %>%
    tibble::as_tibble(., rownames = "feature") %>%
    tidyr::pivot_longer(cols = contains(".mzML")) %>%
    dplyr::left_join(
      x = .,
      y = dplyr::select(meta, sample, group),
      by = c("name" = "sample"),
      # y = tibble::as_tibble(x = meta, rownames = "file") %>%
      # dplyr::select(file, group),
      # by = c("name" = "file")
    )

  p2_data_signif <- full_limma %>%
    dplyr::select(feature, adj.P.Val, contrast) %>%
    dplyr::mutate(
      group1 = stringr::str_split_i(contrast, "-", 1),
      group2 = stringr::str_split_i(contrast, "-", 2),
    ) %>%
    dplyr::select(-contrast) %>%
    dplyr::filter(feature %in% unique(p2_data$feature)) %>%
    rstatix::add_significance(p.col = "adj.P.Val") %>%
    find_y_position(
      test_df = .,
      df = p2_data,
      formula = "value ~ group",
      fun_data = "max"
    )

  p2_data_signif_only <- p2_data_signif %>%
    dplyr::filter(!adj.P.Val.signif %in% c("ns"))

  p2 <- p2_data %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = group,
        y = value
      )
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = group),
      outliers = FALSE
    ) +
    ggplot2::geom_point(
      position = ggplot2::position_jitter(width = 0.15),
      size = 2,
      pch = 21,
      color = "black"
    ) +
    ggplot2::scale_y_continuous(expand =  ggplot2::expansion(c(0.1, 0.1))) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
      y = paste0("Peak area (", value, ")"),
      caption = paste0(
        "method: ", method,
        ", value: ", value,
        ", filled: ", filled,
        ", missing: ", missing
      )
    ) + {
    if (nrow(p2_data_signif_only) > 0) {
      ggpubr::geom_bracket(
        data = p2_data_signif_only,
        ggplot2::aes(
          xmin = group1,
          xmax = group2,
          label = adj.P.Val.signif,
          y.position = y.pos * 1.005
        ),
        step.increase = 0.15,
        vjust = 0.1
      )
    } else {
      NULL
    }
  }

  p3 <- cowplot::plot_grid(
    p1,
    p2,
    ncol = 1,
    align = "v",
    axis = "bt",
    rel_heights = c(1.5, 1)
  )

  if (!is.null(save_loc)) {
    ggplot2::ggsave(
      filename = file_nm,
      plot = p3,
      device = device,
      height = 6,
      width = 6,
      units = "in"
    )
    cli::cli_alert_success(
      "Saved {.path {file_nm}}"
    )
  }

  data_list <- list(
    "combined" = p3,
    "chromatogram" = p1,
    "boxplot" = p2,
    "boxplot_data" = p2_data,
    "boxplot_signif" = p2_data_signif
  )

  return(data_list)
}

plot_feature_pairs <- function(
  feature_chrom = NULL,
  filt_match_row = NULL,
  method = NULL,
  value = NULL,
  filled = FALSE,
  missing = NULL,
  ms_level = 1,
  save_pairs_loc = NULL,
  device = "pdf",
  overwrite = FALSE
) {

  ft_pair <- filt_match_row[["pair"]][[1]]
  file_nm <- paste0(
    save_pairs_loc,
    ft_pair[1], "_", ft_pair[2], ".", device
  )

  if (interactive() && file.exists(file_nm) && !overwrite) {
    cli::cli_alert_info(
      "{.path {file_nm}} already exists, skipping"
    )
    return(invisible(NULL))
  }

  feat1 <- plot_feat_chrom_int(
    feature_chrom = feature_chrom,
    feature = ft_pair[1],
    method = method,
    value = value,
    filled = filled,
    missing = missing,
    ms_level = ms_level,
    save_loc = NULL,
    device = NULL,
    feat_pairs = TRUE
  )

  feat2 <- plot_feat_chrom_int(
    feature_chrom = feature_chrom,
    feature = ft_pair[2],
    method = method,
    value = value,
    filled = filled,
    missing = missing,
    ms_level = ms_level,
    save_loc = NULL,
    device = NULL,
    feat_pairs = TRUE
  )

  ft_pair_p <- ((
    cowplot::ggdraw(feat1$chromatogram) /
      (feat1$boxplot + ggplot2::labs(caption = ""))
  ) |
    (
      cowplot::ggdraw(feat2$chromatogram) /
        (feat2$boxplot + ggplot2::labs(caption = ""))
    )
  ) +
    patchwork::plot_layout(
      guides = "collect",
      axes = "collect",
      heights = c(1.3, 0.7)
    ) +
    patchwork::plot_annotation(
      title = paste0(
        ft_pair[1], " & ",
        ft_pair[2], "\n",
        "Potential ", filt_match_row$name, " ", filt_match_row$chem_change
      )
    )

  if (!is.null(save_pairs_loc)) {
    ggplot2::ggsave(
      filename = file_nm,
      plot = ft_pair_p,
      device = device,
      height = 7,
      width = 7,
      units = "in"
    )
    cli::cli_alert_success(
      "Saved {.path {file_nm}}"
    )
  }

  return(ft_pair_p)
}

feat_to_idx <- function(feature_idx = NULL) {
  clean_idx <- as.numeric(gsub("[A-Za-z]", "", feature_idx))
  return(clean_idx)
}

register_parallel <- function(workers = NULL) {
  if (is.null(workers) || workers < 0) {
    workers <- parallel::detectCores() - 1
  }
  sys <- Sys.info()["sysname"]
  if (sys == "Windows") {
    bp <<- BiocParallel::SnowParam(
      workers = workers,
      type = "SOCK"
    )
  } else if (sys %in% c("Linux", "Darwin")) {
    bp <<- BiocParallel::MulticoreParam(
      workers = workers
    )
  }

  BiocParallel::register(bp, default = TRUE)
}

filt_features_old <- function(
  object = NULL,
  beta_cor_threshold = 0.3,
  beta_snr_threshold = 6,
  sn_threshold = 10,
  filt_vector = NULL
) {
  filt_chrompeaks_tib <- tibble::as_tibble(
    xcms::chromPeaks(object),
    rownames = "feature"
  ) %>%
    dplyr::filter(!is.na(beta_cor) & !is.na(beta_snr)) %>%
    dplyr::group_by(feature) %>%
    dplyr::filter(
      beta_cor >= beta_cor_threshold &
        beta_snr >= beta_snr_threshold
    ) %>%
    dplyr::filter(sn >= sn_threshold) %>%
    dplyr::mutate(feature2 = as.numeric(gsub("[A-Za-z]", "", feature))) %>%
    dplyr::relocate(feature2, .after = feature)

  filt_features_tib <- tibble::as_tibble(
    xcms::featureDefinitions(object),
    rownames = "feature"
  ) %>%
    tidyr::unnest(peakidx) %>%
    dplyr::filter(peakidx %in% filt_chrompeaks_tib$feature2) %>%
    tidyr::nest(data = peakidx)

  filt_sig_features_tib <- filt_features_tib %>%
    # Added so we don't remove interesting ones
    # that don't have a predicted biotransformation
    dplyr::filter(feature %in% filt_vector)

  filt_sig_features <- filt_sig_features_tib$feature

  # Now only look at features that have at least
  # one significantly different feature
  biot_filt_sig_feature_tib <- matched_diffs %>%
    dplyr::filter(
      dplyr::if_any(
        tidyselect::all_of(c("feat1", "feat2")),
        ~ .x %in% filt_sig_features
      )
    ) %>%
    dplyr:::mutate(
      pair = purrr::map2(feat1, feat2, ~ c(.x, .y)),
      # or ppm global
      mz1_forms = purrr::map(
        mz1, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 0))
      ),
      mz2_forms = purrr::map(
        mz2, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 0))
      )
    )

  biot_filt_sig_features <- unique(
    c(
      biot_filt_sig_feature_tib$feat1,
      biot_filt_sig_feature_tib$feat2
    )
  )

  final_plotting_features <- unique(
    c(
      filt_sig_features,
      biot_filt_sig_features
    )
  )

  filt_list <- list(
    # filtering features below:

    # quality filtered peak tib
    "filt_chrompeaks_tib" = filt_chrompeaks_tib,
    # quality filtered feature tib
    "filt_features_tib" = filt_features_tib,
    # quality + sig filtered feature tib
    "filt_sig_features_tib" = filt_sig_features_tib,
    # quality + sig filtered features
    "filt_sig_features" = filt_sig_features,
    # biotransformation features below:

    # quality + sig filtered biotransf tib
    "biot_filt_sig_features_tib" = biot_filt_sig_feature_tib,
    # quality + sig filtered feature
    "biot_filt_sig_features" = biot_filt_sig_features,
    # final features for plotting: below

    # all feats in biot and in filt_sig_features
    "final_plotting_features" = final_plotting_features
  )

  return(filt_list)

}

plot_pca <- function(
  prcomp_res = NULL,
  metad = NULL,
  x = NULL,
  y = NULL
) {
  scores <- tibble::as_tibble(prcomp_res$x, rownames = "sample")

  scores <- scores %>%
    dplyr::left_join(
      x = .,
      y = tibble::as_tibble(metad, rownames = "sample"),
      by = "sample"
    )

  # % variance explained for axis labels
  var_expl <- (prcomp_res$sdev^2) / sum(prcomp_res$sdev^2) * 100
  names(var_expl) <- paste0("PC", seq_along(var_expl))

  pca_p <- scores %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = {{ x }},
        y = {{ y }}
      )
    ) +
    ggplot2::geom_point(aes(color = group), size = 3) +
    ggplot2::scale_color_manual(values = group_colors) +
    ggplot2::labs(
      x = sprintf(
        "%s (%.2f%%)",
        deparse(substitute(x)),
        var_expl[[deparse(substitute(x))]]
      ),
      y = sprintf(
        "%s (%.2f%%)",
        deparse(substitute(y)),
        var_expl[[deparse(substitute(y))]]
      )
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.title = ggplot2::element_blank())

  return(pca_p)
}

pred_biot <- function(
  data = NULL,
  biotransf_data = NULL,
  tolerance_ppm = NULL,
  tolerance = NULL,
  features_of_interest = NULL,
  parallel = TRUE,
  n_workers = snakemake@config$cores
) {
  if (is.null(tolerance) & is.null(tolerance_ppm)) {
    stop("Both tolerance and tolerance_ppm are NULL")
  } else if (!is.null(tolerance) & !is.null(tolerance_ppm)) {
    stop("Either tolerance or tolerance_ppm need to be set to NULL")
  }

  ## 1. Prepare peaks table (all peaks, sorted by m/z)
  peaks <- data %>%
    dplyr::arrange(mass) %>%
    dplyr::mutate(
      peak_id = dplyr::row_number(),
      feature = feature,
      mass = mass
    )
  n_peaks <- nrow(peaks)

  # convenience vectors
  mz_vec <- peaks$mzmed
  rt_vec <- peaks$rtmed
  adduct_vec <- peaks$adduct
  mass_vec <- peaks$mass
  id_vec <- peaks$peak_id
  feat_vec <- peaks$feature
  n_trans <- nrow(biotransf_data)

  # ============================================================================
  # NEW: indices of features of interest (if provided)
  # ============================================================================
  if (!is.null(features_of_interest)) {
    foi_idx <- which(feat_vec %in% features_of_interest)
    if (length(foi_idx) == 0) {
      stop("None of the features_of_interest found in data")
    }
  }
  # ============================================================================

  # Calculate tolerance
  if (is.null(tolerance)) {
    tol_used <- MsCoreUtils::ppm(mz_vec, tolerance_ppm)
  } else {
    tol_used <- rep(tolerance, n_peaks)
  }

  ## 3. For each biotransformation, find all matching peak pairs
  if (parallel && n_trans > 1) {
    # Set up future plan
    oplan <- future::plan(future::multisession, workers = n_workers)
    on.exit(future::plan(oplan), add = TRUE)

    all_matches <- future.apply::future_lapply(1:n_trans, function(k) {
      delta <- biotransf_data$delta_mass[k]
      this_name <- biotransf_data$name[k]
      this_formula <- biotransf_data$delta_formula[k]

      target_lower <- mass_vec + delta - tol_used
      target_upper <- mass_vec + delta + tol_used

      idx_start <- findInterval(target_lower, mass_vec) + 1L
      idx_end <- findInterval(target_upper, mass_vec)

      res_list <- vector("list", n_peaks)

      # Only iterate over features of interest if provided, otherwise all peaks
      iter_idx <- if (!is.null(features_of_interest)) {
        foi_idx
      } else {
        1:n_peaks
      }

      for (i in iter_idx) { # 1:n_peaks
        start_idx <- idx_start[i]
        end_idx <- idx_end[i]

        if (start_idx > end_idx) next

        partner_idx <- seq.int(start_idx, end_idx)

        # """"""""""""""""""""""""""""""
        # Also exclude partners with same feature name (same feature,
        # different adduct)
        partner_idx <- if (!is.null(features_of_interest)) {
          partner_idx[partner_idx != i & feat_vec[partner_idx] != feat_vec[i]]
        } else {
          partner_idx[partner_idx > i & feat_vec[partner_idx] != feat_vec[i]]
        }
        # """"""""""""""""""""""""""""""

        if (length(partner_idx) == 0) next

        res_list[[i]] <- tibble::tibble(
          name = this_name,
          chem_change = this_formula,
          feat1 = feat_vec[i],
          feat2 = feat_vec[partner_idx],
          mz1 = mz_vec[i],
          mz2 = mz_vec[partner_idx],
          adduct1 = adduct_vec[i],
          adduct2 = adduct_vec[partner_idx],
          mass1 = mass_vec[i],
          mass2 = mass_vec[partner_idx],
          rt1 = rt_vec[i],
          rt2 = rt_vec[partner_idx],
          delta_mass = delta,
          obs_delta_mass = mass_vec[partner_idx] - mass_vec[i],
          peak1_id = id_vec[i],
          peak2_id = id_vec[partner_idx]
        )
      }
      dplyr::bind_rows(res_list)
    })
  } else {
    # Sequential version
    all_matches <- vector("list", n_trans)

    for (k in 1:n_trans) {
      delta <- biotransf_data$delta_mass[k]
      this_name <- biotransf_data$name[k]
      this_formula <- biotransf_data$delta_formula[k]

      target_lower <- mass_vec + delta - tol_used
      target_upper <- mass_vec + delta + tol_used

      idx_start <- findInterval(target_lower, mass_vec) + 1L
      idx_end <- findInterval(target_upper, mass_vec)

      res_list <- vector("list", n_peaks)

      # Only iterate over features of interest if provided, otherwise all peaks
      iter_idx <- if (!is.null(features_of_interest)) {
        foi_idx
      } else {
        1:n_peaks
      }

      for (i in iter_idx) { # 1:n_peaks
        start_idx <- idx_start[i]
        end_idx <- idx_end[i]

        if (start_idx > end_idx) next

        partner_idx <- seq.int(start_idx, end_idx)

        # """"""""""""""""""""""""""""""
        # Also exclude partners with same feature name (same feature,
        # different adduct)
        partner_idx <- if (!is.null(features_of_interest)) {
          partner_idx[partner_idx != i & feat_vec[partner_idx] != feat_vec[i]]
        } else {
          partner_idx[partner_idx > i & feat_vec[partner_idx] != feat_vec[i]]
        }
        # """"""""""""""""""""""""""""""

        if (length(partner_idx) == 0) next

        res_list[[i]] <- tibble::tibble(
          name = this_name,
          chem_change = this_formula,
          feat1 = feat_vec[i],
          feat2 = feat_vec[partner_idx],
          mz1 = mz_vec[i],
          mz2 = mz_vec[partner_idx],
          adduct1 = adduct_vec[i],
          adduct2 = adduct_vec[partner_idx],
          mass1 = mass_vec[i],
          mass2 = mass_vec[partner_idx],
          rt1 = rt_vec[i],
          rt2 = rt_vec[partner_idx],
          delta_mass = delta,
          obs_delta_mass = mass_vec[partner_idx] - mass_vec[i],
          peak1_id = id_vec[i],
          peak2_id = id_vec[partner_idx]
        )
      }
      all_matches[[k]] <- dplyr::bind_rows(res_list)
    }
  }

  # 4. Final table of all matched pairs
  matched_diffs <- dplyr::bind_rows(all_matches)
  return(matched_diffs)
}

filt_features <- function(
  object = NULL,
  beta_cor_threshold = 0.8,
  beta_snr_threshold = 3,
  sn_threshold = 10
) {
  filt_chrompeaks_tib <- tibble::as_tibble(
    xcms::chromPeaks(object),
    rownames = "feature"
  ) %>%
    dplyr::filter(!is.na(beta_cor) & !is.na(beta_snr)) %>%
    dplyr::group_by(feature) %>%
    dplyr::filter(
      beta_cor >= beta_cor_threshold &
        beta_snr >= beta_snr_threshold
    ) %>%
    dplyr::filter(sn >= sn_threshold) %>%
    dplyr::mutate(feature2 = as.numeric(gsub("[A-Za-z]", "", feature))) %>%
    dplyr::relocate(feature2, .after = feature)

  filt_features_tib <- tibble::as_tibble(
    xcms::featureDefinitions(object),
    rownames = "feature"
  ) %>%
    tidyr::unnest(peakidx) %>%
    dplyr::filter(peakidx %in% filt_chrompeaks_tib$feature2) %>%
    tidyr::nest(data = peakidx)

  return(filt_features_tib)

}

start_log <- function(output_folder = NULL) {
  Sys.setlocale("LC_CTYPE", ".utf8")
  options(
    cli.ansi = FALSE,
    cli.unicode = TRUE,
    cli.width = 80,
    width = 80,
    crayon.enabled = FALSE,
    readr.show_progress = FALSE
  )

  if (!interactive()) {
    log_file <- file.path(
      output_folder,
      "logs",
      paste0(
        basename(output_folder),
        "_pipeline",
        ".log"
      )
    )
    con <<- file(log_file, open = "at", encoding = "UTF-8")
    sink(con, type = "output", append = TRUE)
    sink(con, type = "message", append = TRUE)
  }
}

script_header <- function() {
  rule_name <- gsub(
    "\\d+_|\\.rds$", "",
    basename(snakemake@output[[1]])
  )
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cli::cli_h1("{rule_name} {.emph {timestamp}}")
}


end_log <- function() {
  # cli::cli_alert_success(
  #   "Pipeline finished on {.time {format(Sys.time())}}"
  # )

  if (!interactive()) {
    sink(type = "message")
    sink(type = "output")
    close(con)
  }

  options(
    cli.width = 220,
    width = 220
  )
}

start_pipeline_msg <- function() {
  block_rule(col = col_cyan)
  cli::cli_rule(center = "Metabolite Disco")
  block_rule(col = col_cyan)
  cli::cli_bullets(
    c(
      " " = "",
      "i" = "Creator: {.emph Wilhelm Sjöland}",
      "i" = "Email: {.email wilhelm.sjoland@wlab.gu.se}",
      "i" = "Version: {.val {0.1}}",
      "i" = "Started pipeline: {.time {format(Sys.time())}}"
      # "i" = "A URL: {.url https://acme.com}"
    )
  )
}

block_rule <- function(col = col_white) {
  rule <- paste(
    rep(
      cli::symbol$full_block,
      cli::console_width()
    ),
    collapse = ""
  )
  cli::cli_text("{col(rule)}")
}

# kind of confusing I'll admit
param_msg <- function(process_history = NULL) {
  purrr::walk(
    .x = slotNames(process_history),
    .f = ~ cli::cli_bullets(
      c(
        "i" = paste0(
          .x, ": ",
          paste( # needed to not repeat the names twice for > 1 vectors
            slot(process_history, .x),
            collapse = ", "
          )
        )
      )
    )
  )
}

median_scale_tidy <- function(
  res_obj = NULL,
  assay = NULL
) {
  # Compute median and generate normalization factor
  mdns <- purrr::map_dbl(
    .x = as.data.frame(SummarizedExperiment::assay(res_obj, assay)),
    .f = ~ {
      median(.x, na.rm = TRUE)
    }
  )

  nf_mdn <- mdns / median(mdns)

  # Dividing dataset by median of median and creating a new assay
  med_scaled <- as.data.frame(
    SummarizedExperiment::assay(
      res_obj,
      assay
    )
  ) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::everything(),
        .fns = ~ . / nf_mdn[dplyr::cur_column()]
      )
    ) %>%
    as.matrix()

  cli::cli_alert_success(
    paste0(
      "assay ({.val {assay}}) median scaled"
    )
  )
  return(med_scaled)
}

median_scale_base <- function(
  res_obj = NULL,
  assay = NULL
) {
  # Compute median and generate normalization factor
  mdns <- base::apply(
    X = SummarizedExperiment::assay(res_obj, assay),
    MARGIN = 2,
    FUN = function(x) {
      median(x, na.rm = TRUE)
    }
  )
  nf_mdn <- mdns / median(mdns)

  # Dividing dataset by median of median and creating a new assay
  med_scaled <- base::sweep(
    x = SummarizedExperiment::assay(res_obj, assay),
    MARGIN = 2,
    STATS = nf_mdn,
    FUN = "/"
  )

  cli::cli_alert_success(
    paste0(
      "assay ({.val {assay}}) median scaled"
    )
  )
  return(med_scaled)
}

produce_complex_upset <- function(
  input = upset_intersect,
  comps = comparisons,
  qvalue = snakemake@params$qvalue
) {
  tmp_p <- ComplexHeatmap::UpSet(
    input,
    set_order = comps,
    comb_order = order(comb_size(input), decreasing = TRUE),
    column_title = paste0("Features with p adjusted < ", qvalue),
    row_names_max_width = ComplexHeatmap::max_text_width(
      ComplexHeatmap::set_name(
        input
      )
    ),
    # Numbers on the intersection size bars (top)
    top_annotation = ComplexHeatmap::HeatmapAnnotation(
      "Intersection size" = ComplexHeatmap::anno_barplot(
        comb_size(input),
        add_numbers = TRUE,
        border = FALSE,
        height = unit(14, "cm")
      ),
      annotation_name_side = "left"
    ),
    # Set size bars on the right with numbers
    left_annotation = ComplexHeatmap::upset_left_annotation(
      input,
      add_numbers = TRUE
    )
  )
  draw(tmp_p)
}

extract_upset_id <- function(
  upset = NULL,
  comps = NULL
) {
  set_comp <- ComplexHeatmap::set_name(upset)
  paste0(as.integer(set_comp %in% comps), collapse = "")
}

extract_upset_comps <- function(
  upset = NULL
) {
  combs <- ComplexHeatmap::comb_name(upset, readable = FALSE)
  purrr::map(
    .x = magrittr::set_names(combs, combs),
    .f = ~ ComplexHeatmap::extract_comb(upset, .x)
  )
}

mol_similarity <- function(
  query_smiles = NULL,
  target_smiles = NULL,
  kekulise = TRUE,
  omit_nulls = TRUE,
  fingerprint = "circular",
  circular_type = "ECFP6",
  method = "tanimoto"
) {

  query_mol <- rcdk::parse.smiles(
    smiles = query_smiles,
    kekulise = kekulise,
    omit.nulls = omit_nulls
  )[[1]]

  target_mols <- rcdk::parse.smiles(
    smiles = target_smiles,
    kekulise = kekulise,
    omit.nulls = omit_nulls
  )

  query_fp <- rcdk::get.fingerprint(
    molecule = query_mol,
    type = fingerprint,
    circular.type = circular_type
  )

  target_fps <- lapply(
    X = target_mols,
    FUN = function(x) {
      rcdk::get.fingerprint(
        molecule = x,
        type = fingerprint,
        circular.type = circular_type
      )
    }
  )

  sims <- target_fps %>%
    purrr::map2(
      .x = .,
      .y = names(.),
      .f = ~ tibble::tibble(
        feature = .y,
        sim = fingerprint::distance(
          fp1 = .x,
          fp2 = query_fp,
          method = "tanimoto"
        )
      )
    ) %>%
    purrr::list_rbind() %>%
    dplyr::arrange(dplyr::desc(sim))

  return(sims)
}

# TODO
# check later -> hclust of similarity
# fps <- lapply(mols, get.fingerprint, type='circular')
# fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')
# fp.dist <- 1 - fp.sim
# cls <- hclust(as.dist(fp.dist))
# plot(cls, labels=FALSE)

run_biotransformer <- function(
  bt_dir = "biotransformer3.0jar",
  smiles = NULL,
  b_type = "superbio",
  k_task = "pred", # pred for prediction, or cid for compound identification
  output_file = "prediction",
  results_path = snakemake@params$output
) {
  old_wd <- getwd()
  new_wd <- normalizePath(bt_dir, mustWork = FALSE)
  biot_output_loc <- normalizePath(
    file.path(results_path, "tables"),
    mustWork = FALSE
  )
  clean_nm <- gsub("\\..*$", "", output_file)

  setwd(new_wd)
  biot_output <- system2(
    command = "java",
    args = c(
      "-jar", "BioTransformer3.0_20230525.jar",
      "-k", k_task,
      "-b", b_type,
      "-ismi", smiles,
      "-ocsv", paste0(biot_output_loc, "/", clean_nm, ".csv"),
      "-a"
    ),
    stdout = TRUE,
    stderr = TRUE
  )
  cat(biot_output, sep = "\n")
  setwd(old_wd)
}

script_header <- function() {
  rule_name <- sub("^[^.]*\\.", "", snakemake@rule)
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cli::cli_h1(
    "{.file {rule_name}} {.emph {timestamp}}"
  )
}

get_chr_data <- function(chrs, feature_name, meta) {
  feature_idxs <- match(feature_name, rownames(xcms::featureDefinitions(chrs)))
  purrr::map_dfr(
    seq_along(feature_idxs),
    .f = \(i) {
      purrr::map_dfr(
        seq_len(ncol(chrs)),
        .f = \(j) {
          chr <- chrs[feature_idxs[i], j]
          tibble::tibble(
            rtime = chr@rtime,
            intensity = chr@intensity,
            mzmin = chr@mz[1],
            mzmax = chr@mz[2],
            sample = colnames(chrs)[j],
            feature = feature_name[i]
          )
        }
      )
    }
  )
}

plot_feature <- function(
  feature_chrom,
  feature,
  meta,
  limma_results,
  method = "sum",
  value = "into",
  filled = TRUE,
  missing = 0,
  ms_level = 1L,
  save_loc = NULL,
  device = "pdf",
  overwrite = FALSE
) {

  if (!is.null(save_loc)) {
    file_nm <- paste0(save_loc, feature, ".", device)
    if (interactive() && file.exists(file_nm) && !overwrite) {
      cli::cli_alert_info("{.path {file_nm}} already exists, skipping")
      return(invisible(NULL))
    }
  }

  tmp <- get_chr_data(feature_chrom, feature, meta) %>%
    dplyr::left_join(
      x = .,
      y = dplyr::select(meta, sample, group),
      by = "sample"
    ) %>%
    dplyr::mutate(intensity = tidyr::replace_na(intensity, 0)) %>%
    dplyr::arrange(group, sample) %>%
    dplyr::mutate(sample = forcats::fct_inorder(sample))

  p1 <- tmp %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = rtime,
        y = intensity,
        group = sample,
        color = group
      )
    ) +
    ggplot2::geom_line(linewidth = 1.5) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.title = ggplot2::element_blank()) +
    ggplot2::labs(x = "Retention time (s)")

  p2 <- tmp %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = rtime,
        y = sample,
        height = intensity,
        color = group
      )
    ) +
    ggridges::geom_ridgeline(scale = 3 / max(tmp$intensity), fill = NA) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      legend.position = "none"
    ) +
    ggplot2::labs(x = "Retention time (s)")

  boxplot_data <- xcms::featureValues(
    feature_chrom,
    method = method,
    value = value,
    intensity = value,
    filled = filled,
    missing = missing,
    ms_level = ms_level
  ) %>%
    tibble::as_tibble(., rownames = "feature") %>%
    dplyr::filter(feature == {{ feature }}) %>%
    tidyr::pivot_longer(cols = dplyr::contains(".mzML")) %>%
    dplyr::left_join(
      x = .,
      y = meta,
      by = c("name" = "sample")
    )

  signif_data <- limma_results %>%
    dplyr::filter(feature == {{ feature }}) %>%
    dplyr::select(feature, adj.P.Val, contrast) %>%
    dplyr::mutate(
      group1 = stringr::str_split_i(contrast, "-", 1),
      group2 = stringr::str_split_i(contrast, "-", 2),
    ) %>%
    dplyr::select(-contrast) %>%
    dplyr::filter(feature %in% unique(boxplot_data$feature)) %>%
    rstatix::add_significance(p.col = "adj.P.Val") %>%
    find_y_position(
      test_df = .,
      df = boxplot_data,
      formula = "value ~ group",
      fun_data = "max"
    )

  signif_only <- signif_data %>%
    dplyr::filter(adj.P.Val.signif != "ns")

  p3 <- boxplot_data %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = group,
        y = value
      )
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(fill = group),
      show.legend = FALSE,
      outliers = FALSE
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = group),
      position = ggplot2::position_jitter(width = 0.3),
      size = 3,
      show.legend = FALSE,
      color = "black",
      pch = 21
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank()
    ) +
    ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(c(0.1, 0.1))
    ) +
    ggplot2::labs(
      y = "Peak area"
    ) +
    {
      if (nrow(signif_only) > 0) {
        ggpubr::geom_bracket(
          data = signif_only,
          ggplot2::aes(
            xmin = group1,
            xmax = group2,
            label = adj.P.Val.signif,
            y.position = y.pos * 1.1
          ),
          step.increase = 0.15,
          vjust = 0.1
        )
      } else {
        NULL
      }
    }

  feat_defs <- xcms::featureDefinitions(feature_chrom) %>%
    tibble::as_tibble(., rownames = "feature") %>%
    dplyr::filter(feature == {{ feature }})

  combined <- (patchwork::free(p2) | (p1 / p3)) +
    patchwork::plot_layout(
      axes = "collect",
      guides = "collect"
    ) +
    patchwork::plot_annotation(
      title = paste0(
        feature, ", ",
        "m/z: ", round(feat_defs$mzmed, 2), ", ",
        "retention time: ", round(feat_defs$rtmed, 2)
      )
    )

  half <- (p1 / p3) +
    patchwork::plot_layout(
      axes = "collect",
      guides = "collect"
    ) +
    patchwork::plot_annotation(
      title = paste0(
        feature, ", ",
        "m/z: ", round(feat_defs$mzmed, 2), ", ",
        "retention time: ", round(feat_defs$rtmed, 2)
      )
    )

  if (!is.null(save_loc)) {
    ggplot2::ggsave(
      filename = file_nm,
      plot = p_final,
      device = device,
      height = 7,
      width = 10,
      units = "in"
    )
    cli::cli_alert_success("Saved {.path {file_nm}}")
  }

  list(
    full = combined,
    half = half,
    overlay = p1,
    ridgeline = p2,
    boxplot = p3,
    boxplot_data = boxplot_data,
    signif_data = signif_data
  )
}
