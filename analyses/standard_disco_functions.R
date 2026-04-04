read_mzml_memory_spectra <- function(
  path = NULL,
  file_name = basename(path)
) {
  ms_file <- mzR::openMSfile(path)
  on.exit(mzR::close(ms_file))

  header <- mzR::header(ms_file)
  is_empty_scan <- is.na(header$peaksCount) | header$peaksCount < 1
  peak_data <- vector("list", nrow(header))
  skipped_nonempty <- list()

  for (i in seq_len(nrow(header))) {
    if (is_empty_scan[i]) {
      peak_data[[i]] <- list(mz = numeric(), intensity = numeric())
      next
    }

    peak_matrix <- tryCatch(
      mzR::peaks(ms_file, scans = i),
      error = function(e) e
    )

    if (inherits(peak_matrix, "error")) {
      skipped_nonempty[[length(skipped_nonempty) + 1]] <- tibble::tibble(
        file_name = file_name,
        path = path,
        scan_index = i,
        acquisitionNum = header$acquisitionNum[i],
        retentionTime = header$retentionTime[i],
        peaksCount = header$peaksCount[i],
        msLevel = header$msLevel[i],
        error = conditionMessage(peak_matrix)
      )
      peak_data[[i]] <- list(mz = numeric(), intensity = numeric())
      next
    }

    if (is.null(peak_matrix) || !length(peak_matrix) || is.null(dim(peak_matrix))) {
      skipped_nonempty[[length(skipped_nonempty) + 1]] <- tibble::tibble(
        file_name = file_name,
        path = path,
        scan_index = i,
        acquisitionNum = header$acquisitionNum[i],
        retentionTime = header$retentionTime[i],
        peaksCount = header$peaksCount[i],
        msLevel = header$msLevel[i],
        error = "NULL/invalid peak matrix"
      )
      peak_data[[i]] <- list(mz = numeric(), intensity = numeric())
      next
    }

    peak_data[[i]] <- list(
      mz = as.numeric(peak_matrix[, 1]),
      intensity = as.numeric(peak_matrix[, 2])
    )
  }

  spectra_df <- data.frame(
    msLevel = as.integer(header$msLevel),
    rtime = as.numeric(header$retentionTime),
    acquisitionNum = as.integer(header$acquisitionNum),
    scanIndex = seq_len(nrow(header)),
    dataOrigin = rep(path, nrow(header)),
    centroided = as.logical(header$centroided),
    smoothed = rep(NA, nrow(header)),
    polarity = as.integer(header$polarity),
    precScanNum = as.integer(header$precursorScanNum),
    precursorMz = as.numeric(header$precursorMZ),
    precursorIntensity = as.numeric(header$precursorIntensity),
    precursorCharge = as.integer(header$precursorCharge),
    collisionEnergy = as.numeric(header$collisionEnergy),
    isolationWindowLowerMz = as.numeric(
      header$isolationWindowTargetMZ - header$isolationWindowLowerOffset
    ),
    isolationWindowTargetMz = as.numeric(header$isolationWindowTargetMZ),
    isolationWindowUpperMz = as.numeric(
      header$isolationWindowTargetMZ + header$isolationWindowUpperOffset
    ),
    peaksCount = as.integer(header$peaksCount),
    totIonCurrent = as.numeric(header$totIonCurrent),
    basePeakMZ = as.numeric(header$basePeakMZ),
    basePeakIntensity = as.numeric(header$basePeakIntensity),
    electronBeamEnergy = as.numeric(header$electronBeamEnergy),
    ionisationEnergy = as.numeric(header$ionisationEnergy),
    lowMZ = as.numeric(header$lowMZ),
    highMZ = as.numeric(header$highMZ),
    mergedScan = as.integer(header$mergedScan),
    mergedResultScanNum = as.integer(header$mergedResultScanNum),
    mergedResultStartScanNum = as.integer(header$mergedResultStartScanNum),
    mergedResultEndScanNum = as.integer(header$mergedResultEndScanNum),
    injectionTime = as.numeric(header$injectionTime),
    filterString = as.character(header$filterString),
    spectrumId = as.character(header$spectrumId),
    ionMobilityDriftTime = as.numeric(header$ionMobilityDriftTime),
    scanWindowLowerLimit = as.numeric(header$scanWindowLowerLimit),
    scanWindowUpperLimit = as.numeric(header$scanWindowUpperLimit),
    file_name = rep(file_name, nrow(header)),
    stringsAsFactors = FALSE
  )
  spectra_df$mz <- purrr::map(peak_data, "mz")
  spectra_df$intensity <- purrr::map(peak_data, "intensity")

  backend <- Spectra::backendInitialize(
    Spectra::MsBackendMemory(),
    data = spectra_df,
    peaksVariables = c("mz", "intensity")
  )

  list(
    spectra = Spectra::Spectra(backend),
    empty_scan_count = sum(is_empty_scan, na.rm = TRUE),
    skipped_nonempty = dplyr::bind_rows(skipped_nonempty)
  )
}

import_msexperiment_memory <- function(
  stds = NULL,
  failure_log_path = NULL,
  stop_on_failed_nonempty = FALSE
) {
  stds_df <- as.data.frame(stds)

  if (is.null(rownames(stds_df)) || !length(rownames(stds_df))) {
    rownames(stds_df) <- stds_df$file_name
  }
  if (!"file_name" %in% colnames(stds_df)) {
    stds_df$file_name <- basename(stds_df$path)
  }

  spectra_imports <- purrr::map2(
    .x = stds_df$path,
    .y = stds_df$file_name,
    .f = read_mzml_memory_spectra
  )
  spectra_obj <- do.call(c, purrr::map(spectra_imports, "spectra"))
  skipped_nonempty <- purrr::map(spectra_imports, "skipped_nonempty") %>%
    dplyr::bind_rows()

  if (!is.null(failure_log_path)) {
    dir.create(dirname(failure_log_path), recursive = TRUE, showWarnings = FALSE)
    if (nrow(skipped_nonempty)) {
      readr::write_csv(skipped_nonempty, failure_log_path)
    } else if (file.exists(failure_log_path)) {
      unlink(failure_log_path)
    }
  }

  if (nrow(skipped_nonempty)) {
    warn_msg <- paste0(
      "Skipped ", nrow(skipped_nonempty),
      " non-empty scan(s) during mzML import"
    )
    if (!is.null(failure_log_path)) {
      warn_msg <- paste0(warn_msg, ". Details written to ", failure_log_path)
    }

    if (isTRUE(stop_on_failed_nonempty)) {
      stop(warn_msg, call. = FALSE)
    }
    warning(warn_msg, call. = FALSE)
  }

  experiment <- MsExperiment::MsExperiment(
    spectra = spectra_obj,
    sampleData = S4Vectors::DataFrame(
      stds_df,
      row.names = rownames(stds_df)
    )
  )

  MsExperiment::linkSampleData(
    experiment,
    with = "sampleData.file_name=spectra.file_name"
  )
}

get_target_adduct_priority <- function() {
  c(
    "[M-H]-" = 1,
    "[M+Cl]-" = 2,
    "[M+CHO2]-" = 3,
    "[M+C2H3O2]-" = 4,
    "[M+Na-2H]-" = 5,
    "[M+K-2H]-" = 6,
    "[2M-H]-" = 7,
    "[M-H+HCOONa]-" = 8
  )
}

get_target_peak_param <- function(ppm = 25) {
  xcms::CentWaveParam(
    ppm = ppm,
    peakwidth = c(1, 80),
    prefilter = c(1, 1),
    snthresh = 1,
    mzCenterFun = "wMean",
    mzdiff = -0.001,
    integrate = 2,
    noise = 0,
    fitgauss = FALSE,
    verboseBetaColumns = FALSE
  )
}

safe_median <- function(x) {
  if (length(x) && any(is.finite(x))) {
    stats::median(x, na.rm = TRUE)
  } else {
    0
  }
}

safe_max <- function(x) {
  if (length(x) && any(is.finite(x))) {
    max(x, na.rm = TRUE)
  } else {
    0
  }
}

safe_q90 <- function(x) {
  if (length(x) && any(is.finite(x))) {
    as.numeric(stats::quantile(x, 0.9, na.rm = TRUE, names = FALSE))
  } else {
    0
  }
}

safe_mad <- function(x) {
  if (sum(is.finite(x)) >= 2) {
    stats::mad(x, na.rm = TRUE)
  } else {
    NA_real_
  }
}

summarize_eic_window <- function(
  chrom_trace = NULL,
  rt_window = NULL
) {
  trace_rt <- chrom_trace@rtime
  trace_intensity <- tidyr::replace_na(chrom_trace@intensity, replace = 0)

  if (!length(trace_intensity)) {
    return(tibble::tibble(
      apex_rt = NA_real_,
      apex_intensity = 0,
      window_area = 0
    ))
  }

  if (!is.null(rt_window)) {
    keep_idx <- trace_rt >= rt_window[1] & trace_rt <= rt_window[2]
    trace_rt <- trace_rt[keep_idx]
    trace_intensity <- trace_intensity[keep_idx]
  }

  if (!length(trace_intensity)) {
    return(tibble::tibble(
      apex_rt = NA_real_,
      apex_intensity = 0,
      window_area = 0
    ))
  }

  apex_idx <- which.max(trace_intensity)[1]
  apex_intensity <- safe_max(trace_intensity)
  apex_rt <- if (apex_intensity > 0) trace_rt[apex_idx] else NA_real_

  tibble::tibble(
    apex_rt = apex_rt,
    apex_intensity = apex_intensity,
    window_area = sum(trace_intensity, na.rm = TRUE)
  )
}

find_target_chrom_peaks <- function(
  target_chrom = NULL,
  positive_samples = NULL,
  mz_theory = NULL,
  ppm = 25,
  peak_param = get_target_peak_param(ppm)
) {
  mz_window <- ppm_to_num(mz_theory, ppm)
  pos_idx <- which(colnames(target_chrom) %in% positive_samples)

  if (!length(pos_idx)) {
    return(tibble::tibble())
  }

  pos_chrom <- target_chrom[, pos_idx]

  peak_obj <- tryCatch(
    suppressWarnings(
      xcms::findChromPeaks(
        object = pos_chrom,
        param = peak_param
      )
    ),
    error = function(e) NULL
  )

  if (is.null(peak_obj) || nrow(xcms::chromPeaks(peak_obj)) == 0) {
    return(tibble::tibble())
  }

  sample_map <- tibble::tibble(
    column = seq_along(colnames(pos_chrom)),
    filename = colnames(pos_chrom)
  )

  xcms::chromPeaks(peak_obj) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("peak_id") %>%
    tibble::as_tibble() %>%
    dplyr::left_join(sample_map, by = "column") %>%
    dplyr::mutate(
      mz_theory = mz_theory,
      mz_window = mz_window
    )
}

cluster_target_peaks <- function(
  peak_table = NULL,
  rt_gap = 2
) {
  if (is.null(peak_table) || nrow(peak_table) == 0) {
    return(tibble::tibble())
  }

  peak_table <- peak_table %>%
    dplyr::arrange(rtmin, rtmax, rt)

  feature_index <- integer(nrow(peak_table))
  feature_index[1] <- 1L
  current_feature <- 1L
  current_rtmax <- peak_table$rtmax[1]

  if (nrow(peak_table) > 1) {
    for (i in 2:nrow(peak_table)) {
      if (peak_table$rtmin[i] <= current_rtmax + rt_gap) {
        feature_index[i] <- current_feature
        current_rtmax <- max(current_rtmax, peak_table$rtmax[i], na.rm = TRUE)
      } else {
        current_feature <- current_feature + 1L
        feature_index[i] <- current_feature
        current_rtmax <- peak_table$rtmax[i]
      }
    }
  }

  peak_table %>%
    dplyr::mutate(feature_index = feature_index)
}

make_empty_target_feature_row <- function(
  group = NULL,
  formula = NULL,
  adduct = NULL,
  mz_theory = NULL,
  adduct_rank = NULL,
  positive_samples = NULL
) {
  tibble::tibble(
    group = group,
    formula = formula,
    adduct = adduct,
    adduct_rank = adduct_rank,
    mz_theory = mz_theory,
    feature_id = NA_character_,
    feature_found = FALSE,
    n_expected_samples = length(positive_samples),
    n_pos_samples_detected = 0L,
    prop_pos_samples_detected = 0,
    npeaks = 0L,
    mzmed = NA_real_,
    mzmin = NA_real_,
    mzmax = NA_real_,
    rtmed = NA_real_,
    rtmin = NA_real_,
    rtmax = NA_real_,
    peak_rt_mad = NA_real_,
    median_peak_apex = 0,
    max_peak_apex = 0,
    median_peak_area = 0,
    sum_peak_area = 0,
    median_peak_sn = 0,
    positive_peak_files = NA_character_,
    signal_cutoff_area = 0,
    n_pos_signal_area = 0L,
    prop_pos_signal_area = 0,
    median_pos_window_area = 0,
    median_neg_window_area = 0,
    enrichment_area = 1,
    pos_samples_signal_area = NA_character_,
    signal_cutoff_apex = 0,
    n_pos_signal_apex = 0L,
    prop_pos_signal_apex = 0,
    median_pos_window_apex = 0,
    median_neg_window_apex = 0,
    enrichment_apex = 1,
    pos_samples_signal_apex = NA_character_
  )
}

summarize_target_features <- function(
  peak_table = NULL,
  target_chrom = NULL,
  positive_samples = NULL,
  group = NULL,
  formula = NULL,
  adduct = NULL,
  mz_theory = NULL,
  adduct_rank = NULL,
  rt_gap = 2,
  BPPARAM = BiocParallel::MulticoreParam(opt$cores),
  chunkSize = opt$cores
) {
  if (is.null(peak_table) || nrow(peak_table) == 0) {
    return(
      make_empty_target_feature_row(
        group = group,
        formula = formula,
        adduct = adduct,
        mz_theory = mz_theory,
        adduct_rank = adduct_rank,
        positive_samples = positive_samples
      )
    )
  }

  clustered_peaks <- cluster_target_peaks(
    peak_table = peak_table,
    rt_gap = rt_gap
  )

  purrr::map_dfr(
    .x = sort(unique(clustered_peaks$feature_index)),
    .f = function(i) {
      tmp_peaks <- clustered_peaks %>%
        dplyr::filter(feature_index == i)

      feature_row <- tibble::tibble(
        group = group,
        formula = formula,
        adduct = adduct,
        adduct_rank = adduct_rank,
        mz_theory = mz_theory,
        feature_id = paste0("SF", stringr::str_pad(i, 3, pad = "0")),
        feature_found = TRUE,
        n_expected_samples = length(positive_samples),
        n_pos_samples_detected = dplyr::n_distinct(tmp_peaks$filename),
        prop_pos_samples_detected =
          dplyr::n_distinct(tmp_peaks$filename) / max(1, length(positive_samples)),
        npeaks = nrow(tmp_peaks),
        mzmed = stats::median(tmp_peaks$mz, na.rm = TRUE),
        mzmin = min(tmp_peaks$mzmin, na.rm = TRUE),
        mzmax = max(tmp_peaks$mzmax, na.rm = TRUE),
        rtmed = stats::median(tmp_peaks$rt, na.rm = TRUE),
        rtmin = min(tmp_peaks$rtmin, na.rm = TRUE),
        rtmax = max(tmp_peaks$rtmax, na.rm = TRUE),
        peak_rt_mad = safe_mad(tmp_peaks$rt),
        median_peak_apex = safe_median(tmp_peaks$maxo),
        max_peak_apex = safe_max(tmp_peaks$maxo),
        median_peak_area = safe_median(tmp_peaks$into),
        sum_peak_area = sum(tmp_peaks$into, na.rm = TRUE),
        median_peak_sn = safe_median(tmp_peaks$sn),
        positive_peak_files = paste(sort(unique(tmp_peaks$filename)), collapse = "; ")
      )

      window_stats <- purrr::map_dfr(
        .x = seq_len(ncol(target_chrom)),
        .f = function(j) {
          summarize_eic_window(
            chrom_trace = target_chrom[1, j],
            rt_window = c(feature_row$rtmin, feature_row$rtmax)
          ) %>%
            dplyr::mutate(sample = colnames(target_chrom)[j])
        }
      ) %>%
        dplyr::mutate(is_pos = sample %in% positive_samples)

      pos_stats <- window_stats %>%
        dplyr::filter(is_pos)
      neg_stats <- window_stats %>%
        dplyr::filter(!is_pos)

      neg_area_med <- safe_median(neg_stats$window_area)
      neg_apex_med <- safe_median(neg_stats$apex_intensity)

      signal_cutoff_area <- max(
        1000,
        safe_q90(neg_stats$window_area),
        neg_area_med * 3
      )
      signal_cutoff_apex <- max(
        1000,
        safe_q90(neg_stats$apex_intensity),
        neg_apex_med * 3
      )

      pos_signal_area <- pos_stats %>%
        dplyr::filter(window_area > signal_cutoff_area)
      pos_signal_apex <- pos_stats %>%
        dplyr::filter(apex_intensity > signal_cutoff_apex)

      feature_row %>%
        dplyr::mutate(
          signal_cutoff_area = signal_cutoff_area,
          n_pos_signal_area = nrow(pos_signal_area),
          prop_pos_signal_area = nrow(pos_signal_area) / max(1, nrow(pos_stats)),
          median_pos_window_area = safe_median(pos_stats$window_area),
          median_neg_window_area = neg_area_med,
          enrichment_area =
            (median_pos_window_area + 1) / (median_neg_window_area + 1),
          pos_samples_signal_area = paste(pos_signal_area$sample, collapse = "; "),
          signal_cutoff_apex = signal_cutoff_apex,
          n_pos_signal_apex = nrow(pos_signal_apex),
          prop_pos_signal_apex = nrow(pos_signal_apex) / max(1, nrow(pos_stats)),
          median_pos_window_apex = safe_median(pos_stats$apex_intensity),
          median_neg_window_apex = neg_apex_med,
          enrichment_apex =
            (median_pos_window_apex + 1) / (median_neg_window_apex + 1),
          pos_samples_signal_apex = paste(pos_signal_apex$sample, collapse = "; ")
        )
    }
  ) %>%
    dplyr::mutate(
      feature_id = paste(group, adduct, feature_id, sep = "__")
    )
}

build_targeted_feature_table <- function(
  experiment = NULL,
  sample_targets = NULL,
  ppm = 25,
  BPPARAM = BiocParallel::MulticoreParam(opt$cores),
  chunkSize = opt$cores,
  rt_gap = 2,
  adduct_priority = get_target_adduct_priority()
) {
  unique_targets <- sample_targets %>%
    dplyr::distinct(group, formula) %>%
    dplyr::arrange(group)

  purrr::map_dfr(
    .x = seq_len(nrow(unique_targets)),
    .f = function(i) {
      target_group <- unique_targets$group[i]
      target_formula <- unique_targets$formula[i]
      positive_samples <- sample_targets %>%
        dplyr::filter(group == target_group) %>%
        dplyr::distinct(filename) %>%
        dplyr::pull(filename)

      target_mass <- Rdisop::getMonoisotopic(
        Rdisop::getMolecule(target_formula)
      )

      target_mz <- MetaboCoreUtils::mass2mz(
        x = setNames(target_mass, target_formula),
        adduct = names(adduct_priority)
      ) %>%
        as.numeric()

      purrr::map_dfr(
        .x = seq_along(target_mz),
        .f = function(j) {
          tmp_adduct <- names(adduct_priority)[j]
          tmp_mz <- target_mz[j]
          tmp_mz_window <- ppm_to_num(tmp_mz, ppm)
          tmp_chrom <- suppressMessages(
            xcms::chromatogram(
              object = experiment,
              mz = c(tmp_mz - tmp_mz_window, tmp_mz + tmp_mz_window),
              aggregationFun = "sum",
              BPPARAM = BPPARAM,
              chunkSize = chunkSize
            )
          )
          tmp_peaks <- find_target_chrom_peaks(
            target_chrom = tmp_chrom,
            positive_samples = positive_samples,
            mz_theory = tmp_mz,
            ppm = ppm
          )

          summarize_target_features(
            peak_table = tmp_peaks,
            target_chrom = tmp_chrom,
            positive_samples = positive_samples,
            group = target_group,
            formula = target_formula,
            adduct = tmp_adduct,
            mz_theory = tmp_mz,
            adduct_rank = unname(adduct_priority[tmp_adduct]),
            rt_gap = rt_gap,
            BPPARAM = BPPARAM,
            chunkSize = chunkSize
          )
        }
      )
    }
  )
}
