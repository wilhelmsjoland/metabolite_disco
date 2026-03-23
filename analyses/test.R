
# getMonoisotopic(getMolecule("C15H10O5"))

# library(MetaboCoreUtils)
# library(MsCoreUtils)

# test <- getMonoisotopic(getMolecule("C15H10O5"))

# unname(mass2mz(test, "[M-H]-"))[,1]

# unname(mass2mz(getMass(getMolecule("C26H28O14")), "[M-H]-"))[, 1]


subset_matched_diffs <- pred_biot(
    data = possible_adducts_signif,
    biotransf_data = bio_transf2,
    tolerance_ppm = 5,
    features_of_interest = pot_glycosides,
    parallel = TRUE
  ) 
  # %>%
  #   dplyr::mutate(
  #     pair = purrr::map2(
  #       .x = feat1,
  #       .y = feat2,
  #       .f = c
  #     ),
  #     obs_diff = abs(obs_delta_mass - delta_mass),
  #     ppm_diff = num_to_ppm(mz = delta_mass, diff = obs_diff),
  #     ppm_mz = num_to_ppm(mz = pmax(mz1, mz2), diff = obs_diff)
  #   )

getMass(getMolecule("C21H20O10"))

ms_exp <- MsExperiment::readMsExperiment(
  spectraFiles = meta$path,
  sampleData = meta
)

xcms::findChromPeaks(
  object = ms_exp,
  BPPARAM = SerialParam(),
  return.type = "XCMSnExp",
  ms_level = 1L,
  # chunkSize = 10,
  param = xcms::CentWaveParam(
    ppm = 25,
    peakwidth = c(min_peak_width, max_peak_width),
    snthresh = 10,
    prefilter = c(3, 1000), # k pks (left) over intens (right) # c(4, 1000)
    mzCenterFun = "wMean",
    integrate = 2,
    mzdiff = 0.01, # 0.001
    fitgauss = TRUE,
    noise = 1000,
    verboseColumns = TRUE,
    roiList = list(),
    firstBaselineCheck = TRUE,
    roiScales = numeric(),
    extendLengthMSW = FALSE,
    verboseBetaColumns = TRUE
  )
)

is_chr <- xcms::chromatogram(
    BPPARAM = BiocParallel::bpparam(),
    chunkSize = 10,
    object = ms_exp,
    mz = mz_range,
    aggregationFun = "sum",
    
  )


is_chr <- xcms::chromatogram(
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 2L,
    object = ms_exp,
    mz = mz_range,
    aggregationFun = "sum",
)


 xcms::adjustRtime(
    object = xchr5,
    chunkSize = 10,
    BPPARAM = BiocParallel::bpparam(),
    param = xcms::PeakGroupsParam(
      minFraction = 0.9, # 0.8
      extraPeaks = 0, # 0
      smooth = "loess",
      peakGroupsMatrix = pgm_filt,
      span = 0.6, # 0.6 # 0.8
      family = "gaussian",
      # peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
      subset = integer(),
      subsetAdjust = c("average", "previous")
    )
  )

xcms::fillChromPeaks(
    object = xchr7,
    BPPARAM = SnowParam(workers = 4, type = "SOCK"),
    chunkSize = 4,
    param = xcms::ChromPeakAreaParam()
  )

xcms::fillChromPeaks(
    object = xchr7,
    BPPARAM = SerialParam(),
    chunkSize = 1L,
    param = xcms::ChromPeakAreaParam()
  )

xcms::findChromPeaks(
    object = is_chr_wide,
    BPPARAM = BiocParallel::bpparam(),
    chunkSize = 10,
    param = xcms::CentWaveParam(
      ppm = 25,
      peakwidth = c(2, 20),
      prefilter = c(1, 1),
      snthresh = 10, # 10
      mzCenterFun = "wMean",
      mzdiff = 0.001,
      integrate = 2,
      noise = 1000,
      verboseBetaColumns = TRUE
    ),
    ms_level = 1
  )

xcms::findChromPeaks(
    object = is_chr_wide,
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1,
    param = xcms::CentWaveParam(
      ppm = 25,
      peakwidth = c(2, 20),
      prefilter = c(1, 1),
      snthresh = 10, # 10
      mzCenterFun = "wMean",
      mzdiff = 0.001,
      integrate = 2,
      noise = 1000,
      verboseBetaColumns = TRUE
    ),
    ms_level = 1
  )

xcms::findChromPeaks(
    object = ms_exp,
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = snakemake@params$cores,
    mslevel = 1L,
    param = xcms::CentWaveParam(
      ppm = snakemake@params$ppm_global,
      peakwidth = c(min_peak_width, max_peak_width),
      snthresh = snakemake@params$sn_threshold,
      prefilter = c(3, 1000), # k pks (left) over intens (right) # c(4, 1000)
      mzCenterFun = "wMean",
      integrate = 2,
      mzdiff = snakemake@params$mzdiff, # 0.001
      fitgauss = FALSE, # TRUE
      noise = 1000,
      verboseColumns = TRUE,
      roiList = list(),
      firstBaselineCheck = TRUE,
      roiScales = numeric(),
      extendLengthMSW = FALSE,
      verboseBetaColumns = TRUE
    )
  )

xcms::findChromPeaks(
    object = ms_exp,
    BPPARAM = BiocParallel::bpparam(),
    chunkSize = 10L,
    mslevel = 1L,
    param = xcms::CentWaveParam(
      ppm = 25,
      peakwidth = c(min_peak_width, max_peak_width),
      snthresh = 10,
      prefilter = c(4, 1000), # k pks (left) over intens (right) # c(4, 1000)
      mzCenterFun = "wMean",
      integrate = 1,
      mzdiff = 0.01, # 0.001
      fitgauss = FALSE, # TRUE
      noise = 1000,
      verboseColumns = FALSE,
      roiList = list(),
      firstBaselineCheck = TRUE,
      roiScales = numeric(),
      extendLengthMSW = FALSE,
      verboseBetaColumns = TRUE
  )
)


test <- rcdk::parse.smiles(glycone_pairs_metadata$aglycone_SMILES)

test

# TODO
# FIX SO EACH SCRIPT CHOOSES ITS OWN PARAMETERS

test <- run_biotransformer(
  bt_dir = snakemake@config$biot_dir,
  smiles = snakemake@config$smiles,
  b_type = "superbio", #superbio
  k_task = "pred",
  output_file = "prediction",
  results_path = "C:/Users/wilhelm/Documents/MEGA/01_juniper/01_arbete/01_projekt/03_psm/output/test/"
)

test2 <- readr::read_csv(
  prediction_path
)

# 356 rows in allHuman - time 105475 with 3 iterations

# - time 284886



test <- xcms::chromPeaks(xchr) %>% 
  tibble::as_tibble(., rownames = "peaks")

test %>%
  dplyr::mutate(peaks2 = as.numeric(gsub("[A-Za-z]", "", peaks))) %>%
  dplyr::mutate(
    param_filter = dplyr::if_else(
      sn > 1000,
      "low",
      "high"
    )
  ) %>%
  ggplot(
    aes(
      x = sn,
      fill = param_filter
    )
  ) +
  geom_histogram() +
  scale_x_continuous(transform = "log10", labels = scales::comma)

test %>%
  ggplot(
    aes(
      x = maxo,
    )
  ) +
  geom_histogram() +
  scale_x_continuous(transform = "log10", labels = scales::comma)


test2 <- xcms::featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(feature %in% all_sig_diff) %>%
  unnest(cols = "peakidx")


test %>%
  dplyr::mutate(peaks2 = as.numeric(gsub("[A-Za-z]", "", peaks))) %>%
  dplyr::filter(peaks2 %in% test2$peakidx) %>%
  ggplot(
    aes(
      x = sn
    )
  ) +
  geom_histogram() +
  scale_x_continuous(transform = "log10", labels = scales::comma)

test3 <- xcms::findChromPeaks(
    object = ms_exp[1],
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1L,
    mslevel = 1L,
    param = xcms::CentWaveParam(
      ppm = 25,
      peakwidth = c(min_peak_width, max_peak_width),
      snthresh = 10,
      prefilter = c(3, 1000), # k pks (left) over intens (right) # c(4, 1000)
      mzCenterFun = "wMean",
      integrate = 2,
      mzdiff = 0.01, # 0.001
      fitgauss = FALSE, # TRUE
      noise = 1000,
      verboseColumns = FALSE,
      roiList = list(),
      firstBaselineCheck = TRUE,
      roiScales = numeric(),
      extendLengthMSW = FALSE,
      verboseBetaColumns = TRUE # TRUE
    )
  )

pks <- as.data.frame(chromPeaks(xchr))

library(MetaboCoreUtils)

xchr_spectra <- xcms::chromPeakSpectra(
  object = xchr,
  return.type = "Spectra"
)

deisotopeSpectra(
  x,
  substDefinition = isotopicSubstitutionMatrix("HMDB_NEUTRAL"),
  tolerance = 0,
  ppm = 20,
  charge = 1
)


iso_groups <- isotopologues(
  x = xchr_spectra,
  ppm = 25
)

test %>%
  dplyr::arrange(maxo) %>%
  dplyr::arrange(sn)
  dplyr::filter(maxo > 1000)


predicted_feats_path <- file.path(
  snakemake@config$output,
  "tables",
  "predicted_annotated_feats.csv"
)


test <- paste0(
  "/Users/wilhelm/Documents/from_ssd/aglycone_release_100um_24h/output",
  "/afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
  "snakemake_objects", "19_molecular_similarity.rds"
)

test
