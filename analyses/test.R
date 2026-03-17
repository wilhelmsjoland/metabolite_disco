
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
    BPPARAM = BiocParallel::bpparam(),
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
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1L,
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
