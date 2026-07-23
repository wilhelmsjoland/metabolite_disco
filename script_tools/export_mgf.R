################################################################################
# Export mgfs ------------------------------------------------------------------
################################################################################

# Export path
mgf_export_path <- file.path(input_path, "sirius")

# Extract MS1 spectra at feature apexes
feat_spectra <- xcms::featureSpectra(
  xchr9,
  msLevel = 1L,
  method = "closest_rt", # Makes sense for MS1 but not for MS2
  BPPARAM = BiocParallel::SerialParam(),
  features = feature_levels
)

feat_spectra$precursorCharge <- -1L

feat_export <- Spectra::selectSpectraVariables(
  feat_spectra,
  c(
    "scanIndex",
    "mz",
    "intensity",
    "rtime",
    "precursorMz",
    "precursorCharge",
    "msLevel"
  )
)

map <- c(
  feature_id = "TITLE",
  Spectra::spectraVariableMapping(MsBackendMgf::MsBackendMgf())
)

dir.create(
  mgf_export_path,
  recursive = TRUE,
  showWarnings = FALSE
)

# Export into the output folder
Spectra::export(
  object = feat_export,
  BPPARAM = BiocParallel::SerialParam(),
  backend = MsBackendMgf::MsBackendMgf(),
  mapping = map,
  file = file.path(mgf_export_path, "ms1_features.mgf")
)