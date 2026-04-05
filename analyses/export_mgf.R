library(Spectra)
library(MsBackendMgf)
library(BiocParallel)

# import from export_mgf.R
feature_levels

# # Extract MS1 spectra at feature apexes
feat_spectra <- xcms::featureSpectra(
  xchr9,
  msLevel = 1L,
  method = "closest_rt", # Makes sense for MS1 but not for MS2
  BPPARAM = SerialParam(),
  features = feature_levels
)

feat_spectra %>% spectraVariables()

feat_spectra$precursorCharge <- -1L

feat_export <- selectSpectraVariables(
  feat_spectra,
  c(
    "scanIndex",
    "mz",
    "intensity",
    "rtime",
    "precursorMz",
    "precursorCharge",
    "msLevel",
  )
)

spectraVariableMapping(MsBackendMgf())

map <- c(feature_id = "TITLE", spectraVariableMapping(MsBackendMgf()))

# Export
Spectra::export(
  object = feat_export,
  BPPARAM = BiocParallel::MulticoreParam(),
  backend = MsBackendMgf(),
  mapping = map,
  file = file.path("output", "features.mgf")
)

file.path(
  "/Volumes/bluecub/aglycone_release_100um_24h/output/experiment_mgfs",
  "features.mgf"
)
