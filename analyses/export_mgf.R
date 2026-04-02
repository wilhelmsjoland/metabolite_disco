library(Spectra)
library(MsBackendMgf)
# library(BiocParallel)
# register(SerialParam())

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
  BPPARAM = SerialParam(),
  backend = MsBackendMgf(),
  mapping = map,
  file = file.path("output", "features.mgf")
)
