library(Spectra)
library(MsBackendMgf)
# library(BiocParallel)
# register(SerialParam())

pot_glycosides

# # Extract MS1 spectra at feature apexes
feat_spectra <- xcms::featureSpectra(
  xchr9,
  msLevel = 1L,
  BPPARAM = SerialParam(),
  features = pot_glycosides
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
    "feature_id"
  )
)

spectraVariableMapping(MsBackendMgf())

map <- c(feature_id = "TITLE", spectraVariableMapping(MsBackendMgf()))
# map

# Export
Spectra::export(
  feat_export,
  backend = MsBackendMgf(),
  mapping = "map",
  file = file.path("output", "features.mgf")
)
