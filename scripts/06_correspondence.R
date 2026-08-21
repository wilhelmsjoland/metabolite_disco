# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
source("scripts/functions.R")
start_log(snakemake@params$output)
script_header()

set.seed(snakemake@params$seed)
register_parallel(snakemake@params$cores)
suppressWarnings(
  suppressPackageStartupMessages({
    library(cli)
    library(BiocParallel)
    library(xcms)
    library(MsExperiment)
    library(RSQLite) # for sql backend
    library(MsBackendSql) # for sql backend
  })
)

alignment <- readRDS(snakemake@input[["alignment"]])
is_data <- readRDS(snakemake@input[["internal_std"]])

xchr6 <- alignment$xchr6
ranges <- is_data$ranges

# ==============================================================================
# Generating chromatograms for simulated bandwiths -----------------------------
# before correspondence for internal standard ----------------------------------
# ==============================================================================
if (is.null(ranges)) {
  cli::cli_alert_info(
    "No internal standard - skipping simulated bandwidth plots"
  )
} else {
cli::cli_alert_info(
  "Producing simulated bandwidth plots for internal standard"
)
bw_chr_1_path <- file.path(
  snakemake@params$output,
  "objects",
  "bw_chr_1.rds"
)
if (interactive() && file.exists(bw_chr_1_path)) {
  bw_chr_1 <- readRDS(file = bw_chr_1_path)
  cli::cli_alert_success(
    paste0(
      "Imported internal standard chromatogram for bandwidth simulation ",
      "for internal standard from: ",
      "{.path {bw_chr_1_path}}"
    )
  )
} else {
  cli::cli_alert_info(
    "Generating chromatogram for internal standard for bandwidth simulation"
  )
  bw_chr_1 <- xcms::chromatogram(
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1,
    object = xchr6,
    mz = ranges$mz_range,
    rt = ranges$rt_range + c(-16, 16)
  )
  saveRDS(object = bw_chr_1, file = bw_chr_1_path)
  cli::cli_alert_success(
    paste0(
      "Saved internal standard chromatogram for bandwidth simulation ",
      "for internal standard to: ",
      "{.path {bw_chr_1_path}}"
    )
  )
}

# ==============================================================================
# Plotting simulated bandwiths before correspondence for internal standard -----
# ==============================================================================
is_simul_first_grouping_path <- file.path(
  snakemake@params$output,
  "graphs",
  "internal_standard",
  "is_simul_first_grouping.pdf"
)
pdf(is_simul_first_grouping_path)
xcms::plotChromPeakDensity(
  object = bw_chr_1,
  param = xcms::PeakDensityParam(
    sampleGroups = MsExperiment::sampleData(xchr6)$group,
    bw = snakemake@params$bw_first_grouping
  )
)
invisible(dev.off())
cli::cli_alert_success(
  paste0(
    "Saved simulated bandwidth plot for first grouping ",
    "(bw_first_grouping: {.val {snakemake@params$bw_first_grouping}}) to: ",
    "{.path {is_simul_first_grouping_path}}"
  )
)

is_simul_second_grouping_path <- file.path(
  snakemake@params$output,
  "graphs",
  "internal_standard",
  "is_simul_second_grouping.pdf"
)
pdf(is_simul_second_grouping_path)
xcms::plotChromPeakDensity(
  object = bw_chr_1,
  param = xcms::PeakDensityParam(
    sampleGroups = MsExperiment::sampleData(xchr6)$group,
    bw = snakemake@params$bw_second_grouping
  )
)
invisible(dev.off())
cli::cli_alert_success(
  paste0(
    "Saved simulated bandwidth plot for second grouping ",
    "(bw_second_grouping: {.val {snakemake@params$bw_second_grouping}}) to: ",
    "{.path {is_simul_second_grouping_path}}"
  )
)
}

# ==============================================================================
# Correspondence ---------------------------------------------------------------
# ==============================================================================
cli::cli_h3("Performing correspondence")

xchr7_path <- file.path(
  snakemake@params$output,
  "objects",
  "xchr7.rds"
)
if (interactive() && file.exists(xchr7_path)) {
  xchr7 <- readRDS(file = xchr7_path)
  cli::cli_alert_success(
    paste0(
      "Imported saved correspondence object from: ",
      "{.path {xchr7_path}}"
    )
  )
} else {
  cli::cli_alert_info("Performing correspondence")
  xchr7 <- xcms::groupChromPeaks(
    object = xchr6,
    param = xcms::PeakDensityParam(
      sampleGroups = MsExperiment::sampleData(xchr6)$group,
      bw = snakemake@params$bw_second_grouping, # 0.5
      minFraction = 0.5, # 0.7
      binSize = 0.01,
      maxFeatures = 1000, # 200
      ppm = snakemake@params$ppm_global,
      minSamples = 2 # 1
    )
  )
  saveRDS(object = xchr7, file = xchr7_path)
  cli::cli_alert_success(
    paste0(
      "Saved correspondence object to: ",
      "{.path {xchr7_path}}"
    )
  )
}

# ==============================================================================
# Generating internal standard chromatograms for -------------------------------
# bandwith checking after correspondence ---------------------------------------
# ==============================================================================
if (is.null(ranges)) {
  cli::cli_alert_info(
    "No internal standard - skipping post-correspondence bandwidth plots"
  )
} else {
cli::cli_alert_info(
  paste0(
    "Producing bandwidth plot for correspondence parameters",
    "for internal standard"
  )
)
bw_chr_2_path <- file.path(
  snakemake@params$output,
  "objects",
  "bw_chr_2.rds"
)
if (interactive() && file.exists(bw_chr_2_path)) {
  bw_chr_2 <- readRDS(file = bw_chr_2_path)
  cli::cli_alert_success(
    paste0(
      "Imported saved chromatogram for internal standards from: ",
      "{.path {bw_chr_2_path}}"
    )
  )
} else {
  cli::cli_alert_info(
    "Generating chromatogram of internal standards for bandwidth ploting"
  )
  bw_chr_2 <- xcms::chromatogram(
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1,
    object = xchr7,
    mz = ranges$mz_range,
    rt = ranges$rt_range + c(-16, 16),
    aggregationFun = "max"
  )
  saveRDS(object = bw_chr_2, file = bw_chr_2_path)
  cli::cli_alert_success(
    paste0(
      "Saved chromatogram for internal standards to: ",
      "{.path {bw_chr_2_path}}"
    )
  )
}

# ==============================================================================
# Bandwith plotting after correspondence ---------------------------------------
# ==============================================================================
is_non_simul_second_group_path <- file.path(
  snakemake@params$output,
  "graphs",
  "internal_standard",
  "is_non_simul_second_grouping.pdf"
)
pdf(is_non_simul_second_group_path)
xcms::plotChromPeakDensity(
  object = bw_chr_2,
  simulate = FALSE
)
invisible(dev.off())
cli::cli_alert_success(
  paste0(
    "Saved bandwidth plot for internals standards for ",
    "parameters used for correspondence to: ",
    "{.path {is_non_simul_second_group_path}}"
  )
)
}

# ==============================================================================
# Snakesave --------------------------------------------------------------------
# ==============================================================================

saveRDS(
  object = list(xchr7 = xchr7),
  file = snakemake@output[[1]]
)

script_footer()
end_log()