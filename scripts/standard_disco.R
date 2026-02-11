source("scripts/functions.R")
source("scripts/standard_disco_args.R")
suppressWarnings(
  suppressPackageStartupMessages({
    library(dplyr)
    # library(MSnbase)
    library(xcms)
    library(MsExperiment)
    library(optparse)
    library(writexl)
    library(BiocParallel)
  })
)

dir.create(file.path(stds_output_path, "peaks"), FALSE, TRUE)

stds <- import_mzml(data_path = data_path, meta_file = meta_file)

if (file.exists(file.path(stds_output_path, "stds_exp.rds"))) {
  stds_exp <- readRDS(file.path(stds_output_path, "stds_exp.rds"))
} else {
  stds_exp <- MsExperiment::readMsExperiment(
    spectraFiles = stds$path,
    sampleData = stds
  )
  saveRDS(stds_exp, file.path(stds_output_path, "stds_exp.rds"))
}

std_mz_theory <- get_theory_mz(chem_form = glycoside_form, adduct = std_adduct)
std_mz_range <- get_short_mz_range(std_mz_theory, mz_window = 0.02)

if (file.exists(file.path(stds_output_path, "std_chr_wide.rds"))) {
  std_chr <- readRDS(file.path(stds_output_path, "std_chr_wide.rds"))
} else {
  std_chr <- xcms::chromatogram(
    object = stds_exp,
    mz = std_mz_range,
    aggregationFun = "max"
  )
  saveRDS(std_chr, file.path(stds_output_path, "std_chr_wide.rds"))
}

# open a plotting window when running via Rscript/Rterm
if (!interactive()) {
  sys <- Sys.info()[["sysname"]]
  if (sys == "Windows") {
    grDevices::windows()
  } else if (sys == "Darwin") {
    grDevices::quartz()
  } else {
    grDevices::x11()   # requires X11; on SSH you need X-forwarding
  }
}

# Individual IS XICs
full_rm_samp <- c()
for (i in seq_along(std_chr)) {
  samp_idx <- which(colnames(std_chr)[i] == rownames(stds))
  plot(
    x = std_chr[, i],
    lwd = 3,
    main = paste0(
      stds$group[samp_idx], "\n",
      colnames(std_chr)[i]
    )
  )
  if (interactive()) {
    rm_samp <- readline("Selection: Good [1] Bad [2] ")
  } else if (!interactive()) {
    cat("Selection: Good [1] Bad [2] ")
    rm_samp <- readLines("stdin", n = 1)
  }

  rm_samp2 <- setNames(colnames(std_chr)[i], rm_samp)

  full_rm_samp <- c(full_rm_samp, rm_samp2)
}
invisible(dev.off())

good_samps <- full_rm_samp[which(names(full_rm_samp) == "1")]
good_chr <- std_chr[, which(colnames(std_chr) %in% good_samps)]

std_ranges <- get_rt_mz_range(chromatogram = good_chr, rt_window = 0.02)

# Get the IS XIC
if (file.exists(file.path(stds_output_path, "std_eic_wide.rds"))) {
  std_eic_wide <- readRDS(file.path(stds_output_path, "std_eic_wide.rds"))
} else {
  std_eic_wide <- xcms::chromatogram(
    stds_exp,
    mz = std_ranges$mz_range + c(-expand_mz, expand_mz),
    rt = std_ranges$rt_range + c(-expand_rt, expand_rt),
    aggregationFun = "sum"
  )
  saveRDS(std_eic_wide, file.path(stds_output_path, "std_eic_wide.rds"))
}

if (file.exists(file.path(stds_output_path, "std_chroms.rds"))) {
  std_chr2 <- readRDS(file.path(stds_output_path, "std_chroms.rds"))
} else {
  # Run peak detection on the EIC
  std_chr2 <- findChromPeaks(
    object = std_eic_wide,
    param = CentWaveParam(
      ppm = std_ppm,
      peakwidth = c(1, 50),
      prefilter = c(1, 1),
      snthresh = 1,
      mzCenterFun = "wMean",
      mzdiff = 0.001, # -0.001,
      integrate = 2,
      noise = 0,
      verboseBetaColumns = TRUE
    ),
    msLevel = 1
  )
  saveRDS(std_chr2, file.path(stds_output_path, "std_chroms.rds"))
}

for (i in seq_along(std_chr2)) {
  file_nm <- colnames(std_chr2)[i]
  samp_idx <- which(file_nm == rownames(stds))
  group_nm <- stds$group[samp_idx]
  clean_file_nm <- gsub(".mzML", "", file_nm)

  save_file <- file.path(
    stds_output_path, "peaks", paste0(group_nm, "_", clean_file_nm, ".pdf")
  )
  pdf(save_file, width = 4, height = 4)
  plot(
    x = std_chr2[, i],
    lwd = 3,
    main = paste0(
      group_nm, "\n",
      file_nm
    )
  )
  invisible(dev.off())
}

# TODO Double-check so this is correct later
std_chr_peaks <- tibble::as_tibble(
  chromPeaks(std_chr2),
  rownames = "peak_id"
) %>%
  dplyr::mutate(file = Biobase::sampleNames(std_chr2)[column]) %>%
  dplyr::left_join(
    x = .,
    y = as_tibble(stds, rownames = "file"),
    by = "file"
  )

write_xlsx(
  x = std_chr_peaks,
  path = file.path(stds_output_path, "std_chr_peaks.xlsx"),
  col_names = TRUE,
  format_headers = TRUE,
  use_zip64 = FALSE
)
message("Done.")