source("R/functions.R")
source("R/standard_disco_args.R")
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

dir.create(file.path(stds.output.path, "peaks"), FALSE, TRUE)

stds <- importFiles(data.path = data.path, meta.file = meta.file)

if (file.exists(file.path(stds.output.path, "stds_exp.rds"))) {
    stds.exp <- readRDS(file.path(stds.output.path, "stds_exp.rds"))
} else {
    stds.exp <- MsExperiment::readMsExperiment(spectraFiles = stds$path, sampleData = stds)
    saveRDS(stds.exp, file.path(stds.output.path, "stds_exp.rds"))
}

std.mz.theory <- getTheoryMz(chem_form = glycoside.form, adduct = std.adduct)
std.mz.range <- getShortMzRange(std.mz.theory, mz.window = 0.02)

if (file.exists(file.path(stds.output.path, "std_chr_wide.rds"))) {
    std.chr <- readRDS(file.path(stds.output.path, "std_chr_wide.rds"))
} else {
    std.chr <- xcms::chromatogram(object = stds.exp, mz = std.mz.range, aggregationFun = "max")
    saveRDS(std.chr, file.path(stds.output.path, "std_chr_wide.rds"))
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
full.rm.samp <- c()
for (i in 1:ncol(std.chr)) {
    samp.idx <- which(colnames(std.chr)[i] == rownames(stds))
    plot(
        x = std.chr[, i],
        lwd = 3,
        main = paste0(
            stds$group[samp.idx], "\n",
            colnames(std.chr)[i]
        )
    )
    if (interactive()) {
        rm.samp <- readline("Selection: Good [1] Bad [2] ")
    } else if (!interactive()) {
        cat("Selection: Good [1] Bad [2] ")
        rm.samp <- readLines("stdin", n = 1)
    }

    rm.samp2 <- setNames(colnames(std.chr)[i], rm.samp)

    full.rm.samp <- c(full.rm.samp, rm.samp2)
}
invisible(dev.off())

good.samps <- full.rm.samp[which(names(full.rm.samp) == "1")]
good.chr <- std.chr[,which(colnames(std.chr) %in% good.samps)]

std.ranges <- getRtMzRange(chromatogram = good.chr, rt_window = 0.02)

# Get the IS XIC
if (file.exists(file.path(stds.output.path, "std_eic_wide.rds"))) {
    std.eic.wide <- readRDS(file.path(stds.output.path, "std_eic_wide.rds"))
} else {
    std.eic.wide <- xcms::chromatogram(
        stds.exp,
        mz = std.ranges$mz.range + c(-expand.mz, expand.mz),
        rt = std.ranges$rt.range + c(-expand.rt, expand.rt),
        aggregationFun = "sum"
    )
    saveRDS(std.eic.wide, file.path(stds.output.path, "std_eic_wide.rds"))
}

if (file.exists(file.path(stds.output.path, "std_chroms.rds"))) {
    std.chr2 <- readRDS(file.path(stds.output.path, "std_chroms.rds"))
} else {
    # Run peak detection on the EIC
    std.chr2 <- findChromPeaks(
        object = std.eic.wide,
        param = CentWaveParam(
            ppm = std.ppm,
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
    saveRDS(std.chr2, file.path(stds.output.path, "std_chroms.rds"))
}

for (i in 1:ncol(std.chr2)) {
    file.nm <- colnames(std.chr2)[i]
    samp.idx <- which(file.nm == rownames(stds))
    group.nm <- stds$group[samp.idx]
    clean.file.nm <- gsub(".mzML", "", file.nm)

    save.file <- file.path(stds.output.path, "peaks", paste0(group.nm, "_", clean.file.nm, ".pdf"))
    pdf(save.file, width = 4, height = 4)
    plot(
        x = std.chr2[, i],
        lwd = 3,
        main = paste0(
            group.nm, "\n",
            file.nm
        )
    )
    invisible(dev.off())
}

# TODO Double-check so this is correct later
std.chr.peaks <- tibble::as_tibble(chromPeaks(std.chr2), rownames = "peak_id") %>%
    dplyr::mutate(file = Biobase::sampleNames(std.chr2)[column]) %>%
    dplyr::left_join(
        x = .,
        y = as_tibble(stds, rownames = "file"),
        by = "file"
    )

write_xlsx(
    x = std.chr.peaks,
    path = file.path(stds.output.path, "std_chr_peaks.xlsx"),
    col_names = TRUE,
    format_headers = TRUE,
    use_zip64 = FALSE
)