source("R/functions.R")
source("R/met_disco_args.R")
suppressWarnings(
    suppressPackageStartupMessages({
        library(tidyverse)
        library(MSnbase)
        library(xcms)
        library(MsExperiment)
        library(RforMassSpectrometry)
        library(Spectra)
        library(Chromatograms)
        library(RColorBrewer)
        library(pheatmap)
        library(QFeatures)
        library(Rdisop)
        library(limma)
        library(BiocParallel)
        library(ggrepel)
        library(ggfortify)
        library(ComplexUpset)
        library(MetaboAnnotation)
        library(CompoundDb)
        library(MetaboCoreUtils)
        library(curl)
        library(gt)
        library(patchwork)
        library(AnnotationHub)
        library(optparse)
        })
    )

# BiocParallel::register(SerialParam())
# bpstop()

################################################################################
# Create output folders
################################################################################
dir.create(res.folder, FALSE, TRUE)
dir.create(file.path(res.folder, "objects"), FALSE, TRUE)
dir.create(file.path(res.folder, "tables"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "bpc"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "internal_standard"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "volcano"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "upset"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "feature_boxplot"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "rtime"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "filled_peaks"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "pca"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "feature_chromatogram_intensity"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "per_sample_peaks"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "features"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "feature_pairs"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "glycoside"), FALSE, TRUE)
dir.create(file.path(res.folder, "graphs", "glycoside_feature_pairs"), FALSE, TRUE)
dir.create(file.path("annotation_databases"), FALSE, TRUE)

################################################################################
# Import metadata
################################################################################
message("Importing metadata...")
meta <- importFiles(data.path, meta.file)
if (check_saved("ms_exp.rds")) {
    ms.exp <- readRDS(file = file.path(res.folder, "objects/ms_exp.rds"))
} else {
    ms.exp <- MsExperiment::readMsExperiment(
        spectraFiles = meta$path,
        sampleData = meta
    )
    saveRDS(object = ms.exp, file = file.path(res.folder, "objects/ms_exp.rds"))
}

################################################################################
# Set colors for groups
################################################################################
message("Setting colors for groups...")
groups.to.use <- unique(MsExperiment::sampleData(ms.exp)$group)
group.colors <- paste0(brewer.pal(n = length(groups.to.use), "Set1")[1:length(groups.to.use)])
group.colors <- setNames(group.colors, groups.to.use)
 
################################################################################
# Create and plot base peak chromatograms
################################################################################
message("Creating base peak chromatograms...")
if (check_saved("bpcs.rds")) {
    bpcs <- readRDS(file = file.path(res.folder, "objects/bpcs.rds"))
} else {
    bpcs <- xcms::chromatogram(ms.exp, aggregationFun = "max")
    saveRDS(object = bpcs, file = file.path(res.folder, "objects/bpcs.rds"))
}

message("Plotting base peak chromatograms...")
pdf(file.path(res.folder, "graphs/bpc/raw_bpc.pdf"))
par(mar = c(4, 4, 3, 2))
plot(
    x = bpcs,
    col = group.colors[MsExperiment::sampleData(ms.exp)$group],
    main = "Base peak chromatogram"
)
invisible(dev.off())

################################################################################
# Create a heatmap of base peak intensities
################################################################################
## Calculate correlation on the log2 transformed base peak intensities
bpcs_bin <- bin(bpcs, binSize = 1)
cormat <- cor(
    log2(
        do.call(
            cbind,
            lapply(bpcs_bin, intensity)
        )
    ),
    use = "complete.obs" # Because NAs
)

cormat.rownames <- basename(Biobase::pData(xcms::phenoData(bpcs))$path)
colnames(cormat) <- rownames(cormat) <- cormat.rownames

## Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = bpcs_bin$group)
rownames(ann) <- cormat.rownames

## Perform the cluster analysis

cormat.p <- pheatmap::pheatmap(
    cormat,
    annotation = ann,
    annotation_color = list(group = group.colors)
)

ggplot2::ggsave(
    filename = paste0(res.folder, "/graphs/bpc/raw_bpc_hmp.pdf"),
    plot = cormat.p,
    device = "pdf",
    height = 10,
    width = 10,
    units = "in"
)

################################################################################
# Part 1 - inspect IS prior to peak-calling
# Inspect the internal standard
# Define the rt and m/z range of the peak area
################################################################################

message("Inspecting internal standard peaks prior to peak-calling...")
mz.theory <- getTheoryMz(chem_form = internal_standard, adduct = adduct)
mz.range <- getShortMzRange(mz.theory, mz.window = 0.02)
if (check_saved("is_chr.rds")) {
    is.chr <- readRDS(file = paste0(res.folder, "/objects/is_chr.rds"))
} else {
    is.chr <- xcms::chromatogram(object = ms.exp, mz = mz.range, aggregationFun = "sum")
    saveRDS(object = is.chr, file = paste0(res.folder, "/objects/is_chr.rds"))
}
ranges <- getRtMzRange(chromatogram = is.chr, rt_window = 0.02)

# Wide IS chromatogram
pdf(paste0(res.folder, "/graphs/internal_standard/all_is_wide.pdf"))
plot(x = is.chr, col = group.colors[is.chr$group], lwd = 3)
legend("topright", legend = names(group.colors), col = group.colors, pch = 16)
invisible(dev.off())

# Get the IS XIC
if (check_saved("is_eic.rds")) {
    is.eic <- readRDS(file = paste0(res.folder, "/objects/is_eic.rds"))
} else {
    is.eic <- xcms::chromatogram(
        object = ms.exp,
        mz = ranges$mz.range,
        rt = ranges$rt.range,
        aggregationFun = "sum"
    )
    saveRDS(object = is.eic, file = paste0(res.folder, "/objects/is_eic.rds"))
}

# All IS EICs together
pdf(paste0(res.folder, "/graphs/internal_standard/all_is.pdf"))
plot(x = is.eic, col = group.colors[is.eic$group],lwd = 3)
legend("topleft", legend = names(group.colors), col = group.colors, pch = 16)
invisible(dev.off())

# Individual IS XICs
# for (i in 1:ncol(is.eic)) {
#     plot(
#         x = is.eic[, i],
#         col = group.colors[names(group.colors) %in% dplyr::filter(meta, rownames(meta) == colnames(bpcs)[i])$group],
#         lwd = 3,
#         main = colnames(bpcs)[i]
#     )
#     readline("Enter for next: ")
# }
# invisible(dev.off())


# TODO FIX
# The estimation of the ppm should be fixed
# try4 <- findMzError(object = ms.exp, ref_peaks_obj = ranges$data)
# max(full.tib$max.mz.error)

################################################################################
# 2. Check if data is centroided
################################################################################

# This filters the amount of speectra not the mz?
# sps <- spectra(ms.exp) |>
#     filterRt(getRtMzRange(chromatogram = is.chr, rt_window = 0.02)$rt.range) |> 
#     filterMzRange(mz = getRtMzRange(chromatogram = is.chr, rt_window = 0.02)$mz.range)
# 
# Spectra::plotSpectraOverlay(x = sps, lwd = 2, xlim = mz.range)
# 
# ms.exp |>
#     filterRt(getRtMzRange(chromatogram = is.chr, rt_window = 0.02)$rt.range) |> 
#     filterMzRange(mz = getRtMzRange(chromatogram = is.chr, rt_window = 0.02)$mz.range) |>
#     plot()
# 
# # Looks pretty good but is the data smoothed - it looks that way?
# # looks like centroid mode
# # TODO Add the ranges used in this example under
# ms.exp |>
#     Spectra::filterRt(rt.range.long) |>
#     Spectra::filterMzRange(c(is.theory.mz - 0.025, is.theory.mz + 0.025)) |>
#     plot()
# invisible(dev.off())
# 
# srn <- ms.exp[1] |>
#     filterRt(rt = rt.range.long) |>
#     filterMzRange(mz = mz.range)
# plot(srn)
# invisible(dev.off())
# 
# srn_1 <- srn[1] |>
#     filterMzRange(c(123.052, 123.040)) |>
#     spectra()
# 
# # Calculate the difference in m/z values between scans
# mz_diff <- srn_1 |>
#     mz() |>
#     unlist() |>
#     diff() |>
#     abs()
# mz_diff
# 
# # Differences in m/z values expressed as ppm
# mz_error <- mz_diff * 1e6 / mean(unlist(mz(srn_1)))
# 
# max_error <- max(mz_error, na.rm = TRUE)
# round_max_error <- round(max_error, -1)
# round_max_error <- 60 # just use this for now

# message(
#     "Calculated max error: ", max_error, " ppm", 
#     "\nUsing rounded max error: ", round_max_error, " ppm",
#     sep = ""
#     )

# shorter range
# ms.exp |>
#     filterRt(rt.range.s) |>
#     filterMzRange(mz.range) |>
#     plot()
# invisible(dev.off())

################################################################################
# Part 3 
# - Ensure IS peak is chosen
# - Determine max and min peak width for CentWaveParam() from the IS.
#  TODO - do this for several peaks and not only the IS
################################################################################

message("Determining minimal and maximal peakwidth based on the internal standard...")

if (check_saved("is_eic_wide.rds")) {
    is.eic.wide <- readRDS(file = paste0(res.folder, "/objects/is_eic_wide.rds"))
} else {
    is.eic.wide <- xcms::chromatogram(
        ms.exp,
        mz = mz.range + c(-0.05, 0.05),
        rt = ranges$rt.range + c(-16, 16),
        aggregationFun = "sum"
    )
    saveRDS(object = is.eic.wide, file = paste0(res.folder, "/objects/is_eic_wide.rds"))
}


# Run peak detection on the EIC
if (check_saved("is_chr2.rds")) {
    is.chr2 <- readRDS(file = paste0(res.folder, "/objects/is_chr2.rds"))
} else {
    is.chr2 <- xcms::findChromPeaks(
        object = is.eic.wide,
        param = xcms::CentWaveParam(
            ppm = ppm.global,
            peakwidth = c(2, 20),
            prefilter = c(1, 1),
            snthresh = sn_threshold, # 10
            mzCenterFun = "wMean",
            mzdiff = 0.001,
            integrate = 2,
            noise = 1000,
            verboseBetaColumns = TRUE
        ),
        msLevel = 1
        )
    saveRDS(object = is.chr2, file = paste0(res.folder, "/objects/is_chr2.rds"))
}

# Calculate peakwidth
is.peaks <- tibble::as_tibble(xcms::chromPeaks(is.chr2), rownames = "rownames") %>%
    dplyr::rowwise() %>%
    dplyr::mutate(delta_rt = rtmax - rtmin) %>%
    dplyr::ungroup()


# Min: to half of some peaks in the datasets
# Max: 2-4x times the average size
is.min.peak.width <- min(is.peaks$delta_rt, na.rm = TRUE)
is.max.peak.width <- max(is.peaks$delta_rt, na.rm = TRUE)

min.peak.width <- quantile(is.peaks$delta_rt, 0.05, na.rm = TRUE) * 0.3 # 0.3
max.peak.width <- quantile(is.peaks$delta_rt, 0.95, na.rm = TRUE) * 4 # 4

message(
    "Internal standard: ", "C7H8O2", ", Theoretical m/z: ", mz.theory,
    "\n\tMinimal peak width of IS: ", is.min.peak.width,
    "\n\tMaximal peak width of IS: ", is.max.peak.width,
    "\nPeak widths used for peak picking: ",
    "\n\t Minimal width: ", min.peak.width,
    "\n\t Maximal width: ", max.peak.width,
    sep = ""
    )

################################################################################
# 4. Call peaks on whole dataset with params
################################################################################

message("Calling peaks...")

if (check_saved("xchr.rds")) {
    xchr <- readRDS(file = paste0(res.folder, "/objects/xchr.rds"))
} else {
    xchr <- xcms::findChromPeaks(
        object = ms.exp,
        BPPARAM = BiocParallel::bpparam(),
        return.type = "XCMSnExp",
        msLevel = 1L,
        param = xcms::CentWaveParam(
            ppm = ppm.global,
            peakwidth = c(min.peak.width, max.peak.width),
            snthresh = sn_threshold, # 10
            prefilter = c(4, 1000), # k peaks (left) over intensity (right) # c(3, 100)
            mzCenterFun = "wMean",
            integrate = 2,
            mzdiff = 0.001,
            fitgauss = FALSE,
            noise = 1000,
            verboseColumns = TRUE,
            roiList = list(),
            firstBaselineCheck = TRUE,
            roiScales = numeric(),
            extendLengthMSW = FALSE,
            verboseBetaColumns = TRUE
            )
        )
    saveRDS(object = xchr, file = paste0(res.folder, "/objects/xchr.rds"))
}

plot(
    bpcs, 
    col = group.colors[MsExperiment::sampleData(xchr)$group],
    main = "Base peak chromatogram after peak picking"
    )
legend(
    "topright",
    col = unique(group.colors[MsExperiment::sampleData(xchr)$group]), 
    legend = unique(names(group.colors[MsExperiment::sampleData(ms.exp)$group])), 
    pch = 16)
invisible(dev.off())


# Use this for the beta-distribution parameters
xchr_data <- tibble::as_tibble(xcms::chromPeaks(xchr), rownames = "peak")

################################################################################
# 5. Inspect peaks
# Extract some peaks here and check quality of peak picking
################################################################################
message("Inspecting called peaks...")

# peaks_to_inspect <- as_tibble(chromPeaks(xchr[1:2]), rownames = "peak") %>%
#     tidyr::drop_na(beta_snr) %>%
#     dplyr::arrange(desc(into)) %>%
#     dplyr::slice(1:10) %>%
#     dplyr::pull(peak) %>%
#     stringr::str_remove_all(., "[A-Za-z]") %>%
#     as.numeric(.)
# 
# for (peak in peaks_to_inspect) {
#     inspect_peaks <- inspectPeak(
#         chromatogram = xchr[1:2],
#         peak_data = as_tibble(chromPeaks(xchr[1:2]), rownames = "peak"),
#         peak.idx = peak,
#         sample.no = FALSE,
#         save.graph = TRUE
#     )
# }


# Poor peaks: beta_cor < 0.5 (or even < 0.2
# Good peaks: beta_snr > 7
# Keep signal to noise at 10, and filter by beta_cor < 0.5, and beta_snr > 7?

# These are the per-sample peak counts
inspectPeakInt(chr_data = xchr, value = into, save.graph = TRUE)

################################################################################
# 5. Refine peaks
# Let these be part of the pipeline as well and as a choice to do or not
################################################################################

if (check_saved("xchr2.rds")) {
    xchr2 <- readRDS(paste0(res.folder, "/objects/xchr2.rds"))
} else {
    xchr2 <- xcms::refineChromPeaks(
        object = xchr,
        param = xcms::CleanPeaksParam(maxPeakwidth = max.peak.width)
    )
    saveRDS(object = xchr2, file = paste0(res.folder, "/objects/xchr2.rds"))
}

if (check_saved("xchr3.rds")) {
    xchr3 <- readRDS(paste0(res.folder, "/objects/xchr3.rds"))
} else {
    xchr3 <- xcms::refineChromPeaks(
        object = xchr2,
        param = xcms::FilterIntensityParam(
            threshold = 0,
            nValues = 1L,
            value = "maxo"
        )
    )
    saveRDS(object = xchr3, file = paste0(res.folder, "/objects/xchr3.rds"))
}

if (check_saved("xchr4.rds")) {
    xchr4 <- readRDS(paste0(res.folder, "/objects/xchr4.rds"))
} else {
    xchr4 <- xcms::refineChromPeaks(
        object = xchr, # xchr3
        param = xcms::MergeNeighboringPeaksParam(
            expandRt = 0.25,
            expandMz = 0,
            ppm =  ppm.global,
            minProp = 0.95 # between 0 & 1
        )
    )
    saveRDS(object = xchr4, file = paste0(res.folder, "/objects/xchr4.rds"))
}

# TODO 
# Make a check for if these exist
clean.removed.peaks <- dplyr::anti_join(
    x = tibble::as_tibble(xcms::chromPeaks(xchr), rownames = "peaks"),
    y = tibble::as_tibble(xcms::chromPeaks(xchr2), rownames = "peaks"),
    by = "peaks"
)
clean.removed.peaks %>% nrow

intensity.removed.peaks <- dplyr::anti_join(
    x = tibble::as_tibble(xcms::chromPeaks(xchr2), rownames = "peaks"),
    y = tibble::as_tibble(xcms::chromPeaks(xchr3), rownames = "peaks"),
    by = "peaks"
)
intensity.removed.peaks %>% nrow

merged.peaks <- dplyr::anti_join(
    x = tibble::as_tibble(xcms::chromPeaks(xchr3), rownames = "peaks"),
    y = tibble::as_tibble(xcms::chromPeaks(xchr4), rownames = "peaks"),
    by = "peaks"
)

################################################################################
# 6. Align retention times
# Correspondence, adjust and correspondence again 
################################################################################
message("Aligning retention times across samples...")

# TODO Change these to more broad so I don't get a million anchor peaks?
if (check_saved("xchr5.rds")) {
    xchr5 <- readRDS(paste0(res.folder, "/objects/xchr5.rds"))
} else {
    xchr5 <- xcms::groupChromPeaks(
        object = xchr, # xchr4
        param = xcms::PeakDensityParam(
            sampleGroups = MsExperiment::sampleData(xchr)$group,
            bw = bw_first_grouping,
            minFraction = 0.5, # 0.7
            binSize = 0.01,
            maxFeatures = 1000, # 200
            ppm = ppm.global,
            minSamples = 2 # 1
            )
        )

    saveRDS(object = xchr5, file = paste0(res.folder, "/objects/xchr5.rds"))
}

# TODO Use this part to check if anchor peaks cover the whole RT
# Get the anchor peaks that would be selected
pgm <- xcms::adjustRtimePeakGroups(
    xchr5, 
    xcms::PeakGroupsParam(minFraction = 0.9)
    )

# Evaluate distribution of anchor peaks' rt in the first sample
# TODO compare this to the full range of retention times from the runs
chrom_peaks5 <- xcms::chromPeaks(xchr5)
max.rt.time <- max(chrom_peaks5[, "rtmax"], na.rm = TRUE)
min.rt.time <- min(chrom_peaks5[, "rtmin"], na.rm = TRUE)

# TODO Add a check here as well!
# quantile(pgm[, 1], na.rm = TRUE) # Remove na.rm = TRUE?
# End of checking for anchor peak rt distribution

# Alignment
if (check_saved("xchr6.rds")) {
    xchr6 <- readRDS(paste0(res.folder, "/objects/xchr6.rds"))
} else {
    xchr6 <- xcms::adjustRtime(
        object = xchr5,
        param = xcms::PeakGroupsParam(
            minFraction = 0.8, # 0.9
            extraPeaks = 0, # 1
            smooth = "loess",
            span = 0.6, # 0.6
            family = "gaussian",
            peakGroupsMatrix = matrix(nrow = 0, ncol = 0),
            subset = integer(),
            subsetAdjust = c("average", "previous")
        )
    )
    saveRDS(object = xchr6, file = paste0(res.folder, "/objects/xchr6.rds"))
}

# Checking for retention drift
# Extract base peak chromatograms
if (check_saved("bpc_after.rds")) {
    bpc_after <- readRDS(file = paste0(res.folder, "/objects/bpc_after.rds"))
} else {
    bpc_after <- xcms::chromatogram(xchr6, aggregationFun = "max", chromPeaks = "none")
    saveRDS(object = bpc_after, file = paste0(res.folder, "/objects/bpc_after.rds"))
}

pdf(paste0(res.folder, "/graphs/rtime/before_after_alignment.pdf"))
par(mfrow = c(2, 1))
# Before retention time alignment
plot(
    bpcs, 
    col = group.colors[MsExperiment::sampleData(xchr6)$group],
    main = "Before retention time alignment"
)

# After retention time alignment
plot(
    bpc_after, 
    col = group.colors[MsExperiment::sampleData(xchr6)$group],
    main = "After retention time alignment"
)
invisible(dev.off())

# Checking for retention drift in IS
if (check_saved("is_drift_check_before.rds")) {
    is.drift.check.before <- readRDS(file = paste0(res.folder, "/objects/is_drift_check_before.rds"))
} else {
    is.drift.check.before <- xchr %>%
        Spectra::filterRt(ranges$rt.range) |>
        Spectra::filterMzRange(ranges$mz.range) |>
        xcms::chromatogram(
            aggregationFun = "max",
            chromPeaks = "none"
        )
    saveRDS(object = is.drift.check.before, file = paste0(res.folder, "/objects/is_drift_check_before.rds"))
}

if (check_saved("is_drift_check_after.rds")) {
    is.drift.check.after <- readRDS(file = paste0(res.folder, "/objects/is_drift_check_after.rds"))
} else {
    is.drift.check.after <- xchr6 %>%
        Spectra::filterRt(ranges$rt.range) |>
        Spectra::filterMzRange(ranges$mz.range) |>
        xcms::chromatogram(
            aggregationFun = "max",
            chromPeaks = "none"
            )
    saveRDS(object = is.drift.check.after, file = paste0(res.folder, "/objects/is_drift_check_after.rds"))
}

# Checking the adjustment in the IS peak
pdf(paste0(res.folder, "/graphs/rtime/is_before_after_alignment.pdf"))
par(mfrow = c(1, 2))
plot(
    is.drift.check.before,
    col = group.colors[MsExperiment::sampleData(xchr6)$group],
    main = "Before:\nRT: 130 - 175 (s)\nM/z range: 226.9 - 228",
    lwd = 3
)

plot(
    is.drift.check.after,
    col = group.colors[MsExperiment::sampleData(xchr6)$group],
    main = "After:\nRT: 130 - 175 (s)\nM/z range: 226.9 - 228",
    lwd = 3
)
invisible(dev.off())

# From Sattely paper
# Retention time correction was performed using the obiwarp method, with a 
# step size of m/z 0.5. Peak alignment was performed with bandwidth 
# of 3 seconds and minimum fraction (minfrac) of samples 
# necessary for a valid group of 0.5.

################################################################################
# 7. Correspondence
# Extract a chromatogram for a m/z range containing internal standard
# Test these settings on the extracted slice
# Do this for several
################################################################################

message("Producing simulated bandwidth plots...")
# Check bandwidth
if (check_saved("chr_1.rds")) {
    chr_1 <- readRDS(file = paste0(res.folder, "/objects/chr_1.rds"))
} else {
    chr_1 <- xcms::chromatogram(
        object = xchr6,
        mz = ranges$mz.range,
        rt = ranges$rt.range + c(-16, 16)
    )
    saveRDS(object = chr_1, file = paste0(res.folder, "/objects/chr_1.rds"))
}

# Test these settings on the extracted slice
pdf(paste0(res.folder, "/graphs/internal_standard/is_simul_first_grouping.pdf"))
density.simul.p <- xcms::plotChromPeakDensity(
    object = chr_1, 
    param = xcms::PeakDensityParam(
        sampleGroups = MsExperiment::sampleData(xchr6)$group, 
        bw = bw_first_grouping
        ))
invisible(dev.off())

pdf(paste0(res.folder, "/graphs/internal_standard/is_simul_second_grouping.pdf"))
density.simul.p <- xcms::plotChromPeakDensity(
    object = chr_1, 
    param = xcms::PeakDensityParam(
        sampleGroups = MsExperiment::sampleData(xchr6)$group, 
        bw = bw_second_grouping
    ))
invisible(dev.off())
# End of checking for bandwidth to use

message("Performing correspondence...")
# Correspondence
if (check_saved("xchr7.rds")) {
    xchr7 <- readRDS(file = paste0(res.folder, "/objects/xchr7.rds"))
} else {
    xchr7 <- xcms::groupChromPeaks(
        object = xchr6,
        param = xcms::PeakDensityParam(
            sampleGroups = MsExperiment::sampleData(xchr6)$group,
            bw = bw_second_grouping, # 0.5
            minFraction = 0.5, # T0.7
            binSize = 0.01,
            maxFeatures = 1000, # 200
            ppm = ppm.global,
            minSamples = 2 # 1
        )
    )
    saveRDS(object = xchr7, file = paste0(res.folder, "/objects/xchr7.rds"))
}

# TODO Check this for several ions as well
# Extract chromatogram including signal for is
if (check_saved("chr_2.rds")) {
    chr_2 <- readRDS(file = paste0(res.folder, "/objects/chr_2.rds"))
} else {
    chr_2 <- xcms::chromatogram(
        object = xchr7,
        mz = ranges$mz.range,
        rt = ranges$rt.range + c(-16, 16),
        aggregationFun = "max"
        )
    saveRDS(object = chr_2, file = paste0(res.folder, "/objects/chr_2.rds"))
}

message("Producing second grouping bandwidth plots...")
# Setting simulate = FALSE to show the actual correspondence results
pdf(paste0(res.folder, "/graphs/internal_standard/is_non_simul_second_grouping.pdf"))
density.non.simul.p <- xcms::plotChromPeakDensity(
    object = chr_2, 
    simulate = FALSE
    )
invisible(dev.off())

################################################################################
# 8. Gap filling
################################################################################
# Checking features
feat.def <- as_tibble(xcms::featureDefinitions(xchr7), rownames = "feature")
feat.val <- as_tibble(xcms::featureValues(xchr7, method = "sum"), rownames = "feature")

# Extract features with nas for peak filling
feat.with.na <- feat.val %>%
    tidyr::pivot_longer(cols = 2:ncol(.)) %>%
    dplyr::filter(is.na(value)) %>%
    dplyr::pull(feature) %>%
    unique(.)

# Filter feature defintions to features with nas
feat.def.nas <- feat.def %>%
    dplyr::filter(feature %in% feat.with.na)

# Create a list for checking the features chromatograms
feat.def.nas.vals <- feat.def.nas %>%
    dplyr::rowwise() %>%
    dplyr::mutate(feat_extract = rbind(
        c(
            mzmed - 0.0015, 
            mzmed + 0.0015,
            rtmin - 2,
            rtmax + 2
            )
        )) %>%
    dplyr::pull(feat_extract)

message("Amount of features with NAs prior to gap filling: ", sum(is.na(featureValues(xchr7))))

# Perform gap filling
if (check_saved("xchr8.rds")) {
    xchr8 <- readRDS(file = paste0(res.folder, "/objects/xchr8.rds"))
} else {
    xchr8 <- xcms::fillChromPeaks(
        object = xchr7,
        param = xcms::ChromPeakAreaParam()
        )
    saveRDS(object = xchr8, file = paste0(res.folder, "/objects/xchr8.rds"))
}

# Number of missing values after gap filling
message("Amount of features with NAs after gap filling: ", sum(is.na(xcms::featureValues(xchr8))))

# Number of filled peaks
message("Number of filled peaks: ", sum(is.na(featureValues(xchr7))) - sum(is.na(xcms::featureValues(xchr8))))

# Plot all non-filled peaks
# Extract the m/z - rt regions for these features
# Extract features with nas for peak filling
feat.with.na.after <- tibble::as_tibble(
    xcms::featureValues(xchr8, method = "sum"), 
    rownames = "feature"
    ) %>%
    tidyr::pivot_longer(cols = 2:ncol(.)) %>%
    dplyr::filter(is.na(value)) %>%
    dplyr::pull(feature) %>%
    unique(.)

chrs.na.feat <- xcms::featureArea(xchr8, features = feat.with.na.after)

# Expand the retention time by 1 second on both sides
chrs.na.feat[, "rtmin"] <- chrs.na.feat[, "rtmin"] - 1
chrs.na.feat[, "rtmax"] <- chrs.na.feat[, "rtmax"] + 1

if (check_saved("chrs_na.rds")) {
    chrs_na <- readRDS(file = paste0(res.folder, "/objects/chrs_na.rds"))
} else {
    chrs_na <- xcms::chromatogram(
        xchr8,
        mz = chrs.na.feat[, c("mzmin", "mzmax")],
        rt = chrs.na.feat[, c("rtmin", "rtmax")] # probably increase this a little
        )
    saveRDS(object = chrs_na, file = paste0(res.folder, "/objects/chrs_na.rds"))
}

################################################################################
# 8. Filtering features and input to SummarizedExperiment
################################################################################

# Filter features that are uninteresting
# TODO https://rformassspectrometry.github.io/Metabonaut/articles/a-end-to-end-untargeted-metabolomics.html#filtering-features-missing-values

message("Filtering features based on missingness...")

group.factor <- MsExperiment::sampleData(xchr8)$group
group.factor <- as.factor(group.factor)

if (check_saved("xchr9.rds")) {
    xchr9 <- readRDS(file = paste0(res.folder, "/objects/xchr9.rds"))
} else {
    xchr9 <- QFeatures::filterFeatures(
        xchr8, 
        xcms::PercentMissingFilter(
            threshold = missing_threshold,
            f = group.factor
        ))
    saveRDS(object = xchr9, file = paste0(res.folder, "/objects/xchr9.rds"))
}

################################################################################
# 9. Median scaling & PCA
################################################################################

message("Producing PCA plots...")

res <- xcms::quantify(
    xchr9,
    method = "sum",
    value = "into",
    filled = FALSE,
    missing = "rowmin_half" # 0 ?
)

SummarizedExperiment::assays(res)$raw_filled <- xcms::featureValues(
    xchr9, 
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = "rowmin_half" # 0 ?
)

# Compute median and generate normalization factor
mdns <- apply(
    SummarizedExperiment::assay(
        res, 
        "raw"
    ), 
    MARGIN = 2,
    median, 
    na.rm = TRUE 
)
nf_mdn <- mdns / median(mdns)

# Dividing dataset by median of median and creating a new assay
SummarizedExperiment::assays(res)$norm <- sweep(
    SummarizedExperiment::assay(res, "raw"), 
    MARGIN = 2, 
    nf_mdn, 
    '/'
)

# Compute median and generate normalization factor
mdns <- apply(
    SummarizedExperiment::assay(
        res, 
        "raw_filled"
    ), 
    MARGIN = 2,
    median, 
    na.rm = TRUE 
)
nf_mdn <- mdns / median(mdns)

# Dividing dataset by median of median and creating a new assay
SummarizedExperiment::assays(res)$norm_filled <- sweep(
    SummarizedExperiment::assay(res, "raw_filled"), 
    MARGIN = 2, 
    nf_mdn, 
    '/'
)

# Log2 transform and scale data
vals <- SummarizedExperiment::assay(res, "raw_filled") |>
    log2() |>
    t() |>
    scale(center = TRUE, scale = TRUE) |>
    as.matrix(.)

pca_res <- prcomp(vals, scale = FALSE, center = FALSE)

# Data before normalization
vals_st <- cbind(vals, phenotype = res$group)
pca_raw <- autoplot(
    pca_res, 
    data = vals_st,
    colour = 'phenotype', 
    scale = 0,
    size = 3
    ) +
    ggplot2::scale_color_manual(values = group.colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Before normalization")

# Data after normalization
vals_norm <- SummarizedExperiment::assay(res, "norm_filled") |>
    log2() |>
    t() |>
    scale(center = TRUE, scale = TRUE) |>
    as.matrix(.)

pca_res_norm <- prcomp(vals_norm, scale = FALSE, center = FALSE)
vals_st_norm <- cbind(vals_norm, phenotype = res$group)
pca_adj <- autoplot(
    pca_res_norm, 
    data = vals_st_norm,
    colour = 'phenotype', 
    scale = 0,
    size = 3
    ) +
    ggplot2::scale_color_manual(values = group.colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "After normalization")

norm.filled.12.pca.p <- pca_raw / pca_adj +
    patchwork::plot_layout(guides = "collect")

ggplot2::ggsave(
    filename = paste0(res.folder, "/graphs/pca/norm_filled_pca_1_2.pdf"),
    plot = norm.filled.12.pca.p,
    device = "pdf",
    height = 6,
    width = 6,
    units = "in"
)

# PC 3-4 before & after normalization
pca_raw <- autoplot(
    pca_res, 
    data = vals_st,
    colour = 'phenotype', 
    x = 3, 
    y = 4, 
    scale = 0,
    size = 3
    ) +
    ggplot2::scale_color_manual(values = group.colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Before normalization")

pca_adj <- autoplot(
    pca_res_norm, 
    data = vals_st_norm,
    colour = 'phenotype', 
    x = 3, 
    y = 4, 
    scale = 0,
    size = 3
    ) +
    ggplot2::scale_color_manual(values = group.colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "After normalization")

norm.filled.34.pca.p  <- pca_raw / pca_adj +
    patchwork::plot_layout(guides = "collect")

ggplot2::ggsave(
    filename = paste0(res.folder, "/graphs/pca/norm_filled_pca_3_4.pdf"),
    plot = norm.filled.34.pca.p,
    device = "pdf",
    height = 6,
    width = 6,
    units = "in"
)

################################################################################
# 9. Limma
################################################################################

message("Running linear models...")

group.used <- factor(meta$group)
design <- model.matrix(~ 0 + group.used)
colnames(design) <- levels(group.used)

comparisons <- combn(
    x = levels(group.used),
    m = 2,
    simplify = TRUE
    ) %>%
    t(.) %>%
    tibble::as_tibble(.) %>%
    dplyr::mutate(comp = paste0(V1, "-", V2)) %>%
    dplyr::pull(comp)

contrasts.mat <- limma::makeContrasts(
    contrasts = comparisons, 
    levels = design
    )

fit <- limma::lmFit(
    log2(SummarizedExperiment::assay(res, "norm")), # don't use imputed
    design = design
    )
fit <- limma::contrasts.fit(fit, contrasts.mat)
fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)

limma_res <- list()
for (i in comparisons) {
    tmp <- limma::topTable(
        fit, 
        coef = i,
        number = Inf, 
        adjust.method = "BH",
        sort.by = "none"
        ) %>%
        tibble::as_tibble(., rownames = "feature") %>%
        dplyr::mutate(contrast = i) %>%
        dplyr::left_join(
            x = .,
            y = SummarizedExperiment::rowData(res) %>%
                tibble::as_tibble(., rownames = "feature"),
            by = "feature"
        )
    
    limma_res[[i]] <- tmp
}

full.limma <- tibble::tibble()
for (i in names(limma_res)) {
    tmp.tib <- limma_res[[i]]
    full.limma <- dplyr::bind_rows(full.limma, tmp.tib)
}

message("Producing volcano plots...")
volc_plot_list <- list()
for (i in names(limma_res)) {
    
    tmp <- limma_res[[i]]
    
    tmp.tib <- tmp %>%
        dplyr::mutate(
            label.p = dplyr::if_else(
                adj.P.Val < p.value.global & abs(logFC) > quantile(abs(logFC), 0.99), # -3 
                mzmed,
                NA
            ),
            direction = dplyr::case_when(
                logFC >= 0 & adj.P.Val < p.value.global ~ "Up",
                logFC < 0 & adj.P.Val < p.value.global ~ "Down",
                adj.P.Val >= p.value.global ~ "ns",
                TRUE ~ as.character("check")
            )
        ) %>%
        dplyr::mutate(direction = forcats::fct_relevel(
            direction, 
            c(
                "Up",
                "ns",
                "Down"
            )
        )) %>%
        tidyr::drop_na(logFC)
    
    tmp.p <- tmp.tib %>%
        ggplot2::ggplot(.,
               ggplot2::aes(
                   x = logFC,
                   y = -log10(adj.P.Val),
                   color = direction
               )) +
        ggplot2::geom_point() +
        ggplot2::scale_color_manual(values = c(
            "Up" = "firebrick",
            "Down" = "cornflowerblue",
            "ns" = "grey"
        )) +
        ggrepel::geom_label_repel(
            data = tidyr::drop_na(tmp.tib, label.p) %>%
                dplyr::arrange(dplyr::desc(abs(logFC))) %>%
                dplyr::slice(1:50),
            ggplot2::aes(label = round(label.p, 2)),
            size = 3,
            max.overlaps = 100,
            box.padding = 0.5,
            color = "black"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.title = ggplot2::element_blank(),
            plot.subtitle = ggplot2::element_text(hjust = 0.5)
        ) +
        ggplot2::labs(
            title = i,
            subtitle = paste0(
                "Rounded mass-to-charge ratios with adj.p < ",
                p.value.global,
                " are labelled"
                ),
            x = "Log2 fold change",
            y = "-Log10 adjusted p-value"
        )
    
    volc_plot_list[[i]] <- tmp.p
    
    ggplot2::ggsave(
        filename = paste0(res.folder, "/graphs/volcano/", i, ".pdf"),
        plot = tmp.p,
        device = "pdf",
        height = 10,
        width = 10,
        units = "in"
    )
}

message("Saving intensity information to tables...")
# Creating tables of all output data and saving to tables
# TODO Fix this dumb logic here or use as witch statement?
assay.names <- names(SummarizedExperiment::assays(res))
for (i in assay.names) {
    if (i %in% c("norm", "norm_filled")) {
        full_data <- dplyr::left_join(
            x = SummarizedExperiment::rowData(res) %>%
                tibble::as_tibble(., rownames = "feature"),
            y = SummarizedExperiment::assay(res, i) %>%
                log2() %>%
                t() %>%
                scale(., center = TRUE, scale = TRUE) %>%
                t() %>%
                tibble::as_tibble(., rownames = "feature"),
            by = "feature"
        )
    } else if (i %in% c("raw", "raw_filled")) {
        full_data <- dplyr::left_join(
            x = SummarizedExperiment::rowData(res) %>%
                tibble::as_tibble(., rownames = "feature"),
            y = SummarizedExperiment::assay(res, i) %>%
                tibble::as_tibble(., rownames = "feature"),
            by = "feature"
        )
    }

    assign(
        x = paste0("full_", i),
        value = full_data,
        envir = .GlobalEnv
    )
    
    readr::write_csv(
        x = full_data,
        file = paste0(res.folder, "/tables/full_", i, ".csv"),
        na = "NA",
        col_names = TRUE,
        append = FALSE
    )
}

################################################################################
# 10. Check intersections
################################################################################

message("Producing upset plots...")

upset.tib <- full.limma %>%
    dplyr::select(feature, contrast, adj.P.Val) %>%
    tidyr::pivot_wider(
        names_from = "contrast",
        values_from = "adj.P.Val"
    ) %>%
    dplyr::mutate(
        dplyr::across(
            .cols = 2:ncol(.),
            .fns = ~ dplyr::if_else(
                . < p.value.global,
                TRUE,
                FALSE
                )
        )
    )

upset.p <- upset.tib %>%
    ComplexUpset::upset(
        intersect = comparisons,
        name = paste0("Features with p adjusted < ", p.value.global),
        width_ratio = 0.15,
        base_annotations = list(
            "Intersecting features" = ComplexUpset::intersection_size()
        )
    )

# upset.tib %>% 
#     dplyr::select(c("feature", all_of(upset.comps))) %>%
#     ComplexUpset::upset(
#         intersect = upset.comps,
#         name = paste0("Features with p adjusted < ", p.value.global),
#         width_ratio = 0.15,
#         base_annotations = list(
#             "Intersecting features" = ComplexUpset::intersection_size()
#         ),
#         min_degree = 3
#     )

ggplot2::ggsave(
    filename = paste0(res.folder, "/graphs/upset/upset.pdf"),
    plot = upset.p,
    device = "pdf",
    height = 7,
    width = 22,
    units = "in"
)

################################################################################
# 11. Find intersecting features
################################################################################

message("Finding intersecting features...")

# TODO Create all interesting comparisons
upset.comp <- c(
    "bu_wt_apiin-bu_wt_control", # 1
    "bu_mutant_control-bu_wt_apiin", # 1
    "bu_mutant_apiin-bu_wt_apiin" # , # 1
    # "bu_mutant_apiin-bu_wt_control" # 2
)

# Add to list
upset.comps <- list()
intersecting.feats <- findIntersectFeat(
    data = upset.tib,
    set = upset.comp,
    full_set = comparisons
)
upset.comps[[stringr::str_flatten(upset.comp, collapse = "*")]] <- intersecting.feats$feature
# TODO END of part that needs fixing

all.int.comps <- upset.tib %>%
    tidyr::pivot_longer(cols = 2:ncol(.)) %>%
    dplyr::filter(value == TRUE) %>%
    dplyr::group_by(feature) %>%
    dplyr::filter(dplyr::n() >= 2) %>%
    dplyr::pull(feature) %>%
    unique(.)

intersect.data <- full_norm_filled %>%
    tidyr::pivot_longer(cols = tidyselect::all_of(rownames(meta))) %>%
    dplyr::filter(feature %in% intersecting.feats$feature) %>%
    dplyr::left_join(
        x = .,
        y = meta %>%
            tibble::as_tibble(., rownames = "sample"),
        by = c("name" = "sample")
    ) %>%
    dplyr::mutate(group = forcats::fct_relevel(
        group,
        c(
            "bu_wt_control",
            "bu_wt_apiin",
            "bu_mutant_control",
            "bu_mutant_apiin"
        )
    ))

limma.p.res <- full.limma %>%
    dplyr::select(feature, adj.P.Val, contrast) %>%
    dplyr::mutate(
        group1 = stringr::str_split_i(contrast, "-", 1),
        group2 = stringr::str_split_i(contrast, "-", 2),
    ) %>%
    dplyr::select(-contrast) %>%
    dplyr::filter(feature %in% intersecting.feats$feature) %>%
    rstatix::add_significance(p.col = "adj.P.Val") %>%
    find_y_position(
        test.df = .,
        df = intersect.data,
        formula = "value ~ group",
        fun.data = "max",
        grouping = "feature"
    )

message("Producing significant intersecting feature boxplots...")
# Intersecting feats individually
for (i in intersecting.feats$feature) {
    tmp.inter.data <- full_norm_filled %>%
        tidyr::pivot_longer(cols = tidyselect::all_of(rownames(meta))) %>%
        dplyr::filter(feature %in% i) %>%
        dplyr::left_join(
            x = .,
            y = meta %>%
                tibble::as_tibble(., rownames = "sample"),
            by = c("name" = "sample")
        ) %>%
        dplyr::mutate(group = forcats::fct_relevel(
            group,
            c(
                "bu_wt_control",
                "bu_wt_apiin",
                "bu_mutant_control",
                "bu_mutant_apiin"
            )
        ))
    
    tmp.inter.limma <- full.limma %>%
        dplyr::select(feature, adj.P.Val, contrast) %>%
        dplyr::mutate(
            group1 = stringr::str_split_i(contrast, "-", 1),
            group2 = stringr::str_split_i(contrast, "-", 2),
        ) %>%
        dplyr::select(-contrast) %>%
        dplyr::filter(feature %in% i) %>%
        rstatix::add_significance(p.col = "adj.P.Val") %>%
        find_y_position(
            test.df = .,
            df = tmp.inter.data,
            formula = "value ~ group",
            fun.data = "max"
        )
    
    tmp.inter.p <- tmp.inter.data %>%
        ggplot2::ggplot(.,
               ggplot2::aes(
                   x = group,
                   y = value
               )) +
        ggplot2::geom_boxplot(
            ggplot2::aes(fill = group),
            outliers = FALSE
        ) +
        ggplot2::geom_point() +
        ggplot2::facet_wrap(
            facets = ~ feature, 
            scales = "free_y"
        ) +
        ggplot2::scale_y_continuous(expand = ggplot2::expansion(c(0.1, 0.15))) +
        ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            strip.background = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            axis.title.x = ggplot2::element_blank(),
            legend.title = ggplot2::element_blank()
        ) +
        ggplot2::labs(
            y = "Log2 median-scaled intensity"
        ) +
        ggpubr::geom_bracket(
            data = tmp.inter.limma %>%
                dplyr::filter(!adj.P.Val.signif %in% c("ns")),
            ggplot2::aes(
                xmin = group1,
                xmax = group2,
                label = adj.P.Val.signif,
                y.position = y.pos * 1.005
            ),
            step.increase = 0.1
        )
    
    ggplot2::ggsave(
        filename = paste0(res.folder, "/graphs/feature_boxplot/", i, ".pdf"),
        plot = tmp.inter.p,
        device = "pdf",
        height = 5,
        width = 5,
        units = "in"
    )
}

################################################################################
# 12. Predictions of m/z
# 1. First find interesting ones to look at (p.val & logFC)
# 2. Match for the ones with a hit in the prediction of biotransformations
# 3. Plot them from the actual data
################################################################################

# TODO Set an optparse argument for:
# polarity
# or by inputting specific adducts
message("Expanding possible adducts...")
xchr9.defs <- xcms::featureDefinitions(xchr9) %>%
    tibble::as_tibble(., rownames = "feature")

xchr9.mzs <- xchr9.defs$mzmed
names(xchr9.mzs) <- xchr9.defs$feature

# TODO Check so correct
possible.adducts <- MetaboCoreUtils::mz2mass(
    xchr9.mzs, # peak.mz
    adduct = adducts(polarity = "negative")
    ) %>%
    tibble::as_tibble(., rownames = "feature") %>%
    tidyr::pivot_longer(
        cols = 2:ncol(.),
        names_to = "adduct",
        values_to = "mass"
    ) %>%
    dplyr::left_join(
        x = .,
        y = xchr9.defs,
        by = "feature"
    )

if (nrow(xchr9.defs) * 17 == nrow(possible.adducts)) {
    message("Adducts expanded")
} else {
    warning("Adducts did not correctly expanded")
}

message("Importing biotransformation file...")
bio.transf <- importBiotransfMeta(file = paste0(data.path, "/", biotransf.file))

message(
    "Predicting potential biotransformations based on:\n\t", 
    "Biotransformation database: ", biotransf.file,
    "\n\tppm: ", ppm.global,
    sep = ""
    )

# TODO
# CHECK IF THIS MAKES SENSE NOW!!!!!
# CHECK THIS FOR SURE
# TODO
# Fix so the observed ppm is added
matched.diffs <- predictBiotransfAdducts(
    data = possible.adducts,
    biotransf.data = bio.transf,
    tolerance_ppm = 1 # ppm.global
)

# TODO This is slow for now
message("Writing predictions to table...")
writexl::write_xlsx(
    x = matched.diffs,
    path = paste0(res.folder, "/tables/matched_diffs.xlsx"),
    col_names = TRUE,
    format_headers = TRUE,
    use_zip64 = FALSE
)

message("
 based on significant features in contrasts:\n\t",
    paste0(upset.comps, "\n\t"),
    sep = ""
    )

filt.match.diffs <- matched.diffs %>%
    dplyr::filter(
        dplyr::if_any(
            dplyr::all_of(c("feat1", "feat2")),
            # ~ .x %in% intersecting.feats$feature
            ~ .x %in% all.int.comps
        )
    ) %>%
    dplyr:::mutate(
        pair = purrr::map2(feat1, feat2, ~ c(.x, .y)),
        mz1_forms = purrr::map(mz1, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 0))), # or ppm global
        mz2_forms = purrr::map(mz2, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 0))) # or ppm global
    ) %>%
    dplyr::filter(grepl("1 x", name))

message("Writing filtered predictions to table...")
writexl::write_xlsx(
    x = filt.match.diffs %>%
        dplyr::select(-c(
            "mz1_forms",
            "mz2_forms",
            "pair"
        )),
    path = paste0(res.folder, "/tables/filt_matched_diffs.xlsx"),
    col_names = TRUE,
    format_headers = TRUE,
    use_zip64 = FALSE
)

################################################################################
# 13. Predictions of m/z vs a database
################################################################################

int.mets <- c(filt.match.diffs$feat1, filt.match.diffs$feat2)

peaks_used <- full_norm_filled %>%
    dplyr::select(feature, mzmed, rtmed) %>%
    dplyr::rename("mz" = "mzmed", "rtime" = "rtmed") %>%
    dplyr::filter(feature %in% int.mets) %>%
    tibble::column_to_rownames(var = "feature")

peaks_used$peak_id <- rownames(peaks_used) # keep XCMS peak IDs

# Load a compound database via CompDb
# The chebi import doesn't work

annotation_hub <- AnnotationHub()
# query(annotation_hub, "CompDb")
cdb <- annotation_hub[["AH119519"]]

# dbname  <- "CompDb.Hsapiens.HMDB.5.0.sqlite"
# db_file <- file.path("annotation_databases", dbname) # temp.dir()
# if (!file.exists(db_file)) {
#     curl::curl_download(
#         url = paste0(
#             "https://github.com/jorainer/MetaboAnnotationTutorials/",
#             "releases/download/2021-11-02/", dbname
#         ),
#         destfile = db_file
#     )
# }
# cdb <- CompoundDb::CompDb(db_file)

target_df <- ProtGenerics::compounds(
    cdb,
    columns = c(
        "compound_id",
        "name",
        "formula",
        "exactmass",
        "smiles",
        "inchi",
        "inchikey",
        "cas",
        "pubchem"
        )
)

# Check potential adducts
# MetaboCoreUtils::adducts(polarity = "negative")
# MetaboCoreUtils::adductNames(polarity = "negative")

# parameters to match by
mz_match_param <- MetaboAnnotation::Mass2MzParam(
    adducts = c(MetaboCoreUtils::adductNames(polarity = "negative")),
    ppm = ppm.global # 10
)

matches <- MetaboAnnotation::matchValues(
    query = peaks_used,
    target = target_df,
    param = mz_match_param
)

anno <- MetaboAnnotation::matchedData(matches) %>%
    tibble::as_tibble(., rownames = "feature") %>%
    dplyr::mutate(abs_score = abs(score)) %>%
    dplyr::arrange(abs_score)

################################################################################
# 14. Filtering chromatogram object
################################################################################

# TODO
# Filter this weaker
# Don't filter as hard, and keep all the significant ones only?
# For the delta m/z comparisons -> use filt.features + all.int.comps
# For the plotting of individuals features -> use only the all.int.comps

message(sprintf(
    "Filtering features with sn: %s, beta_cor: %s, beta_snr: %s", 
    sn_threshold, 
    beta_cor_threshold, 
    beta_snr_threshold
    ))
if (check_saved("xchr9_filt.rds")) {
    xchr9.filt <- readRDS(file = paste0(res.folder, "/objects/xchr9_filt.rds"))
} else {
    xchr9.filt <- filtFeatures(
        object = xchr9,
        sn_threshold = sn_threshold,
        beta_cor_threshold = beta_cor_threshold, 
        beta_snr_threshold = beta_snr_threshold,
        filt_vector = all.int.comps
    )
    saveRDS(object = xchr9.filt, file = paste0(res.folder, "/objects/xchr9_filt.rds"))
}

# Checking specifically for the glycoside anad aglycone m/zs
glycoside <- MetaboCoreUtils::mass2mz(
    Rdisop::getMass(Rdisop::getMolecule(glycoside_form)),
    adduct = MetaboCoreUtils::adducts(polarity = "negative")
    ) %>%
    t() %>%
    tibble::as_tibble(., rownames = "adduct") %>%
    dplyr::rename("glycoside" = V1)

aglycone <- MetaboCoreUtils::mass2mz(
    Rdisop::getMass(Rdisop::getMolecule(aglycone_form)),
    adduct = MetaboCoreUtils::adducts(polarity = "negative")
    ) %>%
    t() %>%
    tibble::as_tibble(., rownames = "adduct") %>%
    dplyr::rename("aglycone" = V1)

range.tol <- ppm_to_num(glycoside_ppm)

gly.agly.adducts <- glycoside %>%
    dplyr::left_join(
        x = .,
        y = aglycone,
        by = "adduct"
    ) %>%
    dplyr::mutate(
        glycoside.min = glycoside - range.tol,
        glycoside.max = glycoside + range.tol,
        aglycone.min = aglycone - range.tol,
        aglycone.max = aglycone + range.tol
    )

gly.agly <- tibble::tibble()
for (i in 1:nrow(gly.agly.adducts)) {
    tmp <- full_raw_filled %>%
        dplyr::filter(
            dplyr::between(mzmed, gly.agly.adducts[i,]$glycoside.min, gly.agly.adducts[i,]$glycoside.max) |
            dplyr::between(mzmed, gly.agly.adducts[i,]$aglycone.min, gly.agly.adducts[i,]$aglycone.max)
        ) %>%
        dplyr::mutate(adduct = gly.agly.adducts[i,]$adduct) %>%
        dplyr::relocate(adduct, .after = "feature")
    
    gly.agly <- bind_rows(gly.agly, tmp)
}

pot.glycosides <- unique(gly.agly$feature)

# TODO Change this to match the Adduct one
matched.diffs2 <- predictBiotransfAdductsSubset(
    data = possible.adducts,
    biotransf.data = bio.transf,
    tolerance_ppm = 15, # glycoside_ppm
    feat_filt = pot.glycosides
)

filt.match.diffs2 <- matched.diffs2 %>%
    dplyr::rowwise() %>%
    dplyr::mutate(pair = list(c(feat1, feat2))) %>%
    dplyr::ungroup()

glycoside.pairs <- unique(c(filt.match.diffs2$feat1, filt.match.diffs2$feat2))

xchr9.filt$final.plotting.features <- unique(c(
    glycoside.pairs, 
    xchr9.filt$final.plotting.features
    ))

################################################################################
# 14. Plotting interesting features
################################################################################

message("Producing feature chromatograms...")
if (check_saved("feature_chrs.rds")) {
    feature.chrs <- readRDS(file = paste0(res.folder, "/objects/feature_chrs.rds"))
} else {
    feature.chrs <- xcms::featureChromatograms(
        object = xchr9,
        expandRt = 0,
        expandMz = 0,
        aggregationFun = "sum",
        filled = TRUE,
        features = xchr9.filt$final.plotting.features,
        missing = 0,
        # features = rownames(xcms::featureDefinitions(xchr9)),
        return.type = "XChromatograms"
    )
    saveRDS(object = feature.chrs, file = paste0(res.folder, "/objects/feature_chrs.rds"))
}

# # Note
# # The EIC data of a feature is extracted from every sample using the 
# # same m/z - rt area. The EIC in a sample does thus not exactly represent the 
# # signal of the actually identified chromatographic peak in that sample. 
# # The chromPeakChromatograms() function would allow to extract the actual EIC 
# # of the chromatographic peak in a specific sample. See also examples below.
 
# message("Writing feature chromatograms to plots...")
# # Features
# plotChromatogramInts(
#     chromatogram = feature.chrs,
#     chrom_object = xchr9,
#     save_loc = "/graphs/features/",
#     # amount = 1,
#     peaks_or_feats = "features"
# )
# 
# message("Writing all gap filled peaks to plots...")
# # Gap filled peaks only
# plotChromatogramInts(
#     chromatogram = chrs_na,
#     chrom_object = xchr8,
#     save_loc = "/graphs/filled_peaks/",
#     # amount = 1,
#     peaks_or_feats = "peaks"
# )

# message("Writing feature chromatograms and intensity boxplots...")
# # feats.to.plot <- rownames(xcms::featureDefinitions(feature.chrs))
# feats.to.plot <- xchr9.filt$filt.sig.features
# for (i in feats.to.plot) {
#     ft.p <- plotFeatChrInt(
#         feature_chrom = feature.chrs,
#         feature = i,
#         method = "sum",
#         value = "into",
#         filled = TRUE,
#         missing = "rowmin_half",
#         msLevel = 1,
#         save_loc = "/graphs/feature_chromatogram_intensity/",
#         device = "pdf",
#         feat_pairs = FALSE
#     )
# }
# 
# message("Plotting feature pairs in filtered biotransformations...")
# for (i in seq_len(nrow(xchr9.filt$biot.filt.sig.features.tib))) {
#     ft.pair.p <- plotFeatPairs(
#         feature_chrom = feature.chrs,
#         filt.match.row = xchr9.filt$biot.filt.sig.features.tib[i,],
#         method = "sum",
#         value = "into",
#         filled = TRUE,
#         missing = 0,
#         msLevel = 1,
#         save_pairs_loc = "/graphs/feature_pairs/",
#         device = "pdf"
#     )
# }
# 
# message("Writing feature chromatograms and intensity boxplots for glycosides/aglycones...")
# for (i in pot.glycosides) {
#     ft.p <- plotFeatChrInt(
#         feature_chrom = feature.chrs,
#         feature = i,
#         method = "sum",
#         value = "into",
#         filled = TRUE,
#         missing = "rowmin_half",
#         msLevel = 1,
#         save_loc = "/graphs/glycoside/",
#         device = "pdf",
#         feat_pairs = FALSE
#     )
# }
# 
# message("Plotting glycoside/aglycone feature pairs biotransformations...")
# for (i in seq_len(nrow(filt.match.diffs2))) {
#     ft.pair.p <- plotFeatPairs(
#         feature_chrom = feature.chrs,
#         filt.match.row = filt.match.diffs2[i,],
#         method = "sum",
#         value = "into",
#         filled = TRUE,
#         missing = 0,
#         msLevel = 1,
#         save_pairs_loc = "/graphs/glycoside_feature_pairs/",
#         device = "pdf"
#     )
# }
