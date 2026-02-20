cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Create and plot base peak chromatograms --------------------------------------
# ==============================================================================
if (check_saved("bpcs.rds")) {
  bpcs <- readRDS(file = file.path(opt$output, "objects/bpcs.rds"))
  cli::cli_alert_success(
    paste0(
      "Imported saved base peak chromatograms from",
      " {.path {file.path(opt$output, 'objects/bpcs.rds')}}"
    )
  )
} else {
  cli::cli_alert_info("Creating base peak chromatograms")
  bpcs <- xcms::chromatogram(ms_exp, aggregationFun = "max")
  saveRDS(object = bpcs, file = file.path(opt$output, "objects/bpcs.rds"))
  cli::cli_alert_success(
    paste0(
      "Saved base peak chromatograms to ",
      "{.path {file.path(opt$output, 'objects/bpcs.rds')}}"
    )
  )
}

pdf(file.path(opt$output, "graphs/bpc/raw_bpc.pdf"))
par(mar = c(4, 4, 3, 2))
plot(
  x = bpcs,
  col = group_colors[MsExperiment::sampleData(ms_exp)$group],
  main = "Base peak chromatogram"
)
invisible(dev.off())
cli::cli_alert_success(
  paste0(
    "Saved base peak chromatogram plot to ",
    "{.path {file.path(opt$output, 'graphs/bpc/raw_bpc.pdf')}}"
  )
)

# ==============================================================================
# Create a heatmap of base peak intensities ------------------------------------
# ==============================================================================
# Calculate correlation on the log2 transformed base peak intensities
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

cormat_rownames <- basename(Biobase::pData(xcms::phenoData(bpcs))$path)
colnames(cormat) <- rownames(cormat) <- cormat_rownames

# Define which phenodata columns should be highlighted in the plot
ann <- data.frame(group = bpcs_bin$group)
rownames(ann) <- cormat_rownames

# Perform the cluster analysis
cormat_p <- pheatmap::pheatmap(
  cormat,
  annotation = ann,
  annotation_color = list(group = group_colors),
  silent = TRUE
)

ggplot2::ggsave(
  filename = paste0(opt$output, "/graphs/bpc/raw_bpc_hmp.pdf"),
  plot = cormat_p,
  device = "pdf",
  height = 10,
  width = 10,
  units = "in"
)
cli::cli_alert_success(
  paste0(
    "Saved base peak intensities heatmap to ",
    "{.path {file.path(opt$output, '/graphs/bpc/raw_bpc_hmp.pdf')}}"
  )
)