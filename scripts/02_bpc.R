# ==============================================================================
# Source functions and minimal startup parameters ------------------------------
# ==============================================================================
source("scripts/functions.R")
start_log(snakemake@params$output)
script_header()

set.seed(snakemake@params$seed)
register_parallel(snakemake@params$cores)
suppressWarnings(
  suppressPackageStartupMessages(
    {
      library(cli)
      library(BiocParallel)
      library(xcms)
      library(MsExperiment)
      library(Biobase)
      library(pheatmap)
      library(ggplot2)
      library(ProtGenerics)
    }
  )
)

# Load previous step
setup <- readRDS(snakemake@input[[1]])
ms_exp <- setup$ms_exp
group_colors <- setup$group_colors

# ==============================================================================
# Create and plot base peak chromatograms --------------------------------------
# ==============================================================================
bpc_path <- file.path(snakemake@params$output, "objects", "bpcs.rds")
if (interactive() && file.exists(bpc_path)) {
  bpcs <- readRDS(file = bpc_path)
  cli::cli_alert_success(
    paste0(
      "Imported saved base peak chromatograms from",
      " {.path {bpc_path}}"
    )
  )
} else {
  cli::cli_alert_info("Creating base peak chromatograms")
  bpcs <- xcms::chromatogram(
    BPPARAM = BiocParallel::SerialParam(),
    chunkSize = 1L,
    object = ms_exp,
    aggregationFun = "max"
  )
  saveRDS(object = bpcs, file = bpc_path)
  cli::cli_alert_success(
    paste0(
      "Saved base peak chromatograms to ",
      "{.path {bpc_path}}"
    )
  )
}

raw_bpc_p_path <- file.path(
  snakemake@params$output,
  "graphs",
  "bpc",
  "raw_bpc.pdf"
)
pdf(raw_bpc_p_path)
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
    "{.path {raw_bpc_p_path}}"
  )
)

# ==============================================================================
# Create a heatmap of base peak intensities ------------------------------------
# ==============================================================================
# Calculate correlation on the log2 transformed base peak intensities
bpcs_bin <- ProtGenerics::bin(bpcs, binSize = 1)
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
raw_bpc_hmp <- pheatmap::pheatmap(
  cormat,
  annotation = ann,
  annotation_color = list(group = group_colors),
  silent = TRUE
)

raw_bpc_hmp_path <- file.path(
  snakemake@params$output,
  "graphs",
  "bpc",
  "raw_bpc_hmp.pdf"
)
ggplot2::ggsave(
  filename = raw_bpc_hmp_path,
  plot = raw_bpc_hmp,
  device = "pdf",
  height = 10,
  width = 10,
  units = "in"
)
cli::cli_alert_success(
  paste0(
    "Saved base peak intensities heatmap to ",
    "{.path {raw_bpc_hmp_path}}"
  )
)

saveRDS(
  object = list(bpcs = bpcs),
  file = snakemake@output[[1]]
)
end_log()