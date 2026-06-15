# TODO
# Check intersection between significantly different
# by mutant and annotated metabolites

source("scripts/functions.R")
library(rcdk)
library(rJava)
library(fingerprint)
# Visualizations
library(cowplot)
library(magick)


apig_mz <- featureDefinitions(xchr9) %>%
  tibble::as_tibble(., rownames = "feature") %>%
  dplyr::filter(grepl("FT02089", feature)) %>%
  pull(mzmed)

anno_peak_ids <- unique(anno$peak_id)

anno_chrs <- xcms::featureChromatograms(
  object = xchr9,
  expandRt = 0,
  expandMz = 0,
  aggregationFun = "sum",
  filled = TRUE,
  features = anno.peak.ids,
  missing = 0,
  # features = rownames(xcms::featureDefinitions(xchr9)),
  return.type = "XChromatograms"
)

# for (i in anno.peak.ids) {
#     tmp.anno.chr <- plotFeatChrInt(
#         feature_chrom = anno.chrs,
#         feature = i,
#         method = "sum",
#         value = "into",
#         filled = TRUE,
#         missing = "rowmin_half",
#         msLevel = 1,
#         save_loc = NULL,
#         device = NULL,
#         feat_pairs = FALSE
#     )
#     print(tmp.anno.chr$combined)
# 
#     readline("Enter for next: ")
# }

# TODO Split more
int_annos <- c(
  "FT02161",
  "FT00292",
  "FT03755",
  "FT00485",
  "FT01696",
  "FT06069",
  "FT01657",
  "FT02948",
  "FT03779",
  "FT02246",
  "FT03341",
  "FT02024",
  "FT09148",
  "FT04929",
  "FT03009",
  "FT01252",
  "FT00081",
  "FT05181",
  "FT04812",
  "FT04750", # int no prod in 2 groups
  "FT01901",
  "FT01345",
  "FT05216",
  "FT01428",
  "FT00099",
  "FT05532",
  "FT01840",
  "FT01655",
  "FT02539",
  "FT05753",
  "FT00762",
  "FT00049",
  "FT02718", # weird
  "FT05496",
  "FT05809",
  "FT02012",
  "FT00307",
  "FT07141",
  "FT03793",
  "FT01327",
  "FT01404",
  "FT01658",
  "FT00183",
  "FT03331", # int
  "FT00157",
  "FT03095",
  "FT00103",
  "FT05777",
  "FT08395",
  "FT01133",
  "FT02027",
  "FT03826",
  "FT06714",
  "FT00372",
  "FT08744",
  "FT06988", # only in controls
  "FT02994",
  "FT07149", # only in mutant
  "FT01879", # super low in mutant apiin
  "FT00475",
  "FT06113",
  "FT04384",
  "FT03101", # int
  "FT03999",
  "FT02109", # almost only in mutant apiin
  "FT06809",
  "FT07476",
  "FT02952",
  "FT04898", # none in mutant apiin
  "FT10427",
  "FT04154",
  "FT06532", # int, strong peak
  "FT02265",
  "FT03040",
  "FT01936",
  "FT02266",
  "FT06508",
  "FT01385",
  "FT00397",
  "FT01640",
  "FT04480",
  "FT06507", # int, big difference
  "FT03435",
  "FT06584", # int - almost only in controls
  "FT07425", # int
  "FT02898",
  "FT05339",
  "FT01151", # only in mutant apiin?
  "FT02134",
  "FT09648",
  "FT08147",
  "FT09207",
  "FT01121",
  "FT04548",
  "FT01150", # weird
  "FT05943",
  "FT01119",
  "FT04168",
  "FT07026", # int only in controls
  "FT05666", # int only in controls
  "FT05782", # int only in controls
  "FT07752",
  "FT07522",
  "FT06533",
  "FT06244", # int only in controls
  "FT05874", # much higher in mutants
  "FT10419", # much higher in only mutant apiin
  "FT09078",
  "FT11262", # int only in controls
  "FT04555",
  "FT06998", # int much higher in controls
  "FT10557", # int much higher in controls
  "FT06635", # int higher in controls
  "FT06935", # int higher in controls
  "FT07834",
  "FT09271", # int only in apiin
  "FT06876", # int only in apiin
  "FT10018", # int only in apiin
  "FT02104",
  "FT07269",
  "FT05914", # int, probably different peaks in mutants and wt
  "FT06738",
  "FT05991",
  "FT10696", # weird
  "FT06263",
  "FT08063",
  "FT09227",
  "FT07139",
  "FT06789",
  "FT03730",
  "FT10864",
  "FT07014",
  "FT07616",
  "FT08733", # int only in controls
  "FT09819", # int only in mutant control
  "FT10162",
  "FT06676",
  "FT10048",
  "FT10885",
  "FT09160",
  "FT09410", # int only in controls
  "FT11592",
  "FT08566",
  "FT09268", # int higher in controls
  "FT10851"
)

anno %>%
  dplyr::filter(peak_id %in% int_annos) %>%
  dplyr::arrange(abs_score)

anno_peak_chr <- plot_feat_chrom_int(
  feature_chrom = anno_chrs,
  feature = "FT09271",
  method = "sum",
  value = "into",
  filled = TRUE,
  missing = "rowmin_half",
  ms_level = 1,
  save_loc = NULL,
  device = NULL,
  feat_pairs = FALSE
)
print(anno_peak_chr$combined)


# library(ChemmineR)
# library(ChemmineOB)
# 
# smiles <- anno$target_smiles
# names(smiles) <- anno$target_name
# 
# # 2) SMILES -> SDFset
# sdfset <- smiles2sdf(smiles)  # uses ChemmineOB/OpenBabel if available :contentReference[oaicite:1]{index=1}
# 
# # (optional) drop invalid structures if any
# sdfset <- sdfset[validSDF(sdfset)]
# 
# # 3) Atom-pair descriptors -> fixed-length binary fingerprints
# apset <- sdf2ap(sdfset)  # :contentReference[oaicite:2]{index=2}
# fpset <- desc2fp(apset, descnames = 1024, type = "FPset")  # atom-pair FPset :contentReference[oaicite:3]{index=3}
# 
# # 4) Full NxN Tanimoto similarity matrix (0..1)
# simMA <- sapply(cid(fpset), function(id) fpSim(fpset[id], fpset, sorted = FALSE, method = "Tanimoto"))
# # fpSim returns similarity coefficients; method can be Tanimoto/Dice/Tversky/etc. :contentReference[oaicite:4]{index=4}
# 
# # 5) Dissimilarity / distance matrix (0..1)
# dissMA <- 1 - simMA
# diag(dissMA) <- 0

smiles <- anno$target_smiles
names(smiles) <- anno$feature

mols <- rcdk::parse.smiles(
  smiles = smiles,
  omit.nulls = TRUE
)

# fps <- lapply(mols, get.fingerprint, type='circular')
# 
# fp.sim <- fingerprint::fp.sim.matrix(fps, method='tanimoto')
# fp.dist <- 1 - fp.sim
# 
# cls <- hclust(as.dist(fp.dist))
# plot(cls, labels=FALSE)

query.mol <- rcdk::parse.smiles("C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O")[[1]]
target.mols <- mols
query.fp <- get.fingerprint(query.mol, type = 'pubchem') # circular
target.fps <- lapply(target.mols, get.fingerprint, type = 'pubchem') # circular
sims <- data.frame(
    sim = do.call(rbind, 
                  lapply(
                      target.fps,
                      fingerprint::distance,
                      fp2 = query.fp, 
                      method = 'tanimoto'
                      )
                  ))

tani <- tibble::as_tibble(sims, rownames = "smiles") %>%
    dplyr::arrange(desc(sim))

tani.comb <- tani %>%
    dplyr::left_join(
        x = .,
        y = anno,
        by = c("smiles" = "feature")
    ) %>%
    dplyr::filter(sim > 0.20)

tani.comb %>%
    dplyr::select(target_name, sim, target_compound_id, smiles, peak_id, adduct)

hmp.p <- tani.comb %>%
    dplyr::mutate(apigenin = "Apigenin") %>%
    dplyr::distinct(target_inchi, .keep_all = TRUE) %>%
    dplyr::mutate(target_name = forcats::fct_reorder(
        .f = target_name,
        .x = sim,
        .fun = "mean",
        .desc = FALSE
    )) %>%
    ggplot(.,
           aes(
               x = apigenin,
               y = target_name,
               fill = sim
           )) +
    geom_tile() +
    scale_y_discrete(expand = expansion(c(0, 0))) +
    scale_x_discrete(expand = expansion(c(0, 0))) +
    scale_fill_gradient() +
    theme_cowplot() +
    theme(
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.7)
    ) +
    labs(
        fill = "Similarity",
        title = "Structural similarity to Apigenin"
    )
hmp.p

tani.comb %>%
    dplyr::filter(peak_id %in% c("FT04480")) %>%
    dplyr::select(smiles, target_name)

# doesnt work
# rcdk::view.molecule.2d(parse.smiles("C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O")[[1]])

mol.2d <- parse.smiles("C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O")[[1]]
img <- view.image.2d(mol.2d)
img.grob <- grid::rasterGrob(img, interpolate = TRUE)
img.p <- wrap_elements(img.grob, clip = TRUE)

img.p / hmp.p +
    plot_layout(
        heights = c(1, 5)
    )

apiin.sim <- findSimilarity(
    smiles.obj = smiles,
    query.smiles = "C1[C@@]([C@H]([C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H](O[C@H]2OC3=CC(=C4C(=C3)OC(=CC4=O)C5=CC=C(C=C5)O)O)CO)O)O)O)(CO)O",
    anno.obj = anno,
    type = "circular",
    method = "tanimoto"
)

apigenin.sim <- findSimilarity(
    smiles.obj = smiles,
    query.smiles = "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O",
    anno.obj = anno,
    type = "circular",
    method = "tanimoto"
)

p1 <- apigenin.sim %>%
    dplyr::filter(sim > 0.2) %>%
    dplyr::distinct(target_inchi, .keep_all = TRUE) %>%
    dplyr::mutate(target_name = forcats::fct_reorder(
        .f = target_name,
        .x = sim,
        .fun = "mean",
        .desc = FALSE
    )) %>%
    ggplot(.,
           aes(
               x = "Apigenin",
               y = target_name,
               fill = sim
           )) +
    geom_tile() +
    scale_y_discrete(expand = expansion(c(0, 0))) +
    scale_x_discrete(expand = expansion(c(0, 0))) +
    scale_fill_gradient(
        limits = c(0.2, 0.7)
    ) +
    theme_cowplot() +
    theme(
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.7)
    ) +
    labs(
        fill = "Tanimoto similarity",
        title = "Structural similarity to Apigenin"
    )

p2 <- apiin.sim %>%
    dplyr::filter(sim > 0.2) %>%
    dplyr::distinct(target_inchi, .keep_all = TRUE) %>%
    dplyr::mutate(target_name = forcats::fct_reorder(
        .f = target_name,
        .x = sim,
        .fun = "mean",
        .desc = FALSE
    )) %>%
    ggplot(.,
           aes(
               x = "Apiin",
               y = target_name,
               fill = sim
           )) +
    geom_tile() +
    scale_y_discrete(expand = expansion(c(0, 0))) +
    scale_x_discrete(expand = expansion(c(0, 0))) +
    scale_fill_gradient(
        limits = c(0.2, 0.7)
    ) +
    theme_cowplot() +
    theme(
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.7)
    ) +
    labs(
        fill = "Tanimoto similarity",
        title = "Structural similarity to Apiin"
    )

p1 + p2 +
    plot_layout(
        guides = "collect",
        axes = "collect"
    )

check.peaks <- unique(
    c(
        apiin.sim %>% 
            dplyr::filter(sim > 0.2) %>%
            pull(peak_id), 
        apigenin.sim %>% 
            dplyr::filter(sim > 0.2) %>%
            pull(peak_id)
        ))

saved.feats <- c()
for (i in check.peaks) {
    anno.peak.chr <- plotFeatChrInt(
        feature_chrom = anno.chrs,
        feature = i,
        method = "sum",
        value = "into",
        filled = TRUE,
        missing = "rowmin_half",
        msLevel = 1,
        save_loc = NULL,
        device = NULL,
        feat_pairs = FALSE
    ) 
    print(anno.peak.chr$combined)
    
    saved.feat <- readline("1 for save, 2 for throw away: ")
    saved.feat <- setNames(saved.feat, unique(anno.peak.chr$p2_data$feature))
    saved.feats <- c(saved.feats, saved.feat)
}

saved.feat.names <- names(saved.feats[saved.feats == "1"])

fp1 <- apigenin.sim %>%
    dplyr::filter(sim > 0.20) %>%
    dplyr::distinct(target_inchi, .keep_all = TRUE) %>%
    dplyr::filter(peak_id %in% saved.feat.names) %>%
    dplyr::mutate(name_feat = paste0(target_name, "_", peak_id)) %>%
    dplyr::mutate(name_feat = forcats::fct_reorder(
        .f = name_feat,
        .x = sim,
        .fun = "mean",
        .desc = FALSE
    )) %>%
    ggplot(.,
           aes(
               x = "Apigenin",
               y = name_feat,
               fill = sim
           )) +
    geom_tile() +
    scale_y_discrete(expand = expansion(c(0, 0))) +
    scale_x_discrete(expand = expansion(c(0, 0))) +
    scale_fill_gradient(
        limits = c(0.2, 0.7)
    ) +
    theme_cowplot() +
    theme(
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.7)
    ) +
    labs(
        fill = "Tanimoto similarity",
        title = "Structural similarity to Apigenin"
    )

fp2 <- apiin.sim %>%
    dplyr::filter(sim > 0.20) %>%
    dplyr::distinct(target_inchi, .keep_all = TRUE) %>%
    dplyr::filter(peak_id %in% saved.feat.names) %>%
    dplyr::mutate(name_feat = paste0(target_name, "_", peak_id)) %>%
    dplyr::mutate(name_feat = forcats::fct_reorder(
        .f = name_feat,
        .x = sim,
        .fun = "mean",
        .desc = FALSE
    )) %>%
    ggplot(.,
           aes(
               x = "Apiin",
               y = name_feat,
               fill = sim
           )) +
    geom_tile() +
    scale_y_discrete(expand = expansion(c(0, 0))) +
    scale_x_discrete(expand = expansion(c(0, 0))) +
    scale_fill_gradient(
        limits = c(0.2, 0.7)
    ) +
    theme_cowplot() +
    theme(
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.7)
    ) +
    labs(
        fill = "Tanimoto similarity",
        title = "Structural similarity to Apiin"
    )

fp1 + fp2 +
    plot_layout(
        guides = "collect",
        axes = "collect"
    )

for (i in saved.feat.names) {
  try <- plotFeatChrInt(
    feature_chrom = anno.chrs,
    feature = i,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = "rowmin_half",
    msLevel = 1,
    save_loc = "/chem_sim/",
    device = "svg",
    feat_pairs = FALSE
  )
}
