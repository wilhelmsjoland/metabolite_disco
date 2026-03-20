# =============================================================================
# Molecular similarity ---------------------------------------------------------
# =============================================================================
# from the heatmap
int_bio_fts <- c(
  "FT11675",
  "FT08181",
  "FT08191",
  "FT04076",
  "FT05019",
  "FT05465",
  "FT02089",
  "FT04760",
  "FT01640",
  "FT06679"
)

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

# missing smiles
anno$peak_id[is.na(anno$target_smiles)]

# 648 rows
# 1,592
anno_filt <- anno %>%
  # TODO
  # DO A more intelligent filtering than this based on
  # how much information is available in all the rows
  dplyr::group_by(adduct) %>% # Think this might have fixed the loss of same
  dplyr::distinct(target_inchikey, .keep_all = TRUE) %>%
  dplyr::filter(peak_id %in% int_annos) # int.bio.feats

smiles <- anno_filt$target_smiles
names(smiles) <- anno_filt$feature

apigenin_smiles <- "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O"
apiin_smiles <- paste0(
  "C1[C@@]([C@H]([C@@H](O1)O[C@@H]2[C@H]([C@@H]([C@H]",
  "(O[C@H]2OC3=CC(=C4C(=C3)OC(=CC4=O)C5=CC=C(C=C5)O)O)CO)O)O)O)(CO)O"
)

apigenin_sims <- mol_similarity(
  query_smiles = apigenin_smiles,
  target_smiles = smiles,
  kekulise = TRUE, # parsing incorrect smiles with electrons
  omit_nulls = TRUE,
  fingerprint = "circular",
  circular_type = "ECFP6",
  method = "tanimoto"
) %>%
  dplyr::left_join(
    x = .,
    y = anno,
    by = "feature"
  ) %>%
  dplyr::filter(sim > 0.15)

# Get names for individual inchikey
inchi_ks <- paste0(apigenin_sims$target_inchikey, collapse = ",")
cmd <- paste0(
  "curl -s -d \"inchikey=", inchi_ks, "\" ",
  "\"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/property/",
  "Title,MolecularFormula,InChIKey,ExactMass,XLogP,TPSA/",
  "CSV\""
)
results <- system(cmd, intern = TRUE)
inchikey_map <- read_csv(file = paste0(results, collapse = "\n"))

apigenin_sims %>%
  dplyr::left_join(
    x = .,
    y = inchikey_map,
    by = c("target_inchikey" = "InChIKey")
  ) %>%
  dplyr::mutate(Title = forcats::fct_reorder( # target_name
    .f = Title, # target_name
    .x = sim,
    .fun = "mean",
    .desc = TRUE
  )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = Title, # target_name
      y = sim
    )
  ) +
  # because there are duplicates - just choose the best one
  ggplot2::geom_col(
    aes(fill = peak_id),
    stat = "summary",
    fun = "max",
    color = "black",
    position = ggplot2::position_dodge()
  ) +
  ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
  ggplot2::scale_y_continuous(
    expand = ggplot2::expansion(c(0, 0)),
    limits = c(0, 0.8)
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(
    axis.title.x = ggplot2::element_blank(),
    axis.title.y = ggplot2::element_text(angle = -90),
    legend.title = ggplot2::element_blank()
  ) +
  ggplot2::labs(y = "Tanimoto similarity")
