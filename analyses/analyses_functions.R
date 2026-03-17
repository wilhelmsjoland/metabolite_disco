read_snakes <- function(
  path = NULL
) {
  list2env(
    readRDS(
      file.path(path)
    ),
    envir = .GlobalEnv
  )
}

mock_snakemake <- function(rule = NULL, config_file = NULL) {

  setClass("Snakemake", slots = list(
    input  = "list",
    output = "list",
    params = "list",
    config = "list",
    rule   = "character"
  ))

  config <- yaml::read_yaml(config_file)
  out <- config[["output"]]
  sm <- function(...) file.path(out, "snakemake_objects", ...)

  # Common params shared by all rules
  common <- list(
    output = out,
    seed   = config[["seed"]],
    cores  = config[["cores"]]
  )

  rules <- list(
    setup = list(
      input  = list(),
      output = list(sm("01_setup.rds")),
      params = c(common, list(
        data_path = config[["data_path"]],
        meta_file = config[["meta_file"]]
      ))
    ),
    bpc = list(
      input  = list(sm("01_setup.rds")),
      output = list(sm("02_bpc.rds")),
      params = common
    ),
    internal_standard = list(
      input  = list(
        setup = sm("01_setup.rds"),
        bpc   = sm("02_bpc.rds")
      ),
      output = list(sm("03_internal_standard.rds")),
      params = c(common, list(
        internal_standard = config[["internal_standard"]],
        is_adduct         = config[["is_adduct"]],
        ppm_global        = config[["ppm_global"]],
        sn_threshold      = config[["sn_threshold"]]
      ))
    ),
    peak_calling = list(
      input  = list(
        setup        = sm("01_setup.rds"),
        internal_std = sm("03_internal_standard.rds")
      ),
      output = list(sm("04_peak_calling.rds")),
      params = c(common, list(
        ppm_global         = config[["ppm_global"]],
        sn_threshold       = config[["sn_threshold"]],
        mzdiff             = config[["mzdiff"]],
        beta_cor_threshold = config[["beta_cor_threshold"]],
        beta_snr_threshold = config[["beta_snr_threshold"]]
      ))
    ),
    alignment = list(
      input  = list(
        setup        = sm("01_setup.rds"),
        bpc          = sm("02_bpc.rds"),
        internal_std = sm("03_internal_standard.rds"),
        peak_calling = sm("04_peak_calling.rds")
      ),
      output = list(sm("05_alignment.rds")),
      params = c(common, list(
        ppm_global         = config[["ppm_global"]],
        bw_first_grouping  = config[["bw_first_grouping"]],
        peak_anchor_sd     = config[["peak_anchor_sd"]],
        min_fraction_align = config[["min_fraction_align"]],
        extra_peaks        = config[["extra_peaks"]],
        span               = config[["span"]]
      ))
    ),
    correspondence = list(
      input  = list(
        alignment    = sm("05_alignment.rds"),
        internal_std = sm("03_internal_standard.rds")
      ),
      output = list(sm("06_correspondence.rds")),
      params = c(common, list(
        ppm_global        = config[["ppm_global"]],
        bw_first_grouping = config[["bw_first_grouping"]],
        bw_second_grouping = config[["bw_second_grouping"]]
      ))
    ),
    gap_filling = list(
      input  = list(
        correspondence = sm("06_correspondence.rds")
      ),
      output = list(sm("07_gap_filling.rds")),
      params = common
    ),
    filter_features = list(
      input  = list(
        gap_filling = sm("07_gap_filling.rds")
      ),
      output = list(sm("08_filter_features.rds")),
      params = c(common, list(
        missingness        = config[["missingness"]],
        sn_threshold       = config[["sn_threshold"]],
        beta_cor_threshold = config[["beta_cor_threshold"]],
        beta_snr_threshold = config[["beta_snr_threshold"]]
      ))
    ),
    scaling = list(
      input  = list(
        filter_features = sm("08_filter_features.rds")
      ),
      output = list(sm("09_scaling.rds")),
      params = common
    ),
    limma = list(
      input  = list(
        scaling = sm("09_scaling.rds"),
        setup   = sm("01_setup.rds")
      ),
      output = list(sm("10_limma.rds")),
      params = c(common, list(
        gap_filling = config[["gap_filling"]]
      ))
    ),
    pca = list(
      input  = list(
        limma = sm("10_limma.rds"),
        setup = sm("01_setup.rds")
      ),
      output = list(sm("11_pca.rds")),
      params = common
    ),
    volcano = list(
      input  = list(
        limma = sm("10_limma.rds")
      ),
      output = list(sm("12_volcano.rds")),
      params = c(common, list(
        qvalue      = config[["qvalue"]],
        gap_filling = config[["gap_filling"]]
      ))
    ),
    upset = list(
      input  = list(
        limma = sm("10_limma.rds")
      ),
      output = list(sm("13_upset.rds")),
      params = c(common, list(
        qvalue      = config[["qvalue"]],
        gap_filling = config[["gap_filling"]]
      ))
    ),
    intersecting_features = list(
      input  = list(
        upset = sm("13_upset.rds")
      ),
      output = list(sm("14_intersecting_features.rds")),
      params = common
    ),
    prep_annotation_biotransformation = list(
      input  = list(
        filter_features = sm("08_filter_features.rds"),
        limma           = sm("10_limma.rds")
      ),
      output = list(sm("15_prep_annotation_biotransformation.rds")),
      params = c(common, list(
        data_path      = config[["data_path"]],
        polarity       = config[["polarity"]],
        biotransf_file = config[["biotransf_file"]],
        rpairs_path    = config[["rpairs_path"]],
        qvalue         = config[["qvalue"]]
      ))
    ),
    mz_predictions = list(
      input  = list(
        limma     = sm("10_limma.rds"),
        prep_biot = sm("15_prep_annotation_biotransformation.rds")
      ),
      output = list(sm("16_mz_predictions.rds")),
      params = c(common, list(
        glycoside      = config[["glycoside"]],
        aglycone       = config[["aglycone"]],
        polarity       = config[["polarity"]],
        biotransf_file = config[["biotransf_file"]],
        ppm_match      = config[["ppm_match"]],
        all_vs_all     = config[["all_vs_all"]]
      ))
    ),
    annotation = list(
      input  = list(
        limma     = sm("10_limma.rds"),
        prep_biot = sm("15_prep_annotation_biotransformation.rds")
      ),
      output = list(sm("17_annotation.rds")),
      params = c(common, list(
        polarity  = config[["polarity"]],
        ppm_match = config[["ppm_match"]]
      ))
    ),
    biotransformer = list(
      input  = list(
        filter_features = sm("08_filter_features.rds")
      ),
      output = list(sm("18_biotransformer.rds")),
      params = c(common, list(
        smiles   = config[["smiles"]],
        biot_dir = config[["biot_dir"]],
        polarity = config[["polarity"]],
        ppm_match = config[["ppm_match"]]
      ))
    ),
    molecular_similarity = list(
      input  = list(
        annotation     = sm("17_annotation.rds"),
        biotransformer = sm("18_biotransformer.rds")
      ),
      output = list(sm("19_molecular_similarity.rds")),
      params = c(common, list(
        smiles = config[["smiles"]]
      ))
    )
  )

  if (!rule %in% names(rules)) {
    stop(
      "Unknown rule: '", rule, "'. Available rules: ",
      paste(names(rules), collapse = ", ")
    )
  }

  r <- rules[[rule]]
  new("Snakemake",
    input  = r$input,
    output = r$output,
    params = r$params,
    config = config,
    rule   = paste0("mock.", rule, ".R")
  )
}
