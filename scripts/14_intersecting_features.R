cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Finding intersecting features ------------------------------------------------
# ==============================================================================
cli::cli_h3("Finding intersecting features")

# Create all comparisons
combinations <- unlist(
  lapply(
    X = seq_along(comparisons),
    FUN = function(x) {
      combn(comparisons, x, simplify = FALSE)
    }
  ),
  recursive = FALSE
)

upset_comps <- list()
for (i in seq_along(combinations)) {
  tmp_intersect_feats <- find_intersect_feat(
    data = upset_tib,
    set = combinations[[i]],
    full_set = comparisons
  )

  upset_comps[[
    stringr::str_flatten(combinations[[i]], collapse = "*")
  ]] <- tmp_intersect_feats$feature
}
cli::cli_alert_success(
  paste0(
    "Generated intersecting features"
  )
)

# Only temporary
upset_comp <- upset_comps[[
  paste0(
    "bu_mutant_apiin-bu_wt_apiin*",
    "bu_mutant_control-bu_wt_apiin*",
    "bu_wt_apiin-bu_wt_control"
  )
]]