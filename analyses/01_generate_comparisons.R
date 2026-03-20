cli::cli_h1(basename(this.path::this.path()))
# ==============================================================================
# Generating comparisons -------------------------------------------------------
# ==============================================================================
# Find the exact index of the groups, regardless of order
int_upset_ids <- purrr::map_vec(
  .x = int_upset_comps,
  .f = ~ extract_upset_id(upset_distinct, .x)
)

upset_distinct_comp <- upset_distinct_comps[c(int_upset_ids)]

int_sig_diff <- unname(unlist(upset_distinct_comp))
