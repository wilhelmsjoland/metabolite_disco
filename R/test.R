test <- filt.match.diffs2 %>%
  dplyr::mutate(obs_ppm = num_to_ppm(abs(delta_mass - obs_delta_mass))) %>% 
  dplyr::filter(feat1 != feat2) %>%
  dplyr::arrange(obs_ppm) %>%
  dplyr::slice(1:100)

test.cols <- unique(c(test$feat1, test$feat2))

# test %>%
#   view()

# TODO
# now left join the rpairs numbers to the rpair metadata
# also filter to only features that are significantly different

test.chrs <- xcms::featureChromatograms(
  object = xchr9,
  expandRt = 0,
  expandMz = 0,
  aggregationFun = "sum",
  filled = TRUE,
  features = test.cols,
  missing = 0,
  return.type = "XChromatograms"
)

for (i in test.cols) {
  tmp.test.chr <- plotFeatChrInt(
    feature_chrom = test.chrs,
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
  print(tmp.test.chr$combined)
  x <- readline("Enter for next, break for stop: ")
  if (x == "break") [
    break
  ]
}
