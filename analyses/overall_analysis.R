################################################################################
# Counting interesting metabolites ---------------------------------------------
################################################################################

# How many are interesting
# This is incorrect because this is done on filtered data for now
xchr9_full_ints <- all_xchr9_data %>%
  # dplyr::left_join(
  #   x = .,
  #   y = dplyr::select(meta, sample, group, path),
  #   by = c("name" = "sample")
  # ) %>%
  # dplyr::relocate(group, .after = "feature") %>%
  dplyr::mutate(
    title = paste0(
      feature, " ",
      round(mzmed, 2), " ",
      round(rtmed, 2)
    )
  )

source("metadata_prep.R")

# # Took the metadata comparisons thing
summarized_exps <- experiments %>%
  dplyr::ungroup() %>%
  dplyr::select(
    c(
      "experiment_id",
      "experiment",
      "strain",
      "condition",
      "unclean_strain",
      "unclean_condition"
    )
  ) %>%
  dplyr::distinct(experiment_id, .keep_all = TRUE)

total_ints <- xchr9_full_ints %>%
  dplyr::group_by(experiment) %>%
  dplyr::summarize(n_features = dplyr::n_distinct(feature)) %>%
  dplyr::arrange(desc(n_features)) %>%
  dplyr::left_join(
    x = .,
    y = summarized_exps,
    by = "experiment"
  ) %>%
  dplyr::left_join(
    x = .,
    y = glycone_pairs_metadata,
    by = c("unclean_condition" = "glycoside")
  )

total_ints %>%
  ggplot(
    aes(
      x = n_features,
      y = strain,
      color = condition
    )
  ) +
  geom_point() +
  scale_x_continuous(transform = "log10")

total_ints %>%
  ggplot(
    aes(
      x = n_features,
      y = strain,
      color = aglycone
    )
  ) +
  geom_point() +
  scale_x_continuous(transform = "log10")

total_ints %>%
  ggplot(
    aes(
      x = n_features,
      y = condition,
      color = strain
    )
  ) +
  geom_point() +
  scale_x_continuous(transform = "log10")

total_ints %>%
  ggplot(
    aes(
      x = n_features,
      y = aglycone,
      color = strain
    )
  ) +
  geom_point() +
  scale_x_continuous(transform = "log10")


