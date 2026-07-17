library(tidyverse)
library(xcms)

exps_for_tic <- list.files(
  path = paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/old/old_windows_output",
    "/experiment"
  ),
  include.dirs = TRUE,
  recursive = TRUE,
  full.names = TRUE,
  pattern = "xchr.rds"
)

summed_tics <- exps_for_tic %>%
  purrr::map(
    .x = .,
    .f = ~ {

      xchr <- readRDS(file = .x)

      xchr_peaks <- tibble::as_tibble(
        xcms::chromPeaks(xchr),
        rownames = "peak"
      )

      xchr_meta <- tibble::as_tibble(
        MsExperiment::sampleData(xchr),
        rownames = "sample_id"
      ) %>%
        dplyr::mutate(sample = dplyr::row_number())

      xchr_data_comb <- dplyr::left_join(
        x = xchr_peaks,
        y = xchr_meta,
        by = "sample"
      ) %>%
        dplyr::mutate(
          dplyr::across(
            .cols = c("into", "intb"),
            .fns = ~ log2(.)
          )
        ) %>%
        dplyr::mutate(comb = paste0(sample, "_", sample_id)) %>%
        dplyr::mutate(
          comb = forcats::fct_reorder(
            .f = comb,
            .x = sample
          )
        )

      xchr_box_summary <- xchr_data_comb %>%
        dplyr::group_by(comb, group) %>%
        dplyr::summarize(
          n = dplyr::n(),
          stats = list(boxplot.stats(into)),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          ymin = purrr::map_dbl(stats, ~ .x$stats[1]), # lower whisker
          lower = purrr::map_dbl(stats, ~ .x$stats[2]), # Q1 (hinge)
          median = purrr::map_dbl(stats, ~ .x$stats[3]),
          upper = purrr::map_dbl(stats, ~ .x$stats[4]), # Q3 (hinge)
          ymax = purrr::map_dbl(stats, ~ .x$stats[5]), # upper whisker
          iqr = upper - lower,
          # outliers = purrr::map(stats, "out") # list, one vector per box
        ) %>%
        dplyr::select(-stats) %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(
          peak_ratio = n / stats::median(n), # centered within its own group
          int_ratio  = median / stats::median(median)
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(
          # Euclidean distance from the "ideal" point (int_ratio = 1,
          # peak_ratio = 1) - a single filterable number combining both
          # ratios, instead of eyeballing the 2D scatter for each sample.
          ratio_dist = sqrt((1 - int_ratio)^2 + (1 - peak_ratio)^2)
        ) %>%
        dplyr::mutate(
          # direct % of own-group median - no spread estimate needed, so it
          # stays stable at n = 4 replicates per group (unlike the MAD/z-score
          # version this replaced, which fell apart at this sample size)
          peak_stat_group = factor(
            dplyr::if_else(
              peak_ratio < 0.65, # | int_ratio < 0.7,
              "bad",
              "good"
            ),
            levels = c("good", "bad")
          )
        ) %>%
        # same scaling geom_boxplot uses for varwidth
        dplyr::mutate(width = 0.9 * sqrt(n) / max(sqrt(n))) %>%
        dplyr::mutate(
          experiment = basename(dirname(dirname(.x))),
          .before = "comb"
        )
    }
  ) %>%
  dplyr::bind_rows()

# Base the cutoff on this
summed_tics %>%
  dplyr::filter(
    experiment %in% basename(dirname(dirname(exps_for_tic)))
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = int_ratio,
      y = peak_ratio
    )
  ) +
  ggplot2::geom_text(
    ggplot2::aes(label = stringr::str_split_i(comb, "_", 1))
  ) +
  ggplot2::scale_y_continuous(breaks = seq(0, 2, by = 0.1)) +
  ggplot2::theme_bw()

# Redundant - the last graph is better for deciding cutoff
# summed_tics %>%
#   dplyr::mutate(row_number = dplyr::row_number()) %>%
#   ggplot2::ggplot(
#     ggplot2::aes(
#       x = row_number,
#       y = ratio_dist
#     )
#   ) +
#   ggplot2::geom_point() +
#   ggplot2::theme_bw() +
#   ggplot2::labs(
#     x = "Sample number",
#     y = "Distance from perfect agreement between samples"
#   )

## graphs - to check each experiment and if filtering was successful
review_flag <- c()
for (i in basename(dirname(dirname(exps_for_tic)))) {
  tmp_p <- summed_tics %>%
    dplyr::filter(experiment %in% i) %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = comb,
        fill = group
      )
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = comb,
        xend = comb,
        color = peak_stat_group,
        y = ymin,
        yend = ymax
      ),
      linewidth = 6
    ) +
    ggplot2::geom_boxplot(
      ggplot2::aes(
        ymin = ymin,
        lower = lower,
        middle = median,
        upper = upper,
        ymax = ymax,
        width = width
      ),
      stat = "identity"
    ) +
    ggplot2::scale_color_manual(values = c(NA, "firebrick")) +
    ggplot2::guides(x = ggplot2::guide_axis(angle = -45)) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::labs(title = i)

  print(tmp_p)

  tmp_review_flag <- readline("Flag for review: 1, Enter for next: ")
  tmp_review_flag <- setNames(tmp_review_flag, i)

  if (tmp_review_flag == "1") {
    review_flag <- c(review_flag, tmp_review_flag)
  }
}

review_flag

# Inspect bad sample tibble
summed_tics %>%
  dplyr::filter(peak_stat_group == "bad")

# Amount of unique samples to remove
summed_tics %>%
  dplyr::filter(peak_ratio < 0.65) %>%
  dplyr::arrange(comb) %>%
  pull(comb) %>%
  stringr::str_split_i(., "_", 2) %>%
  unique() %>%
  length()

bad_samples <- summed_tics %>%
  dplyr::filter(peak_stat_group == "bad") %>%
  dplyr::mutate(comb = stringr::str_split_i(comb, "_", 2)) %>%
  dplyr::pull(comb) %>%
  unique()

for (i in bad_samples) {
  cat(i, "\n")
}
