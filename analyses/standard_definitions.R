library(tidyverse)
source("scripts/functions.R")

for (i in c("glycoside", "aglycone")) {
  top_summary_path <- paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/standard/",
    i,
    "/tables/",
    "top_standard_summary.csv"
  )

  top_summary <- readr::read_csv(
    file = top_summary_path,
    progress = FALSE,
    show_col_types = FALSE
  )

  full_peak_table_path <- paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/standard/",
    i,
    "/tables/",
    "full_peak_table.csv"
  )

  full_peak_table <- readr::read_csv(
    file = full_peak_table_path,
    progress = FALSE,
    show_col_types = FALSE
  )

  chr_vals_path <- paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/standard/",
    i,
    "/tables/",
    "chromatogram_values.csv"
  )

  chr_vals <- readr::read_csv(
    file = chr_vals_path,
    progress = FALSE,
    show_col_types = FALSE
  )

  rt_ref <- top_summary %>%
    dplyr::distinct(group, rt) %>%
    dplyr::rename(rt_ref = rt)

  # looks good
  all_summary <- full_peak_table %>%
    dplyr::filter(maxo > 10000) %>%
    dplyr::left_join(rt_ref, by = "group") %>%
    dplyr::filter(abs(rt - rt_ref) <= 0.5) %>%
    dplyr::select(-rt_ref) %>%
    dplyr::group_by(group, adduct) %>%
    dplyr::summarize(
      median_into = median(into),
      median_maxo = median(maxo),
      rt = median(rt),
      mz = median(mz),
      rtmin = median(rtmin),
      rtmax = median(rtmax),
      mzmin = median(mzmin),
      mzmax = median(mzmax),
      .groups = "drop"
    ) %>%
    dplyr::arrange(group, desc(median_into))

  chr_vals2 <- chr_vals %>%
    dplyr::semi_join(all_summary, by = c("group", "adduct")) %>%
    dplyr::left_join(
      x = .,
      y = all_summary %>%
        dplyr::select(group, adduct, rt_center = rt),
      by = c("group", "adduct")
    ) %>%
    dplyr::filter(abs(rt - rt_center) <= 8) %>%
    dplyr::select(-rt_center) %>%
    dplyr::group_by(group, adduct) %>%
    dplyr::filter(sample == sample[which.max(intensity)]) %>%
    dplyr::ungroup()

  chr_vals2 %>%
    dplyr::filter(group == "Apiin") %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = rt,
        y = intensity,
        group = interaction(sample, adduct),
        color = adduct
      )
    ) +
    ggplot2::geom_line() +
    ggplot2::theme_minimal()

  adduct_colors <- c(
    "[2M-H]-" = "#aab7d9",
    "[2M+C2H3O2]-" = "#85c4bc",
    "[2M+CHO2]-" = "#dec49f",
    "[3M-H]-" = "#cae3e9",
    "[M-2H]2-" = "#a6cca7",
    "[M-H]-" = "#e7b8b7",
    "[M-H+HCOONa]-" = "#8bd0eb",
    "[M]-" = "#ddefce",
    "[M+C2F3O2]-" = "#c2ceab",
    "[M+C2H3N-H]-" = "#a9e7e4",
    "[M+CHO2]-" = "#bec0ac",
    "[M+Cl]-" = "#a6c8ce",
    "[M+K-2H]-" = "#9db7b1",
    "[M+Na-2H]-" = "#d8bfe3"
  )


  save_graph_path <- paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/standard/",
    i,
    "/peaks_manual"
  )
  dir.create(save_graph_path, FALSE, TRUE)
  for (j in unique(chr_vals2$group)) {
    tmp_p <- chr_vals2 %>%
      dplyr::filter(group == j) %>%
      dplyr::mutate(
        adduct = forcats::fct_reorder(adduct, intensity, .fun = max),
        scaled_int = 30 * intensity / max(intensity)
      ) %>%
      ggplot2::ggplot(
        ggplot2::aes(
          x = rt,
          y = adduct,
          height = scaled_int,
          group = interaction(adduct, sample),
          fill = adduct
        )
      ) +
      ggridges::geom_ridgeline(scale = 1, alpha = 0.7) +
      ggplot2::scale_fill_manual(values = adduct_colors, drop = TRUE) +
      ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
      ggplot2::scale_x_continuous(expand = ggplot2::expansion(c(0, 0))) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank()
      ) +
      ggplot2::labs(
        x = "Retention time (s)",
        fill = "Adduct"
      )

    ggplot2::ggsave(
      filename = file.path(save_graph_path, paste0(j, ".svg")),
      plot = tmp_p,
      device = "svg",
      height = 6,
      width = 8,
      units = "in"
    )
  }
}

