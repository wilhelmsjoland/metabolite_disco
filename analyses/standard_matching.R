library(tidyverse)

glycone_feats <- purrr::map(
  .x = c("aglycone", "glycoside"),
  .f = ~ {
    path <- paste0(
      "/Volumes/bluecub/aglycone_release_100um_24h/output/standard/",
      .x,
      "/tables/top_standard_summary.csv"
    )
    readr::read_csv(file = path) %>%
      dplyr::mutate(source = .x, .before = "group")
  }
) %>%
  dplyr::bind_rows()

dirs <- list.dirs(
  path = "/Volumes/bluecub/aglycone_release_100um_24h/output/experiment",
  recursive = FALSE
)

standard_map <- purrr::map(
  .x = dirs,
  .f = ~ {
    path <- paste0(.x, "/objects/xchr9.rds")
    if (!file.exists(path)) return(NULL)
    xchr9 <- readRDS(path)
    query <- S4Vectors::DataFrame(glycone_feats)
    target <- xcms::featureDefinitions(xchr9) %>%
    tibble::rownames_to_column(var = "feature") %>%
    dplyr::select(
      c(
        "feature",
        "mzmed",
        "mzmin",
        "mzmax",
        "rtmed",
        "rtmin",
        "rtmax"
      )
    )

    param <- MetaboAnnotation::MzRtParam(
      ppm = 25,
      toleranceRt = 20
    )

    matches <- MetaboAnnotation::matchValues(
      query = query,
      target = target,
      param = param,
      mzColname = c("mz", "mzmed"),
      rtColname = c("rt", "rtmed")
    )

    matched_data <- MetaboAnnotation::matchedData(matches) %>%
      tibble::as_tibble() %>%
      dplyr::mutate(experiment = basename(.x), .before = "source")
  }
) %>%
  dplyr::bind_rows()

final_standard_map <- standard_map %>%
  tidyr::drop_na(target_feature)


final_standard_vals <- purrr::map(
  .x = dirs,
  .f = ~ {
    path <- paste0(.x, "/objects/xchr9.rds")
    if (!file.exists(path)) return(NULL)
    readRDS(path) %>%
      xcms::featureValues(
        object = .,
        method = "sum",
        value = "into",
        intensity = "into",
        filled = TRUE,
        missing = 0,
        msLevel = 1L
      ) %>%
      tibble::as_tibble(., rownames = "feature") %>%
      dplyr::rename_with(
        .cols = dplyr::contains(".mzML"),
        .fn = ~ sub(".*\\\\", "", .x)
      ) %>%
      tidyr::pivot_longer(cols = dplyr::contains(".mzML")) %>%
      dplyr::left_join(
        x = .,
        y = readr::read_csv(
          file = paste0(.x, "/tables/metadata.csv"),
          progres = FALSE,
          show_col_types = FALSE
        ) %>%
          dplyr::mutate(name = basename(path), .before = "group"),
        by = "name"
      ) %>%
      dplyr::inner_join(
        x = .,
        y = final_standard_map %>%
          dplyr::filter(experiment == basename(.x)) %>%
          dplyr::select(target_feature, group, adduct),
        by = c("feature" = "target_feature")
      )
  }
) %>%
  dplyr::bind_rows()

final_standard_p <- final_standard_vals %>%
  dplyr::filter(condition != "ycfa_glucose") %>%
  dplyr::mutate(
    short_group_x = sub(
      "_[^_]+$",
      "",
      group.x
    )
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    comb = paste0(
      unclean_condition, " -> ", group.y, " ", adduct
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    strain = forcats::fct_relevel(
      strain,
      c(
        "b_uniformis_atcc_8492",
        "bu_gsh_d_ggh_c_gsh_g",
        "b_ovatus_atcc_8483",
        "b_ovatus_atcc_8483_d_operon",
        "b_thetaiotaomicron_vpi_5482",
        "p_copri_i_ak263"
      )
    )
  ) %>%
  dplyr::mutate(
    unclean_strain = forcats::fct_reorder(
      unclean_strain,
      as.numeric(factor(strain))
    )
  ) %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x = unclean_strain,
      y = comb,
      fill = value
    )
  ) +
  ggplot2::geom_tile() +
  guides(x = ggplot2::guide_axis(angle = -45)) +
  ggplot2::theme_minimal() +
  ggplot2::theme(
    axis.title = ggplot2::element_blank()
  ) +
  ggplot2::labs(
    fill = "Peak area"
  )

final_standard_p

ggplot2::ggsave(
  filename = file.path(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/standard_matching/",
    "final_standard_plot.pdf"
  ),
  plot = final_standard_p,
  device = "pdf",
  height = 8,
  width = 8,
  units = "in"
)




glycone_feats %>%
  filter(group == "Afzelin")

featureDefinitions(tester) %>%
  as_tibble(., rownames = "feature") %>%
  dplyr::filter(dplyr::between(rtmed, 248, 249))

# Extract ion chromatogram for mz 431 from the raw data
xcms::chromatogram(tester, mz = c(430.5, 431.5), rt = c(230, 260))

final_standard_vals %>%
  dplyr::filter(condition != "ycfa_glucose") %>%
  dplyr::mutate(
    short_group_x = sub(
      "_[^_]+$",
      "",
      group.x
    )
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    comb = paste0(
      unclean_condition, " -> ", group.y, " ", adduct
    )
  ) %>%
  ungroup() %>%
  dplyr::mutate(
    unclean_condition = stringr::str_to_lower(unclean_condition),
    group.y = stringr::str_to_lower(group.y)
  ) %>%
  dplyr::filter(unclean_condition == group.y)

final_standard_vals %>%
  distinct(unclean_condition, group.y) %>%
  print(n = "all")

featureDefinitions(tester) %>%
  as_tibble(., rownames = "feature") %>%
  dplyr::filter(dplyr::between(mzmed, 431, 432)) %>%
  dplyr::arrange(rtmed) %>%
  print(n = "all")

glycone_feats %>% filter(group == "Afzelin")

tester <- readRDS(
  paste0(
    "/Volumes/bluecub/aglycone_release_100um_24h/output/experiment/",
    "afzelin_b_ovatus_atcc_8483_and_b_ovatus_atcc_8483_d_operon",
    "/objects/xchr9.rds"
  )
)

featureDefinitions(tester) %>%
  as_tibble(., rownames = "feature") %>%
  dplyr::filter(
    dplyr::between(mzmed, 430.5, 431.5),
    dplyr::between(rtmed, 245, 251)
  )

sp <- spectra(tester)
sp@backend@spectraData$dataStorage <- gsub(
  r"(V:\aglycone_release_100um_24h\data\mzml_files\)",
  "/Volumes/bluecub/aglycone_release_100um_24h/data/experiment/mzml/",
  sp@backend@spectraData$dataStorage,
  fixed = TRUE
)
tester@spectra <- sp
tester@spectra <- Spectra::setBackend(
  spectra(tester),
  MsBackendMemory()
)

sd <- MsExperiment::sampleData(tester)
sd$path <- meta$path
rownames(sd) <- meta$sample
sd$spectraOrigin <- gsub(
  r"(V:\aglycone_release_100um_24h\data\mzml_files\)",
  "/Volumes/bluecub/aglycone_release_100um_24h/data/experiment/mzml/",
  sd$spectraOrigin,
  fixed = TRUE
)
MsExperiment::sampleData(tester) <- sd


take_this <- glycone_feats %>%
  dplyr::filter(group == "Afzelin")

take_fin <- featureDefinitions(tester) %>%
  as_tibble(., rownames = "feature") %>%
  dplyr::filter(
    mzmin > (take_this$mzmin - 0.001) &
    mzmax < (take_this$mzmax + 0.001)
  ) %>%
  dplyr::filter(
    rtmin > (take_this$rtmin - 10) &
    rtmax < (take_this$rtmax + 10)
  ) %>%
  print(n = "all")

take_find

# Now extract the chromatogram
chr <- xcms::chromatogram(
  tester,
  mz = c(take_fin$mzmin - 0.01, take_fin$mzmax + 0.01),
  rt = c(take_fin$rtmin - 10, take_fin$rtmax + 10)
)
plot(chr)
