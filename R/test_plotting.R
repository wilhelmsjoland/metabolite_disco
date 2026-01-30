
lone_feat <- feature.chrs[130,]
full.tib <- tibble::tibble()

for (i in 1:ncol(lone_feat)) {
    
    lone_samp <- lone_feat[[1, i]]
    
    tmp.tib <- tibble::tibble(
        rt = Spectra::rtime(lone_samp),
        int = Spectra::intensity(lone_samp),
        file_idx = MSnbase::fromFile(lone_samp),
        file = colnames(lone_feat)[i],
        mz = stringr::str_flatten(paste0(round(Spectra::mz(lone_samp), 3)), " - ")
    )
    
    full.tib <- bind_rows(full.tib, tmp.tib)
}

p1.data <- full.tib %>%
    dplyr::left_join(
        x = .,
        y = tibble::as_tibble(meta, rownames = "file") %>%
            dplyr::select(file, group),
        by = c("file" = "file")
    )

p1 <- p1.data %>%
    ggplot2::ggplot(.,
                    ggplot2::aes(
                        x = rt,
                        y = int,
                        group = file,
                        color = group
                    )) +
    ggplot2::geom_line(lwd = 1) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
        title = paste0(
            "Feature: ", 130, 
            ", M/z: ", unique(p1.data$mz),
            ", RT: ", 
            round(min(p1.data$rt), 3),
            " - ",
            round(max(p1.data$rt), 3)
        ),
        y = "Intensity",
        x = "Retention time"
    ) +
    ggplot2::scale_color_manual(values = group.colors)

p2.data <- xcms::featureValues(
    lone_feat,
    method = "sum",
    value = "into",
    filled = TRUE,
    missing = "rowmin_half", # 0 # NA
    msLevel = 1
) %>%
    tibble::as_tibble(., rownames = "feature") %>%
    dplyr::mutate(feature2 = stringr::str_remove_all(
        feature, 
        stringr::regex("[A-Za-z]", ignore_case = FALSE))
    ) %>%
    dplyr::mutate(feature2 = as.numeric(feature2)) %>%
    dplyr::relocate(feature2, .after = feature) %>%
    tidyr::pivot_longer(cols = contains(".mzML")) %>%
    dplyr::left_join(
        x = .,
        y = tibble::as_tibble(meta, rownames = "file") %>%
            dplyr::select(file, group),
        by = c("name" = "file")
    ) 

p2 <- p2.data %>%
    ggplot2::ggplot(.,
                    ggplot2::aes(
                        x = group, # name
                        y = value,
                        fill = group
                    )) +
    ggplot2::geom_boxplot(outliers = FALSE) +
    ggplot2::geom_point(
        position =  ggplot2::position_jitter(width = 0.15),
        size = 2,
        pch = 21,
        color = "black"
    ) +
    ggplot2::scale_y_continuous(expand =  ggplot2::expansion(c(0.1, 0.1))) +
    ggplot2::scale_fill_manual(values = group.colors) +
    ggplot2::guides(x =  ggplot2::guide_axis(angle = -45)) +
    ggplot2::theme_bw() +
    ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        axis.title.x = ggplot2::element_blank(),
        legend.title = ggplot2::element_blank()
    ) +
    ggplot2::labs(
        y = paste0("Peak area (", value, ")"),
        caption = paste0(
            "method: ", method,
            ", value: ", value,
            ", filled: ", filled,
            ", missing: ", missing
        )
    )

p3 <- p1 / p2 +
    patchwork::plot_layout(
        axes = "collect",
        guides = "collect",
        heights = c(
            0.55,
            0.45
        )
    )

file.nm <- paste0(res.folder, save_loc,"FT_", feat.idx, ".", device)

if (!is.null(save_loc)) {
    ggplot2::ggsave(
        filename = file.nm,
        plot = p3,
        device = device,
        height = 6,
        width = 6,
        units = "in"
    )
}

data_list <- list(
    "combined" = p3,
    "chromatogram" = p1,
    "boxplot" = p2,
    "p1_data" = p1.data,
    "p2_data" = p2.data
)