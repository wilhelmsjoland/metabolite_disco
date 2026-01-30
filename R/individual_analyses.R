apiin.apigenin <- matched.diffs %>%
    dplyr::filter(
        dplyr::if_any(
            dplyr::all_of(c("mz1", "mz2")),
            ~ between(.x, 268, 271)
        )
    ) %>%
    print(n = 100)

matched.diffs %>%
    dplyr::filter(
        dplyr::if_any(
            dplyr::all_of(c("mz1", "mz2")),
            ~ stringr::str_detect(., "1757")
        )
    )

full.data %>%
    dplyr::filter(
        dplyr::if_any(
            dplyr::all_of(c("mzmed")),
            ~ between(.x, getTheoryMz("C15H10O5") - 0.005, getTheoryMz("C15H10O5") + 0.005)
        )
    ) %>%
    print(n = 100)


matched.diffs %>%
    dplyr::filter(stringr::str_detect(name, "apiofurano")) %>%
    dplyr::filter(stringr::str_detect(name, "1 x")) %>%
    dplyr::arrange(desc(mz1)) %>%
    print(n = 1000)

full.data %>%
    dplyr::filter(between(mzmed, 268, 271)) %>%
    dplyr::filter(between(rtmed, (4.7*60), (4.9*60))) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ft_mz_rt = paste0(
        feature, 
        "; m/z: ", round(mzmed, 3), 
        "; rtime: ", round(rtmed, 3)
    )) %>%
    tidyr::pivot_longer(cols = all_of(rownames(meta))) %>%
    dplyr::left_join(
        x = .,
        y = as_tibble(meta, rownames = "name"),
        by = "name"
    ) %>%
    dplyr::mutate(group = forcats::fct_relevel(
        group,
        c(
            "bu_wt_control",
            "bu_wt_apiin",
            "bu_mutant_control",
            "bu_mutant_apiin"
        )
    )) %>%
    ggplot(.,
           aes(
               x = group,
               y = value,
               fill = group
           )) +
    geom_boxplot(outliers = FALSE) +
    geom_point() +
    facet_wrap(~ft_mz_rt, scales = "free_y", strip.position = "left") +
    guides(x = guide_axis(angle = -45)) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()
    )



# TODO
apiin.apigenin

# FT02117 
# FT07320

# Apiin
getMass(getMolecule("C26H28O14")) - getMass(getMolecule("H"))
# Apigenin
getMass(getMolecule("C15H10O5")) - getMass(getMolecule("H"))

full.data %>%
    dplyr::filter(feature %in% c("FT02117", "FT07320")) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(ft_mz_rt = paste0(
        feature, 
        "; m/z: ", round(mzmed, 3), 
        "; rtime: ", round(rtmed, 3)
    )) %>%
    tidyr::pivot_longer(cols = all_of(rownames(meta))) %>%
    dplyr::left_join(
        x = .,
        y = as_tibble(meta, rownames = "name"),
        by = "name"
    ) %>%
    dplyr::mutate(group = forcats::fct_relevel(
        group,
        c(
            "bu_wt_control",
            "bu_wt_apiin",
            "bu_mutant_control",
            "bu_mutant_apiin"
        )
    )) %>%
    ggplot(.,
           aes(
               x = group,
               y = value,
               fill = group
           )) +
    geom_boxplot(outliers = FALSE) +
    geom_point() +
    facet_wrap(~ft_mz_rt, scales = "free_y", strip.position = "left") +
    guides(x = guide_axis(angle = -45)) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank()
    )

getTheoryMz("C26H28O14")
getTheoryMz("C15H10O5")

full.data %>%
    dplyr::filter(
        dplyr::if_any(
            dplyr::all_of(c("mzmed")),
            ~ between(.x, 269.0455 - 0.01, 269.0455 + 0.01) | between(.x, 563.1406 - 0.01, 563.1406 + 0.01)
        )
    ) %>%
    dplyr::arrange(desc(rtmed)) %>%
    print(n = 100)

full.data %>%
    dplyr::filter(
        dplyr::if_any(
            dplyr::all_of(c("mzmed")),
            ~ between(.x, 269.0455 - 0.1, 269.0455 + 0.1) | between(.x, 563.1406 - 0.1, 563.1406 + 0.1)
        )
    ) %>%
    dplyr::arrange(desc(rtmed))
    print(n = 100)

matched.diffs %>%
    dplyr::filter(
        dplyr::if_any(
            dplyr::all_of(c("mz1", "mz2")),
            ~ between(.x, 269.0455 - 0.1, 269.0455 + 0.1) | between(.x, 563.1406 - 0.1, 563.1406 + 0.1)
        )
    ) %>%
    dplyr::arrange(desc(mz1)) %>%
    print(n = 1000)

intersecting.feats.p <- intersect.data %>%
    ggplot(.,
           aes(
               x = group,
               y = value
           )) +
    geom_boxplot(
        aes(fill = group),
        outliers = FALSE
    ) +
    geom_point() +
    facet_wrap(
        facets = ~ feature, 
        scales = "free_y"
    ) +
    scale_y_continuous(expand = expansion(c(0.1, 0.15))) +
    guides(x = guide_axis(angle = -45)) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        legend.title = element_blank()
    ) +
    labs(
        y = "Log2 median-scaled intensity"
    ) +
    ggpubr::geom_bracket(
        data = limma.p.res %>%
            dplyr::filter(!adj.P.Val.signif %in% c("ns")),
        aes(
            xmin = group1,
            xmax = group2,
            label = adj.P.Val.signif,
            y.position = y.pos
        ),
        step.increase = 0.015
    )
ggsave(
    filename = paste0(res.folder, "/graphs/bar/", "all_intersecting.pdf"),
    plot = intersecting.feats.p,
    device = "pdf",
    height = 10,
    width = 10,
    units = "in"
)
