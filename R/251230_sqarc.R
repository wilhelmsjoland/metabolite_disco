# TODO 
# Do this for the plotting functions, mabye faster?
df <- purrr::map_dfr(seq_len(ncol(std.chr2)), \(i) {
    x <- std.chr2[, i]
    plot.samp.idx <- which(colnames(std.chr2)[i] == meta$sample)
    
    tibble::tibble(
        sample = colnames(std.chr2)[i],
        rt = x@rtime,
        intensity = x@intensity,
        group = meta$group[plot.samp.idx]
    )
})

df <- df %>%
    dplyr::mutate(group2 = case_when(
        group == "bo_atcc_apiin" ~ "B. ovatus ATCC 8483 + Apiin",
        group == "bo_delta_apiin" ~ "B. ovatus ATCC 8483 ΔOperon + Apiin",
        group == "bt_apiin" ~ "B. thetaiotaomicron VPI-5482 + Apiin",
        group == "bu_atcc_apiin" ~ "B. uniformis ATCC 8492 + Apiin",
        group == "bu_tko_apiin" ~ "B. uniformis gshD,gghC,gshG + Apiin",
        group == "pc_apiin" ~ "P.copri iAK263 + Apiin",
    )) %>%
    dplyr::mutate(group2 = forcats::fct_relevel(
        group2,
        c(
            "B. ovatus ATCC 8483 + Apiin",
            "B. ovatus ATCC 8483 ΔOperon + Apiin",
            "B. thetaiotaomicron VPI-5482 + Apiin",
            "B. uniformis ATCC 8492 + Apiin",
            "B. uniformis gshD,gghC,gshG + Apiin",
            "P.copri iAK263 + Apiin"
        )
    )) %>%
    dplyr::arrange(group2) %>%
    dplyr::mutate(sample = factor(sample, levels = unique(sample)))


apiin.standards.p <- df %>%
    ggplot(.,
           aes(
               x = rt,
               y = intensity,
               color = group2
           )) +
    geom_line(linewidth = 1) +
    facet_wrap(~ sample, ncol = 4) +
    theme_bw() +
    theme(
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        strip.text = element_blank()
    ) +
    labs(
        y = "Intensity",
        x = "Retention time (s)",
        title = "Chromatograms of Apiin"
    )

ggsave(
    filename = file.path(stds.output.path, "apiin_standards.svg"),
    plot = apiin.standards.p,
    device = "svg",
    height = 5,
    width = 10,
    units = "in"
)
