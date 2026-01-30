# apiin
glycone <- getTheoryMz("C26H28O14", adduct = "[M-H]-")
# aglycone
aglycone <- getTheoryMz("C15H10O5", adduct = "[M-H]-")

# Why is this so insanely high?
range.tol <- ppm_to_num(2000)

# range to use
glycone.range <- glycone + c(-range.tol, +range.tol)
aglycone.range <- aglycone + c(-range.tol, +range.tol)


gly.agly <- full_raw_filled %>%
    dplyr::filter(
        dplyr::between(mzmed, min(glycone.range), max(glycone.range)) |
        dplyr::between(mzmed, min(aglycone.range), max(aglycone.range))
    )

# gly.agly[, "ppm_error"] <- num_to_ppm(abs(gly.agly[,"mzmed"] - aglycone))
# pot.glycosides <- as.numeric(gsub("[A-Za-z]", "", unique(gly.agly$feature)))

pot.glycosides <- unique(gly.agly$feature)

# Fix so this function also returns the peak area and the peak chromatogram
glycoside.data <- list()
for (i in pot.glycosides) {
    tmp.glycoside <- plotFeatChrInt(
        feature_chrom = feature.chrs,
        feature = i,
        method = "sum",
        value = "into",
        filled = TRUE,
        missing = "rowmin_half",
        msLevel = 1,
        save_loc = "/graphs/glycosylation/",
        device = "pdf",
        feat_pairs = FALSE
    )
    
    # list.name <- paste0("FT_", i)
    list.name <- paste0(i)
    glycoside.data[[list.name]] <- tmp.glycoside
}

ft.for.hmp <- tibble::tibble()
for (i in names(glycoside.data)) {
    tmp.ft.for.hmp <- glycoside.data[[i]]$p2_data
    ft.for.hmp <- bind_rows(ft.for.hmp, tmp.ft.for.hmp)
}

# Could use what i already have scaled as well
ft.for.hmp %>%
    dplyr::group_by(feature) %>%
    dplyr::mutate(value = scale(value, center = TRUE, scale = TRUE)[, 1]) %>%
    dplyr::group_by(feature, group) %>%
    dplyr::summarize(value = mean(value)) %>%
    ggplot2::ggplot(.,
           aes(
               x = group,
               y = feature,
               fill = value
           )) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(
        low = "cornflowerblue",
        mid = "white",
        high = "firebrick"
    ) +
    ggplot2::scale_y_discrete(expand = expansion(c(0, 0))) +
    ggplot2::scale_x_discrete(expand = expansion(c(0, 0))) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title = element_blank())