
read_csv(biotransf.file) %>%
    gt::gt() %>%
    gt::opt_stylize(style = 1, color = "gray", add_row_striping = TRUE) %>%
    gt::cols_align(align = "center") %>%
    gtsave(filename = paste0(res.folder, "/tables/orig_biotransf.png"))


bio.transf %>%
    gt::gt() %>%
    gt::opt_stylize(style = 1, color = "gray", add_row_striping = TRUE) %>%
    gt::cols_align(align = "center") %>%
    gtsave(filename = paste0(res.folder, "/tables/2nd_biotransf.html"))

meta %>%
    gt::gt() %>%
    gt::opt_stylize(style = 1, color = "gray", add_row_striping = TRUE) %>%
    gt::cols_align(align = "center") %>%
    gtsave(filename = paste0(res.folder, "/tables/2nd_meta.png"))

read_csv(paste0(data.path, "/", meta.file)) %>%
    gt::gt() %>%
    gt::opt_stylize(style = 1, color = "gray", add_row_striping = TRUE) %>%
    gt::cols_align(align = "center") %>%
    gtsave(filename = paste0(res.folder, "/tables/orig_meta.png"))

upset.tib.t <- full.limma %>%
    dplyr::select(feature, contrast, adj.P.Val) %>%
    tidyr::pivot_wider(
        names_from = "contrast",
        values_from = "adj.P.Val"
    ) %>%
    dplyr::mutate(
        across(
            .cols = 2:ncol(.),
            .fns = ~ if_else(
                . < 0.05,
                TRUE,
                FALSE
            )
        )
    )

upset.p2 <- upset.tib %>%
    ComplexUpset::upset(
        intersect = comparisons,
        name = "Features with p adjusted < 0.05",
        width_ratio = 0.15,
        base_annotations = list(
            "Intersecting features" = intersection_size()
        ),
        min_degree = 3,
        max_degree = 3
    )

upset.p2

ggsave(
    filename = paste0(res.folder, "/graphs/upset/upset_pres.pdf"),
    plot = upset.p2,
    device = "pdf",
    height = 5,
    width = 13,
    units = "in"
)
