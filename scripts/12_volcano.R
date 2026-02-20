message("Producing volcano plots...")
volc_plot_list <- list()
for (i in names(limma_res)) {

  tmp <- limma_res[[i]]

  tmp_tib <- tmp %>%
    dplyr::mutate(
      label.p = dplyr::if_else(
        adj.P.Val < opt$qvalue & abs(logFC) > quantile(abs(logFC), 0.99),
        mzmed,
        NA
      ),
      direction = dplyr::case_when(
        logFC >= 0 & adj.P.Val < opt$qvalue ~ "Up",
        logFC < 0 & adj.P.Val < opt$qvalue ~ "Down",
        adj.P.Val >= opt$qvalue ~ "ns",
        TRUE ~ as.character("check")
      )
    ) %>%
    dplyr::mutate(direction = forcats::fct_relevel(
      direction,
      c(
        "Up",
        "ns",
        "Down"
      )
    )) %>%
    tidyr::drop_na(logFC)

  tmp_p <- tmp_tib %>%
    ggplot2::ggplot(
      ggplot2::aes(
        x = logFC,
        y = -log10(adj.P.Val),
        color = direction
      )
    ) +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = c(
      "Up" = "firebrick",
      "Down" = "cornflowerblue",
      "ns" = "grey"
    )) +
    ggrepel::geom_label_repel(
      data = tidyr::drop_na(tmp_tib, label.p) %>%
        dplyr::arrange(dplyr::desc(abs(logFC))) %>%
        dplyr::slice(1:50),
      ggplot2::aes(label = round(label.p, 2)),
      size = 3,
      max.overlaps = 100,
      box.padding = 0.5,
      color = "black"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      legend.title = ggplot2::element_blank(),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    ) +
    ggplot2::labs(
      title = i,
      subtitle = paste0(
        "Rounded mass-to-charge ratios with adj.p < ",
        opt$qvalue,
        " are labelled"
      ),
      x = "Log2 fold change",
      y = "-Log10 adjusted p-value"
    )

  volc_plot_list[[i]] <- tmp_p

  ggplot2::ggsave(
    filename = paste0(opt$output, "/graphs/volcano/", i, ".pdf"),
    plot = tmp_p,
    device = "pdf",
    height = 10,
    width = 10,
    units = "in"
  )
}