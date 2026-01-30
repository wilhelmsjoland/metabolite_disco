importFiles <- function(
        data.path = NULL,
        meta.file = NULL
        ) {
    mzml.files <- list.files(
        data.path, 
        # pattern = ".test",
        pattern = "\\.mzML$|\\.mzml$", 
        full.names = TRUE, 
        recursive = TRUE, 
        ignore.case = TRUE
    )
    
    if (length(mzml.files) == 0) {
        stop("No .mzML files found in path")
    }
    
    mzml.tib <- data.frame(
        sample = basename(mzml.files),
        path = mzml.files
    )
    
    meta <<- readr::read_csv(
        file = meta.file,
        show_col_types = FALSE
        )
    # %>%
    #     bind_rows(tibble(sample = "x.mzML", group = "y"))
    
    meta.matched <- dplyr::left_join(
        x = meta,
        y = mzml.tib,
        na_matches = "never",
        relationship = "one-to-one",
        by = join_by(sample)
    )
    
    samp.lacking.file <- meta.matched$sample[!meta.matched$sample %in% basename(mzml.files)]
    if(length(samp.lacking.file) > 0) {
        stop(paste0("Samples in metadata lacking files:", "\n", samp.lacking.file))
    }
    non.used.files <- mzml.files[!mzml.files %in% meta.matched$path]
    if (length(non.used.files) > 0) {
        warning("\nFiles in path not in metadata:", paste0("\n", non.used.files))
    }
    
    meta.matched.filt <- meta.matched %>%
        dplyr::filter(!sample %in% samp.lacking.file) %>%
        tibble::column_to_rownames(var = "sample")
    
    return(meta.matched.filt)
}

multChem <- function(chem_formula = NULL, multiply_by = 1) {
    split.form <- stringr::str_extract_all(chem_formula, "[A-Z][a-z]?[0-9]*")[[1]]
    split.form.append <- stringr::str_replace_all(split.form, "([A-Za-z])$", "\\11")
    multiply.by <- multiply_by
    for (i in 1:length(split.form.append)) {
        chem.letter <- stringr::str_extract_all(split.form.append[i], "[A-Za-z]+")[[1]]
        chem.number <- as.integer(stringr::str_extract_all(split.form.append[i], "[0-9]+$")[[1]])
        split.form.append[i] <- paste0(chem.letter, chem.number * multiply.by)
    }
    split.form.append <- stringr::str_flatten(split.form.append)
    clean.form <- Rdisop::getFormula(Rdisop::getMolecule(split.form.append))
    return(clean.form)
}

findIntersectFeat <- function(data = NULL, set = NULL, full_set = NULL) {
    filt.data <- data %>%
        dplyr::filter(
            dplyr::if_all(
                .cols = tidyselect::all_of(full_set[full_set %in% set]),
                .fns = ~ .x == TRUE
            ),
            dplyr::if_all(
                .cols = tidyselect::all_of(full_set[!full_set %in% set]),
                .fns = ~ .x == FALSE
            )
        )
    return(filt.data)
}

find_y_position <- function(
        test.df, 
        df, 
        formula, 
        fun.data, 
        grouping = NULL
        ) {
    
    if (!fun.data %in% c(
        "max",
        "mean",
        "min",
        "mean_sem",
        "max_sem",
        "lower_bound",
        "upper_bound"
    )) {
        stop("Incorrect specification of fun.data")
    }
    
    num.col <- stringr::str_trim(stringr::str_split_i(formula, "~", 1))
    group.col <- stringr::str_trim(stringr::str_split_i(formula, "~", 2))
    
    sum.df <- df %>%
        dplyr::group_by(.data[[group.col]]) %>%
        dplyr::summarize(
            max = base::max(.data[[num.col]], na.rm = TRUE),
            mean = base::mean(.data[[num.col]], na.rm = TRUE),
            min = base::min(.data[[num.col]], na.rm = TRUE),
            mean_sem = mean + stats::sd(.data[[num.col]])/base::sqrt(base::length(.data[[num.col]])),
            max_sem = max + stats::sd(.data[[num.col]])/base::sqrt(base::length(.data[[num.col]])),
            Q1 = stats::quantile(.data[[num.col]], 0.25, na.rm = TRUE),
            Q3 = stats::quantile(.data[[num.col]], 0.75, na.rm = TRUE),
            IQR.val = stats::IQR(.data[[num.col]], na.rm = TRUE),
            lower_bound = Q1 - 1.5 * IQR.val,
            upper_bound = Q3 + 1.5 * IQR.val
        )
    
    sum.df <- sum.df %>%
        dplyr::select(.data[[group.col]], tidyselect::all_of(fun.data)) %>%
        dplyr::rename("y.pos" = .data[[fun.data]])
    
    group1.names <- setNames(group.col, "group1")
    group2.names <- setNames(group.col, "group2")
    
    final.df <- test.df %>%
        dplyr::left_join(
            x = .,
            y = sum.df,
            by = group1.names
        ) %>%
        dplyr::left_join(
            x = .,
            y = sum.df,
            by = group2.names
        ) %>%
        {
            if (!is.null(grouping)) {
                dplyr::group_by(., dplyr::across(dplyr::all_of(grouping)))
            } else {
                .
            }
        } %>%
        # dplyr::group_by(comparisons) %>%
        dplyr::mutate(y.pos = base::max(y.pos.x, y.pos.y, na.rm = TRUE)) %>%
        dplyr::select(-c("y.pos.x", "y.pos.y"))
    
    return(final.df)
}

predictBiotransf <- function(
        data = NULL,
        biotransf.data = NULL,
        tolerance_ppm = NULL,
        tolerance = NULL
) {
    
    if (is.null(tolerance) & is.null(tolerance_ppm)) {
        stop("Both tolerance and tolerance_ppm are NULL")
    } else if (!is.null(tolerance) & !is.null(tolerance_ppm)) {
        stop("Either tolerance or tolerance_ppm need to be set to NULL")
    }
    
    if (is.null(tolerance)) {
        tol_used = ppm_to_num(tolerance_ppm)
    }
    if (is.null(tolerance_ppm)) {
        tol_used = tolerance
    }
    
    ## 1. Prepare peaks table (all peaks, sorted by m/z)
    peaks <- data %>%
        # needs to be sorted for findInterval() indexing
        dplyr::arrange(mzmed) %>%
        dplyr::mutate(
            peak_id = dplyr::row_number(),    # simple integer ID
            feature,
            mzmed
        )
    n_peaks <- nrow(peaks)

    # convenience vectors so we don't keep indexing peaks$...
    mz_vec <- peaks$mzmed
    rt_vec <- peaks$rtmed
    id_vec <- peaks$peak_id
    feat_vec <- peaks$feature
    n_trans <- nrow(biotransf.data)

    ## 3. For each biotransformation, find all matching peak pairs
    all_matches <- vector("list", n_trans)

    for (k in 1:n_trans) {

        # current transformation
        delta <- biotransf.data$delta_mass[k]
        this_name <- biotransf.data$name[k]
        this_formula <- biotransf.data$chem_formula[k]

        # for each peak i, valid partners j must have:
        # mz_vec[j] in [mz_vec[i] + delta - tol, mz_vec[i] + delta + tol]
        target_lower <- mz_vec + delta - tol_used
        target_upper <- mz_vec + delta + tol_used

        # find, for each i, the index range [start_i, end_i] in the sorted mz_vec
        # that lies within [target_lower[i], target_upper[i]]
        idx_start <- findInterval(target_lower, mz_vec) + 1L
        idx_end <- findInterval(target_upper, mz_vec)

        # collect matches for this transformation here
        res_list <- vector("list", n_peaks)

        for (i in 1:n_peaks) {
            start_idx <- idx_start[i]
            end_idx <- idx_end[i]

            # no overlap for this i -> skip
            if (start_idx > end_idx) next

            # candidate partner indices
            partner_idx <- seq.int(start_idx, end_idx)

            # avoid self-pairs and symmetric duplicates (i,j) vs (j,i)
            partner_idx <- partner_idx[partner_idx > i]
            if (length(partner_idx) == 0) next

            # build rows for all partners of peak i
            res_list[[i]] <- tibble::tibble(
                name = this_name,
                chem_change = this_formula,
                delta_mass = delta,
                feat1 = feat_vec[i],
                feat2 = feat_vec[partner_idx],
                mz1 = mz_vec[i],
                mz2 = mz_vec[partner_idx],
                rt1 = rt_vec[i],
                rt2 = rt_vec[partner_idx],
                obs_delta_mass = mz_vec[partner_idx] - mz_vec[i],
                peak1_id = id_vec[i],
                peak2_id = id_vec[partner_idx]
            )
        }

        # bind all i-level results for this transformation
        all_matches[[k]] <- dplyr::bind_rows(res_list)
    }

    # 4. Final table of all matched pairs
    matched.diffs <- dplyr::bind_rows(all_matches)

    return(matched.diffs)
}

importBiotransfMeta <- function(
        file = NULL
        ) {
    biotransf.meta <- readr::read_csv(
        file = file,
        comment = "#",
        show_col_types = FALSE
        ) %>%
        tidyr::uncount(., weights = allowed_n, .id = "multiplier", .remove = FALSE) %>%
        dplyr::mutate(name = paste0(multiplier, " x ", name)) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(chem_formula = multChem(chem_formula, multiplier)) %>%
        dplyr::mutate(delta_mass = getTheoryDeltaMass(chem_formula)) %>%
        dplyr::mutate(chem_formula = paste0("± ", chem_formula)) %>%
        dplyr::ungroup()
    
    return(biotransf.meta)
}

num_to_ppm <- function(number = NULL) {
    ppm.val <- number * 1e6
    return(ppm.val)
}

ppm_to_num <- function(number = NULL) {
    ppm.val <- number / 1e6
    return(ppm.val)
}

getTheoryMz <- function(
        chem_form = NULL,
        adduct = "[M-H]-"
        ) {
    chem.mass <- MetaboCoreUtils::calculateMass(chem_form)
    chem.theory.mz <- MetaboCoreUtils::mass2mz(chem.mass, adduct)
    chem.theory.mz <- chem.theory.mz[1, 1]
    
    return(chem.theory.mz)
}

getTheoryDeltaMass <- function(
        chem_form = NULL
        ) {
    chem.mass <- MetaboCoreUtils::calculateMass(chem_form)
    chem.theory.mass <- unname(chem.mass)
    
    return(chem.theory.mass)
}

getShortMzRange <- function(
    chem.theory.mz = NULL,
    mz.window = 0.01
    ) {
    chem.mz.range <- c(chem.theory.mz - mz.window, chem.theory.mz + mz.window)
    return(chem.mz.range)
}

getRtMzRange <- function(
        chromatogram = NULL,
        rt_window = 0.02 # percentage increase of window
) {
    ref.tib <- tibble::tibble()
    for (samp.n in 1:length(chromatogram)) {
        ref.peak.intensity <- max(intensity(chromatogram[1, samp.n]), na.rm = TRUE)
        max.intensity.idx <- which(intensity(chromatogram[1, samp.n]) == ref.peak.intensity)
        ref.peak.rt <- Spectra::rtime(chromatogram[1, samp.n])[[max.intensity.idx]]
        mzmin <- min(Spectra::mz(chromatogram[1, samp.n]), na.rm = TRUE)
        mzmax <- max(Spectra::mz(chromatogram[1, samp.n]), na.rm = TRUE)
        rt.window.min <- c((1 - rt_window) * ref.peak.rt)
        rt.window.max <- c((1 + rt_window) * ref.peak.rt)
        file <- colnames(chromatogram)[samp.n]
        
        tmp.tib <- tibble::tibble(
            file,
            ref.peak.intensity, 
            ref.peak.rt,
            rt.window.min,
            rt.window.max,
            max.intensity.idx, 
            mzmin,
            mzmax,
            samp.n
            )
        
        ref.tib <- dplyr::bind_rows(ref.tib, tmp.tib)
    }
    
    rt.range <- c(min(ref.tib$rt.window.min), max(ref.tib$rt.window.max))
    mz.range <- c(min(ref.tib$mzmin), max(ref.tib$mzmax))
    
    ranges <- list(
        "mz.range" = mz.range, 
        "rt.range" = rt.range,
        "data" = ref.tib
        )
    
    return(ranges)
}

# TODO
# Rewrite this
# Could loop over all samples to identify which
# intensities for a rt are most consistent later
# findMzError <- function(
#         object = NULL,
#         ref_peaks_obj = NULL
#         ) {
#     
#     rt.range <- c(min(ref_peaks_obj$rt.window.min), max(ref_peaks_obj$rt.window.max))
#     mz.range <- c(min(ref_peaks_obj$mzmin, na.rm = TRUE), max(ref_peaks_obj$mzmax, na.rm = TRUE))
#     
#     srn <- object |> 
#         Spectra::filterRt(rt = rt.range) |>
#         Spectra::filterMzRange(mz.range)
#     
#     print(plot(srn))
#     
#     min.filter <- readline("Write min: ")
#     max.filter <- readline("Write max: ")
#     
#     full.tmp.tib <- tibble()
#     for (i in 1:length(srn)) {
#         srn_1 <- srn[i] |> 
#             Spectra::filterRt(rt = rt.range) |>
#             Spectra::filterMzRange(c(min.filter, max.filter)) |>
#             spectra()
#         
#         try3 <- srn_1 |>
#             Spectra::mz() |>
#             unlist() |>
#             diff() |>
#             abs()
#         
#         # Differences in m/z values expressed as ppm
#         mz.error <- try3 * 1e6 / mean(unlist(mz(srn_1)))
#         
#         max.mz.error <- max(mz.error, na.rm = TRUE)
#         min.mz.error <- max(mz.error, na.rm = TRUE)
#         
#         spectra.data <- Spectra::spectraData(srn_1)
#         spectra.origin <- unique(spectra.data$dataOrigin)
#         spectra.basename <- basename(spectra.origin)
#         
#         tmp.tib <- tibble(
#             max.mz.error,
#             min.mz.error,
#             i,
#             spectra.basename,
#             min.filter,
#             max.filter
#         )
#         
#         full.tmp.tib <- bind_rows(full.tmp.tib, tmp.tib)
#     }
#     
#     return(full.tmp.tib)
# }

# TODO 
# Add colors to the plotting as well
inspectPeak <- function(
        chromatogram = NULL,
        peak_data = NULL,
        peak.idx = TRUE, # if true pick a random peak
        sample.no = FALSE,
        save.graph = TRUE
        ) {
    
    if (isTRUE(peak.idx)) {
        pk.chk <- sample(nrow(peak_data), 1)
    } else {
        pk.chk <- peak.idx
    }
    
    if (isFALSE(sample.no)) {
        test_peak <- chromatogram %>%
            Spectra::filterMzRange(c(peak_data[pk.chk, ]$mzmin, peak_data[pk.chk, ]$mzmax)) %>%
            Spectra::filterRt(c(peak_data[pk.chk, ]$rtmin, peak_data[pk.chk, ]$rtmax)) |>
            xcms::chromatogram(aggregationFun = "sum")
    } else {
        test_peak <- chromatogram[sample.no] %>%
            Spectra::filterMzRange(c(peak_data[pk.chk, ]$mzmin, peak_data[pk.chk, ]$mzmax)) %>%
            Spectra::filterRt(c(peak_data[pk.chk, ]$rtmin, peak_data[pk.chk, ]$rtmax)) |>
            xcms::chromatogram(aggregationFun = "sum")
    }
    if (save.graph == TRUE) {
        pdf(paste0(res.folder, "/graphs/quality_control/", pk.chk, ".pdf"))
        plot(
            x = test_peak,
            main = paste0(
                "Peak: ", peak_data[pk.chk, ]$peak, 
                "\n m/z: ", round(peak_data[pk.chk, ]$mzmin, 3), 
                " - ", round(peak_data[pk.chk, ]$mzmax, 3),
                "\n RT: ", round(peak_data[pk.chk, ]$rtmin, 3), 
                " - ", round(peak_data[pk.chk, ]$rtmax, 3)
            ),
            sub = paste0(
                "Samples: ",
                stringr::str_flatten(colnames(inspect_peaks), ", ")
            ),
            peakBg = NA,
            lwd = 3
        )
        dev.off()
    }

    return(test_peak)
}

inspectPeakInt <- function(
        chr_data = NULL,
        value = "into", # into, intb, maxo
        save.graph = TRUE
        ) {
    # Extract a list of per-sample peak intensities (in log2 scale)
    # TODO Do the samples match the colors?
    xchr_peaks <- tibble::as_tibble(xcms::chromPeaks(chr_data), rownames = "peak")
    xchr_meta <- tibble::as_tibble(MsExperiment::sampleData(chr_data), rownames = "sample_id") %>%
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
        )
    
    # lower number of detected peaks = the smaller width of the boxes
    xchr_data_p <- xchr_data_comb %>%
        dplyr::mutate(comb = paste0(sample, "_", sample_id)) %>%
        dplyr::mutate(comb = forcats::fct_reorder(
            .f = comb,
            .x = sample
        )) %>%
        ggplot2::ggplot(.,
                        ggplot2::aes(
                   x = comb,
                   y = {{ value }},
                   fill = group
               )) +
        ggplot2::geom_boxplot(varwidth = TRUE) +
        ggplot2::scale_fill_manual(values = group.colors) +
        ggplot2::guides(x = guide_axis(angle = -45)) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            legend.title = ggplot2::element_blank()
        ) +
        ggplot2::labs(
            y = "Log2 intensity",
            subtitle = paste0(
                deparse(substitute(value)), 
                ", box width: number of detected peaks"
                )
        )
    
    file.nm <- paste0(
        res.folder, 
        "/graphs/per_sample_peaks/", 
        deparse(substitute(xchr)), 
        "_detected_peaks.pdf"
        )
    
    if (save.graph == TRUE) {
        ggplot2::ggsave(
            filename = file.nm,
            plot = xchr_data_p,
            device = "pdf",
            height = 6,
            width = 6,
            units = "in"
        )
    }
    
    return(xchr_data_p)
}

plotChromatogramInts <- function(
        chromatogram = NULL,
        chrom_object = NULL,
        save_loc = NULL,
        amount = NULL,
        peaks_or_feats = c("features", "peaks")
) {
    seq.plots <- seq.int(from = 0, to = nrow(chromatogram), by = 20)
    # TODO This doesn't work when divisible by 20
    # chrs.na.seq.slices <- c(seq.plots, nrow(chromatogram)) # add until last value
    chrs.na.seq.slices <- pmin(seq.plots + 20 - 1, seq.plots)
    
    if (is.null(amount)) {
        length.chrs.na.seq.slices <- length(chrs.na.seq.slices)-1
    } else {
        length.chrs.na.seq.slices <- amount
    }
    
    for (i in 1:length.chrs.na.seq.slices) {
        # choose slices e.g. 1 : 20, 21:40, 41:60
        # remove first value of end to not plot same value twice
        first.slice <- chrs.na.seq.slices[i]+1 # because starts from 0
        second.slice <- chrs.na.seq.slices[i+1] # because access next part of slice
        
        sliced.chr <- chromatogram[first.slice:second.slice,]
        
        if (peaks_or_feats == "features") {
            slice.chr.peaks <- as.data.frame(xcms::featureDefinitions(sliced.chr))
        } else if (peaks_or_feats == "peaks") {
            slice.chr.peaks <- as.data.frame(xcms::chromPeaks(sliced.chr))
        } else {
            stop("'peaks_or_feats' needs to be either 'features' or 'peaks")
        }
        
        min.mz <- min(slice.chr.peaks$mzmin, na.rm = TRUE)
        max.mz <- max(slice.chr.peaks$mzmax, na.rm = TRUE)
        
        file.name <- paste0(round(min.mz, 3), "-", round(max.mz, 3))
        
        pdf(
            file = paste0(res.folder, save_loc, file.name, ".pdf"),
            width = 12,
            height = 10,
        )

        plot(
            sliced.chr,
            # TODO triple check if this is correct
            # Probably do it the same way i do it for the plotFeatChrInt()
            col = group.colors[MsExperiment::sampleData(chrom_object)$group],
            peakType = "polygon",
            peakBg = NA,
            lwd = 2
        )
        invisible(dev.off())
    }

}

plotFeatChrInt <- function(
        feature_chrom = NULL,
        feature = NULL,
        method = NULL,
        value = NULL,
        filled = FALSE,
        missing = NULL,
        msLevel = 1,
        save_loc = NULL,
        device = "pdf",
        feat_pairs = FALSE
        ) {
    
    feat.tib <<- tibble::as_tibble(xcms::featureDefinitions(feature_chrom), rownames = "feature")
    feat.idx <- which(feat.tib$feature == feature)
    lone_feat <- feature_chrom[feat.idx,]
    lone_feat_def <- tibble::as_tibble(xcms::featureDefinitions(lone_feat), rownames = "feature")
    
    # Only take the ones with peaks otherwise the coloring in the base plot is wrong
    keep_peaks <- xcms::hasChromPeaks(lone_feat)
    lone_feat_peaks <- xcms::hasChromPeaks(lone_feat)[, keep_peaks]
    base.group.colors <- meta$group[match(names(lone_feat_peaks), rownames(meta))]
    
    # full.tib <- tibble::tibble()
    # for (i in 1:ncol(lone_feat)) {
    #     
    #     lone_samp <- lone_feat[[1, i]]
    #     
    #     tmp.tib <- tibble::tibble(
    #         rt = Spectra::rtime(lone_samp),
    #         int = Spectra::intensity(lone_samp),
    #         file_idx = MSnbase::fromFile(lone_samp),
    #         file = colnames(lone_feat)[i],
    #         mz = stringr::str_flatten(paste0(round(Spectra::mz(lone_samp), 3)), " - ")
    #     )
    #     
    #     full.tib <- bind_rows(full.tib, tmp.tib)
    # }
    # 
    # p1.data <- full.tib %>%
    #     dplyr::left_join(
    #         x = .,
    #         y = tibble::as_tibble(meta, rownames = "file") %>%
    #             dplyr::select(file, group),
    #         by = c("file" = "file")
    #     )
    # 
    # p1 <- p1.data %>%
    #     ggplot2::ggplot(.,
    #            ggplot2::aes(
    #                x = rt,
    #                y = int,
    #                group = file,
    #                color = group
    #            )) +
    #     ggplot2::geom_line(lwd = 1) +
    #     ggplot2::theme_bw() +
    #     ggplot2::theme(
    #         legend.title = ggplot2::element_blank()
    #     ) +
    #     ggplot2::labs(
    #         title = paste0(
    #             "Feature: ", feature, 
    #             ", M/z: ", unique(p1.data$mz),
    #             ", RT: ", 
    #             round(min(p1.data$rt), 3),
    #             " - ",
    #             round(max(p1.data$rt), 3)
    #         ),
    #         y = "Intensity",
    #         x = "Retention time"
    #     ) +
    #     ggplot2::scale_color_manual(values = group.colors)

    if (isTRUE(feat_pairs)) {
        p1 <- function() {
            par(
                mar = c(bottom = 2.3, left = 0, top = 0, right = 0),
                mgp = c(1.5, 0.5, 0),
                tck = -0.02
            )
            plot(
                lone_feat, 
                peakType = "polygon",
                peakCol = group.colors[base.group.colors],
                peakBg = NA,
                lwd = 3,
                main = "",
                ylab = "",
                xlab = "",
                xaxt = "n",
                yaxt = "n"
            )
            axis(1, las = 0, cex.axis = 0.8)
            axis(2, las = 2, cex.axis = 0.8)
            mtext("Intensity", side = 2, line = 2.8)
            mtext("Retention time", side = 1, line = 1.5)
        }
    } else {
        p1 <- function() {
            par(mar = c(5.1, 4.1, 4.1, 8.6))
            plot(
                lone_feat, 
                peakType = "polygon",
                peakCol = group.colors[base.group.colors],
                peakBg = NA,
                lwd = 3,
                main = paste0(
                    lone_feat_def$feature, ", ",
                    "m/z: ",
                    round(lone_feat_def$mzmin, 3),
                    " - ",
                    round(lone_feat_def$mzmax, 3),
                    ", RT: ",
                    round(lone_feat_def$rtmin, 2),
                    " - ",
                    round(lone_feat_def$rtmax, 2)
                ),
                ylab = "",
                xlab = "",
                xaxt = "n",
                yaxt = "n"
            )
            axis(1, las = 0, cex.axis = 0.8)
            axis(2, las = 2, cex.axis = 0.8)
            mtext("Intensity", side = 2, line = 3.2)
            mtext("Retention time (s)", side = 1, line = 2.5)
        }
    }
    

    
    p2.data <- xcms::featureValues(
        lone_feat,
        method = method,
        value = value,
        filled = filled,
        missing = missing, # 0 # NA
        msLevel = msLevel
    ) %>%
        tibble::as_tibble(., rownames = "feature") %>%
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
            position = ggplot2::position_jitter(width = 0.15),
            size = 2,
            pch = 21,
            color = "black"
            ) +
        ggplot2::scale_y_continuous(expand =  ggplot2::expansion(c(0.1, 0.1))) +
        ggplot2::scale_fill_manual(values = group.colors) +
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
    
    # p3 <- p1 / p2 +
    #     patchwork::plot_layout(
    #         axes = "collect",
    #         guides = "collect",
    #         heights = c(
    #             0.55,
    #             0.45
    #         )
    #     )
    
    # TODO Change to patchwork?
    p3 <- cowplot::plot_grid(p1, p2, ncol = 1, align = "v", axis = "bt", rel_heights = c(1.5, 1))
    
    file.nm <- paste0(res.folder, save_loc, feature, ".", device)
    
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
        # "p1_data" = p1.data,
        "p2_data" = p2.data
    )
    
    return(data_list)
}

plotFeatPairs <- function(
        feature_chrom = NULL,
        filt.match.row = NULL,
        method = NULL,
        value = NULL,
        filled = FALSE,
        missing = NULL,
        msLevel = 1,
        save_pairs_loc = NULL,
        device = "pdf"
) {
    
    ft.pair <- filt.match.row[["pair"]][[1]]
    
    feat1 <- plotFeatChrInt(
        feature_chrom = feature_chrom,
        feature = ft.pair[1],
        method = method,
        value = value,
        filled = filled,
        missing = missing,
        msLevel = msLevel,
        save_loc = NULL,
        device = NULL,
        feat_pairs = TRUE
    )
    
    feat2 <- plotFeatChrInt(
        feature_chrom = feature_chrom,
        feature = ft.pair[2],
        method = method,
        value = value,
        filled = filled,
        missing = missing,
        msLevel = msLevel,
        save_loc = NULL,
        device = NULL,
        feat_pairs = TRUE
    )
    
    # This is for the old ggplot + ggplot version
    # ft.pair.p <- (
    #     feat1$chromatogram + ggplot2::labs(title = "") | 
    #         feat2$chromatogram + ggplot2::labs(title = "")
    # ) /
    #     (
    #         feat1$boxplot + ggplot2::labs(caption = "") | 
    #             feat2$boxplot + ggplot2::labs(caption = "")
    #     ) +
    #     patchwork::plot_layout(guides = "collect") +
    #     patchwork::plot_annotation(
    #         title = paste0(
    #             ft.pair[1], " & ", 
    #             ft.pair[2], "\n",
    #             "Potential ", filt.match.row$name, " ", filt.match.row$chem_change
    #         )
    #     )
    
    # This is for the base + ggplot version
    ft.pair.p <- ((
        cowplot::ggdraw(feat1$chromatogram) /
            (feat1$boxplot + ggplot2::labs(caption = ""))
    ) |
            (
                cowplot::ggdraw(feat2$chromatogram) /
                    (feat2$boxplot + ggplot2::labs(caption = ""))
            )) +
        patchwork::plot_layout(
            guides = "collect",
            axes = "collect",
            heights = c(1.3, 0.7)
        ) +
        patchwork::plot_annotation(
            title = paste0(
                ft.pair[1], " & ", 
                ft.pair[2], "\n",
                "Potential ", filt.match.row$name, " ", filt.match.row$chem_change
            )
        )
    
    file.nm <- paste0(res.folder, save_pairs_loc, ft.pair[1], "_", ft.pair[2],".", device)
    
    if (!is.null(save_pairs_loc)) {
        ggplot2::ggsave(
            filename = file.nm,
            plot = ft.pair.p,
            device = device,
            height = 7,
            width = 7,
            units = "in"
        )
    }
    
    return(ft.pair.p)
}

feat_to_idx <- function(feature.idx = NULL) {
    clean.idx <- as.numeric(gsub("[A-Za-z]", "", feature.idx))
    return(clean.idx)
}

# TODO FIX This so the path isn't hardcoded
check_saved <- function(
        filename = NULL
        ) {
    object.bool <- file.exists(file.path(res.folder, "objects", filename))
    return(object.bool)
}

registerParallel <- function(
        workers = NULL
        ) {
    sys <- Sys.info()["sysname"]
    if (sys == "Windows") {
        bp <<- BiocParallel::SnowParam(
            workers = workers,
            type = "SOCK"
        )
    } else if (sys %in% c("Linux", "Darwin")) {
        bp <<- BiocParallel::MulticoreParam(
            workers = workers
        )
    }
    
    BiocParallel::register(bp, default = TRUE)
}

filtFeatures <- function(
        object = NULL,
        beta_cor_threshold = 0.3,
        beta_snr_threshold = 6,
        sn_threshold = 10,
        filt_vector = NULL
        ) {
    
    # Poor peaks: beta_cor < 0.5 (or even < 0.2
    # Good peaks: beta_snr > 7
    # Keep signal to noise at 10, and filter by beta_cor < 0.5, and beta_snr > 7?
    
    filt.chrompeaks.tib <- tibble::as_tibble(xcms::chromPeaks(object), rownames = "feature") %>% 
        dplyr::filter(!is.na(beta_cor) & !is.na(beta_snr)) %>%
        dplyr::group_by(feature) %>%
        dplyr::filter(beta_cor >= beta_cor_threshold & beta_snr >= beta_snr_threshold) %>%
        dplyr::filter(sn >= sn_threshold) %>%
        dplyr::mutate(feature2 = as.numeric(gsub("[A-Za-z]", "", feature))) %>%
        dplyr::relocate(feature2, .after = feature)
    
    filt.features.tib <- tibble::as_tibble(xcms::featureDefinitions(object), rownames = "feature") %>%
        tidyr::unnest(peakidx) %>%
        dplyr::filter(peakidx %in% filt.chrompeaks.tib$feature2) %>%
        tidyr::nest(data = peakidx) 
    
    filt.sig.features.tib <- filt.features.tib %>%
        # Added so we don't remove interesting ones 
        # that don't have a predicted biotransformation
        dplyr::filter(feature %in% filt_vector)
    
    filt.sig.features <- filt.sig.features.tib$feature
    
    # Now only look at features that have at least one significantly different feature
    biot.filt.sig.feature.tib <- matched.diffs %>%
        dplyr::filter(
            dplyr::if_any( # if_any for if only one is significant
                dplyr::all_of(c("feat1", "feat2")),
                ~ .x %in% filt.sig.features
            )
        ) %>%
        dplyr:::mutate(
            pair = purrr::map2(feat1, feat2, ~ c(.x, .y)),
            # or ppm global
            mz1_forms = purrr::map(mz1, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 0))),
            mz2_forms = purrr::map(mz2, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 0)))
        ) # %>%
    # dplyr::filter(grepl("1 x", name))
    
    # TODO OLD version where I filder all.int.comps in the biotransformation tibble
    # filt.sig.features <- matched.diffs %>%
    #     dplyr::filter(
    #         dplyr::if_any(
    #             dplyr::all_of(c("feat1", "feat2")),
    #             ~ .x %in% filt.features$feature
    #         )
    #     ) %>%
    #     dplyr::filter(
    #         dplyr::if_any(
    #             dplyr::all_of(c("feat1", "feat2")),
    #             ~ .x %in% all.int.comps
    #         )
    #     ) %>%
    #     dplyr:::mutate(
    #         pair = purrr::map2(feat1, feat2, ~ c(.x, .y)),
    #         # or ppm global
    #         mz1_forms = purrr::map(mz1, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 0))),
    #         mz2_forms = purrr::map(mz2, ~ Rdisop::getFormula(Rdisop::decomposeMass(.x, ppm = 0)))
    #     ) # %>%
    # # dplyr::filter(grepl("1 x", name))
    
    biot.filt.sig.features <- unique(c(
        biot.filt.sig.feature.tib$feat1, 
        biot.filt.sig.feature.tib$feat2
        ))
    
    final.plotting.features <- unique(c(filt.sig.features, biot.filt.sig.features))
    
    filt.list <- list(
        # filtering features
        "filt.chrompeaks.tib" = filt.chrompeaks.tib, # quality filtered peak tib
        "filt.features.tib" = filt.features.tib, # quality filtered feature tib
        "filt.sig.features.tib" = filt.sig.features.tib, # quality + sig filtered feature tib
        "filt.sig.features" = filt.sig.features, # quality + sig filtered features
        # biotransformation features
        "biot.filt.sig.features.tib" = biot.filt.sig.feature.tib, # quality + sig filtered biotransf tib
        "biot.filt.sig.features" = biot.filt.sig.features, # quality + sig filtered feature
        # final features for plotting
        "final.plotting.features" = final.plotting.features # all feats in biot and in filt.sig.features
        
    )
    
    return(filt.list)
    
}

predictBiotransfSubset <- function(
        data = NULL,
        biotransf.data = NULL,
        tolerance_ppm = NULL,
        tolerance = NULL,
        feat_filt = NULL
) {
    
    if (is.null(tolerance) & is.null(tolerance_ppm)) {
        stop("Both tolerance and tolerance_ppm are NULL")
    } else if (!is.null(tolerance) & !is.null(tolerance_ppm)) {
        stop("Either tolerance or tolerance_ppm need to be set to NULL")
    }
    
    if (is.null(tolerance)) {
        tol_used = ppm_to_num(tolerance_ppm)
    }
    if (is.null(tolerance_ppm)) {
        tol_used = tolerance
    }
    
    ## 1. Prepare peaks table (all peaks, sorted by m/z)
    peaks <- data %>%
        # needs to be sorted for findInterval() indexing
        dplyr::arrange(mzmed) %>%
        dplyr::mutate(
            peak_id = dplyr::row_number(),    # simple integer ID
            feature,
            mzmed
        )
    n_peaks <- nrow(peaks)
    
    # convenience vectors so we don't keep indexing peaks$...
    mz_vec <- peaks$mzmed
    rt_vec <- peaks$rtmed
    id_vec <- peaks$peak_id
    feat_vec <- peaks$feature
    n_trans <- nrow(biotransf.data)
    
    ## ONLY do comparisons for FT02089
    anchor_i <- which(feat_vec %in% feat_filt)
    
    ## 3. For each biotransformation, find all matching peak pairs
    all_matches <- vector("list", n_trans)
    
    for (k in 1:n_trans) {
        
        # current transformation
        delta <- biotransf.data$delta_mass[k]
        this_name <- biotransf.data$name[k]
        this_formula <- biotransf.data$chem_formula[k]
        
        # ONLY compute windows for anchor peaks (not all peaks)
        target_lower <- mz_vec[anchor_i] + delta - tol_used
        target_upper <- mz_vec[anchor_i] + delta + tol_used
        
        # ONLY compute interval bounds for anchor peaks
        idx_start <- findInterval(target_lower, mz_vec) + 1L
        idx_end <- findInterval(target_upper, mz_vec)
        
        # ONLY allocate list for anchors
        res_list <- vector("list", length(anchor_i))
        
        for (a in seq_along(anchor_i)) {
            
            i <- anchor_i[a]
            start_idx <- idx_start[a]
            end_idx <- idx_end[a]
            
            # no overlap for this i -> skip
            if (start_idx > end_idx) next
            
            # candidate partner indices
            partner_idx <- seq.int(start_idx, end_idx)
            
            # keep your original duplicate-avoidance logic
            partner_idx <- partner_idx[partner_idx > i]
            if (length(partner_idx) == 0) next
            
            res_list[[a]] <- tibble::tibble(
                name = this_name,
                chem_change = this_formula,
                delta_mass = delta,
                feat1 = feat_vec[i],
                feat2 = feat_vec[partner_idx],
                mz1 = mz_vec[i],
                mz2 = mz_vec[partner_idx],
                rt1 = rt_vec[i],
                rt2 = rt_vec[partner_idx],
                obs_delta_mass = mz_vec[partner_idx] - mz_vec[i],
                peak1_id = id_vec[i],
                peak2_id = id_vec[partner_idx]
            )
        }
        
        all_matches[[k]] <- dplyr::bind_rows(res_list)
    }
    
    matched.diffs <- dplyr::bind_rows(all_matches)
    
    return(matched.diffs)
}

predictBiotransfAdductsSubset <- function(
        data = NULL,
        biotransf.data = NULL,
        tolerance_ppm = NULL,
        tolerance = NULL,
        feat_filt = NULL
) {
    
    if (is.null(tolerance) & is.null(tolerance_ppm)) {
        stop("Both tolerance and tolerance_ppm are NULL")
    } else if (!is.null(tolerance) & !is.null(tolerance_ppm)) {
        stop("Either tolerance or tolerance_ppm need to be set to NULL")
    }
    
    if (is.null(tolerance)) {
        tol_used = ppm_to_num(tolerance_ppm)
    }
    if (is.null(tolerance_ppm)) {
        tol_used = tolerance
    }
    
    ## 1. Prepare peaks table (all peaks, sorted by m/z)
    peaks <- data %>%
        # needs to be sorted for findInterval() indexing
        dplyr::arrange(mass) %>%
        dplyr::mutate(
            peak_id = dplyr::row_number(),    # simple integer ID
            feature,
            mass
        )
    n_peaks <- nrow(peaks)
    
    # convenience vectors so we don't keep indexing peaks$...
    mz_vec <- peaks$mzmed
    rt_vec <- peaks$rtmed
    adduct_vec <- peaks$adduct
    mass_vec <- peaks$mass
    id_vec <- peaks$peak_id
    feat_vec <- peaks$feature
    n_trans <- nrow(biotransf.data)
    
    ## ONLY do comparisons for feat_filt
    anchor_i <- which(feat_vec %in% feat_filt)
    
    ## 3. For each biotransformation, find all matching peak pairs
    all_matches <- vector("list", n_trans)
    
    for (k in 1:n_trans) {
        
        # current transformation
        delta <- biotransf.data$delta_mass[k]
        this_name <- biotransf.data$name[k]
        this_formula <- biotransf.data$chem_formula[k]
        
        # for each peak i, valid partners j must have:
        # mz_vec[j] in [mz_vec[i] + delta - tol, mz_vec[i] + delta + tol]
        target_lower <- mass_vec[anchor_i] + delta - tol_used
        target_upper <- mass_vec[anchor_i] + delta + tol_used
        
        # find, for each i, the index range [start_i, end_i] in the sorted mz_vec
        # that lies within [target_lower[i], target_upper[i]]
        idx_start <- findInterval(target_lower, mass_vec) + 1L
        idx_end <- findInterval(target_upper, mass_vec)
        
        # ONLY allocate list for anchors
        res_list <- vector("list", length(anchor_i))
        
        for (a in seq_along(anchor_i)) {
            
            i <- anchor_i[a]
            start_idx <- idx_start[a]
            end_idx <- idx_end[a]
            
            # no overlap for this i -> skip
            if (start_idx > end_idx) next
            
            # candidate partner indices
            partner_idx <- seq.int(start_idx, end_idx)
            
            # keep your original duplicate-avoidance logic
            partner_idx <- partner_idx[partner_idx > i]
            if (length(partner_idx) == 0) next
            
            res_list[[a]] <- tibble::tibble(
                name = this_name,
                chem_change = this_formula,
                feat1 = feat_vec[i],
                feat2 = feat_vec[partner_idx],
                mz1 = mz_vec[i],
                mz2 = mz_vec[partner_idx],
                
                adduct1 = adduct_vec[i],
                adduct2 = adduct_vec[partner_idx],
                mass1 = mass_vec[i],
                mass2 = mass_vec[partner_idx],
                
                rt1 = rt_vec[i],
                rt2 = rt_vec[partner_idx],
                delta_mass = delta,
                obs_delta_mass = mass_vec[partner_idx] - mass_vec[i],
                peak1_id = id_vec[i],
                peak2_id = id_vec[partner_idx]
            )
        }
        
        # bind all i-level results for this transformation
        all_matches[[k]] <- dplyr::bind_rows(res_list)
    }
    
    # 4. Final table of all matched pairs
    matched.diffs <- dplyr::bind_rows(all_matches)
    
    return(matched.diffs)
}

# predictBiotransfAdductsSubset <- function(
#         data = NULL,
#         biotransf.data = NULL,
#         tolerance_ppm = NULL,
#         tolerance = NULL,
#         feat_filt = NULL
# ) {
    
#     if (is.null(tolerance) & is.null(tolerance_ppm)) {
#         stop("Both tolerance and tolerance_ppm are NULL")
#     } else if (!is.null(tolerance) & !is.null(tolerance_ppm)) {
#         stop("Either tolerance or tolerance_ppm need to be set to NULL")
#     }
    
#     if (is.null(tolerance)) {
#         tol_used = ppm_to_num(tolerance_ppm)
#     }
#     if (is.null(tolerance_ppm)) {
#         tol_used = tolerance
#     }
    
#     ## 1. Prepare peaks table (all peaks, sorted by m/z)
#     peaks <- data %>%
#         # needs to be sorted for findInterval() indexing
#         dplyr::arrange(mass) %>%
#         dplyr::mutate(
#             peak_id = dplyr::row_number(),    # simple integer ID
#             feature,
#             mass
#         )
#     n_peaks <- nrow(peaks)
    
#     # convenience vectors so we don't keep indexing peaks$...
#     mz_vec <- peaks$mzmed
#     rt_vec <- peaks$rtmed
#     adduct_vec <- peaks$adduct
#     mass_vec <- peaks$mass
#     id_vec <- peaks$peak_id
#     feat_vec <- peaks$feature
#     n_trans <- nrow(biotransf.data)
    
#     ## ONLY do comparisons for feat_filt
#     anchor_i <- which(feat_vec %in% feat_filt)
    
#     ## 3. For each biotransformation, find all matching peak pairs
#     all_matches <- vector("list", n_trans)
    
#     for (k in 1:n_trans) {
        
#         # current transformation
#         delta <- biotransf.data$delta_mass[k]
#         this_name <- biotransf.data$name[k]
#         this_formula <- biotransf.data$chem_formula[k]
        
#         # for each peak i, valid partners j must have:
#         # mz_vec[j] in [mz_vec[i] + delta - tol, mz_vec[i] + delta + tol]
#         target_lower <- mass_vec[anchor_i] + delta - tol_used
#         target_upper <- mass_vec[anchor_i] + delta + tol_used
        
#         # find, for each i, the index range [start_i, end_i] in the sorted mz_vec
#         # that lies within [target_lower[i], target_upper[i]]
#         idx_start <- findInterval(target_lower, mass_vec) + 1L
#         idx_end <- findInterval(target_upper, mass_vec)
        
#         # ONLY allocate list for anchors
#         res_list <- vector("list", length(anchor_i))
        
#         for (a in seq_along(anchor_i)) {
            
#             i <- anchor_i[a]
#             start_idx <- idx_start[a]
#             end_idx <- idx_end[a]
            
#             # no overlap for this i -> skip
#             if (start_idx > end_idx) next
            
#             # candidate partner indices
#             partner_idx <- seq.int(start_idx, end_idx)
            
#             # keep your original duplicate-avoidance logic
#             partner_idx <- partner_idx[partner_idx > i]
#             if (length(partner_idx) == 0) next
            
#             res_list[[a]] <- tibble::tibble(
#                 name = this_name,
#                 chem_change = this_formula,
#                 feat1 = feat_vec[i],
#                 feat2 = feat_vec[partner_idx],
#                 mz1 = mz_vec[i],
#                 mz2 = mz_vec[partner_idx],
                
#                 adduct1 = adduct_vec[i],
#                 adduct2 = adduct_vec[partner_idx],
#                 mass1 = mass_vec[i],
#                 mass2 = mass_vec[partner_idx],
                
#                 rt1 = rt_vec[i],
#                 rt2 = rt_vec[partner_idx],
#                 delta_mass = delta,
#                 obs_delta_mass = mass_vec[partner_idx] - mass_vec[i],
#                 peak1_id = id_vec[i],
#                 peak2_id = id_vec[partner_idx]
#             )
#         }
        
#         # bind all i-level results for this transformation
#         all_matches[[k]] <- dplyr::bind_rows(res_list)
#     }
    
#     # 4. Final table of all matched pairs
#     matched.diffs <- dplyr::bind_rows(all_matches)
    
#     return(matched.diffs)
# }

predictBiotransfAdducts <- function(
        data = NULL,
        biotransf.data = NULL,
        tolerance_ppm = NULL,
        tolerance = NULL
) {
    
    if (is.null(tolerance) & is.null(tolerance_ppm)) {
        stop("Both tolerance and tolerance_ppm are NULL")
    } else if (!is.null(tolerance) & !is.null(tolerance_ppm)) {
        stop("Either tolerance or tolerance_ppm need to be set to NULL")
    }
    
    if (is.null(tolerance)) {
        tol_used = ppm_to_num(tolerance_ppm)
    }
    if (is.null(tolerance_ppm)) {
        tol_used = tolerance
    }
    
    ## 1. Prepare peaks table (all peaks, sorted by m/z)
    peaks <- data %>%
        # needs to be sorted for findInterval() indexing
        dplyr::arrange(mass) %>%
        dplyr::mutate(
            peak_id = dplyr::row_number(),    # simple integer ID
            feature,
            mass
        )
    n_peaks <- nrow(peaks)

    # convenience vectors so we don't keep indexing peaks$...
    mz_vec <- peaks$mzmed
    rt_vec <- peaks$rtmed
    adduct_vec <- peaks$adduct
    mass_vec <- peaks$mass
    id_vec <- peaks$peak_id
    feat_vec <- peaks$feature
    n_trans <- nrow(biotransf.data)

    ## 3. For each biotransformation, find all matching peak pairs
    all_matches <- vector("list", n_trans)

    for (k in 1:n_trans) {

        # current transformation
        delta <- biotransf.data$delta_mass[k]
        this_name <- biotransf.data$name[k]
        this_formula <- biotransf.data$chem_formula[k]

        # for each peak i, valid partners j must have:
        # mass_vec[j] in [mass_vec[i] + delta - tol, mass_vec[i] + delta + tol]
        target_lower <- mass_vec + delta - tol_used
        target_upper <- mass_vec + delta + tol_used

        # find, for each i, the index range [start_i, end_i] in the sorted mass_vec
        # that lies within [target_lower[i], target_upper[i]]
        idx_start <- findInterval(target_lower, mass_vec) + 1L
        idx_end <- findInterval(target_upper, mass_vec)

        # collect matches for this transformation here
        res_list <- vector("list", n_peaks)

        for (i in 1:n_peaks) {
            start_idx <- idx_start[i]
            end_idx <- idx_end[i]

            # no overlap for this i -> skip
            if (start_idx > end_idx) next

            # candidate partner indices
            partner_idx <- seq.int(start_idx, end_idx)

            # avoid self-pairs and symmetric duplicates (i,j) vs (j,i)
            partner_idx <- partner_idx[partner_idx > i]
            if (length(partner_idx) == 0) next

            # build rows for all partners of peak i
            res_list[[i]] <- tibble::tibble(
                name = this_name,
                chem_change = this_formula,
                feat1 = feat_vec[i],
                feat2 = feat_vec[partner_idx],
                mz1 = mz_vec[i],
                mz2 = mz_vec[partner_idx],

                adduct1 = adduct_vec[i],
                adduct2 = adduct_vec[partner_idx],
                mass1 = mass_vec[i],
                mass2 = mass_vec[partner_idx],

                rt1 = rt_vec[i],
                rt2 = rt_vec[partner_idx],
                delta_mass = delta,
                obs_delta_mass = mass_vec[partner_idx] - mass_vec[i],
                peak1_id = id_vec[i],
                peak2_id = id_vec[partner_idx]
            )
        }

        # bind all i-level results for this transformation
        all_matches[[k]] <- dplyr::bind_rows(res_list)
    }

    # 4. Final table of all matched pairs
    matched.diffs <- dplyr::bind_rows(all_matches)

    return(matched.diffs)
}