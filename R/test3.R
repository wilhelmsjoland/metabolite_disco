predictBiotransfAdduct <- function(
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