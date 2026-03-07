# =============================================================================
# 08_biomarker_analysis.R
# TNBC QSP Model - Biomarker Analysis
# Computes Response Probability and RIS for each biomarker
# Reproduces Data_S2.xlsx results
# =============================================================================

# =============================================================================
# compute_response_probability()
# For a given biomarker, threshold, and direction, compute response probability
# =============================================================================
compute_response_probability <- function(bm_values, responder, cutoff, direction = "greater") {
  if (direction == "greater") {
    in_group <- bm_values >= cutoff
  } else {
    in_group <- bm_values <= cutoff
  }
  if (sum(in_group) == 0) return(list(rp = NA, n_sub = 0, n_resp = 0, n_nonresp = 0))
  n_resp    <- sum(responder[in_group])
  n_nonresp <- sum(!responder[in_group])
  rp <- n_resp / sum(in_group)
  list(rp = rp, n_sub = sum(in_group), n_resp = n_resp, n_nonresp = n_nonresp)
}

# =============================================================================
# compute_RIS()
# Responder Inclusion Score = (resp_in_sub / total_resp) - (nonresp_in_sub / total_nonresp)
# =============================================================================
compute_RIS <- function(bm_values, responder, cutoff, direction = "greater") {
  if (direction == "greater") {
    in_group <- bm_values >= cutoff
  } else {
    in_group <- bm_values <= cutoff
  }
  total_resp    <- sum(responder)
  total_nonresp <- sum(!responder)
  if (total_resp == 0 || total_nonresp == 0) return(NA)
  n_resp_sub    <- sum(responder & in_group)
  n_nonresp_sub <- sum(!responder & in_group)
  (n_resp_sub / total_resp) - (n_nonresp_sub / total_nonresp)
}

# =============================================================================
# analyze_single_biomarker()
# Scans thresholds for optimal Response Probability and RIS
# =============================================================================
analyze_single_biomarker <- function(bm_values, responder, bm_name,
                                     n_thresholds = 100) {
  valid <- !is.na(bm_values)
  bm_v  <- bm_values[valid]
  resp  <- responder[valid]
  if (length(bm_v) < 10 || var(bm_v) == 0) {
    return(data.frame(
      biomarker = bm_name, cutoff_rp = NA, direction_rp = NA,
      response_prob = NA, RIS_at_rp_cutoff = NA,
      cutoff_ris = NA, direction_ris = NA,
      RIS = NA, resp_prob_at_ris_cutoff = NA,
      n_total = sum(valid), n_resp = sum(resp), n_nonresp = sum(!resp)
    ))
  }

  quantiles <- quantile(bm_v, probs = seq(0.05, 0.95, length.out = n_thresholds))
  thresholds <- unique(quantiles)

  best_rp  <- -Inf;  best_rp_cut <- NA;  best_rp_dir <- NA
  best_ris <- -Inf;  best_ris_cut <- NA; best_ris_dir <- NA
  rp_at_best_ris <- NA;  ris_at_best_rp <- NA

  for (thr in thresholds) {
    for (dir in c("greater","less")) {
      rp_res <- compute_response_probability(bm_v, resp, thr, dir)
      ris    <- compute_RIS(bm_v, resp, thr, dir)

      if (!is.na(rp_res$rp) && rp_res$rp > best_rp) {
        best_rp      <- rp_res$rp
        best_rp_cut  <- thr
        best_rp_dir  <- dir
        ris_at_best_rp <- ris
      }
      if (!is.na(ris) && ris > best_ris) {
        best_ris      <- ris
        best_ris_cut  <- thr
        best_ris_dir  <- dir
        rp_at_best_ris <- compute_response_probability(bm_v, resp, thr, dir)$rp
      }
    }
  }

  data.frame(
    biomarker         = bm_name,
    cutoff_rp         = best_rp_cut,
    direction_rp      = best_rp_dir,
    response_prob     = best_rp,
    RIS_at_rp_cutoff  = ris_at_best_rp,
    cutoff_ris        = best_ris_cut,
    direction_ris     = best_ris_dir,
    RIS               = best_ris,
    resp_prob_at_ris_cutoff = rp_at_best_ris,
    n_total           = sum(valid),
    n_resp            = sum(resp),
    n_nonresp         = sum(!resp),
    stringsAsFactors  = FALSE
  )
}

# =============================================================================
# analyze_all_biomarkers()
# Runs single-biomarker analysis on all pre-defined biomarker columns
# =============================================================================
analyze_all_biomarkers <- function(vp_results_df) {
  # vp_results_df: each row = one VP; columns include biomarkers + "responder"
  responder <- vp_results_df$responder

  # Define biomarker columns to analyze (matching Data_S2 names)
  bm_cols <- c(
    # APC
    "APC_dens_LN",
    # Cancer clone diversity
    "clone_richness_T", "clone_richness_Ln1",
    "clone_shannon_T",  "clone_evenness_T",
    # T cell fractions
    "Tcyt_frac_C", "Tcyt_frac_LN", "Tcyt_frac_T",
    "Texh_frac_T",
    # Treg
    "Treg_frac_LN", "Treg_dens_LN",
    # Immune fractions
    "immune_frac_T", "immune_frac_Ln1",
    # T cell diversity
    "Tcyt_richness_LN", "Tcyt_shannon_LN", "Tcyt_evenness_LN",
    # Macrophages
    "M2M1_T", "M2M1_Ln1",
    # Tumor size
    "d_T", "d_Ln1", "SLD",
    # PD-L1
    "PDL1_tumor"
  )

  bm_cols_valid <- bm_cols[bm_cols %in% names(vp_results_df)]
  if (length(bm_cols_valid) == 0) {
    message("No valid biomarker columns found in vp_results_df")
    return(NULL)
  }

  results <- lapply(bm_cols_valid, function(col) {
    analyze_single_biomarker(vp_results_df[[col]], responder, col)
  })

  res_df <- do.call(rbind, results)
  res_df <- res_df[order(-res_df$response_prob), ]
  rownames(res_df) <- NULL
  res_df$rank_rp  <- seq_len(nrow(res_df))
  res_df_ris <- res_df[order(-res_df$RIS), ]
  res_df_ris$rank_ris <- seq_len(nrow(res_df_ris))

  list(by_resp_prob = res_df,
       by_ris       = res_df_ris[, c("biomarker","rank_ris","RIS","resp_prob_at_ris_cutoff","cutoff_ris","direction_ris")])
}

# =============================================================================
# analyze_biomarker_combinations()
# Pairwise biomarker combination analysis (AND logic)
# =============================================================================
analyze_biomarker_combinations <- function(vp_results_df,
                                           top_n_single = 20) {
  responder <- vp_results_df$responder

  # Use top N biomarkers from single analysis
  single_res <- analyze_all_biomarkers(vp_results_df)
  if (is.null(single_res)) return(NULL)

  top_bms <- single_res$by_resp_prob$biomarker[1:min(top_n_single, nrow(single_res$by_resp_prob))]
  top_cuts <- single_res$by_resp_prob[1:length(top_bms), c("biomarker","cutoff_rp","direction_rp")]

  combo_results <- list()
  pairs <- combn(length(top_bms), 2)

  for (k in seq_len(ncol(pairs))) {
    i <- pairs[1, k];  j <- pairs[2, k]
    bm1 <- top_bms[i]; bm2 <- top_bms[j]
    cut1 <- top_cuts$cutoff_rp[i];   dir1 <- top_cuts$direction_rp[i]
    cut2 <- top_cuts$cutoff_rp[j];   dir2 <- top_cuts$direction_rp[j]

    v1 <- vp_results_df[[bm1]];  v2 <- vp_results_df[[bm2]]
    valid <- !is.na(v1) & !is.na(v2)
    if (sum(valid) < 5) next

    in1 <- if (dir1 == "greater") v1 >= cut1 else v1 <= cut1
    in2 <- if (dir2 == "greater") v2 >= cut2 else v2 <= cut2
    in_both <- in1 & in2 & valid

    n_sub  <- sum(in_both)
    if (n_sub == 0) next
    n_resp_sub <- sum(responder[in_both])
    rp   <- n_resp_sub / n_sub
    ris  <- (sum(responder & in_both) / sum(responder)) -
            (sum(!responder & in_both) / sum(!responder))

    combo_results[[k]] <- data.frame(
      bm1 = bm1, cut1 = cut1, dir1 = dir1,
      bm2 = bm2, cut2 = cut2, dir2 = dir2,
      response_prob = rp, RIS = ris,
      n_subgroup = n_sub,
      n_resp = n_resp_sub,
      stringsAsFactors = FALSE
    )
  }

  if (length(combo_results) == 0) return(NULL)
  combo_df <- do.call(rbind, combo_results[!sapply(combo_results, is.null)])
  combo_df <- combo_df[order(-combo_df$response_prob), ]
  rownames(combo_df) <- NULL
  combo_df
}
