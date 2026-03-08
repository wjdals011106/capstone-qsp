# =============================================================================
# 04_rules.R
# TNBC QSP Model - Algebraic Rules & Derived Quantities
# These are computed from the state variables at each time point
# =============================================================================

compute_rules <- function(state, parms) {
  s <- as.list(state)
  p <- as.list(parms)

  rules <- list()

  # ---------------------------------------------------------------------------
  # TOTAL CELL COUNTS
  # ---------------------------------------------------------------------------
  rules$C_total_T   <- s$C1_T   + s$C2_T   + s$C3_T   + s$C4_T   + s$C5_T
  rules$C_total_Ln1 <- s$C1_Ln1 + s$C2_Ln1 + s$C3_Ln1 + s$C4_Ln1 + s$C5_Ln1
  rules$C_total_Ln2 <- s$C1_Ln2 + s$C2_Ln2 + s$C3_Ln2 + s$C4_Ln2 + s$C5_Ln2
  rules$C_total_oth <- s$C1_oth + s$C2_oth + s$C3_oth + s$C4_oth + s$C5_oth

  Tj_T <- sapply(1:8, function(j) s[[paste0("T",j,"_T")]])
  Tj_Ln1 <- sapply(1:8, function(j) s[[paste0("T",j,"_Ln1")]])
  Tj_Ln2 <- sapply(1:8, function(j) s[[paste0("T",j,"_Ln2")]])
  Tj_oth <- sapply(1:8, function(j) s[[paste0("T",j,"_oth")]])
  Tj_C   <- sapply(1:8, function(j) s[[paste0("T",j,"_C")]])
  Tj_LN  <- sapply(1:8, function(j) s[[paste0("T",j,"_LN")]])

  rules$Tcyt_total_T   <- sum(Tj_T)
  rules$Tcyt_total_Ln1 <- sum(Tj_Ln1)
  rules$Tcyt_total_Ln2 <- sum(Tj_Ln2)
  rules$Tcyt_total_oth <- sum(Tj_oth)
  rules$Tcyt_total_C   <- sum(Tj_C)
  rules$Tcyt_total_LN  <- sum(Tj_LN)

  rules$T_total_T   <- s$T0_T   + s$Th_T   + rules$Tcyt_total_T
  rules$T_total_Ln1 <- s$T0_Ln1 + s$Th_Ln1 + rules$Tcyt_total_Ln1
  rules$T_total_Ln2 <- s$T0_Ln2 + s$Th_Ln2 + rules$Tcyt_total_Ln2
  rules$T_total_oth <- s$T0_oth  + s$Th_oth  + rules$Tcyt_total_oth
  rules$T_total_C   <- s$T0_C   + s$Th_C   + rules$Tcyt_total_C
  rules$T_total_LN  <- s$T0_LN  + s$Th_LN  + rules$Tcyt_total_LN

  rules$M_total_T   <- s$Mac_M1_T   + s$Mac_M2_T
  rules$M_total_Ln1 <- s$Mac_M1_Ln1 + s$Mac_M2_Ln1
  rules$M_total_Ln2 <- s$Mac_M1_Ln2 + s$Mac_M2_Ln2
  rules$M_total_oth <- s$Mac_M1_oth  + s$Mac_M2_oth

  # ---------------------------------------------------------------------------
  # TUMOR VOLUMES (mL) and DIAMETERS (cm)
  # ---------------------------------------------------------------------------
  compute_vol <- function(C_tot, T_tot, M_tot) {
    max(p$V_Tmin,
        p$V_Tmin + (C_tot * p$vol_cell + T_tot * p$vol_Tcell + M_tot * p$vol_Mcell) / p$Ve_T * 1e-12)
  }

  rules$V_T_vol   <- compute_vol(rules$C_total_T,   rules$T_total_T,   rules$M_total_T)
  rules$V_Ln1_vol <- compute_vol(rules$C_total_Ln1, rules$T_total_Ln1, rules$M_total_Ln1)
  rules$V_Ln2_vol <- compute_vol(rules$C_total_Ln2, rules$T_total_Ln2, rules$M_total_Ln2)
  rules$V_oth_vol <- compute_vol(rules$C_total_oth,  rules$T_total_oth,  rules$M_total_oth)

  vol_to_diam <- function(v) (6 * v / pi)^(1/3)
  rules$d_T   <- vol_to_diam(rules$V_T_vol)
  rules$d_Ln1 <- vol_to_diam(rules$V_Ln1_vol)
  rules$d_Ln2 <- vol_to_diam(rules$V_Ln2_vol)
  rules$d_oth <- vol_to_diam(rules$V_oth_vol)

  # ---------------------------------------------------------------------------
  # RECIST: Sum of longest diameters (SLD)
  # Active metastatic lesions only
  # ---------------------------------------------------------------------------
  start_Ln1   <- if (!is.null(s$start_Ln1)   && is.finite(s$start_Ln1))   s$start_Ln1   else 0
  start_Ln2   <- if (!is.null(s$start_Ln2)   && is.finite(s$start_Ln2))   s$start_Ln2   else 0
  start_other <- if (!is.null(s$start_other) && is.finite(s$start_other)) s$start_other else 0

  active_diams <- c(
    if (isTRUE(start_Ln1   == 1)) rules$d_Ln1 else 0,
    if (isTRUE(start_Ln2   == 1)) rules$d_Ln2 else 0,
    if (isTRUE(start_other == 1)) rules$d_oth  else 0
  )
  rules$SLD <- rules$d_T + sum(active_diams)  # primary + active mets (cm)

  # ---------------------------------------------------------------------------
  # IMMUNE CELL FRACTIONS (for biomarker analysis)
  # ---------------------------------------------------------------------------
  # Immune cell fraction in tumor = (T_total + M_total + MDSC) / (C_total + T_total + M_total + MDSC)
  immune_T   <- rules$T_total_T   + rules$M_total_T   + s$MDSC_T
  immune_Ln1 <- rules$T_total_Ln1 + rules$M_total_Ln1 + s$MDSC_Ln1
  immune_Ln2 <- rules$T_total_Ln2 + rules$M_total_Ln2 + s$MDSC_Ln2
  immune_oth <- rules$T_total_oth  + rules$M_total_oth  + s$MDSC_oth

  total_T   <- rules$C_total_T   + immune_T
  total_Ln1 <- rules$C_total_Ln1 + immune_Ln1

  rules$immune_frac_T   <- if (total_T   > 0) immune_T   / total_T   else 0
  rules$immune_frac_Ln1 <- if (total_Ln1 > 0) immune_Ln1 / total_Ln1 else 0

  # Tcyt / T cells (fraction)
  rules$Tcyt_frac_T   <- if (rules$T_total_T   > 0) rules$Tcyt_total_T   / rules$T_total_T   else 0
  rules$Tcyt_frac_Ln1 <- if (rules$T_total_Ln1 > 0) rules$Tcyt_total_Ln1 / rules$T_total_Ln1 else 0
  rules$Tcyt_frac_C   <- if (rules$T_total_C   > 0) rules$Tcyt_total_C   / rules$T_total_C   else 0
  rules$Tcyt_frac_LN  <- if (rules$T_total_LN  > 0) rules$Tcyt_total_LN  / rules$T_total_LN  else 0

  # Treg fraction in LN
  T_LN_all <- s$T0_LN + s$T0_LNl + rules$T_total_LN
  rules$Treg_frac_LN  <- if (T_LN_all > 0) (s$T0_LN + s$T0_LNl) / T_LN_all else 0
  rules$Treg_dens_LN  <- (s$T0_LN + s$T0_LNl) / (p$V_LN * 1e-3)  # cells/mm^3

  # APC density in LN (cells/mm^3)
  rules$APC_dens_LN  <- (s$mAPC_LN + s$mAPC_LNl) / (p$V_LN * 2 * 1e-3)  # two LNs

  # M2/M1 ratio in tumor
  rules$M2M1_T   <- if (s$Mac_M1_T   > 0) s$Mac_M2_T   / s$Mac_M1_T   else NA
  rules$M2M1_Ln1 <- if (s$Mac_M1_Ln1 > 0) s$Mac_M2_Ln1 / s$Mac_M1_Ln1 else NA

  # ---------------------------------------------------------------------------
  # CANCER CLONE RICHNESS AND DIVERSITY (Shannon Index)
  # ---------------------------------------------------------------------------
  shannon_index <- function(counts) {
    counts <- counts[counts > 0]
    if (length(counts) == 0) return(list(H = 0, S = 0, J = 0))
    total <- sum(counts)
    pk    <- counts / total
    H     <- -sum(pk * log(pk))
    S     <- length(counts)
    J     <- if (S > 1) H / log(S) else 0
    list(H = H, S = S, J = J)
  }

  # Cancer clone diversity in tumor
  cancer_T <- c(s$C1_T, s$C2_T, s$C3_T, s$C4_T, s$C5_T)
  div_C_T  <- shannon_index(cancer_T)
  rules$clone_richness_T  <- div_C_T$S
  rules$clone_shannon_T   <- div_C_T$H
  rules$clone_evenness_T  <- div_C_T$J

  cancer_Ln1 <- c(s$C1_Ln1, s$C2_Ln1, s$C3_Ln1, s$C4_Ln1, s$C5_Ln1)
  div_C_Ln1  <- shannon_index(cancer_Ln1)
  rules$clone_richness_Ln1  <- div_C_Ln1$S
  rules$clone_shannon_Ln1   <- div_C_Ln1$H
  rules$clone_evenness_Ln1  <- div_C_Ln1$J

  # Tcyt diversity in LN (across 8 neo-antigen specificities)
  tcyt_LN  <- sapply(1:8, function(j) s[[paste0("T",j,"_LN")]] + s[[paste0("T",j,"_LNl")]])
  div_Tcyt <- shannon_index(tcyt_LN)
  rules$Tcyt_richness_LN  <- div_Tcyt$S
  rules$Tcyt_shannon_LN   <- div_Tcyt$H
  rules$Tcyt_evenness_LN  <- div_Tcyt$J

  # Exhausted Tcyt fraction in tumor
  Texh_T <- sum(sapply(1:8, function(j) s[[paste0("T",j,"_T")]]))  # exhausted tracked by T_tumor itself
  rules$Texh_frac_T <- if (rules$Tcyt_total_T + Texh_T > 0) Texh_T / (rules$Tcyt_total_T + Texh_T) else 0

  # ---------------------------------------------------------------------------
  # PD-L1 EXPRESSION IN TUMOR
  # ---------------------------------------------------------------------------
  IFNg_T <- s$IFNg_T
  rules$PDL1_tumor <- p$C1_PDL1_base * (1 + p$r_PDL1_IFNg * IFNg_T / (p$IFNg_50_ind * 1e-3 + IFNg_T))

  return(rules)
}


# =============================================================================
# RECIST v1.1 Response Classification
# =============================================================================
classify_recist <- function(SLD_baseline, SLD_t, all_tumors_zero = FALSE) {
  if (all_tumors_zero || SLD_t == 0)       return("CR")
  ratio <- SLD_t / SLD_baseline
  if (ratio <= 0.70)                        return("PR")   # >=30% decrease
  if (ratio >= 1.20)                        return("PD")   # >=20% increase
  return("SD")
}

# Determine best overall response from time series
best_overall_response <- function(SLD_vec, times, SLD_baseline,
                                  min_SD_weeks = 24) {
  responses <- character(length(SLD_vec))
  for (i in seq_along(SLD_vec)) {
    all_zero <- SLD_vec[i] == 0
    responses[i] <- classify_recist(SLD_baseline, SLD_vec[i], all_zero)
  }

  # Check if SD is sustained >= 24 weeks (168 days)
  sd_idx <- which(responses == "SD")
  if (length(sd_idx) > 0) {
    sd_times <- times[sd_idx]
    first_sd <- min(sd_times)
    last_sd  <- max(sd_times)
    if ((last_sd - first_sd) < min_SD_weeks * 7) {
      # Short SD - reclassify as PD if preceded by progression
      # (simplified: only count SD if duration > min_SD_weeks)
      responses[sd_idx] <- "PD"
    }
  }

  # Best response priority: CR > PR > SD > PD
  if ("CR" %in% responses) return("CR")
  if ("PR" %in% responses) return("PR")
  if ("SD" %in% responses) return("SD")
  return("PD")
}
