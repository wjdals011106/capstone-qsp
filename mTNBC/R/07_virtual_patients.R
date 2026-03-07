# =============================================================================
# 07_virtual_patients.R
# TNBC QSP Model - Virtual Patient (VP) Generation
# Latin Hypercube Sampling (LHS) for 127 uncertain parameters
# Based on: adg0289_Data_S1.xlsx - "Parameter distributions for VPs" sheet
# =============================================================================

library(lhs)

source("R/01_parameters.R")
source("R/05_initial_conditions.R")

# =============================================================================
# define_vp_distributions()
# Returns a list of 127 parameter distributions for VP generation
# Format: list(name, distribution, median_or_range, sd_or_NA, unit)
# =============================================================================
define_vp_distributions <- function() {
  dists <- list()

  add_lnorm <- function(name, med, gsd) {
    dists[[length(dists)+1]] <<- list(name=name, dist="lnorm", median=med, gsd=gsd)
  }
  add_lunif <- function(name, lo, hi) {
    dists[[length(dists)+1]] <<- list(name=name, dist="lunif", lo=lo, hi=hi)
  }
  add_unif <- function(name, lo, hi) {
    dists[[length(dists)+1]] <<- list(name=name, dist="unif", lo=lo, hi=hi)
  }

  # 1. Initial met diameter (log-normal)
  add_lnorm("initial_met_diameter", 1.65, 0.3)

  # 2-41. Neoantigen concentrations per clone (5 clones x 8 neoantigens)
  for (ci in 1:5) {
    for (j in 1:8) {
      add_lnorm(paste0("agconc_C",ci,"_",j), 6.75e-14, 1)
    }
  }

  # 42-61. Initial cancer cell counts per clone per compartment
  for (ci in 1:5) {
    add_lnorm(paste0("ncells_C",ci), 9.4e5, 0.3)
  }
  for (comp in c("Ln1","Ln2","other")) {
    for (ci in 1:5) {
      add_lnorm(paste0("ncells_C",ci,"_",comp), 9.4e5, 0.3)
    }
  }

  # 62-64. Seeding times
  add_lnorm("delay_Ln1",   2100, 1)
  add_lnorm("delay_Ln2",   2100, 1)
  add_lnorm("delay_other", 2100, 1)

  # 65-84. Cancer clone growth rates (all compartments)
  for (ci in 1:5) {
    add_lnorm(paste0("k_C",ci,"_growth"), 0.0065, 0.7)
  }
  for (comp in c("Ln1","Ln2","other")) {
    for (ci in 1:5) {
      add_lnorm(paste0("k_C",ci,"_growth_",comp), 0.0065, 0.7)
    }
  }

  # 85. Cancer cell killing rate by T cells
  add_lnorm("k_C_T1", 0.95, 1)

  # 86-93. Neoantigen-MHC Kd (k_P1_d1 ... k_P8_d1)
  for (j in 1:8) {
    add_lnorm(paste0("k_P",j,"_d1"), 27e-9, 1)
  }

  # 94-101. Number of T cell clones per neoantigen
  clones_median <- c(189, 190, 191, 192, 193, 194, 195, 196)
  for (j in 1:8) {
    add_lnorm(paste0("n_T",j,"_clones"), clones_median[j], 0.7)
  }
  # 102. Self-antigen T cell clones
  add_lnorm("n_T0_clones", 63, 0.7)

  # 103-104. IL-2 division numbers
  add_unif("N_IL2_CD8", 10, 12)
  add_unif("N_IL2_CD4", 7, 10)

  # 105. PD-1 on T cells (log-uniform)
  add_lunif("T_PD1_total", 1e4, 6.2e4)

  # 106-107. PDL1 baseline
  add_lunif("C1_PDL1_base", 5e4, 2.6e5)
  add_lunif("APC_PDL1_base", 5e4, 2.6e6)

  # 108. PD1_50
  add_lnorm("PD1_50", 6, 1)

  # 109-110. PDL2 ratios
  add_lunif("r_PDL2C1", 0.01, 0.1)
  add_lunif("r_PDL2APC", 0.01, 0.1)

  # 111. Th->Treg differentiation
  add_lnorm("k_Th_Treg", 0.022, 1)

  # 112. Vasculature growth rate
  add_unif("k_K_g", 2.12, 6.12)

  # 113-116. MDSC recruitment
  add_lnorm("k_MDSC_mig",       1.1e4, 0.6)
  add_lnorm("k_MDSC_mig_Ln1",   1.1e4, 0.6)
  add_lnorm("k_MDSC_mig_Ln2",   1.1e4, 0.6)
  add_lnorm("k_MDSC_mig_other", 6.5e4, 0.6)

  # 117-120. Macrophage recruitment
  add_lnorm("k_Mac_mig",       1.7e5, 0.6)
  add_lnorm("k_Mac_mig_Ln1",   3.5e5, 0.6)
  add_lnorm("k_Mac_mig_Ln2",   3.5e5, 0.6)
  add_lnorm("k_Mac_mig_other", 3.5e5, 0.6)

  # 121. CCL2 secretion
  add_lnorm("k_CCL2_sec", 1.7e-12, 0.6)

  # 122-123. CD47/SIRPa expression
  add_unif("C_CD47", 100, 700)
  add_unif("M_SIRPa", 20, 180)

  # 124. PD-1 on macrophages
  add_lunif("M_PD1_total", 1.5e3, 62e3)

  # 125. M1 phagocytosis rate
  add_lnorm("k_M1_phago", 0.33, 1)

  # 126. SIRPa half-maximal
  add_lnorm("SIRPa_50", 37, 0.5)

  # 127. Macrophage-cancer association
  add_lnorm("K_Mac_C", 2, 1)

  return(dists)
}

# =============================================================================
# sample_from_distribution()
# Converts LHS uniform [0,1] sample to parameter value
# =============================================================================
sample_from_distribution <- function(u, dist_info) {
  d <- dist_info
  if (d$dist == "lnorm") {
    # log-normal: GSD parameterization (geometric SD)
    # ln(X) ~ N(ln(median), ln(gsd)^2) if gsd > 0 is log-SD
    # Here gsd is the log-SD of ln(X)
    if (d$gsd == 0 || d$gsd == 1) {
      return(d$median)  # degenerate: no variability
    }
    mu_log <- log(d$median)
    s_log  <- d$gsd  # treating gsd column as log-scale SD
    return(exp(qnorm(u, mu_log, s_log)))
  } else if (d$dist == "lunif") {
    # log-uniform
    return(exp(qunif(u, log(d$lo), log(d$hi))))
  } else if (d$dist == "unif") {
    # uniform
    return(qunif(u, d$lo, d$hi))
  }
}

# =============================================================================
# generate_virtual_patients()
# Returns a data.frame where each row is a VP parameter set
# =============================================================================
generate_virtual_patients <- function(n_vp = 1000, seed = 42) {

  set.seed(seed)
  dists <- define_vp_distributions()
  n_params <- length(dists)
  param_names <- sapply(dists, function(d) d$name)

  # Latin Hypercube Sampling: n_vp rows x n_params columns, uniform [0,1]
  lhs_mat <- randomLHS(n_vp, n_params)

  # Convert to parameter values
  vp_mat <- matrix(NA, nrow = n_vp, ncol = n_params,
                   dimnames = list(NULL, param_names))
  for (j in seq_len(n_params)) {
    vp_mat[, j] <- sapply(lhs_mat[, j], function(u) sample_from_distribution(u, dists[[j]]))
  }

  as.data.frame(vp_mat)
}

# =============================================================================
# build_vp_params()
# Merges base parameters with VP-specific overrides
# =============================================================================
build_vp_params <- function(vp_row) {
  p <- define_parameters()
  for (nm in names(vp_row)) {
    val <- vp_row[[nm]]
    if (!is.na(val) && !is.null(val)) p[[nm]] <- as.numeric(val)
  }
  # Propagate k_Ci_growth to all compartments (overridden if explicitly set)
  for (ci in 1:5) {
    if (!paste0("k_C",ci,"_growth_Ln1") %in% names(vp_row))
      p[[paste0("k_C",ci,"_growth_Ln1")]] <- p[[paste0("k_C",ci,"_growth")]]
    if (!paste0("k_C",ci,"_growth_Ln2") %in% names(vp_row))
      p[[paste0("k_C",ci,"_growth_Ln2")]] <- p[[paste0("k_C",ci,"_growth")]]
    if (!paste0("k_C",ci,"_growth_other") %in% names(vp_row))
      p[[paste0("k_C",ci,"_growth_other")]] <- p[[paste0("k_C",ci,"_growth")]]
  }
  return(p)
}
