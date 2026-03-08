# =============================================================================
# 05_initial_conditions.R
# TNBC QSP Model - Initial Conditions
# Based on: adg0289_Data_S1.xlsx - Species sheet
# All metastatic tumors start at near-zero (seeded via events at delay times)
# =============================================================================

define_initial_conditions <- function(p) {

  # Helper: named zero vector
  IC <- c(
    # -------------------------------------------------------------------------
    # PEMBROLIZUMAB (anti-PD1) - all start at 0
    # -------------------------------------------------------------------------
    aPD1_C   = 0,  aPD1_P   = 0,
    aPD1_T   = 0,  aPD1_Ln1 = 0,  aPD1_Ln2 = 0,  aPD1_oth = 0,
    aPD1_LN  = 0,  aPD1_LNl = 0,  aPD1_LNo = 0,

    # -------------------------------------------------------------------------
    # CARRYING CAPACITY K (initial from Desai 2006: ~10300 cells at start)
    # -------------------------------------------------------------------------
    K_T   = 10300,
    K_Ln1 = 10300,
    K_Ln2 = 10300,
    K_oth = 10300,

    # -------------------------------------------------------------------------
    # ANGIOGENIC FACTOR (starts at 0, builds up from cancer secretion)
    # -------------------------------------------------------------------------
    cvas_T   = 0, cvas_Ln1 = 0, cvas_Ln2 = 0, cvas_oth = 0,

    # -------------------------------------------------------------------------
    # CANCER CLONES (PRIMARY TUMOR - from Data_S1 Species sheet)
    # -------------------------------------------------------------------------
    C1_T = 1982278.03,
    C2_T = 45658.47,
    C3_T = 1927.60,
    C4_T = 2488714.56,
    C5_T = 181421.35,

    # METASTATIC TUMORS - start at 1e-6 cells (seeded by events)
    C1_Ln1 = 1e-6, C2_Ln1 = 1e-6, C3_Ln1 = 1e-6, C4_Ln1 = 1e-6, C5_Ln1 = 1e-6,
    C1_Ln2 = 1e-6, C2_Ln2 = 1e-6, C3_Ln2 = 1e-6, C4_Ln2 = 1e-6, C5_Ln2 = 1e-6,
    C1_oth = 1e-6, C2_oth = 1e-6, C3_oth = 1e-6, C4_oth = 1e-6, C5_oth = 1e-6,

    # -------------------------------------------------------------------------
    # CCL2 / MDSC / NO / ArgI - all 0
    # -------------------------------------------------------------------------
    CCL2_T=0, CCL2_Ln1=0, CCL2_Ln2=0, CCL2_oth=0,
    MDSC_T=0, MDSC_Ln1=0, MDSC_Ln2=0, MDSC_oth=0,
    NO_T=0,   NO_Ln1=0,   NO_Ln2=0,   NO_oth=0,
    ArgI_T=0, ArgI_Ln1=0, ArgI_Ln2=0, ArgI_oth=0,

    # -------------------------------------------------------------------------
    # MACROPHAGES - start at 0 (recruited via CCL2)
    # -------------------------------------------------------------------------
    Mac_M1_T=0, Mac_M2_T=0,
    Mac_M1_Ln1=0, Mac_M2_Ln1=0,
    Mac_M1_Ln2=0, Mac_M2_Ln2=0,
    Mac_M1_oth=0, Mac_M2_oth=0,

    # -------------------------------------------------------------------------
    # CYTOKINES
    # -------------------------------------------------------------------------
    # TGFb: baseline in breast tumor (Ivanovic 2003)
    TGFb_T=p$TGFbase, TGFb_Ln1=p$TGFbase, TGFb_Ln2=p$TGFbase, TGFb_oth=p$TGFbase,

    IFNg_T=0, IFNg_Ln1=0, IFNg_Ln2=0, IFNg_oth=0,
    IL12_T=0, IL12_Ln1=0, IL12_Ln2=0, IL12_oth=0,
    IL10_T=0, IL10_Ln1=0, IL10_Ln2=0, IL10_oth=0,

    # -------------------------------------------------------------------------
    # APC MATURATION CYTOKINE c - starts at 0
    # -------------------------------------------------------------------------
    c_T=0, c_Ln1=0, c_Ln2=0, c_oth=0,

    # -------------------------------------------------------------------------
    # APC (immature) - starts at steady-state density (APC0)
    # -------------------------------------------------------------------------
    APC_T   = p$APC0_T      * p$V_Tmin,   # cells (density * volume)
    APC_Ln1 = p$APC0_T_Ln1 * p$V_Tmin,
    APC_Ln2 = p$APC0_T_Ln2 * p$V_Tmin,
    APC_oth = p$APC0_T_other * p$V_Tmin,

    # mAPC - starts at 0
    mAPC_T=0, mAPC_Ln1=0, mAPC_Ln2=0, mAPC_oth=0,
    mAPC_LN=0, mAPC_LNl=0, mAPC_LNo=0,

    # -------------------------------------------------------------------------
    # NAIVE T CELLS (from Data_S1 Species sheet)
    # -------------------------------------------------------------------------
    # Central: nT0=3706.9, nT1=2274.8
    nT0_C = 3706.90,
    nT1_C = 2274.77,

    # Peripheral: nT0=185344.8, nT1=113738.7
    nT0_P = 185344.83,
    nT1_P = 113738.74,

    # Lymph nodes: nT0=155.7, nT1=113.7
    nT0_LN  = 155.69,  nT0_LNl = 155.69,  nT0_LNo = 155.69,
    nT1_LN  = 113.74,  nT1_LNl = 113.74,  nT1_LNo = 113.74,

    # -------------------------------------------------------------------------
    # ACTIVATED T CELLS IN LN - all start at 0
    # -------------------------------------------------------------------------
    aT0_LN=0, aT0_LNl=0,
    aTh_LN=0, aTh_LNl=0,
    T0_LN=0,  T0_LNl=0, T0_LNo=0,
    Th_LN=0,  Th_LNl=0, Th_LNo=0,

    # IL-2 in LN: ~0.00019 nM from Data_S1
    IL2_LN  = 0.00019,
    IL2_LNl = 0.00019,

    # -------------------------------------------------------------------------
    # TREG / TH IN CENTRAL AND PERIPHERAL - start at 0
    # -------------------------------------------------------------------------
    T0_C=0, T0_P=0,
    Th_C=0, Th_P=0,

    # Treg / Th in tumors - start at 0
    T0_T=0, T0_Ln1=0, T0_Ln2=0, T0_oth=0,
    Th_T=0, Th_Ln1=0, Th_Ln2=0, Th_oth=0,

    # -------------------------------------------------------------------------
    # SIMULATION CONTROL FLAGS
    # -------------------------------------------------------------------------
    start      = 1,    # primary tumor active
    start_Ln1  = 0,    # lung met 1: off until seeding
    start_Ln2  = 0,    # lung met 2: off until seeding
    start_other= 0     # other met: off until seeding
  )

  # --- Per-neoantigen CD8 T cells (j = 1..8): all start at 0 ---
  for (j in 1:8) {
    IC[paste0("aT",j,"_LN")]  <- 0
    IC[paste0("aT",j,"_LNl")] <- 0
    IC[paste0("T",j,"_LN")]   <- 0
    IC[paste0("T",j,"_LNl")]  <- 0
    IC[paste0("T",j,"_LNo")]  <- 0
    IC[paste0("T",j,"_C")]    <- 0
    IC[paste0("T",j,"_P")]    <- 0
    IC[paste0("T",j,"_T")]    <- 0
    IC[paste0("T",j,"_Ln1")]  <- 0
    IC[paste0("T",j,"_Ln2")]  <- 0
    IC[paste0("T",j,"_oth")]  <- 0
  }

  # --- Initial ncells for metastatic compartments (used in events) ---
  # These are stored in parameters; defaults here match Data_S1
  if (!is.null(p$ncells_C1_Ln1)) {
    # Already set in parameters - no action needed
  } else {
    # Fallback defaults
    for (comp in c("Ln1","Ln2","other")) {
      for (ci in 1:5) {
        pname <- paste0("ncells_C",ci,"_",comp)
        if (is.null(p[[pname]])) p[[pname]] <- 1e5
      }
    }
    # Also for primary tumor
    for (ci in 1:5) {
      pname <- paste0("ncells_C",ci)
      if (is.null(p[[pname]])) p[[pname]] <- IC[paste0("C",ci,"_T")]
    }
  }

  return(IC)
}

# =============================================================================
# Scale initial conditions from p$initial_met_diameter
# Used during virtual patient generation (VP)
# =============================================================================
scale_IC_by_diameter <- function(IC, p) {
  # Target diameter -> volume -> cancer cell count at K
  d_target <- p$initial_met_diameter  # cm
  V_target <- pi / 6 * d_target^3    # cm^3 = mL

  # Rough estimate: C_total ~ V_target * Ve_T / vol_cell (in mL/um^3 -> need conversion)
  # vol_cell in um^3; 1 mL = 1e12 um^3
  C_est <- V_target * p$Ve_T * 1e12 / p$vol_cell
  # Distribute evenly across 5 clones
  for (comp in c("Ln1","Ln2","oth")) {
    for (ci in 1:5) {
      p[[paste0("ncells_C",ci,"_",ifelse(comp=="oth","other",comp))]] <- C_est / 5
    }
  }
  return(p)
}
