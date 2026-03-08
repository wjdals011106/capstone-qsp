# =============================================================================
# 05_initial_conditions.R
# TNBC QSP Model - Initial Conditions
# Based on: adg0289_Data_S1.xlsx - Species sheet
# All metastatic tumors start at near-zero (seeded via events at delay times)
# =============================================================================

define_initial_conditions <- function(p) {

  # =========================================================================
  # Pre-compute quasi-steady-state values for fast T cell activation compartments.
  # These variables start at 0 in the original model but fill up within ~1 day.
  # Setting them at their SS eliminates the stiffness from relative rates -> infinity.
  # =========================================================================

  # Naive T cell counts in LN compartments (from Data_S1 Species sheet)
  nT0_LN0    <- 155.69    # nT0 per LN compartment
  nT1_LN0    <- 113.74    # nT1 per LN compartment

  # f_act at t=0: mAPC_LN = APC0_LN (baseline mature APC in lymph node)
  f_act0     <- p$APC0_LN / (p$mAPC_50 + p$APC0_LN)

  # TGFb Hill at baseline
  H_TGFb0   <- p$TGFbase / (p$TGFb_50 + p$TGFbase)

  # neoantigen probabilities
  p_vec      <- c(p$p_1, p$p_2, p$p_3, p$p_4, p$p_5, p$p_6, p$p_7, p$p_8)

  # --- LN activation intermediates (fast, time scale ~ 1/(k_pro+k_death) ~ 1 day) ---
  aTh_LN_ss  <- p$k_Th_act * f_act0 * nT0_LN0  / (p$k_T0_pro + p$k_T0_death)
  aT0_LN_ss  <- p$k_T0_act * f_act0 * H_TGFb0 * nT0_LN0 / (p$k_T0_pro + p$k_T0_death)

  # --- LN effector pools (time scale ~ 1/(q_out+k_death) ~ 0.04 day) ---
  Th_LN_ss   <- p$k_T0_pro * aTh_LN_ss * p$N_aTh / (p$q_T0_LN_out + p$k_T0_death)
  Th_LNo_ss  <- p$k_Th_act * f_act0 * nT0_LN0   / (p$q_T0_LN_out + p$k_T0_death)
  T0_LN_ss   <- p$k_T0_pro * aT0_LN_ss * p$N_aT0 / (p$q_T0_LN_out + p$k_T0_death)
  T0_LNo_ss  <- p$k_T0_act * f_act0 * H_TGFb0 * nT0_LN0 / (p$q_T0_LN_out + p$k_T0_death)

  # Per-neoantigen CD8 T in LN (time scale ~ 0.04 day)
  aT_LN_ss   <- p$k_T1_act * f_act0 * nT1_LN0 * p_vec / (p$k_T1_pro + p$k_T1_death)
  T_LN_ss    <- p$k_T1_pro * aT_LN_ss * p$N_aT1 / (p$q_T1_LN_out + p$k_T1_death)
  T_LNo_ss   <- p$k_T1_act * f_act0 * nT1_LN0  * p_vec / (p$q_T1_LN_out + p$k_T1_death)

  # --- Central T cell pools ---
  # SS: q_LN_out * 3 * T_LN_ss = (k_death + q_P_in - q_P_out*q_P_in/(q_P_out+k_death)) * T_C
  # eff_loss = k_death * (1 + q_P_in/(q_P_out+k_death))
  eff_loss_T0C <- p$k_T0_death * (1 + p$q_T0_P_in * 1440 / (p$q_T0_P_out + p$k_T0_death))
  eff_loss_T1C <- p$k_T1_death * (1 + p$q_T1_P_in / (p$q_T1_P_out + p$k_T1_death))

  T0_C_ss    <- p$q_T0_LN_out * (T0_LN_ss + T0_LN_ss + T0_LNo_ss) / eff_loss_T0C
  Th_C_ss    <- p$q_T0_LN_out * (Th_LN_ss + Th_LN_ss + Th_LNo_ss) / eff_loss_T0C
  T_C_ss     <- p$q_T1_LN_out * (T_LN_ss  + T_LN_ss  + T_LNo_ss)  / eff_loss_T1C

  # --- Peripheral T cell pools at SS ---
  # dT0_P = q_T0_P_in * 1440 * T0_C - (q_T0_P_out + k_T0_death) * T0_P = 0
  T0_P_ss    <- p$q_T0_P_in * 1440 * T0_C_ss / (p$q_T0_P_out + p$k_T0_death)
  Th_P_ss    <- p$q_T0_P_in * 1440 * Th_C_ss / (p$q_T0_P_out + p$k_T0_death)
  # dT1_P = q_T1_P_in * T1_C - (q_T1_P_out + k_T1_death) * T1_P = 0
  T_P_ss     <- p$q_T1_P_in * T_C_ss / (p$q_T1_P_out + p$k_T1_death)

  # --- Tumor-infiltrating T cell pools at SS ---
  # Tumor volume at IC (primary tumor = 4.7e6 cells)
  C_total_T0 <- 1982278.03 + 45658.47 + 1927.60 + 2488714.56 + 181421.35  # = 4.7e6
  V_T0_vol   <- max(p$V_Tmin, C_total_T0 * p$vol_cell / p$Ve_T * 1e-12)
  d_T0       <- (6 * V_T0_vol / pi)^(1/3)
  v_T0       <- d_T0^3 * pi / 6
  # Met volumes are V_Tmin (near-zero cancer cells)
  v_met0     <- 1e-6  # ~V_Tmin in cm^3

  # dT_T = q_T_T * v_T * T_C - k_death * T_T - kill_terms (= 0 at t=0, no effectors)
  T_T_ss     <- p$q_T1_T_in     * 1440 * v_T0  * T_C_ss / p$k_T1_death
  T_Ln1_ss   <- p$q_T1_T_in_Ln1 * 1440 * v_met0 * T_C_ss / p$k_T1_death
  T_Ln2_ss   <- p$q_T1_T_in_Ln2 * 1440 * v_met0 * T_C_ss / p$k_T1_death
  T_oth_ss   <- p$q_T1_T_in_other * 1440 * v_met0 * T_C_ss / p$k_T1_death

  # T0/Th in primary tumor at SS
  # dT0_T = q_T0_P_in * 1440 * T0_C * V_T/V_C + k_Th_Treg * H_TGFb_T * Th_T - k_T0_death * T0_T
  # At t=0 with Th_T ~ 0: T0_T_ss = q_T0_P_in * 1440 * T0_C_ss * V_T / V_C / k_T0_death
  # But V_C is central volume in liters - check parameter
  T0_T_ss    <- p$q_T0_P_in * 1440 * T0_C_ss * V_T0_vol / p$V_C / p$k_T0_death
  Th_T_ss    <- p$q_T0_P_in * 1440 * Th_C_ss * V_T0_vol / p$V_C / p$k_T0_death

  # Helper: named zero vector
  IC <- c(
    # -------------------------------------------------------------------------
    # PEMBROLIZUMAB (anti-PD1) - all start at 0
    # -------------------------------------------------------------------------
    aPD1_C   = 0,  aPD1_P   = 0,
    aPD1_T   = 0,  aPD1_Ln1 = 0,  aPD1_Ln2 = 0,  aPD1_oth = 0,
    aPD1_LN  = 0,  aPD1_LNl = 0,  aPD1_LNo = 0,

    # -------------------------------------------------------------------------
    # CARRYING CAPACITY K: initialized at quasi-SS of angiogenesis model
    # K_target = C_total * (1 + f_vas * k_K_g/k_K_d)
    # -------------------------------------------------------------------------
    K_T   = max(4.7e6 * (1 + (p$k_vas_Csec*4.7e6/p$k_vas_deg)/(p$c_vas_50 + p$k_vas_Csec*4.7e6/p$k_vas_deg) * p$k_K_g/p$k_K_d), 1e7),
    K_Ln1 = 1e3,
    K_Ln2 = 1e3,
    K_oth = 1e3,

    # -------------------------------------------------------------------------
    # ANGIOGENIC FACTOR: steady state = k_vas_Csec * C_total / k_vas_deg
    # Primary tumor C_total ~ 4.7e6 cells
    # -------------------------------------------------------------------------
    cvas_T   = p$k_vas_Csec * 4.7e6 / p$k_vas_deg,
    cvas_Ln1 = 0, cvas_Ln2 = 0, cvas_oth = 0,

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
    # APC (immature) - steady state: APC = APC0 (source = k_APC_death * APC0,
    #   loss = k_APC_death * APC => dAPC = 0 when APC = APC0)
    # -------------------------------------------------------------------------
    APC_T   = p$APC0_T,
    APC_Ln1 = p$APC0_T_Ln1,
    APC_Ln2 = p$APC0_T_Ln2,
    APC_oth = p$APC0_T_other,

    # mAPC: steady state = APC0_LN (source = k_mAPC_death*APC0_LN, loss = k_mAPC_death*mAPC)
    mAPC_T=0, mAPC_Ln1=0, mAPC_Ln2=0, mAPC_oth=0,
    mAPC_LN=p$APC0_LN, mAPC_LNl=p$APC0_LN, mAPC_LNo=p$APC0_LN,

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
    # ACTIVATED T CELLS IN LN - initialized at quasi-SS to avoid stiffness
    # (were 0 in original but fill within ~1 day from APC0_LN baseline)
    # -------------------------------------------------------------------------
    aT0_LN=aT0_LN_ss, aT0_LNl=aT0_LN_ss,
    aTh_LN=aTh_LN_ss, aTh_LNl=aTh_LN_ss,
    T0_LN=T0_LN_ss,   T0_LNl=T0_LN_ss,   T0_LNo=T0_LNo_ss,
    Th_LN=Th_LN_ss,   Th_LNl=Th_LN_ss,   Th_LNo=Th_LNo_ss,

    # IL-2 in LN: ~0.00019 nM from Data_S1
    IL2_LN  = 0.00019,
    IL2_LNl = 0.00019,

    # -------------------------------------------------------------------------
    # TREG / TH IN CENTRAL AND PERIPHERAL - initialized at quasi-SS
    # -------------------------------------------------------------------------
    T0_C=T0_C_ss, T0_P=T0_P_ss,
    Th_C=Th_C_ss, Th_P=Th_P_ss,

    # Treg / Th in tumors - initialized at quasi-SS
    T0_T=T0_T_ss, T0_Ln1=0, T0_Ln2=0, T0_oth=0,
    Th_T=Th_T_ss, Th_Ln1=0, Th_Ln2=0, Th_oth=0,

    # -------------------------------------------------------------------------
    # SIMULATION CONTROL FLAGS
    # -------------------------------------------------------------------------
    start      = 1,    # primary tumor active
    start_Ln1  = 0,    # lung met 1: off until seeding
    start_Ln2  = 0,    # lung met 2: off until seeding
    start_other= 0     # other met: off until seeding
  )

  # --- Per-neoantigen CD8 T cells (j = 1..8): initialized at quasi-SS ---
  # All pools (LN, central, peripheral, tumor) at quasi-SS for numerical stability.
  # Met compartments start at 0 (no active metastatic disease at t=0).
  for (j in 1:8) {
    IC[paste0("aT",j,"_LN")]  <- aT_LN_ss[j]
    IC[paste0("aT",j,"_LNl")] <- aT_LN_ss[j]
    IC[paste0("T",j,"_LN")]   <- T_LN_ss[j]
    IC[paste0("T",j,"_LNl")]  <- T_LN_ss[j]
    IC[paste0("T",j,"_LNo")]  <- T_LNo_ss[j]
    IC[paste0("T",j,"_C")]    <- T_C_ss[j]
    IC[paste0("T",j,"_P")]    <- T_P_ss[j]
    IC[paste0("T",j,"_T")]    <- T_T_ss[j]
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
