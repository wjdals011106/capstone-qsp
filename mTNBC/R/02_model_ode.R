# =============================================================================
# 02_model_ode.R
# TNBC QSP Model - ODE System
# Based on: Arulraj et al., Sci. Adv. 9, eadg0289 (2023)
# =============================================================================
# Compartment abbreviations:
#   _T   = primary tumor
#   _Ln1 = lung metastatic tumor 1
#   _Ln2 = lung metastatic tumor 2
#   _oth = other metastatic tumor
#   _C   = central (blood)
#   _P   = peripheral
#   _LN  = lymph node draining primary tumor
#   _LNl = LN draining lung mets
#   _LNo = LN draining other met
# =============================================================================

tnbc_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # =========================================================================
    # ALGEBRAIC RULES (computed quantities - updated every time step)
    # =========================================================================

    # --- Tumor volumes (mL) ---
    C_total_T   <- C1_T + C2_T + C3_T + C4_T + C5_T
    C_total_Ln1 <- C1_Ln1 + C2_Ln1 + C3_Ln1 + C4_Ln1 + C5_Ln1
    C_total_Ln2 <- C1_Ln2 + C2_Ln2 + C3_Ln2 + C4_Ln2 + C5_Ln2
    C_total_oth <- C1_oth + C2_oth + C3_oth + C4_oth + C5_oth

    T_total_T   <- T0_T + Th_T + sum(sapply(1:8, function(j) get(paste0("T",j,"_T"))))
    T_total_Ln1 <- T0_Ln1 + Th_Ln1 + sum(sapply(1:8, function(j) get(paste0("T",j,"_Ln1"))))
    T_total_Ln2 <- T0_Ln2 + Th_Ln2 + sum(sapply(1:8, function(j) get(paste0("T",j,"_Ln2"))))
    T_total_oth <- T0_oth + Th_oth + sum(sapply(1:8, function(j) get(paste0("T",j,"_oth"))))

    M_total_T   <- Mac_M1_T + Mac_M2_T
    M_total_Ln1 <- Mac_M1_Ln1 + Mac_M2_Ln1
    M_total_Ln2 <- Mac_M1_Ln2 + Mac_M2_Ln2
    M_total_oth <- Mac_M1_oth + Mac_M2_oth

    V_T_vol   <- max(V_Tmin, V_Tmin + (C_total_T   * vol_cell + T_total_T   * vol_Tcell + M_total_T   * vol_Mcell) / Ve_T * 1e-12)
    V_Ln1_vol <- max(V_Tmin, V_Tmin + (C_total_Ln1 * vol_cell + T_total_Ln1 * vol_Tcell + M_total_Ln1 * vol_Mcell) / Ve_T * 1e-12)
    V_Ln2_vol <- max(V_Tmin, V_Tmin + (C_total_Ln2 * vol_cell + T_total_Ln2 * vol_Tcell + M_total_Ln2 * vol_Mcell) / Ve_T * 1e-12)
    V_oth_vol <- max(V_Tmin, V_Tmin + (C_total_oth * vol_cell + T_total_oth * vol_Tcell + M_total_oth * vol_Mcell) / Ve_T * 1e-12)

    # Tumor diameter (cm): V [mL] = pi/6 * d^3 [cm^3]
    d_T   <- (6 * V_T_vol   / pi)^(1/3)
    d_Ln1 <- (6 * V_Ln1_vol / pi)^(1/3)
    d_Ln2 <- (6 * V_Ln2_vol / pi)^(1/3)
    d_oth <- (6 * V_oth_vol / pi)^(1/3)

    # --- Cytotoxic T cells in each tumor ---
    Tcyt_total_T   <- sum(sapply(1:8, function(j) get(paste0("T",j,"_T"))))
    Tcyt_total_Ln1 <- sum(sapply(1:8, function(j) get(paste0("T",j,"_Ln1"))))
    Tcyt_total_Ln2 <- sum(sapply(1:8, function(j) get(paste0("T",j,"_Ln2"))))
    Tcyt_total_oth <- sum(sapply(1:8, function(j) get(paste0("T",j,"_oth"))))

    # --- TGFb Hill functions ---
    H_TGFb_T   <- TGFb_T   / (TGFb_50 + TGFb_T)
    H_TGFb_Ln1 <- TGFb_Ln1 / (TGFb_50 + TGFb_Ln1)
    H_TGFb_Ln2 <- TGFb_Ln2 / (TGFb_50 + TGFb_Ln2)
    H_TGFb_oth <- TGFb_oth / (TGFb_50 + TGFb_oth)

    H_TGFb_Teff_T   <- TGFb_T   / (TGFb_50_Teff + TGFb_T)
    H_TGFb_Teff_Ln1 <- TGFb_Ln1 / (TGFb_50_Teff + TGFb_Ln1)
    H_TGFb_Teff_Ln2 <- TGFb_Ln2 / (TGFb_50_Teff + TGFb_Ln2)
    H_TGFb_Teff_oth <- TGFb_oth / (TGFb_50_Teff + TGFb_oth)

    # --- MDSC Hill functions ---
    H_MDSC_T   <- (ArgI_T / ArgI_50_Teff + NO_T / NO_50_Teff) / (1 + ArgI_T / ArgI_50_Teff + NO_T / NO_50_Teff)
    H_MDSC_Ln1 <- (ArgI_Ln1 / ArgI_50_Teff + NO_Ln1 / NO_50_Teff) / (1 + ArgI_Ln1 / ArgI_50_Teff + NO_Ln1 / NO_50_Teff)
    H_MDSC_Ln2 <- (ArgI_Ln2 / ArgI_50_Teff + NO_Ln2 / NO_50_Teff) / (1 + ArgI_Ln2 / ArgI_50_Teff + NO_Ln2 / NO_50_Teff)
    H_MDSC_oth <- (ArgI_oth / ArgI_50_Teff + NO_oth / NO_50_Teff) / (1 + ArgI_oth / ArgI_50_Teff + NO_oth / NO_50_Teff)

    # --- PD-1 Hill functions (simplified - fraction of PD-1 ligated) ---
    # H_PD1 ≈ PD1_bound^n / (PD1_50^n + PD1_bound^n)
    # For simplicity, we track occupancy via equilibrium approximation:
    PDL1_T   <- C1_PDL1_base * (1 + r_PDL1_IFNg * IFNg_T / (IFNg_50_ind * 1e-3 + IFNg_T))
    PDL1_Ln1 <- C1_PDL1_base * (1 + r_PDL1_IFNg * IFNg_Ln1 / (IFNg_50_ind * 1e-3 + IFNg_Ln1))
    PDL1_Ln2 <- C1_PDL1_base * (1 + r_PDL1_IFNg * IFNg_Ln2 / (IFNg_50_ind * 1e-3 + IFNg_Ln2))
    PDL1_oth <- C1_PDL1_base * (1 + r_PDL1_IFNg * IFNg_oth / (IFNg_50_ind * 1e-3 + IFNg_oth))

    # PD-1 ligation fraction (Hill, pembrolizumab competitively blocks)
    # free_PDL1 ∝ total_PDL1 / (1 + aPD1/Kd), simplified:
    f_PD1_block_T   <- aPD1_T   / (koff_PD1_aPD1 / kon_PD1_aPD1 * 1e9 + aPD1_T)
    f_PD1_block_Ln1 <- aPD1_Ln1 / (koff_PD1_aPD1 / kon_PD1_aPD1 * 1e9 + aPD1_Ln1)
    f_PD1_block_Ln2 <- aPD1_Ln2 / (koff_PD1_aPD1 / kon_PD1_aPD1 * 1e9 + aPD1_Ln2)
    f_PD1_block_oth <- aPD1_oth / (koff_PD1_aPD1 / kon_PD1_aPD1 * 1e9 + aPD1_oth)

    PDL1_free_T   <- PDL1_T   * (1 - f_PD1_block_T)   / A_cell
    PDL1_free_Ln1 <- PDL1_Ln1 * (1 - f_PD1_block_Ln1) / A_cell
    PDL1_free_Ln2 <- PDL1_Ln2 * (1 - f_PD1_block_Ln2) / A_cell
    PDL1_free_oth <- PDL1_oth * (1 - f_PD1_block_oth)  / A_cell

    H_PD1_T   <- PDL1_free_T^n_PD1   / (PD1_50^n_PD1 + PDL1_free_T^n_PD1)
    H_PD1_Ln1 <- PDL1_free_Ln1^n_PD1 / (PD1_50^n_PD1 + PDL1_free_Ln1^n_PD1)
    H_PD1_Ln2 <- PDL1_free_Ln2^n_PD1 / (PD1_50^n_PD1 + PDL1_free_Ln2^n_PD1)
    H_PD1_oth <- PDL1_free_oth^n_PD1  / (PD1_50^n_PD1 + PDL1_free_oth^n_PD1)

    # --- Effective T cell killing rate (inhibited by PD-1, TGFb, MDSC) ---
    k_kill_T   <- k_C_T1 * (1 - H_PD1_T)   * (1 - H_TGFb_Teff_T)   * (1 - H_MDSC_T)
    k_kill_Ln1 <- k_C_T1 * (1 - H_PD1_Ln1) * (1 - H_TGFb_Teff_Ln1) * (1 - H_MDSC_Ln1)
    k_kill_Ln2 <- k_C_T1 * (1 - H_PD1_Ln2) * (1 - H_TGFb_Teff_Ln2) * (1 - H_MDSC_Ln2)
    k_kill_oth <- k_C_T1 * (1 - H_PD1_oth) * (1 - H_TGFb_Teff_oth)  * (1 - H_MDSC_oth)

    # --- IL-10 Hill function for macrophage polarization / phagocytosis ---
    H_IL10_T   <- IL10_T^2   / (IL10_50^2 + IL10_T^2)   # use pM units
    H_IL10_Ln1 <- IL10_Ln1^2 / (IL10_50^2 + IL10_Ln1^2)
    H_IL10_Ln2 <- IL10_Ln2^2 / (IL10_50^2 + IL10_Ln2^2)
    H_IL10_oth <- IL10_oth^2 / (IL10_50^2 + IL10_oth^2)

    H_IL12_T   <- IL12_T^2   / (IL12_50^2 + IL12_T^2)
    H_IL12_Ln1 <- IL12_Ln1^2 / (IL12_50^2 + IL12_Ln1^2)
    H_IL12_Ln2 <- IL12_Ln2^2 / (IL12_50^2 + IL12_Ln2^2)
    H_IL12_oth <- IL12_oth^2 / (IL12_50^2 + IL12_oth^2)

    H_IFNg_T   <- IFNg_T^2   / (IFNg_50^2 + IFNg_T^2)
    H_IFNg_Ln1 <- IFNg_Ln1^2 / (IFNg_50^2 + IFNg_Ln1^2)
    H_IFNg_Ln2 <- IFNg_Ln2^2 / (IFNg_50^2 + IFNg_Ln2^2)
    H_IFNg_oth <- IFNg_oth^2 / (IFNg_50^2 + IFNg_oth^2)

    # --- Phagocytosis inhibition by IL-10 and CD47/SIRPa ---
    H_IL10_phago_T   <- IL10_T   / (IL10_50_phago + IL10_T)
    H_IL10_phago_Ln1 <- IL10_Ln1 / (IL10_50_phago + IL10_Ln1)
    H_IL10_phago_Ln2 <- IL10_Ln2 / (IL10_50_phago + IL10_Ln2)
    H_IL10_phago_oth <- IL10_oth / (IL10_50_phago + IL10_oth)

    SIRPa_syn_T   <- M_SIRPa * A_syn / A_Mcell
    H_SIRPa_T   <- SIRPa_syn_T^n_SIRPa   / (SIRPa_50^n_SIRPa + SIRPa_syn_T^n_SIRPa)
    H_SIRPa_Ln1 <- H_SIRPa_T
    H_SIRPa_Ln2 <- H_SIRPa_T
    H_SIRPa_oth <- H_SIRPa_T

    H_PD1_M_T   <- H_PD1_T
    H_PD1_M_Ln1 <- H_PD1_Ln1
    H_PD1_M_Ln2 <- H_PD1_Ln2
    H_PD1_M_oth <- H_PD1_oth

    k_phago_T   <- k_M1_phago * (1 - H_IL10_phago_T)   * (1 - H_SIRPa_T)   * (1 - H_PD1_M_T)
    k_phago_Ln1 <- k_M1_phago * (1 - H_IL10_phago_Ln1) * (1 - H_SIRPa_Ln1) * (1 - H_PD1_M_Ln1)
    k_phago_Ln2 <- k_M1_phago * (1 - H_IL10_phago_Ln2) * (1 - H_SIRPa_Ln2) * (1 - H_PD1_M_Ln2)
    k_phago_oth <- k_M1_phago * (1 - H_IL10_phago_oth)  * (1 - H_SIRPa_oth) * (1 - H_PD1_M_oth)

    # =========================================================================
    # MODULE 1: PEMBROLIZUMAB PK (2-compartment)
    # =========================================================================
    # Convert dose: 200 mg / 146700 g/mol = 1.363e-3 mol -> in nM given V_C=5L
    # concentration = moles / V_C(L) * 1e9 nM/M
    # Infusion modeled as bolus events (handled in 03_events.R)
    CL_aPD1   <- k_cl_aPD1 / V_C           # clearance per volume [1/day]

    daPD1_C <- - CL_aPD1 * aPD1_C \
               - q_P_aPD1 * 86400 / (V_C * 1000) * (aPD1_C - aPD1_P) \
               - q_T_aPD1 * 86400 / (V_C * 1e6)  * (aPD1_C - aPD1_T) * gamma_T_aPD1 \
               - q_LN_aPD1 * 86400 / (V_C * 1e6) * (aPD1_C - aPD1_LN) * gamma_LN_aPD1

    daPD1_P <- q_P_aPD1 * 86400 / (V_P * 1000) * (aPD1_C - aPD1_P)

    daPD1_T   <- q_T_aPD1   * 86400 / (V_T_vol   * 1e6) * (aPD1_C - aPD1_T)   * gamma_T_aPD1
    daPD1_Ln1 <- q_T_aPD1   * 86400 / (V_Ln1_vol * 1e6) * (aPD1_C - aPD1_Ln1) * gamma_T_aPD1
    daPD1_Ln2 <- q_T_aPD1   * 86400 / (V_Ln2_vol * 1e6) * (aPD1_C - aPD1_Ln2) * gamma_T_aPD1
    daPD1_oth <- q_T_aPD1   * 86400 / (V_oth_vol * 1e6) * (aPD1_C - aPD1_oth) * gamma_T_aPD1
    daPD1_LN  <- q_LN_aPD1  * 86400 / (V_LN * 1e-6) * (aPD1_C - aPD1_LN) * gamma_LN_aPD1 \
                 - q_LD_aPD1 * 1440 * aPD1_LN   # lymphatic drainage [1/min->1/day]
    daPD1_LNl <- q_LN_aPD1  * 86400 / (V_LN * 1e-6) * (aPD1_C - aPD1_LNl) * gamma_LN_aPD1 \
                 - q_LD_aPD1 * 1440 * aPD1_LNl
    daPD1_LNo <- q_LN_aPD1  * 86400 / (V_LN * 1e-6) * (aPD1_C - aPD1_LNo) * gamma_LN_aPD1 \
                 - q_LD_aPD1 * 1440 * aPD1_LNo

    # =========================================================================
    # MODULE 2: ANGIOGENESIS / CARRYING CAPACITY
    # =========================================================================
    # Helper: angiogenic factor
    dcvas_T   <- k_vas_Csec * C_total_T   + k_vas_Msec * Mac_M2_T   - k_vas_deg * cvas_T
    dcvas_Ln1 <- k_vas_Csec * C_total_Ln1 + k_vas_Msec * Mac_M2_Ln1 - k_vas_deg * cvas_Ln1
    dcvas_Ln2 <- k_vas_Csec * C_total_Ln2 + k_vas_Msec * Mac_M2_Ln2 - k_vas_deg * cvas_Ln2
    dcvas_oth <- k_vas_Csec * C_total_oth  + k_vas_Msec * Mac_M2_oth  - k_vas_deg * cvas_oth

    # Carrying capacity K (cells)
    f_vas_T   <- cvas_T   / (c_vas_50 + cvas_T)
    f_vas_Ln1 <- cvas_Ln1 / (c_vas_50 + cvas_Ln1)
    f_vas_Ln2 <- cvas_Ln2 / (c_vas_50 + cvas_Ln2)
    f_vas_oth <- cvas_oth / (c_vas_50 + cvas_oth)

    dK_T   <- k_K_g * f_vas_T   * K_T   - k_K_d * K_T
    dK_Ln1 <- k_K_g * f_vas_Ln1 * K_Ln1 - k_K_d * K_Ln1
    dK_Ln2 <- k_K_g * f_vas_Ln2 * K_Ln2 - k_K_d * K_Ln2
    dK_oth <- k_K_g * f_vas_oth * K_oth  - k_K_d * K_oth

    # =========================================================================
    # MODULE 3: CANCER CELL DYNAMICS (Modified Gompertz)
    # =========================================================================
    # Helper: cancer cell killing rate per Tcyt (Teff-to-C ratio dependent)
    # R_kill_i = k_kill * Tcyt / (K_T_C * C_i + Tcyt) (saturation in effector:target ratio)

    # Per-clone ODE for one compartment
    cancer_clone_ode <- function(Ci, C_total, Tcyt, K, k_growth, k_death,
                                 k_kill, K_T_C, Mac_M1, k_phago) {
      grow  <- k_growth * Ci * log(max(K, 1) / max(C_total + 1, 1)) * (C_total > 0)
      death <- k_death  * Ci
      kill  <- k_kill   * Tcyt / (K_T_C * C_total + Tcyt + 1e-10) * Ci
      phago <- k_phago  * Mac_M1 / (K_Mac_C * C_total + Mac_M1 + 1e-10) * Ci
      grow - death - kill - phago
    }

    # PRIMARY TUMOR
    dC1_T <- start * cancer_clone_ode(C1_T, C_total_T, Tcyt_total_T, K_T, k_C1_growth, k_C1_death, k_kill_T, K_T_C, Mac_M1_T, k_phago_T)
    dC2_T <- start * cancer_clone_ode(C2_T, C_total_T, Tcyt_total_T, K_T, k_C2_growth, k_C2_death, k_kill_T, K_T_C, Mac_M1_T, k_phago_T)
    dC3_T <- start * cancer_clone_ode(C3_T, C_total_T, Tcyt_total_T, K_T, k_C3_growth, k_C3_death, k_kill_T, K_T_C, Mac_M1_T, k_phago_T)
    dC4_T <- start * cancer_clone_ode(C4_T, C_total_T, Tcyt_total_T, K_T, k_C4_growth, k_C4_death, k_kill_T, K_T_C, Mac_M1_T, k_phago_T)
    dC5_T <- start * cancer_clone_ode(C5_T, C_total_T, Tcyt_total_T, K_T, k_C5_growth, k_C5_death, k_kill_T, K_T_C, Mac_M1_T, k_phago_T)

    # LUNG MET 1
    dC1_Ln1 <- start_Ln1 * cancer_clone_ode(C1_Ln1, C_total_Ln1, Tcyt_total_Ln1, K_Ln1, k_C1_growth_Ln1, k_C1_death_Ln1, k_kill_Ln1, K_T_C, Mac_M1_Ln1, k_phago_Ln1)
    dC2_Ln1 <- start_Ln1 * cancer_clone_ode(C2_Ln1, C_total_Ln1, Tcyt_total_Ln1, K_Ln1, k_C2_growth_Ln1, k_C2_death_Ln1, k_kill_Ln1, K_T_C, Mac_M1_Ln1, k_phago_Ln1)
    dC3_Ln1 <- start_Ln1 * cancer_clone_ode(C3_Ln1, C_total_Ln1, Tcyt_total_Ln1, K_Ln1, k_C3_growth_Ln1, k_C3_death_Ln1, k_kill_Ln1, K_T_C, Mac_M1_Ln1, k_phago_Ln1)
    dC4_Ln1 <- start_Ln1 * cancer_clone_ode(C4_Ln1, C_total_Ln1, Tcyt_total_Ln1, K_Ln1, k_C4_growth_Ln1, k_C4_death_Ln1, k_kill_Ln1, K_T_C, Mac_M1_Ln1, k_phago_Ln1)
    dC5_Ln1 <- start_Ln1 * cancer_clone_ode(C5_Ln1, C_total_Ln1, Tcyt_total_Ln1, K_Ln1, k_C5_growth_Ln1, k_C5_death_Ln1, k_kill_Ln1, K_T_C, Mac_M1_Ln1, k_phago_Ln1)

    # LUNG MET 2
    dC1_Ln2 <- start_Ln2 * cancer_clone_ode(C1_Ln2, C_total_Ln2, Tcyt_total_Ln2, K_Ln2, k_C1_growth_Ln2, k_C1_death_Ln2, k_kill_Ln2, K_T_C, Mac_M1_Ln2, k_phago_Ln2)
    dC2_Ln2 <- start_Ln2 * cancer_clone_ode(C2_Ln2, C_total_Ln2, Tcyt_total_Ln2, K_Ln2, k_C2_growth_Ln2, k_C2_death_Ln2, k_kill_Ln2, K_T_C, Mac_M1_Ln2, k_phago_Ln2)
    dC3_Ln2 <- start_Ln2 * cancer_clone_ode(C3_Ln2, C_total_Ln2, Tcyt_total_Ln2, K_Ln2, k_C3_growth_Ln2, k_C3_death_Ln2, k_kill_Ln2, K_T_C, Mac_M1_Ln2, k_phago_Ln2)
    dC4_Ln2 <- start_Ln2 * cancer_clone_ode(C4_Ln2, C_total_Ln2, Tcyt_total_Ln2, K_Ln2, k_C4_growth_Ln2, k_C4_death_Ln2, k_kill_Ln2, K_T_C, Mac_M1_Ln2, k_phago_Ln2)
    dC5_Ln2 <- start_Ln2 * cancer_clone_ode(C5_Ln2, C_total_Ln2, Tcyt_total_Ln2, K_Ln2, k_C5_growth_Ln2, k_C5_death_Ln2, k_kill_Ln2, K_T_C, Mac_M1_Ln2, k_phago_Ln2)

    # OTHER MET
    dC1_oth <- start_other * cancer_clone_ode(C1_oth, C_total_oth, Tcyt_total_oth, K_oth, k_C1_growth_other, k_C1_death_other, k_kill_oth, K_T_C, Mac_M1_oth, k_phago_oth)
    dC2_oth <- start_other * cancer_clone_ode(C2_oth, C_total_oth, Tcyt_total_oth, K_oth, k_C2_growth_other, k_C2_death_other, k_kill_oth, K_T_C, Mac_M1_oth, k_phago_oth)
    dC3_oth <- start_other * cancer_clone_ode(C3_oth, C_total_oth, Tcyt_total_oth, K_oth, k_C3_growth_other, k_C3_death_other, k_kill_oth, K_T_C, Mac_M1_oth, k_phago_oth)
    dC4_oth <- start_other * cancer_clone_ode(C4_oth, C_total_oth, Tcyt_total_oth, K_oth, k_C4_growth_other, k_C4_death_other, k_kill_oth, K_T_C, Mac_M1_oth, k_phago_oth)
    dC5_oth <- start_other * cancer_clone_ode(C5_oth, C_total_oth, Tcyt_total_oth, K_oth, k_C5_growth_other, k_C5_death_other, k_kill_oth, K_T_C, Mac_M1_oth, k_phago_oth)

    # =========================================================================
    # MODULE 4: CCL2, MDSC, NO, ArgI
    # =========================================================================
    k_CCL2_deg_day <- k_CCL2_deg * 24  # convert 1/hr -> 1/day

    dCCL2_T   <- k_CCL2_sec * C_total_T   - k_CCL2_deg_day * CCL2_T
    dCCL2_Ln1 <- k_CCL2_sec * C_total_Ln1 - k_CCL2_deg_day * CCL2_Ln1
    dCCL2_Ln2 <- k_CCL2_sec * C_total_Ln2 - k_CCL2_deg_day * CCL2_Ln2
    dCCL2_oth <- k_CCL2_sec * C_total_oth  - k_CCL2_deg_day * CCL2_oth

    f_CCL2_T   <- CCL2_T   / (CCL2_50 + CCL2_T)
    f_CCL2_Ln1 <- CCL2_Ln1 / (CCL2_50 + CCL2_Ln1)
    f_CCL2_Ln2 <- CCL2_Ln2 / (CCL2_50 + CCL2_Ln2)
    f_CCL2_oth <- CCL2_oth / (CCL2_50 + CCL2_oth)

    dMDSC_T   <- k_MDSC_mig        * f_CCL2_T   - k_MDSC_death * MDSC_T
    dMDSC_Ln1 <- k_MDSC_mig_Ln1    * f_CCL2_Ln1 - k_MDSC_death * MDSC_Ln1
    dMDSC_Ln2 <- k_MDSC_mig_Ln2    * f_CCL2_Ln2 - k_MDSC_death * MDSC_Ln2
    dMDSC_oth <- k_MDSC_mig_other   * f_CCL2_oth  - k_MDSC_death * MDSC_oth

    dNO_T   <- k_NO_sec   * MDSC_T   - k_NO_deg   * NO_T
    dNO_Ln1 <- k_NO_sec   * MDSC_Ln1 - k_NO_deg   * NO_Ln1
    dNO_Ln2 <- k_NO_sec   * MDSC_Ln2 - k_NO_deg   * NO_Ln2
    dNO_oth <- k_NO_sec   * MDSC_oth  - k_NO_deg   * NO_oth

    dArgI_T   <- k_ArgI_sec * MDSC_T   - k_ArgI_deg * ArgI_T
    dArgI_Ln1 <- k_ArgI_sec * MDSC_Ln1 - k_ArgI_deg * ArgI_Ln1
    dArgI_Ln2 <- k_ArgI_sec * MDSC_Ln2 - k_ArgI_deg * ArgI_Ln2
    dArgI_oth <- k_ArgI_sec * MDSC_oth  - k_ArgI_deg * ArgI_oth

    # =========================================================================
    # MODULE 5: MACROPHAGE DYNAMICS (M1 / M2)
    # =========================================================================
    # M1 recruitment (CCL2-dependent), polarization M1<->M2
    dMac_M1_T <- k_Mac_mig * f_CCL2_T \
                 + k_M1_pol * (H_IL12_T + H_IFNg_T) * Mac_M2_T \
                 - k_M2_pol * H_IL10_T * Mac_M1_T \
                 - k_Mac_death * Mac_M1_T

    dMac_M2_T <- k_M2_pol * H_IL10_T * Mac_M1_T \
                 - k_M1_pol * (H_IL12_T + H_IFNg_T) * Mac_M2_T \
                 - k_Mac_death * Mac_M2_T

    dMac_M1_Ln1 <- k_Mac_mig_Ln1 * f_CCL2_Ln1 \
                   + k_M1_pol * (H_IL12_Ln1 + H_IFNg_Ln1) * Mac_M2_Ln1 \
                   - k_M2_pol_Ln1 * H_IL10_Ln1 * Mac_M1_Ln1 \
                   - k_Mac_death * Mac_M1_Ln1

    dMac_M2_Ln1 <- k_M2_pol_Ln1 * H_IL10_Ln1 * Mac_M1_Ln1 \
                   - k_M1_pol * (H_IL12_Ln1 + H_IFNg_Ln1) * Mac_M2_Ln1 \
                   - k_Mac_death * Mac_M2_Ln1

    dMac_M1_Ln2 <- k_Mac_mig_Ln2 * f_CCL2_Ln2 \
                   + k_M1_pol * (H_IL12_Ln2 + H_IFNg_Ln2) * Mac_M2_Ln2 \
                   - k_M2_pol_Ln2 * H_IL10_Ln2 * Mac_M1_Ln2 \
                   - k_Mac_death * Mac_M1_Ln2

    dMac_M2_Ln2 <- k_M2_pol_Ln2 * H_IL10_Ln2 * Mac_M1_Ln2 \
                   - k_M1_pol * (H_IL12_Ln2 + H_IFNg_Ln2) * Mac_M2_Ln2 \
                   - k_Mac_death * Mac_M2_Ln2

    dMac_M1_oth <- k_Mac_mig_other * f_CCL2_oth \
                   + k_M1_pol * (H_IL12_oth + H_IFNg_oth) * Mac_M2_oth \
                   - k_M2_pol_other * H_IL10_oth * Mac_M1_oth \
                   - k_Mac_death * Mac_M1_oth

    dMac_M2_oth <- k_M2_pol_other * H_IL10_oth * Mac_M1_oth \
                   - k_M1_pol * (H_IL12_oth + H_IFNg_oth) * Mac_M2_oth \
                   - k_Mac_death * Mac_M2_oth

    # =========================================================================
    # MODULE 6: CYTOKINES (TGFb, IFNg, IL-12, IL-10)
    # =========================================================================
    # TGFb: secreted by Tregs and M2 macrophages
    dTGFb_T   <- k_TGFb_Tsec * T0_T   + k_TGFb_Msec * Mac_M2_T   - k_TGFb_deg * TGFb_T   + TGFbase * k_TGFb_deg
    dTGFb_Ln1 <- k_TGFb_Tsec * T0_Ln1 + k_TGFb_Msec * Mac_M2_Ln1 - k_TGFb_deg * TGFb_Ln1 + TGFbase * k_TGFb_deg
    dTGFb_Ln2 <- k_TGFb_Tsec * T0_Ln2 + k_TGFb_Msec * Mac_M2_Ln2 - k_TGFb_deg * TGFb_Ln2 + TGFbase * k_TGFb_deg
    dTGFb_oth <- k_TGFb_Tsec * T0_oth  + k_TGFb_Msec * Mac_M2_oth  - k_TGFb_deg * TGFb_oth  + TGFbase * k_TGFb_deg

    # IFNg: secreted by Th cells
    dIFNg_T   <- k_IFNg_sec * Th_T   - k_IFNg_deg * IFNg_T
    dIFNg_Ln1 <- k_IFNg_sec * Th_Ln1 - k_IFNg_deg * IFNg_Ln1
    dIFNg_Ln2 <- k_IFNg_sec * Th_Ln2 - k_IFNg_deg * IFNg_Ln2
    dIFNg_oth <- k_IFNg_sec * Th_oth  - k_IFNg_deg * IFNg_oth

    k_IL12_deg_day <- k_IL12_deg * 24  # 1/hr -> 1/day
    dIL12_T   <- k_IL12_sec * mAPC_T   + k_IL12_Msec * Mac_M1_T   - k_IL12_deg_day * IL12_T
    dIL12_Ln1 <- k_IL12_sec * mAPC_Ln1 + k_IL12_Msec * Mac_M1_Ln1 - k_IL12_deg_day * IL12_Ln1
    dIL12_Ln2 <- k_IL12_sec * mAPC_Ln2 + k_IL12_Msec * Mac_M1_Ln2 - k_IL12_deg_day * IL12_Ln2
    dIL12_oth <- k_IL12_sec * mAPC_oth  + k_IL12_Msec * Mac_M1_oth  - k_IL12_deg_day * IL12_oth

    dIL10_T   <- k_IL10_sec * Mac_M2_T   - k_IL10_deg * IL10_T
    dIL10_Ln1 <- k_IL10_sec * Mac_M2_Ln1 - k_IL10_deg * IL10_Ln1
    dIL10_Ln2 <- k_IL10_sec * Mac_M2_Ln2 - k_IL10_deg * IL10_Ln2
    dIL10_oth <- k_IL10_sec * Mac_M2_oth  - k_IL10_deg * IL10_oth

    # =========================================================================
    # MODULE 7: APC DYNAMICS
    # =========================================================================
    # Maturation cytokine c (DAMPs from dying cancer cells)
    # Rate of cancer cell death -> DAMPs release
    death_rate_T   <- (k_C1_death * C1_T + k_C2_death * C2_T + k_C3_death * C3_T + k_C4_death * C4_T + k_C5_death * C5_T)
    death_rate_Ln1 <- (k_C1_death_Ln1 * C1_Ln1 + k_C2_death_Ln1 * C2_Ln1 + k_C3_death_Ln1 * C3_Ln1 + k_C4_death_Ln1 * C4_Ln1 + k_C5_death_Ln1 * C5_Ln1)
    death_rate_Ln2 <- (k_C1_death_Ln2 * C1_Ln2 + k_C2_death_Ln2 * C2_Ln2 + k_C3_death_Ln2 * C3_Ln2 + k_C4_death_Ln2 * C4_Ln2 + k_C5_death_Ln2 * C5_Ln2)
    death_rate_oth <- (k_C1_death_other * C1_oth + k_C2_death_other * C2_oth + k_C3_death_other * C3_oth + k_C4_death_other * C4_oth + k_C5_death_other * C5_oth)

    dc_T   <- DAMPs * death_rate_T   / V_T_vol   - k_c * c_T
    dc_Ln1 <- DAMPs * death_rate_Ln1 / V_Ln1_vol - k_c * c_Ln1
    dc_Ln2 <- DAMPs * death_rate_Ln2 / V_Ln2_vol - k_c * c_Ln2
    dc_oth <- DAMPs * death_rate_oth  / V_oth_vol  - k_c * c_oth

    f_c_T   <- c_T   / (c50 + c_T)
    f_c_Ln1 <- c_Ln1 / (c50 + c_Ln1)
    f_c_Ln2 <- c_Ln2 / (c50 + c_Ln2)
    f_c_oth <- c_oth / (c50 + c_oth)

    # APC (immature -> mature by DAMPs, IL-12; inhibited by IL-10)
    f_mat_T   <- k_APC_mat * f_c_T   * (1 - H_IL10_T)
    f_mat_Ln1 <- k_APC_mat * f_c_Ln1 * (1 - H_IL10_Ln1)
    f_mat_Ln2 <- k_APC_mat * f_c_Ln2 * (1 - H_IL10_Ln2)
    f_mat_oth <- k_APC_mat * f_c_oth  * (1 - H_IL10_oth)

    dAPC_T   <- APC0_T   * k_APC_death - f_mat_T   * APC_T   - k_APC_death * APC_T
    dAPC_Ln1 <- APC0_T_Ln1 * k_APC_death - f_mat_Ln1 * APC_Ln1 - k_APC_death * APC_Ln1
    dAPC_Ln2 <- APC0_T_Ln2 * k_APC_death - f_mat_Ln2 * APC_Ln2 - k_APC_death * APC_Ln2
    dAPC_oth <- APC0_T_other * k_APC_death - f_mat_oth * APC_oth - k_APC_death * APC_oth

    # mAPC: matured in tumor, migrate to LN
    dmAPC_T   <- f_mat_T   * APC_T   - k_APC_mig * mAPC_T   - k_mAPC_death * mAPC_T
    dmAPC_Ln1 <- f_mat_Ln1 * APC_Ln1 - k_APC_mig * mAPC_Ln1 - k_mAPC_death * mAPC_Ln1
    dmAPC_Ln2 <- f_mat_Ln2 * APC_Ln2 - k_APC_mig * mAPC_Ln2 - k_mAPC_death * mAPC_Ln2
    dmAPC_oth <- f_mat_oth  * APC_oth  - k_APC_mig * mAPC_oth  - k_mAPC_death * mAPC_oth

    # mAPC in LN: arrive from tumor (k_APC_mig * mAPC_T * V_T / V_LN)
    dmAPC_LN  <- k_APC_mig * mAPC_T   * V_T_vol   / (V_LN * 1e-6) + APC0_LN * k_mAPC_death - k_mAPC_death * mAPC_LN
    dmAPC_LNl <- k_APC_mig * mAPC_Ln1 * V_Ln1_vol / (V_LN * 1e-6) + APC0_LN * k_mAPC_death - k_mAPC_death * mAPC_LNl
    dmAPC_LNo <- k_APC_mig * mAPC_oth  * V_oth_vol  / (V_LN * 1e-6) + APC0_LN * k_mAPC_death - k_mAPC_death * mAPC_LNo

    # =========================================================================
    # MODULE 8: NAIVE T CELL TRAFFICKING (Central <-> Peripheral <-> LN)
    # =========================================================================
    # Convert flow rates to 1/day where needed (q in 1/min -> * 1440 min/day)
    q_nT0_P_in_day  <- q_nT0_P_in  * 1440
    q_nT0_P_out_day <- q_nT0_P_out
    q_nT1_P_in_day  <- q_nT1_P_in  * 1440
    q_nT1_P_out_day <- q_nT1_P_out

    # Thymic export + peripheral homeostatic proliferation
    thymic_nT0 <- Q_nT0_thym
    thymic_nT1 <- Q_nT1_thym
    homeo_nT0  <- k_nT0_pro * nT0_C / (K_nT0_pro + nT0_C)
    homeo_nT1  <- k_nT1_pro * nT1_C / (K_nT1_pro + nT1_C)

    dnT0_C <- thymic_nT0 + homeo_nT0 \
              - k_nT0_death * nT0_C \
              - q_nT0_P_in_day * nT0_C \
              + q_nT0_P_out_day * nT0_P \
              - q_nT0_LN_in * nT0_C \
              + q_nT0_LN_out * nT0_LN \
              + q_nT0_LN_out * nT0_LNl

    dnT0_P <- q_nT0_P_in_day * nT0_C - q_nT0_P_out_day * nT0_P - k_nT0_death * nT0_P

    dnT0_LN  <- q_nT0_LN_in * nT0_C - q_nT0_LN_out * nT0_LN  - k_nT0_death * nT0_LN
    dnT0_LNl <- q_nT0_LN_in * nT0_C - q_nT0_LN_out * nT0_LNl - k_nT0_death * nT0_LNl
    dnT0_LNo <- q_nT0_LN_in * nT0_C - q_nT0_LN_out * nT0_LNo - k_nT0_death * nT0_LNo

    dnT1_C <- thymic_nT1 + homeo_nT1 \
              - k_nT1_death * nT1_C \
              - q_nT1_P_in_day * nT1_C \
              + q_nT1_P_out_day * nT1_P \
              - q_nT1_LN_in * nT1_C \
              + q_nT1_LN_out * nT1_LN

    dnT1_P <- q_nT1_P_in_day * nT1_C - q_nT1_P_out_day * nT1_P - k_nT1_death * nT1_P

    dnT1_LN  <- q_nT1_LN_in * nT1_C - q_nT1_LN_out * nT1_LN  - k_nT1_death * nT1_LN
    dnT1_LNl <- q_nT1_LN_in * nT1_C - q_nT1_LN_out * nT1_LNl - k_nT1_death * nT1_LNl
    dnT1_LNo <- q_nT1_LN_in * nT1_C - q_nT1_LN_out * nT1_LNo - k_nT1_death * nT1_LNo

    # =========================================================================
    # MODULE 9: T CELL ACTIVATION IN LN (per neo-antigen j = 1..8 and self j=0)
    # =========================================================================
    # Activation signal: TCR ligation ~ pMHC_j / (p_j_50 + pMHC_j)
    # Simplified: we track active T cells directly with lumped activation rates

    # IL-2 in LN (drives proliferation)
    T_total_LN_all  <- sum(sapply(1:8, function(j) get(paste0("aT",j,"_LN")))) + T0_LN + Th_LN + aT0_LN
    T_total_LNl_all <- sum(sapply(1:8, function(j) get(paste0("aT",j,"_LNl")))) + T0_LNl + Th_LNl + aT0_LNl

    dIL2_LN <- k_IL2_sec * T_total_LN_all  - k_IL2_deg * 1440 * IL2_LN  - k_IL2_cons * T_total_LN_all  * IL2_LN
    dIL2_LNl<- k_IL2_sec * T_total_LNl_all - k_IL2_deg * 1440 * IL2_LNl - k_IL2_cons * T_total_LNl_all * IL2_LNl

    # IL-2 signal for division
    f_IL2_LN  <- IL2_LN  / (IL2_50 + IL2_LN)
    f_IL2_LNl <- IL2_LNl / (IL2_50 + IL2_LNl)

    # Activation + expansion per CD8 T cell clone j in each LN
    # mAPC drives activation (Hill function); N_IL2_CD8 * f_IL2 drives expansion
    f_act_LN  <- mAPC_LN  / (n_sites_APC * APC0_LN + mAPC_LN)   # simplified pMHC Hill
    f_act_LNl <- mAPC_LNl / (n_sites_APC * APC0_LN + mAPC_LNl)

    # T Treg (CD4 Treg) activation in LN
    # Uses self-antigen presentation (H_APC ~ self-pMHC) and TGFb
    # Simplified: activation proportional to mAPC * TGFb signal

    # CD4 Treg (T0) in LN
    daT0_LN  <- k_T0_act * f_act_LN  * H_TGFb_T   * nT0_LN  - k_T0_pro * aT0_LN  - k_T0_death * aT0_LN
    daT0_LNl <- k_T0_act * f_act_LNl * H_TGFb_Ln1 * nT0_LNl - k_T0_pro * aT0_LNl - k_T0_death * aT0_LNl

    dT0_LN  <- k_T0_pro * aT0_LN  * N_aT0 - q_T0_LN_out * T0_LN  - k_T0_death * T0_LN
    dT0_LNl <- k_T0_pro * aT0_LNl * N_aT0 - q_T0_LN_out * T0_LNl - k_T0_death * T0_LNl
    dT0_LNo <- k_T0_act * f_act_LNl * H_TGFb_oth * nT0_LNo - q_T0_LN_out * T0_LNo - k_T0_death * T0_LNo

    # Helper T cell (Th) activation in LN
    daTh_LN  <- k_Th_act * f_act_LN  * nT0_LN  - k_T0_pro * aTh_LN  - k_T0_death * aTh_LN
    daTh_LNl <- k_Th_act * f_act_LNl * nT0_LNl - k_T0_pro * aTh_LNl - k_T0_death * aTh_LNl

    dTh_LN  <- k_T0_pro * aTh_LN  * N_aTh - q_T0_LN_out * Th_LN  - k_T0_death * Th_LN
    dTh_LNl <- k_T0_pro * aTh_LNl * N_aTh - q_T0_LN_out * Th_LNl - k_T0_death * Th_LNl
    dTh_LNo <- k_Th_act * f_act_LNl * nT0_LNo - q_T0_LN_out * Th_LNo - k_T0_death * Th_LNo

    # CD8 cytotoxic T cells per neoantigen j = 1..8
    daT_list_LN  <- vector("numeric", 8)
    daT_list_LNl <- vector("numeric", 8)
    dT_LN_list   <- vector("numeric", 8)
    dT_LNl_list  <- vector("numeric", 8)
    dT_LNo_list  <- vector("numeric", 8)
    dT_C_list    <- vector("numeric", 8)
    dT_P_list    <- vector("numeric", 8)

    for (j in 1:8) {
      nj   <- get(paste0("n_T",j,"_clones"))
      aT_j_LN  <- get(paste0("aT",j,"_LN"))
      aT_j_LNl <- get(paste0("aT",j,"_LNl"))
      T_j_LN   <- get(paste0("T",j,"_LN"))
      T_j_LNl  <- get(paste0("T",j,"_LNl"))
      T_j_LNo  <- get(paste0("T",j,"_LNo"))
      T_j_C    <- get(paste0("T",j,"_C"))
      T_j_P    <- get(paste0("T",j,"_P"))

      # Proliferation boost from IL-2
      N_div_j_LN  <- N0 + N_costim + N_IL2_CD8 * f_IL2_LN
      N_div_j_LNl <- N0 + N_costim + N_IL2_CD8 * f_IL2_LNl

      daT_list_LN[j]  <- k_T1_act * f_act_LN  * nj * nT1_LN  / div_T1 - k_T1_pro * aT_j_LN  - k_T1_death * aT_j_LN
      daT_list_LNl[j] <- k_T1_act * f_act_LNl * nj * nT1_LNl / div_T1 - k_T1_pro * aT_j_LNl - k_T1_death * aT_j_LNl

      dT_LN_list[j]   <- k_T1_pro * aT_j_LN  * N_div_j_LN  - q_T1_LN_out * T_j_LN  - k_T1_death * T_j_LN
      dT_LNl_list[j]  <- k_T1_pro * aT_j_LNl * N_div_j_LNl - q_T1_LN_out * T_j_LNl - k_T1_death * T_j_LNl
      dT_LNo_list[j]  <- k_T1_act * f_act_LNl * nj * nT1_LNo / div_T1 - q_T1_LN_out * T_j_LNo - k_T1_death * T_j_LNo

      # T cells exported from LN -> Central
      dT_C_list[j] <- q_T1_LN_out * T_j_LN + q_T1_LN_out * T_j_LNl + q_T1_LN_out * T_j_LNo \
                      - k_T1_death * T_j_C \
                      - q_T1_T_in * 1440 * d_T^3 * pi / 6 * T_j_C \
                      - q_T1_T_in_Ln1 * 1440 * d_Ln1^3 * pi / 6 * T_j_C \
                      - q_T1_T_in_Ln2 * 1440 * d_Ln2^3 * pi / 6 * T_j_C \
                      - q_T1_T_in_other * 1440 * d_oth^3 * pi / 6 * T_j_C

      dT_P_list[j] <- q_T1_P_in * 1440 * T_j_C - q_T1_P_out * T_j_P - k_T1_death * T_j_P
    }

    # =========================================================================
    # MODULE 10: T CELL DYNAMICS IN TUMOR (infiltration + exhaustion + Treg killing)
    # =========================================================================
    dT_tumor <- function(T_j_T, T_j_C, Tcyt_T, C_tot_T, T0_T, d_t, q_T1_T_in_loc) {
      infiltration   <- q_T1_T_in_loc * 1440 * d_t^3 * pi / 6 * T_j_C
      exhaustion     <- k_T1 * H_PD1_T * T_j_T
      treg_kill      <- k_Treg * T0_T / (T0_T + 1) * T_j_T
      death          <- k_T1_death * T_j_T
      infiltration - exhaustion - treg_kill - death
    }

    dT_oth_tumor <- function(T_j_T, T_j_C, H_PD1, H_TGFb_e, H_MDSC, T0, k_T1_T_in) {
      infiltration <- k_T1_T_in * 1440 * d_oth^3 * pi / 6 * T_j_C
      exhaustion   <- k_T1 * H_PD1 * T_j_T
      treg_kill    <- k_Treg * T0 / (T0 + 1) * T_j_T
      infiltration - exhaustion - treg_kill - k_T1_death * T_j_T
    }

    dT_Ln1_list <- vector("numeric", 8)
    dT_Ln2_list <- vector("numeric", 8)
    dT_oth_list <- vector("numeric", 8)
    dT_T_list   <- vector("numeric", 8)

    for (j in 1:8) {
      T_j_C    <- get(paste0("T",j,"_C"))
      T_j_T    <- get(paste0("T",j,"_T"))
      T_j_Ln1  <- get(paste0("T",j,"_Ln1"))
      T_j_Ln2  <- get(paste0("T",j,"_Ln2"))
      T_j_oth  <- get(paste0("T",j,"_oth"))

      dT_T_list[j]   <- dT_tumor(T_j_T,   T_j_C, Tcyt_total_T,   C_total_T,   T0_T,   d_T,   q_T1_T_in)
      dT_Ln1_list[j] <- dT_tumor(T_j_Ln1, T_j_C, Tcyt_total_Ln1, C_total_Ln1, T0_Ln1, d_Ln1, q_T1_T_in_Ln1)
      dT_Ln2_list[j] <- dT_tumor(T_j_Ln2, T_j_C, Tcyt_total_Ln2, C_total_Ln2, T0_Ln2, d_Ln2, q_T1_T_in_Ln2)
      dT_oth_list[j] <- dT_oth_tumor(T_j_oth, T_j_C, H_PD1_oth, H_TGFb_Teff_oth, H_MDSC_oth, T0_oth, q_T1_T_in_other)
    }

    # Treg (T0) and Th in tumors (infiltrate from central)
    dT0_T <- q_T0_P_in * 1440 * T0_C * V_T_vol / V_C \
             + k_Th_Treg * H_TGFb_T * Th_T \
             - k_T0_death * T0_T
    dT0_Ln1 <- q_T0_P_in * 1440 * T0_C * V_Ln1_vol / V_C + k_Th_Treg * H_TGFb_Ln1 * Th_Ln1 - k_T0_death * T0_Ln1
    dT0_Ln2 <- q_T0_P_in * 1440 * T0_C * V_Ln2_vol / V_C + k_Th_Treg * H_TGFb_Ln2 * Th_Ln2 - k_T0_death * T0_Ln2
    dT0_oth <- q_T0_P_in * 1440 * T0_C * V_oth_vol  / V_C + k_Th_Treg * H_TGFb_oth * Th_oth  - k_T0_death * T0_oth

    dTh_T <- q_T0_P_in * 1440 * Th_C * V_T_vol   / V_C - k_T0_death * Th_T   - k_Th_Treg * H_TGFb_T   * Th_T
    dTh_Ln1 <- q_T0_P_in * 1440 * Th_C * V_Ln1_vol / V_C - k_T0_death * Th_Ln1 - k_Th_Treg * H_TGFb_Ln1 * Th_Ln1
    dTh_Ln2 <- q_T0_P_in * 1440 * Th_C * V_Ln2_vol / V_C - k_T0_death * Th_Ln2 - k_Th_Treg * H_TGFb_Ln2 * Th_Ln2
    dTh_oth <- q_T0_P_in * 1440 * Th_C * V_oth_vol  / V_C - k_T0_death * Th_oth  - k_Th_Treg * H_TGFb_oth * Th_oth

    # Treg / Th in Central and Peripheral
    dT0_C  <- q_T0_LN_out * T0_LN + q_T0_LN_out * T0_LNl + q_T0_LN_out * T0_LNo \
              - k_T0_death * T0_C
    dTh_C  <- q_T0_LN_out * Th_LN + q_T0_LN_out * Th_LNl + q_T0_LN_out * Th_LNo \
              - k_T0_death * Th_C
    dT0_P  <- q_T0_P_in * 1440 * T0_C - q_T0_P_out * T0_P - k_T0_death * T0_P
    dTh_P  <- q_T0_P_in * 1440 * Th_C - q_T0_P_out * Th_P - k_T0_death * Th_P

    # =========================================================================
    # ASSEMBLE DERIVATIVES
    # =========================================================================
    derivs <- c(
      # Antibody
      daPD1_C = daPD1_C, daPD1_P = daPD1_P,
      daPD1_T = daPD1_T, daPD1_Ln1 = daPD1_Ln1, daPD1_Ln2 = daPD1_Ln2, daPD1_oth = daPD1_oth,
      daPD1_LN = daPD1_LN, daPD1_LNl = daPD1_LNl, daPD1_LNo = daPD1_LNo,
      # Carrying capacity
      dK_T = dK_T, dK_Ln1 = dK_Ln1, dK_Ln2 = dK_Ln2, dK_oth = dK_oth,
      # Angiogenic factor
      dcvas_T = dcvas_T, dcvas_Ln1 = dcvas_Ln1, dcvas_Ln2 = dcvas_Ln2, dcvas_oth = dcvas_oth,
      # Cancer clones
      dC1_T=dC1_T, dC2_T=dC2_T, dC3_T=dC3_T, dC4_T=dC4_T, dC5_T=dC5_T,
      dC1_Ln1=dC1_Ln1, dC2_Ln1=dC2_Ln1, dC3_Ln1=dC3_Ln1, dC4_Ln1=dC4_Ln1, dC5_Ln1=dC5_Ln1,
      dC1_Ln2=dC1_Ln2, dC2_Ln2=dC2_Ln2, dC3_Ln2=dC3_Ln2, dC4_Ln2=dC4_Ln2, dC5_Ln2=dC5_Ln2,
      dC1_oth=dC1_oth, dC2_oth=dC2_oth, dC3_oth=dC3_oth, dC4_oth=dC4_oth, dC5_oth=dC5_oth,
      # CCL2 / MDSC / NO / ArgI
      dCCL2_T=dCCL2_T, dCCL2_Ln1=dCCL2_Ln1, dCCL2_Ln2=dCCL2_Ln2, dCCL2_oth=dCCL2_oth,
      dMDSC_T=dMDSC_T, dMDSC_Ln1=dMDSC_Ln1, dMDSC_Ln2=dMDSC_Ln2, dMDSC_oth=dMDSC_oth,
      dNO_T=dNO_T, dNO_Ln1=dNO_Ln1, dNO_Ln2=dNO_Ln2, dNO_oth=dNO_oth,
      dArgI_T=dArgI_T, dArgI_Ln1=dArgI_Ln1, dArgI_Ln2=dArgI_Ln2, dArgI_oth=dArgI_oth,
      # Macrophages
      dMac_M1_T=dMac_M1_T, dMac_M2_T=dMac_M2_T,
      dMac_M1_Ln1=dMac_M1_Ln1, dMac_M2_Ln1=dMac_M2_Ln1,
      dMac_M1_Ln2=dMac_M1_Ln2, dMac_M2_Ln2=dMac_M2_Ln2,
      dMac_M1_oth=dMac_M1_oth, dMac_M2_oth=dMac_M2_oth,
      # Cytokines
      dTGFb_T=dTGFb_T, dTGFb_Ln1=dTGFb_Ln1, dTGFb_Ln2=dTGFb_Ln2, dTGFb_oth=dTGFb_oth,
      dIFNg_T=dIFNg_T, dIFNg_Ln1=dIFNg_Ln1, dIFNg_Ln2=dIFNg_Ln2, dIFNg_oth=dIFNg_oth,
      dIL12_T=dIL12_T, dIL12_Ln1=dIL12_Ln1, dIL12_Ln2=dIL12_Ln2, dIL12_oth=dIL12_oth,
      dIL10_T=dIL10_T, dIL10_Ln1=dIL10_Ln1, dIL10_Ln2=dIL10_Ln2, dIL10_oth=dIL10_oth,
      # APC
      dc_T=dc_T, dc_Ln1=dc_Ln1, dc_Ln2=dc_Ln2, dc_oth=dc_oth,
      dAPC_T=dAPC_T, dAPC_Ln1=dAPC_Ln1, dAPC_Ln2=dAPC_Ln2, dAPC_oth=dAPC_oth,
      dmAPC_T=dmAPC_T, dmAPC_Ln1=dmAPC_Ln1, dmAPC_Ln2=dmAPC_Ln2, dmAPC_oth=dmAPC_oth,
      dmAPC_LN=dmAPC_LN, dmAPC_LNl=dmAPC_LNl, dmAPC_LNo=dmAPC_LNo,
      # Naive T cells
      dnT0_C=dnT0_C, dnT0_P=dnT0_P, dnT0_LN=dnT0_LN, dnT0_LNl=dnT0_LNl, dnT0_LNo=dnT0_LNo,
      dnT1_C=dnT1_C, dnT1_P=dnT1_P, dnT1_LN=dnT1_LN, dnT1_LNl=dnT1_LNl, dnT1_LNo=dnT1_LNo,
      # Activated T in LN
      daT0_LN=daT0_LN, daT0_LNl=daT0_LNl,
      daTh_LN=daTh_LN, daTh_LNl=daTh_LNl,
      dT0_LN=dT0_LN, dT0_LNl=dT0_LNl, dT0_LNo=dT0_LNo,
      dTh_LN=dTh_LN, dTh_LNl=dTh_LNl, dTh_LNo=dTh_LNo,
      dIL2_LN=dIL2_LN, dIL2_LNl=dIL2_LNl,
      # Treg/Th in central/peripheral/tumor
      dT0_C=dT0_C, dT0_P=dT0_P,
      dTh_C=dTh_C, dTh_P=dTh_P,
      dT0_T=dT0_T, dT0_Ln1=dT0_Ln1, dT0_Ln2=dT0_Ln2, dT0_oth=dT0_oth,
      dTh_T=dTh_T, dTh_Ln1=dTh_Ln1, dTh_Ln2=dTh_Ln2, dTh_oth=dTh_oth
    )

    # CD8 T cells per neo-antigen
    for (j in 1:8) {
      derivs[paste0("daT",j,"_LN")]  <- daT_list_LN[j]
      derivs[paste0("daT",j,"_LNl")] <- daT_list_LNl[j]
      derivs[paste0("dT",j,"_LN")]   <- dT_LN_list[j]
      derivs[paste0("dT",j,"_LNl")]  <- dT_LNl_list[j]
      derivs[paste0("dT",j,"_LNo")]  <- dT_LNo_list[j]
      derivs[paste0("dT",j,"_C")]    <- dT_C_list[j]
      derivs[paste0("dT",j,"_P")]    <- dT_P_list[j]
      derivs[paste0("dT",j,"_T")]    <- dT_T_list[j]
      derivs[paste0("dT",j,"_Ln1")]  <- dT_Ln1_list[j]
      derivs[paste0("dT",j,"_Ln2")]  <- dT_Ln2_list[j]
      derivs[paste0("dT",j,"_oth")]  <- dT_oth_list[j]
    }

    list(derivs)
  })
}
