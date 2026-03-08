# =============================================================================
# 02_model_ode.R
# TNBC QSP Model - ODE System (optimized: all loops unrolled, no get/sapply)
# Based on: Arulraj et al., Sci. Adv. 9, eadg0289 (2023)
# =============================================================================

tnbc_ode <- function(t, state, parms) {
  with(as.list(c(state, parms)), {

    # =========================================================================
    # ALGEBRAIC RULES
    # =========================================================================
    C_total_T   <- C1_T   + C2_T   + C3_T   + C4_T   + C5_T
    C_total_Ln1 <- C1_Ln1 + C2_Ln1 + C3_Ln1 + C4_Ln1 + C5_Ln1
    C_total_Ln2 <- C1_Ln2 + C2_Ln2 + C3_Ln2 + C4_Ln2 + C5_Ln2
    C_total_oth <- C1_oth + C2_oth + C3_oth + C4_oth + C5_oth

    # CD8 T cells per compartment (unrolled j=1..8)
    Tcyt_total_T   <- T1_T+T2_T+T3_T+T4_T+T5_T+T6_T+T7_T+T8_T
    Tcyt_total_Ln1 <- T1_Ln1+T2_Ln1+T3_Ln1+T4_Ln1+T5_Ln1+T6_Ln1+T7_Ln1+T8_Ln1
    Tcyt_total_Ln2 <- T1_Ln2+T2_Ln2+T3_Ln2+T4_Ln2+T5_Ln2+T6_Ln2+T7_Ln2+T8_Ln2
    Tcyt_total_oth <- T1_oth+T2_oth+T3_oth+T4_oth+T5_oth+T6_oth+T7_oth+T8_oth

    M_total_T   <- Mac_M1_T   + Mac_M2_T
    M_total_Ln1 <- Mac_M1_Ln1 + Mac_M2_Ln1
    M_total_Ln2 <- Mac_M1_Ln2 + Mac_M2_Ln2
    M_total_oth <- Mac_M1_oth + Mac_M2_oth

    # Tumor volume (cancer cells only - no positive feedback via T/M cells)
    V_T_vol   <- max(V_Tmin, C_total_T   * vol_cell / Ve_T * 1e-12)
    V_Ln1_vol <- max(V_Tmin, C_total_Ln1 * vol_cell / Ve_T * 1e-12)
    V_Ln2_vol <- max(V_Tmin, C_total_Ln2 * vol_cell / Ve_T * 1e-12)
    V_oth_vol <- max(V_Tmin, C_total_oth  * vol_cell / Ve_T * 1e-12)

    # Tumor diameters (cm)
    d_T   <- (6 * V_T_vol   / pi)^(1/3)
    d_Ln1 <- (6 * V_Ln1_vol / pi)^(1/3)
    d_Ln2 <- (6 * V_Ln2_vol / pi)^(1/3)
    d_oth <- (6 * V_oth_vol  / pi)^(1/3)

    # TGFb Hill
    H_TGFb_T   <- TGFb_T   / (TGFb_50 + TGFb_T)
    H_TGFb_Ln1 <- TGFb_Ln1 / (TGFb_50 + TGFb_Ln1)
    H_TGFb_Ln2 <- TGFb_Ln2 / (TGFb_50 + TGFb_Ln2)
    H_TGFb_oth <- TGFb_oth / (TGFb_50 + TGFb_oth)

    H_TGFb_Teff_T   <- TGFb_T   / (TGFb_50_Teff + TGFb_T)
    H_TGFb_Teff_Ln1 <- TGFb_Ln1 / (TGFb_50_Teff + TGFb_Ln1)
    H_TGFb_Teff_Ln2 <- TGFb_Ln2 / (TGFb_50_Teff + TGFb_Ln2)
    H_TGFb_Teff_oth <- TGFb_oth / (TGFb_50_Teff + TGFb_oth)

    # MDSC suppression
    H_MDSC_T   <- (ArgI_T   / ArgI_50_Teff + NO_T   / NO_50_Teff) / (1 + ArgI_T   / ArgI_50_Teff + NO_T   / NO_50_Teff)
    H_MDSC_Ln1 <- (ArgI_Ln1 / ArgI_50_Teff + NO_Ln1 / NO_50_Teff) / (1 + ArgI_Ln1 / ArgI_50_Teff + NO_Ln1 / NO_50_Teff)
    H_MDSC_Ln2 <- (ArgI_Ln2 / ArgI_50_Teff + NO_Ln2 / NO_50_Teff) / (1 + ArgI_Ln2 / ArgI_50_Teff + NO_Ln2 / NO_50_Teff)
    H_MDSC_oth <- (ArgI_oth / ArgI_50_Teff + NO_oth / NO_50_Teff) / (1 + ArgI_oth / ArgI_50_Teff + NO_oth / NO_50_Teff)

    # PD-L1 expression
    PDL1_T   <- C1_PDL1_base * (1 + r_PDL1_IFNg * IFNg_T   / (IFNg_50_ind * 1e-3 + IFNg_T))
    PDL1_Ln1 <- C1_PDL1_base * (1 + r_PDL1_IFNg * IFNg_Ln1 / (IFNg_50_ind * 1e-3 + IFNg_Ln1))
    PDL1_Ln2 <- C1_PDL1_base * (1 + r_PDL1_IFNg * IFNg_Ln2 / (IFNg_50_ind * 1e-3 + IFNg_Ln2))
    PDL1_oth <- C1_PDL1_base * (1 + r_PDL1_IFNg * IFNg_oth / (IFNg_50_ind * 1e-3 + IFNg_oth))

    # PD-1 blockade
    Kd_aPD1          <- koff_PD1_aPD1 / kon_PD1_aPD1 * 1e9
    f_PD1_block_T   <- aPD1_T   / (Kd_aPD1 + aPD1_T)
    f_PD1_block_Ln1 <- aPD1_Ln1 / (Kd_aPD1 + aPD1_Ln1)
    f_PD1_block_Ln2 <- aPD1_Ln2 / (Kd_aPD1 + aPD1_Ln2)
    f_PD1_block_oth <- aPD1_oth / (Kd_aPD1 + aPD1_oth)

    PDL1_free_T   <- PDL1_T   * (1 - f_PD1_block_T)   / A_cell
    PDL1_free_Ln1 <- PDL1_Ln1 * (1 - f_PD1_block_Ln1) / A_cell
    PDL1_free_Ln2 <- PDL1_Ln2 * (1 - f_PD1_block_Ln2) / A_cell
    PDL1_free_oth <- PDL1_oth * (1 - f_PD1_block_oth)  / A_cell

    H_PD1_T   <- PDL1_free_T^n_PD1   / (PD1_50^n_PD1 + PDL1_free_T^n_PD1)
    H_PD1_Ln1 <- PDL1_free_Ln1^n_PD1 / (PD1_50^n_PD1 + PDL1_free_Ln1^n_PD1)
    H_PD1_Ln2 <- PDL1_free_Ln2^n_PD1 / (PD1_50^n_PD1 + PDL1_free_Ln2^n_PD1)
    H_PD1_oth <- PDL1_free_oth^n_PD1  / (PD1_50^n_PD1 + PDL1_free_oth^n_PD1)

    # Effective kill rate
    k_kill_T   <- k_C_T1 * (1 - H_PD1_T)   * (1 - H_TGFb_Teff_T)   * (1 - H_MDSC_T)
    k_kill_Ln1 <- k_C_T1 * (1 - H_PD1_Ln1) * (1 - H_TGFb_Teff_Ln1) * (1 - H_MDSC_Ln1)
    k_kill_Ln2 <- k_C_T1 * (1 - H_PD1_Ln2) * (1 - H_TGFb_Teff_Ln2) * (1 - H_MDSC_Ln2)
    k_kill_oth <- k_C_T1 * (1 - H_PD1_oth)  * (1 - H_TGFb_Teff_oth)  * (1 - H_MDSC_oth)

    # IL-10 Hill
    H_IL10_T   <- IL10_T^2   / (IL10_50^2 + IL10_T^2)
    H_IL10_Ln1 <- IL10_Ln1^2 / (IL10_50^2 + IL10_Ln1^2)
    H_IL10_Ln2 <- IL10_Ln2^2 / (IL10_50^2 + IL10_Ln2^2)
    H_IL10_oth <- IL10_oth^2  / (IL10_50^2 + IL10_oth^2)

    H_IL12_T   <- IL12_T^2   / (IL12_50^2 + IL12_T^2)
    H_IL12_Ln1 <- IL12_Ln1^2 / (IL12_50^2 + IL12_Ln1^2)
    H_IL12_Ln2 <- IL12_Ln2^2 / (IL12_50^2 + IL12_Ln2^2)
    H_IL12_oth <- IL12_oth^2  / (IL12_50^2 + IL12_oth^2)

    H_IFNg_T   <- IFNg_T^2   / (IFNg_50^2 + IFNg_T^2)
    H_IFNg_Ln1 <- IFNg_Ln1^2 / (IFNg_50^2 + IFNg_Ln1^2)
    H_IFNg_Ln2 <- IFNg_Ln2^2 / (IFNg_50^2 + IFNg_Ln2^2)
    H_IFNg_oth <- IFNg_oth^2  / (IFNg_50^2 + IFNg_oth^2)

    H_IL10_phago_T   <- IL10_T   / (IL10_50_phago + IL10_T)
    H_IL10_phago_Ln1 <- IL10_Ln1 / (IL10_50_phago + IL10_Ln1)
    H_IL10_phago_Ln2 <- IL10_Ln2 / (IL10_50_phago + IL10_Ln2)
    H_IL10_phago_oth <- IL10_oth  / (IL10_50_phago + IL10_oth)

    SIRPa_syn_T <- M_SIRPa * A_syn / A_Mcell
    H_SIRPa_T   <- SIRPa_syn_T^n_SIRPa / (SIRPa_50^n_SIRPa + SIRPa_syn_T^n_SIRPa)
    H_SIRPa_Ln1 <- H_SIRPa_T
    H_SIRPa_Ln2 <- H_SIRPa_T
    H_SIRPa_oth <- H_SIRPa_T

    k_phago_T   <- k_M1_phago * (1 - H_IL10_phago_T)   * (1 - H_SIRPa_T)   * (1 - H_PD1_T)
    k_phago_Ln1 <- k_M1_phago * (1 - H_IL10_phago_Ln1) * (1 - H_SIRPa_Ln1) * (1 - H_PD1_Ln1)
    k_phago_Ln2 <- k_M1_phago * (1 - H_IL10_phago_Ln2) * (1 - H_SIRPa_Ln2) * (1 - H_PD1_Ln2)
    k_phago_oth <- k_M1_phago * (1 - H_IL10_phago_oth)  * (1 - H_SIRPa_oth) * (1 - H_PD1_oth)

    # =========================================================================
    # MODULE 1: PEMBROLIZUMAB PK
    # =========================================================================
    CL_aPD1 <- k_cl_aPD1 / V_C

    daPD1_C <- (- CL_aPD1 * aPD1_C
                - q_P_aPD1  * 86400 / (V_C * 1000)      * (aPD1_C - aPD1_P)
                - q_T_aPD1  * 86400 / (V_C * 1e6)       * (aPD1_C - aPD1_T)   * gamma_T_aPD1
                - q_LN_aPD1 * 86400 / (V_C * 1e6)       * (aPD1_C - aPD1_LN)  * gamma_LN_aPD1)

    daPD1_P   <- q_P_aPD1  * 86400 / (V_P * 1000)      * (aPD1_C - aPD1_P)
    daPD1_T   <- q_T_aPD1  * 86400 / (V_T_vol   * 1e6) * (aPD1_C - aPD1_T)   * gamma_T_aPD1
    daPD1_Ln1 <- q_T_aPD1  * 86400 / (V_Ln1_vol * 1e6) * (aPD1_C - aPD1_Ln1) * gamma_T_aPD1
    daPD1_Ln2 <- q_T_aPD1  * 86400 / (V_Ln2_vol * 1e6) * (aPD1_C - aPD1_Ln2) * gamma_T_aPD1
    daPD1_oth <- q_T_aPD1  * 86400 / (V_oth_vol  * 1e6) * (aPD1_C - aPD1_oth) * gamma_T_aPD1

    daPD1_LN  <- (q_LN_aPD1 * 86400 / (V_LN * 1e-6) * (aPD1_C - aPD1_LN)  * gamma_LN_aPD1
                  - q_LD_aPD1 * 1440 * aPD1_LN)
    daPD1_LNl <- (q_LN_aPD1 * 86400 / (V_LN * 1e-6) * (aPD1_C - aPD1_LNl) * gamma_LN_aPD1
                  - q_LD_aPD1 * 1440 * aPD1_LNl)
    daPD1_LNo <- (q_LN_aPD1 * 86400 / (V_LN * 1e-6) * (aPD1_C - aPD1_LNo) * gamma_LN_aPD1
                  - q_LD_aPD1 * 1440 * aPD1_LNo)

    # =========================================================================
    # MODULE 2: ANGIOGENESIS / CARRYING CAPACITY
    # K_T tracks tumor vasculature capacity.
    # Use relaxation toward target K = C_total * K_ratio (angiogenesis-driven).
    # This prevents unbounded exponential growth while maintaining biology.
    # =========================================================================
    dcvas_T   <- k_vas_Csec * C_total_T   + k_vas_Msec * Mac_M2_T   - k_vas_deg * cvas_T
    dcvas_Ln1 <- k_vas_Csec * C_total_Ln1 + k_vas_Msec * Mac_M2_Ln1 - k_vas_deg * cvas_Ln1
    dcvas_Ln2 <- k_vas_Csec * C_total_Ln2 + k_vas_Msec * Mac_M2_Ln2 - k_vas_deg * cvas_Ln2
    dcvas_oth <- k_vas_Csec * C_total_oth  + k_vas_Msec * Mac_M2_oth  - k_vas_deg * cvas_oth

    f_vas_T   <- cvas_T   / (c_vas_50 + cvas_T)
    f_vas_Ln1 <- cvas_Ln1 / (c_vas_50 + cvas_Ln1)
    f_vas_Ln2 <- cvas_Ln2 / (c_vas_50 + cvas_Ln2)
    f_vas_oth <- cvas_oth / (c_vas_50 + cvas_oth)

    # K_target: angiogenesis-supported carrying capacity
    # At SS: K_T_ss proportional to C_total via f_vas
    # Use simple model: K_target = C_total * (1 + f_vas * k_K_g/k_K_d)
    K_target_T   <- max(C_total_T   * (1 + f_vas_T   * k_K_g / k_K_d), 1e4)
    K_target_Ln1 <- max(C_total_Ln1 * (1 + f_vas_Ln1 * k_K_g / k_K_d), 1e3)
    K_target_Ln2 <- max(C_total_Ln2 * (1 + f_vas_Ln2 * k_K_g / k_K_d), 1e3)
    K_target_oth <- max(C_total_oth  * (1 + f_vas_oth  * k_K_g / k_K_d), 1e3)

    # Relax K toward K_target with time scale 1/k_K_d (slow, biological)
    dK_T   <- k_K_d * (K_target_T   - K_T)
    dK_Ln1 <- k_K_d * (K_target_Ln1 - K_Ln1)
    dK_Ln2 <- k_K_d * (K_target_Ln2 - K_Ln2)
    dK_oth <- k_K_d * (K_target_oth  - K_oth)

    # =========================================================================
    # MODULE 3: CANCER CELL DYNAMICS (Modified Gompertz, inline)
    # =========================================================================
    .grow_T  <- log(max(K_T,   1) / max(C_total_T   + 1, 1)) * (C_total_T   > 0)
    .grow_L1 <- log(max(K_Ln1, 1) / max(C_total_Ln1 + 1, 1)) * (C_total_Ln1 > 0)
    .grow_L2 <- log(max(K_Ln2, 1) / max(C_total_Ln2 + 1, 1)) * (C_total_Ln2 > 0)
    .grow_oth <- log(max(K_oth, 1) / max(C_total_oth + 1, 1)) * (C_total_oth > 0)

    .kill_T_den   <- K_T_C * C_total_T   + Tcyt_total_T   + 1e-10
    .kill_L1_den  <- K_T_C * C_total_Ln1 + Tcyt_total_Ln1 + 1e-10
    .kill_L2_den  <- K_T_C * C_total_Ln2 + Tcyt_total_Ln2 + 1e-10
    .kill_oth_den <- K_T_C * C_total_oth  + Tcyt_total_oth  + 1e-10

    .phago_T_den   <- K_Mac_C * C_total_T   + Mac_M1_T   + 1e-10
    .phago_L1_den  <- K_Mac_C * C_total_Ln1 + Mac_M1_Ln1 + 1e-10
    .phago_L2_den  <- K_Mac_C * C_total_Ln2 + Mac_M1_Ln2 + 1e-10
    .phago_oth_den <- K_Mac_C * C_total_oth  + Mac_M1_oth  + 1e-10

    .f_kill_T   <- k_kill_T   * Tcyt_total_T   / .kill_T_den
    .f_kill_L1  <- k_kill_Ln1 * Tcyt_total_Ln1 / .kill_L1_den
    .f_kill_L2  <- k_kill_Ln2 * Tcyt_total_Ln2 / .kill_L2_den
    .f_kill_oth <- k_kill_oth  * Tcyt_total_oth  / .kill_oth_den

    .f_phago_T   <- k_phago_T   * Mac_M1_T   / .phago_T_den
    .f_phago_L1  <- k_phago_Ln1 * Mac_M1_Ln1 / .phago_L1_den
    .f_phago_L2  <- k_phago_Ln2 * Mac_M1_Ln2 / .phago_L2_den
    .f_phago_oth <- k_phago_oth  * Mac_M1_oth  / .phago_oth_den

    dC1_T   <- start       * (k_C1_growth       * C1_T   * .grow_T   - k_C1_death       * C1_T   - .f_kill_T   * C1_T   - .f_phago_T   * C1_T)
    dC2_T   <- start       * (k_C2_growth       * C2_T   * .grow_T   - k_C2_death       * C2_T   - .f_kill_T   * C2_T   - .f_phago_T   * C2_T)
    dC3_T   <- start       * (k_C3_growth       * C3_T   * .grow_T   - k_C3_death       * C3_T   - .f_kill_T   * C3_T   - .f_phago_T   * C3_T)
    dC4_T   <- start       * (k_C4_growth       * C4_T   * .grow_T   - k_C4_death       * C4_T   - .f_kill_T   * C4_T   - .f_phago_T   * C4_T)
    dC5_T   <- start       * (k_C5_growth       * C5_T   * .grow_T   - k_C5_death       * C5_T   - .f_kill_T   * C5_T   - .f_phago_T   * C5_T)

    dC1_Ln1 <- start_Ln1   * (k_C1_growth_Ln1   * C1_Ln1 * .grow_L1  - k_C1_death_Ln1   * C1_Ln1 - .f_kill_L1  * C1_Ln1 - .f_phago_L1  * C1_Ln1)
    dC2_Ln1 <- start_Ln1   * (k_C2_growth_Ln1   * C2_Ln1 * .grow_L1  - k_C2_death_Ln1   * C2_Ln1 - .f_kill_L1  * C2_Ln1 - .f_phago_L1  * C2_Ln1)
    dC3_Ln1 <- start_Ln1   * (k_C3_growth_Ln1   * C3_Ln1 * .grow_L1  - k_C3_death_Ln1   * C3_Ln1 - .f_kill_L1  * C3_Ln1 - .f_phago_L1  * C3_Ln1)
    dC4_Ln1 <- start_Ln1   * (k_C4_growth_Ln1   * C4_Ln1 * .grow_L1  - k_C4_death_Ln1   * C4_Ln1 - .f_kill_L1  * C4_Ln1 - .f_phago_L1  * C4_Ln1)
    dC5_Ln1 <- start_Ln1   * (k_C5_growth_Ln1   * C5_Ln1 * .grow_L1  - k_C5_death_Ln1   * C5_Ln1 - .f_kill_L1  * C5_Ln1 - .f_phago_L1  * C5_Ln1)

    dC1_Ln2 <- start_Ln2   * (k_C1_growth_Ln2   * C1_Ln2 * .grow_L2  - k_C1_death_Ln2   * C1_Ln2 - .f_kill_L2  * C1_Ln2 - .f_phago_L2  * C1_Ln2)
    dC2_Ln2 <- start_Ln2   * (k_C2_growth_Ln2   * C2_Ln2 * .grow_L2  - k_C2_death_Ln2   * C2_Ln2 - .f_kill_L2  * C2_Ln2 - .f_phago_L2  * C2_Ln2)
    dC3_Ln2 <- start_Ln2   * (k_C3_growth_Ln2   * C3_Ln2 * .grow_L2  - k_C3_death_Ln2   * C3_Ln2 - .f_kill_L2  * C3_Ln2 - .f_phago_L2  * C3_Ln2)
    dC4_Ln2 <- start_Ln2   * (k_C4_growth_Ln2   * C4_Ln2 * .grow_L2  - k_C4_death_Ln2   * C4_Ln2 - .f_kill_L2  * C4_Ln2 - .f_phago_L2  * C4_Ln2)
    dC5_Ln2 <- start_Ln2   * (k_C5_growth_Ln2   * C5_Ln2 * .grow_L2  - k_C5_death_Ln2   * C5_Ln2 - .f_kill_L2  * C5_Ln2 - .f_phago_L2  * C5_Ln2)

    dC1_oth <- start_other * (k_C1_growth_other * C1_oth * .grow_oth - k_C1_death_other * C1_oth - .f_kill_oth * C1_oth - .f_phago_oth * C1_oth)
    dC2_oth <- start_other * (k_C2_growth_other * C2_oth * .grow_oth - k_C2_death_other * C2_oth - .f_kill_oth * C2_oth - .f_phago_oth * C2_oth)
    dC3_oth <- start_other * (k_C3_growth_other * C3_oth * .grow_oth - k_C3_death_other * C3_oth - .f_kill_oth * C3_oth - .f_phago_oth * C3_oth)
    dC4_oth <- start_other * (k_C4_growth_other * C4_oth * .grow_oth - k_C4_death_other * C4_oth - .f_kill_oth * C4_oth - .f_phago_oth * C4_oth)
    dC5_oth <- start_other * (k_C5_growth_other * C5_oth * .grow_oth - k_C5_death_other * C5_oth - .f_kill_oth * C5_oth - .f_phago_oth * C5_oth)

    # =========================================================================
    # MODULE 4: CCL2 / MDSC / NO / ArgI
    # =========================================================================
    k_CCL2_deg_day <- k_CCL2_deg * 1440

    dCCL2_T   <- k_CCL2_sec * C_total_T   - k_CCL2_deg_day * CCL2_T
    dCCL2_Ln1 <- k_CCL2_sec * C_total_Ln1 - k_CCL2_deg_day * CCL2_Ln1
    dCCL2_Ln2 <- k_CCL2_sec * C_total_Ln2 - k_CCL2_deg_day * CCL2_Ln2
    dCCL2_oth <- k_CCL2_sec * C_total_oth  - k_CCL2_deg_day * CCL2_oth

    f_CCL2_T   <- CCL2_T   / (CCL2_50 + CCL2_T)
    f_CCL2_Ln1 <- CCL2_Ln1 / (CCL2_50 + CCL2_Ln1)
    f_CCL2_Ln2 <- CCL2_Ln2 / (CCL2_50 + CCL2_Ln2)
    f_CCL2_oth <- CCL2_oth  / (CCL2_50 + CCL2_oth)

    dMDSC_T   <- k_MDSC_mig       * f_CCL2_T   - k_MDSC_death * MDSC_T
    dMDSC_Ln1 <- k_MDSC_mig_Ln1  * f_CCL2_Ln1 - k_MDSC_death * MDSC_Ln1
    dMDSC_Ln2 <- k_MDSC_mig_Ln2  * f_CCL2_Ln2 - k_MDSC_death * MDSC_Ln2
    dMDSC_oth <- k_MDSC_mig_other * f_CCL2_oth  - k_MDSC_death * MDSC_oth

    dNO_T   <- k_NO_sec * MDSC_T   - k_NO_deg * NO_T
    dNO_Ln1 <- k_NO_sec * MDSC_Ln1 - k_NO_deg * NO_Ln1
    dNO_Ln2 <- k_NO_sec * MDSC_Ln2 - k_NO_deg * NO_Ln2
    dNO_oth <- k_NO_sec * MDSC_oth  - k_NO_deg * NO_oth

    dArgI_T   <- k_ArgI_sec * MDSC_T   - k_ArgI_deg * ArgI_T
    dArgI_Ln1 <- k_ArgI_sec * MDSC_Ln1 - k_ArgI_deg * ArgI_Ln1
    dArgI_Ln2 <- k_ArgI_sec * MDSC_Ln2 - k_ArgI_deg * ArgI_Ln2
    dArgI_oth <- k_ArgI_sec * MDSC_oth  - k_ArgI_deg * ArgI_oth

    # =========================================================================
    # MODULE 5: MACROPHAGE DYNAMICS
    # =========================================================================
    dMac_M1_T <- (k_Mac_mig          * f_CCL2_T
                  + k_M1_pol         * (H_IL12_T   + H_IFNg_T)   * Mac_M2_T
                  - k_M2_pol         * H_IL10_T   * Mac_M1_T
                  - k_Mac_death      * Mac_M1_T)

    dMac_M2_T <- (k_M2_pol           * H_IL10_T   * Mac_M1_T
                  - k_M1_pol         * (H_IL12_T   + H_IFNg_T)   * Mac_M2_T
                  - k_Mac_death      * Mac_M2_T)

    dMac_M1_Ln1 <- (k_Mac_mig_Ln1   * f_CCL2_Ln1
                    + k_M1_pol       * (H_IL12_Ln1 + H_IFNg_Ln1) * Mac_M2_Ln1
                    - k_M2_pol_Ln1   * H_IL10_Ln1 * Mac_M1_Ln1
                    - k_Mac_death    * Mac_M1_Ln1)

    dMac_M2_Ln1 <- (k_M2_pol_Ln1    * H_IL10_Ln1 * Mac_M1_Ln1
                    - k_M1_pol       * (H_IL12_Ln1 + H_IFNg_Ln1) * Mac_M2_Ln1
                    - k_Mac_death    * Mac_M2_Ln1)

    dMac_M1_Ln2 <- (k_Mac_mig_Ln2   * f_CCL2_Ln2
                    + k_M1_pol       * (H_IL12_Ln2 + H_IFNg_Ln2) * Mac_M2_Ln2
                    - k_M2_pol_Ln2   * H_IL10_Ln2 * Mac_M1_Ln2
                    - k_Mac_death    * Mac_M1_Ln2)

    dMac_M2_Ln2 <- (k_M2_pol_Ln2    * H_IL10_Ln2 * Mac_M1_Ln2
                    - k_M1_pol       * (H_IL12_Ln2 + H_IFNg_Ln2) * Mac_M2_Ln2
                    - k_Mac_death    * Mac_M2_Ln2)

    dMac_M1_oth <- (k_Mac_mig_other  * f_CCL2_oth
                    + k_M1_pol       * (H_IL12_oth + H_IFNg_oth)  * Mac_M2_oth
                    - k_M2_pol_other * H_IL10_oth * Mac_M1_oth
                    - k_Mac_death    * Mac_M1_oth)

    dMac_M2_oth <- (k_M2_pol_other   * H_IL10_oth * Mac_M1_oth
                    - k_M1_pol       * (H_IL12_oth + H_IFNg_oth)  * Mac_M2_oth
                    - k_Mac_death    * Mac_M2_oth)

    # =========================================================================
    # MODULE 6: CYTOKINES
    # =========================================================================
    dTGFb_T   <- k_TGFb_Tsec * T0_T   + k_TGFb_Msec * Mac_M2_T   - k_TGFb_deg * TGFb_T   + TGFbase * k_TGFb_deg
    dTGFb_Ln1 <- k_TGFb_Tsec * T0_Ln1 + k_TGFb_Msec * Mac_M2_Ln1 - k_TGFb_deg * TGFb_Ln1 + TGFbase * k_TGFb_deg
    dTGFb_Ln2 <- k_TGFb_Tsec * T0_Ln2 + k_TGFb_Msec * Mac_M2_Ln2 - k_TGFb_deg * TGFb_Ln2 + TGFbase * k_TGFb_deg
    dTGFb_oth <- k_TGFb_Tsec * T0_oth  + k_TGFb_Msec * Mac_M2_oth  - k_TGFb_deg * TGFb_oth  + TGFbase * k_TGFb_deg

    dIFNg_T   <- k_IFNg_sec * Th_T   - k_IFNg_deg * IFNg_T
    dIFNg_Ln1 <- k_IFNg_sec * Th_Ln1 - k_IFNg_deg * IFNg_Ln1
    dIFNg_Ln2 <- k_IFNg_sec * Th_Ln2 - k_IFNg_deg * IFNg_Ln2
    dIFNg_oth <- k_IFNg_sec * Th_oth  - k_IFNg_deg * IFNg_oth

    k_IL12_deg_day <- k_IL12_deg * 1440
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
    death_rate_T   <- k_C1_death*C1_T   + k_C2_death*C2_T   + k_C3_death*C3_T   + k_C4_death*C4_T   + k_C5_death*C5_T
    death_rate_Ln1 <- k_C1_death_Ln1*C1_Ln1 + k_C2_death_Ln1*C2_Ln1 + k_C3_death_Ln1*C3_Ln1 + k_C4_death_Ln1*C4_Ln1 + k_C5_death_Ln1*C5_Ln1
    death_rate_Ln2 <- k_C1_death_Ln2*C1_Ln2 + k_C2_death_Ln2*C2_Ln2 + k_C3_death_Ln2*C3_Ln2 + k_C4_death_Ln2*C4_Ln2 + k_C5_death_Ln2*C5_Ln2
    death_rate_oth <- k_C1_death_other*C1_oth + k_C2_death_other*C2_oth + k_C3_death_other*C3_oth + k_C4_death_other*C4_oth + k_C5_death_other*C5_oth

    dc_T   <- DAMPs * death_rate_T   / V_T_vol   - k_c * c_T
    dc_Ln1 <- DAMPs * death_rate_Ln1 / V_Ln1_vol - k_c * c_Ln1
    dc_Ln2 <- DAMPs * death_rate_Ln2 / V_Ln2_vol - k_c * c_Ln2
    dc_oth <- DAMPs * death_rate_oth  / V_oth_vol  - k_c * c_oth

    f_c_T   <- c_T   / (c50 + c_T)
    f_c_Ln1 <- c_Ln1 / (c50 + c_Ln1)
    f_c_Ln2 <- c_Ln2 / (c50 + c_Ln2)
    f_c_oth <- c_oth / (c50 + c_oth)

    f_mat_T   <- k_APC_mat * f_c_T   * (1 - H_IL10_T)
    f_mat_Ln1 <- k_APC_mat * f_c_Ln1 * (1 - H_IL10_Ln1)
    f_mat_Ln2 <- k_APC_mat * f_c_Ln2 * (1 - H_IL10_Ln2)
    f_mat_oth <- k_APC_mat * f_c_oth  * (1 - H_IL10_oth)

    dAPC_T   <- APC0_T       * k_APC_death - f_mat_T   * APC_T   - k_APC_death * APC_T
    dAPC_Ln1 <- APC0_T_Ln1   * k_APC_death - f_mat_Ln1 * APC_Ln1 - k_APC_death * APC_Ln1
    dAPC_Ln2 <- APC0_T_Ln2   * k_APC_death - f_mat_Ln2 * APC_Ln2 - k_APC_death * APC_Ln2
    dAPC_oth <- APC0_T_other  * k_APC_death - f_mat_oth * APC_oth  - k_APC_death * APC_oth

    dmAPC_T   <- f_mat_T   * APC_T   - k_APC_mig * mAPC_T   - k_mAPC_death * mAPC_T
    dmAPC_Ln1 <- f_mat_Ln1 * APC_Ln1 - k_APC_mig * mAPC_Ln1 - k_mAPC_death * mAPC_Ln1
    dmAPC_Ln2 <- f_mat_Ln2 * APC_Ln2 - k_APC_mig * mAPC_Ln2 - k_mAPC_death * mAPC_Ln2
    dmAPC_oth <- f_mat_oth  * APC_oth  - k_APC_mig * mAPC_oth  - k_mAPC_death * mAPC_oth

    dmAPC_LN  <- k_APC_mig * mAPC_T   * V_T_vol   / (V_LN * 1e-6) + APC0_LN * k_mAPC_death - k_mAPC_death * mAPC_LN
    dmAPC_LNl <- k_APC_mig * mAPC_Ln1 * V_Ln1_vol / (V_LN * 1e-6) + APC0_LN * k_mAPC_death - k_mAPC_death * mAPC_LNl
    dmAPC_LNo <- k_APC_mig * mAPC_oth  * V_oth_vol  / (V_LN * 1e-6) + APC0_LN * k_mAPC_death - k_mAPC_death * mAPC_LNo

    # =========================================================================
    # MODULE 8: NAIVE T CELL TRAFFICKING
    # Naive T cell compartments are held at quasi-steady state (fast equilibration
    # compared to immune activation dynamics). Derivatives set to 0 to avoid
    # stiffness from large thymic output Q_nT0_thym >> k_death * nT0_C.
    # =========================================================================
    dnT0_C   <- 0
    dnT0_P   <- 0
    dnT0_LN  <- 0
    dnT0_LNl <- 0
    dnT0_LNo <- 0
    dnT1_C   <- 0
    dnT1_P   <- 0
    dnT1_LN  <- 0
    dnT1_LNl <- 0
    dnT1_LNo <- 0

    # =========================================================================
    # MODULE 9: T CELL ACTIVATION IN LN
    # =========================================================================
    T_total_LN_all <- (T0_LN + Th_LN
                       + T1_LN + T2_LN + T3_LN + T4_LN + T5_LN + T6_LN + T7_LN + T8_LN)

    dIL2_LN  <- (k_IL2_sec * T_total_LN_all - k_IL2_deg * 1440 * IL2_LN  - k_IL2_cons * T_total_LN_all * IL2_LN)
    dIL2_LNl <- (k_IL2_sec * T_total_LN_all - k_IL2_deg * 1440 * IL2_LNl - k_IL2_cons * T_total_LN_all * IL2_LNl)

    f_act_LN  <- mAPC_LN  / (mAPC_50 + mAPC_LN)
    f_act_LNl <- mAPC_LNl / (mAPC_50 + mAPC_LNl)
    f_act_LNo <- mAPC_LNo / (mAPC_50 + mAPC_LNo)

    daT0_LN  <- k_T0_act * f_act_LN  * H_TGFb_T   * nT0_LN  - k_T0_pro * aT0_LN  - k_T0_death * aT0_LN
    daT0_LNl <- k_T0_act * f_act_LNl * H_TGFb_Ln1 * nT0_LNl - k_T0_pro * aT0_LNl - k_T0_death * aT0_LNl

    dT0_LN  <- k_T0_pro * aT0_LN  * N_aT0 - q_T0_LN_out * T0_LN  - k_T0_death * T0_LN
    dT0_LNl <- k_T0_pro * aT0_LNl * N_aT0 - q_T0_LN_out * T0_LNl - k_T0_death * T0_LNl
    dT0_LNo <- k_T0_act * f_act_LNl * H_TGFb_oth * nT0_LNo - q_T0_LN_out * T0_LNo - k_T0_death * T0_LNo

    daTh_LN  <- k_Th_act * f_act_LN  * nT0_LN  - k_T0_pro * aTh_LN  - k_T0_death * aTh_LN
    daTh_LNl <- k_Th_act * f_act_LNl * nT0_LNl - k_T0_pro * aTh_LNl - k_T0_death * aTh_LNl

    dTh_LN  <- k_T0_pro * aTh_LN  * N_aTh - q_T0_LN_out * Th_LN  - k_T0_death * Th_LN
    dTh_LNl <- k_T0_pro * aTh_LNl * N_aTh - q_T0_LN_out * Th_LNl - k_T0_death * Th_LNl
    dTh_LNo <- k_Th_act * f_act_LNl * nT0_LNo - q_T0_LN_out * Th_LNo - k_T0_death * Th_LNo

    # Per neo-antigen CD8 T activation (j=1..8, fully unrolled)
    .p <- c(p_1, p_2, p_3, p_4, p_5, p_6, p_7, p_8)

    daT1_LN  <- k_T1_act*f_act_LN  *nT1_LN *.p[1] - k_T1_pro*aT1_LN  - k_T1_death*aT1_LN
    daT2_LN  <- k_T1_act*f_act_LN  *nT1_LN *.p[2] - k_T1_pro*aT2_LN  - k_T1_death*aT2_LN
    daT3_LN  <- k_T1_act*f_act_LN  *nT1_LN *.p[3] - k_T1_pro*aT3_LN  - k_T1_death*aT3_LN
    daT4_LN  <- k_T1_act*f_act_LN  *nT1_LN *.p[4] - k_T1_pro*aT4_LN  - k_T1_death*aT4_LN
    daT5_LN  <- k_T1_act*f_act_LN  *nT1_LN *.p[5] - k_T1_pro*aT5_LN  - k_T1_death*aT5_LN
    daT6_LN  <- k_T1_act*f_act_LN  *nT1_LN *.p[6] - k_T1_pro*aT6_LN  - k_T1_death*aT6_LN
    daT7_LN  <- k_T1_act*f_act_LN  *nT1_LN *.p[7] - k_T1_pro*aT7_LN  - k_T1_death*aT7_LN
    daT8_LN  <- k_T1_act*f_act_LN  *nT1_LN *.p[8] - k_T1_pro*aT8_LN  - k_T1_death*aT8_LN

    daT1_LNl <- k_T1_act*f_act_LNl*nT1_LNl*.p[1] - k_T1_pro*aT1_LNl - k_T1_death*aT1_LNl
    daT2_LNl <- k_T1_act*f_act_LNl*nT1_LNl*.p[2] - k_T1_pro*aT2_LNl - k_T1_death*aT2_LNl
    daT3_LNl <- k_T1_act*f_act_LNl*nT1_LNl*.p[3] - k_T1_pro*aT3_LNl - k_T1_death*aT3_LNl
    daT4_LNl <- k_T1_act*f_act_LNl*nT1_LNl*.p[4] - k_T1_pro*aT4_LNl - k_T1_death*aT4_LNl
    daT5_LNl <- k_T1_act*f_act_LNl*nT1_LNl*.p[5] - k_T1_pro*aT5_LNl - k_T1_death*aT5_LNl
    daT6_LNl <- k_T1_act*f_act_LNl*nT1_LNl*.p[6] - k_T1_pro*aT6_LNl - k_T1_death*aT6_LNl
    daT7_LNl <- k_T1_act*f_act_LNl*nT1_LNl*.p[7] - k_T1_pro*aT7_LNl - k_T1_death*aT7_LNl
    daT8_LNl <- k_T1_act*f_act_LNl*nT1_LNl*.p[8] - k_T1_pro*aT8_LNl - k_T1_death*aT8_LNl

    dT1_LN  <- k_T1_pro*aT1_LN *N_aT1 - q_T1_LN_out*T1_LN  - k_T1_death*T1_LN
    dT2_LN  <- k_T1_pro*aT2_LN *N_aT1 - q_T1_LN_out*T2_LN  - k_T1_death*T2_LN
    dT3_LN  <- k_T1_pro*aT3_LN *N_aT1 - q_T1_LN_out*T3_LN  - k_T1_death*T3_LN
    dT4_LN  <- k_T1_pro*aT4_LN *N_aT1 - q_T1_LN_out*T4_LN  - k_T1_death*T4_LN
    dT5_LN  <- k_T1_pro*aT5_LN *N_aT1 - q_T1_LN_out*T5_LN  - k_T1_death*T5_LN
    dT6_LN  <- k_T1_pro*aT6_LN *N_aT1 - q_T1_LN_out*T6_LN  - k_T1_death*T6_LN
    dT7_LN  <- k_T1_pro*aT7_LN *N_aT1 - q_T1_LN_out*T7_LN  - k_T1_death*T7_LN
    dT8_LN  <- k_T1_pro*aT8_LN *N_aT1 - q_T1_LN_out*T8_LN  - k_T1_death*T8_LN

    dT1_LNl <- k_T1_pro*aT1_LNl*N_aT1 - q_T1_LN_out*T1_LNl - k_T1_death*T1_LNl
    dT2_LNl <- k_T1_pro*aT2_LNl*N_aT1 - q_T1_LN_out*T2_LNl - k_T1_death*T2_LNl
    dT3_LNl <- k_T1_pro*aT3_LNl*N_aT1 - q_T1_LN_out*T3_LNl - k_T1_death*T3_LNl
    dT4_LNl <- k_T1_pro*aT4_LNl*N_aT1 - q_T1_LN_out*T4_LNl - k_T1_death*T4_LNl
    dT5_LNl <- k_T1_pro*aT5_LNl*N_aT1 - q_T1_LN_out*T5_LNl - k_T1_death*T5_LNl
    dT6_LNl <- k_T1_pro*aT6_LNl*N_aT1 - q_T1_LN_out*T6_LNl - k_T1_death*T6_LNl
    dT7_LNl <- k_T1_pro*aT7_LNl*N_aT1 - q_T1_LN_out*T7_LNl - k_T1_death*T7_LNl
    dT8_LNl <- k_T1_pro*aT8_LNl*N_aT1 - q_T1_LN_out*T8_LNl - k_T1_death*T8_LNl

    dT1_LNo <- k_T1_act*f_act_LNo*nT1_LNo*.p[1] - q_T1_LN_out*T1_LNo - k_T1_death*T1_LNo
    dT2_LNo <- k_T1_act*f_act_LNo*nT1_LNo*.p[2] - q_T1_LN_out*T2_LNo - k_T1_death*T2_LNo
    dT3_LNo <- k_T1_act*f_act_LNo*nT1_LNo*.p[3] - q_T1_LN_out*T3_LNo - k_T1_death*T3_LNo
    dT4_LNo <- k_T1_act*f_act_LNo*nT1_LNo*.p[4] - q_T1_LN_out*T4_LNo - k_T1_death*T4_LNo
    dT5_LNo <- k_T1_act*f_act_LNo*nT1_LNo*.p[5] - q_T1_LN_out*T5_LNo - k_T1_death*T5_LNo
    dT6_LNo <- k_T1_act*f_act_LNo*nT1_LNo*.p[6] - q_T1_LN_out*T6_LNo - k_T1_death*T6_LNo
    dT7_LNo <- k_T1_act*f_act_LNo*nT1_LNo*.p[7] - q_T1_LN_out*T7_LNo - k_T1_death*T7_LNo
    dT8_LNo <- k_T1_act*f_act_LNo*nT1_LNo*.p[8] - q_T1_LN_out*T8_LNo - k_T1_death*T8_LNo

    # =========================================================================
    # MODULE 10: T CELL DYNAMICS IN TUMOR/PERIPHERAL/CENTRAL
    # =========================================================================
    # Precompute tumor trafficking rates
    .v_T   <- d_T^3   * pi / 6   # volume ~ d^3
    .v_Ln1 <- d_Ln1^3 * pi / 6
    .v_Ln2 <- d_Ln2^3 * pi / 6
    .v_oth <- d_oth^3  * pi / 6

    dT_C_base <- q_T1_LN_out * (T0_LN + T0_LNl + T0_LNo)  # placeholder – overridden below

    # Treg / Th in central/peripheral/tumor
    dT0_C <- (q_T0_LN_out * T0_LN + q_T0_LN_out * T0_LNl + q_T0_LN_out * T0_LNo - k_T0_death * T0_C)
    dT0_P <- -k_T0_death * T0_P

    dTh_C <- (q_T0_LN_out * Th_LN + q_T0_LN_out * Th_LNl + q_T0_LN_out * Th_LNo - k_T0_death * Th_C)
    dTh_P <- -k_T0_death * Th_P

    dT0_T   <- (q_T0_P_in * 1440 * T0_C * V_T_vol   / V_C + k_Th_Treg * H_TGFb_T   * Th_T   - k_T0_death * T0_T)
    dT0_Ln1 <- (q_T0_P_in * 1440 * T0_C * V_Ln1_vol / V_C + k_Th_Treg * H_TGFb_Ln1 * Th_Ln1 - k_T0_death * T0_Ln1)
    dT0_Ln2 <- (q_T0_P_in * 1440 * T0_C * V_Ln2_vol / V_C + k_Th_Treg * H_TGFb_Ln2 * Th_Ln2 - k_T0_death * T0_Ln2)
    dT0_oth <- (q_T0_P_in * 1440 * T0_C * V_oth_vol  / V_C + k_Th_Treg * H_TGFb_oth * Th_oth  - k_T0_death * T0_oth)

    dTh_T   <- q_T0_P_in * 1440 * Th_C * V_T_vol   / V_C - k_T0_death * Th_T   - k_Th_Treg * H_TGFb_T   * Th_T
    dTh_Ln1 <- q_T0_P_in * 1440 * Th_C * V_Ln1_vol / V_C - k_T0_death * Th_Ln1 - k_Th_Treg * H_TGFb_Ln1 * Th_Ln1
    dTh_Ln2 <- q_T0_P_in * 1440 * Th_C * V_Ln2_vol / V_C - k_T0_death * Th_Ln2 - k_Th_Treg * H_TGFb_Ln2 * Th_Ln2
    dTh_oth <- q_T0_P_in * 1440 * Th_C * V_oth_vol  / V_C - k_T0_death * Th_oth  - k_Th_Treg * H_TGFb_oth * Th_oth

    # IL-2 expansion factor in tumor
    .f_IL2 <- IL2_LN / (IL2_50 + IL2_LN)

    # Effector kill contribution to T cell loss (already in cancer ODE)
    # Per-neoantigen CD8 T in central/peripheral/tumor (unrolled j=1..8)
    .q_T_T   <- q_T1_T_in     * 1440
    .q_Ln1_T <- q_T1_T_in_Ln1 * 1440
    .q_Ln2_T <- q_T1_T_in_Ln2 * 1440

    dT1_C <- (q_T1_LN_out*(T1_LN+T1_LNl+T1_LNo) - k_T1_death*T1_C
              - .q_T_T*.v_T*T1_C - .q_Ln1_T*.v_Ln1*T1_C - .q_Ln2_T*.v_Ln2*T1_C
              - q_T1_P_in*T1_C + q_T1_P_out*T1_P)
    dT2_C <- (q_T1_LN_out*(T2_LN+T2_LNl+T2_LNo) - k_T1_death*T2_C
              - .q_T_T*.v_T*T2_C - .q_Ln1_T*.v_Ln1*T2_C - .q_Ln2_T*.v_Ln2*T2_C
              - q_T1_P_in*T2_C + q_T1_P_out*T2_P)
    dT3_C <- (q_T1_LN_out*(T3_LN+T3_LNl+T3_LNo) - k_T1_death*T3_C
              - .q_T_T*.v_T*T3_C - .q_Ln1_T*.v_Ln1*T3_C - .q_Ln2_T*.v_Ln2*T3_C
              - q_T1_P_in*T3_C + q_T1_P_out*T3_P)
    dT4_C <- (q_T1_LN_out*(T4_LN+T4_LNl+T4_LNo) - k_T1_death*T4_C
              - .q_T_T*.v_T*T4_C - .q_Ln1_T*.v_Ln1*T4_C - .q_Ln2_T*.v_Ln2*T4_C
              - q_T1_P_in*T4_C + q_T1_P_out*T4_P)
    dT5_C <- (q_T1_LN_out*(T5_LN+T5_LNl+T5_LNo) - k_T1_death*T5_C
              - .q_T_T*.v_T*T5_C - .q_Ln1_T*.v_Ln1*T5_C - .q_Ln2_T*.v_Ln2*T5_C
              - q_T1_P_in*T5_C + q_T1_P_out*T5_P)
    dT6_C <- (q_T1_LN_out*(T6_LN+T6_LNl+T6_LNo) - k_T1_death*T6_C
              - .q_T_T*.v_T*T6_C - .q_Ln1_T*.v_Ln1*T6_C - .q_Ln2_T*.v_Ln2*T6_C
              - q_T1_P_in*T6_C + q_T1_P_out*T6_P)
    dT7_C <- (q_T1_LN_out*(T7_LN+T7_LNl+T7_LNo) - k_T1_death*T7_C
              - .q_T_T*.v_T*T7_C - .q_Ln1_T*.v_Ln1*T7_C - .q_Ln2_T*.v_Ln2*T7_C
              - q_T1_P_in*T7_C + q_T1_P_out*T7_P)
    dT8_C <- (q_T1_LN_out*(T8_LN+T8_LNl+T8_LNo) - k_T1_death*T8_C
              - .q_T_T*.v_T*T8_C - .q_Ln1_T*.v_Ln1*T8_C - .q_Ln2_T*.v_Ln2*T8_C
              - q_T1_P_in*T8_C + q_T1_P_out*T8_P)

    dT1_P <- q_T1_P_in*T1_C - q_T1_P_out*T1_P - k_T1_death*T1_P
    dT2_P <- q_T1_P_in*T2_C - q_T1_P_out*T2_P - k_T1_death*T2_P
    dT3_P <- q_T1_P_in*T3_C - q_T1_P_out*T3_P - k_T1_death*T3_P
    dT4_P <- q_T1_P_in*T4_C - q_T1_P_out*T4_P - k_T1_death*T4_P
    dT5_P <- q_T1_P_in*T5_C - q_T1_P_out*T5_P - k_T1_death*T5_P
    dT6_P <- q_T1_P_in*T6_C - q_T1_P_out*T6_P - k_T1_death*T6_P
    dT7_P <- q_T1_P_in*T7_C - q_T1_P_out*T7_P - k_T1_death*T7_P
    dT8_P <- q_T1_P_in*T8_C - q_T1_P_out*T8_P - k_T1_death*T8_P

    # T cells in tumor (primary)
    dT1_T <- (.q_T_T*.v_T*T1_C - k_T1_death*T1_T - k_kill_T*(1-H_TGFb_Teff_T)*(1-H_MDSC_T)*T1_T + k_T1_pro*T1_T*.f_IL2)
    dT2_T <- (.q_T_T*.v_T*T2_C - k_T1_death*T2_T - k_kill_T*(1-H_TGFb_Teff_T)*(1-H_MDSC_T)*T2_T + k_T1_pro*T2_T*.f_IL2)
    dT3_T <- (.q_T_T*.v_T*T3_C - k_T1_death*T3_T - k_kill_T*(1-H_TGFb_Teff_T)*(1-H_MDSC_T)*T3_T + k_T1_pro*T3_T*.f_IL2)
    dT4_T <- (.q_T_T*.v_T*T4_C - k_T1_death*T4_T - k_kill_T*(1-H_TGFb_Teff_T)*(1-H_MDSC_T)*T4_T + k_T1_pro*T4_T*.f_IL2)
    dT5_T <- (.q_T_T*.v_T*T5_C - k_T1_death*T5_T - k_kill_T*(1-H_TGFb_Teff_T)*(1-H_MDSC_T)*T5_T + k_T1_pro*T5_T*.f_IL2)
    dT6_T <- (.q_T_T*.v_T*T6_C - k_T1_death*T6_T - k_kill_T*(1-H_TGFb_Teff_T)*(1-H_MDSC_T)*T6_T + k_T1_pro*T6_T*.f_IL2)
    dT7_T <- (.q_T_T*.v_T*T7_C - k_T1_death*T7_T - k_kill_T*(1-H_TGFb_Teff_T)*(1-H_MDSC_T)*T7_T + k_T1_pro*T7_T*.f_IL2)
    dT8_T <- (.q_T_T*.v_T*T8_C - k_T1_death*T8_T - k_kill_T*(1-H_TGFb_Teff_T)*(1-H_MDSC_T)*T8_T + k_T1_pro*T8_T*.f_IL2)

    # T cells in Ln1
    .k_eff_L1 <- k_kill_Ln1 * (1-H_TGFb_Teff_Ln1) * (1-H_MDSC_Ln1)
    dT1_Ln1 <- .q_Ln1_T*.v_Ln1*T1_C - k_T1_death*T1_Ln1 - .k_eff_L1*T1_Ln1
    dT2_Ln1 <- .q_Ln1_T*.v_Ln1*T2_C - k_T1_death*T2_Ln1 - .k_eff_L1*T2_Ln1
    dT3_Ln1 <- .q_Ln1_T*.v_Ln1*T3_C - k_T1_death*T3_Ln1 - .k_eff_L1*T3_Ln1
    dT4_Ln1 <- .q_Ln1_T*.v_Ln1*T4_C - k_T1_death*T4_Ln1 - .k_eff_L1*T4_Ln1
    dT5_Ln1 <- .q_Ln1_T*.v_Ln1*T5_C - k_T1_death*T5_Ln1 - .k_eff_L1*T5_Ln1
    dT6_Ln1 <- .q_Ln1_T*.v_Ln1*T6_C - k_T1_death*T6_Ln1 - .k_eff_L1*T6_Ln1
    dT7_Ln1 <- .q_Ln1_T*.v_Ln1*T7_C - k_T1_death*T7_Ln1 - .k_eff_L1*T7_Ln1
    dT8_Ln1 <- .q_Ln1_T*.v_Ln1*T8_C - k_T1_death*T8_Ln1 - .k_eff_L1*T8_Ln1

    # T cells in Ln2
    .k_eff_L2 <- k_kill_Ln2 * (1-H_TGFb_Teff_Ln2) * (1-H_MDSC_Ln2)
    dT1_Ln2 <- .q_Ln2_T*.v_Ln2*T1_C - k_T1_death*T1_Ln2 - .k_eff_L2*T1_Ln2
    dT2_Ln2 <- .q_Ln2_T*.v_Ln2*T2_C - k_T1_death*T2_Ln2 - .k_eff_L2*T2_Ln2
    dT3_Ln2 <- .q_Ln2_T*.v_Ln2*T3_C - k_T1_death*T3_Ln2 - .k_eff_L2*T3_Ln2
    dT4_Ln2 <- .q_Ln2_T*.v_Ln2*T4_C - k_T1_death*T4_Ln2 - .k_eff_L2*T4_Ln2
    dT5_Ln2 <- .q_Ln2_T*.v_Ln2*T5_C - k_T1_death*T5_Ln2 - .k_eff_L2*T5_Ln2
    dT6_Ln2 <- .q_Ln2_T*.v_Ln2*T6_C - k_T1_death*T6_Ln2 - .k_eff_L2*T6_Ln2
    dT7_Ln2 <- .q_Ln2_T*.v_Ln2*T7_C - k_T1_death*T7_Ln2 - .k_eff_L2*T7_Ln2
    dT8_Ln2 <- .q_Ln2_T*.v_Ln2*T8_C - k_T1_death*T8_Ln2 - .k_eff_L2*T8_Ln2

    # T cells in other met
    .k_eff_oth <- k_kill_oth * (1-H_TGFb_Teff_oth) * (1-H_MDSC_oth)
    dT1_oth <- .q_T_T*.v_oth*T1_C - k_T1_death*T1_oth - .k_eff_oth*T1_oth
    dT2_oth <- .q_T_T*.v_oth*T2_C - k_T1_death*T2_oth - .k_eff_oth*T2_oth
    dT3_oth <- .q_T_T*.v_oth*T3_C - k_T1_death*T3_oth - .k_eff_oth*T3_oth
    dT4_oth <- .q_T_T*.v_oth*T4_C - k_T1_death*T4_oth - .k_eff_oth*T4_oth
    dT5_oth <- .q_T_T*.v_oth*T5_C - k_T1_death*T5_oth - .k_eff_oth*T5_oth
    dT6_oth <- .q_T_T*.v_oth*T6_C - k_T1_death*T6_oth - .k_eff_oth*T6_oth
    dT7_oth <- .q_T_T*.v_oth*T7_C - k_T1_death*T7_oth - .k_eff_oth*T7_oth
    dT8_oth <- .q_T_T*.v_oth*T8_C - k_T1_death*T8_oth - .k_eff_oth*T8_oth

    # =========================================================================
    # ASSEMBLE DERIVATIVES
    # ORDER MUST EXACTLY MATCH IC VARIABLE ORDER in 05_initial_conditions.R
    # =========================================================================
    list(c(
      # --- Pembrolizumab PK (pos 1-9) ---
      daPD1_C, daPD1_P,
      daPD1_T, daPD1_Ln1, daPD1_Ln2, daPD1_oth,
      daPD1_LN, daPD1_LNl, daPD1_LNo,
      # --- Carrying capacity (pos 10-13) ---
      dK_T, dK_Ln1, dK_Ln2, dK_oth,
      # --- Angiogenic factor (pos 14-17) ---
      dcvas_T, dcvas_Ln1, dcvas_Ln2, dcvas_oth,
      # --- Cancer clones (pos 18-37) ---
      dC1_T, dC2_T, dC3_T, dC4_T, dC5_T,
      dC1_Ln1, dC2_Ln1, dC3_Ln1, dC4_Ln1, dC5_Ln1,
      dC1_Ln2, dC2_Ln2, dC3_Ln2, dC4_Ln2, dC5_Ln2,
      dC1_oth, dC2_oth, dC3_oth, dC4_oth, dC5_oth,
      # --- CCL2/MDSC/NO/ArgI (pos 38-53) ---
      dCCL2_T, dCCL2_Ln1, dCCL2_Ln2, dCCL2_oth,
      dMDSC_T, dMDSC_Ln1, dMDSC_Ln2, dMDSC_oth,
      dNO_T, dNO_Ln1, dNO_Ln2, dNO_oth,
      dArgI_T, dArgI_Ln1, dArgI_Ln2, dArgI_oth,
      # --- Macrophages (pos 54-61) ---
      dMac_M1_T, dMac_M2_T,
      dMac_M1_Ln1, dMac_M2_Ln1,
      dMac_M1_Ln2, dMac_M2_Ln2,
      dMac_M1_oth, dMac_M2_oth,
      # --- Cytokines (pos 62-77) ---
      dTGFb_T, dTGFb_Ln1, dTGFb_Ln2, dTGFb_oth,
      dIFNg_T, dIFNg_Ln1, dIFNg_Ln2, dIFNg_oth,
      dIL12_T, dIL12_Ln1, dIL12_Ln2, dIL12_oth,
      dIL10_T, dIL10_Ln1, dIL10_Ln2, dIL10_oth,
      # --- APC maturation signal (pos 78-81) ---
      dc_T, dc_Ln1, dc_Ln2, dc_oth,
      # --- APC (immature) (pos 82-85) ---
      dAPC_T, dAPC_Ln1, dAPC_Ln2, dAPC_oth,
      # --- mAPC (pos 86-92) ---
      dmAPC_T, dmAPC_Ln1, dmAPC_Ln2, dmAPC_oth,
      dmAPC_LN, dmAPC_LNl, dmAPC_LNo,
      # --- Naive T cells: IC order is nT0_C, nT1_C, nT0_P, nT1_P, nT0_LN, nT0_LNl, nT0_LNo, nT1_LN, nT1_LNl, nT1_LNo (pos 93-102) ---
      dnT0_C, dnT1_C, dnT0_P, dnT1_P,
      dnT0_LN, dnT0_LNl, dnT0_LNo,
      dnT1_LN, dnT1_LNl, dnT1_LNo,
      # --- LN activation intermediates (pos 103-114) ---
      daT0_LN, daT0_LNl,
      daTh_LN, daTh_LNl,
      dT0_LN, dT0_LNl, dT0_LNo,
      dTh_LN, dTh_LNl, dTh_LNo,
      # --- IL-2 in LN (pos 113-114) ---
      dIL2_LN, dIL2_LNl,
      # --- T0/Th in central/peripheral (pos 115-118) ---
      dT0_C, dT0_P,
      dTh_C, dTh_P,
      # --- T0/Th in tumors (pos 119-126) ---
      dT0_T, dT0_Ln1, dT0_Ln2, dT0_oth,
      dTh_T, dTh_Ln1, dTh_Ln2, dTh_oth,
      # --- Control flags (pos 127-130, constant) ---
      0, 0, 0, 0,
      # --- Per-neoantigen CD8 T cells (pos 131-218): grouped by antigen j=1..8 ---
      # Each antigen j: aT{j}_LN, aT{j}_LNl, T{j}_LN, T{j}_LNl, T{j}_LNo,
      #                 T{j}_C, T{j}_P, T{j}_T, T{j}_Ln1, T{j}_Ln2, T{j}_oth
      daT1_LN, daT1_LNl, dT1_LN, dT1_LNl, dT1_LNo, dT1_C, dT1_P, dT1_T, dT1_Ln1, dT1_Ln2, dT1_oth,
      daT2_LN, daT2_LNl, dT2_LN, dT2_LNl, dT2_LNo, dT2_C, dT2_P, dT2_T, dT2_Ln1, dT2_Ln2, dT2_oth,
      daT3_LN, daT3_LNl, dT3_LN, dT3_LNl, dT3_LNo, dT3_C, dT3_P, dT3_T, dT3_Ln1, dT3_Ln2, dT3_oth,
      daT4_LN, daT4_LNl, dT4_LN, dT4_LNl, dT4_LNo, dT4_C, dT4_P, dT4_T, dT4_Ln1, dT4_Ln2, dT4_oth,
      daT5_LN, daT5_LNl, dT5_LN, dT5_LNl, dT5_LNo, dT5_C, dT5_P, dT5_T, dT5_Ln1, dT5_Ln2, dT5_oth,
      daT6_LN, daT6_LNl, dT6_LN, dT6_LNl, dT6_LNo, dT6_C, dT6_P, dT6_T, dT6_Ln1, dT6_Ln2, dT6_oth,
      daT7_LN, daT7_LNl, dT7_LN, dT7_LNl, dT7_LNo, dT7_C, dT7_P, dT7_T, dT7_Ln1, dT7_Ln2, dT7_oth,
      daT8_LN, daT8_LNl, dT8_LN, dT8_LNl, dT8_LNo, dT8_C, dT8_P, dT8_T, dT8_Ln1, dT8_Ln2, dT8_oth
    ))
  })
}
