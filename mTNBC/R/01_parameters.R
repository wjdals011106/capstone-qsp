# =============================================================================
# 01_parameters.R
# TNBC QSP Model - Parameter Definitions
# Based on: Arulraj et al., Sci. Adv. 9, eadg0289 (2023)
# Data source: adg0289_Data_S1.xlsx - Parameters sheet (737 parameters)
# =============================================================================

define_parameters <- function() {
  p <- list()

  # ---------------------------------------------------------------------------
  # UNIT / STRUCTURAL PARAMETERS
  # ---------------------------------------------------------------------------
  p$cell        <- 1        # cell
  p$day         <- 1        # day
  p$V_Tmin      <- 1e-6     # mL - minimum tumor compartment volume
  p$vol_cell    <- 2572.44  # um^3/cell - cancer cell volume (Abramczyk 2015)
  p$vol_Tcell   <- 175.016  # um^3/cell - T cell volume (Chapman 1981)
  p$vol_Mcell   <- 4849.048 # um^3/cell - macrophage volume (Krombach 1997)
  p$k_cell_clear<- 0.1      # 1/day - dead cell clearance rate
  p$Ve_T        <- 0.37     # dimensionless - tumor void fraction (Finley 2012)

  # ---------------------------------------------------------------------------
  # COMPARTMENT VOLUMES
  # ---------------------------------------------------------------------------
  p$V_C         <- 5        # L - central compartment
  p$V_P         <- 60       # L - peripheral compartment
  p$V_LN        <- 1112.647 # mm^3 - lymph node (Schmidt 2007)
  p$A_syn       <- 37.8     # um^2 - synapse surface area (Jansson 2005)
  p$A_Tcell     <- 151.310  # um^2 - T cell surface area (Chapman 1981)
  p$A_cell      <- 907.920  # um^2 - cancer cell surface area (Abramczyk 2015)
  p$A_APC       <- 900      # um^2 - APC surface area (Agrawal & Linderman 1996)
  p$A_Mcell     <- 1385.442 # um^2 - macrophage surface area (Krombach 1997)
  p$V_e         <- 4e-16    # L - APC endosomal compartment (Agrawal 1996)
  p$A_e         <- 15       # um^2 - APC endosomal surface
  p$A_s         <- 900      # um^2 - APC cell surface

  # ---------------------------------------------------------------------------
  # SIMULATION CONTROL FLAGS
  # ---------------------------------------------------------------------------
  p$seeding     <- 1        # turn on/off metastatic tumor seeding
  p$start       <- 1        # turn on/off tumor evaluation
  p$do_surgery  <- 0        # turn on/off surgical resection of primary tumor
  p$surgery_time<- 200      # day - time of surgery
  p$start_Ln1   <- 0
  p$start_Ln2   <- 0
  p$start_other <- 0
  p$delay_Ln1   <- 2100     # day - seeding time of lung met 1
  p$delay_Ln2   <- 2100     # day - seeding time of lung met 2
  p$delay_other <- 2100     # day - seeding time of other met

  # ---------------------------------------------------------------------------
  # METASTATIC SEEDING: initial cancer cell counts per clone (default 1e5 each)
  # These can be overridden by scale_IC_by_diameter() for VP generation
  # ---------------------------------------------------------------------------
  for (comp in c("Ln1", "Ln2", "other")) {
    for (ci in 1:5) {
      p[[ paste0("ncells_C", ci, "_", comp) ]] <- 1e5
    }
  }

  # ---------------------------------------------------------------------------
  # ANTIGEN CONCENTRATIONS (self-antigen + neoantigens per clone)
  # ---------------------------------------------------------------------------
  ag0 <- 5.4e-13  # mol/cell - self-antigen
  ag_neo <- 6.75e-14  # mol/cell - neoantigen per clone

  for (ci in 1:5) {
    p[[paste0("agconc_C", ci, "_0")]] <- ag0
    for (j in 1:8) {
      p[[paste0("agconc_C", ci, "_", j)]] <- ag_neo
    }
    p[[paste0("totalag_C", ci)]] <- ag0
  }
  # Self-antigen per cell for each clone (P0)
  for (ci in 1:5) {
    p[[paste0("P0_C", ci)]] <- ag0
  }

  # ---------------------------------------------------------------------------
  # CANCER CELL GROWTH / DEATH (all compartments)
  # ---------------------------------------------------------------------------
  p$k_C1_growth       <- 0.0065  # 1/day (Ryu 2014; Desai 2006)
  p$k_C1_death        <- 0.0001  # 1/day (Palsson 2013)
  p$k_C2_growth       <- 0.0065
  p$k_C2_death        <- 0.0001
  p$k_C3_growth       <- 0.0065
  p$k_C3_death        <- 0.0001
  p$k_C4_growth       <- 0.0065
  p$k_C4_death        <- 0.0001
  p$k_C5_growth       <- 0.0065
  p$k_C5_death        <- 0.0001

  for (ci in 1:5) {
    for (comp in c("Ln1","Ln2","other")) {
      p[[paste0("k_C",ci,"_growth_",comp)]] <- 0.0065
      p[[paste0("k_C",ci,"_death_",comp)]]  <- 0.0001
    }
  }

  p$C_max             <- 27000   # cells - cancer cell carrying capacity (Desai 2006)
  p$C_max_Ln1         <- 27000
  p$C_max_Ln2         <- 27000
  p$C_max_other       <- 27000
  p$initial_tumour_diameter <- 2.5  # cm
  p$initial_met_diameter    <- 1.6  # cm (Huang 2021)

  # ---------------------------------------------------------------------------
  # TUMOR VASCULATURE / CARRYING CAPACITY
  # ---------------------------------------------------------------------------
  p$k_K_g       <- 4.12     # 1/day - vasculature growth rate (Desai 2006)
  p$k_K_d       <- 0.0034   # 1/day - vasculature inhibition (Desai 2006; Hahnfeldt 1999)
  p$k_vas_Csec  <- 0.00011  # pg/cell/day - angiogenic factor secretion by cancer cells (Volk 2008)
  p$k_vas_deg   <- 16.6     # 1/day - angiogenic factor degradation (Finley 2011)
  p$c_vas_50    <- 1070     # pg/mL - half-maximal concentration (Desai 2006)
  p$k_vas_Msec  <- 1.7e-5   # pg/cell/day - secretion by macrophage (Wu 2010)

  # ---------------------------------------------------------------------------
  # NAIVE CD4 T CELLS (nT0)
  # ---------------------------------------------------------------------------
  p$div_T0      <- 1160000  # T cell diversity (Robins 2009)
  p$n_T0_clones <- 63       # number of T cell clones
  p$q_nT0_LN_in <- 0.1      # 1/day - naive T cell into LN (Autissier 2010)
  p$q_T0_LN_out <- 24       # 1/day - activated T cell out of LN (De Boer 1995)
  p$q_nT0_LN_out<- 2.88     # 1/day - naive T cell out of LN (Mandl 2012)
  p$k_T0_act    <- 5        # 1/day - T0 activation rate (De Boer 1995)
  p$k_T0_pro    <- 1        # 1/day - T0 proliferation rate (Marchingo 2014)
  p$k_T0_death  <- 0.01     # 1/day
  p$q_T0_P_in   <- 2.436    # 1/min - T0 into peripheral (Zhu 1996)
  p$q_T0_P_out  <- 24       # 1/day - T0 out of peripheral
  p$q_T0_T_in   <- 5.8e-5   # 1/(cm^3*min) - T0 into tumor (Zhu 1996)
  p$q_nT0_P_in  <- 0.1764   # 1/min
  p$q_nT0_P_out <- 5.1      # 1/day (Autissier 2010)
  p$Q_nT0_thym  <- 70e6     # cell/day - thymic output CD4 (Bains 2009)
  p$k_nT0_pro   <- 320e6    # cell/day - naive T0 proliferation (Braber 2012)
  p$K_nT0_pro   <- 1e9      # cell - half-maximal proliferation
  p$k_nT0_death <- 0.002    # 1/day (Braber 2012)

  # ---------------------------------------------------------------------------
  # NAIVE CD8 T CELLS (nT1) - 8 neo-antigen specificities
  # ---------------------------------------------------------------------------
  p$div_T1      <- 1110000  # T cell diversity (Robins 2009)
  p$n_T1_clones <- 0.01991  # Number of T cell clones specific to neoantigen 1
  p$n_T2_clones <- 0.13964
  p$n_T3_clones <- 0.00298
  p$n_T4_clones <- 0.00768
  p$n_T5_clones <- 0.05259
  p$n_T6_clones <- 51.3276
  p$n_T7_clones <- 10.5728
  p$n_T8_clones <- 0.87684
  p$q_nT1_LN_in <- 0.076    # 1/day (Autissier 2010)
  p$q_T1_LN_out <- 24       # 1/day
  p$q_nT1_LN_out<- 1.8      # 1/day (Mandl 2012)
  p$k_T1_act    <- 23       # 1/day (De Boer 1995)
  p$k_T1_pro    <- 1        # 1/day (Marchingo 2014)
  p$k_T1_death  <- 0.01     # 1/day
  p$q_T1_T_in   <- 5.8e-5   # 1/(cm^3*min) - T1 into primary tumor (Zhu 1996)
  p$q_T1_T_in_Ln1   <- 3.5e-5
  p$q_T1_T_in_Ln2   <- 3.5e-5
  p$q_T1_T_in_other <- 5.0e-5
  p$q_T1_P_in   <- 2.436    # 1/min
  p$q_T1_P_out  <- 24       # 1/day
  p$q_nT1_P_in  <- 0.1764   # 1/min
  p$q_nT1_P_out <- 5.1      # 1/day (Autissier 2010)
  p$Q_nT1_thym  <- 35e6     # cell/day (Bains 2009)
  p$k_nT1_pro   <- 320e6    # cell/day (Braber 2012)
  p$K_nT1_pro   <- 1e9      # cell
  p$k_nT1_death <- 0.002    # 1/day (Braber 2012)

  # ---------------------------------------------------------------------------
  # T CELL ACTIVATION / KILLING
  # ---------------------------------------------------------------------------
  p$k_T1        <- 0.1      # 1/day - T cell exhaustion rate (estimated)
  p$k_C_T1      <- 0.95     # 1/day - cancer cell killing by T cells (Robertson-Tessi 2012)
  p$k_Treg      <- 0.05     # 1/day - T cell death by Tregs (estimated)
  p$K_T_C       <- 1.2      # dimensionless - Teff/C ratio dependence (Robertson-Tessi 2012)
  p$K_T_Treg    <- 11       # dimensionless - Teff/Treg ratio dependence
  p$Kc_rec      <- 20.2e6   # cell^2 - half-maximal C for T cell recruitment (Pillis 2006)
  p$N_aT        <- 1        # Number of Activated CD8+ T Cell Generations
  p$N_aT0       <- 1        # Number of Activated Treg Generations

  # ---------------------------------------------------------------------------
  # T CELL DIVISION (Marchingo 2014)
  # ---------------------------------------------------------------------------
  p$N0          <- 2        # divisions by TCR signaling only
  p$N_costim    <- 3        # divisions by co-stimulatory signaling only
  p$N_IL2_CD8   <- 11       # max divisions due to IL-2 (CD8)
  p$N_IL2_CD4   <- 8.5      # max divisions due to IL-2 (CD4)
  p$N_aTh       <- 1        # Activated T Helper Cell Generations
  # LN variants (same values)
  for (comp in c("Ln1","Ln2")) {
    p[[paste0("N0_",comp)]]        <- 2
    p[[paste0("N_costim_",comp)]]  <- 3
    p[[paste0("N_IL2_CD8_",comp)]] <- 11
    p[[paste0("N_IL2_CD4_",comp)]] <- 8.5
    p[[paste0("N_aT_",comp)]]      <- 1
    p[[paste0("N_aT0_",comp)]]     <- 1
    p[[paste0("N_aTh_",comp)]]     <- 1
  }

  # ---------------------------------------------------------------------------
  # IL-2
  # ---------------------------------------------------------------------------
  p$k_IL2_deg   <- 0.2      # 1/min - IL-2 degradation (Lotze 1985)
  p$k_IL2_cons  <- 6e-6     # nmol/cell/hr - IL-2 consumption (Lotze 1985)
  p$k_IL2_sec   <- 3e-5     # nmol/cell/hr - IL-2 secretion (Han 2012)
  p$IL2_50      <- 0.32     # nM - half-maximal IL-2 for T cell activation (Marchingo 2014)
  p$IL2_50_Treg <- 0.32     # nM (Wang & Smith 1987)

  # ---------------------------------------------------------------------------
  # APC DYNAMICS
  # ---------------------------------------------------------------------------
  p$k_APC_mat   <- 1.5      # 1/day - APC maturation rate (Chen 2014)
  p$k_APC_mig   <- 4        # 1/day - APC migration rate (Russo 2016)
  p$k_APC_death <- 0.01     # 1/day - APC death rate (Marino & Kirschner 2004)
  p$k_mAPC_death<- 0.02     # 1/day - mAPC death rate
  p$APC0_T      <- 400000   # cell/mL - APC in primary tumor (Lavin 2017)
  p$APC0_T_Ln1  <- 100000   # cell/mL - APC in met tumor (Siegel 2018)
  p$APC0_T_Ln2  <- 100000
  p$APC0_T_other<- 100000
  p$APC0_LN     <- 1200000  # cell/mL - APC in LN (Catron 2004; Lavin 2017)
  p$n_sites_APC <- 10       # max T cells per APC (De Boer 1995)

  # ---------------------------------------------------------------------------
  # ANTIGEN PROCESSING (MHC, endosome - same for P0 self-antigen)
  # ---------------------------------------------------------------------------
  p$k_c         <- 2        # 1/day - cytokine rate constant (Chen 2014)
  p$c0          <- 1e-9     # M - baseline cytokine concentration (Chen 2014)
  p$c50         <- 1e-9     # M - half-maximal cytokine conc (Chen 2014)
  p$DAMPs       <- 1.34e-14 # mol/cell - DAMPs released per dying cancer cell (Milo 2013)
  p$kin         <- 14.4     # 1/day - MHC internalization (Chen 2014)
  p$kout        <- 28.8     # 1/day - MHC externalization (Chen 2014)
  p$k_P0_up     <- 14.4     # 1/day/cell - antigen uptake (Chen 2014)
  p$k_xP0_deg   <- 2        # 1/day - extracellular antigen degradation (Palsson 2013)
  p$k_P0_deg    <- 17.28    # 1/day - endosomal antigen degradation (Chen 2014)
  p$k_p0_deg    <- 144      # 1/day - endosomal epitope degradation (Chen 2014)
  p$k_P0_on     <- 144000   # 1/(day*M) - antigen-MHC binding (Agrawal 1996)
  p$k_P0_d1     <- 1e-7     # M - antigen-MHC Kd (Stone 2015)
  p$p0_50       <- 2.6455e-5# molecule/um^2 - half-maximal MHC for T cell activation (Kimachi 1997)

  # Neoantigens P1-P8 (same kinetics, different Kd)
  for (j in 1:8) {
    p[[paste0("k_P",j,"_up")]]   <- 14.4
    p[[paste0("k_xP",j,"_deg")]] <- 2
    p[[paste0("k_P",j,"_deg")]]  <- 17.28
    p[[paste0("k_p",j,"_deg")]]  <- 144
    p[[paste0("k_P",j,"_on")]]   <- 144000
    p[[paste0("k_P",j,"_d1")]]   <- 1e-8     # neoantigen Kd (different from self-antigen)
    p[[paste0("p",j,"_50")]]     <- 2.6455e-5
  }

  # TCR kinetics (Lever 2014)
  p$k_M1p0_TCR_on  <- 1    # 1/(s*molecule/um^2)
  p$k_M1p0_TCR_off <- 1    # 1/s
  p$k_M1p0_TCR_p   <- 1    # 1/s
  p$phi_M1p0_TCR   <- 0.09 # 1/s
  p$N_M1p0_TCR     <- 10   # modification steps
  p$TCR_p0_tot     <- 103.813 # molecule/um^2 - TCR per naive T cell (Lever 2014)

  for (j in 1:8) {
    p[[paste0("k_M1p",j,"_TCR_on")]]  <- 1
    p[[paste0("k_M1p",j,"_TCR_off")]] <- 1
    p[[paste0("k_M1p",j,"_TCR_p")]]   <- 1
    p[[paste0("phi_M1p",j,"_TCR")]]   <- 0.09
    p[[paste0("N_M1p",j,"_TCR")]]     <- 10
    p[[paste0("TCR_p",j,"_tot")]]     <- 103.813
  }

  # ---------------------------------------------------------------------------
  # CHECKPOINT MOLECULES - PD-1 / PD-L1 / PD-L2
  # ---------------------------------------------------------------------------
  p$kon_PD1_PDL1    <- 0.06     # 1/(uM*nm*s) (Cheng 2013)
  p$koff_PD1_PDL1   <- 1.476    # 1/s
  p$kon_PD1_PDL2    <- 0.08     # 1/(uM*nm*s) (Cheng 2013)
  p$koff_PD1_PDL2   <- 0.552    # 1/s
  p$kon_PD1_aPD1    <- 670000   # 1/(M*s) (Wang 2014; Brown 2020)
  p$koff_PD1_aPD1   <- 0.00268  # 1/s
  p$kon_PDL1_aPDL1  <- 430000   # 1/(M*s) (Wang 2014)
  p$koff_PDL1_aPDL1 <- 0.000172 # 1/s
  p$Chi_PD1_aPD1    <- 33.333   # 1/nm - cross-arm binding efficiency (Wang 2014)
  p$Chi_PDL1_aPDL1  <- 33.333
  p$Chi_CTLA4_aCTLA4<- 33.333
  p$PD1_50          <- 6        # molecule/um^2 - half-maximal for T cell inactivation (Jafarnejad 2019)
  p$n_PD1           <- 2        # Hill coefficient (estimated)
  p$T_PD1_total     <- 27900    # molecules - PD-1 on T cells (Cheng 2013)
  p$M_PD1_total     <- 12500    # molecules - PD-1 on macrophages (Gordon 2017)

  # PD-L1 on cancer cells
  p$k_out_PDL1      <- 50000    # molecule/day - expression rate (Mimura 2018)
  p$k_in_PDL1       <- 1        # 1/day - degradation rate (Hsu 2018)
  p$r_PDL1_IFNg     <- 6        # fold increase by IFNg (Shin 2017)
  p$C1_PDL1_base    <- 120000   # molecules baseline per cancer cell (Cheng 2013)
  p$r_PDL2C1        <- 0.07     # PDL2/PDL1 ratio (Cheng 2013)
  p$APC_PDL1_base   <- 120000   # molecules baseline per APC
  p$r_PDL2APC       <- 0.07

  # CD28 / CTLA-4 / CD80 / CD86
  p$kon_CD28_CD80   <- 0.1333   # 1/(uM*nm*s) (van der Merwe 1997)
  p$koff_CD28_CD80  <- 1.6      # 1/s
  p$kon_CD28_CD86   <- 0.4667   # 1/(uM*nm*s) (Collins 2002)
  p$koff_CD28_CD86  <- 28       # 1/s
  p$kon_CTLA4_CD80  <- 0.3333   # 1/(uM*nm*s) (van der Merwe 1997)
  p$koff_CTLA4_CD80 <- 0.42     # 1/s
  p$kon_CTLA4_CD86  <- 0.6667   # 1/(uM*nm*s) (Collins 2002)
  p$koff_CTLA4_CD86 <- 5.2      # 1/s
  p$kon_CD80_PDL1   <- 0.1067   # 1/(uM*nm*s) (Cheng 2013)
  p$koff_CD80_PDL1  <- 6.016    # 1/s
  p$kon_CTLA4_aCTLA4<- 383000   # 1/(M*s) (Wang 2014)
  p$koff_CTLA4_aCTLA4<- 0.006971
  p$kon_CD80_CD80   <- 1.9667e-4 # 1/(uM*nm*s) - CD80 self-association (Ikemizu 2000)
  p$koff_CD80_CD80  <- 0.01003  # 1/s
  p$CD28_CD8X_50    <- 200      # molecule/um^2 - half-maximal CD28 co-stim (Parry 2005)
  p$n_CD28_CD8X     <- 2        # Hill coefficient
  p$T_CD28_total    <- 184000   # molecules - CD28 on T cells (Jansson 2005)
  p$T_CTLA4_syn     <- 8000     # molecules - CTLA-4 on T cells (Jansson 2005)
  p$T_PDL1_total    <- 83700    # molecules - PDL1 on T cells
  p$C1_CD80_total   <- 80000    # molecules - CD80 on cancer cells (Jansson 2005)
  p$C1_CD86_total   <- 860000   # molecules - CD86 on cancer cells
  p$APC_CD80_total  <- 80000    # molecules - CD80 on APCs
  p$APC_CD86_total  <- 860000   # molecules - CD86 on APCs

  # ---------------------------------------------------------------------------
  # PEMBROLIZUMAB (anti-PD1) PK
  # ---------------------------------------------------------------------------
  p$q_P_aPD1    <- 9.16e-6  # L/s - central<->peripheral flow (Ahamadi 2016)
  p$q_T_aPD1    <- 8.52e-5  # mL/s - central<->tumor flow (Wang 2019)
  p$q_LN_aPD1   <- 3.25e-6  # mL/s - central<->LN flow (Meijer 2017)
  p$q_LD_aPD1   <- 0.0015   # 1/min - lymphatic drainage TDLN->central (Zhu 1996)
  p$k_cl_aPD1   <- 0.315    # L/day - clearance from central (Ahamadi 2016)
  p$gamma_C_aPD1  <- 0.698  # dimensionless - volume fraction central (Ahamadi 2016)
  p$gamma_P_aPD1  <- 0.0789 # peripheral
  p$gamma_T_aPD1  <- 0.522  # tumor (Coghlin 2010; Finley 2015)
  p$gamma_LN_aPD1 <- 0.2    # LN (Wang 2019)

  # anti-PDL1 PK
  p$q_P_aPDL1   <- 8.7e-6
  p$q_T_aPDL1   <- 8.52e-5
  p$q_LN_aPDL1  <- 3.25e-6
  p$q_LD_aPDL1  <- 0.0015
  p$k_cl_aPDL1  <- 0.324    # L/day (Stroh 2017)
  p$gamma_C_aPDL1  <- 0.61
  p$gamma_P_aPDL1  <- 0.068
  p$gamma_T_aPDL1  <- 0.522
  p$gamma_LN_aPDL1 <- 0.2

  # anti-CTLA4 PK
  p$q_P_aCTLA4  <- 5.331e-6
  p$q_T_aCTLA4  <- 8.52e-5
  p$q_LN_aCTLA4 <- 3.25e-6
  p$q_LD_aCTLA4 <- 0.0015
  p$k_cl_aCTLA4 <- 0.3276   # L/day (Wang 2014)
  p$gamma_C_aCTLA4  <- 0.794
  p$gamma_P_aCTLA4  <- 0.062
  p$gamma_T_aCTLA4  <- 0.522
  p$gamma_LN_aCTLA4 <- 0.2

  # ---------------------------------------------------------------------------
  # HELPER T CELLS (Th) / TREG
  # ---------------------------------------------------------------------------
  p$k_Th_act    <- 10       # 1/day - Th activation rate (Robertson-Tessi 2012)
  p$k_Th_Treg   <- 0.022    # 1/day - Th -> Treg differentiation (Robertson-Tessi 2012)

  # ---------------------------------------------------------------------------
  # TGF-beta
  # ---------------------------------------------------------------------------
  p$k_TGFb_Tsec  <- 1.2e-10  # nmol/cell/day - TGFb secretion by Treg (Liyanage 2002)
  p$k_TGFb_Msec  <- 2e-11    # nmol/cell/day - TGFb secretion by macrophage (Matsuo 2000)
  p$k_TGFb_deg   <- 14.3     # 1/day - TGFb degradation (Robertson-Tessi 2012)
  p$TGFb_50      <- 0.07     # nM - half-maximal TGFb for Th->Treg (Robertson-Tessi 2012)
  p$TGFb_50_Teff <- 0.14     # nM - half-maximal TGFb for CD8 inhibition
  p$TGFbase      <- 0.016    # nM - baseline TGFb in breast tumor (Ivanovic 2003)

  # ---------------------------------------------------------------------------
  # IFN-gamma
  # ---------------------------------------------------------------------------
  p$k_IFNg_sec   <- 5e-13    # nmol/cell/day - IFNg secretion by Th (Autenshlyus 2021)
  p$k_IFNg_deg   <- 7.68     # 1/day - IFNg degradation (Hofstra 1998)
  p$IFNg_50_ind  <- 2.96     # pM - half-maximal IFNg for PDL1 induction (Shin 2017)
  p$IFNg_50      <- 2.9      # pM - half-maximal IFNg for M2->M1 (Altan-Bonnet 2019)

  # ---------------------------------------------------------------------------
  # CCL2
  # ---------------------------------------------------------------------------
  p$k_CCL2_sec   <- 1.7e-12  # nmol/cell/day - CCL2 secretion (Huang 2007; Dutta 2018)
  p$k_CCL2_deg   <- 0.06     # 1/hr - CCL2 degradation (Tanimoto 2007)
  p$CCL2_50      <- 0.23     # nM - half-maximal CCL2 for MDSC recruitment (Opalek 2007)

  # ---------------------------------------------------------------------------
  # MDSC
  # ---------------------------------------------------------------------------
  p$k_MDSC_mig        <- 11000   # cell/(mL*day) - MDSC recruitment (Diaz-Montero 2009)
  p$k_MDSC_mig_Ln1    <- 11000
  p$k_MDSC_mig_Ln2    <- 11000
  p$k_MDSC_mig_other  <- 65000   # (Siegel 2018)
  p$k_MDSC_death      <- 0.015   # 1/day - MDSC death rate (Lai 2018)
  p$k_NO_deg          <- 135     # 1/day - NO degradation (Hakim 1996)
  p$k_ArgI_deg        <- 0.173   # 1/day - ArgI degradation (Schimke 1964)
  p$k_NO_sec          <- 4.8e-7  # nmol/cell/day - NO secretion by MDSC (Serafini 2008)
  p$k_ArgI_sec        <- 0.014   # mU*uL/cell/day - ArgI secretion (Serafini 2008)
  p$ArgI_50_Teff      <- 61.7    # mU - ArgI for T cell inhibition (Serafini 2008)
  p$NO_50_Teff        <- 0.75    # nM - NO for T cell inhibition (Serafini 2008)
  p$ArgI_50_Treg      <- 22.1    # mU - ArgI for Treg expansion (estimated)

  # ---------------------------------------------------------------------------
  # MACROPHAGES
  # ---------------------------------------------------------------------------
  p$k_Mac_mig         <- 170000  # cell/(mL*day) - macrophage recruitment (Yang 2018)
  p$k_Mac_mig_Ln1     <- 350000  # (Siegel 2018)
  p$k_Mac_mig_Ln2     <- 350000
  p$k_Mac_mig_other   <- 300000
  p$k_M2_pol          <- 0.25    # 1/day - M1->M2 polarization (estimated)
  p$k_M2_pol_Ln1      <- 1.1
  p$k_M2_pol_Ln2      <- 1.1
  p$k_M2_pol_other    <- 0.4
  p$k_M1_pol          <- 0.045   # 1/day - M2->M1 polarization
  p$k_Mac_death       <- 0.02    # 1/day - macrophage death (Ginhoux 2016)
  p$k_M1_phago        <- 0.33    # 1/day - M1 phagocytosis of cancer (Pan 2020)
  p$K_Mac_C           <- 2       # dimensionless - macrophage-cancer association (estimated)

  # Macrophage checkpoint/cytokine
  p$kon_CD47_SIRPa    <- 0.0625  # 1/(uM*min*nm) (Hayes 2020)
  p$koff_CD47_SIRPa   <- 0.3     # 1/min
  p$SIRPa_50          <- 37      # molecule/um^2 (Willingham 2012)
  p$n_SIRPa           <- 2       # Hill coefficient (estimated)
  p$C_CD47            <- 400     # molecule/um^2 - CD47 on cancer cells (Morrissey 2020)
  p$M_SIRPa           <- 100     # molecule/um^2 - SIRPa on macrophages (Subramanian 2006)
  p$IL10_50_phago     <- 270     # pM - IL-10 half-max for phagocytosis inhibition (Bian 2016)

  # ---------------------------------------------------------------------------
  # IL-12 / IL-10
  # ---------------------------------------------------------------------------
  p$k_IL12_sec        <- 8.5e-12  # nmol/cell/day - IL-12 secretion by mAPC (Heufler 1996)
  p$k_IL12_Msec       <- 5e-13    # nmol/cell/day - IL-12 secretion by macrophage (Clough 2007)
  p$k_IL12_deg        <- 0.0231   # 1/hr - IL-12 degradation (Carreno 2000)
  p$k_IL10_sec        <- 3e-12    # nmol/cell/day - IL-10 secretion by macrophage (Autenshlyus 2021)
  p$k_IL10_deg        <- 4        # 1/day - IL-10 degradation (Saxena 2015)
  p$IL10_50           <- 8        # pM - IL-10 half-max for M1->M2 (Altan-Bonnet 2019)
  p$IL12_50           <- 0.14     # pM - IL-12 half-max for M2->M1 (Altan-Bonnet 2019)

  # ---------------------------------------------------------------------------
  # DOSING (Pembrolizumab - KEYNOTE-119)
  # ---------------------------------------------------------------------------
  p$dose_pembrolizumab <- 200   # mg - dose per administration
  p$MW_pembrolizumab   <- 146700 # Da (g/mol) - molecular weight of pembrolizumab
  p$dose_interval      <- 21    # days - q3w
  p$n_doses            <- 35    # maximum number of doses (~2 years)

  # ---------------------------------------------------------------------------
  # PARAMETERS USED IN ODE BUT DEFINED SEPARATELY
  # ---------------------------------------------------------------------------
  # mAPC half-saturation for T cell activation in LN
  p$mAPC_50    <- 100       # cells - mAPC half-max for T act (Luber 2010)

  # Per-neoantigen TCR frequency (8 neoantigens, sum to ~0.04 of naive T cells)
  # Equal distribution by default: each neoantigen ~0.005 of nT1 pool
  p$p_1 <- 0.125; p$p_2 <- 0.125; p$p_3 <- 0.125; p$p_4 <- 0.125
  p$p_5 <- 0.125; p$p_6 <- 0.125; p$p_7 <- 0.125; p$p_8 <- 0.125

  # Expansion ratio: activated T cell -> effector T cells (per neoantigen)
  p$N_aT1   <- 2000   # fold expansion (Gattinoni 2009)

  # Initial tumor diameters (used as algebraic variables in ODE, but also stored)
  # These are computed dynamically from state - provide default initialization
  p$d_T_init   <- 3.0  # cm - primary tumor initial diameter
  p$d_Ln1_init <- 0.1  # cm
  p$d_Ln2_init <- 0.1  # cm
  p$d_oth_init <- 0.1  # cm

  return(p)
}
