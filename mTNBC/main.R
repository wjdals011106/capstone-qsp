# =============================================================================
# main.R
# TNBC QSP Model - Main Execution Script
# Arulraj et al., Sci. Adv. 9, eadg0289 (2023)
# =============================================================================
# USAGE:
#   source("main.R")            # Run full pipeline
#   source("main.R"); step1()   # Run individual steps
# =============================================================================

# --- Working directory: project root (capstone-qsp/mTNBC) ---
if (!file.exists("R/01_parameters.R")) {
  stop("Run this script from the mTNBC/ directory. Use setwd('path/to/mTNBC')")
}

# --- Load all modules ---
source("R/01_parameters.R")
source("R/02_model_ode.R")
source("R/03_events.R")
source("R/04_rules.R")
source("R/05_initial_conditions.R")
source("R/06_simulation.R")
source("R/07_virtual_patients.R")
source("R/08_biomarker_analysis.R")
source("R/09_visualization.R")
source("R/utils.R")

# --- Create output directories ---
dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("output/results", showWarnings = FALSE, recursive = TRUE)

cat("======================================================\n")
cat("  TNBC QSP Model - Arulraj et al. (2023)\n")
cat("======================================================\n\n")

# =============================================================================
# STEP 1: Validate parameters against Data_S1.xlsx
# =============================================================================
step1_validate <- function() {
  cat("[Step 1] Validating parameters...\n")
  p <- define_parameters()
  validate_parameters(p, "data/adg0289_Data_S1.xlsx")
  cat("  Total parameters defined:", length(p), "\n\n")
  invisible(p)
}

# =============================================================================
# STEP 2: Single patient simulation (quick test)
# =============================================================================
step2_single_patient <- function(t_end = 500, save_plots = TRUE) {
  cat("[Step 2] Running single patient simulation...\n")

  # Set seeding times close for quick test
  p_test <- list(
    delay_Ln1   = 50,
    delay_Ln2   = 80,
    delay_other = 100,
    start_Ln1   = 0,
    start_Ln2   = 0,
    start_other = 0
  )

  result <- run_simulation(p_override = p_test, t_treat_end = t_end, verbose = TRUE)

  if (result$success) {
    cat("  BOR:", result$BOR, "\n")
    cat("  Responder:", result$responder, "\n")
    cat("  SLD baseline:", round(result$SLD_baseline, 3), "cm\n")

    if (save_plots) {
      p1 <- plot_tumor_growth(result, "Single Patient: Tumor Diameters")
      save_figure(p1, "step2_tumor_growth.png")

      p2 <- plot_cytokine_timeseries(result)
      save_figure(p2, "step2_cytokines.png")
    }
  } else {
    cat("  Simulation failed.\n")
  }
  invisible(result)
}

# =============================================================================
# STEP 3: Generate virtual patients (1,000 VPs)
# =============================================================================
step3_generate_vps <- function(n_vp = 1000, seed = 42) {
  cat("[Step 3] Generating", n_vp, "virtual patients (LHS)...\n")
  vp_df <- generate_virtual_patients(n_vp, seed)
  cat("  VP dimensions:", nrow(vp_df), "x", ncol(vp_df), "\n")
  saveRDS(vp_df, "output/results/vp_parameters.rds")
  cat("  Saved: output/results/vp_parameters.rds\n\n")
  invisible(vp_df)
}

# =============================================================================
# STEP 4: Run virtual clinical trial (parallel)
# =============================================================================
step4_virtual_trial <- function(vp_df = NULL, n_cores = 1, t_treat_end = 700) {
  cat("[Step 4] Running virtual clinical trial...\n")

  if (is.null(vp_df)) {
    if (file.exists("output/results/vp_parameters.rds")) {
      vp_df <- readRDS("output/results/vp_parameters.rds")
    } else {
      vp_df <- step3_generate_vps()
    }
  }

  cat("  Running", nrow(vp_df), "simulations on", n_cores, "core(s)...\n")
  t0 <- proc.time()
  results <- run_vp_cohort(vp_df, t_treat_end = t_treat_end, n_cores = n_cores, verbose = TRUE)
  elapsed <- (proc.time() - t0)["elapsed"]
  cat("  Elapsed:", round(elapsed / 60, 1), "minutes\n")

  # Filter successful simulations
  success_idx <- which(sapply(results, function(r) isTRUE(r$success)))
  cat("  Successful simulations:", length(success_idx), "/", length(results), "\n")

  results_ok <- results[success_idx]
  summary    <- summarize_cohort(results_ok)

  cat("\n  --- Response Rate Summary ---\n")
  cat("  CR:", round(summary$CR * 100, 1), "%\n")
  cat("  PR:", round(summary$PR * 100, 1), "%\n")
  cat("  SD:", round(summary$SD * 100, 1), "%\n")
  cat("  PD:", round(summary$PD * 100, 1), "%\n")
  cat("  ORR (CR+PR):", round(summary$ORR * 100, 1), "%\n")
  cat("  Target ORR (KEYNOTE-119): ~11%\n\n")

  saveRDS(results_ok, "output/results/vp_results.rds")
  cat("  Saved: output/results/vp_results.rds\n")

  # Waterfall plot
  vp_summary_df <- data.frame(
    vp_id        = seq_along(results_ok),
    BOR          = sapply(results_ok, function(r) r$BOR),
    SLD_baseline = sapply(results_ok, function(r) r$SLD_baseline),
    SLD_min      = sapply(results_ok, get_SLD_min)
  )
  p_wf <- plot_waterfall(vp_summary_df)
  save_figure(p_wf, "step4_waterfall.png")

  # Response rate comparison (KEYNOTE-119)
  keynote119 <- list(CR = 0.0397, PR = 0.0685, SD = 0.2238, PD = 0.6678)
  p_rr <- plot_response_rates(summary, keynote119)
  save_figure(p_rr, "step4_response_rates.png")

  invisible(list(results = results_ok, summary = summary, vp_df = vp_summary_df))
}

# =============================================================================
# STEP 5: Biomarker analysis
# =============================================================================
step5_biomarker <- function(results_list = NULL, save_output = TRUE) {
  cat("[Step 5] Running biomarker analysis...\n")

  if (is.null(results_list)) {
    if (file.exists("output/results/vp_results.rds")) {
      results_list <- readRDS("output/results/vp_results.rds")
    } else {
      stop("Run step4_virtual_trial() first.")
    }
  }

  # Build biomarker data.frame
  bm_df <- extract_biomarkers_df(results_list)
  cat("  Biomarker data.frame:", nrow(bm_df), "VPs x", ncol(bm_df), "columns\n")

  # Single biomarker analysis
  bm_results <- analyze_all_biomarkers(bm_df)
  cat("\n  --- Top 10 Biomarkers (Response Probability) ---\n")
  print(head(bm_results$by_resp_prob[, c("biomarker","response_prob","RIS","cutoff_rp","direction_rp")], 10))

  cat("\n  --- Top 10 Biomarkers (RIS) ---\n")
  print(head(bm_results$by_ris[, c("biomarker","RIS","cutoff_ris","direction_ris")], 10))

  # Combination analysis
  cat("\n  Running pairwise biomarker combinations...\n")
  combo_df <- analyze_biomarker_combinations(bm_df)
  if (!is.null(combo_df)) {
    cat("  --- Top 5 Biomarker Combinations ---\n")
    print(head(combo_df[, c("bm1","bm2","response_prob","RIS","n_subgroup")], 5))
  }

  if (save_output) {
    saveRDS(bm_results, "output/results/biomarker_single.rds")
    if (!is.null(combo_df)) saveRDS(combo_df, "output/results/biomarker_combinations.rds")
    write.csv(bm_results$by_resp_prob, "output/results/biomarker_ranking.csv", row.names = FALSE)

    p_bm <- plot_biomarker_ranking(bm_results, top_n = 15)
    save_figure(p_bm, "step5_biomarker_ranking.png")
  }

  invisible(list(single = bm_results, combinations = combo_df))
}

# =============================================================================
# FULL PIPELINE
# =============================================================================
run_full_pipeline <- function(n_vp = 1000, n_cores = 1, t_treat_end = 700) {
  cat("Starting full TNBC QSP pipeline...\n\n")
  step1_validate()
  step2_single_patient(save_plots = TRUE)
  vp_df    <- step3_generate_vps(n_vp = n_vp)
  trial    <- step4_virtual_trial(vp_df = vp_df, n_cores = n_cores, t_treat_end = t_treat_end)
  bm_res   <- step5_biomarker(trial$results)
  cat("\n======================================================\n")
  cat("  Pipeline complete!\n")
  cat("  Results saved in output/\n")
  cat("======================================================\n")
  invisible(list(trial = trial, biomarkers = bm_res))
}

# =============================================================================
# AUTO-RUN (comment out to run manually)
# =============================================================================
# Uncomment the line below to run full pipeline automatically:
# run_full_pipeline(n_vp = 100, n_cores = 1, t_treat_end = 500)

# Quick single-patient test:
cat("Ready. To run:\n")
cat("  step1_validate()               # Validate parameters\n")
cat("  step2_single_patient()         # Quick single-patient test\n")
cat("  vp_df <- step3_generate_vps()  # Generate VPs\n")
cat("  step4_virtual_trial(vp_df)     # Run clinical trial\n")
cat("  step5_biomarker()              # Biomarker analysis\n")
cat("  run_full_pipeline(n_vp=100)    # Full run (100 VPs)\n\n")
