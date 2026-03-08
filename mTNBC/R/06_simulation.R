# =============================================================================
# 06_simulation.R
# TNBC QSP Model - Simulation Runner
# =============================================================================

library(deSolve)

source("R/01_parameters.R")
source("R/02_model_ode.R")
source("R/03_events.R")
source("R/04_rules.R")
source("R/05_initial_conditions.R")
source("R/utils.R")

# =============================================================================
# run_simulation()
# Runs a single patient simulation (pre-treatment growth + treatment)
# Returns a list with time-series and summary metrics
# =============================================================================
run_simulation <- function(p_override = list(),
                            t_treat_start = 0,
                            t_treat_end   = 700,
                            t_step        = 1,
                            solver        = "lsoda",
                            verbose       = FALSE) {

  # --- Build parameters ---
  p <- define_parameters()
  for (nm in names(p_override)) p[[nm]] <- p_override[[nm]]

  # --- Build initial conditions ---
  IC <- define_initial_conditions(p)
  parms <- unlist(p)

  # --- Time grid ---
  times <- seq(t_treat_start, t_treat_end, by = t_step)
  # Add event times to ensure they are hit
  ev_times <- get_event_times(p, t_treat_start, t_treat_end)
  times    <- sort(unique(c(times, ev_times)))

  # --- Event function ---
  eventfun <- make_event_function(p, t_treat_start, t_treat_end)

  # --- Run ODE ---
  out_df <- NULL
  tryCatch({
    out <- ode(
      y       = IC,
      times   = times,
      func    = tnbc_ode,
      parms   = parms,
      method  = "vode",     # Adams/BDF, well-suited for stiff large systems
      events  = list(func = eventfun, time = ev_times),
      atol    = 1e-6,       # absolute tolerance
      rtol    = 1e-4,       # relative tolerance
      hmax    = 1.0,        # max step = 1 day
      maxsteps= 500000
    )
    out_df <- as.data.frame(out)
  }, error = function(e) {
    if (verbose) message("ODE solver error: ", e$message)
  })

  if (is.null(out_df)) return(list(success = FALSE))

  # --- Compute derived rules at each time point ---
  state_names <- names(IC)
  rules_list  <- lapply(seq_len(nrow(out_df)), function(i) {
    state_i <- setNames(as.numeric(out_df[i, state_names]), state_names)
    compute_rules(state_i, p)
  })

  # Convert rules list to data.frame
  rules_df <- do.call(rbind, lapply(rules_list, function(r) {
    as.data.frame(lapply(r, function(x) if (length(x) == 1) x else NA))
  }))

  # Combine
  result_df <- cbind(out_df["time"], rules_df)

  # --- Baseline SLD (at t=0, or first evaluation) ---
  SLD_baseline <- result_df$SLD[1]
  if (SLD_baseline <= 0) SLD_baseline <- 0.1  # avoid division by zero

  # --- RECIST assessment every 9 weeks (63 days) ---
  eval_times <- seq(t_treat_start + 63, t_treat_end, by = 63)
  eval_idx   <- sapply(eval_times, function(te) which.min(abs(result_df$time - te)))
  SLD_series <- result_df$SLD[eval_idx]
  t_series   <- result_df$time[eval_idx]

  BOR <- best_overall_response(SLD_series, t_series, SLD_baseline)

  # --- Baseline biomarkers (at t_treat_start) ---
  bm <- as.list(rules_df[1, ])

  # --- Summary output ---
  list(
    success      = TRUE,
    timeseries   = result_df,
    BOR          = BOR,
    responder    = BOR %in% c("CR","PR","SD"),
    SLD_baseline = SLD_baseline,
    SLD_series   = data.frame(time = t_series, SLD = SLD_series),
    biomarkers   = bm,
    parameters   = p
  )
}


# =============================================================================
# run_single_patient_example()
# Quick single-patient test run
# =============================================================================
run_single_patient_example <- function() {
  cat("Running single patient simulation...\n")

  p_over <- list(
    delay_Ln1   = 100,  # seed lung met 1 at day 100
    delay_Ln2   = 150,  # seed lung met 2 at day 150
    delay_other = 200,  # seed other met at day 200
    start_Ln1   = 0,
    start_Ln2   = 0,
    start_other = 0
  )

  # Precompute ncells from initial_met_diameter
  p_base <- define_parameters()
  for (nm in names(p_over)) p_base[[nm]] <- p_over[[nm]]
  p_base <- scale_IC_by_diameter(define_initial_conditions(p_base), p_base)

  result <- run_simulation(
    p_override    = p_over,
    t_treat_start = 0,
    t_treat_end   = 500,
    t_step        = 1,
    verbose       = TRUE
  )

  if (result$success) {
    cat("  Best Overall Response:", result$BOR, "\n")
    cat("  Responder:", result$responder, "\n")
    cat("  SLD baseline:", round(result$SLD_baseline, 3), "cm\n")
    cat("  SLD at last eval:", round(tail(result$SLD_series$SLD, 1), 3), "cm\n")

    # Quick tumor growth plot
    ts <- result$timeseries
    cat("  Primary tumor diameter range:", round(range(ts$d_T, na.rm=TRUE), 3), "\n")
    cat("  Lung met 1 diameter range:   ", round(range(ts$d_Ln1, na.rm=TRUE), 3), "\n")
  } else {
    cat("  Simulation failed.\n")
  }

  invisible(result)
}
