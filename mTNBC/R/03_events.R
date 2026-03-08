# =============================================================================
# 03_events.R
# TNBC QSP Model - Discrete Events
# - Pembrolizumab dosing (200 mg IV q3w)
# - Metastatic tumor seeding (at delay_Ln1, delay_Ln2, delay_other)
# - Primary tumor surgery (at surgery_time, if do_surgery=1)
# =============================================================================

# -----------------------------------------------------------------------------
# build_event_data()
# Returns a data.frame of events for deSolve::lsoda (eventfun or event$data)
# -----------------------------------------------------------------------------
build_event_data <- function(p, t_start = 0, t_end = 700) {

  events <- list()

  # ---------------------------------------------------------------------------
  # 1. PEMBROLIZUMAB DOSING (200 mg IV q3w)
  # ---------------------------------------------------------------------------
  # Convert dose to nM in central compartment (V_C = 5 L)
  # 200 mg / 146700 g/mol = 1.3633e-3 mmol = 1.3633e-6 mol
  # [nM] = mol / V_C(L) * 1e9 = 1.3633e-6 / 5 * 1e9 = 272.65 nM
  dose_nM <- (p$dose_pembrolizumab / p$MW_pembrolizumab * 1e3) / p$V_C * 1e9  # nM

  dose_times <- seq(t_start, t_end, by = p$dose_interval)
  dose_times <- dose_times[dose_times <= t_end]

  for (td in dose_times) {
    events[[ length(events) + 1 ]] <- data.frame(
      var   = "aPD1_C",
      time  = td,
      value = dose_nM,
      method = "add"
    )
  }

  # ---------------------------------------------------------------------------
  # 2. METASTATIC TUMOR SEEDING
  # ---------------------------------------------------------------------------
  # At seeding time, set start flag and add initial cancer cells
  # (start_Ln1, start_Ln2, start_other switch to 1)
  seeding_vars <- list(
    list(delay = p$delay_Ln1,   flag = "start_Ln1",   clones = c("C1_Ln1","C2_Ln1","C3_Ln1","C4_Ln1","C5_Ln1"),
         ncells = c(p$ncells_C1_Ln1, p$ncells_C2_Ln1, p$ncells_C3_Ln1, p$ncells_C4_Ln1, p$ncells_C5_Ln1)),
    list(delay = p$delay_Ln2,   flag = "start_Ln2",   clones = c("C1_Ln2","C2_Ln2","C3_Ln2","C4_Ln2","C5_Ln2"),
         ncells = c(p$ncells_C1_Ln2, p$ncells_C2_Ln2, p$ncells_C3_Ln2, p$ncells_C4_Ln2, p$ncells_C5_Ln2)),
    list(delay = p$delay_other, flag = "start_other", clones = c("C1_oth","C2_oth","C3_oth","C4_oth","C5_oth"),
         ncells = c(p$ncells_C1_other, p$ncells_C2_other, p$ncells_C3_other, p$ncells_C4_other, p$ncells_C5_other))
  )

  for (sv in seeding_vars) {
    if (sv$delay >= t_start && sv$delay <= t_end) {
      # Turn on start flag
      events[[ length(events) + 1 ]] <- data.frame(
        var   = sv$flag,
        time  = sv$delay,
        value = 1,
        method = "replace"
      )
      # Seed cancer cells
      for (k in seq_along(sv$clones)) {
        events[[ length(events) + 1 ]] <- data.frame(
          var   = sv$clones[k],
          time  = sv$delay,
          value = sv$ncells[k],
          method = "replace"
        )
      }
    }
  }

  # ---------------------------------------------------------------------------
  # 3. PRIMARY TUMOR SURGERY
  # ---------------------------------------------------------------------------
  if (p$do_surgery == 1 && p$surgery_time >= t_start && p$surgery_time <= t_end) {
    for (ci in c("C1_T","C2_T","C3_T","C4_T","C5_T")) {
      events[[ length(events) + 1 ]] <- data.frame(
        var   = ci,
        time  = p$surgery_time,
        value = 0,
        method = "replace"
      )
    }
    events[[ length(events) + 1 ]] <- data.frame(
      var   = "start",
      time  = p$surgery_time,
      value = 0,
      method = "replace"
    )
  }

  # Combine all events
  if (length(events) == 0) return(NULL)
  ev_df <- do.call(rbind, events)
  ev_df <- ev_df[order(ev_df$time), ]
  rownames(ev_df) <- NULL
  return(ev_df)
}

# -----------------------------------------------------------------------------
# event_function() for deSolve (alternative approach using eventfun)
# This is called by ode() at each event time
# -----------------------------------------------------------------------------
make_event_function <- function(p, t_start = 0, t_end = 700) {

  dose_nM    <- (p$dose_pembrolizumab / p$MW_pembrolizumab * 1e3) / p$V_C * 1e9
  dose_times <- seq(t_start, t_end, by = p$dose_interval)

  # Capture all needed values from p into the closure (avoids parms$ on atomic vector)
  delay_Ln1    <- p$delay_Ln1
  delay_Ln2    <- p$delay_Ln2
  delay_other  <- p$delay_other
  do_surgery   <- p$do_surgery
  surgery_time <- p$surgery_time
  seeding      <- p$seeding

  ncells_C1_Ln1   <- p$ncells_C1_Ln1
  ncells_C2_Ln1   <- p$ncells_C2_Ln1
  ncells_C3_Ln1   <- p$ncells_C3_Ln1
  ncells_C4_Ln1   <- p$ncells_C4_Ln1
  ncells_C5_Ln1   <- p$ncells_C5_Ln1

  ncells_C1_Ln2   <- p$ncells_C1_Ln2
  ncells_C2_Ln2   <- p$ncells_C2_Ln2
  ncells_C3_Ln2   <- p$ncells_C3_Ln2
  ncells_C4_Ln2   <- p$ncells_C4_Ln2
  ncells_C5_Ln2   <- p$ncells_C5_Ln2

  ncells_C1_other <- p$ncells_C1_other
  ncells_C2_other <- p$ncells_C2_other
  ncells_C3_other <- p$ncells_C3_other
  ncells_C4_other <- p$ncells_C4_other
  ncells_C5_other <- p$ncells_C5_other

  function(t, y, parms) {
    # Pembrolizumab dose
    if (any(abs(t - dose_times) < 0.01)) {
      y["aPD1_C"] <- y["aPD1_C"] + dose_nM
    }
    # Metastatic seeding (Ln1)
    if (abs(t - delay_Ln1) < 0.01 && seeding == 1) {
      y["start_Ln1"] <- 1
      y["C1_Ln1"] <- ncells_C1_Ln1
      y["C2_Ln1"] <- ncells_C2_Ln1
      y["C3_Ln1"] <- ncells_C3_Ln1
      y["C4_Ln1"] <- ncells_C4_Ln1
      y["C5_Ln1"] <- ncells_C5_Ln1
    }
    # Metastatic seeding (Ln2)
    if (abs(t - delay_Ln2) < 0.01 && seeding == 1) {
      y["start_Ln2"] <- 1
      y["C1_Ln2"] <- ncells_C1_Ln2
      y["C2_Ln2"] <- ncells_C2_Ln2
      y["C3_Ln2"] <- ncells_C3_Ln2
      y["C4_Ln2"] <- ncells_C4_Ln2
      y["C5_Ln2"] <- ncells_C5_Ln2
    }
    # Metastatic seeding (other)
    if (abs(t - delay_other) < 0.01 && seeding == 1) {
      y["start_other"] <- 1
      y["C1_oth"] <- ncells_C1_other
      y["C2_oth"] <- ncells_C2_other
      y["C3_oth"] <- ncells_C3_other
      y["C4_oth"] <- ncells_C4_other
      y["C5_oth"] <- ncells_C5_other
    }
    # Surgery
    if (do_surgery == 1 && abs(t - surgery_time) < 0.01) {
      y["C1_T"] <- 0; y["C2_T"] <- 0; y["C3_T"] <- 0
      y["C4_T"] <- 0; y["C5_T"] <- 0; y["start"] <- 0
    }
    return(y)
  }
}

# Helper: get all event times for ode() times argument
get_event_times <- function(p, t_start = 0, t_end = 700) {
  dose_times <- seq(t_start, t_end, by = p$dose_interval)
  seed_times <- c(p$delay_Ln1, p$delay_Ln2, p$delay_other)
  seed_times <- seed_times[seed_times >= t_start & seed_times <= t_end]
  surg_time  <- if (p$do_surgery == 1) p$surgery_time else numeric(0)
  all_times  <- sort(unique(c(dose_times, seed_times, surg_time)))
  all_times[all_times >= t_start & all_times <= t_end]
}
