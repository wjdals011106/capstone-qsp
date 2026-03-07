# =============================================================================
# utils.R
# TNBC QSP Model - Utility Functions
# =============================================================================

# =============================================================================
# Shannon diversity index
# =============================================================================
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

# =============================================================================
# Safe division (returns 0 if denominator is 0 or NA)
# =============================================================================
safe_div <- function(a, b, fill = 0) {
  ifelse(is.na(b) | b == 0, fill, a / b)
}

# =============================================================================
# Hill function  f = x^n / (k^n + x^n)
# =============================================================================
hill <- function(x, k, n = 1) {
  x  <- pmax(x, 0)
  xn <- x^n
  kn <- k^n
  xn / (kn + xn + 1e-300)
}

# =============================================================================
# Summarize VP cohort results
# =============================================================================
summarize_cohort <- function(results_list) {
  bor_table <- table(sapply(results_list, function(r) r$BOR))
  n_total   <- length(results_list)

  rates <- list(
    CR = as.numeric(bor_table["CR"]) / n_total,
    PR = as.numeric(bor_table["PR"]) / n_total,
    SD = as.numeric(bor_table["SD"]) / n_total,
    PD = as.numeric(bor_table["PD"]) / n_total
  )
  rates[sapply(rates, is.na)] <- 0
  rates$ORR <- rates$CR + rates$PR
  rates$n_total <- n_total
  rates$bor_table <- bor_table
  rates
}

# =============================================================================
# Extract biomarker data.frame from results list
# =============================================================================
extract_biomarkers_df <- function(results_list) {
  bm_list <- lapply(results_list, function(r) {
    bm <- r$biomarkers
    bm$BOR       <- r$BOR
    bm$responder <- r$responder
    as.data.frame(bm)
  })
  df <- do.call(plyr::rbind.fill, bm_list)
  df
}

# =============================================================================
# Compute SLD minimum from time series
# =============================================================================
get_SLD_min <- function(result) {
  if (!result$success) return(NA)
  min(result$SLD_series$SLD, na.rm = TRUE)
}

# =============================================================================
# Load parameters from Data_S1 Excel (validation helper)
# =============================================================================
load_data_s1_params <- function(filepath) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' required. Install with: install.packages('readxl')")
  }
  readxl::read_excel(filepath, sheet = "Parameters")
}

# =============================================================================
# Compare simulated vs expected parameter values
# =============================================================================
validate_parameters <- function(p, data_s1_path = "data/adg0289_Data_S1.xlsx") {
  if (!file.exists(data_s1_path)) {
    message("Data_S1 not found at: ", data_s1_path)
    return(invisible(NULL))
  }
  excel_params <- load_data_s1_params(data_s1_path)
  colnames(excel_params) <- c("Name","Value","Unit","Description")

  mismatches <- c()
  for (i in seq_len(nrow(excel_params))) {
    nm  <- excel_params$Name[i]
    val <- suppressWarnings(as.numeric(excel_params$Value[i]))
    if (is.na(val)) next
    if (!is.null(p[[nm]])) {
      rel_diff <- abs(p[[nm]] - val) / (abs(val) + 1e-300)
      if (rel_diff > 0.01) {
        mismatches <- c(mismatches, sprintf("  %s: model=%.4g, excel=%.4g", nm, p[[nm]], val))
      }
    }
  }

  if (length(mismatches) > 0) {
    message("Parameter mismatches (>1%):\n", paste(mismatches, collapse="\n"))
  } else {
    message("All parameters match Data_S1 within 1%.")
  }
  invisible(mismatches)
}

# =============================================================================
# Parallel VP simulation runner
# =============================================================================
run_vp_cohort <- function(vp_df, t_treat_end = 700, n_cores = 1,
                           verbose = FALSE) {
  source("R/06_simulation.R", local = TRUE)

  run_one <- function(i) {
    vp_row <- as.list(vp_df[i, ])
    p_vp   <- build_vp_params(vp_row)
    tryCatch(
      run_simulation(p_override    = p_vp,
                     t_treat_start = 0,
                     t_treat_end   = t_treat_end,
                     t_step        = 1,
                     verbose       = verbose),
      error = function(e) list(success = FALSE, error = e$message)
    )
  }

  if (n_cores > 1 && requireNamespace("parallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, c("vp_df","t_treat_end","verbose"), envir = environment())
    parallel::clusterEvalQ(cl, {
      source("R/01_parameters.R"); source("R/02_model_ode.R")
      source("R/03_events.R");     source("R/04_rules.R")
      source("R/05_initial_conditions.R"); source("R/06_simulation.R")
      source("R/07_virtual_patients.R");   source("R/utils.R")
    })
    results <- parallel::parLapply(cl, seq_len(nrow(vp_df)), run_one)
  } else {
    results <- lapply(seq_len(nrow(vp_df)), function(i) {
      if (verbose && i %% 50 == 0) cat("VP", i, "of", nrow(vp_df), "\n")
      run_one(i)
    })
  }
  results
}
