# =============================================================================
# 09_visualization.R
# TNBC QSP Model - Visualization
# Reproduces key figures from Arulraj et al. (2023)
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

# Default theme
theme_tnbc <- function() {
  theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      legend.position  = "right"
    )
}

# =============================================================================
# plot_tumor_growth()
# Time-series of tumor diameters (single patient)
# =============================================================================
plot_tumor_growth <- function(result, title = "Tumor Growth") {
  ts <- result$timeseries
  df <- ts %>%
    select(time, d_T, d_Ln1, d_Ln2, d_oth) %>%
    pivot_longer(-time, names_to = "compartment", values_to = "diameter") %>%
    mutate(compartment = recode(compartment,
      d_T   = "Primary tumor",
      d_Ln1 = "Lung Met 1",
      d_Ln2 = "Lung Met 2",
      d_oth = "Other Met"))

  ggplot(df, aes(x = time, y = diameter, color = compartment)) +
    geom_line(size = 1) +
    labs(title = title, x = "Time (days)", y = "Diameter (cm)",
         color = "Compartment") +
    theme_tnbc()
}

# =============================================================================
# plot_immune_timeseries()
# Time-series of key immune species
# =============================================================================
plot_immune_timeseries <- function(result, compartment = "T") {
  ts <- result$timeseries
  suffix <- paste0("_", compartment)

  df <- ts %>%
    select(time,
           matches(paste0("Tcyt_total", suffix)),
           matches(paste0("T0", suffix)),
           matches(paste0("Mac_M1", suffix)),
           matches(paste0("Mac_M2", suffix)),
           matches(paste0("MDSC", suffix))) %>%
    pivot_longer(-time, names_to = "species", values_to = "count")

  ggplot(df, aes(x = time, y = count, color = species)) +
    geom_line(size = 1) +
    scale_y_log10() +
    labs(title = paste("Immune cells in", compartment, "compartment"),
         x = "Time (days)", y = "Cell count (log scale)",
         color = "Species") +
    theme_tnbc()
}

# =============================================================================
# plot_cytokine_timeseries()
# Time-series of cytokines in tumor
# =============================================================================
plot_cytokine_timeseries <- function(result) {
  ts <- result$timeseries

  df <- ts %>%
    select(time, TGFb_T, IFNg_T, IL12_T, IL10_T, CCL2_T) %>%
    pivot_longer(-time, names_to = "cytokine", values_to = "concentration") %>%
    mutate(cytokine = recode(cytokine,
      TGFb_T  = "TGF-β", IFNg_T = "IFN-γ",
      IL12_T  = "IL-12", IL10_T = "IL-10", CCL2_T = "CCL2"))

  ggplot(df, aes(x = time, y = concentration, color = cytokine)) +
    geom_line(size = 1) +
    labs(title = "Cytokines in primary tumor",
         x = "Time (days)", y = "Concentration (nM)",
         color = "Cytokine") +
    theme_tnbc()
}

# =============================================================================
# plot_waterfall()
# Waterfall plot of best % change in SLD (Fig. 3 style)
# =============================================================================
plot_waterfall <- function(vp_df) {
  # vp_df: rows = VPs, cols include SLD_baseline, SLD_min, BOR
  vp_sorted <- vp_df %>%
    mutate(pct_change = (SLD_min - SLD_baseline) / SLD_baseline * 100) %>%
    arrange(pct_change) %>%
    mutate(vp_id = row_number())

  ggplot(vp_sorted, aes(x = vp_id, y = pct_change, fill = BOR)) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -30, linetype = "dashed", color = "blue") +
    geom_hline(yintercept =  20, linetype = "dashed", color = "red") +
    scale_fill_manual(values = c(CR = "#2E86AB", PR = "#A0CFA8", SD = "#F6AE2D", PD = "#F26419")) +
    labs(title = "Waterfall Plot - Best % Change in SLD",
         x = "Virtual Patient", y = "Best % Change from Baseline",
         fill = "Response") +
    theme_tnbc()
}

# =============================================================================
# plot_biomarker_ranking()
# Bar chart of top biomarkers by response probability (Fig. 5 style)
# =============================================================================
plot_biomarker_ranking <- function(bm_results, top_n = 20, metric = "response_prob") {
  df <- bm_results$by_resp_prob[1:min(top_n, nrow(bm_results$by_resp_prob)), ] %>%
    arrange(.data[[metric]]) %>%
    mutate(biomarker = factor(biomarker, levels = biomarker))

  ggplot(df, aes(x = biomarker, y = .data[[metric]])) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "red") +
    labs(title = paste("Top", top_n, "Biomarkers by", metric),
         x = "Biomarker", y = metric) +
    theme_tnbc()
}

# =============================================================================
# plot_response_rates()
# Bar chart of CR/PR/SD/PD rates vs KEYNOTE-119 (Fig. 3A style)
# =============================================================================
plot_response_rates <- function(sim_rates, clinical_rates = NULL) {
  df_sim <- data.frame(
    Response = c("CR","PR","SD","PD"),
    Rate     = c(sim_rates$CR, sim_rates$PR, sim_rates$SD, sim_rates$PD),
    Source   = "Simulation"
  )

  if (!is.null(clinical_rates)) {
    df_clin <- data.frame(
      Response = c("CR","PR","SD","PD"),
      Rate     = c(clinical_rates$CR, clinical_rates$PR, clinical_rates$SD, clinical_rates$PD),
      Source   = "KEYNOTE-119"
    )
    df <- rbind(df_sim, df_clin)
  } else {
    df <- df_sim
  }

  df$Response <- factor(df$Response, levels = c("CR","PR","SD","PD"))

  ggplot(df, aes(x = Response, y = Rate * 100, fill = Source)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c(Simulation = "steelblue", "KEYNOTE-119" = "coral")) +
    labs(title = "Response Rate Comparison",
         x = "Response", y = "Rate (%)", fill = "Source") +
    theme_tnbc()
}

# =============================================================================
# plot_spider()
# Spider plot of SLD over time per VP (Fig. 3 style)
# =============================================================================
plot_spider <- function(vp_ts_list, max_vps = 100) {
  n <- min(length(vp_ts_list), max_vps)
  df_list <- lapply(seq_len(n), function(i) {
    ts <- vp_ts_list[[i]]
    data.frame(
      time       = ts$time,
      pct_change = (ts$SLD - ts$SLD[1]) / ts$SLD[1] * 100,
      vp_id      = i,
      BOR        = attr(ts, "BOR")
    )
  })
  df <- do.call(rbind, df_list)
  df$BOR <- factor(df$BOR, levels = c("CR","PR","SD","PD"))

  ggplot(df, aes(x = time, y = pct_change, group = vp_id, color = BOR)) +
    geom_line(alpha = 0.5, size = 0.5) +
    geom_hline(yintercept = -30, linetype = "dashed", color = "blue") +
    geom_hline(yintercept =  20, linetype = "dashed", color = "red") +
    scale_color_manual(values = c(CR = "#2E86AB", PR = "#A0CFA8", SD = "#F6AE2D", PD = "#F26419")) +
    labs(title = "Spider Plot - SLD Change Over Time",
         x = "Time (days)", y = "% Change from Baseline SLD",
         color = "BOR") +
    ylim(-100, 100) +
    theme_tnbc()
}

# =============================================================================
# save_figure()
# Utility to save ggplot figures
# =============================================================================
save_figure <- function(plt, filename, width = 8, height = 6, dpi = 300) {
  dir.create("output/figures", showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path("output/figures", filename), plot = plt,
         width = width, height = height, dpi = dpi)
  message("Saved: output/figures/", filename)
}
