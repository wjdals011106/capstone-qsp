# =============================================================================
# app.R
# TNBC QSP Model — Shiny Application
# =============================================================================

library(shiny)
library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)
library(shinydashboard)

# Source model files (relative to app.R location)
# Works both from RStudio (Run App) and from Rscript
app_dir <- tryCatch(
  dirname(rstudioapi::getActiveDocumentContext()$path),
  error = function(e) getwd()
)
setwd(app_dir)
source("R/01_parameters.R")
source("R/02_model_ode.R")
source("R/03_events.R")
source("R/04_rules.R")
source("R/05_initial_conditions.R")
source("R/utils.R")
source("R/06_simulation.R")

# =============================================================================
# UI
# =============================================================================
ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "TNBC QSP Simulator"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("시뮬레이션 설정", tabName = "tab_setup",   icon = icon("sliders")),
      menuItem("종양 성장",       tabName = "tab_tumor",   icon = icon("chart-line")),
      menuItem("면역 반응",       tabName = "tab_immune",  icon = icon("shield-halved")),
      menuItem("RECIST 평가",     tabName = "tab_recist",  icon = icon("ruler")),
      menuItem("바이오마커",       tabName = "tab_biomarker", icon = icon("dna"))
    )
  ),

  dashboardBody(
    # Custom CSS
    tags$head(tags$style(HTML("
      .info-box { min-height: 70px; }
      .info-box-icon { height: 70px; line-height: 70px; }
      .info-box-content { padding-top: 8px; padding-bottom: 8px; }
      .small-box { min-height: 90px; }
    "))),

    tabItems(

      # ===========================================================
      # TAB 1: 시뮬레이션 설정
      # ===========================================================
      tabItem(tabName = "tab_setup",
        fluidRow(
          # ── 치료 설정 ──────────────────────────────────────────
          box(title = "치료 설정", status = "primary", solidHeader = TRUE, width = 4,
            numericInput("dose_pembro", "Pembrolizumab 용량 (mg)", value = 200, min = 0, max = 500, step = 10),
            numericInput("dose_interval", "투여 간격 (일)", value = 21, min = 7, max = 42, step = 7),
            numericInput("t_treat_end", "시뮬레이션 기간 (일)", value = 300, min = 50, max = 700, step = 50),
            hr(),
            checkboxInput("do_surgery", "원발 종양 수술", value = FALSE),
            conditionalPanel(
              condition = "input.do_surgery == true",
              numericInput("surgery_time", "수술 시점 (일)", value = 90, min = 0, max = 300, step = 10)
            )
          ),

          # ── 종양 생물학 ────────────────────────────────────────
          box(title = "종양 생물학", status = "warning", solidHeader = TRUE, width = 4,
            sliderInput("k_C1_growth", "원발 종양 성장률 (/day)", min = 0.001, max = 0.02,  value = 0.0065, step = 0.0005),
            sliderInput("k_C1_death", "기저 사멸률 (/day)",       min = 0.0001, max = 0.005, value = 0.001,  step = 0.0001),
            sliderInput("n_imm",      "신생항원 수",               min = 1, max = 8, value = 3, step = 1),
            hr(),
            checkboxInput("do_seeding", "전이 병소 시뮬레이션", value = TRUE),
            conditionalPanel(
              condition = "input.do_seeding == true",
              numericInput("delay_Ln1",   "폐 전이 1 시작 (일)", value = 100, min = 0, max = 400, step = 10),
              numericInput("delay_Ln2",   "폐 전이 2 시작 (일)", value = 150, min = 0, max = 400, step = 10),
              numericInput("delay_other", "기타 전이 시작 (일)",  value = 200, min = 0, max = 400, step = 10)
            )
          ),

          # ── 면역 환경 ──────────────────────────────────────────
          box(title = "면역 환경 파라미터", status = "success", solidHeader = TRUE, width = 4,
            sliderInput("k_T1_kill",   "CD8 T세포 살해율 (nM·/day)",  min = 0.1, max = 5,    value = 1.0,  step = 0.1),
            sliderInput("k_MDSC_sup",  "MDSC 억제 강도",              min = 0.1, max = 2,    value = 0.5,  step = 0.05),
            sliderInput("TGFbase",     "기저 TGF-β (nM)",             min = 0.01, max = 1.0, value = 0.1,  step = 0.01),
            sliderInput("IFNg_50_ind", "IFNγ EC50 (pM)",              min = 1,    max = 100,  value = 10,   step = 1),
            hr(),
            sliderInput("k_PD1_block", "PD-1 차단 효능 (fold)",       min = 0.5, max = 5.0,  value = 2.0,  step = 0.1)
          )
        ),

        fluidRow(
          column(12, align = "center",
            actionButton("btn_run", "▶  시뮬레이션 실행", class = "btn btn-primary btn-lg",
                         style = "margin: 10px; padding: 12px 40px; font-size: 16px;"),
            tags$br(),
            uiOutput("run_status")
          )
        ),

        fluidRow(
          # Summary info boxes
          infoBoxOutput("ibox_bor",      width = 3),
          infoBoxOutput("ibox_sld",      width = 3),
          infoBoxOutput("ibox_responder", width = 3),
          infoBoxOutput("ibox_time",     width = 3)
        )
      ),

      # ===========================================================
      # TAB 2: 종양 성장
      # ===========================================================
      tabItem(tabName = "tab_tumor",
        fluidRow(
          box(title = "종양 직경 변화 (RECIST)", status = "primary", solidHeader = TRUE, width = 8,
            plotOutput("plot_diameter", height = "380px")
          ),
          box(title = "표시 옵션", status = "info", width = 4,
            checkboxGroupInput("show_lesions", "병소 선택",
              choices  = list("원발 종양 (T)" = "d_T",
                              "폐 전이 1 (Ln1)" = "d_Ln1",
                              "폐 전이 2 (Ln2)" = "d_Ln2",
                              "기타 전이 (oth)" = "d_oth"),
              selected = c("d_T", "d_Ln1", "d_Ln2")
            ),
            hr(),
            checkboxInput("show_recist_lines", "RECIST 기준선 표시", value = TRUE),
            checkboxInput("show_dose_marks",   "투약 시점 표시",      value = TRUE),
            hr(),
            sliderInput("plot_t_max", "표시 기간 (일)", min = 50, max = 700, value = 300, step = 10)
          )
        ),
        fluidRow(
          box(title = "SLD (Sum of Longest Diameters)", status = "info", solidHeader = TRUE, width = 6,
            plotOutput("plot_sld", height = "280px")
          ),
          box(title = "종양 세포 클론 구성", status = "warning", solidHeader = TRUE, width = 6,
            plotOutput("plot_clones", height = "280px")
          )
        )
      ),

      # ===========================================================
      # TAB 3: 면역 반응
      # ===========================================================
      tabItem(tabName = "tab_immune",
        fluidRow(
          box(title = "종양 내 면역세포 동태", status = "success", solidHeader = TRUE, width = 8,
            plotOutput("plot_immune_T", height = "380px")
          ),
          box(title = "면역세포 선택", status = "info", width = 4,
            checkboxGroupInput("show_immune_cells", "세포 종류",
              choices = list(
                "CD8 T세포 (Tcyt)" = "Tcyt_total_T",
                "Treg (T0)" = "T0_T",
                "Th 세포" = "Th_T",
                "M1 대식세포" = "Mac_M1_T",
                "M2 대식세포" = "Mac_M2_T",
                "MDSC" = "MDSC_T"
              ),
              selected = c("Tcyt_total_T", "T0_T", "Mac_M1_T", "Mac_M2_T")
            ),
            hr(),
            radioButtons("immune_compartment", "구획 선택",
              choices = list("원발 종양" = "T", "폐 전이 1" = "Ln1", "폐 전이 2" = "Ln2"),
              selected = "T"
            )
          )
        ),
        fluidRow(
          box(title = "LN CD8 T세포 (항원별)", status = "primary", solidHeader = TRUE, width = 6,
            plotOutput("plot_Tcyt_LN", height = "280px")
          ),
          box(title = "종양 내 사이토카인", status = "warning", solidHeader = TRUE, width = 6,
            plotOutput("plot_cytokines", height = "280px")
          )
        )
      ),

      # ===========================================================
      # TAB 4: RECIST 평가
      # ===========================================================
      tabItem(tabName = "tab_recist",
        fluidRow(
          box(title = "RECIST v1.1 워터폴 플롯", status = "primary", solidHeader = TRUE, width = 8,
            plotOutput("plot_waterfall", height = "350px")
          ),
          box(title = "평가 결과", status = "info", width = 4,
            h4("Best Overall Response"),
            uiOutput("recist_bor_display"),
            hr(),
            h4("SLD 변화율 (% from baseline)"),
            tableOutput("recist_table")
          )
        ),
        fluidRow(
          box(title = "스파이더 플롯 (SLD 변화율 시계열)", status = "success", solidHeader = TRUE, width = 12,
            plotOutput("plot_spider", height = "280px")
          )
        )
      ),

      # ===========================================================
      # TAB 5: 바이오마커
      # ===========================================================
      tabItem(tabName = "tab_biomarker",
        fluidRow(
          box(title = "면역 분율 (Immune Fraction)", status = "primary", solidHeader = TRUE, width = 6,
            plotOutput("plot_immune_frac", height = "300px")
          ),
          box(title = "클론 다양성 (Shannon Index)", status = "warning", solidHeader = TRUE, width = 6,
            plotOutput("plot_diversity", height = "300px")
          )
        ),
        fluidRow(
          box(title = "M2/M1 비율 & PD-L1 발현", status = "danger", solidHeader = TRUE, width = 6,
            plotOutput("plot_m2m1", height = "280px")
          ),
          box(title = "바이오마커 요약 (t=0 vs 최종)", status = "success", solidHeader = TRUE, width = 6,
            tableOutput("bm_summary_table")
          )
        )
      )
    ) # end tabItems
  ) # end dashboardBody
) # end dashboardPage


# =============================================================================
# SERVER
# =============================================================================
server <- function(input, output, session) {

  # ── Reactive: run simulation on button click ──────────────────────────────
  sim_result <- eventReactive(input$btn_run, {
    showModal(modalDialog(
      title = "시뮬레이션 실행 중...",
      "ODE 솔버를 실행하고 있습니다. 약 1~3분 소요됩니다.",
      easyClose = FALSE, footer = NULL
    ))
    on.exit(removeModal())

    # Build parameter overrides from UI
    p_over <- list(
      dose_pembrolizumab = input$dose_pembro,
      dose_interval      = input$dose_interval,
      k_C1_growth        = input$k_C1_growth,
      k_C2_growth        = input$k_C1_growth,
      k_C3_growth        = input$k_C1_growth * 1.1,
      k_C4_growth        = input$k_C1_growth * 0.9,
      k_C5_growth        = input$k_C1_growth * 1.05,
      k_C1_death         = input$k_C1_death,
      k_T1_kill          = input$k_T1_kill,
      k_MDSC_sup         = input$k_MDSC_sup,
      TGFbase            = input$TGFbase,
      IFNg_50_ind        = input$IFNg_50_ind,
      do_surgery         = as.integer(input$do_surgery),
      surgery_time       = if (input$do_surgery) input$surgery_time else 9999,
      seeding            = as.integer(input$do_seeding),
      delay_Ln1          = input$delay_Ln1,
      delay_Ln2          = input$delay_Ln2,
      delay_other        = input$delay_other
    )

    result <- tryCatch(
      run_simulation(
        p_override    = p_over,
        t_treat_start = 0,
        t_treat_end   = input$t_treat_end,
        t_step        = 1,
        verbose       = FALSE
      ),
      error = function(e) {
        list(success = FALSE, error_msg = e$message)
      }
    )
    result
  })

  # ── Status text ───────────────────────────────────────────────────────────
  output$run_status <- renderUI({
    if (is.null(sim_result())) return(NULL)
    r <- sim_result()
    if (r$success) {
      tags$span(style = "color:green; font-weight:bold;", "✔ 시뮬레이션 완료")
    } else {
      tags$span(style = "color:red; font-weight:bold;",
                paste("✘ 실패:", r$error_msg))
    }
  })

  # ── Info boxes ────────────────────────────────────────────────────────────
  output$ibox_bor <- renderInfoBox({
    r <- sim_result()
    if (is.null(r) || !r$success) {
      infoBox("BOR", "—", icon = icon("question"), color = "gray", fill = TRUE)
    } else {
      col <- switch(r$BOR, CR = "green", PR = "light-blue", SD = "yellow", PD = "red", "gray")
      infoBox("Best Overall Response", r$BOR, icon = icon("star"), color = col, fill = TRUE)
    }
  })

  output$ibox_sld <- renderInfoBox({
    r <- sim_result()
    if (is.null(r) || !r$success) {
      infoBox("SLD 감소", "—", icon = icon("ruler"), color = "gray", fill = TRUE)
    } else {
      bl  <- r$SLD_baseline
      mn  <- min(r$SLD_series$SLD, na.rm = TRUE)
      pct <- round((mn - bl) / bl * 100, 1)
      col <- if (pct < -30) "light-blue" else if (pct < 20) "yellow" else "red"
      infoBox("최대 SLD 변화", paste0(pct, "%"), icon = icon("ruler"), color = col, fill = TRUE)
    }
  })

  output$ibox_responder <- renderInfoBox({
    r <- sim_result()
    if (is.null(r) || !r$success) {
      infoBox("반응 여부", "—", icon = icon("user"), color = "gray", fill = TRUE)
    } else {
      col <- if (r$responder) "green" else "red"
      lbl <- if (r$responder) "Responder" else "Non-Responder"
      infoBox("반응 여부", lbl, icon = icon("user"), color = col, fill = TRUE)
    }
  })

  output$ibox_time <- renderInfoBox({
    r <- sim_result()
    if (is.null(r) || !r$success) {
      infoBox("기간", "—", icon = icon("clock"), color = "gray", fill = TRUE)
    } else {
      n  <- nrow(r$timeseries)
      mx <- max(r$timeseries$time, na.rm = TRUE)
      infoBox("시뮬레이션 기간", paste0(round(mx), "일 / ", n, "행"),
              icon = icon("clock"), color = "blue", fill = TRUE)
    }
  })

  # ── Helper: get timeseries safely ─────────────────────────────────────────
  ts_data <- reactive({
    r <- sim_result()
    if (is.null(r) || !r$success) return(NULL)
    r$timeseries
  })

  raw_data <- reactive({
    r <- sim_result()
    if (is.null(r) || !r$success) return(NULL)
    # rebuild full state timeseries from ODE output columns
    r$timeseries
  })

  # ==========================================================================
  # TAB 2: 종양 성장
  # ==========================================================================

  output$plot_diameter <- renderPlot({
    ts  <- ts_data()
    if (is.null(ts)) return(plot_placeholder("시뮬레이션을 먼저 실행하세요"))

    t_max  <- input$plot_t_max
    chosen <- input$show_lesions

    # Build long-form data
    cols   <- intersect(chosen, colnames(ts))
    if (length(cols) == 0) return(plot_placeholder("병소를 선택하세요"))

    labels <- c(d_T = "원발 종양", d_Ln1 = "폐 전이 1", d_Ln2 = "폐 전이 2", d_oth = "기타 전이")
    df_long <- ts[ts$time <= t_max, c("time", cols)] %>%
      pivot_longer(-time, names_to = "lesion", values_to = "diameter") %>%
      mutate(lesion = recode(lesion, !!!labels))

    p <- ggplot(df_long, aes(x = time, y = diameter, color = lesion)) +
      geom_line(linewidth = 1.1) +
      labs(title = "종양 직경 변화", x = "시간 (일)", y = "직경 (cm)", color = "병소") +
      theme_bw(base_size = 13) +
      theme(legend.position = "top")

    if (input$show_recist_lines) {
      r    <- sim_result()
      bl_T <- if (!is.null(r) && r$success) r$SLD_baseline else NA
      if (!is.na(bl_T)) {
        p <- p +
          geom_hline(yintercept = bl_T * 0.70, linetype = "dashed", color = "steelblue", alpha = 0.7) +
          geom_hline(yintercept = bl_T * 1.20, linetype = "dashed", color = "tomato",    alpha = 0.7) +
          annotate("text", x = t_max * 0.85, y = bl_T * 0.70, label = "PR (-30%)",
                   color = "steelblue", size = 3.5, vjust = -0.5) +
          annotate("text", x = t_max * 0.85, y = bl_T * 1.20, label = "PD (+20%)",
                   color = "tomato",    size = 3.5, vjust = -0.5)
      }
    }

    if (input$show_dose_marks) {
      p_base <- define_parameters()
      dose_t <- seq(0, t_max, by = p_base$dose_interval)
      p <- p + geom_vline(xintercept = dose_t, linetype = "dotted", color = "gray50", alpha = 0.4)
    }
    p
  })

  output$plot_sld <- renderPlot({
    r  <- sim_result()
    if (is.null(r) || !r$success) return(plot_placeholder())
    ts <- ts_data()
    t_max <- input$plot_t_max

    df <- data.frame(
      time = r$SLD_series$time,
      pct  = (r$SLD_series$SLD - r$SLD_baseline) / r$SLD_baseline * 100
    )
    df <- df[df$time <= t_max, ]

    ggplot(df, aes(x = time, y = pct)) +
      geom_hline(yintercept = c(-30, 20), linetype = "dashed",
                 color = c("steelblue", "tomato")) +
      geom_hline(yintercept = 0, color = "gray40") +
      geom_line(color = "darkgreen", linewidth = 1.2) +
      geom_point(size = 3, color = "darkgreen") +
      labs(title = "SLD % 변화 (from baseline)",
           x = "시간 (일)", y = "SLD 변화율 (%)") +
      theme_bw(base_size = 13)
  })

  output$plot_clones <- renderPlot({
    ts <- ts_data()
    if (is.null(ts)) return(plot_placeholder())
    t_max <- input$plot_t_max

    clone_cols <- c("C1_T","C2_T","C3_T","C4_T","C5_T")
    avail <- intersect(clone_cols, colnames(ts))
    if (length(avail) == 0) return(plot_placeholder("클론 정보 없음"))

    df <- ts[ts$time <= t_max, c("time", avail)] %>%
      pivot_longer(-time, names_to = "clone", values_to = "cells") %>%
      group_by(time) %>%
      mutate(frac = cells / sum(cells + 1e-10))

    ggplot(df, aes(x = time, y = frac, fill = clone)) +
      geom_area(alpha = 0.85) +
      scale_fill_brewer(palette = "Set2") +
      labs(title = "원발 종양 클론 구성 비율",
           x = "시간 (일)", y = "클론 비율", fill = "클론") +
      theme_bw(base_size = 13) +
      theme(legend.position = "top")
  })

  # ==========================================================================
  # TAB 3: 면역 반응
  # ==========================================================================

  output$plot_immune_T <- renderPlot({
    ts <- ts_data()
    if (is.null(ts)) return(plot_placeholder())

    comp   <- input$immune_compartment
    chosen <- input$show_immune_cells
    t_max  <- input$plot_t_max

    # Map column names to compartment
    col_map <- list(
      Tcyt_total_T = paste0("Tcyt_total_", comp),
      T0_T         = paste0("T0_", comp),
      Th_T         = paste0("Th_", comp),
      Mac_M1_T     = paste0("Mac_M1_", comp),
      Mac_M2_T     = paste0("Mac_M2_", comp),
      MDSC_T       = paste0("MDSC_", comp)
    )
    label_map <- c(
      Tcyt_total_T = "CD8 Tcyt",
      T0_T         = "Treg (T0)",
      Th_T         = "Th 세포",
      Mac_M1_T     = "M1 대식세포",
      Mac_M2_T     = "M2 대식세포",
      MDSC_T       = "MDSC"
    )

    cols_actual <- unlist(col_map[chosen])
    cols_actual <- intersect(cols_actual, colnames(ts))
    if (length(cols_actual) == 0) return(plot_placeholder("데이터 없음"))

    df <- ts[ts$time <= t_max, c("time", cols_actual)]
    colnames(df)[-1] <- label_map[chosen[chosen %in% names(col_map)]]

    df_long <- pivot_longer(df, -time, names_to = "cell", values_to = "count")

    ggplot(df_long, aes(x = time, y = count, color = cell)) +
      geom_line(linewidth = 1.0) +
      scale_y_log10(labels = scales::comma) +
      labs(title = paste0("면역세포 동태 — ", comp, " 구획"),
           x = "시간 (일)", y = "세포 수 (log)", color = "세포 종류") +
      theme_bw(base_size = 13) +
      theme(legend.position = "right")
  })

  output$plot_Tcyt_LN <- renderPlot({
    ts <- ts_data()
    if (is.null(ts)) return(plot_placeholder())

    t_max <- input$plot_t_max
    cols  <- paste0("T", 1:8, "_LN")
    avail <- intersect(cols, colnames(ts))
    if (length(avail) == 0) return(plot_placeholder("LN T세포 없음"))

    df <- ts[ts$time <= t_max, c("time", avail)] %>%
      pivot_longer(-time, names_to = "spec", values_to = "count") %>%
      mutate(spec = sub("_LN", "", spec))

    ggplot(df, aes(x = time, y = count, color = spec)) +
      geom_line(linewidth = 0.9) +
      labs(title = "LN CD8 T세포 (항원 특이성별)",
           x = "시간 (일)", y = "세포 수", color = "항원") +
      theme_bw(base_size = 13) +
      theme(legend.position = "top")
  })

  output$plot_cytokines <- renderPlot({
    ts <- ts_data()
    if (is.null(ts)) return(plot_placeholder())

    t_max  <- input$plot_t_max
    cyto_cols <- intersect(c("IFNg_T","IL2_T","TGFb_T","IL10_T"), colnames(ts))
    if (length(cyto_cols) == 0) return(plot_placeholder("사이토카인 없음"))

    df <- ts[ts$time <= t_max, c("time", cyto_cols)] %>%
      pivot_longer(-time, names_to = "cyto", values_to = "conc")

    ggplot(df, aes(x = time, y = conc, color = cyto)) +
      geom_line(linewidth = 0.9) +
      labs(title = "종양 내 사이토카인 농도",
           x = "시간 (일)", y = "농도 (nM)", color = "사이토카인") +
      theme_bw(base_size = 13) +
      theme(legend.position = "top")
  })

  # ==========================================================================
  # TAB 4: RECIST 평가
  # ==========================================================================

  output$plot_waterfall <- renderPlot({
    r <- sim_result()
    if (is.null(r) || !r$success) return(plot_placeholder())

    df <- data.frame(
      eval = paste0("평가 ", seq_along(r$SLD_series$SLD)),
      pct  = (r$SLD_series$SLD - r$SLD_baseline) / r$SLD_baseline * 100
    )
    df$col <- ifelse(df$pct <= -30, "PR/CR",
               ifelse(df$pct >= 20, "PD", "SD"))

    ggplot(df, aes(x = eval, y = pct, fill = col)) +
      geom_col(width = 0.6) +
      geom_hline(yintercept = c(-30, 20), linetype = "dashed",
                 color = c("steelblue","tomato"), linewidth = 0.8) +
      scale_fill_manual(values = c("PR/CR" = "#2196F3", "SD" = "#FFC107", "PD" = "#F44336")) +
      labs(title = "RECIST 워터폴 플롯",
           x = "평가 시점 (9주 간격)", y = "SLD 변화율 (%)", fill = "반응") +
      theme_bw(base_size = 13)
  })

  output$recist_bor_display <- renderUI({
    r <- sim_result()
    if (is.null(r) || !r$success) return(tags$p("—"))
    col_map <- c(CR = "#00C853", PR = "#2196F3", SD = "#FFC107", PD = "#F44336")
    col <- col_map[r$BOR]
    tags$div(
      tags$span(r$BOR,
                style = paste0("font-size:36px; font-weight:bold; color:", col, ";")),
      tags$p(if (r$responder) "✔ 반응자 (Responder)" else "✘ 비반응자 (Non-Responder)")
    )
  })

  output$recist_table <- renderTable({
    r <- sim_result()
    if (is.null(r) || !r$success) return(NULL)
    data.frame(
      `평가 시점 (일)` = round(r$SLD_series$time, 1),
      `SLD (cm)`       = round(r$SLD_series$SLD, 3),
      `변화율 (%)`     = round((r$SLD_series$SLD - r$SLD_baseline) / r$SLD_baseline * 100, 1),
      `RECIST`         = sapply(r$SLD_series$SLD, function(s)
                           classify_recist(r$SLD_baseline, s)),
      check.names = FALSE
    )
  }, bordered = TRUE, striped = TRUE)

  output$plot_spider <- renderPlot({
    r <- sim_result()
    if (is.null(r) || !r$success) return(plot_placeholder())

    df <- data.frame(
      time = r$SLD_series$time,
      pct  = (r$SLD_series$SLD - r$SLD_baseline) / r$SLD_baseline * 100,
      patient = "Patient 1"
    )

    ggplot(df, aes(x = time, y = pct, group = patient, color = patient)) +
      geom_line(linewidth = 1.2) +
      geom_point(size = 3) +
      geom_hline(yintercept = c(-30, 20, 0), linetype = "dashed",
                 color = c("steelblue","tomato","gray50"), alpha = 0.7) +
      labs(title = "스파이더 플롯 — SLD 변화율",
           x = "시간 (일)", y = "SLD 변화율 (%)") +
      theme_bw(base_size = 13) +
      theme(legend.position = "none")
  })

  # ==========================================================================
  # TAB 5: 바이오마커
  # ==========================================================================

  output$plot_immune_frac <- renderPlot({
    ts <- ts_data()
    if (is.null(ts)) return(plot_placeholder())

    t_max <- input$plot_t_max
    cols  <- intersect(c("immune_frac_T","Tcyt_frac_T","Treg_frac_LN"), colnames(ts))
    if (length(cols) == 0) return(plot_placeholder("바이오마커 없음"))

    labels <- c(immune_frac_T = "면역세포 분율 (종양)",
                Tcyt_frac_T  = "CD8 Tcyt 분율 (종양 T)",
                Treg_frac_LN = "Treg 분율 (LN)")

    df <- ts[ts$time <= t_max, c("time", cols)] %>%
      pivot_longer(-time, names_to = "bm", values_to = "val") %>%
      mutate(bm = recode(bm, !!!labels))

    ggplot(df, aes(x = time, y = val, color = bm)) +
      geom_line(linewidth = 1.0) +
      scale_y_continuous(limits = c(0, 1), labels = scales::percent) +
      labs(title = "면역세포 분율 변화",
           x = "시간 (일)", y = "분율", color = "바이오마커") +
      theme_bw(base_size = 13) +
      theme(legend.position = "top")
  })

  output$plot_diversity <- renderPlot({
    ts <- ts_data()
    if (is.null(ts)) return(plot_placeholder())

    t_max <- input$plot_t_max
    cols  <- intersect(c("clone_shannon_T","Tcyt_shannon_LN","clone_evenness_T"), colnames(ts))
    if (length(cols) == 0) return(plot_placeholder("다양성 지수 없음"))

    labels <- c(clone_shannon_T  = "암 클론 Shannon (종양)",
                Tcyt_shannon_LN  = "Tcyt Shannon (LN)",
                clone_evenness_T = "암 클론 Evenness (종양)")

    df <- ts[ts$time <= t_max, c("time", cols)] %>%
      pivot_longer(-time, names_to = "idx", values_to = "val") %>%
      mutate(idx = recode(idx, !!!labels))

    ggplot(df, aes(x = time, y = val, color = idx)) +
      geom_line(linewidth = 1.0) +
      labs(title = "클론/면역 다양성 지수",
           x = "시간 (일)", y = "Shannon Index / Evenness", color = "지수") +
      theme_bw(base_size = 13) +
      theme(legend.position = "top")
  })

  output$plot_m2m1 <- renderPlot({
    ts <- ts_data()
    if (is.null(ts)) return(plot_placeholder())

    t_max <- input$plot_t_max
    cols  <- intersect(c("M2M1_T","PDL1_tumor"), colnames(ts))
    if (length(cols) == 0) return(plot_placeholder("없음"))

    labels <- c(M2M1_T = "M2/M1 비율 (종양)", PDL1_tumor = "PD-L1 발현")
    df <- ts[ts$time <= t_max, c("time", cols)] %>%
      pivot_longer(-time, names_to = "bm", values_to = "val") %>%
      mutate(bm = recode(bm, !!!labels))

    ggplot(df, aes(x = time, y = val, color = bm)) +
      geom_line(linewidth = 1.0) +
      facet_wrap(~bm, scales = "free_y") +
      labs(title = "M2/M1 비율 & PD-L1 발현",
           x = "시간 (일)", y = "값") +
      theme_bw(base_size = 13) +
      theme(legend.position = "none")
  })

  output$bm_summary_table <- renderTable({
    ts <- ts_data()
    if (is.null(ts)) return(NULL)

    bm_cols <- c("immune_frac_T","Tcyt_frac_T","Treg_frac_LN",
                 "clone_shannon_T","Tcyt_shannon_LN","M2M1_T","PDL1_tumor")
    avail <- intersect(bm_cols, colnames(ts))
    if (length(avail) == 0) return(NULL)

    row0   <- as.numeric(ts[1,                avail])
    row_end <- as.numeric(ts[nrow(ts),        avail])
    labels  <- c(immune_frac_T = "면역세포 분율",
                 Tcyt_frac_T   = "CD8 Tcyt 분율",
                 Treg_frac_LN  = "Treg 분율 LN",
                 clone_shannon_T = "클론 Shannon",
                 Tcyt_shannon_LN = "Tcyt Shannon",
                 M2M1_T         = "M2/M1 비율",
                 PDL1_tumor      = "PD-L1 발현")

    data.frame(
      바이오마커 = labels[avail],
      `t=0`      = round(row0,    4),
      `최종`     = round(row_end, 4),
      `변화`     = round(row_end - row0, 4),
      check.names = FALSE
    )
  }, bordered = TRUE, striped = TRUE)

} # end server


# =============================================================================
# Helper: placeholder plot when no data
# =============================================================================
plot_placeholder <- function(msg = "시뮬레이션을 실행하세요") {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = msg, size = 7, color = "gray50") +
    theme_void()
}

# =============================================================================
# RUN APP
# =============================================================================
shinyApp(ui = ui, server = server)
