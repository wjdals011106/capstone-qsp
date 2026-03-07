# QSP Model Reproduction Plan
## Metastatic TNBC - Pembrolizumab Monotherapy
### Based on: Arulraj et al., Sci. Adv. 9, eadg0289 (2023)

---

## 1. Project Overview

**Goal:** Arulraj et al. (2023)의 전이성 삼중음성유방암(mTNBC) QSP 모델을 R로 재현한다.
- 원 논문은 MATLAB SimBiology로 구현됨
- R의 `deSolve`, `mrgsolve` 등의 ODE solver 패키지를 활용하여 재현

**Model Scale:**
- 746 species (상태변수)
- 737 parameters
- 1,428 reactions (ODE 정의)
- 314 algebraic rules
- 65 discrete events
- 127 virtual patient parameter distributions

---

## 2. Model Architecture (Compartments)

### 2.1 Compartment Structure (9개 주요 구획)

| Compartment | Description | Initial Volume |
|-------------|-------------|----------------|
| V_C | Central (Blood) | 5 L |
| V_P | Peripheral | 60 L |
| V_T | Primary Tumor | Dynamic (starts ~0.000001 mL) |
| V_T_Ln1 | Lung Metastatic Tumor 1 | Dynamic |
| V_T_Ln2 | Lung Metastatic Tumor 2 | Dynamic |
| V_T_other | Other Metastatic Tumor | Dynamic |
| V_LN | Primary Tumor-draining LN | 1112.6 mm³ |
| V_LN_Ln | Lung Tumor-draining LN | 1112.6 mm³ |
| V_LN_other | Other Tumor-draining LN | 1112.6 mm³ |

### 2.2 Compartment Connectivity
```
Central (Blood) <---> Peripheral
    |
    |---> Primary Tumor <---> Primary LN
    |---> Lung Met Tumor 1 <---> Lung LN
    |---> Lung Met Tumor 2 <---> Lung LN
    |---> Other Met Tumor <---> Other LN
```

---

## 3. Model Species (주요 세포 및 분자)

### 3.1 Cell Types
| Category | Species | Description |
|----------|---------|-------------|
| Cancer | C1-C5 | 5개 cancer clones (각기 다른 성장률, neo-epitope 발현) |
| T cells | nT0 (naive CD4), nT1 (naive CD8) | Naive T cells |
| | T0 (Treg), Th (Helper T) | CD4+ T cell subsets |
| | T1-T8 | Neo-antigen specific cytotoxic T cells (8 specificities) |
| | T1_exh-T8_exh | Exhausted cytotoxic T cells |
| APCs | APC, mAPC | Immature/Mature antigen-presenting cells |
| Macrophages | M1, M2 | M1 (antitumor) / M2 (protumor) macrophages |
| MDSCs | MDSC | Myeloid-derived suppressor cells |

Initial amounts (Central compartment):
- `nT0` (naive CD4) = 3706.9 cells
- `nT1` (naive CD8) = 2274.8 cells
- `nT0` (peripheral) = 185,344.8 cells

### 3.2 Cytokines & Soluble Factors
| Species | Role |
|---------|------|
| IL-2 | T cell proliferation |
| IL-10 | Immunosuppressive |
| IL-12 | APC maturation, M2->M1 polarization |
| IFN-gamma | Antitumor, PD-L1 upregulation |
| TGF-beta | Immunosuppressive, Treg induction |
| CCL2 | MDSC/Macrophage recruitment |
| Angiogenic factors | Tumor vasculature growth |

### 3.3 Checkpoint Molecules
| Molecule | Expression |
|----------|-----------|
| PD-1 | Cytotoxic T cells, Macrophages |
| PD-L1 | Cancer cells (upregulated by IFN-gamma) |
| PD-L2 | Cancer cells |
| CD28 | T cells |
| CTLA-4 | T cells |
| CD80/CD86 | APCs |
| CD47/SIRPalpha | Cancer cells / Macrophages |

---

## 4. Core Model Equations

### 4.1 Cancer Cell Dynamics (Modified Gompertzian Growth)

클론 i (i = 1~5)에 대해:
```
dC_i/dt = k_C_growth_i * C_i * ln(K / C_total)
          - k_C_death * C_i
          - k_C_T * f(T_cyt, C) * C_i
          - k_C_M1 * f(M1, C) * C_i
```
- C_i: i번째 cancer clone
- K: carrying capacity (종양 혈관에 의해 동적으로 변화)
- k_C_T: cytotoxic T cell에 의한 cancer cell killing rate
- k_C_M1: M1 macrophage에 의한 phagocytosis rate

### 4.2 Tumor Carrying Capacity (Angiogenesis)
```
dK/dt = k_K_g * AF / (AF + c_vas_50) - k_K_d * K
d(AF)/dt = k_vas_Csec * C_total - k_vas_deg * AF
```
- AF: angiogenic factor concentration
- Angiogenic factors secreted by cancer cells and M2 macrophages

### 4.3 Tumor Volume (Algebraic Rule)
```
V_T = V_Tmin + (C_total * vol_cell + T_total * vol_Tcell) / Ve_T
              + M_total * vol_Mcell / Ve_T
```

### 4.4 Naive T Cell Dynamics
```
dnT/dt = k_nT_source - k_nT_death * nT
         - k_nT_mig * (nT - nT_peripheral)
         - k_nT_LN * nT   # trafficking to LN
```

### 4.5 T Cell Activation (in Lymph Nodes)
TCR-pMHC ligation 기반 활성화:
```
activation_rate = f_Hill(pMHC_on_mAPC) = pMHC / (p0_50 + pMHC)
```
Number of divisions (linear sum):
```
N_div = N0 * TCR_signal
      + N_costim * CD28_signal / (1 + CTLA4_competition)
      + N_IL2 * f(IL2)
```
- N0 = 2, N_costim = 3
- N_IL2_CD8 = 11, N_IL2_CD4 = 8.5

### 4.6 Cytotoxic T Cell Dynamics (in Tumor)
침윤:
```
infiltration = k_T_infiltration * V_T * f_vas * T_central
```
암세포 살상 (PD-1, TGF-β, MDSC 억제 포함):
```
k_eff = k_C_T1 * (1 - f_PD1) * (1 - f_TGFb) * (1 - f_MDSC)
killing = k_eff * Tcyt / (Tcyt + C_total * K_T_C) * C_i
```
소진(Exhaustion):
```
exhaustion_rate = k_T1 * f(PD1-PDL1) + k_T_IL10 * f(IL10)
```

### 4.7 Macrophage Dynamics
Recruitment (CCL2 의존):
```
recruitment = k_Mac_recruit * CCL2 / (CCL2 + CCL2_50)
```
Polarization:
- M1 -> M2: promoted by IL-10, TGF-beta
- M2 -> M1: promoted by IL-12, IFN-gamma

Phagocytosis: inhibited by IL-10 and CD47-SIRPalpha interaction

### 4.8 MDSC Dynamics
- Recruitment: CCL2-dependent
- Effects: Release Arg-I and NO -> inhibit cytotoxic T cell activity

### 4.9 Cytokine Dynamics
```
dCytokine/dt = secretion_terms - k_deg * Cytokine
```

### 4.10 PD-L1 Expression & Checkpoint
PD-L1 발현 (IFN-γ 유도):
```
PD-L1 = PD-L1_base * (1 + r_PDL1_IFNg * IFNg / (IFNg + IFNg_50))
```
PD-1/PD-L1 binding → T cell killing 및 macrophage phagocytosis 억제 (Hill function)

### 4.11 Pembrolizumab PK/PD
PK (2-compartment model):
```
dA_central/dt = -CL * A_central/V_C - Q * (A_central/V_C - A_peripheral/V_P) + dose_rate
dA_peripheral/dt = Q * (A_central/V_C - A_peripheral/V_P)
```
- Dose: 200 mg IV every 3 weeks (KEYNOTE-119)

PD:
- Anti-PD-1이 T cells 및 macrophage의 PD-1에 결합
- PD-1/PD-L1, PD-1/PD-L2 상호작용 차단
- 억제 신호 감소 → killing/phagocytosis 강화

### 4.12 Biomarker Diversity Indices
Shannon Index:
```
H = -sum(p_k * ln(p_k))    (k = 1..S)
```
Evenness:
```
J = H / ln(S)
```
Responder Inclusion Score (RIS):
```
RIS = [Responders in subgroup / Total responders]
    - [Non-responders in subgroup / Total non-responders]
```
범위: -1 ~ +1

---

## 5. Key Parameters Summary

### 5.1 Physical / Physiological Parameters
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| V_C | 5 | L | Central compartment volume |
| V_P | 60 | L | Peripheral compartment volume |
| V_LN | 1112.6 | mm³ | Lymph node volume |
| vol_cell | 2572.44 | μm³/cell | Cancer cell volume |
| vol_Tcell | 175.02 | μm³/cell | T cell volume |
| Ve_T | 0.37 | - | Tumor void fraction |
| k_cell_clear | 0.1 | 1/day | Dead cell clearance rate |

### 5.2 Cancer Cell Growth Parameters
| Parameter | Value | Unit | Source |
|-----------|-------|------|--------|
| k_C1_growth | 0.0065 | 1/day | Ryu 2014; Desai 2006 |
| k_C1_death | 0.0001 | 1/day | Palsson 2013 |
| k_K_g | 4.12 | 1/day | Desai 2006 |
| k_K_d | 0.0034 | 1/day | Desai 2006; Hahnfeldt 1999 |
| k_vas_Csec | 0.00011 | pg/cell/day | Volk 2008 |
| k_vas_deg | 16.6 | 1/day | Finley 2011 |
| c_vas_50 | 1070 | pg/mL | Desai 2006 |
| C_max | 27000 | cells | Desai 2006 |
| initial_met_diameter | 1.6 | cm | Huang 2021 |

### 5.3 T Cell Parameters
| Parameter | Value | Unit | Source |
|-----------|-------|------|--------|
| div_T0 | 1,160,000 | - | Robins 2009 |
| Q_nT0_thym | 70,000,000 | cell/day | Bains 2009 |
| Q_nT1_thym | 35,000,000 | cell/day | Bains 2009 |
| k_T0_act | 5 | 1/day | De Boer & Perelson 1995 |
| k_T1_act | 23 | 1/day | De Boer & Perelson 1995 |
| k_C_T1 | 0.95 | 1/day | Robertson-Tessi 2012 |
| k_T1 | 0.1 | 1/day | 추정 |
| k_Treg | 0.05 | 1/day | 추정 |
| K_T_C | 1.2 | - | Robertson-Tessi 2012 |
| N0 | 2 | - | Marchingo 2014 |
| N_costim | 3 | - | Marchingo 2014 |
| N_IL2_CD8 | 11 | - | Marchingo 2014 |
| N_IL2_CD4 | 8.5 | - | Marchingo 2014 |
| IL2_50 | 0.32 | nM | Marchingo 2014 |
| k_IL2_deg | 0.2 | 1/min | Lotze 1985 |
| k_IL2_sec | 3e-5 | nmol/cell/hr | Han 2012 |
| k_IL2_cons | 6e-6 | nmol/cell/hr | Lotze 1985 |

### 5.4 APC Parameters
| Parameter | Value | Unit | Source |
|-----------|-------|------|--------|
| k_APC_mat | 1.5 | 1/day | Chen 2014 |
| k_APC_mig | 4 | 1/day | Russo 2016 |
| k_APC_death | 0.01 | 1/day | Marino & Kirschner 2004 |
| APC0_T | 400,000 | cells/mL | Lavin 2017 |
| APC0_LN | 1,200,000 | cells/mL | Catron 2004 |
| p0_50 | 2.6455e-5 | molecules/μm² | Kimachi 1997 |
| DAMPs | 1.34e-14 | mol/cell | Milo 2013 |

### 5.5 Macrophage Parameters
| Parameter | Value | Unit | Description |
|-----------|-------|------|-------------|
| k_Mac_death | 0.01 | 1/day | Macrophage death rate |
| CCL2_50 | (sampled) | nM | Half-max CCL2 for MDSC recruitment |

### 5.6 Virtual Patient Parameters (127 varied parameters)
| Parameter | Distribution | Median (SD) | Unit |
|-----------|-------------|-------------|------|
| initial_met_diameter | Log-normal | 1.65 (0.3) | cm |
| agconc_Ci_j (neo-antigen) | Log-normal | 6.75e-14 (1) | mol/cell |
| seeding_time | Sampled | Variable | day |
| div_T0 | Log-normal | 1,160,000 | - |
| k_C_growth (clone별) | Log-normal | - | 1/day |

### 5.7 Sensitive Parameters (ORR에 영향이 큰 파라미터들)
- `k_C_T1`: Cancer cell killing rate by cytotoxic T cells (제1위)
- `CCL2_50`: Half-maximal CCL2 for MDSC recruitment (증가 시 ORR 증가)
- `div_T0`: T cell diversity (증가 시 ORR 증가)
- `k_CCL2_deg`: CCL2 degradation rate (증가 시 ORR 증가)
- `k_Treg`: T cell death rate by Tregs (증가 시 ORR 감소)
- `agconc_C_0`: Self-antigen concentration (증가 시 ORR 감소)

---

## 6. Model Calibration Targets

### 6.1 KEYNOTE-119 임상시험 보정 목표
| Response | KEYNOTE-119 | 시뮬레이션 목표 |
|----------|-------------|----------------|
| CR | 3.97% | ~2% |
| PR | 6.85% | ~9% |
| SD | 22.38% | ~22% |
| PD | 66.78% | ~67% |
| 폐 전이 환자 비율 | 65% | 65% |

### 6.2 Baseline Species 범위 (Fig. S3/S4 기준)
| Species | 기저값 (중앙값) | 단위 |
|---------|---------------|------|
| T cells in central | 1~2 × 10⁴ | cells/mL |
| Tcyt in tumor | 5~10 × 10⁷ | cells/mL |
| Treg in tumor | 2~5 × 10⁷ | cells/mL |
| M2/M1 Mac ratio | 5~25 | - |
| MDSC in tumor | 2~5 × 10⁵ | cells/mL |
| IL-2 in LN | ~0.47 | nM |
| TGF-β in tumor | 0.15~0.30 | nM |
| IFN-γ in tumor | 1~1.5 × 10⁻³ | nM |
| Tumor diameter | 1.5~4.5 | cm |
| PD-L1 in tumor | 260~380 | molecules/μm² |
| Carrying capacity | 2~4 × 10⁹ | cells |

---

## 7. RECIST v1.1 Response Evaluation

```
# 종양 직경 합 (최대 2개 폐 전이 + 1개 기타 전이)
SLD_baseline = sum(diameter_Ln1, diameter_Ln2, diameter_other)

# 반응 분류
CR:  모든 종양 소실
PR:  SLD_t / SLD_baseline <= 0.70  (30% 이상 감소)
PD:  SLD_t / SLD_baseline >= 1.20  (20% 이상 증가)
SD:  그 외

# Responder = CR/PR 또는 SD (최소 24주 지속)
# Non-responder = PD
```

---

## 8. Top Predictive Biomarkers (Data_S2 기준)

### 8.1 단일 바이오마커 (반응 확률 순)
| Rank | Biomarker | Cutoff | Direction | Response Prob. | RIS |
|------|-----------|--------|-----------|----------------|-----|
| 1 | APC density in LNs | 233 cells/mm³ | > | 0.857 | 0.101 |
| 2 | Cancer clone richness | 2.714 | < | 0.750 | 0.046 |
| 3 | Tcyt/T cells in LNs | 0.924 | > | 0.741 | 0.182 |
| 4 | Tcyt/T cells in central | 0.926 | > | 0.691 | 0.160 |
| 5 | Exhausted Tcyt fraction | 0.205 | > | 0.680 | 0.048 |

RIS 기준 1위: Treg density in LNs (RIS = 0.351)

### 8.2 최우수 복합 바이오마커
| Combination | Response Prob. |
|-------------|----------------|
| Tcyt/T cells in central + APCs in LNs | **1.000** |
| Tcyt/T cells in LNs + APCs in LNs | **1.000** |
| Immune cell fraction + Treg density in LNs | RIS = **0.43** |

---

## 9. Implementation Plan (R Code)

### Directory Structure
```
tnbc/
├── plan.md                    # 이 파일
├── R/
│   ├── 01_parameters.R        # 파라미터 정의 (737개)
│   ├── 02_model_ode.R         # ODE 시스템 정의
│   ├── 03_events.R            # Discrete events (dosing, seeding)
│   ├── 04_rules.R             # Algebraic rules
│   ├── 05_initial_conditions.R # 초기조건 설정
│   ├── 06_simulation.R        # 시뮬레이션 실행
│   ├── 07_virtual_patients.R  # Virtual patient 생성
│   ├── 08_biomarker_analysis.R # 바이오마커 분석
│   ├── 09_visualization.R     # 결과 시각화
│   └── utils.R                # 유틸리티 함수
├── data/
│   ├── adg0289_Data_S1.xlsx   # 모델 파라미터 원본
│   └── adg0289_Data_S2.xlsx   # 바이오마커 분석 원본
├── output/
│   ├── figures/
│   └── results/
└── main.R                     # 메인 실행 스크립트
```

### Required R Packages
```r
deSolve    # ODE solver (lsoda - stiff system 대응)
mrgsolve   # C++ 기반 대규모 ODE (대안, 권장)
readxl     # Excel 파일 읽기
lhs        # Latin Hypercube Sampling (virtual patients)
ggplot2    # 시각화
dplyr, tidyr  # 데이터 처리
parallel   # 병렬 연산 (VP 시뮬레이션)
```

### Phase-by-Phase Implementation

| Phase | Task | Week | Deliverable |
|-------|------|------|-------------|
| 1 | Infrastructure Setup | 1 | 프로젝트 구조, 패키지 설치 |
| 2 | Parameters & Species (01, 05) | 1-2 | 파라미터/초기조건 로딩 |
| 3 | ODE System (02) - 모듈별 | 3-5 | 핵심 ODE 함수 |
| 4 | Rules + Events (03, 04) | 6 | Algebraic rules, events |
| 5 | Single Patient Simulation (06) | 7-8 | 단일 환자 검증 |
| 6 | Virtual Patients (07) | 9-10 | 1000명 VP 생성 |
| 7 | Clinical Trial Simulation | 10 | KEYNOTE-119 보정 |
| 8 | Biomarker Analysis (08) | 11 | 바이오마커 ranking |
| 9 | Visualization (09) | 12 | 논문 Figure 재현 |

### ODE Module Implementation Order
1. Cancer cell dynamics (Gompertzian growth + killing)
2. Naive T cell trafficking
3. APC dynamics & antigen processing
4. T cell activation in LN
5. T cell infiltration & effector functions in tumor
6. T cell exhaustion
7. Macrophage dynamics (recruitment, polarization, phagocytosis)
8. MDSC dynamics
9. Cytokine dynamics (IL-2, IL-10, IL-12, IFN-γ, TGF-β, CCL2)
10. Checkpoint molecules (PD-1, PD-L1, CD28, CTLA-4, CD47)
11. Angiogenesis / carrying capacity
12. Pembrolizumab PK/PD

---

## 10. Technical Considerations

### ODE Solver Selection
- **mrgsolve 권장:** C++ 기반으로 대규모 ODE에 적합, event handling 우수
- **대안:** deSolve의 `lsoda` (stiff system 대응 가능)
- 746개 상태변수 → stiff system일 가능성 높음

### Computational Considerations
- 1000 VP × (pre-treatment + treatment) 시뮬레이션 = 높은 연산 비용
- `parallel` 또는 `future` 패키지로 병렬 처리 필수

### Validation Strategy
1. **단위 검증:** 모든 수식의 단위 일관성 확인
2. **단일 환자 검증:** 논문 Fig. S3, S4의 time series와 비교
3. **VP 검증:** 논문 Fig. 3의 response rates, waterfall plot과 비교
4. **바이오마커 검증:** Data_S2.xlsx 결과와 비교

### Key Challenges
1. **모델 규모:** 746 ODE + 314 rules → 디버깅 난이도 높음
2. **Stiffness:** 세포 역학과 분자 역학의 시간 스케일 차이
3. **Event handling:** Metastatic seeding, drug dosing
4. **Parameter uncertainty:** 87개 파라미터에 강한 실험적 근거 부족
5. **MATLAB → R 변환:** SimBiology reaction format을 R ODE로 변환

---

## 11. Data S1 Excel Structure (참조)

| Sheet | Rows | Description |
|-------|------|-------------|
| Compartments | 131 | 구획 정의 (이름, 용량, 단위) |
| Species | 747 | 상태변수 정의 (초기값, 단위, 위치) |
| Parameters | 738 | 파라미터 (값, 단위, 설명) |
| Uncertain parameters | 87 | 불확실 파라미터 |
| Reactions | 1429 | ODE 반응식 |
| Rules | 315 | 대수 규칙 |
| Events | 66 | 이벤트 (seeding, dosing) |
| Parameter distributions | 128 | VP 생성용 분포 |

---

## 12. Reference Resources

- **논문:** Arulraj et al., Sci. Adv. 9, eadg0289 (2023)
- **MATLAB Code:** http://dx.doi.org/10.17632/r46rk4vwdv.1
- **Base Model (Wang et al.):** iScience 25, 104702 (2022)
- **Data S1:** Model reactions and parameters (adg0289_Data_S1.xlsx)
- **Data S2:** Biomarker analysis results (adg0289_Data_S2.xlsx)
- **Clinical Trial:** KEYNOTE-119 (Winer et al., Lancet Oncol. 2021)

---

*Last updated: 2026-03-07*
