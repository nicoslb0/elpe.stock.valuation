
# ============================================================
#  HELLENiQ ENERGY — FINANCIAL ANALYSIS
#  Single-file script  |  Academic & Professional Grade
#  Last updated: May 2026
# ============================================================


# ---- 0. WORKING DIRECTORY --------------------------------
setwd("/Users/nikoslamprou/Desktop/MScBusinessEcon/Corporate.finance/final.elpe")
message("Working directory: ", getwd())
# ---- 1. LIBRARIES ----------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(lmtest)
  library(sandwich)
  library(car)
  library(tseries)
  library(quantmod)
  library(factoextra)
  library(shiny)
  library(shinydashboard)
  library(DT)
  library(plotly)
  library(depmixS4)
  library(lubridate)
  library(zoo)
  library(rugarch)
  library(vars)
  library(strucchange)
  library(urca)
})
options(stringsAsFactors = FALSE)
dir.create("output", showWarnings = FALSE)

# ============================================================
# SECTION A — MARKET DATA & ECONOMETRICS
# ============================================================

# ---- A.1  CSV DOWNLOAD (only if missing) -----------------
get_csv <- function(symbol, file, from, to) {
  if (!file.exists(file)) {
    message("Downloading ", symbol, " -> ", file)
    obj <- getSymbols(symbol, src = "yahoo",
                      from = from, to = to,
                      auto.assign = FALSE)
    write.csv(as.data.frame(na.omit(obj)), file)
  } else {
    message("Found: ", file)
  }
}

get_csv("ELPE.AT", "ELPE_AT_5yr.csv", "2021-04-07", "2026-04-07")
get_csv("GD.AT",   "GD_AT_5yr.csv",   "2021-04-07", "2026-04-07")
get_csv("BZ=F",    "Brent_5y.csv",    "2019-01-01",
        as.character(Sys.Date()))

# ---- A.2  LOAD & CLEAN PRICES ----------------------------
read_price_csv <- function(file, close_col) {
  df <- readr::read_csv(file, show_col_types = FALSE)
  df <- df[, c(1, which(names(df) == close_col))]
  names(df) <- c("Date", "price")
  df$Date   <- as.Date(df$Date)
  df[order(df$Date), ]
}

elpe_raw <- read_price_csv("ELPE_AT_5yr.csv", "ELPE.AT.Close") |>
  rename(elpe_close = price)

gd_raw   <- read_price_csv("GD_AT_5yr.csv",   "GD.AT.Close")   |>
  rename(mkt_close = price)

returns <- elpe_raw |>
  inner_join(gd_raw, by = "Date") |>
  mutate(
    elpe_ret = log(elpe_close / dplyr::lag(elpe_close)),
    mkt_ret  = log(mkt_close  / dplyr::lag(mkt_close))
  ) |>
  filter(is.finite(elpe_ret), is.finite(mkt_ret)) |>
  as.data.frame()

# Monthly returns for CAPM regression — matches Excel 5Y monthly OLS (beta ≈ 0.74)
# Daily OLS under-estimates systematic risk due to non-synchronous trading & micro-structure noise
returns_monthly <- elpe_raw |>
  inner_join(gd_raw, by = "Date") |>
  mutate(ym = format(Date, "%Y-%m")) |>
  group_by(ym) |>
  slice_tail(n = 1) |>
  ungroup() |>
  arrange(Date) |>
  mutate(
    elpe_ret = log(elpe_close / dplyr::lag(elpe_close)),
    mkt_ret  = log(mkt_close  / dplyr::lag(mkt_close))
  ) |>
  filter(is.finite(elpe_ret), is.finite(mkt_ret)) |>
  as.data.frame()

message(sprintf("ELPE: %d rows  |  GD: %d rows", nrow(elpe_raw), nrow(gd_raw)))
message(sprintf("Daily observations: %d  |  Monthly observations: %d",
                nrow(returns), nrow(returns_monthly)))

# ---- A.3  CAPM OLS + DIAGNOSTICS -------------------------
# Use monthly returns — consistent with Excel WACC sheet (5Y monthly OLS vs ASE GD, β ≈ 0.74)
model_capm <- lm(elpe_ret ~ mkt_ret, data = returns_monthly,
                 model = TRUE, x = TRUE, y = TRUE)

beta_ols  <- unname(coef(model_capm)["mkt_ret"])
alpha_ols <- unname(coef(model_capm)["(Intercept)"])
r2_capm   <- summary(model_capm)$r.squared

# Robust standard errors
robust_hc3 <- lmtest::coeftest(model_capm,
                               vcov = sandwich::vcovHC(model_capm, type = "HC3"))
robust_nw  <- lmtest::coeftest(model_capm,
                               vcov = sandwich::NeweyWest(model_capm, lag = 5))

beta_hc3 <- unname(robust_hc3["mkt_ret", "Estimate"])
beta_nw  <- unname(robust_nw["mkt_ret",  "Estimate"])

# Diagnostic tests
bp  <- lmtest::bptest(model_capm)              # Breusch-Pagan
dw  <- lmtest::dwtest(model_capm)              # Durbin-Watson
lb  <- Box.test(resid(model_capm),
                type = "Ljung-Box", lag = 10)  # Ljung-Box
jb  <- tseries::jarque.bera.test(resid(model_capm))  # Jarque-Bera

diag_tbl <- tibble(
  Test      = c("Breusch-Pagan",
                "Durbin-Watson",
                "Ljung-Box (lag 10)",
                "Jarque-Bera"),
  Statistic = round(c(bp$statistic, dw$statistic,
                      lb$statistic, jb$statistic), 4),
  `P-value` = round(c(bp$p.value,  dw$p.value,
                      lb$p.value,  jb$p.value),  4),
  Conclusion = c(
    ifelse(bp$p.value < 0.05,
           "Heteroskedasticity present — use robust SE",
           "No heteroskedasticity"),
    ifelse(dw$p.value < 0.05,
           "Autocorrelation present — use NW SE",
           "No autocorrelation"),
    ifelse(lb$p.value < 0.05,
           "Serial correlation in residuals",
           "No serial correlation"),
    ifelse(jb$p.value < 0.05,
           "Non-normal residuals (fat tails expected)",
           "Residuals approx. normal")
  )
)

message(sprintf(
  "Beta: %.4f  Alpha: %.6f  R2: %.3f",
  beta_ols, alpha_ols, r2_capm))
message(sprintf(
  "HC3 Beta: %.4f  NW Beta: %.4f",
  beta_hc3, beta_nw))

# ---- A.4  ROLLING BETA (reusable function) ---------------
roll_beta <- function(y, x, window = 250L) {
  n  <- length(y)
  br <- rep(NA_real_, n)
  for (i in window:n) {
    idx   <- (i - window + 1L):i
    br[i] <- coef(lm(y[idx] ~ x[idx]))[2L]
  }
  br
}

returns$beta_1y <- roll_beta(returns$elpe_ret, returns$mkt_ret)
mean_roll       <- mean(returns$beta_1y, na.rm = TRUE)

# Blume (1975) shrinkage toward 1
beta_blume <- 0.67 * beta_ols + 0.33 * 1.0

# Hamada re-levered beta (target D/E = 0.79x, tax = 22%)
# Use Blume-adjusted beta as input to Hamada — matches Excel WACC sheet logic:
#   βU = βL_blume / [1 + (1-t)×D/E]  →  βL* = βU × [1 + (1-t)×target_D/E]
#   With monthly βOLS ≈ 0.74: βBlume ≈ 0.83, βU ≈ 0.445, βL* ≈ 0.72
current_DE <- 1.099  # FY2025: Total Financial Debt / Equity (ELPE Valuation.xlsx, sheet 11)
target_DE  <- 0.79
tax_rate <- 0.22
beta_unlevered <- beta_blume / (1 + (1 - tax_rate) * current_DE)
beta_relevered <- beta_unlevered * (1 + (1 - tax_rate) * target_DE)
rolling_beta_csv <- returns |> filter(!is.na(beta_1y))

message(sprintf(
  "Rolling beta — Static: %.4f  Mean 1Y: %.4f  Blume: %.4f  Re-levered: %.4f",
  beta_ols, mean_roll, beta_blume, beta_relevered))

# ---- A.5  MOH COMPARISON ---------------------------------
moh_xts <- getSymbols("MOH.AT", src = "yahoo",
                      from = "2021-04-07", to = "2026-04-07",
                      auto.assign = FALSE)

moh_df <- tibble(
  Date      = index(moh_xts),
  moh_close = as.numeric(Cl(moh_xts))
) |>
  inner_join(rename(gd_raw, gd_close = mkt_close), by = "Date") |>
  arrange(Date) |>
  mutate(
    moh_ret = log(moh_close / dplyr::lag(moh_close)),
    mkt_ret = log(gd_close  / dplyr::lag(gd_close))
  ) |>
  filter(is.finite(moh_ret), is.finite(mkt_ret)) |>
  as.data.frame()

moh_df$beta_1y  <- roll_beta(moh_df$moh_ret, moh_df$mkt_ret)
beta_moh_static <- unname(coef(lm(moh_ret ~ mkt_ret,
                                  data = moh_df))["mkt_ret"])

df_compare <- bind_rows(
  tibble(Date    = as.Date(returns$Date[!is.na(returns$beta_1y)]),
         Company = "ELPE",
         beta_roll = returns$beta_1y[!is.na(returns$beta_1y)]),
  tibble(Date    = as.Date(moh_df$Date[!is.na(moh_df$beta_1y)]),
         Company = "MOH",
         beta_roll = moh_df$beta_1y[!is.na(moh_df$beta_1y)])
)

message(sprintf("MOH static beta: %.4f", beta_moh_static))

# ---- A.6  MULTI-FACTOR MODEL (Market + Brent) ------------
brent_raw <- readr::read_csv("Brent_5y.csv",
                             show_col_types = FALSE) |>
  rename(Date = 1)

brent_raw$Date <- as.Date(brent_raw$Date)

# Dynamically find adjusted/close column
oil_col <- names(brent_raw)[
  grep("Adjusted|Close", names(brent_raw), ignore.case = TRUE)[1]]

brent <- brent_raw |>
  transmute(Date,
            brent_close = as.numeric(.data[[oil_col]])) |>
  filter(is.finite(brent_close)) |>
  mutate(oil_ret = log(brent_close / dplyr::lag(brent_close))) |>
  dplyr::select(Date, oil_ret) |>
  filter(is.finite(oil_ret))

mf_data <- returns |>
  transmute(Date = as.Date(Date),
            elpe_ret = as.numeric(elpe_ret),
            mkt_ret  = as.numeric(mkt_ret)) |>
  inner_join(brent, by = "Date") |>
  filter(is.finite(elpe_ret),
         is.finite(mkt_ret),
         is.finite(oil_ret)) |>
  as.data.frame()

model_mf <- lm(elpe_ret ~ mkt_ret + oil_ret,
               data = mf_data,
               model = TRUE, x = TRUE, y = TRUE)

beta_oil <- unname(coef(model_mf)["oil_ret"])
r2_mf    <- summary(model_mf)$r.squared

# Nested model F-test: does Brent add explanatory power?
capm_nested <- lm(elpe_ret ~ mkt_ret,
                  data = mf_data,
                  model = TRUE, x = TRUE, y = TRUE)
nested_test <- anova(capm_nested, model_mf)

message(sprintf("MF R2: %.4f  Oil beta: %.4f", r2_mf, beta_oil))
message(sprintf("Nested F-test p-value: %.4f",
                nested_test[["Pr(>F)"]][2]))

# ============================================================
# SECTION B — COST OF CAPITAL & VALUATION
# ============================================================

# ---- B.0  SHARED INPUTS ----------------------------------
rf        <- 0.03081   # ECB Benchmark AAA 10Y, Apr 2026 (ELPE Valuation.xlsx, sheet 11)
erp       <- 0.0752    # Damodaran Apr-2026: US base 4.67% + Greece CRP 2.85%
kd_pretax <- 0.04668   # Blended: 65% synthetic (Rf+CDS 4.63%) + 35% effective (4.74%)
tax_rate  <- 0.22
we        <- 0.5842    # Market-value weight: €3,006M / (€3,006M + €2,140M) — Apr 2026
wd        <- 1 - we
kd        <- kd_pretax * (1 - tax_rate)

# ---- B.1  BOTTOM-UP BETA (Damodaran Pure Play — Base Case) --------
# Revenue-weighted composite unlevered beta, re-levered at actual D/E (Hamada)
# Source: Damodaran Global Jan-2026 sector betas (globals.tern.edu.au); ELPE revenue FY2025
segment_betas <- tibble::tibble(
  Segment  = c("Refining & Supply", "Distribution/EKO",
               "Renewables (RES)", "Power & Gas", "Chemicals"),
  beta_u   = c(0.63, 0.46, 0.57, 1.00, 0.95),
  revenue  = c(5699186, 4916209, 640510, 71000, 284116),
  source   = c("Damodaran Jan-2026 Oil/Gas Integrated",
               "Damodaran Jan-2026 Oil/Gas Distribution",
               "Damodaran Jan-2026 Green & Renewable Energy",
               "Damodaran Jan-2026 Power",
               "Damodaran Jan-2026 Chemical (Basic)")
) |> dplyr::mutate(
  weight  = revenue / sum(revenue),
  contrib = weight * beta_u
)

beta_u_bu <- sum(segment_betas$contrib)                              # 0.5219
beta_bu   <- beta_u_bu * (1 + (1 - tax_rate) * current_DE)          # 0.9692 at D/E=1.099

message(sprintf("Bottom-up βU=%.4f  βL*=%.4f (re-levered at D/E=%.3f)",
                beta_u_bu, beta_bu, current_DE))

# ---- B.2  TWO WACC ESTIMATES --------------------------------

# WACC 1 — OLS-based (robustness check; market-implied, thin ASE trading biases β down)
ke_ols   <- rf + beta_relevered * erp
wacc_ols <- we * ke_ols + wd * kd

# WACC 2 — Bottom-up (base case; Damodaran pure-play segments, preferred for multi-segment conglomerate)
ke_bu   <- rf + beta_bu * erp
wacc_bu <- we * ke_bu + wd * kd

# Base case for all DCF / EVA / Reverse DCF
ke      <- ke_bu
wacc    <- wacc_bu

message(sprintf("OLS WACC  : Ke=%.2f%%  WACC=%.2f%%  (β_rel=%.4f)",
                ke_ols * 100, wacc_ols * 100, beta_relevered))
message(sprintf("Bottom-up WACC (BASE): Ke=%.2f%%  WACC=%.2f%%  (β_bu=%.4f)",
                ke_bu * 100, wacc_bu * 100, beta_bu))

wacc_comparison_tbl <- tibble::tibble(
  Parameter = c("Beta methodology", "Levered beta",
                "Cost of Equity (Ke)", "WACC",
                "WACC–g spread (g=1.8%)", "Academic basis"),
  `OLS (Robustness)`   = c(
    "5Y monthly OLS + Blume adj. + Hamada re-lever",
    sprintf("%.4f", beta_relevered),
    sprintf("%.2f%%", ke_ols * 100),
    sprintf("%.2f%%", wacc_ols * 100),
    sprintf("%.2f pp", (wacc_ols - 0.018) * 100),
    "CAPM; downward-biased on thin ASE"
  ),
  `Bottom-up (Base Case)` = c(
    "Damodaran Jan-2026 pure-play segments, revenue-weighted",
    sprintf("%.4f", beta_bu),
    sprintf("%.2f%%", ke_bu * 100),
    sprintf("%.2f%%", wacc_bu * 100),
    sprintf("%.2f pp", (wacc_bu - 0.018) * 100),
    "Damodaran (2012); preferred for conglomerates"
  )
)

wacc_summary <- tibble::tibble(
  Parameter = c(
    "Risk-Free Rate (Rf)",
    "Equity Risk Premium (ERP)",
    "─── OLS Beta path ───",
    "  Raw OLS Beta (monthly, 5Y)",
    "  Blume-Adjusted Beta",
    "  Re-levered Beta (Hamada, D/E=0.79x)",
    "  Cost of Equity — OLS (Ke_OLS)",
    "  WACC — OLS (robustness check)",
    "─── Bottom-up Beta path ───",
    "  Composite Unlevered Beta (βU, revenue-weighted)",
    "  Re-levered Bottom-up Beta (βL*, D/E=1.099x)",
    "  Cost of Equity — Bottom-up (Ke_BU)",
    "  WACC — Bottom-up ★ BASE CASE ★",
    "─── Debt / Capital Structure ───",
    "Pre-Tax Cost of Debt (Kd)",
    "After-Tax Cost of Debt",
    "Market Cap (EUR M)",
    "Net Debt (EUR M)",
    "Weight of Equity (We)",
    "Weight of Debt (Wd)"
  ),
  Value = c(
    sprintf("%.3f%%", rf * 100),
    "7.52%",
    "",
    sprintf("%.4f", beta_ols),
    sprintf("%.4f", beta_blume),
    sprintf("%.4f", beta_relevered),
    sprintf("%.2f%%", ke_ols * 100),
    sprintf("%.2f%%", wacc_ols * 100),
    "",
    sprintf("%.4f", beta_u_bu),
    sprintf("%.4f", beta_bu),
    sprintf("%.2f%%", ke_bu * 100),
    sprintf("%.2f%%", wacc_bu * 100),
    "",
    sprintf("%.2f%%", kd_pretax * 100),
    sprintf("%.2f%%", kd * 100),
    "3,006", "2,140",
    sprintf("%.1f%%", we * 100),
    sprintf("%.1f%%", wd * 100)
  ),
  Source = c(
    "ECB Benchmark AAA 10Y, Apr 2026",
    "Damodaran Apr-2026 (US base 4.67% + Greece CRP 2.85%)",
    "", "5Y monthly OLS vs ASE GD",
    "Blume (1975): β×2/3 + 1×1/3",
    "Hamada (1972): Blume-βU × [1+(1-t)×target D/E=0.79x]",
    "Rf + β_OLS_relevered × ERP",
    "We×Ke_OLS + Wd×Kd×(1-t)",
    "", "Damodaran Jan-2026 sector βU, revenue-weighted",
    "βU × [1+(1-t)×D/E] at actual D/E=1.099x",
    "Rf + β_BU × ERP",
    "We×Ke_BU + Wd×Kd×(1-t) — used in all DCF/EVA",
    "",
    "65% synthetic (Rf+CDS) + 35% effective rate",
    "Kd × (1 − 22%)",
    "~€9.82 × 306.3M shares, Apr 2026",
    "FY2025 Balance Sheet",
    "Market-value weight", "Market-value weight"
  )
)

# Scenario analysis — use bottom-up beta as base
scenario_wacc <- tibble::tibble(
  Scenario = c("Bear (high β)", "Base — Bottom-up", "OLS Robustness"),
  Beta     = c(1.200,  beta_bu,   beta_relevered),
  Ke       = round(rf + c(1.200, beta_bu, beta_relevered) * erp, 4),
  `Kd(at)` = rep(kd, 3),
  WACC     = round(we * (rf + c(1.200, beta_bu, beta_relevered) * erp) + wd * kd, 4)
) |>
  dplyr::mutate(dplyr::across(c(Ke, `Kd(at)`, WACC),
                              ~ paste0(round(. * 100, 2), "%")))

# Sensitivity table: WACC vs g
sens_tbl <- expand_grid(
  WACC = c(0.045, 0.050, 0.057, 0.060, 0.065, 0.070),
  g    = c(0.005, 0.010, 0.015, 0.018, 0.020, 0.025)
) |>
  mutate(
    spread    = round((WACC - g) * 100, 2),
    WACC_label = paste0(round(WACC * 100, 2), "%")
  ) |>
  dplyr::select(WACC_label, g, spread) |>
  pivot_wider(names_from = g, values_from = spread,
              names_prefix = "g=")


# ============================================================
# C. FUNDAMENTALS / COMPARATIVES / PCA / SENTIMENT
# ============================================================

# ---- C.1 ALTMAN Z-SCORE verified from Excel sheet 7 --------
years <- 2021:2025

altmanclean <- tibble::tibble(
  Year = years,
  ZELPE = c(1.892, 2.899, 2.670, 2.733, 2.352),
  ZMOH  = c(3.625, 4.257, 3.623, 3.623, 3.517)
)

altmanlong <- altmanclean %>%
  tidyr::pivot_longer(
    cols = c(ZELPE, ZMOH),
    names_to = "Company",
    values_to = "ZScore"
  ) %>%
  dplyr::mutate(
    Company = dplyr::case_match(
      Company,
      "ZELPE" ~ "ELPE",
      "ZMOH"  ~ "MOH",
      .default = Company
    )
  )

# ---- C.2 DUPONT verified from Excel sheet 6 ---------------
elpedupont <- tibble::tibble(
  Company = "ELPE",
  Year = years,
  Margin = c(0.03699, 0.06169, 0.03754, 0.00488, 0.01524),
  Turnover = c(1.2627, 1.7699, 1.5360, 1.6098, 1.4233),
  Leverage = c(3.6720, 3.3758, 2.9381, 2.7787, 2.9728)
)

mohdupont <- tibble::tibble(
  Company = "MOH",
  Year = years,
  Margin = c(0.01971, 0.05816, 0.06058, 0.02357, 0.05668),
  Turnover = c(2.1240, 2.3100, 1.7580, 1.6740, 1.4280),
  Leverage = c(4.0580, 3.3660, 2.7320, 2.6380, 2.3950)
)

dupontall <- dplyr::bind_rows(elpedupont, mohdupont) %>%
  dplyr::mutate(ROE = Margin * Turnover * Leverage)

dupont2025 <- dupontall %>%
  dplyr::filter(Year == 2025) %>%
  tidyr::pivot_longer(
    cols = c(Margin, Turnover, Leverage),
    names_to = "Component",
    values_to = "Value"
  )

# ---- C.3 FINANCIAL RATIOS verified from Excel -------------
financialratios <- tibble::tibble(
  Metric = c(
    "Current Ratio",
    "Quick Ratio",
    "Cash Ratio",
    "Debt/Equity (loans only)",
    "Net Debt / Adj. EBITDA",
    "Interest Coverage (EBITDA / Interest)",
    "Adj. EBITDA Margin (%)",
    "Pre-Tax ROA (%)",
    "Pre-Tax ROE (%)"
  ),
  ELPE2024 = c(1.333, 0.741, 0.279, 0.873, 2.072, 6.089, 6.8, 4.1, 11.4),
  MOH2024  = c(2.000, 1.363, 0.751, 0.948, 1.357, 4.710, 9.0, 8.6, 23.2),
  ELPE2025 = c(1.407, 0.861, 0.359, 1.099, 1.890, 8.185, 9.7, 3.1, 9.3),
  MOH2025  = c(1.656, 1.203, 0.685, 0.805, 1.269, 3.967, 9.1, 10.8, 27.1)
)

# ---- C.4 PCA verified values ------------------------------
pcadata <- tibble::tibble(
  Company = rep(c("ELPE", "MOH"), each = 5),
  Year = rep(2021:2025, 2),
  Z = c(1.89, 2.90, 2.67, 2.73, 2.35, 3.63, 4.26, 3.62, 3.62, 3.52),
  ROE = c(0.17, 0.39, 0.17, 0.02, 0.07, 0.16, 0.45, 0.29, 0.10, 0.39),
  EBITDA = c(0.037, 0.062, 0.038, 0.005, 0.015, 0.020, 0.058, 0.061, 0.024, 0.057),
  DE = c(2.67, 2.38, 1.94, 1.78, 1.97, 3.06, 2.37, 1.73, 1.64, 1.40),
  CR = c(1.15, 1.28, 1.25, 1.33, 1.41, 1.89, 2.12, 1.95, 2.00, 1.66),
  NDE = c(2.45, 2.10, 1.85, 2.07, 1.89, 1.62, 1.48, 1.55, 1.36, 1.27)
)

pcadf <- pcadata %>% dplyr::select(-Year)

pcaresult <- stats::prcomp(
  pcadf %>% dplyr::select(-Company),
  scale. = TRUE
)

# ---- C.5 SENTIMENT ----------------------------------------
sentimentdf <- tibble::tibble(
  Category = c("Positive", "Negative"),
  NetScore = c(12, -8)
)

# ============================================================
# D. CRACK SPREAD / HMM / REGRESSION
# ============================================================

# ---- D.1 HMM ON CRACK SPREAD -------------------------------
quantmod::getSymbols("CL=F", from = "2015-01-01", auto.assign = TRUE, warnings = FALSE)
quantmod::getSymbols("HO=F", from = "2015-01-01", auto.assign = TRUE, warnings = FALSE)

crack <- na.omit(merge(Cl(`HO=F`), Cl(`CL=F`)))
crack <- data.frame(
  date = zoo::index(crack),
  value = as.numeric(crack[, 1]) - as.numeric(crack[, 2])
)

mod <- depmixS4::depmix(value ~ 1, data = crack, nstates = 3, family = gaussian())
fit <- depmixS4::fit(mod, verbose = FALSE)

trans <- matrix(
  c(
    0.988, 0.006, 0.006,
    0.007, 0.992, 0.001,
    0.007, 0.001, 0.992
  ),
  nrow = 3, byrow = TRUE
)

statemeans <- c(24.245, 46.457, 15.094)  # S1 Normal, S2 Bull, S3 Bear
statesds   <- c(2.245, 14.372, 3.364)

# ---- D.2 CRACK SPREAD REGRESSION ---------------------------

possibleexcel <- c(
  "ELPE Valuation.xlsx",
  "./ELPE Valuation.xlsx",
  "ELPE-Historical-Analysis.xlsx",
  "./ELPE-Historical-Analysis.xlsx"
)

excelfile <- possibleexcel[file.exists(possibleexcel)][1]
if (is.na(excelfile)) {
  excelfile <- NA_character_
}

# Annual crack spread 2021-2025
crackyearly <- crack %>%
  dplyr::mutate(year = lubridate::year(date)) %>%
  dplyr::group_by(year) %>%
  dplyr::summarise(avgspread = mean(value, na.rm = TRUE), .groups = "drop") %>%
  dplyr::filter(year %in% 2021:2025)

# ------------------------------------------------------------
# EBITDA margin history
# Try Excel first; if it fails, use verified fallback values
# ------------------------------------------------------------

ebitdamarginhist <- NULL

if (!is.na(excelfile) && file.exists(excelfile)) {
  marginsraw <- readxl::read_excel(excelfile, sheet = "4_MARGINS", col_names = FALSE)
  
  adjrows <- which(
    apply(marginsraw, 1, function(r) {
      any(grepl("Adj\\.?\\s*EBITDA", as.character(r), ignore.case = TRUE))
    })
  )
  
  if (length(adjrows) > 0) {
    adjebitdarow <- adjrows[1]
    rowvals <- as.character(unlist(marginsraw[adjebitdarow, ]))
    rowvals <- gsub("%", "", rowvals)
    rowvals <- gsub(",", ".", rowvals)
    numvals <- suppressWarnings(as.numeric(rowvals))
    numvals <- numvals[!is.na(numvals)]
    
    if (length(numvals) >= 5) {
      ebitdamarginhist <- numvals[1:5] / 100
    }
  }
}

# Verified fallback if Excel parsing fails
if (is.null(ebitdamarginhist) || length(ebitdamarginhist) != 5 || any(is.na(ebitdamarginhist))) {
  ebitdamarginhist <- c(0.068, 0.097, 0.037, 0.068, 0.097)
  message("Using fallback EBITDA margin series for 2021-2025.")
}

# Safety check for crack spread too
if (nrow(crackyearly) != 5 || any(is.na(crackyearly$avgspread))) {
  crackyearly <- tibble::tibble(
    year = 2021:2025,
    avgspread = c(18.0, 22.5, 18.2, 12.8, 14.5)
  )
  message("Using fallback crack spread series for 2021-2025.")
}

# EBITDA margin regression
regdata <- data.frame(
  year = 2021:2025,
  margin = ebitdamarginhist,
  spread = crackyearly$avgspread
)

print(regdata)

modelreg <- lm(margin ~ spread, data = regdata)
alphareg <- coef(modelreg)[1]
betareg  <- coef(modelreg)[2]

# Revenue growth regression
revenuehist <- c(9222235, 14508068, 12803061, 12767894, 11614643)
revgrowthhist <- diff(revenuehist) / head(revenuehist, -1)

spreadforgrowth <- crackyearly %>%
  dplyr::filter(year %in% 2022:2025) %>%
  dplyr::pull(avgspread)

if (length(spreadforgrowth) != 4 || any(is.na(spreadforgrowth))) {
  spreadforgrowth <- c(22.5, 18.2, 12.8, 14.5)
  message("Using fallback spread series for revenue growth regression.")
}

reggrowthdata <- data.frame(
  year = 2022:2025,
  revgrowth = revgrowthhist,
  spread = spreadforgrowth
)

modelgrowth <- lm(revgrowth ~ spread, data = reggrowthdata)
alphagrowth <- coef(modelgrowth)[1]
betagrowth  <- coef(modelgrowth)[2]
r2growth    <- summary(modelgrowth)$r.squared

predpts <- data.frame(
  spread = statemeans,
  label = c("Normal", "Bull", "Bear"),
  growth = alphagrowth + betagrowth * statemeans
)

predmarginpts <- data.frame(
  spread = statemeans,
  label = c("Normal", "Bull", "Bear"),
  margin = alphareg + betareg * statemeans
)

crackdf <- dplyr::bind_rows(
  reggrowthdata %>% dplyr::mutate(Type = "Revenue Growth"),
  regdata %>% dplyr::mutate(Type = "EBITDA Margin")
)

message(sprintf("Margin regression alpha = %.4f | beta = %.6f | R2 = %.3f",
                alphareg, betareg, summary(modelreg)$r.squared))
message(sprintf("Growth regression alpha = %.4f | beta = %.6f | R2 = %.3f",
                alphagrowth, betagrowth, r2growth))

# ============================================================
# E. MARKOV DCF - BASE CASE
# ============================================================

nsim <- 1000L
nquarters <- 20L
yearssim <- 5L
nruns <- 200L

basefcff    <- c(284000, 312000, 335000, 358000, 374000)
baserevenue <- c(12331000, 12214000, 11946000, 11555000, 11141000)
baseebitda  <- c(1309000, 1412000, 1400000, 1378000, 1359000)
basenopat   <- c(684000, 752400, 745200, 727200, 717120)

regimemultmean <- c("1" = 1.00, "2" = 1.20, "3" = 0.75)
regimemultsd   <- c("1" = 0.05, "2" = 0.08, "3" = 0.10)

waccval        <- wacc_bu   # Base case: bottom-up Damodaran WACC (~7.57%)
growthterminal <- 0.018
exitmultiple   <- 7.5
wgordon        <- 0.40
wevebitda      <- 0.60
netdebtlatest  <- 2140000
shares         <- 306300000
currentprice   <- 9.50
payoutratio    <- 0.60

discfactors <- sapply(seq_len(yearssim), function(t) 1 / (1 + waccval)^t)

callvalue <- function(S, K, T, r, sigma) {
  d1 <- (log(S / K) + (r + sigma^2 / 2) * T) / (sigma * sqrt(T))
  d2 <- d1 - sigma * sqrt(T)
  S * pnorm(d1) - K * exp(-r * T) * pnorm(d2)
}

Sres <- 800000
Kres <- 750000
Tres <- 3
sigmares <- 0.35
totaloption <- callvalue(Sres, Kres, Tres, waccval, sigmares)

probsmat <- posterior(fit, type = "viterbi")

runmarkovdcf <- function(seedvalue,
                         nsim,
                         nquarters,
                         yearssim,
                         probsmat,
                         trans,
                         regimemultmean,
                         regimemultsd,
                         basefcff,
                         baserevenue,
                         basenopat,
                         baseebitda,
                         discfactors,
                         growthterminal,
                         waccval,
                         exitmultiple,
                         wgordon,
                         wevebitda,
                         totaloption,
                         netdebtlatest,
                         shares,
                         currentprice) {
  
  set.seed(seedvalue)
  
  statessim <- matrix(NA_integer_, nrow = nsim, ncol = nquarters)
  initprobs <- as.numeric(probsmat[nrow(probsmat), 2:4])
  initprobs <- initprobs / sum(initprobs)
  
  for (i in seq_len(nsim)) {
    statessim[i, 1] <- sample(1:3, 1, prob = initprobs)
    for (t in 2:nquarters) {
      statessim[i, t] <- sample(1:3, 1, prob = trans[statessim[i, t - 1], ])
    }
  }
  
  fcfsim     <- matrix(NA_real_, nsim, yearssim)
  revenuesim <- matrix(NA_real_, nsim, yearssim)
  nopatsim   <- matrix(NA_real_, nsim, yearssim)
  ebitdasim  <- matrix(NA_real_, nsim, yearssim)
  
  for (i in seq_len(nsim)) {
    quarterlystates <- statessim[i, ]
    for (y in seq_len(yearssim)) {
      qs <- quarterlystates[((y - 1) * 4 + 1):(y * 4)]
      qm <- sapply(as.character(qs), function(s) {
        rnorm(1, mean = regimemultmean[s], sd = regimemultsd[s])
      })
      annualmult <- max(0.40, mean(qm))
      fcfsim[i, y]     <- basefcff[y]    * annualmult
      revenuesim[i, y] <- baserevenue[y] * annualmult
      nopatsim[i, y]   <- basenopat[y]   * annualmult
      ebitdasim[i, y]  <- baseebitda[y]  * annualmult
    }
  }
  
  pvfcf <- fcfsim %*% discfactors
  tvgordon <- fcfsim[, yearssim] * (1 + growthterminal) /
    (waccval - growthterminal) / ((1 + waccval)^yearssim)
  tvevebitda <- ebitdasim[, yearssim] * exitmultiple /
    ((1 + waccval)^yearssim)
  
  terminalvalue <- wgordon * tvgordon + wevebitda * tvevebitda
  evsim <- as.numeric(pvfcf) + terminalvalue + totaloption
  pricepershare <- (evsim - netdebtlatest) * 1000 / shares
  
  data.frame(
    seed = seedvalue,
    meanprice = mean(pricepershare, na.rm = TRUE),
    medianprice = median(pricepershare, na.rm = TRUE),
    p05 = as.numeric(quantile(pricepershare, 0.05, na.rm = TRUE)),
    p25 = as.numeric(quantile(pricepershare, 0.25, na.rm = TRUE)),
    p75 = as.numeric(quantile(pricepershare, 0.75, na.rm = TRUE)),
    p95 = as.numeric(quantile(pricepershare, 0.95, na.rm = TRUE)),
    probupside = mean(pricepershare > currentprice, na.rm = TRUE)
  )
}

seedgrid <- 1001:(1000 + nruns)

runresults <- do.call(
  rbind,
  lapply(seedgrid, function(s) {
    runmarkovdcf(
      seedvalue = s,
      nsim = nsim,
      nquarters = nquarters,
      yearssim = yearssim,
      probsmat = probsmat,
      trans = trans,
      regimemultmean = regimemultmean,
      regimemultsd = regimemultsd,
      basefcff = basefcff,
      baserevenue = baserevenue,
      basenopat = basenopat,
      baseebitda = baseebitda,
      discfactors = discfactors,
      growthterminal = growthterminal,
      waccval = waccval,
      exitmultiple = exitmultiple,
      wgordon = wgordon,
      wevebitda = wevebitda,
      totaloption = totaloption,
      netdebtlatest = netdebtlatest,
      shares = shares,
      currentprice = currentprice
    )
  })
)

meanofmeans   <- mean(runresults$meanprice, na.rm = TRUE)
cimeanlow     <- quantile(runresults$meanprice, 0.025, na.rm = TRUE)
cimeanhigh    <- quantile(runresults$meanprice, 0.975, na.rm = TRUE)
meanofmedians <- mean(runresults$medianprice, na.rm = TRUE)
cimedianlow   <- quantile(runresults$medianprice, 0.025, na.rm = TRUE)
cimedianhigh  <- quantile(runresults$medianprice, 0.975, na.rm = TRUE)
meanprobupside <- mean(runresults$probupside, na.rm = TRUE)
ciupsidelow    <- quantile(runresults$probupside, 0.025, na.rm = TRUE)
ciupsidehigh   <- quantile(runresults$probupside, 0.975, na.rm = TRUE)
avgp05 <- mean(runresults$p05, na.rm = TRUE)
avgp25 <- mean(runresults$p25, na.rm = TRUE)
avgp75 <- mean(runresults$p75, na.rm = TRUE)
avgp95 <- mean(runresults$p95, na.rm = TRUE)

rec <- if (meanprobupside >= 0.60) "BUY" else if (meanprobupside >= 0.40) "HOLD" else "SELL"

meanprice   <- meanofmeans
medianprice <- meanofmedians
probupside  <- meanprobupside

comparisontbl <- tibble::tibble(
  Model = "Markov DCF base transitions (multi-run)",
  `Mean EUR` = round(meanofmeans, 2),
  `95% CI Mean` = sprintf("%.2f to %.2f", cimeanlow, cimeanhigh),
  `Median EUR` = round(meanofmedians, 2),
  `95% CI Median` = sprintf("%.2f to %.2f", cimedianlow, cimedianhigh),
  `Avg Upside %` = round((meanofmeans - currentprice) / currentprice * 100, 1),
  `Prob(Upside) %` = round(meanprobupside * 100, 1),
  Recommendation = rec
)

markovsummarytbl <- tibble::tibble(
  Metric = c(
    "Outer runs", "Paths per run", "Average of run means (EUR/share)",
    "95% CI of run means - low", "95% CI of run means - high",
    "Average of run medians (EUR/share)",
    "95% CI of run medians - low", "95% CI of run medians - high",
    "Average P05 across runs", "Average P25 across runs",
    "Average P75 across runs", "Average P95 across runs",
    "Average probability of upside",
    "95% CI probability of upside - low",
    "95% CI probability of upside - high",
    "Recommendation"
  ),
  Value = c(
    nruns, nsim, round(meanofmeans, 4),
    round(cimeanlow, 4), round(cimeanhigh, 4),
    round(meanofmedians, 4),
    round(cimedianlow, 4), round(cimedianhigh, 4),
    round(avgp05, 4), round(avgp25, 4), round(avgp75, 4), round(avgp95, 4),
    round(meanprobupside, 4),
    round(ciupsidelow, 4), round(ciupsidehigh, 4),
    rec
  )
)

# ============================================================
# E.2 MARKOV DCF - MOBILE / HIGHER VOLATILITY
# ============================================================

nrunsmobile <- 200L
transmobile <- matrix(
  c(
    0.85, 0.10, 0.05,
    0.10, 0.85, 0.05,
    0.05, 0.05, 0.90
  ),
  nrow = 3, byrow = TRUE
)

runmobilemarkovdcf <- function(seedvalue,
                               nsim,
                               nquarters,
                               yearssim,
                               transmobile,
                               regimemultmean,
                               regimemultsd,
                               basefcff,
                               baseebitda,
                               discfactors,
                               growthterminal,
                               waccval,
                               exitmultiple,
                               wgordon,
                               wevebitda,
                               totaloption,
                               netdebtlatest,
                               shares,
                               currentprice) {
  
  set.seed(seedvalue)
  
  statesmobile <- matrix(NA_integer_, nrow = nsim, ncol = nquarters)
  
  for (i in seq_len(nsim)) {
    statesmobile[i, 1] <- sample(1:3, 1, prob = c(0.2, 0.6, 0.2))
    for (t in 2:nquarters) {
      statesmobile[i, t] <- sample(1:3, 1, prob = transmobile[statesmobile[i, t - 1], ])
    }
  }
  
  fcfmobile    <- matrix(NA_real_, nsim, yearssim)
  ebitdamobile <- matrix(NA_real_, nsim, yearssim)
  
  for (i in seq_len(nsim)) {
    qsall <- statesmobile[i, ]
    for (y in seq_len(yearssim)) {
      qs <- qsall[((y - 1) * 4 + 1):(y * 4)]
      qm <- sapply(as.character(qs), function(s) {
        rnorm(1, mean = regimemultmean[s], sd = regimemultsd[s])
      })
      annualmult <- max(0.40, mean(qm))
      fcfmobile[i, y]    <- basefcff[y]   * annualmult
      ebitdamobile[i, y] <- baseebitda[y] * annualmult
    }
  }
  
  pvmobile <- fcfmobile %*% discfactors
  tvmobg <- fcfmobile[, yearssim] * (1 + growthterminal) /
    (waccval - growthterminal) / ((1 + waccval)^yearssim)
  tvmobev <- ebitdamobile[, yearssim] * exitmultiple /
    ((1 + waccval)^yearssim)
  
  termmobile <- wgordon * tvmobg + wevebitda * tvmobev
  exmobile <- rowSums(fcfmobile) > sum(basefcff)
  evmobile <- as.numeric(pvmobile) + termmobile + ifelse(exmobile, totaloption, 0)
  pricemobile <- (evmobile - netdebtlatest) * 1000 / shares
  
  data.frame(
    seed = seedvalue,
    meanpricemobile = mean(pricemobile, na.rm = TRUE),
    medianpricemobile = median(pricemobile, na.rm = TRUE),
    p05mobile = as.numeric(quantile(pricemobile, 0.05, na.rm = TRUE)),
    p25mobile = as.numeric(quantile(pricemobile, 0.25, na.rm = TRUE)),
    p75mobile = as.numeric(quantile(pricemobile, 0.75, na.rm = TRUE)),
    p95mobile = as.numeric(quantile(pricemobile, 0.95, na.rm = TRUE)),
    probmobile = mean(pricemobile > currentprice, na.rm = TRUE)
  )
}

seedgridmobile <- 5001:(5000 + nrunsmobile)

mobilerunresults <- do.call(
  rbind,
  lapply(seedgridmobile, function(s) {
    runmobilemarkovdcf(
      seedvalue = s,
      nsim = nsim,
      nquarters = nquarters,
      yearssim = yearssim,
      transmobile = transmobile,
      regimemultmean = regimemultmean,
      regimemultsd = regimemultsd,
      basefcff = basefcff,
      baseebitda = baseebitda,
      discfactors = discfactors,
      growthterminal = growthterminal,
      waccval = waccval,
      exitmultiple = exitmultiple,
      wgordon = wgordon,
      wevebitda = wevebitda,
      totaloption = totaloption,
      netdebtlatest = netdebtlatest,
      shares = shares,
      currentprice = currentprice
    )
  })
)


meanpricemobileavg <- mean(mobilerunresults$meanpricemobile, na.rm = TRUE)
cimobilelow <- quantile(mobilerunresults$meanpricemobile, 0.025, na.rm = TRUE)
cimobilehigh <- quantile(mobilerunresults$meanpricemobile, 0.975, na.rm = TRUE)
medianmobileavg <- mean(mobilerunresults$medianpricemobile, na.rm = TRUE)
cimedmobilelow <- quantile(mobilerunresults$medianpricemobile, 0.025, na.rm = TRUE)
cimedmobilehigh <- quantile(mobilerunresults$medianpricemobile, 0.975, na.rm = TRUE)
probmobileavg <- mean(mobilerunresults$probmobile, na.rm = TRUE)
ciprobmobilelow <- quantile(mobilerunresults$probmobile, 0.025, na.rm = TRUE)
ciprobmobilehigh <- quantile(mobilerunresults$probmobile, 0.975, na.rm = TRUE)
avgp05mobile <- mean(mobilerunresults$p05mobile, na.rm = TRUE)
avgp25mobile <- mean(mobilerunresults$p25mobile, na.rm = TRUE)
avgp75mobile <- mean(mobilerunresults$p75mobile, na.rm = TRUE)
avgp95mobile <- mean(mobilerunresults$p95mobile, na.rm = TRUE)

recmobile <- if (probmobileavg >= 0.60) "BUY" else if (probmobileavg >= 0.40) "HOLD" else "SELL"

meanpricemobile <- meanpricemobileavg
probmobile <- probmobileavg

mobilemarkovsummarytbl <- tibble::tibble(
  Metric = c(
    "Outer runs", "Paths per run", "Average of run means (EUR/share)",
    "95% CI of run means - low", "95% CI of run means - high",
    "Average of run medians (EUR/share)",
    "95% CI of run medians - low", "95% CI of run medians - high",
    "Average P05 across runs", "Average P25 across runs",
    "Average P75 across runs", "Average P95 across runs",
    "Average probability of upside",
    "95% CI probability of upside - low",
    "95% CI probability of upside - high",
    "Recommendation"
  ),
  Value = c(
    nrunsmobile, nsim, round(meanpricemobileavg, 4),
    round(cimobilelow, 4), round(cimobilehigh, 4),
    round(medianmobileavg, 4),
    round(cimedmobilelow, 4), round(cimedmobilehigh, 4),
    round(avgp05mobile, 4), round(avgp25mobile, 4), round(avgp75mobile, 4), round(avgp95mobile, 4),
    round(probmobileavg, 4), round(ciprobmobilelow, 4), round(ciprobmobilehigh, 4),
    recmobile
  )
)

# ============================================================
# F. AFN / COVENANT STRESS TEST + GBM
# ============================================================

set.seed(123)

statessimcompat <- matrix(NA_integer_, nrow = nsim, ncol = nquarters)
initprobs <- as.numeric(probsmat[nrow(probsmat), 2:4])
initprobs <- initprobs / sum(initprobs)

for (i in seq_len(nsim)) {
  statessimcompat[i, 1] <- sample(1:3, 1, prob = initprobs)
  for (t in 2:nquarters) {
    statessimcompat[i, t] <- sample(1:3, 1, prob = trans[statessimcompat[i, t - 1], ])
  }
}

fcfsim <- matrix(NA_real_, nsim, yearssim)
revenuesim <- matrix(NA_real_, nsim, yearssim)
nopatsim <- matrix(NA_real_, nsim, yearssim)
ebitdasim <- matrix(NA_real_, nsim, yearssim)

for (i in seq_len(nsim)) {
  quarterlystates <- statessimcompat[i, ]
  for (y in seq_len(yearssim)) {
    qs <- quarterlystates[((y - 1) * 4 + 1):(y * 4)]
    qm <- sapply(as.character(qs), function(s) {
      rnorm(1, mean = regimemultmean[s], sd = regimemultsd[s])
    })
    annualmult <- max(0.40, mean(qm))
    fcfsim[i, y]     <- basefcff[y]    * annualmult
    revenuesim[i, y] <- baserevenue[y] * annualmult
    nopatsim[i, y]   <- basenopat[y]   * annualmult
    ebitdasim[i, y]  <- baseebitda[y]  * annualmult
  }
}

revenuehist <- c(9700000, 10800000, 12331000, 12214000, 11946000)

pvfcfcompat <- fcfsim %*% discfactors
tvgordoncompat <- fcfsim[, yearssim] * (1 + growthterminal) /
  (waccval - growthterminal) / ((1 + waccval)^yearssim)
tvevebitdacompat <- ebitdasim[, yearssim] * exitmultiple /
  ((1 + waccval)^yearssim)
terminalvaluecompat <- wgordon * tvgordoncompat + wevebitda * tvevebitdacompat
evsimcompat <- as.numeric(pvfcfcompat) + terminalvaluecompat + totaloption
pricepershare <- (evsimcompat - netdebtlatest) * 1000 / shares

initialequity <- 2727908
initialdebt   <- 3272837

debtseries    <- matrix(NA_real_, nsim, yearssim)
equityseries  <- matrix(NA_real_, nsim, yearssim)
afnsim        <- matrix(NA_real_, nsim, yearssim)
dscrsim       <- matrix(NA_real_, nsim, yearssim)
ndebitdasim   <- matrix(NA_real_, nsim, yearssim)
interestsim   <- matrix(NA_real_, nsim, yearssim)

debtseries[, 1]   <- initialdebt
equityseries[, 1] <- initialequity

for (i in seq_len(nsim)) {
  for (t in seq_len(yearssim)) {
    currentdebt   <- debtseries[i, t]
    currentequity <- equityseries[i, t]
    interestt     <- currentdebt * kd_pretax
    netincomet    <- nopatsim[i, t] - interestt
    retainedt     <- max(netincomet, 0) * (1 - payoutratio)
    reqassets     <- revenuesim[i, t] / revenuehist[5] * (initialdebt + initialequity)
    afnt          <- reqassets - currentequity - retainedt - currentdebt
    
    interestsim[i, t]  <- interestt
    dscrsim[i, t]      <- fcfsim[i, t] / (interestt + currentdebt * 0.05)
    ndebitdasim[i, t]  <- currentdebt / ebitdasim[i, t]
    afnsim[i, t]       <- afnt
    
    if (t < yearssim) {
      debtseries[i, t + 1]   <- currentdebt + max(afnt, 0) * 0.50
      equityseries[i, t + 1] <- currentequity + retainedt + max(afnt, 0) * 0.50
    }
  }
}

dematrix <- debtseries / equityseries
criterionA <- apply(dscrsim, 1, function(r) all(r >= 1.3, na.rm = TRUE))
criterionB <- apply(dematrix, 1, function(r) all(r <= 0.79, na.rm = TRUE))
criterionC <- apply(ndebitdasim, 1, function(r) all(r <= 2.5, na.rm = TRUE))

# ---- GBM ---------------------------------------------------
mktprice <- currentprice

S0gbm <- 9.50
mugbm <- 0.118
sigmagbm <- 0.288
stepsgbm <- 252
dtgbm <- 1 / stepsgbm
nsimgbm <- 1000

set.seed(42)
gbmpaths <- matrix(NA_real_, nrow = nsimgbm, ncol = stepsgbm + 1)
gbmpaths[, 1] <- S0gbm

for (t in 2:(stepsgbm + 1)) {
  z <- rnorm(nsimgbm)
  gbmpaths[, t] <- gbmpaths[, t - 1] * exp(
    (mugbm - 0.5 * sigmagbm^2) * dtgbm +
      sigmagbm * sqrt(dtgbm) * z
  )
}

gbmfinal <- gbmpaths[, stepsgbm + 1]

# ============================================================
# H. ADVANCED ECONOMETRICS & RISK METRICS
# ============================================================

# ---- H.1 STATIONARITY (ADF + KPSS) -------------------------
adf_elpe   <- tseries::adf.test(returns$elpe_ret)
adf_mkt    <- tseries::adf.test(returns$mkt_ret)
kpss_elpe  <- tseries::kpss.test(returns$elpe_ret, null = "Level")
kpss_mkt   <- tseries::kpss.test(returns$mkt_ret,  null = "Level")
crack_vals <- na.omit(crack$value)
adf_crack  <- tseries::adf.test(crack_vals)
kpss_crack <- tseries::kpss.test(crack_vals, null = "Level")

stationarity_tbl <- tibble::tibble(
  Series         = c("ELPE log-returns", "GD log-returns", "Crack spread (level)"),
  `ADF stat`     = round(c(adf_elpe$statistic, adf_mkt$statistic, adf_crack$statistic), 4),
  `ADF p-val`    = round(c(adf_elpe$p.value,   adf_mkt$p.value,  adf_crack$p.value),   4),
  `KPSS stat`    = round(c(kpss_elpe$statistic, kpss_mkt$statistic, kpss_crack$statistic), 4),
  `KPSS p-val`   = round(c(kpss_elpe$p.value,   kpss_mkt$p.value,  kpss_crack$p.value),  4),
  `ADF verdict`  = dplyr::case_when(
    c(adf_elpe$p.value, adf_mkt$p.value, adf_crack$p.value) < 0.01 ~ "Stationary ***",
    c(adf_elpe$p.value, adf_mkt$p.value, adf_crack$p.value) < 0.05 ~ "Stationary **",
    c(adf_elpe$p.value, adf_mkt$p.value, adf_crack$p.value) < 0.10 ~ "Stationary *",
    TRUE ~ "Non-stationary"
  ),
  `KPSS verdict` = ifelse(
    c(kpss_elpe$p.value, kpss_mkt$p.value, kpss_crack$p.value) > 0.05,
    "Stationary", "Non-stationary **"
  )
)
message("H.1 Stationarity tests complete.")

# ---- H.2 GARCH(1,1) — Student-t on ELPE returns -----------
garch_spec <- rugarch::ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
garch_fit <- rugarch::ugarchfit(spec = garch_spec,
                                data = returns$elpe_ret,
                                solver = "hybrid")

garch_coefs       <- coef(garch_fit)
garch_alpha1      <- garch_coefs["alpha1"]
garch_beta1       <- garch_coefs["beta1"]
garch_persistence <- garch_alpha1 + garch_beta1
garch_sigma_vec   <- as.numeric(rugarch::sigma(garch_fit))
garch_sigma_last  <- tail(garch_sigma_vec, 1) * sqrt(252)

garch_param_tbl <- tibble::tibble(
  Parameter      = c("ω (omega)", "α₁ (alpha1)", "β₁ (beta1)",
                     "Persistence α+β", "Half-life (days)",
                     "GARCH σ annualised (last obs)", "Static σ (GBM input)"),
  Estimate       = round(c(garch_coefs["omega"], garch_alpha1, garch_beta1,
                            garch_persistence,
                            log(0.5) / log(garch_persistence),
                            garch_sigma_last, sigmagbm), 6),
  Interpretation = c("Long-run variance intercept",
                     "Sensitivity to recent shocks (ARCH)",
                     "Persistence of past variance (GARCH)",
                     "Close to 1 → highly persistent volatility",
                     "Days for vol shock to halve in magnitude",
                     "Latest GARCH-implied annualised volatility",
                     "Constant vol currently used in GBM")
)

garch_vol_df <- tibble::tibble(
  Date     = as.Date(returns$Date),
  cond_vol = garch_sigma_vec * sqrt(252),
  static   = sigmagbm
)
message(sprintf("H.2 GARCH: α=%.4f β=%.4f persistence=%.4f σ_last=%.4f",
                garch_alpha1, garch_beta1, garch_persistence, garch_sigma_last))

# ---- H.3 VAR(p) + GRANGER CAUSALITY + IRF ------------------
var_df    <- mf_data[, c("elpe_ret", "mkt_ret", "oil_ret")]
lag_sel   <- vars::VARselect(var_df, lag.max = 10, type = "const")
opt_lag   <- max(1L, as.integer(lag_sel$selection["AIC(n)"]))
var_model <- vars::VAR(var_df, p = opt_lag, type = "const")

gc_oil_elpe <- vars::causality(var_model, cause = "oil_ret")
gc_mkt_elpe <- vars::causality(var_model, cause = "mkt_ret")
gc_elpe_oil <- vars::causality(var_model, cause = "elpe_ret")

granger_tbl <- tibble::tibble(
  Direction  = c("Oil → ELPE", "Market → ELPE", "ELPE → Oil"),
  `F-stat`   = round(c(gc_oil_elpe$Granger$statistic,
                        gc_mkt_elpe$Granger$statistic,
                        gc_elpe_oil$Granger$statistic), 4),
  `p-value`  = round(c(gc_oil_elpe$Granger$p.value,
                        gc_mkt_elpe$Granger$p.value,
                        gc_elpe_oil$Granger$p.value), 4),
  Conclusion = ifelse(c(gc_oil_elpe$Granger$p.value,
                         gc_mkt_elpe$Granger$p.value,
                         gc_elpe_oil$Granger$p.value) < 0.05,
                      "Granger-causes", "No Granger causality")
)

irf_obj <- vars::irf(var_model, impulse = "oil_ret", response = "elpe_ret",
                     n.ahead = 20, boot = TRUE, ci = 0.95, runs = 200)
irf_df <- tibble::tibble(
  Horizon = 0:20,
  IRF     = as.numeric(irf_obj$irf$oil_ret),
  Lower   = as.numeric(irf_obj$Lower$oil_ret),
  Upper   = as.numeric(irf_obj$Upper$oil_ret)
)
message(sprintf("H.3 VAR(%d) | Oil->ELPE p=%.4f | Mkt->ELPE p=%.4f",
                opt_lag, gc_oil_elpe$Granger$p.value, gc_mkt_elpe$Granger$p.value))

# ---- H.4 STRUCTURAL BREAK (Chow + Bai-Perron) -------------
chow_res <- strucchange::sctest(elpe_ret ~ mkt_ret, type = "Chow",
                                data = returns, point = round(nrow(returns) * 0.5))
bp_res   <- strucchange::breakpoints(elpe_ret ~ mkt_ret, data = returns)
n_breaks <- if (!all(is.na(bp_res$breakpoints))) length(bp_res$breakpoints) else 0L
bp_dates_vec <- if (n_breaks > 0) as.character(returns$Date[bp_res$breakpoints]) else "None detected"

break_tbl <- tibble::tibble(
  Test         = c("Chow (midpoint split)",
                   paste0("Bai-Perron break #", seq_len(max(1L, n_breaks)))),
  `p-value`    = round(c(chow_res$p.value, rep(NA_real_, max(1L, n_breaks))), 4),
  `Break date` = c(NA_character_,
                   if (n_breaks > 0) bp_dates_vec else "None"),
  Conclusion   = c(
    ifelse(chow_res$p.value < 0.05, "Structural break present", "No break at midpoint"),
    rep(if (n_breaks > 0) "Break detected" else "No break", max(1L, n_breaks))
  )
)

# CUSUM time series for chart
cusum_ef <- strucchange::efp(elpe_ret ~ mkt_ret, data = returns, type = "OLS-CUSUM")
cusum_df <- tibble::tibble(
  Date  = as.Date(returns$Date[seq_along(cusum_ef$process)]),
  CUSUM = as.numeric(cusum_ef$process)
)
message(sprintf("H.4 Structural break: Chow p=%.4f | BP breaks=%d",
                chow_res$p.value, n_breaks))

# ---- H.5 COINTEGRATION (Johansen trace test) ---------------
elpe_brent_lvl <- elpe_raw |>
  inner_join(
    brent_raw |> transmute(Date, brent = as.numeric(.data[[oil_col]])),
    by = "Date"
  ) |>
  filter(is.finite(elpe_close), is.finite(brent)) |>
  as.data.frame()

jo_mat  <- as.matrix(elpe_brent_lvl[, c("elpe_close", "brent")])
jo_test <- urca::ca.jo(jo_mat, type = "trace", K = 2, ecdet = "const")

jo_tbl <- tibble::tibble(
  `H₀`        = c("r = 0  (no cointegration)", "r ≤ 1"),
  `Trace stat` = round(rev(jo_test@teststat), 4),
  `CV 10%`     = round(jo_test@cval[, "10pct"], 4),
  `CV 5%`      = round(jo_test@cval[, "5pct"],  4),
  `CV 1%`      = round(jo_test@cval[, "1pct"],  4),
  Conclusion   = ifelse(rev(jo_test@teststat) > jo_test@cval[, "5pct"],
                        "Reject H₀ — cointegrated", "Fail to reject H₀")
)
message("H.5 Johansen cointegration complete.")

# ---- H.6 RISK-ADJUSTED PERFORMANCE -------------------------
ann_elpe_ret <- mean(returns$elpe_ret, na.rm = TRUE) * 252
ann_elpe_sd  <- sd(returns$elpe_ret,   na.rm = TRUE) * sqrt(252)
sharpe_elpe  <- (ann_elpe_ret - rf)  / ann_elpe_sd
treynor_elpe <- (ann_elpe_ret - rf)  / beta_ols
jensen_alpha <- alpha_ols * 252

ann_moh_ret  <- mean(moh_df$moh_ret, na.rm = TRUE) * 252
ann_moh_sd   <- sd(moh_df$moh_ret,   na.rm = TRUE) * sqrt(252)
sharpe_moh   <- (ann_moh_ret - rf) / ann_moh_sd
treynor_moh  <- (ann_moh_ret - rf) / beta_moh_static

perf_tbl <- tibble::tibble(
  Metric = c("Annualised Return (%)", "Annualised Volatility (%)",
             "Sharpe Ratio", "Treynor Ratio", "Jensen's α (ann.)"),
  ELPE   = round(c(ann_elpe_ret * 100, ann_elpe_sd * 100,
                    sharpe_elpe, treynor_elpe, jensen_alpha), 4),
  MOH    = round(c(ann_moh_ret  * 100, ann_moh_sd  * 100,
                    sharpe_moh,  treynor_moh,  NA_real_),  4),
  Formula = c("mean(r)×252", "sd(r)×√252",
              "(Rₐ−Rf)/σₐ", "(Rₐ−Rf)/β",  "α_OLS×252 (ELPE only)")
)

# ---- H.7 VaR & CVaR (Markov DCF + GBM) --------------------
var95_dcf  <- as.numeric(quantile(pricepershare, 0.05, na.rm = TRUE))
cvar95_dcf <- mean(pricepershare[pricepershare <= var95_dcf], na.rm = TRUE)
var99_dcf  <- as.numeric(quantile(pricepershare, 0.01, na.rm = TRUE))
cvar99_dcf <- mean(pricepershare[pricepershare <= var99_dcf], na.rm = TRUE)

var95_gbm  <- as.numeric(quantile(gbmfinal, 0.05, na.rm = TRUE))
cvar95_gbm <- mean(gbmfinal[gbmfinal <= var95_gbm], na.rm = TRUE)
var99_gbm  <- as.numeric(quantile(gbmfinal, 0.01, na.rm = TRUE))
cvar99_gbm <- mean(gbmfinal[gbmfinal <= var99_gbm], na.rm = TRUE)

risk_tbl <- tibble::tibble(
  Metric                = c("VaR 95%", "CVaR / ES 95%", "VaR 99%", "CVaR / ES 99%"),
  `Markov DCF (EUR/sh)` = round(c(var95_dcf, cvar95_dcf, var99_dcf, cvar99_dcf), 2),
  `GBM 1Y (EUR/sh)`     = round(c(var95_gbm, cvar95_gbm, var99_gbm, cvar99_gbm), 2),
  `DCF vs current`      = paste0(round((c(var95_dcf, cvar95_dcf, var99_dcf, cvar99_dcf) /
                                          currentprice - 1) * 100, 1), "%"),
  `GBM vs current`      = paste0(round((c(var95_gbm, cvar95_gbm, var99_gbm, cvar99_gbm) /
                                          currentprice - 1) * 100, 1), "%")
)

# ---- H.8 EVA (Economic Value Added, projected 2026–2030) ---
ic_base <- initialequity + initialdebt
ic_proj <- ic_base * (baserevenue / baserevenue[1])

eva_tbl <- tibble::tibble(
  Year                         = 2026:2030,
  `NOPAT (EUR k)`              = round(basenopat, 0),
  `Inv. Capital (EUR k)`       = round(ic_proj, 0),
  `ROIC (%)`                   = paste0(round(basenopat / ic_proj * 100, 2), "%"),
  `WACC×IC — BU ★`             = round(ic_proj * wacc_bu, 0),
  `EVA — BU ★ (EUR k)`         = round(basenopat - ic_proj * wacc_bu, 0),
  `ROIC–WACC BU (pp)`          = paste0(round((basenopat / ic_proj - wacc_bu) * 100, 2), " pp"),
  `WACC×IC — OLS`              = round(ic_proj * wacc_ols, 0),
  `EVA — OLS (EUR k)`          = round(basenopat - ic_proj * wacc_ols, 0),
  `ROIC–WACC OLS (pp)`         = paste0(round((basenopat / ic_proj - wacc_ols) * 100, 2), " pp")
)

# ---- H.9 REVERSE DCF ----------------------------------------
rev_dcf_fn <- function(target_px, w, fcff_v, nd, sh) {
  obj <- function(g) {
    if (g >= w) return(1e9)
    pv <- sum(fcff_v / (1 + w)^seq_along(fcff_v))
    tv <- fcff_v[length(fcff_v)] * (1 + g) / (w - g) / (1 + w)^length(fcff_v)
    (pv + tv - nd) * 1000 / sh - target_px
  }
  tryCatch(uniroot(obj, c(-0.05, w - 0.001))$root, error = function(e) NA_real_)
}

# Bottom-up WACC (base case)
imp_g_base_bu <- rev_dcf_fn(currentprice,        wacc_bu, basefcff, netdebtlatest, shares)
imp_g_bear_bu <- rev_dcf_fn(currentprice * 0.80, wacc_bu, basefcff, netdebtlatest, shares)
imp_g_bull_bu <- rev_dcf_fn(currentprice * 1.20, wacc_bu, basefcff, netdebtlatest, shares)

# OLS WACC (robustness)
imp_g_base_ols <- rev_dcf_fn(currentprice,        wacc_ols, basefcff, netdebtlatest, shares)
imp_g_bear_ols <- rev_dcf_fn(currentprice * 0.80, wacc_ols, basefcff, netdebtlatest, shares)
imp_g_bull_ols <- rev_dcf_fn(currentprice * 1.20, wacc_ols, basefcff, netdebtlatest, shares)

# Alias for legacy charts
imp_g_base <- imp_g_base_bu
imp_g_bear <- imp_g_bear_bu
imp_g_bull <- imp_g_bull_bu

reverse_tbl <- tibble::tibble(
  Scenario      = c("Current (€9.50)", "Bear (–20%)", "Bull (+20%)"),
  `Price (EUR)` = round(c(currentprice, currentprice * 0.80, currentprice * 1.20), 2),
  `Implied g — BU ★`  = paste0(round(c(imp_g_base_bu,  imp_g_bear_bu,  imp_g_bull_bu)  * 100, 2), "%"),
  `Δg vs 1.8% — BU`   = paste0(round((c(imp_g_base_bu,  imp_g_bear_bu,  imp_g_bull_bu)  - growthterminal) * 100, 2), " pp"),
  `Implied g — OLS`   = paste0(round(c(imp_g_base_ols, imp_g_bear_ols, imp_g_bull_ols) * 100, 2), "%"),
  `Δg vs 1.8% — OLS`  = paste0(round((c(imp_g_base_ols, imp_g_bear_ols, imp_g_bull_ols) - growthterminal) * 100, 2), " pp"),
  Interpretation = c(
    ifelse(!is.na(imp_g_base_bu) && imp_g_base_bu > growthterminal,
           "Market prices in above-base growth (BU WACC)", "Market prices in below-base growth (BU WACC)"),
    "Bear scenario needs lower-than-base g to justify",
    "Bull scenario supported only if higher g is realised"
  )
)
message(sprintf("H.9 Reverse DCF (BU):  implied g at €%.2f = %.2f%%", currentprice,
                ifelse(is.na(imp_g_base_bu),  NA, imp_g_base_bu  * 100)))
message(sprintf("H.9 Reverse DCF (OLS): implied g at €%.2f = %.2f%%", currentprice,
                ifelse(is.na(imp_g_base_ols), NA, imp_g_base_ols * 100)))

# ---- H.10 CORRELATED MONTE CARLO ----------------------------
cor_raw <- matrix(c(
  1.00, 0.60, 0.75, 0.70,
  0.60, 1.00, 0.65, 0.80,
  0.75, 0.65, 1.00, 0.90,
  0.70, 0.80, 0.90, 1.00
), nrow = 4L, byrow = TRUE)
chol_corr <- tryCatch(chol(cor_raw), error = function(e) diag(4L))

set.seed(999L)
nsim_corr   <- 2000L
noise_scale <- 0.12
corr_prices <- numeric(nsim_corr)

for (i in seq_len(nsim_corr)) {
  z     <- matrix(rnorm(yearssim * 4L), nrow = yearssim, ncol = 4L) %*% chol_corr
  fcf_c <- pmax(basefcff   * (1 + z[, 1] * noise_scale), basefcff   * 0.40)
  ebi_c <- pmax(baseebitda * (1 + z[, 4] * noise_scale), baseebitda * 0.40)
  pv_c  <- sum(fcf_c * discfactors)
  tv_g  <- fcf_c[yearssim] * (1 + growthterminal) /
              (waccval - growthterminal) / (1 + waccval)^yearssim
  tv_e  <- ebi_c[yearssim] * exitmultiple / (1 + waccval)^yearssim
  ev_c  <- pv_c + wgordon * tv_g + wevebitda * tv_e + totaloption
  corr_prices[i] <- (ev_c - netdebtlatest) * 1000 / shares
}

corr_var95  <- as.numeric(quantile(corr_prices, 0.05,  na.rm = TRUE))
corr_cvar95 <- mean(corr_prices[corr_prices <= corr_var95], na.rm = TRUE)

corr_compare_tbl <- tibble::tibble(
  Metric            = c("Mean price (EUR)", "Median price (EUR)",
                        "P05", "P95", "VaR 95%", "CVaR 95%", "Prob upside (%)"),
  `Correlated MC`   = round(c(mean(corr_prices), median(corr_prices),
                               as.numeric(quantile(corr_prices, 0.05)),
                               as.numeric(quantile(corr_prices, 0.95)),
                               corr_var95, corr_cvar95,
                               mean(corr_prices > currentprice) * 100), 3),
  `Base Markov DCF` = round(c(meanofmeans, meanofmedians,
                               avgp05, avgp95, var95_dcf, cvar95_dcf,
                               meanprobupside * 100), 3)
)

corr_dist_df <- dplyr::bind_rows(
  tibble::tibble(Model = "Correlated MC",   Price = corr_prices),
  tibble::tibble(Model = "Base Markov DCF", Price = runresults$meanprice)
)
message(sprintf("H.10 Correlated MC: mean=%.2f | Markov mean=%.2f",
                mean(corr_prices), meanofmeans))

# ============================================================
# G. SHINY DASHBOARD
# ============================================================

ui <- shinydashboard::dashboardPage(
  skin = "blue",
  shinydashboard::dashboardHeader(title = "HELLENiQ ENERGY Valuation Dashboard"),
  shinydashboard::dashboardSidebar(
    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("Altman Z-Score", tabName = "altman", icon = shiny::icon("shield-alt")),
      shinydashboard::menuItem("CAPM Diagnostics", tabName = "capm", icon = shiny::icon("chart-line")),
      shinydashboard::menuItem("Rolling Beta", tabName = "rollbeta", icon = shiny::icon("sync")),
      shinydashboard::menuItem("WACC Scenarios", tabName = "wacc", icon = shiny::icon("percentage")),
      shinydashboard::menuItem("WACC Comparison", tabName = "wacc_comp", icon = shiny::icon("balance-scale")),
      shinydashboard::menuItem("DuPont Analysis", tabName = "dupont", icon = shiny::icon("layer-group")),
      shinydashboard::menuItem("Financial Ratios", tabName = "ratios", icon = shiny::icon("table")),
      shinydashboard::menuItem("Sensitivity WACC-g", tabName = "sensitivity", icon = shiny::icon("sliders-h")),
      shinydashboard::menuItem("Sentiment", tabName = "sentiment", icon = shiny::icon("comments")),
      shinydashboard::menuItem("PCA Biplot", tabName = "pca", icon = shiny::icon("project-diagram")),
      shinydashboard::menuItem("Crack Spread", tabName = "crack", icon = shiny::icon("oil-can")),
      shinydashboard::menuItem("Markov DCF Valuation", tabName = "dcf", icon = shiny::icon("calculator")),
      shinydashboard::menuItem("AFN Robustness", tabName = "afn", icon = shiny::icon("shield-alt")),
      shinydashboard::menuItem("Interactive WACC-g", tabName = "inter", icon = shiny::icon("sliders")),
      shinydashboard::menuItem("GBM Monte Carlo", tabName = "gbm", icon = shiny::icon("dice")),
      shinydashboard::menuItem("Downloads", tabName = "downloads", icon = shiny::icon("download")),
      shinydashboard::menuItem("── Advanced Econometrics ──", tabName = "sep1", icon = shiny::icon("flask")),
      shinydashboard::menuItem("Stationarity (ADF/KPSS)", tabName = "stationary", icon = shiny::icon("wave-square")),
      shinydashboard::menuItem("GARCH(1,1) Volatility",    tabName = "garch",     icon = shiny::icon("chart-area")),
      shinydashboard::menuItem("VAR & Granger Causality",  tabName = "vargranger",icon = shiny::icon("exchange-alt")),
      shinydashboard::menuItem("Cointegration (Johansen)", tabName = "coint",     icon = shiny::icon("link")),
      shinydashboard::menuItem("Structural Break",         tabName = "break",     icon = shiny::icon("bolt")),
      shinydashboard::menuItem("Risk-Adjusted Performance",tabName = "perf",      icon = shiny::icon("trophy")),
      shinydashboard::menuItem("VaR & CVaR",               tabName = "risk",      icon = shiny::icon("shield-alt")),
      shinydashboard::menuItem("EVA & Reverse DCF",        tabName = "eva",       icon = shiny::icon("balance-scale")),
      shinydashboard::menuItem("Correlated MC",            tabName = "corrmc",    icon = shiny::icon("random"))
    )
  ),
  shinydashboard::dashboardBody(
    shinydashboard::tabItems(
      shinydashboard::tabItem(tabName = "altman",
                              shiny::fluidRow(
                                shinydashboard::box(
                                  title = "Altman Z-Score (ELPE vs MOH)", width = 12, status = "primary", solidHeader = TRUE,
                                  shiny::p(class = "text-muted", "Source: Excel sheet 7_ALTMAN. Caution zone: 1.81–2.99."),
                                  DT::DTOutput("altmantable"),
                                  shiny::br(),
                                  plotly::plotlyOutput("altmanchart", height = "420px"),
                                  shiny::br(),
                                  shiny::downloadButton("dlaltman2", "Download CSV")
                                )
                              )
      ),
      shinydashboard::tabItem(tabName = "capm",
                              shiny::fluidRow(
                                shinydashboard::box(title = "OLS Regression Diagnostics", width = 12, status = "primary", solidHeader = TRUE,
                                                    DT::DTOutput("diagtable"))
                              ),
                              shiny::fluidRow(
                                shinydashboard::box(title = sprintf("CAPM Scatter: β=%.4f α=%.6f R²=%.3f", beta_ols, alpha_ols, r2_capm),
                                                    width = 8, status = "info", solidHeader = TRUE,
                                                    plotly::plotlyOutput("scatterplot", height = "420px")),
                                shinydashboard::box(title = "Multi-Factor Model (Market + Brent)", width = 4, status = "info", solidHeader = TRUE,
                                                    shiny::verbatimTextOutput("mfsummary"))
                              )
      ),
      shinydashboard::tabItem(tabName = "rollbeta",
                              shiny::fluidRow(
                                shinydashboard::box(title = "ELPE Rolling Beta (250-day)", width = 12, status = "warning", solidHeader = TRUE,
                                                    shiny::p(class = "text-muted", "Dashed = static OLS beta. Dotted = rolling mean."),
                                                    plotly::plotlyOutput("betaelpe", height = "370px"),
                                                    shiny::br(),
                                                    plotly::plotlyOutput("betacompare", height = "370px"))
                              )
      ),
      shinydashboard::tabItem(tabName = "wacc",
        shiny::fluidRow(
          shinydashboard::box(title = "Adjust WACC Inputs", width = 3, status = "primary", solidHeader = TRUE,
            shiny::sliderInput("rf_in",   "Risk-Free Rate Rf (%)",      min = 0.5,  max = 6.0,  value = 3.08, step = 0.01),
            shiny::sliderInput("erp_in",  "Equity Risk Premium ERP (%)",min = 3.0,  max = 12.0, value = 7.52, step = 0.01),
            shiny::sliderInput("beta_in", "Levered Beta",               min = 0.20, max = 2.00, value = round(beta_bu, 3), step = 0.01),
            shiny::sliderInput("kd_in",   "Pre-Tax Cost of Debt Kd (%)",min = 1.0,  max = 8.0,  value = 4.54, step = 0.01),
            shiny::sliderInput("we_in",   "Weight of Equity We (%)",    min = 20.0, max = 80.0, value = 48.8, step = 0.1),
            shiny::sliderInput("tax_in",  "Tax Rate (%)",               min = 10.0, max = 35.0, value = 22.0, step = 0.5)
          ),
          shinydashboard::box(title = "Live WACC Computation", width = 5, status = "success", solidHeader = TRUE,
            DT::DTOutput("wacc_live_tbl"),
            shiny::br(),
            plotly::plotlyOutput("wacc_decomp_chart", height = "260px")
          ),
          shinydashboard::box(title = "Live Sensitivity (WACC–g spread)", width = 4, status = "warning", solidHeader = TRUE,
            DT::DTOutput("wacc_sens_live", height = "360px")
          )
        ),
        shiny::fluidRow(
          shinydashboard::box(title = "Base-Case WACC Summary (Apr-2026)", width = 6, status = "info", solidHeader = TRUE,
            DT::DTOutput("wacctable")),
          shinydashboard::box(title = "Scenario Analysis", width = 6, status = "info", solidHeader = TRUE,
            DT::DTOutput("scenariotable"),
            shiny::br(),
            plotly::plotlyOutput("waccchart", height = "250px"))
        )
      ),
      # ── WACC COMPARISON TAB ────────────────────────────────
      shinydashboard::tabItem(tabName = "wacc_comp",
        shiny::fluidRow(
          shinydashboard::box(
            title = "OLS vs Bottom-Up WACC — Side-by-Side", width = 12,
            status = "primary", solidHeader = TRUE,
            shiny::p(class = "text-muted",
              "Bottom-up (Damodaran pure-play segments) is the base case for all DCF/EVA calculations.
               OLS-Hamada is retained as a robustness check. The ~1.45 pp gap is driven entirely by beta methodology."),
            DT::DTOutput("wacc_comp_tbl"),
            shiny::br(),
            plotly::plotlyOutput("wacc_comp_chart", height = "320px")
          )
        ),
        shiny::fluidRow(
          shinydashboard::box(
            title = "Segment Betas (Damodaran Jan-2026)", width = 7,
            status = "info", solidHeader = TRUE,
            shiny::p(class = "text-muted",
              "Revenue-weighted composite βU = 0.5219. Re-levered at actual D/E = 1.099x → βL* = 0.969."),
            DT::DTOutput("seg_beta_tbl")
          ),
          shinydashboard::box(
            title = "Beta Build-up Summary", width = 5,
            status = "warning", solidHeader = TRUE,
            DT::DTOutput("beta_buildup_tbl")
          )
        )
      ),

      shinydashboard::tabItem(tabName = "dupont",
                              shiny::fluidRow(
                                shinydashboard::box(title = "DuPont ROE (ELPE vs MOH, 2021–2025)", width = 12, status = "primary", solidHeader = TRUE,
                                                    shiny::p(class = "text-muted", "ROE = Net Margin × Turnover × Leverage."),
                                                    plotly::plotlyOutput("dupontroe", height = "370px"))
                              ),
                              shiny::fluidRow(
                                shinydashboard::box(title = "Components Breakdown (FY 2025)", width = 6, status = "info", solidHeader = TRUE,
                                                    plotly::plotlyOutput("dupontcomp", height = "370px")),
                                shinydashboard::box(title = "Decomposition Table", width = 6, status = "info", solidHeader = TRUE,
                                                    DT::DTOutput("duponttable"))
                              )
      ),
      shinydashboard::tabItem(tabName = "ratios",
                              shiny::fluidRow(
                                shinydashboard::box(title = "Financial Ratios", width = 12, status = "primary", solidHeader = TRUE,
                                                    DT::DTOutput("ratiostable"),
                                                    shiny::br(),
                                                    shiny::downloadButton("dlratios", "Download CSV"))
                              )
      ),
      shinydashboard::tabItem(tabName = "sensitivity",
                              shiny::fluidRow(
                                shinydashboard::box(title = "Terminal Value Spread (WACC - g)", width = 12, status = "warning", solidHeader = TRUE,
                                                    shiny::p(class = "text-muted", "Base case: WACC 5.68%, g 1.80%."),
                                                    DT::DTOutput("senstable"))
                              )
      ),
      shinydashboard::tabItem(tabName = "sentiment",
                              shiny::fluidRow(
                                shinydashboard::box(title = "Net Sentiment Score", width = 12, status = "primary", solidHeader = TRUE,
                                                    plotly::plotlyOutput("sentimentchart", height = "300px"))
                              )
      ),
      shinydashboard::tabItem(tabName = "pca",
                              shiny::fluidRow(
                                shinydashboard::box(title = "PCA Biplot (ELPE vs MOH)", width = 12, status = "primary", solidHeader = TRUE,
                                                    shiny::p(class = "text-muted", "Variables: Z, ROE, EBITDA margin, D/E, Current Ratio, Net Debt/EBITDA."),
                                                    shiny::plotOutput("pcaplot", height = "520px"))
                              )
      ),
      shinydashboard::tabItem(tabName = "crack",
                              shiny::fluidRow(
                                shinydashboard::box(title = "Crack Spread Regression", width = 12, status = "warning", solidHeader = TRUE,
                                                    shiny::p(class = "text-muted", "OLS regressions linking HO-WTI crack spread to ELPE fundamentals."),
                                                    shiny::fluidRow(
                                                      shiny::column(6, plotly::plotlyOutput("reggrowthplot", height = "400px")),
                                                      shiny::column(6, plotly::plotlyOutput("regmarginplot", height = "400px"))
                                                    ),
                                                    shiny::br(),
                                                    shiny::fluidRow(
                                                      shiny::column(6, DT::DTOutput("reggrowthtbl")),
                                                      shiny::column(6, DT::DTOutput("regmargintbl"))
                                                    ),
                                                    shiny::br(),
                                                    shiny::downloadButton("dlcrack", "Download CSV"))
                              )
      ),
      shinydashboard::tabItem(tabName = "dcf",
                              shiny::fluidRow(
                                shinydashboard::valueBoxOutput("dcfmeanbox", width = 3),
                                shinydashboard::valueBoxOutput("dcfmedianbox", width = 3),
                                shinydashboard::valueBoxOutput("dcfupsidebox", width = 3),
                                shinydashboard::valueBoxOutput("dcfrecbox", width = 3)
                              ),
                              shiny::fluidRow(
                                shinydashboard::box(title = "Model Comparison", width = 4, status = "primary", solidHeader = TRUE,
                                                    DT::DTOutput("comparisontbl")),
                                shinydashboard::box(title = "Price Distribution", width = 8, status = "info", solidHeader = TRUE,
                                                    plotly::plotlyOutput("pricehist", height = "380px"),
                                                    shiny::br(),
                                                    shiny::downloadButton("dldcf", "Download CSV"))
                              )
      ),
      shinydashboard::tabItem(tabName = "afn",
                              shiny::fluidRow(
                                shinydashboard::box(title = "Covenant Stress Test", width = 12, status = "danger", solidHeader = TRUE,
                                                    shiny::p(class = "text-muted", "1,000 paths over 5 years. Single breach = path fails."),
                                                    shiny::fluidRow(
                                                      shiny::column(6, DT::DTOutput("criteriatbl")),
                                                      shiny::column(6, plotly::plotlyOutput("criteriabar", height = "260px"))
                                                    ),
                                                    shiny::fluidRow(
                                                      shinydashboard::box(title = "DSCR Distribution by Year", width = 6, status = "info", solidHeader = TRUE,
                                                                          plotly::plotlyOutput("afndscr", height = "380px")),
                                                      shinydashboard::box(title = "D/E Ratio Fan Chart (p10-p25-p75-p90)", width = 6, status = "info", solidHeader = TRUE,
                                                                          plotly::plotlyOutput("afnde", height = "380px"))
                                                    ),
                                                    shiny::fluidRow(
                                                      shinydashboard::box(title = "AFN Distribution (EUR M)", width = 6, status = "warning", solidHeader = TRUE,
                                                                          plotly::plotlyOutput("afnafn", height = "380px")),
                                                      shinydashboard::box(title = "Net Debt/EBITDA Fan Chart", width = 6, status = "warning", solidHeader = TRUE,
                                                                          plotly::plotlyOutput("afnnde", height = "380px"))
                                                    )
                                )
                              )
      ),
      shinydashboard::tabItem(tabName = "inter",
                              shiny::fluidRow(
                                shinydashboard::box(title = "Fair Value Sensitivity", width = 3, status = "primary", solidHeader = TRUE,
                                                    shiny::sliderInput("waccinter", "WACC %", min = 4.0, max = 10.0, value = round(wacc_bu * 100, 1), step = 0.1),
                                                    shiny::sliderInput("ginter", "Terminal growth %", min = 0.0, max = 3.0, value = 1.8, step = 0.1),
                                                    shiny::actionButton("calcinter", "Refresh Table", class = "btn-primary btn-block")),
                                shinydashboard::box(title = "Breakdown (EUR M)", width = 5, status = "info", solidHeader = TRUE,
                                                    shiny::tableOutput("intertable")),
                                shinydashboard::box(title = "Fair Value per Share", width = 4, status = "success", solidHeader = TRUE,
                                                    shiny::verbatimTextOutput("interprice"))
                              )
      ),
      shinydashboard::tabItem(tabName = "gbm",
                              shiny::fluidRow(
                                shinydashboard::valueBoxOutput("gbmmeanbox", width = 4),
                                shinydashboard::valueBoxOutput("gbmlobox", width = 4),
                                shinydashboard::valueBoxOutput("gbmhibox", width = 4)
                              ),
                              shiny::fluidRow(
                                shinydashboard::box(title = "GBM Price Paths (1Y, 1,000 paths)", width = 12, status = "success", solidHeader = TRUE,
                                                    plotly::plotlyOutput("gbmchart", height = "400px"))
                              )
      ),
      shinydashboard::tabItem(tabName = "downloads",
        shiny::fluidRow(
          shinydashboard::box(title = "Download All Outputs", width = 12, status = "primary", solidHeader = TRUE,
            shiny::p("Click the buttons below to download CSV files of each analysis."),
            shiny::downloadButton("dlaltman2", "Altman CSV"),
            shiny::downloadButton("dlcapm",    "CAPM Table"),
            shiny::downloadButton("dlratios2", "Financial Ratios"),
            shiny::downloadButton("dldcf2",    "DCF Valuation"),
            shiny::downloadButton("dlwacc",    "WACC Table"),
            shiny::downloadButton("dlcrack2",  "Crack Spread Data"))
        )
      ),

      # ── ADVANCED ECONOMETRICS TABS ──────────────────────────
      shinydashboard::tabItem(tabName = "sep1"),

      shinydashboard::tabItem(tabName = "stationary",
        shiny::fluidRow(
          shinydashboard::box(title = "Stationarity Tests: ADF & KPSS", width = 12,
                              status = "primary", solidHeader = TRUE,
            shiny::p(class = "text-muted",
              "ADF H₀: unit root (non-stationary). KPSS H₀: stationary. *** p<0.01 ** p<0.05 * p<0.10."),
            DT::DTOutput("stat_tbl"),
            shiny::br(),
            shiny::p(class = "text-muted",
              "Stationarity is a prerequisite for valid OLS regressions and VAR estimation.
               If crack spread is non-stationary in levels, proceed to cointegration analysis (Johansen tab).")
          )
        )
      ),

      shinydashboard::tabItem(tabName = "garch",
        shiny::fluidRow(
          shinydashboard::box(title = "GARCH(1,1) — Student-t: Parameter Estimates", width = 5,
                              status = "primary", solidHeader = TRUE,
            DT::DTOutput("garch_param_tbl")),
          shinydashboard::box(title = "GARCH-implied Annualised Conditional Volatility vs. Static GBM σ",
                              width = 7, status = "info", solidHeader = TRUE,
            plotly::plotlyOutput("garch_vol_chart", height = "360px"))
        )
      ),

      shinydashboard::tabItem(tabName = "vargranger",
        shiny::fluidRow(
          shinydashboard::box(title = "Granger Causality Tests", width = 4,
                              status = "primary", solidHeader = TRUE,
            shiny::p(class = "text-muted", "H₀: X does not Granger-cause Y."),
            DT::DTOutput("granger_tbl")),
          shinydashboard::box(title = "Impulse Response: Oil shock → ELPE returns (20-day horizon)",
                              width = 8, status = "info", solidHeader = TRUE,
            shiny::p(class = "text-muted", "95% bootstrap CI. Shaded band = uncertainty around point estimate."),
            plotly::plotlyOutput("irf_chart", height = "360px"))
        )
      ),

      shinydashboard::tabItem(tabName = "coint",
        shiny::fluidRow(
          shinydashboard::box(title = "Johansen Trace Test: ELPE Price vs Brent Price", width = 12,
                              status = "primary", solidHeader = TRUE,
            shiny::p(class = "text-muted",
              "Tests whether ELPE and Brent prices share a long-run equilibrium relationship.
               Reject r=0 at 5% → at least one cointegrating vector exists → ECM is appropriate."),
            DT::DTOutput("jo_tbl"),
            shiny::br(),
            plotly::plotlyOutput("coint_levels_chart", height = "340px"))
        )
      ),

      shinydashboard::tabItem(tabName = "break",
        shiny::fluidRow(
          shinydashboard::box(title = "Structural Break Tests", width = 4,
                              status = "danger", solidHeader = TRUE,
            shiny::p(class = "text-muted", "Chow: midpoint split. Bai-Perron: data-driven."),
            DT::DTOutput("break_tbl")),
          shinydashboard::box(title = "OLS-CUSUM Process (boundary = 5% significance)",
                              width = 8, status = "warning", solidHeader = TRUE,
            shiny::p(class = "text-muted",
              "Empirical fluctuation process. Boundary crossing indicates parameter instability."),
            plotly::plotlyOutput("cusum_chart", height = "360px"))
        )
      ),

      shinydashboard::tabItem(tabName = "perf",
        shiny::fluidRow(
          shinydashboard::box(title = "Risk-Adjusted Performance: ELPE vs MOH", width = 6,
                              status = "primary", solidHeader = TRUE,
            shiny::p(class = "text-muted", "All metrics annualised on 252 trading days. Rf = 3.08% (ECB AAA)."),
            DT::DTOutput("perf_tbl")),
          shinydashboard::box(title = "Sharpe & Treynor Comparison", width = 6,
                              status = "info", solidHeader = TRUE,
            plotly::plotlyOutput("perf_chart", height = "340px"))
        )
      ),

      shinydashboard::tabItem(tabName = "risk",
        shiny::fluidRow(
          shinydashboard::valueBoxOutput("var95_dcf_box",  width = 3),
          shinydashboard::valueBoxOutput("cvar95_dcf_box", width = 3),
          shinydashboard::valueBoxOutput("var95_gbm_box",  width = 3),
          shinydashboard::valueBoxOutput("cvar95_gbm_box", width = 3)
        ),
        shiny::fluidRow(
          shinydashboard::box(title = "VaR & CVaR Table", width = 5,
                              status = "danger", solidHeader = TRUE,
            shiny::p(class = "text-muted",
              "VaR = worst case at confidence level. CVaR = average loss beyond VaR (Expected Shortfall)."),
            DT::DTOutput("risk_tbl")),
          shinydashboard::box(title = "Price Distribution with VaR/CVaR overlays", width = 7,
                              status = "warning", solidHeader = TRUE,
            plotly::plotlyOutput("risk_dist_chart", height = "380px"))
        )
      ),

      shinydashboard::tabItem(tabName = "eva",
        shiny::fluidRow(
          shinydashboard::box(title = "Economic Value Added (EVA) — Projected 2026–2030", width = 6,
                              status = "primary", solidHeader = TRUE,
            shiny::p(class = "text-muted", "EVA = NOPAT − WACC × Invested Capital. Positive → value creation."),
            DT::DTOutput("eva_tbl"),
            shiny::br(),
            plotly::plotlyOutput("eva_chart", height = "300px")),
          shinydashboard::box(title = "Reverse DCF — Implied Terminal Growth Rate", width = 6,
                              status = "info", solidHeader = TRUE,
            shiny::p(class = "text-muted",
              "Back-solves for the growth rate implied by the current market price.
               Compares to base-case assumption of 1.80%."),
            DT::DTOutput("reverse_tbl"),
            shiny::br(),
            plotly::plotlyOutput("reverse_chart", height = "260px"))
        )
      ),

      shinydashboard::tabItem(tabName = "corrmc",
        shiny::fluidRow(
          shinydashboard::box(title = "Correlated Monte Carlo vs Base Markov DCF", width = 5,
                              status = "primary", solidHeader = TRUE,
            shiny::p(class = "text-muted",
              "Correlated MC applies Cholesky-decomposed shocks so FCF, Revenue, NOPAT and EBITDA
               move together (ρ calibrated to refining economics). Removes independence assumption."),
            DT::DTOutput("corr_compare_tbl")),
          shinydashboard::box(title = "Overlaid Price Distributions (Correlated MC vs Base Markov)", width = 7,
                              status = "info", solidHeader = TRUE,
            plotly::plotlyOutput("corr_dist_chart", height = "400px"))
        )
      )
    )
  )
)

server <- function(input, output, session) {
  
  output$altmantable <- DT::renderDT({
    DT::datatable(altmanclean, options = list(dom = "t"), rownames = FALSE)
  })
  
  output$altmanchart <- plotly::renderPlotly({
    p <- ggplot2::ggplot(altmanlong, ggplot2::aes(x = Year, y = ZScore, colour = Company)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::geom_point(size = 3) +
      ggplot2::geom_hline(yintercept = 2.99, linetype = "dashed", colour = "darkgreen") +
      ggplot2::geom_hline(yintercept = 1.81, linetype = "dashed", colour = "red") +
      ggplot2::scale_colour_manual(values = c("ELPE" = "#2E86AB", "MOH" = "#E69F00")) +
      ggplot2::labs(title = "Altman Z-Score (2021–2025)", y = "Z-Score", x = NULL) +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  output$diagtable <- DT::renderDT({
    DT::datatable(diag_tbl, options = list(dom = "t"), rownames = FALSE)
  })

  output$scatterplot <- plotly::renderPlotly({
    p <- ggplot2::ggplot(returns, ggplot2::aes(x = mkt_ret, y = elpe_ret)) +
      ggplot2::geom_point(alpha = 0.45, size = 1.5, colour = "#2E86AB") +
      ggplot2::geom_smooth(method = "lm", se = TRUE, colour = "#E74C3C", linewidth = 1.2) +
      ggplot2::labs(title = sprintf("ELPE vs GD: β=%.4f, R²=%.3f", beta_ols, r2_capm),
                    x = "GD log-return", y = "ELPE log-return") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  output$mfsummary <- shiny::renderPrint({
    cat("Two-factor model (Market + Brent)\n")
    cat(sprintf("β_market = %.4f | β_oil = %.4f | R² = %.4f\n", coef(model_mf)["mkt_ret"], beta_oil, r2_mf))
    cat(sprintf("Nested F-test: Brent added, p = %.4f\n", nested_test$`Pr(>F)`[2]))
  })
  
  output$betaelpe <- plotly::renderPlotly({
    p <- ggplot2::ggplot(rolling_beta_csv, ggplot2::aes(x = as.Date(Date), y = beta_1y)) +
      ggplot2::geom_line(colour = "#2E86AB", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = beta_ols, linetype = "dashed", colour = "#E74C3C") +
      ggplot2::geom_hline(yintercept = mean_roll, linetype = "dotted", colour = "#28B463") +
      ggplot2::labs(title = sprintf("ELPE Rolling 1Y Beta: Static=%.4f, Roll mean=%.4f, Blume adj=%.4f",
                                    beta_ols, mean_roll, beta_blume),
                    x = NULL, y = "Rolling Beta") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  output$betacompare <- plotly::renderPlotly({
    p <- ggplot2::ggplot(df_compare, ggplot2::aes(x = Date, y = beta_roll, colour = Company)) +
      ggplot2::geom_line(linewidth = 0.8) +
      ggplot2::scale_colour_manual(values = c("ELPE" = "#2E86AB", "MOH" = "#E69F00")) +
      ggplot2::labs(title = sprintf("Rolling Beta: ELPE (static %.4f) vs MOH (static %.4f)", beta_ols, beta_moh_static),
                    x = NULL, y = "Rolling Beta") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  output$wacctable <- DT::renderDT({
    DT::datatable(wacc_summary, options = list(dom = "t", pageLength = 15), rownames = FALSE)
  })
  
  output$scenariotable <- DT::renderDT({
    DT::datatable(scenario_wacc, options = list(dom = "t"), rownames = FALSE)
  })
  
  output$waccchart <- plotly::renderPlotly({
    wn <- as.numeric(gsub("%", "", as.data.frame(scenario_wacc)$WACC))
    plotly::plot_ly(x = scenario_wacc$Scenario, y = wn, type = "bar",
                    marker = list(color = c("#E74C3C", "#2E86AB", "#28B463"))) %>%
      plotly::layout(title = "WACC by Scenario", yaxis = list(title = "WACC %"), xaxis = list(title = ""))
  })
  
  output$dupontroe <- plotly::renderPlotly({
    p <- ggplot2::ggplot(dupontall, ggplot2::aes(x = Year, y = ROE * 100, colour = Company)) +
      ggplot2::geom_line(linewidth = 1.2) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_colour_manual(values = c("ELPE" = "#2E86AB", "MOH" = "#E69F00")) +
      ggplot2::labs(title = "DuPont ROE (2021–2025)", y = "ROE (%)", x = NULL) +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  output$dupontcomp <- plotly::renderPlotly({
    p <- ggplot2::ggplot(dupont2025, ggplot2::aes(x = Component, y = Value, fill = Company)) +
      ggplot2::geom_col(position = "dodge", width = 0.6) +
      ggplot2::scale_fill_manual(values = c("ELPE" = "#2E86AB", "MOH" = "#E69F00")) +
      ggplot2::labs(title = "DuPont Components (FY 2025)", x = NULL, y = NULL) +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  output$duponttable <- DT::renderDT({
    DT::datatable(
      dupontall %>%
        dplyr::mutate(
          Margin = paste0(round(Margin * 100, 2), "%"),
          Turnover = round(Turnover, 3),
          Leverage = round(Leverage, 3),
          ROE = paste0(round(ROE * 100, 2), "%")
        ) %>%
        dplyr::arrange(Company, Year),
      options = list(dom = "t", pageLength = 10), rownames = FALSE
    )
  })
  
  output$ratiostable <- DT::renderDT({
    DT::datatable(financialratios, options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })
  
  output$senstable <- DT::renderDT({
    DT::datatable(sens_tbl, options = list(dom = "t", scrollX = TRUE), rownames = FALSE)
  })
  
  output$sentimentchart <- plotly::renderPlotly({
    plotly::plot_ly(sentimentdf, x = ~Category, y = ~NetScore, type = "bar",
                    marker = list(color = c("#2E86AB", "#A9C46A"))) %>%
      plotly::layout(title = "Sentiment Net Score", yaxis = list(title = "Net Score"))
  })
  
  output$pcaplot <- shiny::renderPlot({
    factoextra::fviz_pca_biplot(
      pcaresult,
      label = "var",
      habillage = pcadf$Company,
      addEllipses = FALSE,
      palette = c("#2E86AB", "#E69F00"),
      repel = TRUE,
      title = "PCA Biplot: ELPE vs MOH (2021–2025)"
    )
  })
  
  output$reggrowthplot <- plotly::renderPlotly({
    fitline <- data.frame(spread = seq(min(reggrowthdata$spread) - 2, max(reggrowthdata$spread) + 2, length.out = 80))
    fitline$growth <- alphagrowth + betagrowth * fitline$spread
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(data = fitline, ggplot2::aes(x = spread, y = growth * 100), colour = "#2E86AB", linewidth = 1) +
      ggplot2::geom_point(data = reggrowthdata, ggplot2::aes(x = spread, y = revgrowth * 100, label = year),
                          colour = "#2E86AB", size = 3.5) +
      ggplot2::geom_point(data = predpts, ggplot2::aes(x = spread, y = growth * 100, colour = label), size = 5, shape = 18) +
      ggplot2::geom_text(data = predpts, ggplot2::aes(x = spread, y = growth * 100, label = label), vjust = -1.1, size = 3) +
      ggplot2::scale_colour_manual(values = c("Normal" = "#E69F00", "Bull" = "#28B463", "Bear" = "#E74C3C")) +
      ggplot2::labs(title = sprintf("Revenue Growth vs Crack Spread (R² = %.3f, n = 4)", r2growth),
                    x = "Avg Annual Crack Spread (USD/bbl)", y = "Revenue Growth (%)", colour = "Regime") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p, tooltip = c("x", "y", "label"))
  })
  
  output$regmarginplot <- plotly::renderPlotly({
    fitm <- data.frame(spread = seq(min(regdata$spread) - 2, max(regdata$spread) + 2, length.out = 80))
    fitm$margin <- alphareg + betareg * fitm$spread
    p <- ggplot2::ggplot() +
      ggplot2::geom_line(data = fitm, ggplot2::aes(x = spread, y = margin * 100), colour = "#9B59B6", linewidth = 1) +
      ggplot2::geom_point(data = regdata, ggplot2::aes(x = spread, y = margin * 100, label = year),
                          colour = "#9B59B6", size = 3.5) +
      ggplot2::geom_point(data = predmarginpts, ggplot2::aes(x = spread, y = margin * 100, colour = label), size = 5, shape = 18) +
      ggplot2::geom_text(data = predmarginpts, ggplot2::aes(x = spread, y = margin * 100, label = label), vjust = -1.1, size = 3) +
      ggplot2::scale_colour_manual(values = c("Normal" = "#E69F00", "Bull" = "#28B463", "Bear" = "#E74C3C")) +
      ggplot2::labs(title = sprintf("Adj. EBITDA Margin vs Crack Spread (R² = %.3f, n = 5)", summary(modelreg)$r.squared),
                    x = "Avg Annual Crack Spread (USD/bbl)", y = "Adj. EBITDA Margin (%)", colour = "Regime") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p, tooltip = c("x", "y", "label"))
  })
  
  output$reggrowthtbl <- DT::renderDT({
    tbl <- data.frame(
      Regime = c("Bear", "Normal", "Bull"),
      `Spread (USD/bbl)` = round(statemeans[c(3, 1, 2)], 2),
      `Rev Growth (%)` = round((alphagrowth + betagrowth * statemeans[c(3, 1, 2)]) * 100, 2),
      Note = rep("Illustrative (n = 4)", 3)
    )
    DT::datatable(tbl, options = list(dom = "t"), rownames = FALSE)
  })
  
  output$regmargintbl <- DT::renderDT({
    tbl <- data.frame(
      Regime = c("Bear", "Normal", "Bull"),
      `Spread (USD/bbl)` = round(statemeans[c(3, 1, 2)], 2),
      `EBITDA Margin (%)` = round((alphareg + betareg * statemeans[c(3, 1, 2)]) * 100, 2)
    )
    DT::datatable(tbl, options = list(dom = "t"), rownames = FALSE)
  })
  
  output$comparisontbl <- DT::renderDT({
    DT::datatable(comparisontbl, options = list(dom = "t"), rownames = FALSE)
  })
  
  output$pricehist <- plotly::renderPlotly({
    dfhist <- dplyr::bind_rows(
      tibble::tibble(Model = "Base Markov (run means)", Price = runresults$meanprice),
      tibble::tibble(Model = "Mobile Markov (run means)", Price = mobilerunresults$meanpricemobile)
    )
    p <- ggplot2::ggplot(dfhist, ggplot2::aes(x = Price, fill = Model)) +
      ggplot2::geom_histogram(bins = 35, alpha = 0.7, colour = "white", position = "identity") +
      ggplot2::geom_vline(xintercept = currentprice, colour = "#E74C3C", linewidth = 1.2) +
      ggplot2::geom_vline(xintercept = meanprice, colour = "#2E86AB", linewidth = 1.2, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = meanpricemobile, colour = "#E69F00", linewidth = 1.2, linetype = "dashed") +
      ggplot2::scale_fill_manual(values = c("Base Markov (run means)" = "#2E86AB", "Mobile Markov (run means)" = "#E69F00")) +
      ggplot2::labs(title = "Robust Markov DCF: Distribution Across Simulation Runs", x = "EUR/share", y = "Frequency") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  output$dcfmeanbox <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(paste0("€", round(meanprice, 2)), "Mean DCF Price", icon = shiny::icon("euro-sign"), color = "blue")
  })
  
  output$dcfmedianbox <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(paste0("€", round(medianprice, 2)), "Median DCF Price", icon = shiny::icon("chart-bar"), color = "light-blue")
  })
  
  output$dcfupsidebox <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(paste0(round(probupside * 100, 1), "%"), "Prob. Upside", icon = shiny::icon("arrow-up"),
                             color = if (probupside >= 0.6) "green" else if (probupside >= 0.4) "yellow" else "red")
  })
  
  output$dcfrecbox <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(rec, "Recommendation", icon = shiny::icon("thumbs-up"),
                             color = if (rec == "BUY") "green" else if (rec == "HOLD") "yellow" else "red")
  })
  
  output$criteriatbl <- DT::renderDT({
    critdf <- data.frame(
      Criterion = c("A: DSCR ≥1.3x (all years)", "B: D/E ≤0.79x (all years)",
                    "C: ND/EBITDA ≤2.5x (all years)", "All three satisfied"),
      `Probability (%)` = paste0(round(c(mean(criterionA), mean(criterionB), mean(criterionC),
                                         mean(criterionA & criterionB & criterionC)) * 100, 1), "%")
    )
    DT::datatable(critdf, options = list(dom = "t"), rownames = FALSE)
  })
  
  output$criteriabar <- plotly::renderPlotly({
    critdf <- data.frame(
      Label = c("DSCR ≥1.3x", "D/E ≤0.79x", "ND/EBITDA ≤2.5x", "All three"),
      Prob = round(c(mean(criterionA), mean(criterionB), mean(criterionC),
                     mean(criterionA & criterionB & criterionC)) * 100, 1)
    )
    plotly::plot_ly(critdf, x = ~Label, y = ~Prob, type = "bar",
                    marker = list(color = c("#2E86AB", "#2E86AB", "#28B463", "#E74C3C"))) %>%
      plotly::layout(title = "% of paths passing each covenant",
                     yaxis = list(title = "% of 1,000 paths", range = c(0, 100)),
                     xaxis = list(title = ""))
  })
  
  output$afndscr <- plotly::renderPlotly({
    dscrdf <- as.data.frame(dscrsim)
    colnames(dscrdf) <- paste0("Year ", 1:yearssim)
    dscrlong <- dscrdf %>%
      dplyr::mutate(sim = dplyr::row_number()) %>%
      tidyr::pivot_longer(-sim, names_to = "Year", values_to = "DSCR")
    p <- ggplot2::ggplot(dscrlong, ggplot2::aes(x = Year, y = DSCR)) +
      ggplot2::geom_boxplot(fill = "#2E86AB", alpha = 0.7, outlier.size = 0.4) +
      ggplot2::geom_hline(yintercept = 1.3, linetype = "dashed", colour = "red", linewidth = 0.8) +
      ggplot2::labs(title = "DSCR Distribution by Year (1,000 paths)", x = NULL, y = "DSCR") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  output$afnde <- plotly::renderPlotly({
    dedf <- as.data.frame(dematrix)
    colnames(dedf) <- paste0("Year ", 1:yearssim)
    defan <- dedf %>%
      dplyr::mutate(sim = dplyr::row_number()) %>%
      tidyr::pivot_longer(-sim, names_to = "Year", values_to = "DE") %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(
        p10 = quantile(DE, 0.10, na.rm = TRUE),
        p25 = quantile(DE, 0.25, na.rm = TRUE),
        p50 = median(DE, na.rm = TRUE),
        p75 = quantile(DE, 0.75, na.rm = TRUE),
        p90 = quantile(DE, 0.90, na.rm = TRUE),
        .groups = "drop"
      )
    p <- ggplot2::ggplot(defan, ggplot2::aes(x = Year, group = 1)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = p10, ymax = p90), fill = "#2E86AB", alpha = 0.15) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = p25, ymax = p75), fill = "#2E86AB", alpha = 0.35) +
      ggplot2::geom_line(ggplot2::aes(y = p50), colour = "#2E86AB", linewidth = 1.2) +
      ggplot2::geom_hline(yintercept = 0.79, linetype = "dashed", colour = "red", linewidth = 0.8) +
      ggplot2::labs(title = "D/E Ratio Fan Chart", x = NULL, y = "Debt/Equity") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  output$afnafn <- plotly::renderPlotly({
    afndf <- as.data.frame(afnsim / 1000)
    colnames(afndf) <- paste0("Year ", 1:yearssim)
    afnlong <- afndf %>%
      dplyr::mutate(sim = dplyr::row_number()) %>%
      tidyr::pivot_longer(-sim, names_to = "Year", values_to = "AFNM") %>%
      dplyr::filter(AFNM >= 0)
    
    p <- ggplot2::ggplot(afnlong, ggplot2::aes(x = Year, y = AFNM)) +
      ggplot2::geom_boxplot(fill = "#E69F00", alpha = 0.7, outlier.size = 0.4) +
      ggplot2::labs(title = "AFN Distribution (EUR M)", x = NULL, y = "AFN EUR M") +
      ggplot2::theme_minimal()
    
    plotly::ggplotly(p)
  })
  
  output$afnnde <- plotly::renderPlotly({
    ndedf <- as.data.frame(ndebitdasim)
    colnames(ndedf) <- paste0("Year ", 1:yearssim)
    ndefan <- ndedf %>%
      dplyr::mutate(sim = dplyr::row_number()) %>%
      tidyr::pivot_longer(-sim, names_to = "Year", values_to = "NDE") %>%
      dplyr::group_by(Year) %>%
      dplyr::summarise(
        p10 = quantile(NDE, 0.10, na.rm = TRUE),
        p25 = quantile(NDE, 0.25, na.rm = TRUE),
        p50 = median(NDE, na.rm = TRUE),
        p75 = quantile(NDE, 0.75, na.rm = TRUE),
        p90 = quantile(NDE, 0.90, na.rm = TRUE),
        .groups = "drop"
      )
    p <- ggplot2::ggplot(ndefan, ggplot2::aes(x = Year, group = 1)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = p10, ymax = p90), fill = "#28B463", alpha = 0.15) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = p25, ymax = p75), fill = "#28B463", alpha = 0.35) +
      ggplot2::geom_line(ggplot2::aes(y = p50), colour = "#28B463", linewidth = 1.2) +
      ggplot2::geom_hline(yintercept = 2.5, linetype = "dashed", colour = "red", linewidth = 0.8) +
      ggplot2::labs(title = "Net Debt / EBITDA Fan Chart", x = NULL, y = "Net Debt / EBITDA") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  interdata <- shiny::eventReactive(input$calcinter, {
    w <- input$waccinter / 100
    g <- input$ginter / 100
    basefcfflocal <- c(284000, 312000, 335000, 358000, 374000)
    yrs <- 1:5
    disc <- 1 / (1 + w)^yrs
    pvfcf <- sum(basefcfflocal * disc)
    terminal <- basefcfflocal[5] * (1 + g) / (w - g) / (1 + w)^5
    ev <- pvfcf + terminal
    netdebt <- 2140000
    equityval <- ev - netdebt
    sharestotal <- 306300000
    price <- equityval * 1000 / sharestotal
    list(pvfcf = pvfcf, terminal = terminal, ev = ev, equityval = equityval, price = price)
  }, ignoreNULL = FALSE)
  
  output$intertable <- shiny::renderTable({
    res <- interdata()
    data.frame(
      Component = c("PV of 5Y FCFFs", "Terminal Value Gordon", "Enterprise Value", "Less Net Debt", "Equity Value"),
      `EUR millions` = c(
        round(res$pvfcf / 1000, 1),
        round(res$terminal / 1000, 1),
        round(res$ev / 1000, 1),
        round(-2140000 / 1000, 1),
        round(res$equityval / 1000, 1)
      ),
      check.names = FALSE
    )
  }, striped = TRUE, hover = TRUE)
  
  output$interprice <- shiny::renderPrint({
    res <- interdata()
    cat(sprintf("Fair Value per Share: %.2f EUR", res$price))
  })
  
  output$gbmmeanbox <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(paste0("€", round(mean(gbmfinal), 2)), "GBM Mean Price 1Y", icon = shiny::icon("euro-sign"), color = "blue")
  })
  
  output$gbmlobox <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(paste0("€", round(stats::quantile(gbmfinal, 0.025), 2)), "2.5th Percentile", icon = shiny::icon("arrow-down"), color = "red")
  })
  
  output$gbmhibox <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(paste0("€", round(stats::quantile(gbmfinal, 0.975), 2)), "97.5th Percentile", icon = shiny::icon("arrow-up"), color = "green")
  })
  
  output$gbmchart <- plotly::renderPlotly({
    nshow <- min(100, nrow(gbmpaths))
    pathdf <- as.data.frame(t(gbmpaths[1:nshow, ])) %>%
      dplyr::mutate(Step = 0:(nrow(.) - 1)) %>%
      tidyr::pivot_longer(-Step, names_to = "Sim", values_to = "Price")
    p <- ggplot2::ggplot(pathdf, ggplot2::aes(x = Step, y = Price, group = Sim)) +
      ggplot2::geom_line(alpha = 0.12, colour = "#28B463", linewidth = 0.35) +
      ggplot2::geom_hline(yintercept = mktprice, colour = "black", linetype = "dashed", linewidth = 0.9) +
      ggplot2::labs(title = sprintf("GBM Paths: S0 %.2f, μ %.1f%%, σ %.1f%% | Mean 1Y EUR %.2f",
                                    S0gbm, mugbm * 100, sigmagbm * 100, mean(gbmfinal)),
                    x = "Trading days", y = "Price EUR") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
  
  output$dlaltman2 <- shiny::downloadHandler(
    filename = function() paste0("altman_", Sys.Date(), ".csv"),
    content = function(file) write.csv(altmanclean, file, row.names = FALSE)
  )
  
  output$dlcapm <- shiny::downloadHandler(
    filename = function() paste0("capm_", Sys.Date(), ".csv"),
    content = function(file) write.csv(diag_tbl, file, row.names = FALSE)
  )
  
  output$dlratios <- shiny::downloadHandler(
    filename = function() paste0("ratios_", Sys.Date(), ".csv"),
    content = function(file) write.csv(financialratios, file, row.names = FALSE)
  )
  
  output$dlratios2 <- shiny::downloadHandler(
    filename = function() paste0("ratios_", Sys.Date(), ".csv"),
    content = function(file) write.csv(financialratios, file, row.names = FALSE)
  )
  
  output$dldcf <- shiny::downloadHandler(
    filename = function() paste0("dcf_", Sys.Date(), ".csv"),
    content = function(file) write.csv(comparisontbl, file, row.names = FALSE)
  )
  
  output$dldcf2 <- shiny::downloadHandler(
    filename = function() paste0("dcf_", Sys.Date(), ".csv"),
    content = function(file) write.csv(comparisontbl, file, row.names = FALSE)
  )
  
  output$dlwacc <- shiny::downloadHandler(
    filename = function() paste0("wacc_", Sys.Date(), ".csv"),
    content = function(file) write.csv(wacc_summary, file, row.names = FALSE)
  )
  
  output$dlcrack <- shiny::downloadHandler(
    filename = function() paste0("crack_", Sys.Date(), ".csv"),
    content = function(file) write.csv(crackdf, file, row.names = FALSE)
  )
  
  output$dlcrack2 <- shiny::downloadHandler(
    filename = function() paste0("crack_", Sys.Date(), ".csv"),
    content = function(file) write.csv(crackdf, file, row.names = FALSE)
  )

  # ── REACTIVE WACC COMPUTATION ──────────────────────────────
  rwacc <- shiny::reactive({
    rf_r  <- input$rf_in  / 100
    erp_r <- input$erp_in / 100
    b_r   <- input$beta_in
    kd_r  <- input$kd_in  / 100
    we_r  <- input$we_in  / 100
    t_r   <- input$tax_in / 100
    wd_r  <- 1 - we_r
    ke_r  <- rf_r + b_r * erp_r
    kd_at <- kd_r * (1 - t_r)
    w_r   <- we_r * ke_r + wd_r * kd_at
    list(rf = rf_r, erp = erp_r, beta = b_r, kd_pre = kd_r,
         we = we_r, wd = wd_r, tax = t_r, ke = ke_r, kd = kd_at, wacc = w_r)
  })

  output$wacc_live_tbl <- DT::renderDT({
    p <- rwacc()
    tbl <- tibble::tibble(
      Parameter = c("Risk-Free Rate (Rf)", "Equity Risk Premium (ERP)",
                    "Levered Beta", "Pre-Tax Cost of Debt",
                    "Tax Rate", "Weight of Equity (We)", "Weight of Debt (Wd)",
                    "Cost of Equity (Ke) = Rf + β×ERP",
                    "After-Tax Cost of Debt = Kd×(1−t)",
                    "WACC = We×Ke + Wd×Kd(at)"),
      Value = paste0(round(c(p$rf, p$erp, p$beta, p$kd_pre,
                              p$tax, p$we, p$wd, p$ke, p$kd, p$wacc) * 100, 3), "%")
    )
    DT::datatable(tbl, options = list(dom = "t", pageLength = 12),
                  rownames = FALSE) |>
      DT::formatStyle("Parameter",
                      target = "row",
                      backgroundColor = DT::styleEqual(
                        "WACC = We×Ke + Wd×Kd(at)", "#d4edda"))
  })

  output$wacc_decomp_chart <- plotly::renderPlotly({
    p <- rwacc()
    df <- data.frame(
      Component = c("Ke contribution (We×Ke)",
                    "Kd contribution (Wd×Kd_at)",
                    "WACC"),
      Value     = round(c(p$we * p$ke, p$wd * p$kd, p$wacc) * 100, 3),
      Colour    = c("#2E86AB", "#E69F00", "#28B463")
    )
    plotly::plot_ly(df, x = ~Component, y = ~Value, type = "bar",
                    marker = list(color = df$Colour),
                    text = ~paste0(Value, "%"), textposition = "outside") |>
      plotly::layout(title = sprintf("WACC decomposition — Live: %.2f%%", p$wacc * 100),
                     yaxis = list(title = "%", range = c(0, max(df$Value) * 1.25)),
                     xaxis = list(title = ""))
  })

  output$wacc_sens_live <- DT::renderDT({
    w <- rwacc()$wacc
    g_vals   <- c(0.005, 0.010, 0.015, 0.018, 0.020, 0.025)
    w_vals   <- w + c(-0.015, -0.01, -0.005, 0, 0.005, 0.010)
    mat      <- outer(w_vals, g_vals, function(ww, gg) round((ww - gg) * 100, 2))
    rownames(mat) <- paste0(round(w_vals * 100, 2), "%")
    colnames(mat) <- paste0("g=", round(g_vals * 100, 1), "%")
    DT::datatable(as.data.frame(mat),
                  caption = "WACC–g spread (pp) — rows: WACC ±, cols: g",
                  options = list(dom = "t", scrollX = TRUE)) |>
      DT::formatStyle(colnames(mat),
                      background = DT::styleColorBar(c(0, 8), "#2E86AB"),
                      backgroundSize = "90% 80%",
                      backgroundRepeat = "no-repeat",
                      backgroundPosition = "center")
  })

  # ── STATIONARITY ──────────────────────────────────────────
  output$stat_tbl <- DT::renderDT({
    DT::datatable(stationarity_tbl, options = list(dom = "t"), rownames = FALSE) |>
      DT::formatStyle(c("ADF verdict", "KPSS verdict"),
                      color = DT::styleEqual(
                        c("Stationary ***", "Stationary **", "Stationary *",
                          "Stationary", "Non-stationary", "Non-stationary **", "Borderline *"),
                        c("green","green","orange","green","red","red","orange")),
                      fontWeight = "bold")
  })

  # ── GARCH ─────────────────────────────────────────────────
  output$garch_param_tbl <- DT::renderDT({
    DT::datatable(garch_param_tbl, options = list(dom = "t"), rownames = FALSE)
  })

  output$garch_vol_chart <- plotly::renderPlotly({
    p <- ggplot2::ggplot(garch_vol_df, ggplot2::aes(x = Date)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = pmin(cond_vol, static),
                                         ymax = pmax(cond_vol, static)),
                            fill = "#2E86AB", alpha = 0.15) +
      ggplot2::geom_line(ggplot2::aes(y = cond_vol, colour = "GARCH(1,1) σ"),
                          linewidth = 0.8) +
      ggplot2::geom_hline(ggplot2::aes(yintercept = sigmagbm, colour = "Static GBM σ"),
                           linetype = "dashed", linewidth = 1) +
      ggplot2::scale_colour_manual(
        values = c("GARCH(1,1) σ" = "#2E86AB", "Static GBM σ" = "#E74C3C")) +
      ggplot2::labs(title = sprintf(
        "Conditional vol: GARCH persistence=%.3f | Last σ_ann=%.1f%% | Static=%.1f%%",
        garch_persistence, garch_sigma_last * 100, sigmagbm * 100),
        x = NULL, y = "Annualised volatility", colour = NULL) +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })

  # ── VAR / GRANGER / IRF ───────────────────────────────────
  output$granger_tbl <- DT::renderDT({
    DT::datatable(granger_tbl, options = list(dom = "t"), rownames = FALSE) |>
      DT::formatStyle("Conclusion",
                      color = DT::styleEqual(
                        c("Granger-causes", "No Granger causality"),
                        c("#28B463", "#E74C3C")),
                      fontWeight = "bold")
  })

  output$irf_chart <- plotly::renderPlotly({
    p <- ggplot2::ggplot(irf_df, ggplot2::aes(x = Horizon)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = Lower, ymax = Upper),
                            fill = "#2E86AB", alpha = 0.25) +
      ggplot2::geom_line(ggplot2::aes(y = IRF), colour = "#2E86AB", linewidth = 1.2) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                           colour = "#E74C3C", linewidth = 0.8) +
      ggplot2::labs(
        title = "IRF: 1 s.d. oil shock → ELPE log-return (95% bootstrap CI)",
        x = "Horizon (days)", y = "Response") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })

  # ── COINTEGRATION ─────────────────────────────────────────
  output$jo_tbl <- DT::renderDT({
    DT::datatable(jo_tbl, options = list(dom = "t"), rownames = FALSE) |>
      DT::formatStyle("Conclusion",
                      color = DT::styleEqual(
                        c("Reject H₀ — cointegrated", "Fail to reject H₀"),
                        c("#28B463", "#E74C3C")),
                      fontWeight = "bold")
  })

  output$coint_levels_chart <- plotly::renderPlotly({
    df_lvl <- elpe_brent_lvl |>
      dplyr::mutate(
        elpe_norm  = elpe_close / elpe_close[1],
        brent_norm = brent      / brent[1]
      )
    p <- ggplot2::ggplot(df_lvl, ggplot2::aes(x = as.Date(Date))) +
      ggplot2::geom_line(ggplot2::aes(y = elpe_norm,  colour = "ELPE (normalised)"),  linewidth = 0.9) +
      ggplot2::geom_line(ggplot2::aes(y = brent_norm, colour = "Brent (normalised)"), linewidth = 0.9) +
      ggplot2::scale_colour_manual(values = c("ELPE (normalised)" = "#2E86AB",
                                               "Brent (normalised)" = "#E69F00")) +
      ggplot2::labs(title = "ELPE vs Brent — Normalised Price Levels (base=1)",
                    x = NULL, y = "Normalised price", colour = NULL) +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })

  # ── STRUCTURAL BREAK ──────────────────────────────────────
  output$break_tbl <- DT::renderDT({
    DT::datatable(break_tbl, options = list(dom = "t"), rownames = FALSE)
  })

  output$cusum_chart <- plotly::renderPlotly({
    boundary_val <- 1.36 / sqrt(nrow(cusum_df))
    p <- ggplot2::ggplot(cusum_df, ggplot2::aes(x = Date, y = CUSUM)) +
      ggplot2::geom_line(colour = "#2E86AB", linewidth = 0.9) +
      ggplot2::geom_hline(yintercept =  boundary_val, linetype = "dashed",
                           colour = "red", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = -boundary_val, linetype = "dashed",
                           colour = "red", linewidth = 0.8) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
      ggplot2::labs(
        title = "OLS-CUSUM — red dashed = 5% significance boundary",
        x = NULL, y = "Empirical fluctuation process") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })

  # ── PERFORMANCE METRICS ───────────────────────────────────
  output$perf_tbl <- DT::renderDT({
    DT::datatable(perf_tbl, options = list(dom = "t"), rownames = FALSE)
  })

  output$perf_chart <- plotly::renderPlotly({
    df_perf <- tibble::tibble(
      Metric  = rep(c("Sharpe Ratio", "Treynor Ratio"), 2),
      Company = c("ELPE", "ELPE", "MOH", "MOH"),
      Value   = c(sharpe_elpe, treynor_elpe, sharpe_moh, treynor_moh)
    )
    p <- ggplot2::ggplot(df_perf, ggplot2::aes(x = Metric, y = Value, fill = Company)) +
      ggplot2::geom_col(position = "dodge", width = 0.55) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey40") +
      ggplot2::scale_fill_manual(values = c("ELPE" = "#2E86AB", "MOH" = "#E69F00")) +
      ggplot2::labs(title = "Risk-Adjusted Return Metrics: ELPE vs MOH",
                    x = NULL, y = "Ratio value") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })

  # ── VaR / CVaR ────────────────────────────────────────────
  output$var95_dcf_box <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(
      paste0("€", round(var95_dcf, 2)), "VaR 95% — Markov DCF",
      icon = shiny::icon("exclamation-triangle"), color = "red")
  })
  output$cvar95_dcf_box <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(
      paste0("€", round(cvar95_dcf, 2)), "CVaR 95% — Markov DCF",
      icon = shiny::icon("fire"), color = "orange")
  })
  output$var95_gbm_box <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(
      paste0("€", round(var95_gbm, 2)), "VaR 95% — GBM 1Y",
      icon = shiny::icon("exclamation-triangle"), color = "yellow")
  })
  output$cvar95_gbm_box <- shinydashboard::renderValueBox({
    shinydashboard::valueBox(
      paste0("€", round(cvar95_gbm, 2)), "CVaR 95% — GBM 1Y",
      icon = shiny::icon("fire"), color = "purple")
  })

  output$risk_tbl <- DT::renderDT({
    DT::datatable(risk_tbl, options = list(dom = "t"), rownames = FALSE)
  })

  output$risk_dist_chart <- plotly::renderPlotly({
    df_risk <- dplyr::bind_rows(
      tibble::tibble(Model = "Markov DCF", Price = pricepershare),
      tibble::tibble(Model = "GBM 1Y",    Price = gbmfinal)
    )
    p <- ggplot2::ggplot(df_risk, ggplot2::aes(x = Price, fill = Model)) +
      ggplot2::geom_histogram(bins = 60, alpha = 0.55, position = "identity",
                               colour = "white", linewidth = 0.2) +
      ggplot2::geom_vline(xintercept = var95_dcf,  colour = "#E74C3C",
                           linewidth = 1.1, linetype = "solid") +
      ggplot2::geom_vline(xintercept = cvar95_dcf, colour = "#C0392B",
                           linewidth = 1.1, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = var95_gbm,  colour = "#E69F00",
                           linewidth = 1.1, linetype = "solid") +
      ggplot2::geom_vline(xintercept = cvar95_gbm, colour = "#D35400",
                           linewidth = 1.1, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = currentprice, colour = "black",
                           linewidth = 1.2, linetype = "dotdash") +
      ggplot2::annotate("text", x = var95_dcf,  y = Inf, vjust = 2,
                         label = "VaR95 DCF",  colour = "#E74C3C", size = 3) +
      ggplot2::annotate("text", x = var95_gbm,  y = Inf, vjust = 4,
                         label = "VaR95 GBM",  colour = "#E69F00", size = 3) +
      ggplot2::scale_fill_manual(values = c("Markov DCF" = "#2E86AB",
                                             "GBM 1Y" = "#28B463")) +
      ggplot2::labs(title = "Price distributions with VaR & CVaR (solid = VaR, dashed = CVaR)",
                    x = "EUR/share", y = "Count") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })

  # ── EVA & REVERSE DCF ─────────────────────────────────────
  output$eva_tbl <- DT::renderDT({
    DT::datatable(eva_tbl, options = list(dom = "t", scrollX = TRUE), rownames = FALSE) |>
      DT::formatStyle("EVA — BU ★ (EUR k)",
                      color = DT::styleInterval(0, c("#E74C3C", "#28B463")),
                      fontWeight = "bold")
  })

  output$eva_chart <- plotly::renderPlotly({
    eva_long <- tidyr::pivot_longer(
      eva_tbl,
      cols = c("EVA — BU ★ (EUR k)", "EVA — OLS (EUR k)"),
      names_to = "Method", values_to = "EVA"
    )
    p <- ggplot2::ggplot(eva_long, ggplot2::aes(x = factor(Year), y = EVA, fill = Method)) +
      ggplot2::geom_col(position = "dodge", width = 0.65) +
      ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
      ggplot2::scale_fill_manual(values = c("EVA — BU ★ (EUR k)" = "#2E86AB",
                                             "EVA — OLS (EUR k)"   = "#E69F00")) +
      ggplot2::labs(title = "EVA (EUR k) — Bottom-up vs OLS WACC (2026–2030)",
                    x = NULL, y = "EVA (EUR k)", fill = NULL) +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })

  output$reverse_tbl <- DT::renderDT({
    DT::datatable(reverse_tbl, options = list(dom = "t"), rownames = FALSE)
  })

  output$reverse_chart <- plotly::renderPlotly({
    df_rev <- tibble::tibble(
      Scenario  = reverse_tbl$Scenario,
      Price     = c(currentprice, currentprice * 0.80, currentprice * 1.20),
      implied_g = c(imp_g_base, imp_g_bear, imp_g_bull) * 100
    )
    plotly::plot_ly(df_rev, x = ~Scenario, y = ~implied_g, type = "bar",
                    text = ~paste0(round(implied_g, 2), "%"),
                    textposition = "outside",
                    marker = list(color = c("#2E86AB", "#E74C3C", "#28B463"))) |>
      plotly::add_segments(x = 0.5, xend = 3.5,
                            y = growthterminal * 100, yend = growthterminal * 100,
                            line = list(dash = "dash", color = "black", width = 2),
                            name = "Base g = 1.8%") |>
      plotly::layout(title = "Implied Terminal Growth Rate by Price Scenario",
                     yaxis = list(title = "Implied g (%)"),
                     xaxis = list(title = ""),
                     showlegend = TRUE)
  })

  # ── WACC COMPARISON ───────────────────────────────────────
  output$wacc_comp_tbl <- DT::renderDT({
    DT::datatable(wacc_comparison_tbl,
                  options = list(dom = "t", pageLength = 10), rownames = FALSE) |>
      DT::formatStyle("Parameter", fontWeight = "bold") |>
      DT::formatStyle("Bottom-up (Base Case)",
                      backgroundColor = "#d4edda")
  })

  output$wacc_comp_chart <- plotly::renderPlotly({
    df_wc <- data.frame(
      Method     = c("OLS (robustness)", "Bottom-up ★ BASE"),
      WACC       = round(c(wacc_ols, wacc_bu) * 100, 2),
      Beta       = round(c(beta_relevered, beta_bu), 4),
      Ke         = round(c(ke_ols, ke_bu) * 100, 2)
    )
    plotly::plot_ly(df_wc, x = ~Method,
                    y = ~WACC, type = "bar", name = "WACC (%)",
                    marker = list(color = c("#E69F00", "#2E86AB")),
                    text = ~paste0(WACC, "%"), textposition = "outside") |>
      plotly::add_trace(y = ~Ke, type = "bar", name = "Ke (%)",
                        marker = list(color = c("#F0C27F", "#85C1E9")),
                        text = ~paste0(Ke, "%"), textposition = "outside") |>
      plotly::layout(
        title = sprintf("WACC & Ke: OLS=%.2f%% vs Bottom-up=%.2f%%",
                        wacc_ols * 100, wacc_bu * 100),
        barmode = "group",
        yaxis = list(title = "% p.a.", range = c(0, 14)),
        xaxis = list(title = ""),
        legend = list(orientation = "h", x = 0.3, y = -0.15))
  })

  output$seg_beta_tbl <- DT::renderDT({
    tbl <- segment_betas |>
      dplyr::mutate(
        weight  = paste0(round(weight  * 100, 1), "%"),
        contrib = round(contrib, 4),
        revenue = format(revenue, big.mark = ",")
      ) |>
      dplyr::select(Segment, beta_u, revenue, weight, contrib, source)
    names(tbl) <- c("Segment", "βU (Damodaran)", "Revenue (EUR k)",
                    "Weight", "β Contribution", "Source")
    DT::datatable(tbl, options = list(dom = "t"), rownames = FALSE)
  })

  output$beta_buildup_tbl <- DT::renderDT({
    tbl <- tibble::tibble(
      Step  = c("OLS monthly β (5Y)", "Blume adjustment",
                "Hamada unlever (D/E=1.099)", "Hamada re-lever (D/E=0.79)",
                "──────────────────",
                "Composite βU (rev-wtd)", "Re-levered bottom-up β (D/E=1.099)",
                "── WACC inputs ──",
                "OLS β → Ke", "OLS β → WACC",
                "Bottom-up β → Ke ★", "Bottom-up β → WACC ★"),
      Value = c(sprintf("%.4f", beta_ols),
                sprintf("%.4f", beta_blume),
                sprintf("%.4f", beta_unlevered),
                sprintf("%.4f", beta_relevered),
                "──",
                sprintf("%.4f", beta_u_bu),
                sprintf("%.4f", beta_bu),
                "──",
                sprintf("%.2f%%", ke_ols  * 100),
                sprintf("%.2f%%", wacc_ols * 100),
                sprintf("%.2f%%", ke_bu   * 100),
                sprintf("%.2f%%", wacc_bu  * 100))
    )
    DT::datatable(tbl, options = list(dom = "t", pageLength = 15),
                  rownames = FALSE) |>
      DT::formatStyle("Step",
                      target = "row",
                      backgroundColor = DT::styleEqual(
                        c("Bottom-up β → WACC ★", "Re-levered bottom-up β (D/E=1.099)"),
                        c("#d4edda", "#d4edda")))
  })

  # ── CORRELATED MONTE CARLO ────────────────────────────────
  output$corr_compare_tbl <- DT::renderDT({
    DT::datatable(corr_compare_tbl, options = list(dom = "t"), rownames = FALSE)
  })

  output$corr_dist_chart <- plotly::renderPlotly({
    p <- ggplot2::ggplot(corr_dist_df, ggplot2::aes(x = Price, fill = Model)) +
      ggplot2::geom_density(alpha = 0.55, colour = "white", linewidth = 0.3) +
      ggplot2::geom_vline(xintercept = mean(corr_prices), colour = "#E69F00",
                           linewidth = 1.2, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = meanofmeans, colour = "#2E86AB",
                           linewidth = 1.2, linetype = "dashed") +
      ggplot2::geom_vline(xintercept = corr_var95, colour = "#E74C3C",
                           linewidth = 1.0, linetype = "solid") +
      ggplot2::geom_vline(xintercept = currentprice, colour = "black",
                           linewidth = 1.1, linetype = "dotdash") +
      ggplot2::scale_fill_manual(values = c("Correlated MC"   = "#E69F00",
                                             "Base Markov DCF" = "#2E86AB")) +
      ggplot2::annotate("text", x = corr_var95 - 0.3, y = Inf, vjust = 2,
                         label = "Corr VaR95", colour = "#E74C3C", size = 3, angle = 90) +
      ggplot2::labs(
        title = "Density: Correlated MC vs Base Markov DCF (dashed = mean, solid red = VaR95)",
        x = "EUR/share", y = "Density") +
      ggplot2::theme_minimal()
    plotly::ggplotly(p)
  })
}

message("Clean dashboard section loaded.")
message("Launching Shiny dashboard...")
shiny::shinyApp(ui = ui, server = server)
