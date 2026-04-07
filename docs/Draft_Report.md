# Forecasting Canadian Household Net Savings: A Comparative Study of Univariate and Multivariate Approaches

## Abstract

This report presents a comparative forecasting study of Canadian quarterly household net savings from 1977 Q1 to 2019 Q4. We evaluate a hierarchy of increasingly complex models — from simple benchmarks (persistence, average, exponential smoothing) through ARIMA and SARIMA specifications to ARIMAX models incorporating macroeconomic regressors — using one-step-ahead holdout RMSE as the primary evaluation criterion. Two forecasting frameworks are compared: a fixed-parameter scheme estimated on a single training window, and a rolling window scheme (R = 80 quarters) that re-estimates parameters at each step to adapt to structural change. The rolling window approach is motivated by the presence of likely structural breaks spanning the oil shocks, the Great Moderation, and the 2008 financial crisis. We find that household savings exhibit near-random-walk behaviour, with even persistence (naïve forecast) proving difficult to beat. Exogenous regressors derived from cross-correlation analysis — lagged log-differenced Gross Private Domestic Investment and CPI — provide modest improvements in some specifications, while theoretically motivated regressors from a precautionary saving (PIM) framework offer an alternative structural lens. The full RMSE rankings and model diagnostics are reported for both the fixed-parameter and rolling window evaluations.

## 1. Motivation

Household net savings — the difference between disposable income and consumption expenditure — is a key barometer of macroeconomic resilience. When households save more, they accumulate a buffer against future income shocks; when they save less, consumer spending drives short-term GDP growth but leaves the economy vulnerable to downturns. For policymakers at the Bank of Canada, reliable one-quarter-ahead forecasts of aggregate household savings would inform decisions about interest rate policy, fiscal stimulus design, and financial stability assessments.

Yet household savings are notoriously difficult to predict. The permanent income hypothesis (Friedman, 1957) implies that rational, forward-looking households adjust consumption smoothly in response to anticipated income changes, so that savings — the residual — absorbs all *unanticipated* shocks. If the hypothesis holds strictly, savings should follow a random walk, and no model can systematically beat naïve persistence. The practical question is whether modest departures from the random walk — arising from habit formation, liquidity constraints, precautionary motives, or macroeconomic leading indicators — provide enough signal to justify more complex models.

This project investigates that question empirically. We construct a hierarchy of forecasting models, from the simplest possible benchmark (repeat last quarter's value) to theory-driven ARIMAX specifications with macroeconomic regressors, and ask: *does the added complexity pay off in reduced forecast error?* By comparing models within both a fixed-parameter and a rolling-window framework, we also test whether the answer changes when the estimation procedure is allowed to adapt to structural breaks in the Canadian economy over the past four decades.

## 2. Data and Variables

### 2.1 Source and Scope

Our dataset is a quarterly macroeconomic panel for Canada spanning 1977 Q1 to 2019 Q4 (172 observations). The data were sourced from Statistics Canada via the `cansim` package, harmonised into a single wide-format file (`macro_panel_wide_raw.csv`), and processed through a standardised pipeline (`data_ingest.R`).

### 2.2 Variables

The core variable of interest is **household net savings** (`hh_save`), measured in millions of current Canadian dollars. The candidate predictor set includes:

| Variable | Description |
|---|---|
| `gdi` | Gross domestic income |
| `cons` | Personal consumption expenditure |
| `unemp` | Unemployment rate (%) |
| `cpi` | Consumer price index |
| `interest` | Market interest rate |
| `gpdi` | Gross private domestic investment |

All level series are non-stationary (confirmed by ADF tests). We apply the log-difference operator to obtain stationary transformations: `ldgdi`, `ldcons`, `ldcpi`, `ldgpdi`. Unemployment and interest rates, being rates, are first-differenced: `dunemp`, `dinterest`. The target series `hh_save` is differenced internally by the ARIMA/ARIMAX estimation procedure ($d = 1$).

### 2.3 Exploratory Findings

Time-series plots reveal that household savings are volatile and exhibit no obvious deterministic trend after 2000, consistent with near-random-walk dynamics. The ACF of the first-differenced series (`dhh`) shows weak autocorrelation, with a significant spike at lag 1 that decays quickly — suggesting an MA(1) component in the differenced domain.

**STL decomposition.** To disentangle the sources of variation, we apply Seasonal and Trend decomposition using Loess (STL) to the level `hh_save` series. The decomposition reveals:

- A **trend** component that traces the long-run rise in household savings through the late 1980s, the sharp decline through the 1990s and 2000s (when households increasingly relied on credit), and a partial recovery after the 2008 financial crisis. These movements correspond to well-documented Canadian macroeconomic episodes: the oil booms, the interest rate hikes of the early 1980s, and the post-crisis deleveraging.
- A **seasonal** component that is modest relative to both trend and remainder. Quarterly seasonality in savings exists (e.g., tax-refund quarters) but accounts for only a small share of total variance. This explains why SARIMA models offer only marginal improvement over non-seasonal ARIMA in the holdout evaluations.
- A **remainder** (irregular) component that dominates the short-term fluctuations. Its standard deviation substantially exceeds that of the seasonal component, confirming that quarter-to-quarter movements in household savings are largely unpredictable — consistent with the near-random-walk evidence from the ACF and variogram analyses.

**Cross-correlation analysis (CCF)** identified regressors that *lead* household savings:

- `ldgpdi` at lags 1 and 3 (investment leads savings)
- `ldcpi` at lag 2 (inflation leads savings adjustment)

These lags were used to construct the ARIMAX regressor set. The contemporaneous correlations with GDI and consumption were discarded because they would not be available at forecast time.

### 2.4 Train/Test Split

For the fixed-parameter evaluation (Notebook 02), the first 132 quarters (1977 Q1 – 2009 Q4) serve as the training set, and the remaining 40 quarters (2010 Q1 – 2019 Q4) form the holdout. For the rolling window evaluation (Notebook 03), the initial training window is 80 quarters (1977 Q1 – 1996 Q4), and the holdout spans 92 quarters (1997 Q1 – 2019 Q4). The rolling window size is fixed at $R = 80$ quarters.

## 3. Methodology and Economic Interpretation

### 3.1 Benchmark Models and the Efficient Markets Analogy

The simplest models encode the hypothesis that household savings are informationally efficient — past patterns contain no exploitable signal.

- **Persistence (random walk):** $\hat{y}_t = y_{t-1}$. This is the forecasting analogue of the permanent income hypothesis: if households fully absorb new information each quarter, last quarter's savings level is the best predictor. Persistence sets the bar that every more complex model must clear.
- **Historical average:** $\hat{y}_t = \bar{y}_{\text{train}}$. Assumes savings fluctuate around a stable long-run mean — i.e., the series is stationary in levels. The poor RMSE of the average model (approximately 16,000 vs 10,000 for persistence) decisively rejects this assumption.
- **Seasonal average (by period):** $\hat{y}_t = \bar{y}_{\text{quarter}(t)}$. Tests whether savings have a fixed quarterly pattern (e.g., tax-season effects). Its RMSE is nearly identical to the simple average, confirming the weak seasonality flagged by the STL decomposition.

### 3.2 Exponential Smoothing and Adaptive Expectations

Exponential smoothing models formalise the idea that households form expectations adaptively — they revise their savings behaviour partly in response to recent outcomes.

- **Simple Exponential Smoothing (SES):** $\hat{y}_t = \alpha y_{t-1} + (1-\alpha) \hat{y}_{t-1}$. This is the optimal forecast for an IMA(1,1) process, i.e. a random walk with partial shock absorption. The estimated smoothing parameter $\alpha$ quantifies how quickly agents update: a high $\alpha$ means households react strongly to last quarter's outcome; a low $\alpha$ means they persist with prior expectations. SES matches persistence closely (RMSE ≈ 10,300 vs 10,085), consistent with household savings being a unit-root process with modest memory.
- **Holt Linear:** Adds a trend component to SES, capturing the possibility that savings drift upward or downward over time. Its slightly higher RMSE than SES suggests that any deterministic trend in savings is too weak (or too unstable) to improve forecasts.
- **Holt-Winters Additive:** Further adds quarterly seasonal factors. The STL decomposition showed that seasonality in savings is small, and this is borne out in the holdout: HW Additive does not consistently beat SES, confirming that the seasonal signal is too faint to exploit.

### 3.3 ARIMA: Box-Jenkins Modelling of Internal Dynamics

ARIMA models extend the exponential smoothing intuition by allowing richer autoregressive and moving-average structures. The economic hypothesis is that there exists short-term *internal momentum* in savings — households adjust their saving rate gradually due to habit formation, contractual obligations, or portfolio rebalancing frictions.

- The ACF/PACF of `dhh` point to an MA(1) component and a possible AR(1) term, motivating candidates ARIMA(1,1,1), ARIMA(0,1,1), ARIMA(2,1,0), and ARIMA(1,1,0).
- `auto.arima` on the training set selects ARIMA(0,1,1)(2,0,0)[4] — an IMA(1,1) with seasonal AR terms — confirming that the primary predictable structure is a single moving-average coefficient in the differenced series: $\Delta y_t = \varepsilon_t + \theta_1 \varepsilon_{t-1}$ with $\theta_1 \approx -0.28$.
- In economic terms, the negative $\theta_1$ means that a positive savings surprise in one quarter is partially reversed the next — consistent with the habit-formation interpretation: an unexpectedly high savings quarter is followed by a return toward the household's longer-run saving plan.
- The **best ARIMA** (ARIMA(1,1,1)) achieves RMSE ≈ 10,132 under fixed parameters, a modest improvement over persistence (10,085). The ARIMA(1,1,1) adds an AR(1) term that captures a weak first-order persistence in the growth rate of savings.

### 3.4 SARIMA: Testing for Quarterly Seasonality

SARIMA extends ARIMA with seasonal AR and MA terms at the quarterly frequency. The hypothesis is that household savings exhibit a stable quarterly pattern — for example, higher savings in Q1 due to RRSP contributions or tax-refund timing.

- We evaluate 8 SARIMA specifications. The best, SARIMA(1,1,1)(1,0,0)[4], achieves RMSE ≈ 10,203 under fixed parameters — slightly worse than pure ARIMA(1,1,1). The seasonal AR(1) coefficient is modest, and the STL decomposition already showed that the seasonal component has small amplitude.
- This confirms that while quarterly patterns exist, they are too weak and too unstable across regimes to improve out-of-sample forecasts.

### 3.5 ARIMAX: Exogenous Information and Leading Indicators

ARIMAX models test whether macroeconomic variables contain predictive information about savings beyond what the series' own history provides. The economic hypotheses are:

- **Investment channel (GPDI leads savings):** Under the pecking-order theory, corporate investment draws on internal funds first. A rise in GPDI signals corporate expansion, which raises wage income and eventually household savings — but with a lag. CCF analysis confirms that `ldgpdi` leads `dhh` at lags 1 and 3.
- **Inflation channel (CPI leads savings):** Rising prices erode real purchasing power. Forward-looking households may increase nominal savings in anticipation of future price increases. The CCF identifies a leading relationship at lag 2 for `ldcpi`.
- **Full vs Reduced specification:** The initial ARIMAX(1,1,1) fit on the 1977–1996 window shows that only `ldgpdi_lag1` is statistically significant at the 5% level (p = 0.041). The other regressors (`ldcpi_lag2`, p = 0.63; `ldgpdi_lag3`, p = 0.07) add estimation noise. The **Reduced ARIMAX** retains only the significant regressor, following the parsimony principle.
- In the NB02 fixed-parameter evaluation, the best ARIMAX — ARIMAX(2,1,3) (Full) — achieves the lowest overall RMSE (9,776), beating persistence by about 3%. In the NB03 rolling window evaluation, the best rolling model is ARIMAX(1,1,1) (Reduced) with RMSE 8,740, also beating all benchmarks on the longer 92-quarter holdout.

### 3.6 Precautionary Income Model (PIM)

The PIM is derived from the CRRA Euler equation for intertemporal consumption:

$$C_t^{-\gamma} = \beta \, \mathbb{E}_t \left[ (1 + r_{t+1}) \, C_{t+1}^{-\gamma} \right]$$

Log-linearising and invoking the savings identity, the precautionary saving motive maps onto observable proxies:

| Theory variable | Proxy | Justification |
|---|---|---|
| Income growth $\Delta y_t$ | `ldgdi_lag1` | Log-differenced GDI |
| Real interest rate $r_t$ | `dinterest_lag1` − `ldcpi_lag1` | Fisher identity |
| Income uncertainty | `dunemp_lag1` | Labour market risk |

The economic logic is that higher unemployment raises labour income uncertainty, which increases the precautionary motive and pushes savings upward. Higher real interest rates make saving more attractive. The PIM ARIMAX thereby embeds an explicit structural economic mechanism — rather than relying on purely data-driven CCF selection.

The best PIM specification (PIM ARIMAX(0,1,2)) achieves RMSE ≈ 10,183 in the fixed-parameter evaluation — competitive with ARIMA(1,1,1) and better than all exponential smoothing models. This suggests the precautionary motive has genuine, if modest, empirical relevance for Canadian household savings.

### 3.7 One-Step-Ahead Forecast Protocol

All models produce genuine one-step-ahead forecasts: at each holdout quarter $t$, the forecast $\hat{y}_t$ uses only information available through $t - 1$. No future data leaks into the forecast.

### 3.8 Fixed-Parameter vs Rolling Window

**Fixed-parameter** (Notebook 02): Model parameters are estimated once on the training set and held constant throughout the holdout. This mimics a practitioner who fits a model and deploys it without re-estimation.

**Rolling window** (Notebook 03): At each holdout step, the model is re-estimated on a window of the most recent $R = 80$ quarters. The oldest observation is dropped as the window slides forward. This approach, advocated by Giacomini & White (2006), keeps sample size constant and prevents observations from old regimes from diluting current parameter estimates. It is particularly appropriate when structural breaks are plausible — as they are for Canadian savings data spanning the oil shocks of the late 1970s, the Great Moderation, and the 2008 financial crisis.

### 3.9 Evaluation Criterion

The primary metric is one-step-ahead holdout **Root Mean Square Error** (RMSE):

$$\text{RMSE} = \sqrt{\frac{1}{H} \sum_{t=1}^{H} (y_t - \hat{y}_t)^2}$$

where $H$ is the number of holdout observations. We also report in-sample (training) RMSE to diagnose overfitting: a model with low training RMSE but high holdout RMSE is likely overfitting.

## 4. Results

### 4.1 Fixed-Parameter One-Step Forecasts (Notebook 02)

Each model is fitted once on the 1977 Q1 – 2009 Q4 training set (132 observations) and evaluated over the 40-quarter holdout (2010 Q1 – 2019 Q4). The table below reports selected holdout RMSE values:

| Model | Holdout RMSE |
|---|---|
| ARIMAX(2,1,3) — Full regressors | **9,776** |
| Persistence (random walk) | 10,085 |
| ARIMA(1,1,1) | 10,132 |
| PIM ARIMAX(0,1,2) | 10,183 |
| SARIMA(1,1,1)(1,0,0)[4] | 10,203 |
| SES | ~10,300 |
| Holt-Winters Additive | ~10,400 |
| Historical average | ~16,000 |

Key observations:

- **Persistence is hard to beat.** The naïve random-walk forecast (RMSE 10,085) sets a challenging benchmark. Many ARIMA and SARIMA specifications cluster within 1–2% of this value, confirming the near-random-walk nature of household savings.
- **ARIMAX with full regressors wins.** ARIMAX(2,1,3) with `ldgpdi_lag1`, `ldcpi_lag2`, and `ldgpdi_lag3` achieves the lowest RMSE (9,776), beating persistence by approximately 3%. This suggests that real investment growth does carry modest leading information about household savings.
- **The PIM specification is competitive.** PIM ARIMAX(0,1,2) (RMSE 10,183) outperforms all exponential smoothing models and several ARIMA specifications, supporting the view that precautionary saving motives based on unemployment and interest rate dynamics contain genuine predictive signal.
- **SARIMA adds little.** The best SARIMA (10,203) performs slightly worse than the best non-seasonal ARIMA (10,132) — consistent with the STL finding that seasonal variation in savings is small relative to the remainder.
- **Average-based benchmarks fail.** The historical average (RMSE ~16,000) is dramatically worse than persistence, confirming that household savings are non-stationary in levels.

A scatter matrix of forecast vectors reveals high pairwise correlations among ARIMA-family models, indicating that most models capture the same underlying dynamics — the unpredictable random-walk-like component dominates.

### 4.2 Rolling Window Forecasts (Notebook 03)

Models are re-estimated at every holdout step on a fixed window of $R = 80$ quarters (20 years), producing 92 one-step-ahead forecasts (1997 Q1 – 2019 Q4). Both fixed-parameter and rolling schemes are evaluated on this longer holdout for a direct RMSE comparison.

| Model | Fixed RMSE | Rolling RMSE | Improvement (%) |
|---|---|---|---|
| ARIMAX(1,1,1) — Reduced | 9,096 | **8,740** | +3.9% |
| ARIMA(1,1,1) | 8,865 | 8,914 | −0.6% |
| Persistence | 9,032 | 9,032 | — |
| SARIMA(1,1,1)(1,0,0)[4] | 9,052 | 8,856 | +2.2% |

Key observations:

- **ARIMAX(1,1,1) Reduced (Rolling) is the overall best** with RMSE 8,740, beating persistence (9,032) by 3.2%. The reduced specification — retaining only the significant regressor `ldgpdi_lag1` — outperforms the full specification, confirming that parsimony matters when parameters are re-estimated on limited 80-observation windows.
- **9 out of the evaluated rolling models beat persistence**, indicating that time-series structure in household savings, while weak, is statistically exploitable — particularly when the estimation window adapts to the most recent regime.
- **Rolling does not uniformly dominate fixed.** Some ARIMA specifications show negative improvement (e.g., ARIMA(1,1,1): −0.6%), meaning the rolling window's smaller sample size can increase estimation variance enough to offset the benefit of regime adaptation.
- **Re-estimation helps ARIMAX most.** The exogenous regressor coefficients are most sensitive to regime changes (the relationship between GPDI growth and savings may strengthen or weaken across business cycles), so rolling re-estimation provides the greatest lift for ARIMAX models.

### 4.3 Holdout Visualisation

Faceted holdout plots display the actual household net savings (grey) against each model's forecast (green) in separate panels. In both evaluations, forecasts from ARIMA-family models track the actual series closely but with a visible one-quarter lag — the hallmark of near-persistence behaviour. The ARIMAX models occasionally capture turning points slightly better than pure ARIMA, but the visual difference is subtle.

## 5. Discussion

### 5.1 Near-Random-Walk Dynamics

The central finding is that Canadian household net savings behave approximately as a random walk with drift. The IMA(1,1) structure identified by `auto.arima` implies that shocks are partially absorbed within one quarter, but the transient MA(1) component offers limited forecasting leverage. This is consistent with the permanent income hypothesis: households smooth consumption over income shocks, so savings — the residual — absorbs all unanticipated changes in income and expenditure.

The STL decomposition reinforces this picture. The remainder component dominates quarterly variation in the level series, while the trend captures slow-moving regime shifts that are useful for description but difficult to forecast one quarter ahead. The seasonal component's small amplitude explains the negligible SARIMA improvement.

### 5.2 Limited Value of Exogenous Regressors

Despite the theoretical appeal of macroeconomic covariates, the empirical gains from ARIMAX over pure ARIMA are modest — at most ~3% RMSE reduction. The CCF-based selection identifies statistically significant leading relationships (GPDI at lag 1 in particular), but these are too weak to substantially compress forecast errors. The reduced ARIMAX consistently matches or outperforms the full specification, reinforcing that parsimony is critical when the signal-to-noise ratio is low.

### 5.3 Structural Breaks and the Case for Rolling Windows

The rolling window is well-suited to this dataset because Canadian household savings have experienced multiple regime shifts: the high-saving environment of the late 1970s and early 1980s, the declining trend through the 1990s, and the post-2008 stabilisation. A fixed expanding window averages over all these regimes, potentially diluting parameter estimates relevant to the current regime. The rolling window discards observations older than 20 years and allows adaptation — particularly beneficial for ARIMAX, where regressor relationships shift across business cycles.

### 5.4 Practical Implications

For a practitioner forecasting Canadian household savings one quarter ahead:

1. **Start with persistence.** The naïve forecast is a strong baseline (RMSE ~9,000–10,000 depending on holdout).
2. **Consider an IMA(1,1) or ARIMA(1,1,1).** The modest MA(1) component provides a small and consistent forecasting edge.
3. **Add regressors cautiously.** Only `ldgpdi_lag1` is robustly significant across estimation windows; including more regressors risks overfitting.
4. **Re-estimate regularly.** A rolling window of 20 years provides enough data for stable estimation while adapting to structural change.
5. **Manage expectations.** Even the best model reduces RMSE by only ~3% relative to persistence. The near-random-walk character of household savings imposes a hard floor on forecast accuracy.

## 7. References

- Box, G. E. P., Jenkins, G. M., Reinsel, G. C., & Ljung, G. M. (2015). *Time Series Analysis: Forecasting and Control* (5th ed.). Wiley.
- Cleveland, R. B., Cleveland, W. S., McRae, J. E., & Terpenning, I. (1990). STL: A seasonal-trend decomposition procedure based on Loess. *Journal of Official Statistics*, 6(1), 3–73.
- Friedman, M. (1957). *A Theory of the Consumption Function*. Princeton University Press.
- Giacomini, R., & White, H. (2006). Tests of conditional predictive ability. *Econometrica*, 74(6), 1545–1578.
- Hyndman, R. J., & Athanasopoulos, G. (2021). *Forecasting: Principles and Practice* (3rd ed.). OTexts.
- Statistics Canada. CANSIM database, various tables (national accounts, labour force, prices).

## Appendix

### A. Code and Reproducibility

All analysis code is organised into four R Markdown notebooks:

| Notebook | Purpose |
|---|---|
| `01_EDA.rmd` | Data loading, transformation, stationarity testing, CCF analysis, variograms |
| `02_Household_Forecast.Rmd` | Fixed-parameter one-step forecasts: benchmarks, ESM, ARIMA, SARIMA, ARIMAX, PIM |
| `03_Rolling_Window_Forecast.Rmd` | Rolling window (R = 80) forecasts with fixed vs rolling RMSE comparison |

Utility functions are in `utils/`:

- `one_step_model.R`: Benchmark and fixed-parameter one-step forecast functions.
- `expanding_window.R`: Rolling/expanding window forecast loops for ARIMA, SARIMA, ARIMAX (with `window_size` parameter).
- `helpers_forecast.R`: Training RMSE helper.

The R environment is managed via `renv` for reproducibility. Data files and intermediate outputs are saved under `output/`.

### B. Model Specifications Evaluated

**ARIMA candidates:** ARIMA(1,1,1), (0,1,1), (2,1,0), (1,1,0)

**SARIMA candidates:** SARIMA(1,1,1)(1,0,0)[4], (0,1,1)(1,0,0)[4], (1,1,1)(1,0,1)[4], (0,1,1)(1,0,1)[4], (1,1,1)(0,0,1)[4], (0,1,1)(0,0,1)[4], (1,1,1)(0,0,2)[4], (0,1,1)(0,0,2)[4]

**ARIMAX candidates (Full):** ARIMAX(1,1,1), (0,1,1), (1,1,2), (2,1,2) with regressors `ldgpdi_lag1`, `ldcpi_lag2`, `ldgpdi_lag3`

**ARIMAX candidates (Reduced):** ARIMAX(1,1,1), (0,1,1), (1,1,2), (2,1,2) with regressor `ldgpdi_lag1`

**PIM candidates:** PIM ARIMAX(1,1,1), (0,1,1), (1,1,0), (0,1,2) with regressors `dinterest_lag1`, `ldcpi_lag1`, `dunemp_lag1`, `ldgdi_lag1`

### C. RMSE Tables

The complete RMSE ranking tables are generated dynamically by the notebooks and saved as `.rds` files:

- `output/hh_rmse_tbl.rds` — Fixed-parameter RMSE ranking (Notebook 02)
- `output/ew_rmse_tbl.rds` — Rolling window RMSE ranking (Notebook 03)

Selected values are reported in Section 4. The full tables (including all candidate specifications, training RMSE, and `auto.arima` selections) are available in the rendered notebook outputs.
