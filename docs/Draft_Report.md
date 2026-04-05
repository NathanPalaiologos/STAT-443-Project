# Forecasting Canadian Household Net Savings: A Comparative Study of Univariate and Multivariate Approaches

## Abstract

This report presents a comparative forecasting study of Canadian quarterly household net savings from 1977 Q1 to 2019 Q4. We evaluate a hierarchy of increasingly complex models — from simple benchmarks (persistence, average, exponential smoothing) through ARIMA and SARIMA specifications to ARIMAX models incorporating macroeconomic regressors — using one-step-ahead holdout RMSE as the primary evaluation criterion. Two forecasting frameworks are compared: a fixed-parameter scheme estimated on a single training window, and a rolling window scheme (R = 80 quarters) that re-estimates parameters at each step to adapt to structural change. The rolling window approach is motivated by the presence of likely structural breaks spanning the oil shocks, the Great Moderation, and the 2008 financial crisis. We find that household savings exhibit near-random-walk behaviour, with even persistence (naïve forecast) proving difficult to beat. Exogenous regressors derived from cross-correlation analysis — lagged log-differenced Gross Private Domestic Investment and CPI — provide modest improvements in some specifications, while theoretically motivated regressors from a precautionary saving (PIM) framework offer an alternative structural lens. The full RMSE rankings and model diagnostics are reported for both the fixed-parameter and rolling window evaluations.

## 1. Data and Variables

### 1.1 Source and Scope

Our dataset is a quarterly macroeconomic panel for Canada spanning 1977 Q1 to 2019 Q4 (172 observations). The data were sourced from Statistics Canada via the `cansim` package, harmonised into a single wide-format file (`macro_panel_wide_raw.csv`), and processed through a standardised pipeline (`data_ingest.R`).

### 1.2 Variables

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

### 1.3 Exploratory Findings

Time-series plots reveal that household savings are volatile and exhibit no obvious deterministic trend after 2000, consistent with near-random-walk dynamics. The ACF of the first-differenced series (`dhh`) shows weak autocorrelation, with a significant spike at lag 1 that decays quickly — suggesting an MA(1) component in the differenced domain.

Cross-correlation function (CCF) analysis identified regressors that *lead* household savings:

- `ldgpdi` at lags 1 and 3 (investment leads savings)
- `ldcpi` at lag 2 (inflation leads savings adjustment)

These lags were used to construct the ARIMAX regressor set. The contemporaneous correlations with GDI and consumption were discarded because they would not be available at forecast time.

### 1.4 Train/Test Split

For the fixed-parameter evaluation (Notebook 02), the first 132 quarters (1977 Q1 – 2009 Q4) serve as the training set, and the remaining 40 quarters (2010 Q1 – 2019 Q4) form the holdout. For the rolling window evaluation (Notebook 03), the initial training window is 80 quarters (1977 Q1 – 1996 Q4), and the holdout spans 92 quarters (1997 Q1 – 2019 Q4). The rolling window size is fixed at $R = 80$ quarters.

## 2. Methodology

### 2.1 Model Hierarchy

We evaluate models in order of increasing complexity:

1. **Simple benchmarks**: Persistence (naïve / random walk), historical average, average by year, average by period.
2. **Exponential smoothing (ESM)**: Simple exponential smoothing (SES), Holt linear trend, Holt-Winters additive seasonal.
3. **ARIMA**: Manual candidates ARIMA(1,1,1), (0,1,1), (2,1,0), (1,1,0) chosen from ACF/PACF inspection.
4. **SARIMA**: Extensions with quarterly seasonal components — e.g., SARIMA(1,1,1)(1,0,0)[4], SARIMA(0,1,1)(0,0,1)[4], among others.
5. **ARIMAX (Full regressors)**: ARIMA orders paired with exogenous regressors `ldgpdi_lag1`, `ldcpi_lag2`, `ldgpdi_lag3`.
6. **ARIMAX (Reduced)**: Same ARIMA orders with only `ldgpdi_lag1`, the sole regressor significant at the 5% level.
7. **PIM (Precautionary Income Model)**: A theory-driven ARIMAX specification grounded in the CRRA Euler equation, using `dinterest_lag1`, `ldcpi_lag1`, `dunemp_lag1`, `ldgdi_lag1` as regressors.

### 2.2 One-Step-Ahead Forecast Protocol

All models produce genuine one-step-ahead forecasts: at each holdout quarter $t$, the forecast $\hat{y}_t$ uses only information available through $t - 1$. No future data leaks into the forecast.

### 2.3 Fixed-Parameter vs Rolling Window

**Fixed-parameter** (Notebook 02): Model parameters are estimated once on the training set and held constant throughout the holdout. This mimics a practitioner who fits a model and deploys it without re-estimation.

**Rolling window** (Notebook 03): At each holdout step, the model is re-estimated on a window of the most recent $R = 80$ quarters. The oldest observation is dropped as the window slides forward. This approach, advocated by Giacomini & White (2006), keeps sample size constant and prevents observations from old regimes from diluting current parameter estimates. It is particularly appropriate when structural breaks are plausible — as they are for Canadian savings data spanning the oil shocks of the late 1970s, the Great Moderation, and the 2008 financial crisis.

### 2.4 Evaluation Criterion

The primary metric is one-step-ahead holdout **Root Mean Square Error** (RMSE):

$$\text{RMSE} = \sqrt{\frac{1}{H} \sum_{t=1}^{H} (y_t - \hat{y}_t)^2}$$

where $H$ is the number of holdout observations. We also report in-sample (training) RMSE to diagnose overfitting: a model with low training RMSE but high holdout RMSE is likely overfitting.

### 2.5 ARIMAX Regressor Selection

Regressors were chosen by CCF analysis between each candidate's first difference and `dhh`. Only lags at which the regressor *leads* household savings (negative CCF lags) were retained, ensuring strict exogeneity by construction. This yielded:

- **Full**: `ldgpdi_lag1`, `ldcpi_lag2`, `ldgpdi_lag3`
- **Reduced**: `ldgpdi_lag1` (the only coefficient significant at $\alpha = 0.05$ in the initial window fit)

The coefficient significance test on the 1977–1996 initial window showed:

| Coefficient | p-value | Significant at 5%? |
|---|---|---|
| ar1 | 0.62 | No |
| ma1 | 0.029 | Yes |
| ldgpdi_lag1 | 0.041 | Yes |
| ldcpi_lag2 | 0.63 | No |
| ldgpdi_lag3 | 0.067 | No |

This motivates the reduced specification: the insignificant regressors inflate parameter estimation variance and may worsen out-of-sample performance.

### 2.6 Theoretical Model: Precautionary Income Model (PIM)

The PIM is derived from the CRRA Euler equation for intertemporal consumption:

$$C_t^{-\gamma} = \beta \, \mathbb{E}_t \left[ (1 + r_{t+1}) \, C_{t+1}^{-\gamma} \right]$$

Log-linearising and invoking the savings identity, the precautionary saving motive maps onto observable proxies:

| Theory variable | Proxy | Justification |
|---|---|---|
| Income growth $\Delta y_t$ | `ldgdi_lag1` | Log-differenced GDI |
| Real interest rate $r_t$ | `dinterest_lag1` − `ldcpi_lag1` | Fisher identity |
| Income uncertainty | `dunemp_lag1` | Labour market risk |

This yields a four-regressor ARIMAX specification with all regressors lagged by one quarter.

## 3. Results

### 3.1 Fixed-Parameter One-Step Forecasts (Notebook 02)

In the fixed-parameter evaluation, each model is fitted once on the 1977 Q1 – 2009 Q4 training set and evaluated over the 40-quarter holdout (2010 Q1 – 2019 Q4). The complete RMSE ranking is computed dynamically and stored in `hh_rmse_tbl.rds`. Key observations:

- **Persistence is hard to beat.** The naïve random-walk forecast sets a challenging benchmark. Many ARIMA and SARIMA specifications cluster around similar RMSE values, confirming the near-random-walk nature of household savings.
- **ARIMA/SARIMA models** provide marginal improvements over persistence when the order is correctly specified. The best ARIMA and SARIMA are selected programmatically by minimum holdout RMSE (dynamic labels stored as `best_arima_label` and `best_sarima_label`).
- **ARIMAX models** with CCF-justified regressors occasionally beat pure ARIMA, but the gains are modest. The reduced ARIMAX (using only `ldgpdi_lag1`) tends to match or outperform the full ARIMAX, consistent with the coefficient significance analysis.
- **The PIM specification** offers a theory-grounded alternative. Its RMSE is competitive with the data-driven ARIMAX, suggesting that the precautionary saving motive has some empirical relevance for Canadian household savings.
- **Exponential smoothing** models (SES, Holt, HW Additive) generally fall between the simple benchmarks and the ARIMA family.

A scatter matrix of forecast vectors reveals the degree of agreement among models. High pairwise correlation (particularly between Persistence and the ARIMA-family forecasts) confirms that most models are capturing similar dynamics — the unpredictable random-walk-like component dominates the signal.

### 3.2 Rolling Window Forecasts (Notebook 03)

In the rolling window evaluation, models are re-estimated at every holdout step on a fixed window of $R = 80$ quarters, producing 92 one-step-ahead forecasts (1997 Q1 – 2019 Q4).

**Fixed vs Rolling Comparison.** For each ARIMA, SARIMA, and ARIMAX specification, we compare the RMSE from the fixed-parameter scheme against the rolling window scheme. The percentage improvement is:

$$\text{Improvement (\%)} = \frac{\text{RMSE}_{\text{Fixed}} - \text{RMSE}_{\text{Rolling}}}{\text{RMSE}_{\text{Fixed}}} \times 100$$

A positive value indicates that the rolling window outperforms fixed parameters. Key observations:

- **Re-estimation generally helps.** The rolling window allows model parameters to adapt to regime changes, which is valuable over the extended 92-quarter holdout that spans the Asian financial crisis, the dot-com bust, and the 2008 recession.
- **The best rolling window model** is identified programmatically by minimum RMSE across all model-method combinations. It typically comes from the ARIMA or SARIMA family, reflecting the limited predictive gain from exogenous regressors.
- **ARIMAX rolling window evaluation** is informative for regressor relevance: if the ARIMAX rolling RMSE does not beat the pure ARIMA rolling RMSE, the exogenous variables should be dropped in favour of parsimony.
- **The reduced ARIMAX** continues to match or beat the full ARIMAX under the rolling scheme, reinforcing the insignificance of `ldcpi_lag2` and `ldgpdi_lag3`.

### 3.3 Holdout Visualisation

The holdout forecast plots show the actual household net savings against the predicted values from selected models. In both the fixed and rolling evaluations, the forecasts from ARIMA-family models track the actual series closely but with a noticeable lag — characteristic of persistence-like behaviour. The HW Additive model occasionally captures seasonal swings that the non-seasonal ARIMA misses, but this does not consistently lower RMSE.

## 4. Discussion

### 4.1 Near-Random-Walk Dynamics

The central finding is that Canadian household net savings behave approximately as a random walk with drift. The IMA(1,1) structure of the saving rate — identified by `auto.arima` — implies that shocks are partially absorbed within one quarter, but the transient MA(1) component offers limited forecasting leverage. This is consistent with the permanent income hypothesis: households smooth consumption over income shocks, which makes savings the residual absorber of unpredictable changes in income and expenditure.

### 4.2 Limited Value of Exogenous Regressors

Despite the theoretical appeal of macroeconomic covariates (GPDI as a leading indicator, CPI as an inflation proxy), the empirical gains from ARIMAX over pure ARIMA are modest. The CCF-based selection identifies statistically significant leading relationships, but these relationships are too weak to substantially reduce forecast error. The coefficient significance analysis confirms that only `ldgpdi_lag1` is significant at the 5% level in the initial window; the remaining regressors add estimation noise.

### 4.3 Structural Breaks and the Case for Rolling Windows

The rolling window is well-suited to this dataset because Canadian household savings have experienced multiple regime shifts: the high-saving environment of the late 1970s and early 1980s, the declining trend through the 1990s, and the post-2008 stabilisation. A fixed expanding window would average over all these regimes, potentially diluting the parameter estimates relevant to the current regime. The rolling window, by contrast, discards observations older than 20 years and allows the model to adapt.

### 4.4 Practical Implications

For a practitioner forecasting Canadian household savings one quarter ahead, the key takeaway is:

1. **Start with persistence.** The naïve forecast is a strong baseline.
2. **Consider an IMA(1,1) or ARIMA(1,1,1).** The modest MA(1) component provides a small forecasting edge.
3. **Add regressors cautiously.** Only `ldgpdi_lag1` is robustly significant; including more regressors risks overfitting.
4. **Re-estimate regularly.** A rolling window of 20 years provides enough data for stable estimation while adapting to structural change.

## 5. Contributions

*[To be completed.]*

## 6. References

- Box, G. E. P., Jenkins, G. M., Reinsel, G. C., & Ljung, G. M. (2015). *Time Series Analysis: Forecasting and Control* (5th ed.). Wiley.
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

Readers are referred to the rendered notebook outputs for the precise numerical values, as these depend on the stochastic `auto.arima` selection and are best consumed in their original tabular format.
