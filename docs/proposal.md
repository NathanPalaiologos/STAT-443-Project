# Project Pivot Proposal: Forecasting Canada's Unemployment Rate

## 1. Motivation for the Pivot

Our original target — household net savings — proved to be essentially a random walk. Neither ARIMA, SARIMA, ARIMAX with CCF-justified regressors, nor the precautionary-saving theoretical model could beat simple benchmarks (Persistence, Average) on the 20-quarter holdout. This is actually an economically meaningful null result: household savings absorb idiosyncratic shocks so efficiently that the best one-step-ahead predictor is simply the most recent observation.

We propose to pivot to **forecasting Canada's quarterly unemployment rate**, a series already present in our harmonized dataset (`unemployment_rate` / `unemp`). This pivot maximizes reuse of our existing data pipeline, train/test split, one-step-ahead forecast infrastructure, and — crucially — turns the covariates that *failed* to predict savings into covariates with strong *theoretical and empirical justification* for predicting unemployment.

## 2. Why Unemployment?

### 2.1 Statistical Forecastability

Unlike household net savings, the unemployment rate has several properties that favour time-series modelling:

| Property | hh_save | unemp |
|---|---|---|
| Bounded | No (can be negative) | Yes (0–100%) |
| Mean-reverting | Weak (≈ random walk) | Strong (NAIRU anchor) |
| Autocorrelation in levels | Slow decay (unit root) | Slow decay but *stationary* or near-stationary |
| Autocorrelation after differencing | Near white-noise | Significant AR structure |
| Seasonal pattern | Weak | Moderate (quarterly hiring cycles) |

The unemployment rate is well-known in the forecasting literature to exhibit high persistence (large AR(1) coefficient ≈ 0.9) with a mean-reverting tendency toward the Non-Accelerating Inflation Rate of Unemployment (NAIRU). This means ARMA/ARIMA models can capture meaningful internal dynamics — a prerequisite for beating naïve benchmarks.

### 2.2 Our CCF Evidence Already Supports This

In our EDA notebook (`01_EDA.rmd`), we computed cross-correlation functions between all differenced macro variables. The CCF plots show:

- **`ldgdi` vs `dunemp`**: GDI growth *leads* unemployment changes at lags −1 to −3 (Okun's Law in action).
- **`ldgpdi` vs `dunemp`**: Investment growth leads unemployment changes — firms hire after investing.
- **`dhh` vs `dunemp`**: Household saving changes *lead* unemployment at negative lags. This was an inconvenient finding when savings was the response; now it becomes a *useful predictor*.
- **`ldcons` vs `dunemp`**: Consumption growth leads unemployment — demand-side channel.

In other words, the regressors that showed no power for predicting savings have **demonstrated leading relationships with unemployment** in the very same EDA we already produced.

### 2.3 Rich Theoretical Grounding

The economic theory connecting our covariates to unemployment is deep and well-established:

1. **Okun's Law** ($\Delta u_t = \alpha + \beta \Delta \ln Y_t + \varepsilon_t$): Output growth (GDI) is inversely related to unemployment changes. This is one of the most robust empirical regularities in macroeconomics.

2. **The Phillips Curve** ($\pi_t = \pi_t^e - \gamma(u_t - u^*) + \nu_t$): Inflation and unemployment are jointly determined. Lagged CPI changes carry predictive information about the unemployment trajectory.

3. **Investment–Employment Channel**: Gross Private Domestic Investment (GPDI) leads employment because firms commit capital before hiring labour. A rise in GPDI today signals lower unemployment next quarter.

4. **Monetary Policy Transmission**: Changes in market interest rates affect unemployment with a lag through the credit channel — higher rates slow borrowing, investment, hiring, and eventually raise unemployment.

5. **Precautionary Behaviour Feedback**: The precautionary savings motive we studied earlier implies that household savings rise *in anticipation* of unemployment risk. This means lagged `dhh` is a forward-looking signal of labour market deterioration — a theoretically justified and non-obvious predictor.

## 3. Proposed Empirical Specification

### 3.1 Response Variable

$$y_t = \text{unemp}_t$$

We model the unemployment rate in levels (or first-differenced if unit-root tests indicate non-stationarity on the training sample). Given the quarterly frequency and boundedness, the level series is likely stationary or trend-stationary, which simplifies modelling.

### 3.2 Candidate Model Hierarchy

We follow the same progressive complexity structure as in the household savings notebook:

| Stage | Model | Description |
|---|---|---|
| **Benchmarks** | Average, Persistence, Seasonal Average, Seasonal Persistence | Naïve baselines |
| **Exponential Smoothing** | SES, Holt Linear, HW Additive | Capture level, trend, and seasonal components |
| **Univariate ARIMA** | ARIMA($p, d, q$), SARIMA($p,d,q$)($P,D,Q$)[4] | Internal dynamics of the unemployment series |
| **ARIMAX (data-driven)** | ARIMAX with CCF-selected regressors | Lagged `ldgdi`, `ldgpdi`, `ldcons` at CCF-optimal lags |
| **ARIMAX (Okun's Law)** | $\Delta u_t = \alpha + \beta_1 \Delta \ln \text{GDI}_{t-1} + \beta_2 \Delta r_{t-1} + \phi(B)\varepsilon_t$ | Theory-driven: Okun's law + monetary policy channel |
| **ARIMAX (Phillips + Okun hybrid)** | $\Delta u_t = \alpha + \beta_1 \Delta \ln \text{GDI}_{t-1} + \beta_2 \Delta \ln \text{CPI}_{t-1} + \beta_3 \Delta \text{hh\_save}_{t-1} + \phi(B)\varepsilon_t$ | Combines output gap, inflation expectations, and precautionary savings signal |

### 3.3 Exogenous Variables and Justification

| Regressor | Transformation | Lag | Theoretical Justification |
|---|---|---|---|
| `gdi` | $\Delta \ln \text{GDI}_{t-k}$ | 1–2 | Okun's Law: output growth reduces unemployment |
| `gpdi` | $\Delta \ln \text{GPDI}_{t-k}$ | 1–3 | Investment precedes hiring |
| `interest` | $\Delta r_{t-k}$ | 1–2 | Monetary transmission: rate hikes slow employment |
| `cpi` | $\Delta \ln \text{CPI}_{t-k}$ | 1 | Phillips Curve: inflation–unemployment tradeoff |
| `hh_save` | $\Delta \text{hh\_save}_{t-k}$ | 1–2 | Precautionary saving motive as a leading indicator of labour market stress |
| `cons` | $\Delta \ln \text{cons}_{t-k}$ | 1 | Demand channel: consumption drives employment |

Lag selection will be confirmed via CCF analysis on the training sample, following the same protocol as in `02_Household_Forecast.Rmd`.

## 4. What We Can Reuse

| Asset | Reuse? | Notes |
|---|---|---|
| `data_ingest.R` | ✅ Fully | Same CANSIM tables, same harmonization |
| `macro_panel_wide_raw.csv` | ✅ Fully | Unemployment already in the dataset |
| `01_EDA.rmd` | ✅ Mostly | CCF plots already computed; add unemployment-focused diagnostics |
| Train/test split (2015 Q1 cutoff) | ✅ Fully | Same 20-quarter holdout |
| `one_step_model.R` | ✅ Fully | All benchmark, ESM, ARIMA, ARIMAX functions work on any univariate ts |
| `02_Household_Forecast.Rmd` | 🔄 Template | Clone structure, swap `hh_save` → `unemp` |
| Forecast Plan theoretical framework | 🔄 Adapt | Replace precautionary savings theory with Okun + Phillips |

## 5. Narrative Arc for the Report

1. **Introduction**: Forecasting unemployment is a first-order policy question. The Bank of Canada and federal government rely on unemployment forecasts for monetary and fiscal policy respectively.

2. **EDA & Stationarity**: Demonstrate that `unemp` (or `dunemp`) is stationary and has meaningful autocorrelation structure — unlike household savings.

3. **Benchmarks & Exponential Smoothing**: Establish the bar to beat.

4. **ARIMA / SARIMA**: Exploit the strong internal persistence of unemployment.

5. **ARIMAX — Okun's Law**: Test whether lagged output growth improves forecasts, grounded in one of the most established relationships in macroeconomics.

6. **ARIMAX — Phillips + Precautionary Hybrid**: A novel empirical specification that uses lagged household savings change as a forward-looking fear index. This ties back to our *original* project and reframes the random-walk finding as a feature: savings changes encode forward-looking labour market expectations that contain predictive signal for unemployment.

7. **Model Comparison & Conclusion**: RMSE ranking, discussion of whether theory-guided models outperform pure ARIMA, and economic interpretation.

## 6. Risk Assessment

| Risk | Mitigation |
|---|---|
| `unemp` is also near-random-walk | Unlikely given known persistence; but if so, the Phillips + Okun ARIMAX still provides structural insights even without RMSE gains |
| Small holdout (20 quarters) | Same constraint as before; use one-step-ahead expanding-window evaluation consistently |
| Multicollinearity among regressors | Use CCF to pre-select; compare nested ARIMAX specifications |
| Overfitting with too many lags/regressors | Keep model order low ($p + q \leq 3$); limit exogenous variables to 2–3 per specification |

## 7. Recommended Next Steps

1. **Copy** `02_Household_Forecast.Rmd` → `04_Unemployment_Forecast.Rmd`
2. **Swap** response variable from `hh_save` to `unemp` and re-run stationarity checks
3. **Re-examine CCFs** with `dunemp` as the response to confirm lag structures on the training set
4. **Fit the model hierarchy** (benchmarks → ESM → ARIMA → ARIMAX) using the same one-step-ahead protocol
5. **Write the Okun's Law and Phillips Curve ARIMAX** specifications as the theoretical models
6. **Compare RMSE** and select the best model
