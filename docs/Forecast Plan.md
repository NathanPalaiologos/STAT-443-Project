---
title: "Forecast Plan"
output: html_document
---

# Econometric Forecasting Plan: Household and Corporate Net Savings

This document synthesizes our empirical findings and outlines a structured econometric action plan for forecasting household and corporate net savings. It translates the theoretical discussions regarding information efficiency, habit formation, and simultaneous equation bias into a rigorous predictive framework designed to minimize one-step-ahead out-of-sample error.

## 1. Empirical Documentation and Baseline Hypotheses

Our initial univariate analysis utilizing automated ARIMA selection yielded distinct data generating processes for the household and corporate sectors.

### Household Sector Dynamics
The absolute household net savings series followed an ARIMA(1,1,1)(0,0,2)[4] process, which simplified to an IMA(1,1) process when transformed into a saving rate:
$$\Delta hr_t = \varepsilon_t + \theta_1 \varepsilon_{t-1}$$
where $\theta_1 = -0.4004$. 

**Hypothesis:** The household saving rate exhibits short-term mean reversion. Because household wage income is downwardly sticky and consumption is subject to habit formation, households cannot adjust their permanent consumption paths instantaneously. A transitory macroeconomic shock generates a contemporaneous spike in the saving rate, which is partially corrected in the subsequent quarter as behavioral inertia dissipates.

### Corporate Sector Dynamics
Both the absolute corporate savings and the corporate saving rate were identified as pure random walks, specifically ARIMA(0,1,0):
$$\Delta cr_t = \varepsilon_t$$

**Hypothesis:** Corporate cash flows and retained earnings are highly sensitive to real-time market fluctuations. The corporate sector efficiently and completely absorbs aggregate shocks into its new baseline saving rate within a single quarter. The strict white noise nature of $\Delta cr_t$ implies that internal autoregressive memory is insufficient for predictive modeling.

## 2. Resolving the Corporate Random Walk: Exogenous Integration

Relying on a univariate framework for corporate savings yields a naive forecast. To extract structural signal from the white noise residual $\varepsilon_t$, we must project the first-differenced series onto a predetermined exogenous information set $\mathcal{F}_{t-1}$.



The core issue with predicting corporate savings via macroeconomic covariates like Gross Domestic Income (GDI) and consumption is tautological accounting. Corporate savings is fundamentally linked to capital expenditure. Therefore, the most mathematically sound fix is to utilize Gross Private Domestic Investment (GPDI) and the cost of capital as leading indicators. 

Under the pecking order theory of corporate finance, firms prioritize internal funds (retained earnings) for investment. By lagging the structural regressors, we bypass simultaneous equation bias and establish strict exogeneity. The predictive specification for the corporate saving rate ($cr_t$) is:

$$\Delta cr_t = \alpha + \beta_1 \Delta GPDI_{t-1} + \beta_2 \Delta r_{t-1} + \beta_3 \Delta \pi_{t-1} + u_t$$

* **$\Delta GPDI_{t-1}$:** Changes in prior-quarter aggregate investment act as a proxy for corporate capital deployment cycles.
* **$\Delta r_{t-1}$:** The lagged first difference of market interest rates captures the shifting opportunity cost of retaining cash versus distributing dividends or servicing debt.
* **$\Delta \pi_{t-1}$:** Lagged CPI inflation proxies the nominal adjustment required to maintain real operating margins.

## 3. General ARIMAX Model Specification

To forecast the one-step-ahead household saving rate ($hr_t$) while avoiding the endogeneity inherent in national accounting identities, we specify a dynamic regression model utilizing strictly lagged covariates. 

Let $X_{t-1}$ represent the vector of predetermined macroeconomic state variables: the first differences of GDI, the unemployment rate, market interest rates, and CPI.

The ARIMAX specification to be estimated is:

$$\Delta hr_t = c + \sum_{i=1}^p \phi_i \Delta hr_{t-i} + \sum_{j=1}^k \gamma_j X_{t-1, j} + \varepsilon_t + \sum_{m=1}^q \theta_m \varepsilon_{t-m}$$

This structure guarantees that $E[\varepsilon_t | X_{t-1}] = 0$, satisfying the Gauss-Markov assumptions necessary for consistent parameter estimation. The inclusion of the lagged unemployment rate serves to theoretically ground the forecast in the precautionary savings motive.

## 4. Action Plan for Implementation

To achieve the primary objective of minimizing the one-step-ahead Root Mean Square Error (RMSE), we will execute the following methodological steps.

**Step 1: Stationarity and Integration Testing**
* Conduct Augmented Dickey-Fuller (ADF) and Phillips-Perron (PP) tests on all candidate variables (GDI, CPI, Unemployment, Interest Rates, GPDI).
* Apply the difference operator $\Delta$ to any $I(1)$ variables to ensure the entire system is strictly stationary $I(0)$, preventing spurious regression artifacts.

**Step 2: Feature Engineering and Selection**
* Construct the predetermined regressor matrix $X_{t-1}$ by lagging the stationary macroeconomic variables by exactly one quarter.
* Evaluate the inclusion of strictly forward-looking variables, such as a Consumer Confidence Index or the yield curve spread, to enhance the predictive power of the information set.



**Step 3: Rolling Origin Cross-Validation**
* Partition the quarterly dataset into an initial training window (e.g., 80% of the sample) and a testing holdout set.
* Implement a time-series cross-validation protocol. Train the ARIMAX model on observations $1$ through $T$, forecast $T+1$, and calculate the squared error.
* Expand the training window to $T+1$, re-estimate all $\phi$, $\theta$, and $\gamma$ parameters via Maximum Likelihood Estimation, and forecast $T+2$.

**Step 4: Hyperparameter Optimization**
* Iterate the cross-validation protocol across a constrained grid of autoregressive ($p$) and moving average ($q$) orders.
* Select the final model specification that strictly minimizes the average out-of-sample RMSE across all rolling windows, rather than relying exclusively on in-sample metrics like the Akaike Information Criterion. 

Would you like to proceed with setting up the code for the Augmented Dickey-Fuller tests to verify the integration order of the exogenous variables?