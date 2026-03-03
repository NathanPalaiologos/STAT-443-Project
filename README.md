# STAT 443 Project — Forecasting Sectoral National Net Savings in Canada

**Authors:** Kaiyan Zhang, Junhao Wen, Shisheng Xiao, Jingwen Leng, Boyan Yu

---

## Table of Contents

1. [Motivation & Problem](#1-motivation--problem)
2. [Data](#2-data)
3. [Models and Evaluation Metrics](#3-models-and-evaluation-metrics)
4. [Implementation Plan](#4-implementation-plan)
5. [Repository Structure](#5-repository-structure)
6. [Getting Started — renv Setup](#6-getting-started--renv-setup)
7. [Contributing via Pull Requests](#7-contributing-via-pull-requests)

---

## 1. Motivation & Problem

Net savings — calculated as net disposable income minus final consumption expenditure — represent the funds available for investment and future consumption. Modern economists regard them as reflecting a *precautionary motive* in economic agents' decision-making under future uncertainty. Forecasting sectoral saving patterns therefore provides a window into confidence, risk exposure, and the effectiveness of stabilization policy.

This project focuses on the net savings of **Canadian households and firms** to examine which sector is more sensitive to interest rate shocks, unexpected inflation, and labour market frictions. By comparing predictability across sectors and forecasting rules, we assess:

1. Who builds and draws down buffers over the business cycle;
2. Who is more prone to behaviour changes after policy shocks or structural breaks; and
3. How forward- or backward-looking agents appear based on autocorrelation and cross-correlation patterns.

---

## 2. Data

All data are sourced from **Statistics Canada** and the **Bank of Canada**, harmonized to a **quarterly frequency** over the period **1977 Q1 – 2025 Q3** (191 observations).

| Variable | Source |
|---|---|
| Current & Capital Accounts (net savings) | [Table 36-10-0111-01](https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=3610011101) |
| Gross Domestic Income | [Table 36-10-0105-01](https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=3610010501) |
| Unemployment Rate | [Table 14-10-0287-03](https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1410028703) |
| Market Interest Rates | [Table 10-10-0139-01](https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1010013901) |
| CPI (Inflation Indicator) | [Table 18-10-0004-13](https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=1810000413) |
| Household Consumption Expenditure | [Table 36-10-0107-01](https://www150.statcan.gc.ca/t1/tbl1/en/tv.action?pid=3610010701) |

**Pre-processing steps:**

- Exploratory data analysis (time-series plots, ACF/PACF)
- Log-differencing to achieve stationarity where needed
- Stationarity verification via variogram plots and the Augmented Dickey–Fuller test

---

## 3. Models and Evaluation Metrics

### Benchmark rules

| Rule | Description |
|---|---|
| Average forecast | Forecast equals the historical mean |
| Persistence forecast | Forecast equals the most recent observed value |

### Candidate models

We fit models from the **exponential smoothing** and **ARIMA** families, including:

- Specifications inspired by structural economic models (e.g., the Permanent Income Hypothesis)
- Specifications selected by algorithmic search (Auto ARIMA, heuristic grid search)

### One-step-ahead validation RMSE (primary metric)

Forecasts are evaluated on a held-out validation set $\{T_0+1, \ldots, T\}$:

$$\mathrm{RMSE}_{\mathrm{OOB},1} = \sqrt{\frac{1}{T-T_0}\sum_{t=T_0+1}^{T}\!\left(y_t - \hat{y}_{t\mid t-1}\right)^2}$$

### Rolling-window RMSE (time permitting)

Time-series cross-validation via rolling one-step forecasts over a window of length $w$:

$$\mathrm{RMSE}_{\mathrm{Rolling},1}(w) = \sqrt{\frac{1}{T-w}\sum_{t=w+1}^{T}\!\left(y_t - \hat{y}^{(w)}_{t\mid t-1}\right)^2}$$

---

## 4. Implementation Plan

The project is developed in **RStudio** with version control managed through **GitHub**.

- All data ingestion, cleaning, and transformation steps are implemented as reproducible R scripts with a clear, modular structure.
- Common utilities (plotting themes, metric functions, forecast-rule helpers) are centralized in helper scripts to ensure consistency and reproducibility.
- The full workflow runs end-to-end from raw inputs to final tables and figures with minimal manual intervention.

---

## 5. Repository Structure

```
STAT-443-Project/
├── renv/               # renv package library (auto-generated, not committed)
├── renv.lock           # Lockfile pinning exact package versions
├── .Rprofile           # Auto-activates renv on project load
├── STAT-443-Project.Rproj
└── README.md
```

> Additional folders (e.g., `data/`, `R/`, `output/`) will be added as the project develops.

---

## 6. Getting Started — renv Setup

This project uses [**renv**](https://rstudio.github.io/renv/) to manage R package dependencies, ensuring every collaborator works with identical package versions.

### Prerequisites

- R ≥ 4.5.1
- RStudio (recommended) or any R environment

### Steps

1. **Clone the repository**

   ```bash
   git clone https://github.com/NathanPalaiologos/STAT-443-Project.git
   cd STAT-443-Project
   ```

2. **Open the project in RStudio**

   Double-click `STAT-443-Project.Rproj`.  
   `renv` will be activated automatically via `.Rprofile`.

3. **Restore the package library**

   In the R console, run:

   ```r
   renv::restore()
   ```

   This reads `renv.lock` and installs all packages at the pinned versions into the project-local library.

4. **Verify the environment** (optional)

   ```r
   renv::status()
   ```

   A message like *"The project is already synchronized with the lockfile."* confirms everything is in order.

### Keeping the lockfile up to date

After installing or updating any package during development, run:

```r
renv::snapshot()
```

Commit the updated `renv.lock` so that teammates can synchronize.

---

## 7. Contributing via Pull Requests

We follow a **feature-branch workflow**. Please do **not** push directly to `main`.

### Workflow

1. **Sync your local `main` with the remote**

   ```bash
   git checkout main
   git pull origin main
   ```

2. **Create a feature branch**

   Use a descriptive name, e.g.:

   ```bash
   git checkout -b feature/arima-household
   ```

3. **Make your changes**

   - Write clean, well-commented R code.
   - Keep each commit focused on a single logical change.
   - Run `renv::snapshot()` if you add or update packages, and commit the updated `renv.lock`.

4. **Push your branch**

   ```bash
   git push origin feature/arima-household
   ```

5. **Open a Pull Request on GitHub**

   - Navigate to the repository on GitHub and click **"Compare & pull request"**.
   - Give the PR a clear title and a short description of what was changed and why.
   - Request a review from at least one other team member.

6. **Address review feedback**

   Push additional commits to the same branch; the PR updates automatically.

7. **Merge**

   Once approved and all checks pass, merge using **"Squash and merge"** to keep the history clean, then delete the feature branch.

### Commit message style

```
<type>: <short summary>

Optional longer explanation.
```

Common types: `feat`, `fix`, `data`, `docs`, `refactor`, `test`.

Example: `feat: add Auto ARIMA model for household net savings`
 