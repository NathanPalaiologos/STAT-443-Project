library(tidyverse)
library(tsibble)

dat <- read_csv("data/macro_panel_wide_raw.csv")|>
  mutate(quarterly_date = as.Date(quarterly_date)) |>
  arrange(quarterly_date)|>
  filter(quarterly_date <= as.Date("2019-12-31")) # Filter to the desired date range

# Create variables according to Forecast Plan.md.
# Note this is just my implementtion of the plan.
# It may not be 100% accurate, or someone else may already have done this

dat2 <- dat %>%
  arrange(quarterly_date) %>%
  mutate(
    # first differences of targets (consistent with macro_df_diff.rds naming)
    d_hh = household_net_savings - lag(household_net_savings),
    # first differences of candidate macro variables
    ldgdi      = log(gross_domestic_income) - lag(log(gross_domestic_income)),
    dunemp    = unemployment_rate - lag(unemployment_rate),
    dinterest = market_interest_rates - lag(market_interest_rates),
    ldcpi      = log1p(cpi_inflation_indicator) - log1p(lag(cpi_inflation_indicator)),
    ldgpdi     = log1p(gross_private_domestic_investment) - log1p(lag(gross_private_domestic_investment))
  )

glimpse(dat2)
# differenced variables misses first obs; lagged differenced variables missed first two obs.

#The function below is copied from teacher's R code.

# Function for variogram as defined in 
# Bisgaard and Kulahci (2011). Time Series Analysis and Forecasting by Example, Wiley
variogram = function(y, lagmax=10, iprint=F)
{ G = rep(1,lagmax)
n = length(y)
if(lagmax>n) { lagmax = n-2 }
y1 = y[-1]; y2 = y[-n]
d1 = y1-y2; denom = var(d1)
for(k in 2:lagmax)
{ y1 = y[(k+1):n]; y2 = y[1:(n-k)]
dk = y1-y2
numer = var(dk)
G[k] = numer/denom
}
# H is the variogram assuming a stationary time series
H = rep(1,lagmax)
ac = c(acf(y,plot=F,lag.max=lagmax)$acf)
H = (1-ac[-1])/(1-ac[2])
if(iprint) 
{ print(cbind(G,H)) 
  mx = max(cbind(G,H))
  matplot(1:lagmax,cbind(G,H),ylim=c(0,mx),ylab="variogram", xlab="lag")
  abline(h=1)
}
list(G=G, H=H)
}

# get variable names for all the numeric variables that need variogram check
var_names <- names(dat2)[sapply(dat2, is.numeric)]


vg_results <- vector("list", length(var_names))
names(vg_results) <- var_names

# loop through each variable and make variograms
for (v in var_names) {
  y <- dat2[[v]]
  y_clean <- y[!is.na(y)]
  
  vg_results[[v]] <- variogram(y_clean, lagmax = 10, iprint = FALSE)
}

# -------------------------------------------------------------------
# Separate level (undifferenced) vs differenced variables
# -------------------------------------------------------------------
diff_prefixes <- c("d_", "dhh", "dcorp", "ld", "dunemp", "dinterest")

is_diff <- sapply(var_names, function(v) {
  any(sapply(diff_prefixes, function(p) startsWith(v, p)))
})

level_vars <- var_names[!is_diff]
diff_vars  <- var_names[is_diff]

# Helper to plot a group of variograms into a PNG
plot_variogram_panel <- function(vars, vg_list, filepath, title_prefix,
                                  col_g = "#2166AC", col_h = "#B2182B") {
  n <- length(vars)
  n_cols  <- 3
  n_rows  <- ceiling(n / n_cols)

  png(filepath,
      width  = 1200,
      height = n_rows * 350,
      res    = 150)

  par(mfrow = c(n_rows, n_cols),
      mar   = c(4, 4, 3, 1),
      mgp   = c(2, 0.7, 0),
      cex.main = 0.95,
      cex.lab  = 0.85,
      cex.axis = 0.8)

  for (v in vars) {
    res <- vg_list[[v]]
    lags <- seq_along(res$G)
    mx <- max(c(res$G, res$H), na.rm = TRUE)

    matplot(lags, cbind(res$G, res$H),
            type = "b", pch = c(16, 17), lty = c(1, 2),
            col = c(col_g, col_h),
            xlab = "Lag", ylab = "Variogram",
            ylim = c(0, mx * 1.1),
            main = v)
    abline(h = 1, lty = 3, col = "grey50")
    legend("topleft",
           legend = c("G (empirical)", "H (stationary)"),
           lty = c(1, 2), pch = c(16, 17),
           col = c(col_g, col_h),
           bty = "n", cex = 0.7)
  }

  dev.off()
  message("Saved: ", filepath)
}

# Plot level variables
plot_variogram_panel(level_vars, vg_results,
                     "output/variograms_levels.png",
                     "Level Variables")

# Plot differenced variables
plot_variogram_panel(diff_vars, vg_results,
                     "output/variograms_differenced.png",
                     "Differenced Variables")

# Also keep a combined plot for backward compatibility
all_vars <- c(level_vars, diff_vars)
plot_variogram_panel(all_vars, vg_results,
                     "output/variograms_all.png",
                     "All Variables")