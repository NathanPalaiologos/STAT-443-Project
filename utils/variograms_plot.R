library(tidyverse)
library(tsibble)

dat <- read_csv("data/macro_panel_wide_raw.csv")|>
  mutate(quarterly_date = as.Date(quarterly_date)) |>
  arrange(quarterly_date)|>
  filter(quarterly_date <= as.Date("2014-12-31")) # Filter to the desired date range

# Create variables according to Forecast Plan.md.
# Note this is just my implementtion of the plan.
# It may not be 100% accurate, or someone else may already have done this

dat2 <- dat %>%
  arrange(quarterly_date) %>%
  mutate(
    # saving rates (consistent with eda_df_nodiff.rds naming)
    hh_rate   = household_net_savings / gross_domestic_income,
    corp_rate = corporate_net_savings / gross_domestic_income,
    
    # first differences of targets (consistent with macro_df_diff.rds naming)
    d_hh = household_net_savings - lag(household_net_savings),
    d_corp = corporate_net_savings - lag(corporate_net_savings),
    dhh_rate   = hh_rate - lag(hh_rate),
    dcorp_rate = corp_rate - lag(corp_rate),
    
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

# make plots in a grid.
# Compute grid dimensions dynamically from the number of variables.
# output the plot in png so we can import in EDA.rmd
n_vars <- length(var_names)
n_cols  <- 3
n_rows  <- ceiling(n_vars / n_cols)

png("output/variograms_all.png",
    width  = 1200,
    height = n_rows * 400,
    res    = 150)

par(mfrow = c(n_rows, n_cols))

    
for (v in var_names) {
  res <- vg_results[[v]]
  matplot(seq_along(res$G), cbind(res$G, res$H),
            type = "b", pch = c(1, 2), lty = c(1, 2),
            xlab = "Lag", ylab = "Variogram", main = v)
  abline(h = 1, lty = 3)
  legend("topleft", legend = c("G", "H"),
           lty = c(1, 2), pch = c(1, 2), bty = "n", cex = 0.8)
}

dev.off()