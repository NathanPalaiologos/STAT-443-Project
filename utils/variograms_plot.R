library(tidyverse)
library(tsibble)

dat <- read_csv("data/macro_panel_wide_raw.csv")
glimpse(dat)

# Create variables according to Forecast Plan.md.
# Note this is just my implementtion of the plan.
# It may not be 100% accurate, or someone else may already have done this

dat2 <- dat %>%
  arrange(quarterly_date) %>%
  mutate(
    # saving rates
    hr = household_net_savings / gross_domestic_income,
    cr = corporate_net_savings / gross_domestic_income,
    
    # first differences of targets
    d_hr = hr - lag(hr),
    d_cr = cr - lag(cr),
    
    # first differences of candidate macro variables
    d_gdi   = gross_domestic_income - lag(gross_domestic_income),
    d_unemp = unemployment_rate - lag(unemployment_rate),
    d_rate  = market_interest_rates - lag(market_interest_rates),
    d_cpi   = cpi_inflation_indicator - lag(cpi_inflation_indicator),
    d_gpdi  = gross_private_domestic_investment - lag(gross_private_domestic_investment)
  ) %>%
  mutate(
    # one-quarter lagged predetermined information set
    L1_d_gdi   = lag(d_gdi, 1),
    L1_d_unemp = lag(d_unemp, 1),
    L1_d_rate  = lag(d_rate, 1),
    L1_d_cpi   = lag(d_cpi, 1),
    L1_d_gpdi  = lag(d_gpdi, 1)
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
# We have 24 variograms. So will use 8 * 3 grids
# output the plot in png so we can import in EDA.rmd
png("output/variograms_all.png",
    width = 1200,
    height = 3200,   # ~ 8/3 ratio
    res = 150)

par(mfrow = c(8, 3))

    
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