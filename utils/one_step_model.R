# =========================================================
# one_step_model.R
# Stat 443 one-step-ahead forecast functions (instructor-style convention)
# Note: Return structures may vary by function; see each function for details.
# =========================================================

# 1. SIMPLE BASELINE MODELS

# PERSISTENCE MODEL
# (train, holdout)

persist_fc <- function(train, holdout, iprint = FALSE) {
  n <- length(train)
  nholdout <- length(holdout)
  
  mse <- 0
  fc_vec <- rep(NA, nholdout)
  
  fc <- train[n]
  yt <- holdout[1]
  fc_vec[1] <- fc
  fcerror <- yt - fc
  mse <- mse + fcerror^2
  
  if (nholdout >= 2) {
    for (i in 2:nholdout) {
      yt <- holdout[i]
      fc <- holdout[i - 1]
      fc_vec[i] <- fc
      fcerror <- yt - fc
      mse <- mse + fcerror^2
    }
  }
  
  rmse <- sqrt(mse / nholdout)
  return(list(fc = fc_vec, rmse = rmse))
}

# IID/AVERAGE MODEL
# (train, holdout)  

iid_fc <- function(train, holdout, iprint = FALSE) {
  nholdout <- length(holdout)
  avg <- mean(train)
  
  mse <- 0
  fc_vec <- rep(NA, nholdout)
  
  for (i in 1:nholdout) {
    yt <- holdout[i]
    fc <- avg
    fc_vec[i] <- fc
    fcerror <- yt - fc
    mse <- mse + fcerror^2
  }
  
  rmse <- sqrt(mse / nholdout)
  return(list(fc = fc_vec, rmse = rmse))
}

# Simple LINEAR MODEL
# (train, holdout, bcoef)

linear_fc = function(train, holdout, bcoef, iprint = FALSE) {
  n = length(train)
  nholdout = length(holdout)
  
  b0 = bcoef[1]
  b1 = bcoef[2]
  
  mse = 0
  fc_vec = rep(NA, nholdout)
  
  yt = holdout[1]
  fc = b0 + b1 * train[n]
  fc_vec[1] = fc
  fcerror = yt - fc
  mse = mse + fcerror^2
  
  if (nholdout >= 2) {
    for (i in 2:nholdout) {
      fc = b0 + b1 * holdout[i - 1]
      yt = holdout[i]
      fc_vec[i] = fc
      fcerror = yt - fc
      mse = mse + fcerror^2
    }
  }
  
  rmse = sqrt(mse / nholdout)
  return(list(fc = fc_vec, rmse = rmse))
}




# 2. EXPONENTIAL SMOOTHING

# SIMPLE EXPONENTIAL SMOOTHING (SES)
# (train, holdout, alpha, level)

esm_fc = function(train, holdout, alpha, level, iprint = FALSE) {
  nholdout = length(holdout)
  
  sse = 0
  fc_vec = rep(NA, nholdout)
  
  fc = level
  yt = holdout[1]
  fc_vec[1] = fc
  fcerror = yt - fc
  sse = sse + fcerror^2
  ellprev = level
  
  if (nholdout >= 2) {
    for (i in 2:nholdout) {
      ellnew = alpha * holdout[i - 1] + (1 - alpha) * ellprev
      fc = ellnew
      fc_vec[i] = fc
      yt = holdout[i]
      fcerror = yt - fc
      sse = sse + fcerror^2
      ellprev = ellnew
    }
  }
  
  rmse = sqrt(sse / nholdout)
  return(list(fc = fc_vec, rmse = rmse))
}


# Holt linear exponential smoothing
# (train, holdout, alpha, beta, level, slope)

lholt_fc = function(train, holdout, alpha, beta, level, slope, iprint = FALSE) {
  nholdout = length(holdout)
  
  sse = 0
  fc_vec = rep(NA, nholdout)
  
  fc = level + slope
  yt = holdout[1]
  fc_vec[1] = fc
  fcerror = yt - fc
  sse = sse + fcerror^2
  
  ellprev = level
  bprev = slope
  
  if (nholdout >= 2) {
    for (i in 2:nholdout) {
      ellnew = alpha * holdout[i - 1] + (1 - alpha) * fc
      bnew = beta * (ellnew - ellprev) + (1 - beta) * bprev
      fc = ellnew + bnew
      fc_vec[i] = fc
      yt = holdout[i]
      fcerror = yt - fc
      sse = sse + fcerror^2
      ellprev = ellnew
      bprev = bnew
    }
  }
  
  rmse = sqrt(sse / nholdout)
  return(list(fc = fc_vec, rmse = rmse))
}

# Damped Holt linear exponential smoothing

dholt_fc = function(train, holdout, alpha, beta, level, slope, phi, iprint = FALSE) {
  nholdout = length(holdout)
  
  sse = 0
  fc_vec = rep(NA, nholdout)
  
  fc = level + phi * slope
  yt = holdout[1]
  fc_vec[1] = fc
  fcerror = yt - fc
  sse = sse + fcerror^2
  
  ellprev = level
  bprev = slope
  
  if (nholdout >= 2) {
    for (i in 2:nholdout) {
      ellnew = alpha * holdout[i - 1] + (1 - alpha) * fc
      bnew = beta * (ellnew - ellprev) + (1 - beta) * bprev
      fc = ellnew + phi * bnew
      fc_vec[i] = fc
      yt = holdout[i]
      fcerror = yt - fc
      sse = sse + fcerror^2
      ellprev = ellnew
      bprev = bnew
    }
  }
  
  rmse = sqrt(sse / nholdout)
  return(list(fc = fc_vec, rmse = rmse))
}


# 3. SEASONAL MODELS

# SEASONAL PERSISTENCE
# (train, holdout, d=12)

persistbymonth_fc = function(train, holdout, d = 12, iprint = FALSE) {
  n = length(train)
  nholdout = length(holdout)
  
  if (n < d) stop("train must contain at least one full seasonal cycle")
  
  last_season = train[(n - d + 1):n]
  
  mse = 0
  fc_vec = rep(NA, nholdout)
  
  for (i in 1:nholdout) {
    pos = ((i - 1) %% d) + 1
    
    if (i <= d) {
      fc = last_season[pos]
    } else {
      fc = holdout[i - d]
    }
    
    yt = holdout[i]
    fc_vec[i] = fc
    fcerror = yt - fc
    mse = mse + fcerror^2
  }
  
  rmse = sqrt(mse / nholdout)
  return(list(fc = fc_vec, rmse = rmse))
}

# SEASONAL AVERAGE
# (train, holdout, d=12)

iidbymonth_fc = function(train, holdout, iprint = FALSE, d = 12) {
  n = length(train)
  nholdout = length(holdout)
  
  if (n %% d != 0) warning("train length is not a multiple of d")
  
  season_means = rep(NA, d)
  for (j in 1:d) {
    season_means[j] = mean(train[seq(j, n, by = d)])
  }
  
  mse = 0
  fc_vec = rep(NA, nholdout)
  
  for (i in 1:nholdout) {
    pos = ((i - 1) %% d) + 1
    fc = season_means[pos]
    yt = holdout[i]
    fc_vec[i] = fc
    fcerror = yt - fc
    mse = mse + fcerror^2
  }
  
  rmse = sqrt(mse / nholdout)
  return(list(fc = fc_vec, rmse = rmse))
}

# Calculate Annual Average (Seasonal Average by Year)
iidbyyear_fc = function(train, holdout, periodicity = 12, iprint = FALSE) {
  n = length(train)
  nholdout = length(holdout)
  
  if (n %% periodicity != 0) warning("train length is not a multiple of periodicity")
  
  year_means = rep(NA, n / periodicity)
  for (j in 1:(n / periodicity)) {
    year_means[j] = mean(train[((j - 1) * periodicity + 1):(j * periodicity)])
  }
  
  mse = 0
  fc_vec = rep(NA, nholdout)
  
  for (i in 1:nholdout) {
    pos = ((i - 1) %/% periodicity) + 1
    fc = year_means[pos]
    yt = holdout[i]
    fc_vec[i] = fc
    fcerror = yt - fc
    mse = mse + fcerror^2
  }
  
  rmse = sqrt(mse / nholdout)
  
  if (iprint) print(cbind(holdout = holdout, fc = fc_vec))
  
  return(list(fc = fc_vec, rmse = rmse))
}


# Winters additive seasonal exponential smoothing
# (train, holdout, alpha, beta, gamma, level, slope, season)


aseason_fc = function(train, holdout, alpha, beta, gamma, level, slope, season,
                      iprint = FALSE) {
  
  d = length(season)
  nholdout = length(holdout)
  
  sse = 0
  fc_vec = rep(NA, nholdout)
  
  seas = as.numeric(season)
  
  # first forecast
  fc = level + slope + seas[1]
  yt = holdout[1]
  fc_vec[1] = fc
  fcerror = yt - fc
  sse = sse + fcerror^2
  
  ellprev = level
  bprev = slope
  
  if (nholdout >= 2) {
    for (i in 2:nholdout) {
      
      prev_pos = ((i - 2) %% d) + 1
      curr_pos = ((i - 1) %% d) + 1
      
      sold = seas[prev_pos]
      
      ellnew = alpha * (holdout[i - 1] - sold) + (1 - alpha) * (ellprev + bprev)
      bnew   = beta * (ellnew - ellprev) + (1 - beta) * bprev
      snew   = gamma * (holdout[i - 1] - ellnew) + (1 - gamma) * sold
      
      seas[prev_pos] = snew
      
      fc = ellnew + bnew + seas[curr_pos]
      
      fc_vec[i] = fc
      yt = holdout[i]
      fcerror = yt - fc
      sse = sse + fcerror^2
      
      ellprev = ellnew
      bprev = bnew
    }
  }
  
  rmse = sqrt(sse / nholdout)
  
  if (iprint) print(cbind(holdout = holdout, fc = fc_vec))
  
  return(list(fc = fc_vec, rmse = rmse))
}

# Winters multiplicative seasonal exponential smoothing
# (train, holdout, alpha, beta, gamma, level, slope, season)

mseason_fc = function(train, holdout, alpha, beta, gamma, level, slope, season,
                      iprint = FALSE) {
  d = length(season)
  nholdout = length(holdout)
  
  sse = 0
  fc_vec = rep(NA, nholdout)
  
  seas = as.numeric(season)
  
  fc = (level + slope) * seas[1]
  yt = holdout[1]
  fc_vec[1] = fc
  fcerror = yt - fc
  sse = sse + fcerror^2
  
  ellprev = level
  bprev = slope
  
  if (nholdout >= 2) {
    for (i in 2:nholdout) {
      prev_pos = ((i - 2) %% d) + 1
      curr_pos = ((i - 1) %% d) + 1
      
      sold = seas[prev_pos]
      
      ellnew = alpha * (holdout[i - 1] / sold) + (1 - alpha) * (ellprev + bprev)
      bnew = beta * (ellnew - ellprev) + (1 - beta) * bprev
      snew = gamma * (holdout[i - 1] / ellnew) + (1 - gamma) * sold
      
      seas[prev_pos] = snew
      fc = (ellnew + bnew) * seas[curr_pos]
      
      fc_vec[i] = fc
      yt = holdout[i]
      fcerror = yt - fc
      sse = sse + fcerror^2
      
      ellprev = ellnew
      bprev = bnew
    }
  }
  
  rmse = sqrt(sse / nholdout)
  return(list(fc = fc_vec, rmse = rmse))
}


# 4. AUTOREGRESSIVE MODELS


# AR(p) MODEL
# (train, holdout, arvec, mu)

arp_fc = function(train, holdout, arvec, mu, iprint=F) {
  
  p = length(arvec)
  n = length(train)
  nholdout = length(holdout)
  
  # combine last p of train + holdout
  z = c(train[(n-p+1):n], holdout)
  
  # center
  z = z - mu
  
  sse = 0
  fc_vec = rep(NA, nholdout)
  
  for (i in 1:nholdout) {
    
    # take previous p values
    zprev = z[i:(i+p-1)]
    
    # reverse to match y_t, y_{t-1}, ...
    zprev = rev(zprev)
    
    # forecast
    fc = mu + sum(arvec * zprev)
    
    yt = holdout[i]
    fc_vec[i] = fc
    
    fcerror = yt - fc
    sse = sse + fcerror^2
  }
  
  rmse = sqrt(sse / nholdout)
  
  return(list(fc = fc_vec, rmse = rmse))
}

# MA(q) MODEL
# (train, holdout, mavec, mu)

maq_fc = function(train, holdout, mavec, mu, iprint=F) {
  
  q = length(mavec)
  n = length(train)
  nholdout = length(holdout)
  
  # combine last q of train + holdout
  z = c(train[(n-q+1):n], holdout)
  
  # center
  z = z - mu
  
  sse = 0
  fc_vec = rep(NA, nholdout)
  
  for (i in 1:nholdout) {
    
    # take previous q values
    zprev = z[i:(i+q-1)]
    
    # reverse to match e_t, e_{t-1}, ...
    zprev = rev(zprev)
    
    # forecast
    fc = mu + sum(mavec * zprev)
    
    yt = holdout[i]
    fc_vec[i] = fc
    
    fcerror = yt - fc
    sse = sse + fcerror^2
  }
  
  rmse = sqrt(sse / nholdout)
  
  return(list(fc = fc_vec, rmse = rmse))
}

# General ARIMA (p,d,q) Model

arima_one_step_fc <- function(train, holdout, order = c(0,0,0), 
                              include.mean = TRUE) {
  if (length(train) < 2) stop("train must contain at least 2 observations")
  if (length(holdout) < 1) stop("holdout must contain at least 1 observation")
  if (length(order) != 3) stop("order must be c(p,d,q)")
  
  n_train <- length(train)
  n_holdout <- length(holdout)
  
  # Fit ARIMA model on the training data only
  fit_train <- forecast::Arima(train, order = order, include.mean = include.mean)
  
  # Generate the full series by combining train and holdout (for state updates)
  full_series <- c(train, holdout)
  
  fit_full <- forecast::Arima(full_series, model = fit_train)
  
  # One-step-ahead forecasts for the holdout period.
  # These fitted values at holdout times are conditional on y_{t-1},
  # where y_{t-1} comes from observed data (train + realized holdout history).
  all_fitted <- fitted(fit_full)
  fc_vec <- as.numeric(all_fitted[(n_train + 1):(n_train + n_holdout)])
  
  # Holdout RMSE
  yt <- as.numeric(holdout)
  fcerror <- yt - fc_vec
  rmse <- sqrt(mean(fcerror^2, na.rm = TRUE))
  
  return(list(fc = fc_vec, rmse = rmse))
}


# 5. GENERAL SEASONAL RULES

# GENERAL SEASONAL PERSISTENCE
# (train, holdout, periodicity)

persistbyperiod_fc = function(train, holdout, periodicity, iprint = FALSE) {
  ntrain = length(train)
  nholdout = length(holdout)
  m = periodicity
  
  sse = 0
  fc_vec = numeric(nholdout)
  full_data = c(train, holdout)
  
  for (i in 1:nholdout) {
    yt = holdout[i]
    fc = full_data[(ntrain + i) - m]
    fc_vec[i] = fc
    fcerror = yt - fc
    sse = sse + fcerror^2
  }
  
  rmse = sqrt(sse / nholdout)
  if (iprint) print(cbind(holdout = holdout, fc = fc_vec))
  return(list(rmse = rmse, fc = fc_vec))
}


# Seasonal IID by Period
# (train, holdout, periodicity)

iidbyperiod_fc = function(train, holdout, periodicity, iprint = FALSE) {
  ntrain = length(train)
  nholdout = length(holdout)
  m = periodicity
  
  sse = 0
  fc_vec = numeric(nholdout)
  seasonal_means = numeric(m)
  
  for (j in 1:m) {
    period_indices = seq(j, ntrain, by = m)
    seasonal_means[j] = mean(train[period_indices])
  }
  
  for (i in 1:nholdout) {
    yt = holdout[i]
    period_index = ((ntrain + i - 1) %% m) + 1
    fc = seasonal_means[period_index]
    fc_vec[i] = fc
    fcerror = yt - fc
    sse = sse + fcerror^2
  }
  
  rmse = sqrt(sse / nholdout)
  if (iprint) print(cbind(holdout = holdout, fc = fc_vec))
  return(list(rmse = rmse, fc = fc_vec))
}

sarima_one_step_fc <- function(train, holdout, order = c(0,0,0), 
                               seasonal = c(0,0,0), include.mean = TRUE) {
  if (length(train) < 2) stop("train must contain at least 2 observations")
  if (length(holdout) < 1) stop("holdout must contain at least 1 observation")
  if (length(order) != 3) stop("order must be c(p,d,q)")
  if (length(seasonal) != 3) stop("seasonal must be c(P,D,Q)")
  
  n_train <- length(train)
  n_holdout <- length(holdout)
  
  if (stats::is.ts(train)) {
    freq <- stats::frequency(train)
    st <- stats::start(train)
    full_series <- stats::ts(c(as.numeric(train), as.numeric(holdout)), 
                             start = st, frequency = freq)
  } else {
    full_series <- c(train, holdout)
  }
  
  # Fit SARIMA model on the training data
  # The seasonal argument accepts a numeric vector of length 3: c(P,D,Q)
  fit_train <- forecast::Arima(train, order = order, seasonal = seasonal, 
                               include.mean = include.mean)
  
  # Apply the fitted model to the full series
  fit_full <- forecast::Arima(full_series, model = fit_train)
  
  # Extract one-step-ahead forecasts for the holdout period.
  # These are conditional on observed values up to t-1.
  all_fitted <- stats::fitted(fit_full)
  fc_vec <- as.numeric(all_fitted[(n_train + 1):(n_train + n_holdout)])
  
  # Compute Holdout RMSE
  yt <- as.numeric(holdout)
  fcerror <- yt - fc_vec
  rmse <- sqrt(mean(fcerror^2, na.rm = TRUE))
  
  return(list(fc = fc_vec, rmse = rmse))
}

# 6. ARMAX MODELS

# COSINE-SINE REGRESSION
# (holdout, newx, bvec)

reg_cossin_fc = function(holdout, newx, bvec, iprint = FALSE) {
  fc_vec = bvec[1] + bvec[2] * newx[, "cosine"] + bvec[3] * newx[, "sine"]
  
  rmse = sqrt(mean((holdout - fc_vec)^2))
  if (iprint) print(cbind(holdout = holdout, fc = fc_vec))
  return(list(rmse = rmse, fc = fc_vec))
}


# ARMAX (REG + AR(1) RESIDUAL)
# (train, holdout, newx, newx_prev, bvec, phi)

reg_ar1_cossin_fc = function(train, holdout, newx, newx_prev, bvec, phi,
                             iprint = FALSE) {
  nholdout = length(holdout)
  
  mean_now = bvec[1] + bvec[2] * newx[, "cosine"] + bvec[3] * newx[, "sine"]
  mean_prev = bvec[1] + bvec[2] * newx_prev[, "cosine"] + bvec[3] * newx_prev[, "sine"]
  
  # last training residual for the first holdout forecast
  train_last_mean = mean_prev[1]
  resid_prev = numeric(nholdout)
  resid_prev[1] = train[length(train)] - train_last_mean
  
  if (nholdout >= 2) {
    for (i in 2:nholdout) {
      resid_prev[i] = holdout[i - 1] - mean_prev[i]
    }
  }
  
  fc_vec = mean_now + phi * resid_prev
  rmse = sqrt(mean((holdout - fc_vec)^2))
  
  if (iprint) print(cbind(holdout = holdout, fc = fc_vec))
  return(list(rmse = rmse, fc = fc_vec))
}

# 7. FULL ARIMA / ARMAX SYSTEM

# GENERAL ARIMA FORECAST (ARIMA/ARMA.ARMAX)
# (armaobj, n.ahead)

y_predict = function(armaobj, n.ahead, se.fit = TRUE, alpha = 0.10, iprint = FALSE) {
  predobj = predict(armaobj, n.ahead = n.ahead, se.fit = se.fit)
  
  # If standard errors are requested and available, compute intervals
  if (isTRUE(se.fit) && !is.null(predobj$se)) {
    cv = qnorm(1 - alpha / 2)
    
    lwr = predobj$pred - cv * predobj$se
    upr = predobj$pred + cv * predobj$se
    out = cbind(lwr, predobj$pred, upr)
    colnames(out) = c("lwr", "pred", "upr")
    
    if (iprint) {
      print(out)
      ypred1step = out[1, ]
      cat("1-step ahead forecast interval for y\n")
      print(round(ypred1step, 3))
    } else {
      ypred1step = out[1, ]
    }
    
    return(list(ypred1step = ypred1step, pred = predobj$pred, se = predobj$se))
  }
  
  # When se.fit = FALSE or standard errors are unavailable, skip intervals
  return(list(pred = predobj$pred))
}


# ARIMAX one step-ahead forecast
arimax_one_step_fc <- function(train, holdout, xreg_train, xreg_holdout, 
                               order = c(0,0,0), seasonal = c(0,0,0), 
                               include.mean = TRUE) {
  if (length(train) < 2) stop("train must contain at least 2 observations")
  if (length(holdout) < 1) stop("holdout must contain at least 1 observation")
  if (length(order) != 3) stop("order must be c(p,d,q)")
  if (length(seasonal) != 3) stop("seasonal must be c(P,D,Q)")
  
  # Initialize lengths for train and holdout
  n_train <- length(train)
  n_holdout <- length(holdout)
  
  # State space filtering approach to maintain time series properties and state updates
  xreg_train_mat <- as.matrix(xreg_train)
  xreg_holdout_mat <- as.matrix(xreg_holdout)
  
  if (ncol(xreg_train_mat) != ncol(xreg_holdout_mat)) {
    stop("Error: xreg_train and xreg_holdout must have the same number of columns.")
  }
  if (nrow(xreg_train_mat) != n_train) {
    stop("nrow(xreg_train) must equal length(train)")
  }
  if (nrow(xreg_holdout_mat) != n_holdout) {
    stop("nrow(xreg_holdout) must equal length(holdout)")
  }
  
  # Merge into full matrix
  full_xreg <- rbind(xreg_train_mat, xreg_holdout_mat)
  
  # Keep full series for state updates
  if (stats::is.ts(train)) {
    freq <- stats::frequency(train)
    st <- stats::start(train)
    full_series <- stats::ts(c(as.numeric(train), as.numeric(holdout)), 
                             start = st, frequency = freq)
  } else {
    full_series <- c(train, holdout)
  }
  
  # ARIMAX fit on training data
  fit_train <- forecast::Arima(train, order = order, seasonal = seasonal,
                               xreg = xreg_train_mat, include.mean = include.mean)
  
  # Fit trained coefficients to full series for state updates and forecasting
  fit_full <- forecast::Arima(full_series, model = fit_train, xreg = full_xreg)
  
  # Extract one-step-ahead forecasts for the holdout period.
  # This keeps ARIMA parameters fixed at training estimates while updating
  # state recursively with realized observations.
  all_fitted <- stats::fitted(fit_full)
  fc_vec <- as.numeric(all_fitted[(n_train + 1):(n_train + n_holdout)])
  
  # Calculate holdout RMSE
  yt <- as.numeric(holdout)
  fcerror <- yt - fc_vec
  rmse <- sqrt(mean(fcerror^2, na.rm = TRUE))
  
  return(list(fc = fc_vec, rmse = rmse))
}