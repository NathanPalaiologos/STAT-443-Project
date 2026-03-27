# data_pipeline/helpers_forecast.R

rmse_vec <- function(actual, pred) {
  sqrt(mean((actual - pred)^2, na.rm = TRUE))
}

# Average forecast rule
avg_fc <- function(train, holdout) {
  nh <- length(holdout)
  fc_vec <- rep(mean(train, na.rm = TRUE), nh)
  list(rmse = rmse_vec(holdout, fc_vec), fc = fc_vec)
}

# Persistence forecast rule
persist_fc <- function(train, holdout) {
  ntrain <- length(train)
  nh <- length(holdout)
  full_data <- c(train, holdout)
  fc_vec <- numeric(nh)
  
  for (i in seq_len(nh)) {
    fc_vec[i] <- full_data[ntrain + i - 1]
  }
  
  list(rmse = rmse_vec(holdout, fc_vec), fc = fc_vec)
}

# Seasonal persistence by period (professor-style)
persistbyperiod_fc <- function(train, holdout, periodicity) {
  ntrain <- length(train)
  nholdout <- length(holdout)
  m <- periodicity
  sse <- 0
  fc_vec <- numeric(nholdout)
  full_data <- c(train, holdout)
  
  for (i in seq_len(nholdout)) {
    yt <- holdout[i]
    fc <- full_data[(ntrain + i) - m]
    fc_vec[i] <- fc
    sse <- sse + (yt - fc)^2
  }
  
  list(rmse = sqrt(sse / nholdout), fc = fc_vec)
}

# Seasonal mean by period (professor-style)
iidbyperiod_fc <- function(train, holdout, periodicity) {
  ntrain <- length(train)
  nholdout <- length(holdout)
  m <- periodicity
  sse <- 0
  fc_vec <- numeric(nholdout)
  seasonal_means <- numeric(m)
  
  for (j in seq_len(m)) {
    period_indices <- seq(j, ntrain, by = m)
    seasonal_means[j] <- mean(train[period_indices], na.rm = TRUE)
  }
  
  for (i in seq_len(nholdout)) {
    yt <- holdout[i]
    period_index <- ((ntrain + i - 1) %% m) + 1
    fc <- seasonal_means[period_index]
    fc_vec[i] <- fc
    sse <- sse + (yt - fc)^2
  }
  
  list(rmse = sqrt(sse / nholdout), fc = fc_vec)
}

# Winters additive seasonal exponential smoothing (professor-style)
aseason_fc <- function(train, holdout, alpha, beta, gamma, level, slope, season) {
  nh <- length(holdout)
  m <- length(season)
  sse <- 0
  fc_vec <- numeric(nh)
  
  L <- level
  b <- slope
  s <- season
  
  for (i in seq_len(nh)) {
    idx <- ((i - 1) %% m) + 1
    fc <- L + b + s[idx]
    fc_vec[i] <- fc
    
    yt <- holdout[i]
    fcerror <- yt - fc
    sse <- sse + fcerror^2
    
    L_old <- L
    L <- alpha * (yt - s[idx]) + (1 - alpha) * (L_old + b)
    b <- beta * (L - L_old) + (1 - beta) * b
    s[idx] <- gamma * (yt - L_old - b) + (1 - gamma) * s[idx]
  }
  
  list(rmse = sqrt(sse / nh), fc = fc_vec)
}

# 1-step-ahead holdout RMSE for ARIMA (professor-style)
arima_fc <- function(tsdata,
                     ntrain,
                     order,
                     seasonal = list(order = c(0, 0, 0), period = 1),
                     method,
                     traincoef,
                     include.mean,
                     iprint = FALSE) {
  
  obj <- arima(
    tsdata,
    order = order,
    seasonal = seasonal,
    init = traincoef,
    fixed = traincoef,
    method = method,
    include.mean = include.mean,
    optim.control = list(maxit = 0)
  )
  
  fc <- tsdata - obj$residuals
  ntotal <- length(tsdata)
  holdout_fc <- fc[(ntrain + 1):ntotal]
  holdout <- tsdata[(ntrain + 1):ntotal]
  
  if (iprint) print(cbind(holdout, holdout_fc))
  
  rmse <- sqrt(mean((holdout - holdout_fc)^2))
  list(rmse = rmse, fc = holdout_fc)
}

# 1-step-ahead holdout RMSE for ARMAX / ARX (professor-style)
armax_fc <- function(tsdata,
                     ntrain,
                     order,
                     method,
                     traincoef,
                     include.mean,
                     xreg,
                     iprint = FALSE) {
  
  obj <- arima(
    tsdata,
    order = order,
    init = traincoef,
    fixed = traincoef,
    method = method,
    include.mean = include.mean,
    xreg = xreg,
    optim.control = list(maxit = 0)
  )
  
  fc <- tsdata - obj$residuals
  ntotal <- length(tsdata)
  holdout_fc <- fc[(ntrain + 1):ntotal]
  holdout <- tsdata[(ntrain + 1):ntotal]
  
  if (iprint) print(cbind(holdout, holdout_fc))
  
  rmse <- sqrt(mean((holdout - holdout_fc)^2))
  list(rmse = rmse, fc = holdout_fc)
}