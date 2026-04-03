# =========================================================
# expanding_window.R
# Expanding-window one-step-ahead forecast functions
# Re-estimates the model at each holdout step using an
# expanding training window.
# =========================================================

# ARIMA Expanding Window
# At each holdout time t, fits ARIMA on observations 1...(n_train + t - 1),
# then forecasts one step ahead.

arima_expanding_window_fc <- function(train, holdout, order = c(0,0,0),
                                      include.mean = TRUE) {
  n_train <- length(train)
  n_holdout <- length(holdout)
  full <- c(as.numeric(train), as.numeric(holdout))
  
  fc_vec <- numeric(n_holdout)
  
  for (t in 1:n_holdout) {
    window_end <- n_train + t - 1
    train_window <- full[1:window_end]
    if (stats::is.ts(train)) {
      train_window <- stats::ts(train_window,
                                start = stats::start(train),
                                frequency = stats::frequency(train))
    }
    fit_t <- forecast::Arima(train_window, order = order, include.mean = include.mean)
    fc_vec[t] <- forecast::forecast(fit_t, h = 1)$mean
  }
  
  yt <- as.numeric(holdout)
  fcerror <- yt - fc_vec
  rmse <- sqrt(mean(fcerror^2))
  
  return(list(fc = fc_vec, rmse = rmse))
}

# SARIMA Expanding Window

sarima_expanding_window_fc <- function(train, holdout, order = c(0,0,0),
                                       seasonal = c(0,0,0),
                                       include.mean = TRUE) {
  n_train <- length(train)
  n_holdout <- length(holdout)
  full <- c(as.numeric(train), as.numeric(holdout))
  
  fc_vec <- numeric(n_holdout)
  
  for (t in 1:n_holdout) {
    window_end <- n_train + t - 1
    train_window <- full[1:window_end]
    if (stats::is.ts(train)) {
      train_window <- stats::ts(train_window,
                                start = stats::start(train),
                                frequency = stats::frequency(train))
    }
    fit_t <- forecast::Arima(train_window, order = order, seasonal = seasonal,
                             include.mean = include.mean)
    fc_vec[t] <- forecast::forecast(fit_t, h = 1)$mean
  }
  
  yt <- as.numeric(holdout)
  fcerror <- yt - fc_vec
  rmse <- sqrt(mean(fcerror^2))
  
  return(list(fc = fc_vec, rmse = rmse))
}

# ARIMAX Expanding Window

arimax_expanding_window_fc <- function(train, holdout, xreg_train, xreg_holdout,
                                       order = c(0,0,0), seasonal = c(0,0,0),
                                       include.mean = TRUE) {
  n_train <- length(train)
  n_holdout <- length(holdout)
  full <- c(as.numeric(train), as.numeric(holdout))
  
  xreg_train_mat <- as.matrix(xreg_train)
  xreg_holdout_mat <- as.matrix(xreg_holdout)
  full_xreg <- rbind(xreg_train_mat, xreg_holdout_mat)
  
  fc_vec <- numeric(n_holdout)
  
  for (t in 1:n_holdout) {
    window_end <- n_train + t - 1
    train_window <- full[1:window_end]
    xreg_window <- full_xreg[1:window_end, , drop = FALSE]
    
    if (stats::is.ts(train)) {
      train_window <- stats::ts(train_window,
                                start = stats::start(train),
                                frequency = stats::frequency(train))
    }
    fit_t <- forecast::Arima(train_window, order = order, seasonal = seasonal,
                             xreg = xreg_window, include.mean = include.mean)
    newxreg <- full_xreg[window_end + 1, , drop = FALSE]
    fc_vec[t] <- forecast::forecast(fit_t, h = 1, xreg = newxreg)$mean
  }
  
  yt <- as.numeric(holdout)
  fcerror <- yt - fc_vec
  rmse <- sqrt(mean(fcerror^2))
  
  return(list(fc = fc_vec, rmse = rmse))
}
