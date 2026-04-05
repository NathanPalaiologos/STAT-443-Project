# =========================================================
# expanding_window.R
# Rolling / expanding-window one-step-ahead forecast functions
# Re-estimates the model at each holdout step.
# When window_size is specified, uses a fixed rolling window;
# otherwise the training window expands from the initial split.
# =========================================================

# ARIMA Rolling / Expanding Window

arima_expanding_window_fc <- function(train, holdout, order = c(0,0,0),
                                      include.mean = TRUE,
                                      window_size = NULL) {
  n_train <- length(train)
  n_holdout <- length(holdout)
  full <- c(as.numeric(train), as.numeric(holdout))
  
  fc_vec <- numeric(n_holdout)
  
  for (t in 1:n_holdout) {
    window_end <- n_train + t - 1
    window_start <- if (!is.null(window_size)) max(1, window_end - window_size + 1) else 1
    train_window <- full[window_start:window_end]
    if (stats::is.ts(train)) {
      freq <- stats::frequency(train)
      orig_start <- stats::start(train)
      offset <- window_start - 1
      new_year <- orig_start[1] + (orig_start[2] - 1 + offset) %/% freq
      new_qtr  <- (orig_start[2] - 1 + offset) %% freq + 1
      train_window <- stats::ts(train_window, start = c(new_year, new_qtr),
                                frequency = freq)
    }
    fit_t <- forecast::Arima(train_window, order = order, include.mean = include.mean)
    fc_vec[t] <- forecast::forecast(fit_t, h = 1)$mean
  }
  
  yt <- as.numeric(holdout)
  fcerror <- yt - fc_vec
  rmse <- sqrt(mean(fcerror^2))
  
  return(list(fc = fc_vec, rmse = rmse))
}

# SARIMA Rolling / Expanding Window

sarima_expanding_window_fc <- function(train, holdout, order = c(0,0,0),
                                       seasonal = c(0,0,0),
                                       include.mean = TRUE,
                                       window_size = NULL) {
  n_train <- length(train)
  n_holdout <- length(holdout)
  full <- c(as.numeric(train), as.numeric(holdout))
  
  fc_vec <- numeric(n_holdout)
  
  for (t in 1:n_holdout) {
    window_end <- n_train + t - 1
    window_start <- if (!is.null(window_size)) max(1, window_end - window_size + 1) else 1
    train_window <- full[window_start:window_end]
    if (stats::is.ts(train)) {
      freq <- stats::frequency(train)
      orig_start <- stats::start(train)
      offset <- window_start - 1
      new_year <- orig_start[1] + (orig_start[2] - 1 + offset) %/% freq
      new_qtr  <- (orig_start[2] - 1 + offset) %% freq + 1
      train_window <- stats::ts(train_window, start = c(new_year, new_qtr),
                                frequency = freq)
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

# ARIMAX Rolling / Expanding Window

arimax_expanding_window_fc <- function(train, holdout, xreg_train, xreg_holdout,
                                       order = c(0,0,0), seasonal = c(0,0,0),
                                       include.mean = TRUE,
                                       window_size = NULL) {
  n_train <- length(train)
  n_holdout <- length(holdout)
  full <- c(as.numeric(train), as.numeric(holdout))
  
  xreg_train_mat <- as.matrix(xreg_train)
  xreg_holdout_mat <- as.matrix(xreg_holdout)
  full_xreg <- rbind(xreg_train_mat, xreg_holdout_mat)
  
  fc_vec <- numeric(n_holdout)
  
  for (t in 1:n_holdout) {
    window_end <- n_train + t - 1
    window_start <- if (!is.null(window_size)) max(1, window_end - window_size + 1) else 1
    train_window <- full[window_start:window_end]
    xreg_window <- full_xreg[window_start:window_end, , drop = FALSE]
    
    if (stats::is.ts(train)) {
      freq <- stats::frequency(train)
      orig_start <- stats::start(train)
      offset <- window_start - 1
      new_year <- orig_start[1] + (orig_start[2] - 1 + offset) %/% freq
      new_qtr  <- (orig_start[2] - 1 + offset) %% freq + 1
      train_window <- stats::ts(train_window, start = c(new_year, new_qtr),
                                frequency = freq)
    }
    fit_t <- tryCatch(
      forecast::Arima(train_window, order = order, seasonal = seasonal,
                      xreg = xreg_window, include.mean = include.mean),
      error = function(e) {
        forecast::Arima(train_window, order = order, seasonal = seasonal,
                        xreg = xreg_window, include.mean = include.mean,
                        method = "ML")
      }
    )
    newxreg <- full_xreg[window_end + 1, , drop = FALSE]
    fc_vec[t] <- forecast::forecast(fit_t, h = 1, xreg = newxreg)$mean
  }
  
  yt <- as.numeric(holdout)
  fcerror <- yt - fc_vec
  rmse <- sqrt(mean(fcerror^2))
  
  return(list(fc = fc_vec, rmse = rmse))
}
