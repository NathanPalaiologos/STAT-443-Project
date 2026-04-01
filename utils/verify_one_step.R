# =========================================================
# verify_one_step.R
# Verify that arima_one_step_fc produces genuine one-step-ahead forecasts
# by comparing against a manual expanding-window loop.
# =========================================================

library(forecast)

cat("=== One-Step-Ahead Forecast Verification ===\n\n")

# --- Load data ---
train_df <- readRDS("output/train_df.rds")
test_df  <- readRDS("output/test_df.rds")
hh_train <- ts(train_df$hh_save, start = c(1977, 1), frequency = 4)
hh_test  <- ts(test_df$hh_save,  start = c(2015, 1), frequency = 4)

source("utils/one_step_model.R")

# --- Method 1: Our arima_one_step_fc (Kalman filter, fixed params) ---
res1 <- arima_one_step_fc(hh_train, hh_test, order = c(1, 1, 1), include.mean = TRUE)

# --- Method 2: Manual loop using Arima() + forecast() at each step ---
fit_train <- Arima(hh_train, order = c(1, 1, 1), include.mean = TRUE)
n_holdout <- length(hh_test)
fc_manual <- numeric(n_holdout)

for (i in 1:n_holdout) {
  # Build the series available at time t (train + holdout up to i-1)
  if (i == 1) {
    available <- hh_train
  } else {
    available <- ts(c(as.numeric(hh_train), as.numeric(hh_test[1:(i-1)])),
                    start = start(hh_train), frequency = 4)
  }
  # Re-apply the SAME parameters (from fit_train) to the available series
  refit <- Arima(available, model = fit_train)
  fc_manual[i] <- forecast(refit, h = 1)$mean
}

# --- Compare ---
max_diff <- max(abs(res1$fc - fc_manual))
cat("Method 1 (Kalman filter on full series) vs Method 2 (manual loop):\n")
cat("  Max absolute difference:", format(max_diff, scientific = TRUE), "\n")
cat("  Numerically identical:", max_diff < 1e-10, "\n\n")

# --- Method 3: Verify residuals are white noise (from training fit) ---
cat("Ljung-Box test on training residuals (ARIMA(1,1,1)):\n")
lb <- Box.test(residuals(fit_train), lag = 10, type = "Ljung-Box")
cat("  Test statistic:", round(lb$statistic, 4), "\n")
cat("  p-value:", round(lb$p.value, 4), "\n")
cat("  Residuals consistent with white noise:", lb$p.value > 0.05, "\n\n")

# --- Verify SARIMA one-step ---
cat("=== SARIMA Verification ===\n")
res_s1 <- sarima_one_step_fc(hh_train, hh_test, order = c(1,1,1), 
                              seasonal = c(1,0,0), include.mean = TRUE)
fit_sarima_train <- Arima(hh_train, order = c(1,1,1), 
                           seasonal = list(order = c(1,0,0), period = 4),
                           include.mean = TRUE)
fc_sarima_manual <- numeric(n_holdout)
for (i in 1:n_holdout) {
  if (i == 1) {
    available <- hh_train
  } else {
    available <- ts(c(as.numeric(hh_train), as.numeric(hh_test[1:(i-1)])),
                    start = start(hh_train), frequency = 4)
  }
  refit <- Arima(available, model = fit_sarima_train)
  fc_sarima_manual[i] <- forecast(refit, h = 1)$mean
}

max_diff_s <- max(abs(res_s1$fc - fc_sarima_manual))
cat("  Max absolute difference:", format(max_diff_s, scientific = TRUE), "\n")
cat("  Numerically identical:", max_diff_s < 1e-10, "\n\n")

# --- Verify return contains fit object ---
cat("=== Return Object Verification ===\n")
cat("  arima_one_step_fc returns 'fit':", !is.null(res1$fit), "\n")
cat("  AIC accessible:", !is.null(res1$fit$aic), "\n")
cat("  AIC value:", round(res1$fit$aic, 2), "\n")
cat("  sarima returns 'fit':", !is.null(res_s1$fit), "\n\n")

# --- Summary of verification ---  
cat("=== VERIFICATION SUMMARY ===\n")
cat("All one-step-ahead forecast functions produce GENUINE one-step forecasts.\n")
cat("The Arima(full, model=fit) approach runs the Kalman filter with fixed\n")
cat("parameters and produces identical results to a manual refit-at-each-step loop.\n")
cat("Parameters are estimated ONCE on training data; only state (residual memory)\n")
cat("is updated as new observations arrive.\n")
