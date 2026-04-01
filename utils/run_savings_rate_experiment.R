# =========================================================
# run_savings_rate_experiment.R
# Test whether forecasting hh_rate beats benchmarks
# =========================================================

library(forecast)
library(tseries)
suppressPackageStartupMessages(library(tidyverse))

cat("=== Savings Rate Forecasting Experiment ===\n\n")

# --- Load data ---
train_df <- readRDS("output/train_df.rds")
test_df  <- readRDS("output/test_df.rds")
source("utils/one_step_model.R")

hh_rate_train <- ts(train_df$hh_rate, start = c(1977, 1), frequency = 4)
hh_rate_test  <- ts(test_df$hh_rate,  start = c(2015, 1), frequency = 4)

cat("Saving rate range (train):", round(range(hh_rate_train), 4), "\n")
cat("Saving rate range (test):", round(range(hh_rate_test), 4), "\n\n")

# --- Stationarity ---
cat("=== Stationarity Tests ===\n")
adf_level <- adf.test(na.omit(as.numeric(hh_rate_train)))
cat("ADF on hh_rate (level): p =", round(adf_level$p.value, 4), "\n")
adf_diff <- adf.test(na.omit(diff(as.numeric(hh_rate_train))))
cat("ADF on diff(hh_rate):   p =", round(adf_diff$p.value, 4), "\n\n")

# --- Benchmarks ---
cat("=== Benchmarks ===\n")
res_avg     <- iid_fc(hh_rate_train, hh_rate_test)
res_persist <- persist_fc(hh_rate_train, hh_rate_test)
res_period  <- iidbyperiod_fc(hh_rate_train, hh_rate_test, periodicity = 4)

# SES
hw_ses <- HoltWinters(hh_rate_train, beta = FALSE, gamma = FALSE)
res_ses <- esm_fc(hh_rate_train, hh_rate_test,
                  alpha = hw_ses$alpha,
                  level = as.numeric(hw_ses$coefficients[1]))

# Holt
hw_holt <- HoltWinters(hh_rate_train, gamma = FALSE)
res_holt <- lholt_fc(hh_rate_train, hh_rate_test,
                     alpha = hw_holt$alpha,
                     beta  = hw_holt$beta,
                     level = as.numeric(hw_holt$coefficients[1]),
                     slope = as.numeric(hw_holt$coefficients[2]))

benchmark_tbl <- tibble(
  Model = c("Average", "Persistence", "Avg by Period", "SES", "Holt Linear"),
  RMSE = round(c(res_avg$rmse, res_persist$rmse, res_period$rmse,
                  res_ses$rmse, res_holt$rmse), 6)
)
print(benchmark_tbl)

# --- ARIMA models ---
cat("\n=== auto.arima on hh_rate ===\n")
auto_rate <- auto.arima(hh_rate_train, seasonal = TRUE,
                        stepwise = FALSE, approximation = FALSE)
summary(auto_rate)

cat("\n=== ARIMA / SARIMA Specifications ===\n")
specs <- list(
  list(name = "ARIMA(0,1,1)", order = c(0,1,1)),
  list(name = "ARIMA(1,1,1)", order = c(1,1,1)),
  list(name = "ARIMA(1,1,0)", order = c(1,1,0)),
  list(name = "ARIMA(2,1,0)", order = c(2,1,0))
)

arima_results <- map_dfr(specs, function(s) {
  res <- arima_one_step_fc(hh_rate_train, hh_rate_test,
                            order = s$order, include.mean = TRUE)
  tibble(Model = s$name, AIC = round(res$fit$aic, 2), RMSE = round(res$rmse, 6))
})

# SARIMA specs
sarima_specs <- list(
  list(name = "SARIMA(0,1,1)(0,0,2)[4]", order = c(0,1,1), seasonal = c(0,0,2)),
  list(name = "SARIMA(1,1,1)(1,0,0)[4]", order = c(1,1,1), seasonal = c(1,0,0)),
  list(name = "SARIMA(0,1,1)(1,0,0)[4]", order = c(0,1,1), seasonal = c(1,0,0)),
  list(name = "SARIMA(1,1,1)(0,0,1)[4]", order = c(1,1,1), seasonal = c(0,0,1))
)

sarima_results <- map_dfr(sarima_specs, function(s) {
  res <- sarima_one_step_fc(hh_rate_train, hh_rate_test,
                             order = s$order, seasonal = s$seasonal,
                             include.mean = TRUE)
  tibble(Model = s$name, AIC = round(res$fit$aic, 2), RMSE = round(res$rmse, 6))
})

all_arima <- bind_rows(arima_results, sarima_results) |> arrange(RMSE)
print(all_arima)

# --- Combined RMSE ranking ---
cat("\n=== COMBINED RMSE RANKING ===\n")
best_arima_row <- all_arima[1, ]
combined <- bind_rows(benchmark_tbl, tibble(Model = best_arima_row$Model,
                                             RMSE = best_arima_row$RMSE)) |>
  arrange(RMSE)
print(combined)

persist_rmse <- benchmark_tbl$RMSE[benchmark_tbl$Model == "Persistence"]
best_rmse <- combined$RMSE[1]
best_model <- combined$Model[1]

cat("\nBest model:", best_model, "| RMSE:", best_rmse, "\n")
cat("Persistence RMSE:", persist_rmse, "\n")
cat("Beats Persistence?", ifelse(best_rmse < persist_rmse, "YES", "NO"), "\n")

# --- Also test expanding window on best ---
cat("\n=== Expanding Window on Best ARIMA Spec ===\n")
# Parse best ARIMA order
best_order <- specs[[which(map_chr(specs, "name") == best_arima_row$Model)]]$order
ew_res <- arima_expanding_window_fc(hh_rate_train, hh_rate_test,
                                     order = best_order, include.mean = TRUE)
cat("Expanding window RMSE:", round(ew_res$rmse, 6), "\n")
cat("Fixed-param RMSE:     ", best_arima_row$RMSE, "\n")

cat("\n=== CONCLUSION ===\n")
if (best_rmse < persist_rmse) {
  diff_pct <- (persist_rmse - best_rmse) / persist_rmse * 100
  cat("The saving rate CAN be forecast better than persistence by", 
      round(diff_pct, 2), "% RMSE reduction.\n")
} else {
  cat("The saving rate CANNOT beat the persistence benchmark.\n")
  cat("Like the level series, hh_rate is effectively a random walk.\n")
}
