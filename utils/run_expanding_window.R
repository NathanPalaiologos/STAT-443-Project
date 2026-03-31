# Quick expanding window comparison for hh_save level
library(forecast)
suppressPackageStartupMessages(library(tidyverse))

train_df <- readRDS("output/train_df.rds")
test_df  <- readRDS("output/test_df.rds")
source("utils/one_step_model.R")

hh_train <- ts(train_df$hh_save, start = c(1977, 1), frequency = 4)
hh_test  <- ts(test_df$hh_save,  start = c(2015, 1), frequency = 4)

cat("=== Expanding Window vs Fixed on hh_save ===\n\n")

# Persistence (same either way)
persist <- persist_fc(hh_train, hh_test)
cat("Persistence RMSE:", round(persist$rmse, 4), "\n\n")

# ARIMA(1,1,1) fixed
fixed_111 <- arima_one_step_fc(hh_train, hh_test, order = c(1,1,1), include.mean = TRUE)
cat("ARIMA(1,1,1) fixed-param RMSE:", round(fixed_111$rmse, 4), "\n")

# ARIMA(1,1,1) expanding
ew_111 <- arima_expanding_window_fc(hh_train, hh_test, order = c(1,1,1), include.mean = TRUE)
cat("ARIMA(1,1,1) expanding RMSE:  ", round(ew_111$rmse, 4), "\n\n")

# ARIMA(0,1,1) fixed
fixed_011 <- arima_one_step_fc(hh_train, hh_test, order = c(0,1,1), include.mean = TRUE)
cat("ARIMA(0,1,1) fixed-param RMSE:", round(fixed_011$rmse, 4), "\n")

# ARIMA(0,1,1) expanding
ew_011 <- arima_expanding_window_fc(hh_train, hh_test, order = c(0,1,1), include.mean = TRUE)
cat("ARIMA(0,1,1) expanding RMSE:  ", round(ew_011$rmse, 4), "\n\n")

# SARIMA(1,1,1)(1,0,0)[4] 
fixed_s <- sarima_one_step_fc(hh_train, hh_test, order = c(1,1,1), 
                               seasonal = c(1,0,0), include.mean = TRUE)
cat("SARIMA(1,1,1)(1,0,0)[4] fixed RMSE:", round(fixed_s$rmse, 4), "\n")

ew_s <- sarima_expanding_window_fc(hh_train, hh_test, order = c(1,1,1),
                                    seasonal = c(1,0,0), include.mean = TRUE)
cat("SARIMA(1,1,1)(1,0,0)[4] expanding RMSE:", round(ew_s$rmse, 4), "\n")
