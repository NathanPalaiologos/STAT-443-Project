# Load packages
library(tidyverse)
library(cansim)
library(lubridate)

#=== 1. Define CANSIM Metadata & Target Headline Vectors ===
# Swapped the standard CPI table for the Core Inflation table (18-10-0256-01)
macro_metadata <- tribble(
  ~category,               ~table_id,       ~native_freq, ~headline_vector,
  "net_savings", "36-10-0111-01", "quarterly", "v62305783",  # Net saving
  "gross_domestic_income", "36-10-0105-01", "quarterly", "v62305752",  # Gross domestic income
  "household_consumption", "36-10-0107-01", "quarterly", "v62305732",  # Household final consumption expenditure
  "unemployment_rate", "14-10-0287-03", "monthly", "v2062815",   # Unemployment rate, 15+, Both sexes
  "cpi_trim_yoy", "18-10-0256-01", "monthly", "v108785715", # CPI-trim (Y/Y % change)
  "market_interest_rates", "10-10-0139-01", "daily", "v122530"     # Bank rate
)

#=== 2. Data wrangling and harmonization ===
fetch_and_harmonize <- function(table_id, category_name, native_freq, target_vector, start_date = "1977-01-01") {
  
  message(sprintf("Fetching %s (Table ID: %s) - Frequency: %s", category_name, table_id, native_freq))
  
  # Fetch data
  raw_df <- get_cansim(table_id, refresh = TRUE)
  
  # Filter for the specific headline vector right away to keep data light
  df_clean <- raw_df |>
    filter(VECTOR == target_vector) |>
    mutate(Date = as.Date(Date)) |>
    filter(Date >= as.Date(start_date)) |>
    mutate(quarterly_date = floor_date(Date, "quarter"))
  
  # 3. Frequency-specific handling logic
  if (native_freq == "quarterly") {
    
    harmonized_df <- df_clean |>
      select(VECTOR, quarterly_date, VALUE)
    
  } else if (native_freq == "monthly") {
    
    # Natively monthly: Extract the end-of-quarter value
    harmonized_df <- df_clean |>
      group_by(VECTOR, quarterly_date) |>
      filter(Date == max(Date)) |> 
      ungroup() |> 
      select(VECTOR, quarterly_date, VALUE)
    
  } else if (native_freq == "daily") {
    
    # Natively daily: Average over the quarter
    harmonized_df <- df_clean |>
      group_by(VECTOR, quarterly_date) |>
      summarise(VALUE = mean(VALUE, na.rm = TRUE), .groups = "drop")
    
  }
  
  # 4. Add category tags
  harmonized_df <- harmonized_df |>
    mutate(category = category_name)
  
  return(harmonized_df)
}

#=== 3. Stack and Pivot Data ===

# Stack all the harmonized data into a single long-format panel dataset
macro_panel_long <- pmap_dfr(
  list(macro_metadata$table_id, macro_metadata$category, macro_metadata$native_freq, macro_metadata$headline_vector),
  ~fetch_and_harmonize(..1, ..2, ..3, ..4, start_date = "1977-01-01")
) |>
  select(quarterly_date, VALUE, category) |>
  rename(value = VALUE)

# Pivot to wide format
macro_panel_wide <- macro_panel_long |>
  pivot_wider(
    names_from = category, 
    values_from = value
  ) |>
  mutate(quarterly_date = as.Date(quarterly_date)) |>
  arrange(quarterly_date)

glimpse(macro_panel_wide)

#=== 4. Save harmonized data ===
write.csv(macro_panel_wide, "macro_panel_wide_raw.csv", row.names = FALSE)