# Load packages
library(tidyverse)
library(dplyr)
library(cansim)
library(readr)
library(lubridate)

#=== Download data from CANSIM ===
macro_metadata <- tribble(
  ~category,               ~table_id,       ~native_freq,
  "net_savings",           "36-10-0111-01", "quarterly",
  "gross_domestic_income", "36-10-0105-01", "quarterly",
  "household_consumption", "36-10-0107-01", "quarterly",
  "unemployment_rate",     "14-10-0287-03", "monthly",
  "cpi",                   "18-10-0004-13", "monthly",
  "market_interest_rates", "10-10-0139-01", "daily"
)

#=== Data wrangling and harmonization ===

fetch_and_harmonize <- function(table_id, category_name, native_freq, start_date = "1977-01-01") {
  
  message(sprintf("Fetching %s (Table ID: %s) - Frequency: %s", category_name, table_id, native_freq))
  
  # Fetch data
  raw_df <- get_cansim(table_id, refresh = TRUE)
  
  # Extract dimensional metadata to safely attach later
  metadata <- raw_df |>
    select(-any_of(c("Date", "REF_DATE", "VALUE", "val_norm", "STATUS", "SYMBOL", "TERMINATED"))) |>
    distinct(VECTOR, .keep_all = TRUE)
  
  # Filter for national-level data and align basic date formatting
  df_clean <- raw_df |>
    mutate(Date = as.Date(Date)) |>
    filter(Date >= as.Date(start_date), GEO == "Canada") |>
    mutate(quarterly_date = floor_date(Date, "quarter")) # Establishes the Q1/Q2/Q3/Q4 boundary
  
  # 3. Frequency-specific handling logic
  if (native_freq == "quarterly") {
    # Natively quarterly: No aggregation required. Just pass the value through.
    harmonized_df <- df_clean |>
      select(VECTOR, quarterly_date, VALUE)
    
  } else if (native_freq == "monthly") {
    # Natively monthly: CPI & Unemployment. 
    # Macro standard is to take the 3-month average for the quarter. 
    # (If you prefer end-of-period, you would use `filter(Date == max(Date))` instead)
    harmonized_df <- df_clean |>
      group_by(VECTOR, quarterly_date) |>
      filter(Date == max(Date))
    
  } else if (native_freq == "daily") {
    # Natively daily: Interest rates.
    # Average all daily trading observations over the entire quarter.
    harmonized_df <- df_clean |>
      group_by(VECTOR, quarterly_date) |>
      summarise(VALUE = mean(VALUE, na.rm = TRUE), .groups = "drop")
  }
  
  # 4. Re-attach metadata and add category tags
  harmonized_df <- harmonized_df |>
    left_join(metadata, by = "VECTOR") |>
    mutate(category = category_name)
  
  return(harmonized_df)
}

# Stack all the harmonized data into a single long-format panel dataset
macro_panel_long <- pmap_dfr(
  list(macro_metadata$table_id, macro_metadata$category, macro_metadata$native_freq),
  ~fetch_and_harmonize(..1, ..2, ..3, start_date = "1977-01-01")
) |>
  select(quarterly_date, VALUE, category)|>
  rename(
    value = VALUE
  )

glimpse(macro_panel_long)

macro_panel_wide <- macro_panel_long |>
  pivot_wider(names_from = category, values_from = value)|>
  mutate(quarterly_date = as.Date(quarterly_date)) |>
  

glimpse(macro_panel_wide)

#=== Save harmonized data ===
write.csv(macro_panel_wide, "../data/macro_panel_wide_raw.csv", row.names = FALSE)

