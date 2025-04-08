# Step 2.2: ppm2flux - Cumulative emissions and Global Warming Potential (GWP) ####

# 1.  Cumulative emissions ####
## 1.1. Determining function ####

flux2acc <- function(data) {

  library(dplyr)
  library(zoo)
  library(tidyr)
  library(stringr)
  library(lubridate)

  data_name <- deparse(substitute(data))

  gases <- c("CH4", "N2O", "CO2", "Gas1", "Gas2", "Gas3")

  data <- data %>%
    mutate(Date = parse_date_time(as.character(Date), orders = c("dmy", "mdy", "ymd"))) %>% # deal with date format parsing it into parts
    filter(!is.na(Date))

  present_gases <- gases[gases %in% str_remove(names(data), "_flux_corrected")]

  full_dates_df <- data %>%
    group_by(Plot) %>%
    summarise(
      min_date = as.Date(min(Date)),
      max_date = as.Date(max(Date)),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(dates = list(seq(from = min_date, to = max_date, by = "day"))) %>%
    unnest(dates) %>%
    select(Plot, Date = dates)

  flux_df <- full_dates_df %>%
    left_join(data, by = c("Plot", "Date")) %>%
    arrange(Plot, Date)

  for (gas in present_gases) {
    flux_col <- paste0(gas, "_flux_corrected")
    emission_col <- paste0(gas, "_emission_kg_ha")

    flux_df <- flux_df %>%
      group_by(Plot) %>%
      arrange(Date) %>%
      mutate("{flux_col}" := na.approx(.data[[flux_col]], x = Date, na.rm = FALSE, rule = 2)) %>%
      ungroup()

    flux_df <- flux_df %>%
      mutate("{emission_col}" := .data[[flux_col]] * 0.24)
  }

  kept_emissions <- paste0(present_gases, "_emission_kg_ha")

  cumulative_emissions <- flux_df %>%
    group_by(Plot) %>%
    summarise(across(all_of(kept_emissions), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  flux_df <- flux_df %>%
    select(Date, Plot, ID, any_of(paste0(present_gases, "_flux_corrected")), any_of(kept_emissions))

  assign(paste0("daily_", data_name), flux_df, envir = .GlobalEnv)
  assign(paste0("acc_", data_name), cumulative_emissions, envir = .GlobalEnv)

  invisible(NULL)
}

## 1.3. Tests ####
# Using ppm2flux() outputs for tests

# Notes:
# - Try with input data frames with different formats for the Date column. So far the function works for format: "27-Jun-24"

input_acc_test1 <- flux2acc(flux_df_testA) # ouputs: acc_flux_df_testA (cumulative emissions) and daily_flux_df_testA (linear flux interpolation).
input_acc_test2 <- flux2acc(flux_df_testM) # ouputs: acc_flux_df_testM (cumulative emissions) and daily_flux_df_testM (linear flux interpolation).

# 2. Global Warming Potential (GWP) ####

flux2GWP <- function(data,
                     CH4_eq = 27, # default CO2 equivalents according to IPCC, 2021.
                     N2O_eq = 273) {




}
