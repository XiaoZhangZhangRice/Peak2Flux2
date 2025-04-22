# Step 2.2: ppm2flux - Cumulative emissions and Global Warming Potential (GWP) ####

library(dplyr)
library(zoo)
library(tidyr)
library(stringr)
library(lubridate)

# 1.  Cumulative emissions ####
## 1.1. Determining function ####

flux2acc <- function(data) {

  data_name <- deparse(substitute(data))

  gases <- c("CH4", "N2O", "CO2", "Gas1", "Gas2", "Gas3")

  data <- data %>%
    mutate(Date = parse_date_time(as.character(Date), orders = c("dmy", "mdy", "ymd"))) %>% # deal with date format parsing it into parts
    filter(!is.na(Date))

  corrected_cols <- str_remove(names(data), "_flux_corrected")
  mgm2h_cols <- str_remove(names(data), "_flux_mgm2h")

  present_corrected <- gases[gases %in% corrected_cols]
  present_mgm2h <- gases[gases %in% mgm2h_cols]

  if (length(present_corrected) > 0) {
    suffix <- "_flux_corrected"
    present_gases <- present_corrected
  } else if (length(present_mgm2h) > 0) {
    suffix <- "_flux_mgm2h"
    present_gases <- present_mgm2h
  } else {
    stop("No usable flux columns found: expected *_flux_corrected or *_flux_mgm2h.")
  }

  full_dates_df <- data %>%
    group_by(Plot, Treatment) %>%
    summarise(
      min_date = as.Date(min(Date)),
      max_date = as.Date(max(Date)),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(dates = list(seq(from = min_date, to = max_date, by = "day"))) %>%
    unnest(dates) %>%
    select(Plot, Treatment, Date = dates)

  flux_df <- full_dates_df %>%
    left_join(data, by = c("Plot", "Treatment", "Date")) %>%
    arrange(Plot, Treatment, Date)

  for (gas in present_gases) {

    flux_col <- paste0(gas, suffix)

    emission_col <- paste0(gas, "_emission_kg_ha")

    flux_df <- flux_df %>%
      group_by(Plot, Treatment) %>%
      arrange(Date) %>%
      mutate("{flux_col}" := na.approx(.data[[flux_col]], x = Date, na.rm = FALSE, rule = 2)) %>%
      ungroup()

    flux_df <- flux_df %>%
      mutate("{emission_col}" := .data[[flux_col]] * 0.24)
  }

  kept_fluxes <- paste0(present_gases, suffix)

  kept_emissions <- paste0(present_gases, "_emission_kg_ha")

  cumulative_emissions <- flux_df %>%
    group_by(Plot, Treatment) %>%
    summarise(across(all_of(kept_emissions), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  flux_df <- flux_df %>%
    select(Date, Plot, Treatment, ID, any_of(kept_fluxes), any_of(kept_emissions))

  assign(paste0("daily_", data_name), flux_df, envir = .GlobalEnv)
  assign(paste0("acc_", data_name), cumulative_emissions, envir = .GlobalEnv)

  invisible(NULL)
}

## 1.2. Tests ####
# Using ppm2flux() outputs for tests

# Notes:
# - Try with input data frames with different formats for the Date column. So far the function works for format: "27-Jun-24"

# Test 1: Using as input a dataframe with flux corrections:
flux2acc(flux_df_testA) # ouputs: acc_flux_df_testA (cumulative emissions) and daily_flux_df_testA (linear flux interpolation).

# Test 2: Using as input a dataframe without flux corrections:
flux2acc(flux_df_testM) # ouputs: acc_flux_df_testM (cumulative emissions) and daily_flux_df_testM (linear flux interpolation).

# 2. Global Warming Potential (GWP) ####
## 2.1. Determining function ####

acc2GWP <- function(data,
                     CH4_eq = 27, # default CO2 equivalents according to IPCC, 2021.
                     N2O_eq = 273) {

  data <- data %>%
    mutate(CO2_eq_from_CH4 = CH4_emission_kg_ha * CH4_eq,
           CO2_eq_from_N20 = N2O_emission_kg_ha * N2O_eq,
           Proportion_GWP_CH4 = CO2_eq_from_CH4 / (CO2_eq_from_CH4 + CO2_eq_from_N20),
           Proportion_GWP_N20 = CO2_eq_from_N20 / (CO2_eq_from_CH4 + CO2_eq_from_N20),
           GWP_kg_CO2eq_ha = CO2_eq_from_CH4 + CO2_eq_from_N20)

  return(data)
}

## 2.2. Tests ####
GWP_test1 <- acc2GWP(acc_flux_df_testA)
GWP_test2 <- acc2GWP(acc_flux_df_testM)
