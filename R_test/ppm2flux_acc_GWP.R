# Step 2.2: ppm2flux - Cumulative emissions ####

# Load required packages:
library(dplyr)
library(zoo)
library(tidyr)

# 1. Determining function ####

flux2acc <- function(data) {

  Gases <- c("CH4", "N2O", "CO2", "Gas1", "Gas2", "Gas3")

  data_name <- deparse(substitute(data))
  output_name <- paste0("acc_", data_name) # to call resulting cumulative emissions data frame according to input data frame

  flux_df <- data %>%
    mutate(Date =as.Date(Date, format = "%d-%b-%y"))

  start_date <- min(flux_df$Date)
  end_date <- max(flux_df$Date)
  date_seq <- seq.Date(start_date, end_date, by = "day")

  flux_df <- flux_df %>% # fill missing dates in between sampling events
    group_by(Plot) %>%
    complete(Date = date_seq) %>%
    ungroup()

    ## 1.1. Flux interpolation ####
  kept_gases <- c()
  kept_emissions <- c()

  for(gas in Gases) {
    flux_col <- paste0(gas, "_flux_corrected")
    emission_col <- paste0(gas, "_emission_kg_ha")

      if (flux_col %in% names(flux_df)) {
        flux_df <- flux_df %>%
          mutate("{flux_col}" := na.approx(.data[[flux_col]], rule = 2)) %>% # linear extrapolation for fluxes in between sampling events
          mutate("{emission_col}" := .data[[flux_col]] * 0.24) %>% # from mg m-2 h-1 to kg ha-1 day-1
          ungroup()

        kept_gases <- c(kept_gases, flux_col)
        kept_emissions <- c(kept_emissions, emission_col)

  } # closes if() for flux_col %in% names(flux_df)
  } # closes for(gas in Gases)

  cumulative_emissions <- flux_df %>%
    group_by(Plot) %>%
    summarise(across(all_of(kept_emissions), ~ sum(.x, na.rm = TRUE))) %>%
    ungroup()

  assign(output_name, cumulative_emissions, envir = .GlobalEnv)

  flux_df <- flux_df %>%
    select(Date, Plot, all_of(kept_gases), all_of(kept_emissions))

  return(flux_df)

  } # closes f(x)

# 2. Tests ####
# Using ppm2flux() outputs for tests

# Notes:
# - Try with input data frames with different formats for the Date column. So far "%d-%b-%y" works for format: "27-Jun-24"

input_acc_test1 <- flux2acc(flux_df_testA)
