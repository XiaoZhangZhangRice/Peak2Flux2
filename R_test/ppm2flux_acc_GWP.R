# Step 2.2: ppm2flux - Cumulative emissions ####

# Load required packages:
library(dplyr)
library(zoo)
library(tidyr)

# 1. Determining function ####

flux2acc <- function(data) {

  Gases <- c("CH4", "N2O", "CO2", "Gas1", "Gas2", "Gas3")

  flux_df <- data %>%
    mutate(Date =as.Date(Date, format = "%d-%b-%y"))

  start_date <- min(flux_df$Date)
  end_date <- max(flux_df$Date)
  date_seq <- seq.Date(start_date, end_date, by = "day")

    ## 1.1. Flux interpolation ####

  for(gas in Gases) {
    flux_col <- paste0(gas, "_flux_corrected")

      flux_df <- flux_df %>%
          group_by(Plot) %>%
          complete(Date = date_seq) %>%
          ungroup()

      if (flux_col %in% names(flux_df)) {
        flux_df <- flux_df %>%
          mutate("{flux_col}" := na.approx(.data[[flux_col]], rule = 2)) %>%
          ungroup()

    return(flux_df)

  } # closes if() for flux_col %in% names(flux_df)
  } # closes for(gas in Gases)
  } # closes f(x)

# 2. Tests ####
# Using ppm2flux() outputs for tests

# Notes:
# - Try with input data frames with different formats for the Date column. So far "%d-%b-%y" works for format: "27-Jun-24"

input_acc_test1 <- flux2acc(flux_df_testA)
