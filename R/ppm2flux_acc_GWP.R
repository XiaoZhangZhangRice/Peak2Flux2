#' Converts flux into cumulative emission.
#'
#' Interpolates linearly throughout a time series of flux data to calculate daily flux and cumulative emissions.
#'
#' @param data A data frame containing gas fluxes for a series of sampling events.
#' @return Data frames containing calculated daily fluxes in \eqn{mg \cdot m^{-2} \cdot h^{-1}} (name: "daily_data") and cumulative emissions \eqn{kg \cdot ha^{-1}} (name: "acc_data").
#' @importFrom dplyr %>%
#' @export
#' @examples
#' # Test 1: Using as input a data rame with flux corrections:
#' flux2acc(flux_df_testA) # ouputs: acc_flux_df_testA (cumulative emissions) and daily_flux_df_testA (linear flux interpolation).
#'
#' # Test 2: Using as input a data frame without flux corrections:
#' flux2acc(flux_df_testM) # ouputs: acc_flux_df_testM (cumulative emissions) and daily_flux_df_testM (linear flux interpolation).

flux2acc <- function(data) {

  data_name <- deparse(substitute(data))

  gases <- c("CH4", "N2O", "CO2", "Gas1", "Gas2", "Gas3")

  data <- data %>%
    dplyr::mutate(Date = lubridate::parse_date_time(as.character(Date), orders = c("dmy", "mdy", "ymd"))) %>% # deal with date format parsing it into parts
    dplyr::filter(!is.na(Date))

  corrected_cols <- stringr::str_remove(names(data), "_flux_corrected")
  mgm2h_cols <- stringr::str_remove(names(data), "_flux_mgm2h")

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
    dplyr::group_by(Plot, Treatment) %>%
    dplyr::summarise(
      min_date = as.Date(min(Date)),
      max_date = as.Date(max(Date)),
      .groups = "drop"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(dates = list(seq(from = min_date, to = max_date, by = "day"))) %>%
    tidyr::unnest(dates) %>%
    dplyr::select(Plot, Treatment, Date = dates)

  flux_df <- full_dates_df %>%
    dplyr::left_join(data, by = c("Plot", "Treatment", "Date")) %>%
    dplyr::arrange(Plot, Treatment, Date)

  for (gas in present_gases) {

    flux_col <- paste0(gas, suffix)

    emission_col <- paste0(gas, "_emission_kg_ha")

    flux_df <- flux_df %>%
      dplyr::group_by(Plot, Treatment) %>%
      dplyr::arrange(Date) %>%
      dplyr::mutate("{flux_col}" := zoo::na.approx(.data[[flux_col]], x = Date, na.rm = FALSE, rule = 2)) %>%
      dplyr::ungroup()

    flux_df <- flux_df %>%
      dplyr::mutate("{emission_col}" := .data[[flux_col]] * 0.24)
  }

  kept_fluxes <- paste0(present_gases, suffix)

  kept_emissions <- paste0(present_gases, "_emission_kg_ha")

  cumulative_emissions <- flux_df %>%
    dplyr::group_by(Plot, Treatment) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(kept_emissions), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

  flux_df <- flux_df %>%
    dplyr::select(Date, Plot, Treatment, ID, dplyr::any_of(kept_fluxes), dplyr::any_of(kept_emissions))

  assign(paste0("daily_", data_name), flux_df, envir = .GlobalEnv)
  assign(paste0("acc_", data_name), cumulative_emissions, envir = .GlobalEnv)

  invisible(NULL)
}

## 1.2. Tests ####
# Using ppm2flux() outputs for tests
# Notes:
# - Try with input data frames with different formats for the Date column. So far the function works for format: "27-Jun-24"
# Test 1: Using as input a dataframe with flux corrections:
# flux2acc(flux_df_testA) # ouputs: acc_flux_df_testA (cumulative emissions) and daily_flux_df_testA (linear flux interpolation).
# Test 2: Using as input a dataframe without flux corrections:
# flux2acc(flux_df_testM) # ouputs: acc_flux_df_testM (cumulative emissions) and daily_flux_df_testM (linear flux interpolation).

#' Calculates Global Warming Potential.
#'
#' Global Warming potential in \eqn{CO_2} equivalents is calculated from cumulative emissions data in \eqn{kg \cdot ha^{-1}}.
#'
#' @param data A data frame containing cumulative emissions data in \eqn{kg \cdot ha^{-1}}.
#' @return Data frame containing calculated Global Warming potential in \eqn{CO_{2}} equivalents.
#' @param CH4_eq \eqn{CO_2} equivalents with a time horizon of 100 years for \eqn{CH_4}. Default is 27 \eqn{CO_{2}} equivalents according to IPCC, 2021.
#' @param N2O_eq \eqn{CO_2} equivalents on a 100-years horizon for \eqn{N_2O}. Default is 273 \eqn{CO_{2}} equivalents according to IPCC, 2021.
#' @importFrom dplyr %>%
#' @export
#' @examples
#' # Test 1: Using as input a data rame with flux corrections:
#' flux2acc(flux_df_testA) # ouputs: acc_flux_df_testA (cumulative emissions) and daily_flux_df_testA (linear flux interpolation).
#'
#' # Test 2: Using as input a data frame without flux corrections:
#' flux2acc(flux_df_testM) # ouputs: acc_flux_df_testM (cumulative emissions) and daily_flux_df_testM (linear flux interpolation).

acc2GWP <- function(data,
                     CH4_eq = 27, # default CO2 equivalents according to IPCC, 2021.
                     N2O_eq = 273) {

  data <- data %>%
    dplyr::mutate(CO2_eq_from_CH4 = CH4_emission_kg_ha * CH4_eq,
           CO2_eq_from_N20 = N2O_emission_kg_ha * N2O_eq,
           Proportion_GWP_CH4 = CO2_eq_from_CH4 / (CO2_eq_from_CH4 + CO2_eq_from_N20),
           Proportion_GWP_N20 = CO2_eq_from_N20 / (CO2_eq_from_CH4 + CO2_eq_from_N20),
           GWP_kg_CO2eq_ha = CO2_eq_from_CH4 + CO2_eq_from_N20)

  return(data)
}

## 2.2. Tests ####
# GWP_test1 <- acc2GWP(acc_flux_df_testA)
# GWP_test2 <- acc2GWP(acc_flux_df_testM)
