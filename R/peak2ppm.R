#standard_data <- read.csv("C:/Users/zhang/Documents/GitHub/Peak2Flux2/data/Std_Input.csv")
#str(standard_data)
#sample_data <- read.csv("C:/Users/zhang/Documents/GitHub/Peak2Flux2/data/Input_samples_No_ppm.csv")

# Step 1: Split the dataframe by GC_Run
split_by_GC_Run <- function(data) {
  split(data, data$GC_Run)
}

# Step 2 & 3: Perform linear regression for each gas and extract m, c, and r² with R² alert
calibration_curves <- function(data, r_square_threshold = 0.9) {

  # Initialize an empty list to store results
  results <- list(
    CH4 = list(c = NA, m = NA, r2 = NA),
    CO2 = list(c = NA, m = NA, r2 = NA),
    N2O = list(c = NA, m = NA, r2 = NA),
    Gas1 = list(c = NA, m = NA, r2 = NA),
    Gas2 = list(c = NA, m = NA, r2 = NA),
    Gas3 = list(c = NA, m = NA, r2 = NA)
  )

  alerts <- list()  # To store alerts about low R² values

  # CH4 model (handle NAs)
  if (all(!is.na(data$Std_CH4_PPM)) && all(!is.na(data$Std_CH4_Peak))) {
    CH4model <- lm(formula = Std_CH4_PPM ~ Std_CH4_Peak, data = data)
    results$CH4$c <- CH4model$coefficients[1]
    results$CH4$m <- CH4model$coefficients[2]
    results$CH4$r2 <- summary(CH4model)$r.squared
    if (results$CH4$r2 < r_square_threshold) {
      alerts <- c(alerts, paste("Warning: CH4 R² (", round(results$CH4$r2, 2), ") is below", r_square_threshold))
    }
  }

  # CO2 model (handle NAs)
  if (all(!is.na(data$Std_CO2_PPM)) && all(!is.na(data$Std_CO2_Peak))) {
    CO2model <- lm(formula = Std_CO2_PPM ~ Std_CO2_Peak, data = data)
    results$CO2$c <- CO2model$coefficients[1]
    results$CO2$m <- CO2model$coefficients[2]
    results$CO2$r2 <- summary(CO2model)$r.squared
    if (results$CO2$r2 < r_square_threshold) {
      alerts <- c(alerts, paste("Warning: CO2 R² (", round(results$CO2$r2, 2), ") is below", r_square_threshold))
    }
  }

  # N2O model (handle NAs)
  if (all(!is.na(data$Std_N2O_PPM)) && all(!is.na(data$Std_N2O_Peak))) {
    N2Omodel <- lm(formula = Std_N2O_PPM ~ Std_N2O_Peak, data = data)
    results$N2O$c <- N2Omodel$coefficients[1]
    results$N2O$m <- N2Omodel$coefficients[2]
    results$N2O$r2 <- summary(N2Omodel)$r.squared
    if (results$N2O$r2 < r_square_threshold) {
      alerts <- c(alerts, paste("Warning: N2O R² (", round(results$N2O$r2, 2), ") is below", r_square_threshold))
    }
  }

  # Gas1 model (handle NAs)
  if (all(!is.na(data$Std_Gas1_PPM)) && all(!is.na(data$Std_Gas1_Peak))) {
    Gas1model <- lm(formula = Std_Gas1_PPM ~ Std_Gas1_Peak, data = data)
    results$Gas1$c <- Gas1model$coefficients[1]
    results$Gas1$m <- Gas1model$coefficients[2]
    results$Gas1$r2 <- summary(Gas1model)$r.squared
    if (results$Gas1$r2 < r_square_threshold) {
      alerts <- c(alerts, paste("Warning: Gas1 R² (", round(results$Gas1$r2, 2), ") is below", r_square_threshold))
    }
  }

  # Gas2 model (handle NAs)
  if (all(!is.na(data$Std_Gas2_PPM)) && all(!is.na(data$Std_Gas2_Peak))) {
    Gas2model <- lm(formula = Std_Gas2_PPM ~ Std_Gas2_Peak, data = data)
    results$Gas2$c <- Gas2model$coefficients[1]
    results$Gas2$m <- Gas2model$coefficients[2]
    results$Gas2$r2 <- summary(Gas2model)$r.squared
    if (results$Gas2$r2 < r_square_threshold) {
      alerts <- c(alerts, paste("Warning: Gas2 R² (", round(results$Gas2$r2, 2), ") is below", r_square_threshold))
    }
  }

  # Gas3 model (handle NAs)
  if (all(!is.na(data$Std_Gas3_PPM)) && all(!is.na(data$Std_Gas3_Peak))) {
    Gas3model <- lm(formula = Std_Gas3_PPM ~ Std_Gas3_Peak, data = data)
    results$Gas3$c <- Gas3model$coefficients[1]
    results$Gas3$m <- Gas3model$coefficients[2]
    results$Gas3$r2 <- summary(Gas3model)$r.squared
    if (results$Gas3$r2 < r_square_threshold) {
      alerts <- c(alerts, paste("Warning: Gas3 R² (", round(results$Gas3$r2, 2), ") is below", r_square_threshold))
    }
  }

  # If there are any alerts, print them
  if (length(alerts) > 0) {
    cat(paste(alerts, collapse = "\n"), "\n")
  }

  # Convert results to a dataframe
  calibration_curve_results <- data.frame(
    Gas = c("CH4", "CO2", "N2O", "Gas1", "Gas2", "Gas3"),
    Intercept = sapply(results, function(x) x$c),
    Slope = sapply(results, function(x) x$m),
    R_Squared = sapply(results, function(x) x$r2)
  )

  return(calibration_curve_results)
}

# Step 4: Calculate ppm concentrations for samples using calibration curve (ppm = m * peak + c)
conc_calculator <- function(m, x, c) {
  m * x + c
}

# Main function to process the data
peak2ppm <- function(standard_data, sample_data, r_square_threshold = 0.9) {
  # Step 1: Split standard data by GC_Run
  split_data <- split_by_GC_Run(standard_data)

  # Step 2 & 3: For each GC_Run, perform the calibration and calculate ppm for the sample_data
  results <- lapply(split_data, function(gc_run_data) {
    # Perform calibration for this GC_Run
    calibration_results <- calibration_curves(gc_run_data, r_square_threshold)

    # Filter sample_data for this GC_Run
    matching_samples <- sample_data[sample_data$GC_Run == unique(gc_run_data$GC_Run), ]

    # Calculate CH4, CO2, N2O, and Gas 1, 2, and 3 concentrations using the calibration curves
    matching_samples$Sample_CH4_ppm <- conc_calculator(calibration_results$Slope[calibration_results$Gas == "CH4"], matching_samples$Sample_CH4_Peak, calibration_results$Intercept[calibration_results$Gas == "CH4"])
    matching_samples$Sample_CO2_ppm <- conc_calculator(calibration_results$Slope[calibration_results$Gas == "CO2"], matching_samples$Sample_CO2_Peak, calibration_results$Intercept[calibration_results$Gas == "CO2"])
    matching_samples$Sample_N2O_ppm <- conc_calculator(calibration_results$Slope[calibration_results$Gas == "N2O"], matching_samples$Sample_N2O_Peak, calibration_results$Intercept[calibration_results$Gas == "N2O"])
    matching_samples$Sample_Gas1_ppm <- conc_calculator(calibration_results$Slope[calibration_results$Gas == "Gas1"], matching_samples$Sample_Gas1_Peak, calibration_results$Intercept[calibration_results$Gas == "Gas1"])
    matching_samples$Sample_Gas2_ppm <- conc_calculator(calibration_results$Slope[calibration_results$Gas == "Gas2"], matching_samples$Sample_Gas2_Peak, calibration_results$Intercept[calibration_results$Gas == "Gas2"])
    matching_samples$Sample_Gas3_ppm <- conc_calculator(calibration_results$Slope[calibration_results$Gas == "Gas3"], matching_samples$Sample_Gas3_Peak, calibration_results$Intercept[calibration_results$Gas == "Gas3"])

    # Return the updated sample data with calculated concentrations
    return(matching_samples)
  })

  # Combine results into a single dataframe
  combined_results <- do.call(rbind, results)
  return(combined_results)
}



# Example usage
# final_data <- process_data(standard_data, sample_data)


# Example usage
# final_data <- process_data(standard_data, sample_data)


# Example usage
#final_data <- peak2ppm(standard_data, sample_data)
#calibration_results <- calibration_curves(standard_data)
