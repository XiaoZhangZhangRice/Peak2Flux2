# Step 2: ppm2flux ####

# Load required packages:
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)

# 1. Determining function ####

ppm2flux_test <- function(data,
                          Timesteps = 4,
                          CH4_mass = 0,
                          N2O_mass = 0,
                          CO2_mass = 0,
                          Gas1_mass = 0,
                          Gas2_mass = 0,
                          Gas3_mass = 0,
                          Diagnostics = FALSE,
                          R2_Threshold = 0.7,
                          Neg_Rate = TRUE) {

  Gases <- c("CH4", "N2O", "CO2", "Gas1", "Gas2", "Gas3")

  # data frame for flux calculation
  ppm_df <- data %>%
    mutate( Chamber_Temp_K = data$Chamber_Temp_C + 273,
            ID = paste0(Date,"_", Plot),
            Volume_m3 = ifelse(is.na(Volume_m3), Surface_Area_m2 * Height_m, Volume_m3), # Seba 19.03.25: before modification-> Volume_m3 = data$Surface_Area_m2 * data$Height_m,
            CH4_density_g_m3 = (CH4_mass / (82.0575 * Chamber_Temp_K)) * 1000000,
            N2O_density_g_m3 = (N2O_mass / (82.0575 * Chamber_Temp_K)) * 1000000,
            CO2_density_g_m3 = (CO2_mass / (82.0575 * Chamber_Temp_K)) * 1000000,
            Gas1_density_g_m3 = (Gas1_mass / (82.0575 * Chamber_Temp_K)) * 1000000,
            Gas2_density_g_m3 = (Gas2_mass / (82.0575 * Chamber_Temp_K)) * 1000000,
            Gas3_density_g_m3 = (Gas3_mass / (82.0575 * Chamber_Temp_K)) * 1000000,
            CH4_byMass_mgm3 = (CH4_density_g_m3 * Sample_CH4_ppm) / 1000,
            CH4_byMass_mgm2 = (CH4_byMass_mgm3 * Volume_m3) / Surface_Area_m2,
            N2O_byMass_mgm3 = (N2O_density_g_m3 * Sample_N2O_ppm) / 1000,
            N2O_byMass_mgm2 = (N2O_byMass_mgm3 * Volume_m3) / Surface_Area_m2,
            CO2_byMass_mgm3 = (CO2_density_g_m3 * Sample_CO2_ppm) / 1000,
            CO2_byMass_mgm2 = (CO2_byMass_mgm3 * Volume_m3) / Surface_Area_m2,
            Gas1_byMass_mgm3 = (Gas1_density_g_m3 * Sample_Gas1_ppm) / 1000,
            Gas1_byMass_mgm2 = (Gas1_byMass_mgm3 * Volume_m3) / Surface_Area_m2,
            Gas2_byMass_mgm3 = (Gas2_density_g_m3 * Sample_Gas2_ppm) / 1000,
            Gas2_byMass_mgm2 = (Gas2_byMass_mgm3 * Volume_m3) / Surface_Area_m2,
            Gas3_byMass_mgm3 = (Gas3_density_g_m3 * Sample_Gas3_ppm) / 1000,
            Gas3_byMass_mgm2 = (Gas3_byMass_mgm3 * Volume_m3) / Surface_Area_m2)

  flux_df_test_timesteps <- length(unique(ppm_df$Time_mins)) # checks timesteps within input data frame, used after as logic value to decide if alternative models should be calculated

  flux_df_test <- ppm_df %>%  # data frame for flux calculation
    distinct(ID, Date, Plot, Treatment) # Creates a data frame with unique values for certain columns

  if (Timesteps == 4 & flux_df_test_timesteps == 4) { # in case Timesteps argument is left as default and input data frame has 4 timesteps per sampling event

    flux_df_test <- cbind(empty_column1 = NA, flux_df_test)
    new_columns <- paste0("col", 1:116) ## Adds empty columns
    flux_df_test[new_columns] <- NA

    colnames(flux_df_test) = c("Code_Nr", "ID", "Date", "Plot", "Treatment", "CH4_flux_mgm2h", "R2_CH4",
                               "p_CH4", "N2O_flux_mgm2h", "R2_N2O", "p_N2O", "CO2_flux_mgm2h",
                               "R2_CO2", "p_CO2", "Gas1_flux_mgm2h", "R2_Gas1", "p_Gas1",
                               "Gas2_flux_mgm2h", "R2_Gas2", "p_Gas2", "Gas3_flux_mgm2h", "R2_Gas3",
                               "p_Gas3", "CH4_flux_Alt1", "R2_CH4_Alt1", "p_CH4_Alt1", "CH4_flux_Alt2",
                               "R2_CH4_Alt2", "p_CH4_Alt2", "CH4_flux_Alt3", "R2_CH4_Alt3", "p_CH4_Alt3",
                               "CH4_flux_Alt4", "R2_CH4_Alt4", "p_CH4_Alt4", "CH4_model", "CH4_flux_corrected",
                               "R2_CH4_corrected", "Logic_CH4", "N2O_flux_Alt1", "R2_N2O_Alt1", "p_N2O_Alt1",
                               "N2O_flux_Alt2", "R2_N2O_Alt2", "p_N2O_Alt2", "N2O_flux_Alt3", "R2_N2O_Alt3",
                               "p_N2O_Alt3", "N2O_flux_Alt4", "R2_N2O_Alt4", "p_N2O_Alt4","N2O_model",
                               "N2O_flux_corrected", "R2_N2O_corrected", "Logic_N2O", "CO2_flux_Alt1", "R2_CO2_Alt1",
                               "p_CO2_Alt1", "CO2_flux_Alt2", "R2_CO2_Alt2", "p_CO2_Alt2", "CO2_flux_Alt3",
                               "R2_CO2_Alt3", "p_CO2_Alt3", "CO2_flux_Alt4", "R2_CO2_Alt4", "p_CO2_Alt4",
                               "CO2_model", "CO2_flux_corrected", "R2_CO2_corrected", "Logic_CO2", "CO2_flux_Alt1",
                               "R2_CO2_Alt1", "p_CO2_Alt1", "CO2_flux_Alt2", "R2_Gas1_Alt2", "p_Gas1_Alt2",
                               "Gas1_flux_Alt3", "R2_Gas1_Alt3", "p_Gas1_Alt3", "Gas1_flux_Alt4", "R2_Gas1_Alt4",
                               "p_Gas1_Alt4", "Gas1_model", "Gas1_flux_corrected", "R2_Gas1_corrected", "Logic_Gas1",
                               "Gas1_flux_Alt1", "R2_Gas1_Alt1", "p_Gas1_Alt1", "Gas1_flux_Alt2", "R2_Gas2_Alt2",
                               "p_Gas2_Alt2", "Gas2_flux_Alt3", "R2_Gas2_Alt3", "p_Gas2_Alt3", "Gas2_flux_Alt4",
                               "R2_Gas2_Alt4", "p_Gas2_Alt4", "Gas2_model", "Gas2_flux_corrected", "R2_Gas2_corrected",
                               "Logic_Gas2", "Gas3_flux_Alt1", "R2_Gas3_Alt1", "p_Gas3_Alt1", "Gas3_flux_Alt2",
                               "R2_Gas3_Alt2", "p_Gas3_Alt2", "Gas3_flux_Alt3", "R2_Gas3_Alt3", "p_Gas3_Alt3",
                               "Gas3_flux_Alt4", "R2_Gas3_Alt4", "p_Gas3_Alt4", "Gas3_model", "Gas3_flux_corrected",
                               "R2_Gas3_corrected", "Logic_Gas3")

  } else { # in case Timestep argument is modified (not calculating alternative models)

    flux_df_test <- cbind(empty_column1 = NA, flux_df_test)
    new_columns <- paste0("col", 1:20) ## Adds empty columns
    flux_df_test[new_columns] <- NA

    colnames(flux_df_test) = c("Code_Nr", "ID", "Date", "Plot", "Treatment", "CH4_flux_mgm2h", "R2_CH4",
                               "p_CH4", "N2O_flux_mgm2h", "R2_N2O", "p_N2O", "CO2_flux_mgm2h",
                               "R2_CO2", "p_CO2", "Gas1_flux_mgm2h", "R2_Gas1", "p_Gas1",
                               "Gas2_flux_mgm2h", "R2_Gas2", "p_Gas2", "Gas3_flux_mgm2h", "R2_Gas3",
                               "p_Gas3")
  }

  flux_df_test$Code_Nr <- 1:nrow(flux_df_test)

  ## 1.1. Flux calculation ####

  for(gas in Gases) {

    # Variable names for each gas iteration
    mass_var <- paste0(gas, "_byMass_mgm2")
    flux_var <- paste0(gas, "_flux_mgm2h")
    flux_var_Alt1 <- paste0(gas, "_flux_Alt1")
    flux_var_Alt2 <- paste0(gas, "_flux_Alt2")
    flux_var_Alt3 <- paste0(gas, "_flux_Alt3")
    flux_var_Alt4 <- paste0(gas, "_flux_Alt4")
    r2_var <- paste0("R2_", gas)
    r2_var_Alt1 <- paste0("R2_", gas, "_Alt1")
    r2_var_Alt2 <- paste0("R2_", gas, "_Alt2")
    r2_var_Alt3 <- paste0("R2_", gas, "_Alt3")
    r2_var_Alt4 <- paste0("R2_", gas, "_Alt4")
    p_var <- paste0("p_", gas)
    p_var_Alt1 <- paste0("p_", gas, "_Alt1")
    p_var_Alt2 <- paste0("p_", gas, "_Alt2")
    p_var_Alt3 <- paste0("p_", gas, "_Alt3")
    p_var_Alt4 <- paste0("p_", gas, "_Alt4")
    gas_mass <- get(paste0(gas, "_mass"))

    if (Diagnostics == TRUE & gas_mass != 0 & Timesteps == 4 & flux_df_test_timesteps == 4) { # in case Diagnostics argument changed to TRUE

      pdf_filename <- paste0("outputs/flux_diagnostics_", gas, ".pdf")
      pdf(file = pdf_filename)  # Creates a diagnostics file with plots for each ID and its alternative models

    }

    for (i in 1:length(flux_df_test$ID)) {
      Code_i <- flux_df_test$ID[i]
      Filt_i <- filter(ppm_df, ppm_df$ID == Code_i) # if returned as Time-Series, re-run library(dplyr)

      if (gas_mass != 0) { # Linear model will be calculated only for gases with molecular weight defined in Gas_mass arguments

        ### 1.1.1. Loop section 1: Rate calculation ####
        lm_i <- lm(as.formula(paste0(mass_var, "~Time_mins")), data = Filt_i)
        flux_df_test[[flux_var]][i] <- coef(lm_i)[2]*60 # Returns Gas_flux_mgm2h for each Code.
        flux_df_test[[r2_var]][i] <- summary(lm_i)$r.squared # Returns R2_Gas for each Code.
        lmp_i <- function (modelobject) {  # Function created to call later the model's p-value
          if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
          f <- summary(modelobject)$fstatistic
          p <- pf(f[1], f[2], f[3], lower.tail = F)
          attributes(p) <- NULL
          return(p)}
        flux_df_test[[p_var]][i] <- lmp_i(lm_i) # Returns p-value for each Code.

      }

      if (Timesteps == 4 & gas_mass != 0 & flux_df_test_timesteps == 4) { # in case Timesteps argument is left as default then alternative models are calculated

        ### 1.1.2. Loop section 2: Rate correction ####

        ## Fitting 4 alternative "3-values" models (each one removing one time-step)

        ## Alt_1: excluding concentration from time step T0
        Filt_Alt1i <- if(is.na(Filt_i[[mass_var]][1]) == TRUE) {Filt_i} else {Filt_i[-1,]} # Filters excluding first concentration. In case there are NA
        # it doesn't exclude values.
        lm_Alt1i <- lm(as.formula(paste0(mass_var, "~Time_mins")), data=Filt_Alt1i) # Linear model for these 3 values
        flux_df_test[[flux_var_Alt1]][i] <- coef(lm_Alt1i)[2]*60 # Returns flux_mgm2h for each Code.
        flux_df_test[[r2_var_Alt1]][i] <- summary(lm_Alt1i)$r.squared # Returns R2 for each Code.
        flux_df_test[[p_var_Alt1]][i] <- lmp_i(lm_Alt1i) # Returns p-value for each Code.

        #### Alt_2: excluding concentration from time step T1
        Filt_Alt2i <- if(is.na(Filt_i[[mass_var]][2]) == TRUE) {Filt_i} else {Filt_i[-2,]} # Filters excluding second concentration. In case the value to be
        # excluded is NA it doesn't exclude values
        lm_Alt2i <- lm(as.formula(paste0(mass_var, "~Time_mins")), data=Filt_Alt2i) # Linear model for these 3 values
        flux_df_test[[flux_var_Alt2]][i] <- coef(lm_Alt2i)[2]*60 # Returns flux_mgm2h for each Code.
        flux_df_test[[r2_var_Alt2]][i] <- summary(lm_Alt2i)$r.squared # Returns R2for each Code.
        flux_df_test[[p_var_Alt2]][i] <- lmp_i(lm_Alt2i) # Returns p-value for each Code.

        #### Alt_3: excluding concentration from time step T2
        Filt_Alt3i <- if(is.na(Filt_i[[mass_var]][3]) == TRUE) {Filt_i} else {Filt_i[-3,]} # Filters excluding third concentration. In case the value to be
        # excluded is NA it doesn't exclude values
        lm_Alt3i <- lm(as.formula(paste0(mass_var, "~Time_mins")), data=Filt_Alt3i) # Linear model for these 3 values
        flux_df_test[[flux_var_Alt3]][i] <- coef(lm_Alt3i)[2]*60 # Returns _flux_mgm2h for each Code.
        flux_df_test[[r2_var_Alt3]][i] <- summary(lm_Alt3i)$r.squared # Returns R2 for each Code.
        flux_df_test[[p_var_Alt3]][i] <- lmp_i(lm_Alt3i) # Returns p-value for each Code.

        #### Alt_4: excluding concentration from time step T3
        Filt_Alt4i <- if(is.na(Filt_i[[mass_var]][4]) == TRUE) {Filt_i} else {Filt_i[-4,]} # Filters excluding fourth concentration. In case the value to be
        # excluded is NA it doesn't exclude values
        lm_Alt4i <- lm(as.formula(paste0(mass_var, "~Time_mins")), data=Filt_Alt4i) # Linear model for these 3 values
        flux_df_test[[flux_var_Alt4]][i] <- coef(lm_Alt4i)[2]*60 # Returns flux_mgm2h for each Code.
        flux_df_test[[r2_var_Alt4]][i] <- summary(lm_Alt4i)$r.squared # Returns R2 for each Code.
        flux_df_test[[p_var_Alt4]][i] <- lmp_i(lm_Alt4i) # Returns p-value for each Code.

        ## Loop section 2.2: Including restrictions in the loop with nested if_else:
        # Here an if else must let users define the threshold value (or if they want this restriction at all)

        flux_df_test[[paste0(gas, "_flux_corrected")]][i] <- if((flux_df_test[[r2_var]][i] > R2_Threshold) &
                                                                (Neg_Rate || coef(lm_i)[2] > 0)) {flux_df_test[[flux_var]][i]
        } else if((flux_df_test[[r2_var]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt1]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt2]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt3]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt4]][i] < R2_Threshold)) {0
        } else if((flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt4]][i])) {flux_df_test[[mass_va]][i]
        } else if ((flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt1i)[2] > 0)) {flux_df_test[[flux_var_Alt1]][i]
        } else if ((flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt2i)[2] > 0))  {flux_df_test[[flux_var_Alt2]][i]
        } else if ((flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt3i)[2] > 0)) {flux_df_test[[flux_var_Alt3]][i]
        } else if ((flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (Neg_Rate || coef(lm_Alt4i)[2] > 0)) {flux_df_test[[flux_var_Alt4]][i]
        } else {0}

        # ## Column with chosen model:
        flux_df_test[[paste0(gas, "_model")]][i] <- if((flux_df_test[[r2_var]][i] > R2_Threshold) &
                                                       (Neg_Rate || coef(lm_i)[2] > 0)) {"Complete model"
        } else if((flux_df_test[[r2_var]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt1]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt2]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt3]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt4]][i] < R2_Threshold)) {"No flux"
        } else if((flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt4]][i])){"Complete model"
        } else if ((flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt1i)[2] > 0)) {"Alternative model 1"
        } else if ((flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt2i)[2] > 0)) {"Alternative model 2"
        } else if ((flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt3i)[2] > 0)) {"Alternative model 3"
        } else if ((flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (Neg_Rate || coef(lm_Alt4i)[2] > 0)) {"Alternative model 4"
        } else {"No flux"}

        # ## Loop section 2.3: Calculating R2 according to the applied correction (if any):
        flux_df_test[[paste0("R2_", gas, "_corrected")]][i] <- if((flux_df_test[[r2_var]][i] > R2_Threshold) &
                                                                  (Neg_Rate || coef(lm_i)[2] > 0)) {flux_df_test[[r2_var]][i]
        } else if((flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt4]][i])) {flux_df_test[[r2_var]][i]
        } else if((flux_df_test[[r2_var]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt1]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt2]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt3]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt4]][i] < R2_Threshold)) {0
        } else if ((flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt1i)[2] > 0)) {flux_df_test[[r2_var_Alt1]][i]
        } else if ((flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt2i)[2] > 0)) {flux_df_test[[r2_var_Alt2]][i]
        } else if ((flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt3i)[2] > 0)) {flux_df_test[[r2_var_Alt3]][i]
        } else if ((flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (Neg_Rate || coef(lm_Alt4i)[2] > 0)) {flux_df_test[[r2_var_Alt4]][i]
        } else {0}

        # ## Add column with method and rate selection logic:
        flux_df_test[[paste0("Logic_", gas)]][i] <- if((flux_df_test[[r2_var]][i] > R2_Threshold) &
                                                       (Neg_Rate || coef(lm_i)[2] > 0)) {"Complete model has R2 > R2_Threshold"
        } else if((flux_df_test[[r2_var]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt1]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt2]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt3]][i] < R2_Threshold) &
                  (flux_df_test[[r2_var_Alt4]][i] < R2_Threshold)) {"No model achieves  R2 threshold"
        } else if((flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                  (flux_df_test[[r2_var]][i] > flux_df_test[[r2_var_Alt4]][i]))  {"Complete model achieves the highest R2"
        } else if ((flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (flux_df_test[[r2_var_Alt1]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt1i)[2] > 0)) {"Complete model has R2 < R2_Threshold and Alt. 1 achieves the highest R2 (> R2_Threshold)"
        } else if ((flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (flux_df_test[[r2_var_Alt2]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt2i)[2] > 0)) {"Complete model has R2 < R2_Threshold and Alt. 2 achieves the highest R2 (> R2_Threshold)"
        } else if ((flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt3]][i] > flux_df_test[[r2_var_Alt4]][i]) &
                   (Neg_Rate || coef(lm_Alt3i)[2] > 0)) {"Complete model has R2 < R2_Threshold and Alt. 3 achieves the highest R2 (> R2_Threshold)"
        } else if ((flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt1]][i]) &
                   (flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt2]][i]) &
                   (flux_df_test[[r2_var_Alt4]][i] > flux_df_test[[r2_var_Alt3]][i]) &
                   (Neg_Rate || coef(lm_Alt4i)[2] > 0)) {"Complete model has R2 < R2_Threshold and Alt. 4 achieves the highest R2 (> R2_Threshold)"
        } else {"Alternative model achieves higher R2 but negative rate"}


        ## 1.2. Diagnostic plots ####

        if (Diagnostics == TRUE & gas_mass != 0 & Timesteps == 4 & flux_df_test_timesteps == 4) { # in case Diagnostics argument changed to TRUE: creating diagnostic plots

          ## Plot - Original values (before corrections):
          # This is the plot diagnostic section, either as an independent function or determined by an argument for users to decide if they want this output or not.

          Plot_i <- ggplot(data = Filt_i, aes(x=Time_mins, y=!!sym(paste0(gas, "_byMass_mgm2")))) + # !!sym() to dynamically reference the column in Filt_i
            geom_point() +
            xlab("Sample time (min)") +
            ylab(paste0(gas, " by mass (mgm2)")) +
            ggtitle(paste("ID = ", Filt_i$ID[1], "; Complete")) +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_x_continuous(breaks=c(Filt_i$Time_mins[1], Filt_i$Time_mins[2], Filt_i$Time_mins[3], Filt_i$Time_mins[4])) +
            stat_poly_line() +
            stat_poly_eq() +
            annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(flux_df_test[[flux_var]][i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)

          ## Plot Alt_1:
          Plot_Alt_1 <- ggplot(data = Filt_Alt1i, aes(x=Time_mins, y=!!sym(paste0(gas, "_byMass_mgm2")))) + # !!sym() to dynamically reference the column in Filt_i
            geom_point() +
            xlab("Sample time (min)") +
            ylab(paste0(gas, " by mass (mgm2)")) +
            ggtitle(paste("Alt. Model 1")) +
            theme(plot.title = element_text(hjust = 0.5)) +
            scale_x_continuous(breaks=c(Filt_i$Time_mins[1], Filt_i$Time_mins[2], Filt_i$Time_mins[3], Filt_i$Time_mins[4])) +
            stat_poly_line() +
            stat_poly_eq()  +
            annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(flux_df_test[[flux_var_Alt1]][i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)

          ## Plot Alt_2:
          Plot_Alt_2 <- ggplot(data = Filt_Alt2i, aes(x=Time_mins, y=!!sym(paste0(gas, "_byMass_mgm2")))) + # !!sym() to dynamically reference the column in Filt_i
            geom_point() +
            xlab("Sample time (min)") +
            ylab(paste0(gas, " by mass (mgm2)")) +
            ggtitle(paste("Alt. Model 2")) +
            theme(plot.title = element_text(hjust = 0.5))+
            scale_x_continuous(breaks=c(Filt_i$Time_mins[1], Filt_i$Time_mins[2], Filt_i$Time_mins[3], Filt_i$Time_mins[4])) +
            stat_poly_line() +
            stat_poly_eq()  +
            annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(flux_df_test[[flux_var_Alt2]][i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)

          ## Plot Alt_3:
          Plot_Alt_3 <- ggplot(data = Filt_Alt3i, aes(x=Time_mins, y=!!sym(paste0(gas, "_byMass_mgm2")))) + # !!sym() to dynamically reference the column in Filt_i
            geom_point() +
            xlab("Sample time (min)") +
            ylab(paste0(gas, " by mass (mgm2)")) +
            ggtitle(paste("Alt. Model 3")) +
            theme(plot.title = element_text(hjust = 0.5))+
            scale_x_continuous(breaks=c(Filt_i$Time_mins[1], Filt_i$Time_mins[2], Filt_i$Time_mins[3], Filt_i$Time_mins[4])) +
            stat_poly_line() +
            stat_poly_eq()  +
            annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(flux_df_test[[flux_var_Alt3]][i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)

          ## Plot Alt_4:
          Plot_Alt_4 <- ggplot(data = Filt_Alt4i, aes(x=Time_mins, y=!!sym(paste0(gas, "_byMass_mgm2")))) + # !!sym() to dynamically reference the column in Filt_i
            geom_point() +
            xlab("Sample time (min)") +
            ylab(paste0(gas, " by mass (mgm2)")) +
            ggtitle(paste("Alt. Model 4")) +
            theme(plot.title = element_text(hjust = 0.5))+
            scale_x_continuous(breaks=c(Filt_i$Time_mins[1], Filt_i$Time_mins[2], Filt_i$Time_mins[3], Filt_i$Time_mins[4])) +
            stat_poly_line() +
            stat_poly_eq()  +
            annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(flux_df_test[[flux_var_Alt4]][i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)

          gas_arrange <- ggarrange(Plot_i, Plot_Alt_1, Plot_Alt_2, Plot_Alt_3, Plot_Alt_4, ncol = 2, nrow = 3)

          print(gas_arrange)

        } # closing if() for diagnostic plots (only running if diagnostic plots are activated, i.e. Diagnostics argument modified to TRUE)
      } # closing if() for Timesteps (only running for default Timesteps argument)
    } # closing for() loop for each ID data subset (iterates for each ID (unique(Date&Plot)))

    if (Diagnostics == TRUE & gas_mass != 0 & Timesteps == 4 & flux_df_test_timesteps == 4) { # in case Diagnostics argument modified to TRUE: activates diagnostic plots
      dev.off()
    } # closing if() for diagnostic plots (only running if diagnostic plots are activated, i.e. argument modified to TRUE)

  } # closing for() loop for each Gas (iterates through Gases <- c("CH4", "N2O", "CO2", "Gas1", "Gas2", "Gas3"))

  flux_df_test <- flux_df_test %>% # removing empty columns (i.e. Gases without input data)
    select(where(~ !all(is.na(.))))

  return(flux_df_test)

} # closing function()

# 2. Tests ####

input_test <- read.csv("data/Input_samples_ppm.csv")
input_test2 <- read.csv("data/Input_samples_Nawal.csv")
input_test3 <- read.csv("data/Input_samples_Nawal2.csv") # data input contains Volume_m3 info instead of Surface_Area_m2 and Height_m
input_test4 <- read.csv("data/Input_samples_DASIG_1.csv") # DASIG data set: 20240828_Dasig_Input, from Flooded Rice - Arkansas, USA, 2024 (test version)
input_test5 <- read.csv("data/Input_samples_DASIG_2.csv") # DASIG data set: complete (05_April to 28_Aug, 2024), from Flooded Rice - Arkansas, USA, 2024 (test version)

## Tests with CH4 and N2O ppm inputs (these must be then defined in Gas_mass arguments)

# input: input_test
flux_df_testA <- ppm2flux_test(input_test, CH4_mass = 16, N2O_mass = 44) # Test 1: Keeping all arguments as default - Output with alternative models but no diagnostics
flux_df_testB <- ppm2flux_test(input_test, CH4_mass = 16, N2O_mass = 44, Diagnostics = TRUE) # Test 2: Activating diagnostic plots - Output with alternative models and diagnostics
flux_df_testC <- ppm2flux_test(input_test, CH4_mass = 16, N2O_mass = 44, Timesteps = 5) # Test 3: Testing Timesteps - Output without alternative models nor diagnostics
flux_df_testD <- ppm2flux_test(input_test, CH4_mass = 16, N2O_mass = 44, Diagnostics = TRUE, Timesteps = 5) # Test 4: No Diagnostics and Timesteps conflicts - Output without alternative models nor diagnostics
flux_df_testE <- ppm2flux_test(input_test, CH4_mass = 16, N2O_mass = 44, Neg_Rate = FALSE) # Test 5: Testing Neg_Rate argument, does not consider alternative models resulting in negative rates.
flux_df_testF <- ppm2flux_test(input_test, CH4_mass = 16, N2O_mass = 44, Diagnostics = TRUE, Neg_Rate = FALSE) # Test 6: Tests Diagnostic plots with Neg_Rate modified to FALSE

# input: input_test2
flux_df_testG <- ppm2flux_test(input_test2, CH4_mass = 16, N2O_mass = 44) # Test 7: 2nd test data set. Keeping all arguments as default - Output with alternative models but no diagnostics
flux_df_testH <- ppm2flux_test(input_test2, CH4_mass = 16, N2O_mass = 44, Diagnostics = TRUE) # Test 8: Activating diagnostic plots - Output with alternative models and diagnostics

# input: input_test3
flux_df_testI <- ppm2flux_test(input_test3, CH4_mass = 16, N2O_mass = 44) # Test 8: Keeping all arguments as default - Output with alternative models but no diagnostics

# input: input_test4
flux_df_testJ <- ppm2flux_test(input_test4, CH4_mass = 16, N2O_mass = 44, CO2_mass = 44, Timesteps = 5) # Test 9: 3rd test data set (with 5 timesteps). Keeping all arguments as default - Output with alternative models but no diagnostics
flux_df_testK <- ppm2flux_test(input_test4, CH4_mass = 16, N2O_mass = 44, CO2_mass = 44, Diagnostics = TRUE) # Test 10: Activating diagnostic plots - Output without alternative models nor diagnostics (due to 5 timesteps within input data frame)

# input: input_test5
flux_df_testL <- ppm2flux_test(input_test5, CH4_mass = 16, N2O_mass = 44, CO2_mass = 44, Timesteps = 5) # Test 11: 4th test data set (with 5 timesteps). Keeping all arguments as default - Output with alternative models but no diagnostics
flux_df_testM <- ppm2flux_test(input_test5, CH4_mass = 16, N2O_mass = 44, CO2_mass = 44, Diagnostics = TRUE) # Test 12: Activating diagnostic plots - Output with alternative models and diagnostics

