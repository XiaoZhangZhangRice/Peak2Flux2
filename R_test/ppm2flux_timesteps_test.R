# f(x): ppm2flux ####

# Load required packages:
library(dplyr)
library(ggplot2)
library(ggpmisc)
library(ggpubr)

# 1. Input dataframe ####

input_test <- read.csv("data/Input_samples_ppm.csv")
# input_test <- read.csv("Input_files/Input_samples_ppm_short.csv")

# 2. Determining function ####

ppm2flux <- function(data, Timesteps = 4) { # function(data, method = "lm", Timesteps = "4", ppm_CH4 = "CH4", Volume_L = "NA")

  # Arguments:
  ## Timesteps -> Number of samples (vials) collected from a chamber on each sampling event (Date&Plot).
  ##              Default: 4 - If not changed then the f(x) calculates 4 "drop-one" alternative models; if changed then the f(x) calculates only the complete model.

  # data frame for flux calculation
  ppm_df <- data %>%
    mutate(
      Chamber_Temp_K = data$Chamber_Temp_C + 273,
      ID = paste0(Date, Plot),

      Volume_m3 = data$Surface_Area_m2 * data$Height_m,

      # Here should do an if(), the following is for GCs returning CH4 ppm (default), but the if() should work for ppm_CH4 = "C-CH4"

      CH4_density_g_m3 = (16 / (82.0575 * Chamber_Temp_K)) * 1000000,
      N2O_density_g_m3 = (44 / (82.0575 * Chamber_Temp_K)) * 1000000,
      CO2_density_g_m3 = (44 / (82.0575 * Chamber_Temp_K)) * 1000000,
      CH4_byMass_mgm3 = (CH4_density_g_m3 * Sample_CH4_ppm) / 1000,
      CH4_byMass_mgm2 = (CH4_byMass_mgm3 * Volume_m3) / Surface_Area_m2,
      N2O_byMass_mgm3 = (N2O_density_g_m3 * Sample_N2O_ppm) / 1000,
      N2O_byMass_mgm2 = (N2O_byMass_mgm3 * Volume_m3) / Surface_Area_m2,
      CO2_byMass_mgm3 = (CO2_density_g_m3 * Sample_CO2_ppm) / 1000,
      CO2_byMass_mgm2 = (CO2_byMass_mgm3 * Volume_m3) / Surface_Area_m2)

  # data frame for flux calculation
  flux_df <- ppm_df %>%
    distinct(ID, Date) # Creates data frame with unique values for certain columns

  if (Timesteps == 4) { # in case Timestep argument is left as default

  flux_df <- cbind(empty_column1=NA,flux_df,empty_column2=NA, empty_column3=NA, empty_column4=NA, empty_column5=NA,
                   empty_column6=NA, empty_column7=NA, empty_column8=NA, empty_column9=NA, empty_column10=NA, empty_column11=NA,
                   empty_column12=NA, empty_column13=NA, empty_column14=NA, empty_column15=NA, empty_column16=NA, empty_column17=NA,
                   empty_column18=NA, empty_column19=NA, empty_column20=NA, empty_column21=NA, empty_column22=NA, empty_column23=NA,
                   empty_column24=NA, empty_column25=NA, empty_column26=NA)  ## Adds empty columns

  colnames(flux_df) = c("Code_Nr", "ID", "Date", "CH4_flux_mgm2h", "R2_CH4", "p_CH4", "N2O_flux_mgm2h",
                        "R2_N2O", "p_N2O", "CO2_flux_mgm2h", "R2_CO2", "p_CO2", "CH4_flux_Alt1", "R2_CH4_Alt1", "p_CH4_Alt1", "CH4_flux_Alt2",
                        "R2_CH4_Alt2", "p_CH4_Alt2", "CH4_flux_Alt3", "R2_CH4_Alt3", "p_CH4_Alt3", "CH4_flux_Alt4", "R2_CH4_Alt4", "p_CH4_Alt4",
                        "CH4_model", "CH4_flux_corrected", "R2_CH4_corrected", "Logic_CH4")

  } else { # in case Timestep argument is modified (not calculating alternative models)

  flux_df <- cbind(empty_column1=NA,flux_df,empty_column2=NA, empty_column3=NA, empty_column4=NA, empty_column5=NA,
                   empty_column6=NA, empty_column7=NA, empty_column8=NA, empty_column9=NA, empty_column10=NA, empty_column11=NA,
                   empty_column12=NA)  ## Adds empty columns

  colnames(flux_df) = c("Code_Nr", "ID", "Date", "CH4_flux_mgm2h", "R2_CH4", "p_CH4", "N2O_flux_mgm2h",
                        "R2_N2O", "p_N2O", "CO2_flux_mgm2h", "R2_CO2", "p_CO2")
  }

  flux_df$Code_Nr <- 1:nrow(flux_df)

  # flux calculation

  pdf('outputs/flux_diagnostics.pdf')  # Creates a diagnostics file with plots for each ID and its alternative models

  for (i in 1:length(flux_df$ID)) {
    Code_i <- flux_df$ID[i]
    Filt_i <- filter(ppm_df,ppm_df$ID == Code_i) # if returned as Time-Series, re-run library(dplyr)

    ## Loop section 1: Rate calculation.
    lm_i <- lm(CH4_byMass_mgm2~Time_mins, data=Filt_i)
    flux_df$CH4_flux_mgm2h[i] <- coef(lm_i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
    flux_df$R2_CH4[i] <- summary(lm_i)$r.squared # Returns R2_CH4 for each Code.
    lmp_i <- function (modelobject) {  # Function created to call later the model's p-value
      if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
      f <- summary(modelobject)$fstatistic
      p <- pf(f[1],f[2],f[3],lower.tail=F)
      attributes(p) <- NULL
      return(p)}
    flux_df$p_CH4[i] <- lmp_i(lm_i) # Returns p-value for each Code.


    if (Timesteps == 4) { # in case Timesteps argument is left as default

    ## Loop section 2: Rate correction.

    ## Fitting 4 alternative "3-values" models (each one removing one time-step) and a Log model:

    # Here an in for loop should iterate across gases (e.g. CH4, N2O, CO2), explore if it could be open to do it across any other gas defined in the input data frame.

    ## Alt_1: excluding concentration from time step T0
    Filt_Alt1i <- if(is.na(Filt_i$CH4_byMass_mgm2[1]) == TRUE) {Filt_i} else {Filt_i[-1,]} # Filters excluding one concentration. In case there are NA
    # it doesn't exclude values.
    lm_Alt1i <- lm(CH4_byMass_mgm2~Time_mins, data=Filt_Alt1i) # Linear model of these 3 values
    flux_df$CH4_flux_Alt1[i] <- coef(lm_Alt1i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
    flux_df$R2_CH4_Alt1[i] <- summary(lm_Alt1i)$r.squared # Returns R2_CH4 for each Code.
    flux_df$p_CH4_Alt1[i] <- lmp_i(lm_Alt1i) # Returns p-value for each Code.

    #### Alt_2: excluding concentration from time step T1
    Filt_Alt2i <- if(is.na(Filt_i$CH4_byMass_mgm2[2]) == TRUE) {Filt_i} else {Filt_i[-2,]} # Filters excluding one concentration. In case the value to be
    # excluded is NA it doesn't exclude values
    lm_Alt2i <- lm(CH4_byMass_mgm2~Time_mins, data=Filt_Alt2i) # Linear model of these 3 values
    flux_df$CH4_flux_Alt2[i] <- coef(lm_Alt2i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
    flux_df$R2_CH4_Alt2[i] <- summary(lm_Alt2i)$r.squared # Returns R2_CH4 for each Code.
    flux_df$p_CH4_Alt2[i] <- lmp_i(lm_Alt2i) # Returns p-value for each Code.

    #### Alt_3: excluding concentration from time step T2
    Filt_Alt3i <- if(is.na(Filt_i$CH4_byMass_mgm2[3]) == TRUE) {Filt_i} else {Filt_i[-3,]} # Filters excluding one concentration. In case the value to be
    # excluded is NA it doesn't exclude values
    lm_Alt3i <- lm(CH4_byMass_mgm2~Time_mins, data=Filt_Alt3i) # Linear model of these 3 values
    flux_df$CH4_flux_Alt3[i] <- coef(lm_Alt3i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
    flux_df$R2_CH4_Alt3[i] <- summary(lm_Alt3i)$r.squared # Returns R2_CH4 for each Code.
    flux_df$p_CH4_Alt3[i] <- lmp_i(lm_Alt3i) # Returns p-value for each Code.

    #### Alt_4: excluding concentration from time step T3
    Filt_Alt4i <- if(is.na(Filt_i$CH4_byMass_mgm2[4]) == TRUE) {Filt_i} else {Filt_i[-4,]} # Filters excluding one concentration. In case the value to be
    # excluded is NA it doesn't exclude values
    lm_Alt4i <- lm(CH4_byMass_mgm2~Time_mins, data=Filt_Alt4i) # Linear model of these 3 values
    flux_df$CH4_flux_Alt4[i] <- coef(lm_Alt4i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
    flux_df$R2_CH4_Alt4[i] <- summary(lm_Alt4i)$r.squared # Returns R2_CH4 for each Code.
    flux_df$p_CH4_Alt4[i] <- lmp_i(lm_Alt4i) # Returns p-value for each Code.

    # #### Log Model:
    # Log_mod_i<- lm(CH4_byMass_mgm2~log(Time_mins+1), data=Filt_i)
    # flux_df$CH4_flux_Log[i] <- coef(Log_mod_i)[2]*60 # Returns CH4_flux_mgm2h for each Code.
    # flux_df$R2_CH4_Log[i] <- summary(Log_mod_i)$r.squared # Returns R2_CH4 for each Code.
    # flux_df$p_CH4_Log[i] <- lmp_i(Log_mod_i) # Returns p-value for each Code.

    ## Loop section 2.2: Including restrictions in the loop with nested if_else:
    # Here an if else must let users define the threshold value (or if they want this restriction at all)

    flux_df$CH4_flux_corrected[i] <- if((flux_df$R2_CH4[i] < 0.7) & coef(lm_i)[2] < 0) {0
    } else if(flux_df$R2_CH4[i] > 0.7) {flux_df$CH4_flux_mgm2h[i]
    } else if((flux_df$R2_CH4[i] < 0.7) & (flux_df$R2_CH4_Alt1[i] < 0.7) & (flux_df$R2_CH4_Alt2[i] < 0.7) &
              (flux_df$R2_CH4_Alt3[i] < 0.7) & (flux_df$R2_CH4_Alt4[i] < 0.7)) {0
    } else if((flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt2[i]) &
              (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt3[i]) & (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt4[i])) {flux_df$CH4_flux_mgm2h[i]
    } else if ((flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt2[i]) & (flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt3[i]) &
               (flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt1i)[2] > 0)) {flux_df$CH4_flux_Alt1[i]
    } else if ((flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt3[i]) &
               (flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt2i)[2] > 0)) {flux_df$CH4_flux_Alt2[i]
    } else if ((flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt2[i]) &
               (flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt3i)[2] > 0)) {flux_df$CH4_flux_Alt3[i]
    } else if ((flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt2[i]) &
               (flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt3[i]) & (coef(lm_Alt4i)[2] > 0)) {flux_df$CH4_flux_Alt4[i]
    } else {0}

    # ## Column with chosen model:
    flux_df$CH4_model[i] <- if((flux_df$R2_CH4[i] < 0.7) & coef(lm_i)[2] < 0) {"No flux"
    } else if(flux_df$R2_CH4[i] > 0.7) {"Original"
    } else if((flux_df$R2_CH4[i] < 0.7) & (flux_df$R2_CH4_Alt1[i] < 0.7) & (flux_df$R2_CH4_Alt2[i] < 0.7) &
              (flux_df$R2_CH4_Alt3[i] < 0.7) & (flux_df$R2_CH4_Alt4[i] < 0.7)) {"No flux"
    } else if((flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt2[i]) &
              (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt3[i]) & (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt4[i])){"Original"
    } else if ((flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt2[i]) & (flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt3[i]) &
               (flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt1i)[2] > 0)) {"Alt. 1"
    } else if ((flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt3[i]) &
               (flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt2i)[2] > 0)) {"Alt. 2"
    } else if ((flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt2[i]) &
               (flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt3i)[2] > 0)) {"Alt. 3"
    } else if ((flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt2[i]) &
               (flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt3[i]) & (coef(lm_Alt4i)[2] > 0)) {"Alt. 4"
    } else {"No flux"}

    # ## Loop section 2.3: Calculating R2 according to the applied correction (if any):
    flux_df$R2_CH4_corrected[i] <- if((flux_df$R2_CH4[i] < 0.7) & coef(lm_i)[2] < 0) {0
    } else if(flux_df$R2_CH4[i] > 0.7) {flux_df$R2_CH4[i]
    } else if((flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt2[i]) &
              (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt3[i]) & (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt4[i])) {flux_df$R2_CH4[i]
    } else if((flux_df$R2_CH4[i] < 0.7) & (flux_df$R2_CH4_Alt1[i] < 0.7) & (flux_df$R2_CH4_Alt2[i] < 0.7) &
              (flux_df$R2_CH4_Alt3[i] < 0.7) & (flux_df$R2_CH4_Alt4[i] < 0.7)) {0
    } else if ((flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt2[i]) & (flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt3[i]) &
               (flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt1i)[2] > 0)) {flux_df$R2_CH4_Alt1[i]
    } else if ((flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt3[i]) &
               (flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt2i)[2] > 0)) {flux_df$R2_CH4_Alt2[i]
    } else if ((flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt2[i]) &
               (flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt3i)[2] > 0)) {flux_df$R2_CH4_Alt3[i]
    } else if ((flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt2[i]) &
               (flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt3[i]) & (coef(lm_Alt4i)[2] > 0)) {flux_df$R2_CH4_Alt4[i]
    } else {0}
    #
    # ## Add column with method and rate selection logic:
    flux_df$Logic_CH4[i] <- if((flux_df$R2_CH4[i] < 0.7) & coef(lm_i)[2] < 0) {"Original model has negative rate and R2 < 0.7"
    } else if(flux_df$R2_CH4[i] > 0.7) {"Original model has R2 > 0.7"
    } else if((flux_df$R2_CH4[i] < 0.7) & (flux_df$R2_CH4_Alt1[i] < 0.7) & (flux_df$R2_CH4_Alt2[i] < 0.7) &
              (flux_df$R2_CH4_Alt3[i] < 0.7) & (flux_df$R2_CH4_Alt4[i] < 0.7)) {"No model achieves  R2 threshold"
    } else if((flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt2[i]) &
              (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt3[i]) & (flux_df$R2_CH4[i] > flux_df$R2_CH4_Alt4[i]))  {"Original model achieves the highest R2"
    } else if ((flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt2[i]) & (flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt3[i]) &
               (flux_df$R2_CH4_Alt1[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt1i)[2] > 0)) {"Original model has R2 threshold and Alt. 1 achieves the highest
                  R2 (> 0.7) w/positive rate"
    } else if ((flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt3[i]) &
               (flux_df$R2_CH4_Alt2[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt2i)[2] > 0)) {"Original model has R2 threshold and Alt. 2 achieves the highest
                  R2 (> 0.7) w/positive rate"
    } else if ((flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt2[i]) &
               (flux_df$R2_CH4_Alt3[i] > flux_df$R2_CH4_Alt4[i]) & (coef(lm_Alt3i)[2] > 0)) {"Original model has R2 threshold and Alt. 3 achieves the highest
                  R2 (> 0.7) w/positive rate"
    } else if ((flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt1[i]) & (flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt2[i]) &
               (flux_df$R2_CH4_Alt4[i] > flux_df$R2_CH4_Alt3[i]) & (coef(lm_Alt4i)[2] > 0)) {"Original model has R2 threshold and Alt. 4 achieves the highest
                  R2 (> 0.7) w/positive rate"
    } else {"Alternative achieves higher R2 but negative rate"}

    ## Plot - Original values (before corrections):
    # This is the plot diagnostic section, either as an independent function or determined by an argument for users to decide if they want this output or not.

    Plot_i <- ggplot(data = Filt_i, aes(x=Time_mins, y=CH4_byMass_mgm2)) +
      geom_point() +
      xlab("Sample time (min)") +
      ylab("CH4 by mass (mgm2)") +
      ggtitle(paste("Code = ", i, "; Original")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(breaks=c(0, 10, 20, 30)) +
      stat_poly_line() +
      stat_poly_eq() +
      annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(flux_df$CH4_flux_mgm2h[i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)

    # print(Plot_i)

    ## Plot Alt_1:
    Plot_Alt_1 <- ggplot(data = Filt_Alt1i, aes(x=Time_mins, y=CH4_byMass_mgm2)) +
      geom_point() +
      xlab("Sample time (min)") +
      ylab("CH4 by mass (mgm2)") +
      ggtitle(paste("Alt. Model 1")) +
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_x_continuous(breaks=c(0, 10, 20, 30)) +
      stat_poly_line() +
      stat_poly_eq()  +
      annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(flux_df$CH4_flux_Alt1[i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)

    ## Plot Alt_2:
    Plot_Alt_2 <- ggplot(data = Filt_Alt2i, aes(x=Time_mins, y=CH4_byMass_mgm2)) +
      geom_point() +
      xlab("Sample time (min)") +
      ylab("CH4 by mass (mgm2)") +
      ggtitle(paste("Alt. Model 2")) +
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(breaks=c(0, 10, 20, 30)) +
      stat_poly_line() +
      stat_poly_eq()  +
      annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(flux_df$CH4_flux_Alt2[i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)

    ## Plot Alt_3:
    Plot_Alt_3 <- ggplot(data = Filt_Alt3i, aes(x=Time_mins, y=CH4_byMass_mgm2)) +
      geom_point() +
      xlab("Sample time (min)") +
      ylab("CH4 by mass (mgm2)") +
      ggtitle(paste("Alt. Model 3")) +
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(breaks = c(0, 10, 20, 30)) +
      stat_poly_line() +
      stat_poly_eq()  +
      annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(flux_df$CH4_flux_Alt3[i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)

    ## Plot Alt_4:
    Plot_Alt_4 <- ggplot(data = Filt_Alt4i, aes(x=Time_mins, y=CH4_byMass_mgm2)) +
      geom_point() +
      xlab("Sample time (min)") +
      ylab("CH4 by mass (mgm2)") +
      ggtitle(paste("Alt. Model 4")) +
      theme(plot.title = element_text(hjust = 0.5))+
      scale_x_continuous(breaks=c(0, 10, 20, 30)) +
      stat_poly_line() +
      stat_poly_eq()  +
      annotate(geom="text", -Inf, Inf, label=paste("Rate: ", round(flux_df$CH4_flux_Alt4[i], digits = 4),"mgm2h"), hjust = -0.25, vjust = 13)

    # ## Log Model:
    # lm_eqn <- function(Filt_i){                                                    ## Function defined to include R2 in log fitted plot.
    #   log_coef <- lm(CH4_byMass_mgm2~log(Time_mins+1), data=Filt_i);
    #   eq <- substitute(italic(y) == ~~italic(R)^2~"="~r2,
    #                    list(r2 = format(summary(log_coef)$r.squared, digits = 3)))
    #   as.character(as.expression(eq));
    # }

    # Plot_Log <- ggplot(data = Filt_i, aes(x=Time_mins, y=CH4_byMass_mgm2)) +
    #   geom_point() +
    #   xlab("Sample time (min)") +
    #   ylab("CH4 by mass (mgm2)") +
    #   ggtitle(paste("Log Model: y~log(x+1)")) +
    #   theme(plot.title = element_text(hjust = 0.5)) +
    #   stat_poly_line(data = Filt_i, method="lm",formula=y~log(x+1),fill="red") +
    #   stat_poly_eq(data = Filt_i, method="lm",formula=y~log(x+1))

    ## Arrange plots:
    CH4_arrange <- ggarrange(Plot_i, Plot_Alt_1, Plot_Alt_2, Plot_Alt_3, Plot_Alt_4, ncol = 2, nrow = 3)

    print(CH4_arrange)

    } else { # in case Timesteps argument is modified (not calculating alternative models)
    }

  }

  dev.off()

  return(flux_df)

}

flux_df <- ppm2flux(input_test, Timesteps = 4)

