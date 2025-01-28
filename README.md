# Seba's notes:

## To do:
-  Generate an input data frame with more gases: So far we've only tested with CH4 and N2O
- Create a **ppm2flux_threshold_test.R**:
   - "R2_Threshold" argument: the function is written with 0.7 as threshold. Define default with Zhang and modify the code so it considers alternative values if modified.
- Make the code leaner (e.g. try fitting alternate model calculation and plots into for loops.
- Make the code faster (e.g. try lapply() and switch() approaches instead of for loops).
- Make te output leaner (e.g. remove columns without calculated values (e.g. for gases where Gas_mass = 0).

## Done:
1.  **In R_test**:
  - **ppm2flux_timesteps_test.R**:
    - **Timesteps argument**: by default (4) it already calculates "drop-one" alternative models and generates a diagnostics pdf. When modified, output data frame contains only original model results.
    - **Gas_mass**: Default = 0. ppm2flux funtion will only apply to gases where this default value is changed to the gas molecular weight.
    - **Diagnostics**: Default = FALSE. If changed to TRUE diagnostic plots are activated for gases with Gas_mass argument different to 0. Note: it only has an effect in Timesteps argument is left as default (4).
