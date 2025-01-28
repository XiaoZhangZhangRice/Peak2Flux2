# Seba's notes:

## To do:
- Generate an input data frame with more gases: So far we've only tested with CH4 and N2O
- Create a **ppm2flux_threshold_test.R**:
   - "R2_Threshold" argument: the function is written with 0.7 as threshold. Define default with Zhang and modify the code so it considers alternative values if modified.
   - Note: this test already contains arguments developed in ppm2flux_timesteps_test.R.
   - "Negative_flux" argument: default = TRUE. If left as default then corrections only consider positive rates for corrected fluxes ("Gas_flux_corrected").
- Alternative models's plots:
  - Fix x axis labels, currently fixed in breaks=c(0, 10, 20, 30). Make it according to user's sampling timings (Time_mins of input_test).
  - Title according to Date and Plot instead of Code.
- Make the code leaner (e.g. try fitting alternate model calculation and plots into for loops.
- Make the code faster (e.g. try lapply() and switch() approaches instead of for loops).
- Make the output leaner (e.g. remove columns without calculated values (e.g. for gases where Gas_mass = 0).

## Done:
1.  **In R_test**:
  - **ppm2flux_timesteps_test.R**:
    - **Timesteps argument**: by default (4) it already calculates "drop-one" alternative models and generates a diagnostics pdf. When modified, output data frame contains only original model results.
    - **Gas_mass**: Default = 0. ppm2flux funtion will only apply to gases where this default value is changed to the gas molecular weight.
    - **Diagnostics**: Default = FALSE. If changed to TRUE diagnostic plots are activated for gases with Gas_mass argument different to 0. Note: it only has an effect in Timesteps argument is left as default (4).
    - **R2_Threshold**: Default = 0.7. Defines an R<sup>2</sup> threshold under which corrected fluxes (slopes) from linear models are considered as 0 (no emissions). If the user decides to omit this parameter for flux corrections then it must be changed to 0. If increased, then the more strict the correction. **Note**: only applies to corrected fluxes ("Gas_flux_corrected"), it does not modify the original (complete, "Gas_flux_mgm2h") model flux calculation.
    - **Neg_Rate argument**: Default = TRUE.  When modified to FALSE then it includes an extra restriction to the flux correction that excludes alternative models that result in negative flux (even if resulting R<sup>2</sup> is higher than the original model). 
