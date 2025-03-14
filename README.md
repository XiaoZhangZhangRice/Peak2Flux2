# Peak2Flux R Package for calculating GHG fluxes from Chromatography reads.

# 1. Project Description:
The Project allows GC users to transform GC peak data (time) into gas fluxes (mg m2 h<sup>-1</sup>) through 2 steps: 

### Step 1: peak2ppm
Calculates GHG concentrations (ppm) from chromatography peaks data frames (GC output).
- Required input: GC data (.csv) according to provided template.

### Step 2: ppm2flux
Calculates GHG fluxes (mg m2 h<sup>-1</sup>) from GHG concentrations (ppm) data frame (peak2ppm output)
- Required input: Step 1 output (.csv).
- Allows **flux corrections** in the form of alternative models and restriction settings (only for studies with 4 collected samples (vials) per sampling event):
   - **Alternative models**: Besides the complete model (with all four concentrations per sampling event), 4 "drop-one" alternative models are calculated, each dropping consecutively one concentration point (e.g. Alt. 1 has T0, T1, T2 instead of all T0 to T3). The function output is a data frame with fluxes resulting from complete and alternative linear models, their respective R<sup>2</sup>, p-values and a selected model according to argument settings (R2_Threshold and Neg_Rate arguments, see function description below). Fluxes resulting from these alternative models are defined as corrected fluxes. 
   - **Diagnostic plots**: Plots with complete and alternative linear models for all sampling events. 
   - Restrictions to corrected fluxes:
       - **R<sup>2</sup> threshold**: Sets an R<sup>2</sup> threshold for alternative models to be selected instead of the complete model. If the threshold is not surpassed by any model, then the resulting corrected flux will be 0 as no flux is assumed. 
       - **Negative fluxes**: The user can decide if alternative models resulting in negative fluxes can be considered (or not) as proper corrected fluxes.
  
#### ppm2flux() function:
- **Timesteps argument**
   - Number of samples (vials) collected from a chamber on each sampling event (Date&Plot).
   - Default value: 4.
   - If modified from default then **flux corrections** are not calculated.
- **Gas_mass arguments**: CH4_mass; N2O_mass; CO2_mass; Gas1_mass;  Gas2_mass; Gas3_mass.
   - Default value: 0 for all gases.
   - Note: if the user's GC output is atom instead of molecule ppm data (e.g. CH<sub>4</sub>-C instead of CH<sub>4</sub>), this argument schould be the target atom's molecular weight (e.g. 12 instead of 16). 
- **Diagnostics argument**:
   - Allows to activate/desactivate diagnostic plots.
   - Default value: FALSE (no diagnostic plots, plots are generated when changed to TRUE).
- **R2_Threshold argument**:
   - Defines an R<sup>2</sup> threshold under which corrected fluxes (slopes) from linear models are considered as 0 (no emissions). If the user decides to omit this parameter for flux corrections then it must be changed to 0. If increased, then the more strict the correction. 
   - Default value: 0.7.
   - Note: Only applies to corrected fluxes ("Gas_flux_corrected"), it does not modify the original (complete, "Gas_flux_mgm2h") model flux calculation.
- **Neg_Rate argument**:
   - When modified to FALSE then it includes an extra restriction to the flux correction that excludes alternative models that result in negative flux (even if resulting R<sup>2</sup> is higher than the original model). 
   - Default value: TRUE (does not consider this restriction and negative fluxes are accepted).
   - Note: If complete models result in negative flux (and R<sup>2</sup> is higher than the defined R2_Threshold) then the selected model will still be this complete model with regative flux.

# 2. Project Layout:

- 2.1. /R: developed R scipts
- 2.2. /data: input dataframes for script development
- 2.3. /outputs: output files (i.e., data frames, plots, etc) from developed scripts
- 2.4. /R_test: scripts under ongoing development
