# Peak2Flux R Package for calculating GHG fluxes from Chromatography reads.

# 1. Project Description:
The Project allows GC users to transform GC peak data (time) into gas fluxes (mg m2 h<sup>-1</sup>) through 2 steps: 

### Step 1: peak2ppm
Calculates GHG concentrations (ppm) from chromatography peaks data frames (GC output).
- Required input: GC data (.csv) according to provided template.

### Step 2: ppm2flux
Calculates GHG fluxes (mg m2 h<sup>-1</sup>) from GHG concentrations (ppm) data frame (peak2ppm output)
- Required input: Step 1 output (.csv).
  
#### ppm2flux() function:
- **Timesteps argument**
   - Number of samples (vials) collected from a chamber on each sampling event (Date&Plot).
   - Default value: 4
- **Gas_mass arguments**: CH4_mass; N2O_mass; CO2_mass; Gas1_mass;  Gas2_mass; Gas3_mass.
   - Default value: 0 for all gases.
   - Note: if the user's GC output is atom instead of molecule ppm data (e.g. CH<sub>4</sub>-C instead of CH<sub>4</sub>), this argument schould be the target atom's molecular weight (e.g. 12 instead of 16). 
- **Diagnostics argument**:
   - Allows to activate/desactivate diagnostic plots.
   - Default value: FALSE (no diagnostic plots, plots are generated when changed to TRUE)

# 2. Project Layout:

- 2.1. /R: developed R scipts
- 2.2. /data: input dataframes for script development
- 2.3. /outputs: output files (i.e., data frames, plots, etc) from developed scripts
- 2.4. /R_test: scripts under ongoing development
