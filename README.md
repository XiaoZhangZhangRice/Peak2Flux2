# Peak 2 FluxR Package for calculating GHG fluxes from Chromatography reads.

# 1. Project Description:

### Part 1: peak2ppm
Calculates GHG concentrations (ppm) from chromatography peaks data frames (GC output).

### Part 2: ppm2flux
Calculates GHG fluxes (mg m2 h<sup>-1</sup>) from GHG concentrations (ppm) data frame (peak2ppm output)
#### ppm2flux() function:
- **Gas_mass arguments**: CH4_mass; N2O_mass; CO2_mass; Gas1_mass;  Gas2_mass; Gas3_mass.
   - Default values: CH4_mass = 16, N2O_mass = 44, CO2_mass = 44, Gas1_mass = 0,  Gas2_mass = 0, Gas3_mass = 0 
   - Note: if the user's GC output is atom instead of molecule ppm data (e.g. CH<sub>4</sub>-C instead of CH<sub>4</sub>), this argument schould be the target atom's molecular weight (e.g. 12 instead of 16). 

# 2. Folders and files

- 2.1. /R: developed R scipts
- 2.2. /data: input dataframes for script development
- 2.3. /outputs: output files (i.e., data frames, plots, etc) from developed scripts
- 2.4. /R_test: scripts under ongoing development
