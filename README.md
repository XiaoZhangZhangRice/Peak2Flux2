# Seba's notes:

## To do:
- Generate an input data frame with more gases: So far we've only tested with CH4 and N2O
- Work the Molecule_ppm argument: this would also require a data frame with ppm info target atom ppm (e.g. C-CH4)
- Create "R2_Threshold" argument: the function is written with 0.7 as threshold. Define default with Zhang and modify the code so it considers alternative values if modified.

## Done:
- Timesteps argument: by default (4) it already calculates "drop-one" alternative models and generates a diagnostics pdf. When modified, output data frame contains only original model results.
