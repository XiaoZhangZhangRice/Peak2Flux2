library(Peak2Flux2)

input_test_peak <- read.csv("data/Input_samples_No_ppm.csv")
std_test_peak <- read.csv("data/Std_Input.csv")

input_test_flux <- read.csv("data/Input_samples_ppm.csv")

peak2ppm(std_test_peak, input_test_peak)

ppm2flux(input_test_flux)

