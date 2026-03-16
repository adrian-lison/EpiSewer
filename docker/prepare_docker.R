source("R/utils_stan.R")
write_scalar_vars(modelfolder = "inst/stan", modelname = "EpiSewer_main.stan")
write_scalar_params(modelfolder = "inst/stan", modelname = "EpiSewer_main.stan")
