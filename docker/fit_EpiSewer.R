#!/usr/bin/env Rscript
options(cmdstanr_no_ver_check=TRUE)

cat("--- Running EpiSewer in docker container ---\n")

library(cmdstanr)
source("/opt/utils_warnings.R")
source("/opt/utils_stan.R")

args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop(paste(
    "Please supply a file with the EpiSewer job",
    "and the path for the result file."
    ), call.=FALSE)
}

input_filepath <- args[1]
output_filepath <- args[2]

# check if input exists
if (!file.exists(input_filepath)) {
  stop(paste("File", input_filepath, "does not exist."), call.=FALSE)
}

# check if input is rds file
if (tools::file_ext(input_filepath) != "rds") {
  stop(paste("File", input_filepath, "is not an rds file."), call.=FALSE)
}

# read job file
job <- readRDS(input_filepath)

stan_model <- get_stan_model(use_docker = TRUE)
model_instance <- stan_model$load_model[[1]]()

fit_res <- tryCatch(fit_model(
  job = job,
  model = stan_model,
  model_instance = model_instance,
  run_silent = FALSE
), error = function(e) {
  warning(paste("Error fitting model:", e$message))
  return(list(error = e$message))
})

# save result
saveRDS(fit_res, output_filepath)
