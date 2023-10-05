#' Update meta information in modeldata based on available variables
#'
#' @description This update function is for all meta information that depends on
#'   modeldata from several functions. There are also functions which add meta
#'   information directly, in particular when this meta information cannot be
#'   solely inferred from modeldata. This function is designed such that calling
#'   it will never do harm to the modeldata object and not throw errors if
#'   something is missing in the modeldata object.
modeldata_update_metainfo <- function(modeldata) {
  if (modeldata_check(modeldata,
                      required = c("S", "D", "T"),
                      throw_error = FALSE
  )) {
    modeldata$meta_info$length_shedding <- with(modeldata, S + D + T)
  }
  if (modeldata_check(modeldata,
                      required = c("L", "S", "D", "T"),
                      throw_error = FALSE
  )) {
    modeldata$meta_info$length_I <- with(modeldata, L + S + D + T)
  }
  if (modeldata_check(modeldata,
                      required = c("L", "S", "D", "T", "G"),
                      throw_error = FALSE
  )) {
    modeldata$meta_info$length_R <- with(modeldata, L + S + D + T - G)
  }
  if (modeldata_check(
    modeldata,
    required = c(
      "measured_concentrations",
      "flow",
      "meta_info$composite_window",
      "meta_info$load_per_case"
    ),
    throw_error = FALSE
  )) {
    # crude descriptive estimate of cases at start of time series
    modeldata$meta_info$initial_cases_crude <-
      with(
        modeldata,
        0.1 + measured_concentrations[1] *
          mean(flow[1:meta_info$composite_window]) / meta_info$load_per_case
      )
  }
  return(modeldata)
}

modeldata_descriptions <- function() {
  descriptions <- list(
    "meta_info$composite_window" =
      "window length for composite samples in days",
    "meta_info$length_seeding" =
      "length of seeding phase for infections",
    "meta_info$length_I" =
      "number of days over which infections are modeled",
    "meta_info$length_R" =
      "number of days over which Rt is modeled",
    "meta_info$load_per_case" =
      "assumed overall load shed per individual",
    "meta_info$initial_cases_crude" =
      "empirical estimate for #cases at start of time period"
  )
  return(descriptions)
}

modeldata_var_requirements <- function() {
  requirements <- list(
    "meta_info$initial_cases_crude" = "flows_observe",
    "meta_info$length_seeding" = "generation_dist_assume",
    "meta_info$length_I" = c(
      "incubation_dist_assume",
      "shedding_dist_assume",
      "residence_dist_assume"),
    "meta_info$length_R" = c(
      "incubation_dist_assume",
      "shedding_dist_assume",
      "generation_dist_assume",
      "residence_dist_assume"
    )
  )
  return(requirements)
}

modeldata_defaults <- function() {
  defaults <- list(
    "K" = 0,
    "X" = numeric(0),
    "eta_prior" = numeric(0),
    "init$eta" = numeric(0)
  )
  return(defaults)
}

all_components <- function() {
  components <- c(
    "concentrations",
    "droplets",
    "LOD",
    "load_per_case",
    "load_variation",
    "flows",
    "generation_dist",
    "incubation_dist",
    "shedding_dist",
    "residence_dist",
    "R",
    "seeding",
    "infection_noise",
    "sample_effects",
    "noise"
  )
  return(components)
}
