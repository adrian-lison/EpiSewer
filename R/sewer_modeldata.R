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
    modeldata$.metainfo$length_shedding <- with(modeldata, S + D + T)
  }
  if (modeldata_check(modeldata,
    required = c("L", "S", "D", "T"),
    throw_error = FALSE
  )) {
    modeldata$.metainfo$length_I <- with(modeldata, L + S + D + T)
  }
  if (modeldata_check(modeldata,
    required = c("L", "S", "D", "T", "G"),
    throw_error = FALSE
  )) {
    modeldata$.metainfo$length_R <- with(modeldata, L + S + D + T - G)
  }
  if (modeldata_check(modeldata,
                      required = c("residence_dist",
                                   "shedding_dist",
                                   "incubation_dist"),
                      throw_error = FALSE
  )) {
    inc_shed_dist <- with(
      modeldata, convolve(incubation_dist, rev(shedding_dist), type = "o")
      )
    total_dist <- with(
      modeldata, convolve(inc_shed_dist, rev(residence_dist), type = "o")
    )
    modeldata$.metainfo$partial_window <- which(cumsum(total_dist)>0.9)[1]-1
  }
  if (modeldata_check(
    modeldata,
    required = c(
      "measured_concentrations",
      "flow",
      ".metainfo$composite_window",
      ".metainfo$load_per_case"
    ),
    throw_error = FALSE
  )) {
    # crude descriptive estimate of cases at start of time series
    modeldata$.metainfo$initial_cases_crude <-
      with(
        modeldata,
        0.1 + measured_concentrations[1] *
          mean(flow[1:.metainfo$composite_window]) / .metainfo$load_per_case
      )
  }
  return(modeldata)
}

modeldata_descriptions <- function() {
  descriptions <- list(
    ".metainfo$composite_window" =
      "window length for composite samples in days",
    ".metainfo$length_seeding" =
      "length of seeding phase for infections",
    ".metainfo$length_I" =
      "number of days over which infections are modeled",
    ".metainfo$length_R" =
      "number of days over which Rt is modeled",
    ".metainfo$load_per_case" =
      "assumed overall load shed per individual",
    ".metainfo$initial_cases_crude" =
      "empirical estimate for #cases at start of time period"
  )
  return(descriptions)
}

modeldata_var_requirements <- function() {
  requirements <- list(
    ".metainfo$initial_cases_crude" = "flows_observe",
    ".metainfo$length_seeding" = "generation_dist_assume",
    ".metainfo$length_I" = c(
      "incubation_dist_assume",
      "shedding_dist_assume",
      "residence_dist_assume"
    ),
    ".metainfo$length_R" = c(
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
    "partitions",
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
    "noise",
    "horizon"
  )
  return(components)
}

all_parameters <- function(print = FALSE) {
  params <- as.data.frame(matrix(c(
    'measurement_noise_cv','nu_upsilon_a','Coefficient of variation (measurement noise)',1,identity,
    'dPCR_total_partitions','nu_upsilon_b_mu','Average total number of partitions in dPCR',1e4,identity,
    'dPCR_partition_variation','nu_upsilon_b_cv','Partition number variation in dPCR',1,identity,
    'dPCR_conversion_factor','nu_upsilon_c','Conversion factor in dPCR',1e-5,identity,
    'pre_replicate_cv','nu_psi','Coefficient of variation (pre-PCR noise)',1,identity,
    'load_variation_cv','nu_zeta','Individual-level coefficient of load variation',1,identity,
    'infection_overdispersion','I_xi','Overdispersion of infections',1,identity,
    'seeding_intercept','iota_log_seed_intercept','Initial number of infections',1,exp
  ), byrow = T, ncol = 5, dimnames = list(c(),c('short_name','raw_name','long_name',"scaling","transf"))))
  if (print) {
    return(c(paste(apply(params, 1, function(x) paste0("- `",x["short_name"],"` (",x["raw_name"],"): ",x["long_name"])), collapse = "\n")))
  } else {
    return(params)
  }
}
