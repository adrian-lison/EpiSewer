modeldata_descriptions <- function() {
  descriptions <- list(
    ".metainfo$composite_window" =
      "window length for composite samples in days",
    ".metainfo$length_seeding" =
      "length of seeding phase for infections",
    ".metainfo$length_I" =
      "number of days over which infections are modeled",
    ".metainfo$length_R" =
      "number of days over which Rt is estimated (including extended seeding)",
    ".metainfo$length_R_modeled" =
      "number of days over which Rt is explicitly modeled",
    ".metainfo$load_per_case" =
      "assumed overall load shed per individual",
    ".metainfo$initial_cases_crude" =
      "empirical estimate for #cases at start of time period"
  )
  return(descriptions)
}

modeldata_defaults <- function() {
  defaults <- list(
    "K" = 0,
    "X" = numeric(0),
    "eta_prior" = numeric(0),
    ".init$eta" = numeric(0)
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
    "outliers",
    "sample_effects",
    "noise",
    "horizon",
    "damping"
  )
  return(components)
}

all_parameters <- function(print = FALSE) {
  params <- as.data.frame(matrix(c(
    'measurement_noise_cv','nu_upsilon_a','Coefficient of variation (measurement noise)',1,identity,
    'dPCR_maximum_partitions','max_partitions','Maximum number of partitions in dPCR',1e4,identity,
    'dPCR_partition_loss_mean','partition_loss_mu','Mean relative partition loss in dPCR',function(x) x$job$data$partition_loss_max,function(x) plogis(x),
    'dPCR_partition_loss_variation','partition_loss_sigma','Partition number variation in dPCR',1,identity,
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
