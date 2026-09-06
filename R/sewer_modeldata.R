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
