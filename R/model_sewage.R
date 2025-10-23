#' Model the sewage process
#'
#' @description This module function is used to specify the components of the
#'  `sewage` module in `EpiSewer`.
#'
#' @description Each component can be specified using one or several helper
#'  functions (see available options below). See the documentation of the
#'  individual helper functions to adjust model priors and further settings.
#'
#' @param flows Daily flow volumes at the sampling site. The flow can change due
#'   to rainfall or industrial discharge, and directly influences pathogen
#'   concentrations in the wastewater. Modeling options:
#' `r component_functions_("flows")`
#' @param residence_dist Sewer residence time distribution for pathogen
#'   particles. By default, `EpiSewer` assumes that particles arrive at the
#'   sampling site within the day of shedding. However, for larger sewage
#'   systems, particles may travel longer than a day depending on where and
#'   when they were shed into the wastewater. Modeling options:
#' `r component_functions_("residence_dist")`
#'
#' @return A `modeldata` object containing the data and specifications of the
#'   `sewage` module.
#' @export
#' @family {module functions}
model_sewage <- function(
    flows = flows_observe(),
    residence_dist = residence_dist_assume()) {
  verify_is_modeldata(flows, "flows")
  verify_is_modeldata(residence_dist, "residence_dist")
  return(modeldata_combine(flows, residence_dist))
}

#' Assume a constant wastewater flow
#'
#' @description This option assumes a constant flow of wastewater at the
#'   sampling site. Can be used as an approximation if no regular flow
#'   measurements are available.
#'
#' @param flow_constant Daily wastewater flow volume, assumed to be the same for
#'   all days.
#'
#' @details The flow volume unit should be the same as for the concentration
#'   measurements, e.g. if concentrations are measured in gc/mL, then the flow
#'   should be in mL as well.
#'
#' @details Note that when the flow is unknown and therefore some arbitrary
#'   value (e.g. `flow_constant=1`) is assumed, this must also be accounted for
#'   in the assumed `load_per_case`. For this purpose, the function
#'   [suggest_load_per_case()] offers a `flow_constant` argument.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {flow models}
flows_assume <- function(
    flow_constant,
    modeldata = modeldata_init()) {
  modeldata <- tbc(
    "flow_data",
    {
      all_dates <-
        seq.Date(
          modeldata$.metainfo$T_start_date,
          modeldata$.metainfo$T_end_date + modeldata$.metainfo$forecast_horizon,
          by = "1 day"
        )
      modeldata$flow <- rep(flow_constant, length(all_dates))
    },
    required = c(
      ".metainfo$T_start_date",
      ".metainfo$T_end_date",
      ".metainfo$forecast_horizon"
      ),
    modeldata = modeldata
  )

  modeldata$.str$sewage[["flows"]] <- list(
    flows_assume = c()
  )

  return(modeldata)
}

#' Observe wastewater flows
#'
#' @description This option accounts for daily wastewater flow volumes measured
#'   at the sampling site. The flow can change due to rainfall or industrial
#'   discharge, and directly influences pathogen concentrations in the
#'   wastewater.
#'
#' @details The flow volume unit should be the same as for the concentration
#'   measurements, e.g. if concentrations are measured in gc/mL, then the flow
#'   should be in mL as well.
#'
#' @param flows A `data.frame` with each row representing one day. Must have at
#'   least a column with dates and a column with flow measurements.
#' @param date_col Name of the column containing the dates.
#' @param flow_col Name of the column containing the flows.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#' @family {flow models}
flows_observe <-
  function(flows = NULL,
           date_col = "date",
           flow_col = "flow",
           modeldata = modeldata_init()) {
    modeldata <- tbp("flows_observe",
      {
        required_data_cols <- c(date_col, flow_col)
        if (!all(required_data_cols %in% names(flows))) {
          cli::cli_abort(
            c(paste(
              "The following columns must be present",
              "in the provided flow `data.frame`:",
              paste(required_data_cols, collapse = ", ")
            ),
            paste("Please adjust the `data.frame` or specify the right column",
                  "names via the `_col` arguments of this function."))
          )
        }

        flows = as.data.table(flows)[, .SD, .SDcols = required_data_cols]
        flows <- setnames(flows, old = c(date_col, flow_col), new = c("date", "flow"))
        modeldata$.metainfo$flows_cols <- list(
          date_col = date_col,
          flow_col = flow_col
        )

        flows <- check_date_column(flows, date_col)

        if (any(duplicated(flows, by = "date"))) {
          flows <- unique(flows, by = c("date", "flow"))
          if (any(duplicated(flows, by = "date"))) { # if still duplicates
            cli::cli_abort(
              "Flow data is ambigious, duplicate dates with different values found."
              )
          }
        }

        modeldata <- tbc(
          "flow_data",
          {
            all_dates <-
              seq.Date(
                modeldata$.metainfo$T_start_date,
                modeldata$.metainfo$T_end_date + modeldata$.metainfo$forecast_horizon,
                by = "1 day"
              )
            missing_flow_dates <-
              lubridate::as_date(
                setdiff(all_dates, flows[!is.na(flow), date])
              )
            if (length(missing_flow_dates) > 0) {
              if (any(missing_flow_dates %in% modeldata$.metainfo$measured_dates)) {
              cli::cli_abort(paste(
                "Missing flow values for the following sampled dates:",
                paste(lubridate::as_date(intersect(
                  missing_flow_dates,
                  modeldata$.metainfo$measured_dates
                  )), collapse = ", ")
              ))
              } else {
                n_estimate = sum(missing_flow_dates <= modeldata$.metainfo$T_end_date)
                n_forecast = sum(missing_flow_dates > modeldata$.metainfo$T_end_date)
                if (n_forecast == 0) {
                  text_dates <- paste(n_estimate, "unsampled dates")
                } else if (n_estimate == 0) {
                  text_dates <- paste(n_forecast, "forecasting dates")
                } else {
                  text_dates <- paste(
                    n_estimate, "unsampled and", n_forecast, "forecasting dates"
                    )
                }
                cli::cli_inform(c(
                  "!" = paste(
                    "EpiSewer will impute missing flow data for", text_dates,
                    "using the median flow value."),
                  "*" = paste(
                     "This has no impact on model fitting and Rt estimation -",
                    "but means that absolute concentrations on these",
                    "unobserved dates may not be accurately predicted."),
                  "*" = paste(
                    "To disable this warning, please impute the missing",
                    "flow data manually before running EpiSewer.")
                ))

                flows <- flows[
                  CJ(date = all_dates, unique=TRUE),
                  on=.(date)
                ]
                median_flow <- median(flows$flow, na.rm = TRUE)
                setnafill(flows, fill = median_flow, cols = "flow")
              }
            }
            if (any(flows$flow==0)) {
              cli::cli_abort("Flow data must not contain zero flows.")
            }
            flows <-
              flows[
                date >= modeldata$.metainfo$T_start_date &
                date <= modeldata$.metainfo$T_end_date + modeldata$.metainfo$forecast_horizon,
                ]
            flows <- setorderv(flows, cols = "date")
            modeldata$flow <- flows$flow
          },
          required = c(
            ".metainfo$T_start_date",
            ".metainfo$T_end_date",
            ".metainfo$measured_dates",
            ".metainfo$forecast_horizon"
            ),
          modeldata = modeldata
        )

        return(modeldata)
      },
      required_data = "flows",
      modeldata = modeldata
    )

    modeldata$.str$sewage[["flows"]] <- list(
      flows_observe = c()
    )

    return(modeldata)
  }

#' Assume a sewer residence time distribution
#'
#' @description This option assumes a fixed residence time distribution for
#'   pathogen particles. By default, `EpiSewer` assumes that particles arrive at
#'   the sampling site within the day of shedding. However, for larger sewage
#'   systems, particles may travel longer than a day depending on where and
#'   when they were shed into the wastewater.
#'
#' @param residence_dist A numeric vector representing a discrete residence time
#'   distribution, with elements describing the share of load that takes 0 days,
#'   1 day, 2 days, and so on to arrive at the sampling site.
#'
#' @inheritParams template_model_helpers
#' @inherit modeldata_init return
#' @export
#'
#' @examples
#' # Particles arrive within the same day
#' residence_dist_assume(residence_dist = c(1))
#'
#' # Particles always arrive after one day
#' residence_dist_assume(residence_dist = c(0, 1))
#'
#' # 1/4 of particles only arrives after one day
#' residence_dist_assume(residence_dist = c(0.75, 0.25))
residence_dist_assume <-
  function(residence_dist = NULL, modeldata = modeldata_init()) {
    modeldata <- tbp("residence_dist_assume",
      {
        residence_dist <- check_dist(residence_dist, "residence time distribution")
        modeldata$D <- length(residence_dist) - 1
        modeldata$residence_dist <- residence_dist
        return(modeldata)
      },
      required_assumptions = "residence_dist",
      modeldata = modeldata
    )

    modeldata$.str$sewage[["residence_dist"]] <- list(
      residence_dist_assume = c()
    )

    return(modeldata)
  }
