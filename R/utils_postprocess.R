extract_dims_df <- function(df, col) {
  max_ndims <- max(stringr::str_count(df[[col]], ",") + 1)
  for (i in max_ndims:1) {
    df <- cbind(as.integer(stringr::str_extract(
      df[[col]],
      paste0(
        "(?<=\\[",
        paste0(rep("\\d{1,999}", i - 1), collapse = ","),
        ifelse(i > 1, ",)", ")"),
        "\\d+(?=[,\\d]*\\])"
      )
    )), df)
  }
  df <- cbind(stringr::str_remove(df[[col]], "\\[.+\\]"), df)
  names(df)[1:(max_ndims + 1)] <-
    c("var_name", paste0("dim", 1:max_ndims))
  return(df)
}

map_dates_1d_df <- function(fit_summary, date_mapping) {
  fit_summary <- extract_dims_df(fit_summary, "variable")
  fit_summary <- cbind(date = date_mapping[fit_summary$dim1], fit_summary)
  fit_summary <- subset(fit_summary, select = -c(var_name, dim1, variable))
  setDT(fit_summary)
  return(fit_summary)
}

# T_shift: how much does the variable lead or lag the time from 1:T?
get_summary_1d_date <- function(fit, var, T_shift, .metainfo,
                                var_forecast = NULL,
                                intervals = c(0.5, 0.95)) {
  vars <- list(
    var = var
  )
  date_mappings <- list(
    date_mapping = seq.Date(
      .metainfo$T_start_date - T_shift, .metainfo$T_end_date,
      by = "1 day"
    ))
  types <- list(
    type = factor("estimate", levels = c("estimate", "forecast"))
  )

  if (.metainfo$forecast_horizon > 0 && !is.null(var_forecast)) {
    vars$var_forecast <- var_forecast
    date_mappings$date_mapping_forecast <- seq.Date(
      .metainfo$T_end_date + 1,
      .metainfo$T_end_date + .metainfo$forecast_horizon,
      by = "1 day"
    )
    types$type_forecast <- factor("forecast", levels = c("estimate", "forecast"))
  }

  var_summary <- rbindlist(mapply(function(variable, dates, type) {
    if (is.null(variable)) {return(NULL)}
    var_summary <- fit$summary(
      variable,
      mean = function(x) mean(x, na.rm = TRUE),
      median = function(x) median(x, na.rm = TRUE),
      function(x) {
        setNames(
          quantile(x, (1 - rev(intervals)) / 2, na.rm = TRUE), paste0("lower_", rev(intervals))
        )
      },
      function(x) {
        setNames(
          quantile(x, (1 + intervals) / 2, na.rm = TRUE), paste0("upper_", intervals)
        )
      }
    )
    var_summary <- map_dates_1d_df(var_summary, dates)
    var_summary[, type := type]
  },
  variable = vars,
  dates = date_mappings,
  type = types,
  SIMPLIFY = FALSE
  ))
  setDT(var_summary)
  return(var_summary)
}

# T_shift: how much does the variable lead or lag the time from 1:T?
get_summary_1d_date_log <- function(fit, var, T_shift, .metainfo,
                                    var_forecast = NULL,
                                    intervals = c(0.5, 0.95)) {
  vars <- list(
    var = var
    )
  date_mappings <- list(
    date_mapping = seq.Date(
    .metainfo$T_start_date - T_shift, .metainfo$T_end_date,
    by = "1 day"
  ))
  types <- list(
    type = factor("estimate", levels = c("estimate", "forecast"))
    )

  if (.metainfo$forecast_horizon > 0) {
    vars$var_forecast <- var_forecast
    date_mappings$date_mapping_forecast <- seq.Date(
      .metainfo$T_end_date + 1,
      .metainfo$T_end_date + .metainfo$forecast_horizon,
      by = "1 day"
    )
    types$type_forecast <- factor("forecast", levels = c("estimate", "forecast"))
  }

  var_summary <- rbindlist(mapply(function(variable, dates, type) {
    if (is.null(variable)) {return(NULL)}
    var_summary <- fit$summary(
      variable,
      mean = function(x) mean(exp(x), na.rm = TRUE),
      median = function(x) median(exp(x), na.rm = TRUE),
      function(x) {
        setNames(
          quantile(exp(x), (1 - rev(intervals)) / 2, na.rm = TRUE), paste0("lower_", rev(intervals))
        )
      },
      function(x) {
        setNames(
          quantile(exp(x), (1 + intervals) / 2, na.rm = TRUE), paste0("upper_", intervals)
        )
      }
    )
    var_summary <- map_dates_1d_df(var_summary, dates)
    var_summary[, type := type]
  },
  variable = vars,
  dates = date_mappings,
  type = types,
  SIMPLIFY = FALSE
  ))
  setDT(var_summary)
  return(var_summary)
}

# T_shift: how much does the variable lead or lag the time from 1:T?
get_diagnostics_1d_date <- function(fit, var, T_shift, .metainfo,
                                var_forecast = NULL,
                                intervals = c(0.5, 0.95)) {
  vars <- list(
    var = var
  )
  date_mappings <- list(
    date_mapping = seq.Date(
      .metainfo$T_start_date - T_shift, .metainfo$T_end_date,
      by = "1 day"
    ))
  types <- list(
    type = factor("estimate", levels = c("estimate", "forecast"))
  )

  if (.metainfo$forecast_horizon > 0 && !is.null(var_forecast)) {
    vars$var_forecast <- var_forecast
    date_mappings$date_mapping_forecast <- seq.Date(
      .metainfo$T_end_date + 1,
      .metainfo$T_end_date + .metainfo$forecast_horizon,
      by = "1 day"
    )
    types$type_forecast <- factor("forecast", levels = c("estimate", "forecast"))
  }

  var_summary <- rbindlist(mapply(function(variable, dates, type) {
    if (is.null(variable)) {return(NULL)}
    var_summary <- fit$summary(
      variable,
      posterior::default_convergence_measures()
    )
    var_summary <- map_dates_1d_df(var_summary, dates)
    var_summary[, type := type]
  },
  variable = vars,
  dates = date_mappings,
  type = types,
  SIMPLIFY = FALSE
  ))
  setDT(var_summary)
  return(var_summary)
}

get_draws_1d_date <- function(fit, variable, ndraws = NULL, draw_ids = NULL) {
  fit_draws <- fit$draws(variable, format = "df")
  setDT(fit_draws)
  if (!is.null(ndraws)) {
    draw_ids <- sample.int(max(fit_draws$.draw), ndraws, replace = FALSE)
    fit_draws <- fit_draws[.draw %in% draw_ids, ]
  } else if (!is.null(draw_ids)) {
    fit_draws <- fit_draws[.draw %in% draw_ids, ]
  }
  fit_draws <- melt(fit_draws,
    id.vars = c(".chain", ".iteration", ".draw"),
    value.name = variable
  )
  fit_draws[, date := stringr::str_extract(variable, "(?<=\\[)\\d+(?=\\])")]
  return(fit_draws[])
}

get_latent_trajectories <- function(fit, var, T_shift, .metainfo,
                                    var_forecast = NULL,
                                    ndraws = NULL, draw_ids = NULL) {

  fit_draws <- get_draws_1d_date(fit, var, ndraws = ndraws, draw_ids = draw_ids)
  date_mapping <- seq.Date(
    .metainfo$T_start_date - T_shift, .metainfo$T_end_date,
    by = "1 day"
  )
  fit_draws[, date := date_mapping[as.integer(date)]]
  fit_draws[, c(".iteration", "variable") := NULL]
  fit_draws[, type := "estimate"]

  if (.metainfo$forecast_horizon > 0) {
    fit_draws_forecast <- get_draws_1d_date(
      fit, var_forecast, draw_ids = unique(fit_draws$.draw)
    )
    setnames(fit_draws_forecast, var_forecast, var)
    date_mapping_forecast <- seq.Date(
      .metainfo$T_end_date + 1,
      .metainfo$T_end_date + .metainfo$forecast_horizon,
      by = "1 day"
    )
    fit_draws_forecast[, date := date_mapping_forecast[as.integer(date)]]
    fit_draws_forecast[, c(".iteration", "variable") := NULL]
    fit_draws_forecast[, type := "forecast"]
  } else {
    fit_draws_forecast <- NULL
  }

  return(rbindlist(list(fit_draws, fit_draws_forecast))[])
}

get_summary_vector <- function(fit, var, varnames = NULL,
                               intervals = c(0.5, 0.95)) {
  fsummary <- fit$summary(var,
    mean = function(x) mean(x, na.rm = TRUE),
    median = function(x) median(x, na.rm = TRUE),
    function(x) {
      setNames(
        quantile(x, (1 - rev(intervals)) / 2, na.rm = TRUE), paste0("lower_", rev(intervals))
      )
    },
    function(x) {
      setNames(
        quantile(x, (1 + intervals) / 2, na.rm = TRUE), paste0("upper_", intervals)
      )
    }
  )
  if (!is.null(varnames)) {
    if (nrow(fsummary) != length(varnames)) {
      cli::cli_abort("Mismatch between model var and varnames", .internal = TRUE)
    }
    fsummary$variable <- forcats::fct_inorder(varnames, ordered = TRUE)
  }
  setDT(fsummary)
  return(fsummary)
}

get_summary_vector_log <- function(fit, var, varnames = NULL,
                                   intervals = c(0.5, 0.95)) {
  fsummary <- fit$summary(var,
    mean = function(x) mean(exp(x), na.rm = TRUE),
    median = function(x) median(exp(x), na.rm = TRUE),
    function(x) {
      setNames(
        quantile(exp(x), (1 - rev(intervals)) / 2, na.rm = TRUE), paste0("lower_", rev(intervals))
      )
    },
    function(x) {
      setNames(
        quantile(exp(x), (1 + intervals) / 2, na.rm = TRUE), paste0("upper_", intervals)
      )
    }
  )
  if (!is.null(varnames)) {
    if (nrow(fsummary) != length(varnames)) {
      cli::cli_abort("Mismatch between model var and varnames", .internal = TRUE)
    }
    fsummary$variable <- forcats::fct_inorder(varnames, ordered = TRUE)
  }
  setDT(fsummary)
  return(fsummary)
}

# Abort with an informative message if a result has no valid summary, e.g.
# because summarizing the fit failed (stored as a "try-error").
assert_valid_summary <- function(result) {
  summary <- result$summary
  if (
    is.null(summary) || inherits(summary, "try-error") || !is.list(summary)
  ) {
    msg <- "The result has no valid summary."
    if (inherits(summary, "try-error")) {
      msg <- c(msg, "x" = as.character(summary))
    } else if (!is.null(result$errors)) {
      msg <- c(msg, "x" = "The model fit returned errors.")
    }
    cli::cli_abort(msg, call = NULL)
  }
  invisible(TRUE)
}

combine_summaries <- function(result_list, summary_name) {
  lapply(result_list, assert_valid_summary)
  summary_list <- lapply(result_list, function(x) x$summary[[summary_name]])
  combined <- vctrs::vec_rbind(!!!summary_list, .names_to = "model")
  combined$model <- forcats::fct_inorder(
    as.character(combined$model),
    ordered = TRUE
  )
  return(combined)
}

combine_summary_list <- function(summary_list) {
  combined <- vctrs::vec_rbind(!!!summary_list, .names_to = "model")
  combined$model <- forcats::fct_inorder(
    as.character(combined$model),
    ordered = TRUE
  )
  return(combined)
}

combine_samples <- function(
    result_list, variable_name, draws = FALSE, ndraws = NULL) {
  lapply(result_list, assert_valid_summary)
  summary_list <- lapply(result_list, function(x) {
    res <- x$summary$samples[
      , .SD, .SDcols = c("date", "type", "seeding", ".draw", ".chain", variable_name)
      ]
    if (!is.null(ndraws)) {
      draw_ids <- unique(res$.draw)
      if (ndraws > length(draw_ids)) {
        cli::cli_abort(paste0(
          "The maximum number of available draws is ", length(draw_ids), "."
        ))
      }
      selected_ids <- sample(draw_ids, ndraws)
      res <- res[res$.draw %in% selected_ids, ]
    }
    return(res)
  })
  combined <- vctrs::vec_rbind(!!!summary_list, .names_to = "model")
  combined$model <- forcats::fct_inorder(
    as.character(combined$model),
    ordered = TRUE
  )
  setDT(combined)
  return(combined)
}

summarize_samples_dt <- function(dt, index_cols, variable, intervals, cols_end) {
  summary_dt <- dt[, .SD, .SDcols=c(index_cols, variable)]
  summary_functions <- c(
    list(
      function(x) mean(x, na.rm = TRUE),
      function(x) median(x, na.rm = TRUE)
    ),
    lapply(
      rev(intervals),
      function(i) function(x) quantile(x, (1 - i) / 2, na.rm = TRUE)
    ),
    lapply(
      intervals,
      function(i) function(x) quantile(x, (1 + i) / 2, na.rm = TRUE)
    )
  )
  summary_names <- c(
    "mean", "median",
    paste0("lower_", rev(intervals)),
    paste0("upper_", intervals)
  )
  result <- summary_dt[,
                    setNames(lapply(summary_functions, function(f) f(.SD[[1]])), summary_names),
                    by=c("date", "type", "seeding"), .SDcols=c(variable)
  ]
  setcolorder(result, c(1, (max(cols_end)+1):ncol(result), cols_end))
  return(result)
}

# Remove one level of list nesting, keeping the names of the inner elements.
# Only bare lists are spliced; other elements (including data.frames and
# atomic vectors) are kept as they are. Replacement for the deprecated
# rlang::flatten().
flatten_one_level <- function(x) {
  is_bare <- vapply(x, function(e) is.list(e) && !is.object(e), logical(1))
  pieces <- Map(function(e, spliced) if (spliced) e else list(e), x, is_bare)
  inner_names <- Map(function(e, nm, spliced) {
    if (spliced) names(e) else nm
  }, x, names(x), is_bare)
  out <- unlist(unname(pieces), recursive = FALSE)
  if (is.null(out)) out <- list()
  names(out) <- unlist(inner_names)
  out
}
