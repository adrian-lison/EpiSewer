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

get_summary_1d_date <- function(fit, var, T_shift, meta_info, intervals = c(0.95, 0.5)) {
  # T_shift: how much does the variable lead or lag the time from 1:T?
  date_mapping <- seq.Date(
    meta_info$T_start_date - T_shift, meta_info$T_end_date,
    by = "1 day"
  )
  var_summary <- map_dates_1d_df(fit$summary(var,
    mean = mean,
    median = median,
    function(x) setNames(quantile(x, (1-intervals)/2), paste0("lower_", intervals)),
    function(x) setNames(quantile(x, (1+intervals)/2), paste0("upper_", intervals))
  ), date_mapping)
  setDT(var_summary)
  return(var_summary)
}

get_summary_1d_date_log <- function(fit, var, T_shift, meta_info, intervals = c(0.95, 0.5)) {
  # T_shift: how much does the variable lead or lag the time from 1:T?
  date_mapping <- seq.Date(
    meta_info$T_start_date - T_shift, meta_info$T_end_date,
    by = "1 day"
  )
  var_summary <- map_dates_1d_df(fit$summary(var,
    mean = function(x) mean(exp(x)),
    median = function(x) median(exp(x)),
    function(x) setNames(quantile(exp(x), (1-intervals)/2), paste0("lower_", intervals)),
    function(x) setNames(quantile(exp(x), (1+intervals)/2), paste0("upper_", intervals))
  ), date_mapping)
  setDT(var_summary)
  return(var_summary)
}

get_draws_1d_date <- function(fit, variable, ndraws = NULL) {
  fit_draws <- fit$draws(variable, format = "df")
  setDT(fit_draws)
  if (!is.null(ndraws)) {
    draw_ids <- sample.int(max(fit_draws$.draw), ndraws, replace = FALSE)
    fit_draws <- fit_draws[ .draw %in% draw_ids,]
  }
  fit_draws <- melt(fit_draws,
                    id.vars = c(".chain", ".iteration", ".draw"),
                    value.name = variable
  )
  fit_draws[, date := stringr::str_extract(variable, "(?<=\\[)\\d+(?=\\])")]
  return(fit_draws[])
}

get_summary_vector <- function(fit, var, varnames = NULL, intervals = c(0.95, 0.5)) {
  fsummary <- fit$summary(var,
    mean = mean,
    median = median,
    function(x) setNames(quantile(x, (1-intervals)/2), paste0("lower_", intervals)),
    function(x) setNames(quantile(x, (1+intervals)/2), paste0("upper_", intervals))
  )
  if (!is.null(varnames)) {
    if (nrow(fsummary) != length(varnames)) {
      abort("Mismatch between model var and varnames", .internal = TRUE)
    }
    fsummary$variable <- forcats::fct_inorder(varnames, ordered = TRUE)
  }
  setDT(fsummary)
  return(fsummary)
}

get_summary_vector_log <- function(fit, var, varnames = NULL, intervals = c(0.95, 0.5)) {
  fsummary <- fit$summary(var,
    mean = function(x) mean(exp(x)),
    median = function(x) median(exp(x)),
    function(x) setNames(quantile(exp(x), (1-intervals)/2), paste0("lower_", intervals)),
    function(x) setNames(quantile(exp(x), (1+intervals)/2), paste0("upper_", intervals))
  )
  if (!is.null(varnames)) {
    if (nrow(fsummary) != length(varnames)) {
      abort("Mismatch between model var and varnames", .internal = TRUE)
    }
    fsummary$variable <- forcats::fct_inorder(varnames, ordered = TRUE)
  }
  setDT(fsummary)
  return(fsummary)
}

combine_summaries <- function(result_list, summary_name) {
  summary_list <- lapply(result_list, function(x) x$summary[[summary_name]])
  combined <- vctrs::vec_rbind(!!!summary_list, .names_to = "model")
  combined$model <- forcats::fct_inorder(
    as.character(combined$model),
    ordered = TRUE
  )
  return(combined)
}

combine_samples <- function(
    result_list, summary_name, draws = FALSE, ndraws = NULL) {
  summary_list <- lapply(result_list, function(x) {
    res <- x$summary[[summary_name]]
    if (!is.null(ndraws)) {
      draw_ids <- unique(res$.draw)
      if (ndraws > length(draw_ids)) {
        abort(paste(
          "The maximum number of available draws is", length(draw_ids)
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
