#' Load stan model
#'
#' @description Loads a specific stan model for model fitting in `EpiSewer`.
#'
#' @param modeldata A `modeldata` object containing data and specifications of
#'   the model to be fitted. The required stan model is automatically inferred
#'   from the `modeldata`.
#' @param model_filename File name of a specific stan model to load. This is an
#'   alternative to supplying `modeldata`.
#' @param model_folder Path to the folder containing the stan models for
#'   `EpiSewer`.
#' @param profile Should profiling be run during model fitting? Default is
#'   `TRUE`. Disabling profiling can decrease runtime in some cases.
#' @param threads Should multihreading be enabled? Default is `FALSE`, as
#'   `EpiSewer` currently does not support within-chain parallelism.
#' @param force_recompile If `FALSE` (default), the model is only recompiled if
#'   changes to the model code are detected. However, as the change detection is
#'   not fully reliable, it is sometimes necessary to force recompilation after
#'   having made some changes to the stan code.
#'
#' @return A `list` containing the model definition and a link to the compiled
#'   stan model.
#' @export
get_stan_model <- function(
    modeldata = NULL,
    model_filename = NULL,
    model_folder = "stan",
    profile = TRUE,
    threads = FALSE,
    force_recompile = FALSE) {
  model_stan <- list()

  if (is.null(model_filename)) {
    if (is.null(modeldata)) {
      rlang::abort(
        c(
          "Please either provide ",
          "`model_filename` and `model_folder` (i.e. path to a stan model)",
          "or `modeldata` (so that a suitable stan model can be inferred)."
        )
      )
    }
    if (modeldata$.metainfo$R_estimate_approach == "splines") {
      model_filename <- "wastewater_Re_splines.stan"
    } else if (modeldata$.metainfo$R_estimate_approach == "ets") {
      model_filename <- "wastewater_Re.stan"
    } else if (modeldata$.metainfo$R_estimate_approach == "rw") {
      model_filename <- "wastewater_Re.stan"
    } else {
      rlang::abort(
        paste(
          "No suitable model available,",
          "please supply a stan model yourself",
          "using `get_stan_model`, or open an issue."
        )
      )
    }
  }

  model_stan[["model_filename"]] <- model_filename
  model_stan[["model_folder"]] <- model_folder
  model_stan[["force_recompile"]] <- force_recompile
  model_stan[["threads"]] <- threads
  model_stan[["profile"]] <- profile

  model_stan <- update_compiled_stanmodel(model_stan, force_recompile)

  return(model_stan)
}

#' Configure the model fitting
#'
#' @description Sets model fitting options in `EpiSewer`.
#'
#' @param sampler Which type of sampler should be used for model fitting?
#'   Currently, only [sampler_stan_mcmc()] is supported.
#' @param fitted If `TRUE` (default), the fitted model object is also returned,
#'   not only summaries of the model fitting. Note that this makes the results
#'   object substantially larger.
#'
#' @return A `list` with the model fitting options.
#' @export
set_fit_opts <- function(sampler = sampler_stan_mcmc(), fitted = TRUE) {
  opts <- as.list(environment())
  return(opts)
}

#' Use the stan MCMC sampler
#'
#' @description This option will use stan's NUTS sampler via [cmdstanr] for
#'   Markov Chain Monte Carlo (MCMC) sampling of the `EpiSewer` model.
#'
#' @param ... Further arguments to pass to [cmdstanr].
#' @inheritParams cmdstanr::sample
#'
#' @return A `list` with settings for the MCMC sampler.
#' @export
sampler_stan_mcmc <- function(
    chains = 4,
    iter_warmup = 1000,
    iter_sampling = 1000,
    adapt_delta = 0.99,
    max_treedepth = 15,
    step_size = 0.01,
    parallel_chains = NULL,
    threads_per_chain = 1,
    seed = 0,
    refresh = 200,
    show_messages = TRUE,
    show_exceptions = FALSE,
    ...) {
  opts <- c(as.list(environment()), list(...))
  if (opts$threads_per_chain == 1) {
    opts$threads_per_chain <- NULL
  }
  return(opts)
}

update_compiled_stanmodel <- function(model_stan, force_recompile = FALSE) {
  n_models <- length(model_stan$model_filename)

  model_path_list <- lapply(1:n_models, function(i) {
    system.file(model_stan$model_folder,
                model_stan$model_filename[[i]], package = "EpiSewer")
  })
  include_paths_list <- lapply(1:n_models, function(i) {
    system.file(model_stan$model_folder, package = "EpiSewer")
  }) # identical

  if (!model_stan$profile) {
    stan_no_profiling_list <- lapply(1:n_models, function(i) {
      write_stan_files_no_profiling(
        model_path_list[[i]],
        include_paths_list[[i]]
      )
    })
    model_path_list <- lapply(1:n_models, function(i) {
      stan_no_profiling_list[[i]]$model
    })
    include_paths_list <- lapply(1:n_models, function(i) {
      stan_no_profiling_list[[i]]$include_paths
    })
  }

  cpp_options <- list()
  if (model_stan$threads) {
    cpp_options[["stan_threads"]] <- TRUE
  }

  stanmodel_list <- lapply(1:n_models, function(i) {
    try(cmdstanr::cmdstan_model(model_path_list[[i]],
      include_paths = include_paths_list[[i]],
      dir = dirname(model_path_list[[i]]),
      cpp_options = cpp_options,
      force_recompile = force_recompile
    ))
  })

  model_stan[["get_stan_model"]] <- lapply(1:n_models, function(i) {
    function() {
      return(stanmodel_list[[i]])
    }
  })

  return(model_stan)
}

#' Computes a deterministic hash of a stanmodel, i.e. of the code from its main
#' stan file and all included stan files
get_hash_model <- function(stanmodel, only_functions = FALSE) {
  if ("try-error" %in% class(stanmodel)) {
    rlang::abort("There was an error compiling the stan model.")
  }
  stanmodelcode <- paste(stanmodel$code(), collapse = " ")
  includes <- stringr::str_extract_all(
    stanmodelcode,
    "(?<=#include )functions/.*?\\.stan"
  )[[1]]
  include_files <- list.files(
    stanmodel$include_paths(),
    recursive = TRUE,
    full.names = TRUE
  )
  include_files <- include_files[
    sapply(
      include_files,
      function(x) any(stringr::str_detect(x, includes))
    )
  ]
  all_digests <- sapply(c(stanmodel$stan_file(), include_files), function(x) {
    digest::digest(file = x, algo = "md5")
  })
  final_hash <- digest::digest(
    paste0(all_digests, collapse = ""),
    algo = "md5",
    serialize = FALSE
  )
  return(final_hash)
}


