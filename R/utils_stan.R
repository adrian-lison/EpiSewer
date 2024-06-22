#' Load stan model
#'
#' @description Loads a specific stan model for model fitting in `EpiSewer`.
#'
#' @param model_metainfo A `list` containing meta information about the
#'   specified model to be fitted. The required stan model is automatically
#'   inferred from the `model_metainfo`.
#' @param model_filename File name of a specific stan model to load. This is an
#'   alternative to supplying `model_metainfo`.
#' @param model_folder Path to the folder containing the stan models for
#'   `EpiSewer`.
#' @param profile Should profiling be run during model fitting? Default is
#'   `TRUE`. If `FALSE`, will remove all profiling statements from the model
#'   before fitting. This can decrease runtime in some cases.
#' @param threads Should multihreading be enabled? Default is `FALSE`, as
#'   `EpiSewer` currently does not support within-chain parallelism.
#' @param force_recompile If `FALSE` (default), the model is only recompiled if
#'   changes to the model code are detected. However, as the change detection is
#'   not fully reliable, it is sometimes necessary to force recompilation after
#'   having made changes to the stan code.
#' @param package Name of the package in which to search for stan files
#'   (defaults to EpiSewer). If NULL, will search in the normal working
#'   directory.
#'
#' @return A `list` containing the model definition and a link to the compiled
#'   stan model.
#' @export
get_stan_model <- function(
    model_metainfo = NULL,
    model_filename = NULL,
    model_folder = "stan",
    profile = TRUE,
    threads = FALSE,
    force_recompile = FALSE,
    package = "EpiSewer") {
  model_stan <- list()

  if (is.null(model_filename)) {
    if (is.null(model_metainfo)) {
      cli::cli_abort(
        c(
          "Please either provide ",
          "`model_filename` and `model_folder` (i.e. path to a stan model) or",
          "`model_metainfo` (to infer a suitable stan model)."
        )
      )
    }
    if (model_metainfo$R_estimate_approach == "splines") {
      model_filename <- "wastewater_Re_splines.stan"
    } else if (model_metainfo$R_estimate_approach == "ets") {
      model_filename <- "wastewater_Re.stan"
    } else if (model_metainfo$R_estimate_approach == "rw") {
      model_filename <- "wastewater_Re.stan"
    } else if (model_metainfo$R_estimate_approach == "fast") {
      model_filename <- "wastewater_Re_fast.stan"
    } else {
      cli::cli_abort(
        paste(
          "No suitable model available,",
          "please supply a stan model yourself",
          "in `model_stan_opts`, or open an issue."
        )
      )
    }
  }

  model_stan[["model_filename"]] <- model_filename
  model_stan[["model_folder"]] <- model_folder
  model_stan[["force_recompile"]] <- force_recompile
  model_stan[["threads"]] <- threads
  model_stan[["profile"]] <- profile
  model_stan[["package"]] <- package

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
#' @param model Details about the model file to be used, see
#'   [model_stan_opts()]. This also makes it possible for users to supply
#'   customized models.
#'
#' @return A `list` with the model fitting options.
#' @export
set_fit_opts <- function(sampler = sampler_stan_mcmc(), fitted = TRUE, model = model_stan_opts()) {
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
    iter_warmup = 500,
    iter_sampling = 500,
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

#' Specify details of the stan model
#'
#' @param model_filename Name of the stan model file. If NULL (default), the
#'   name of the model file is automatically inferred based on the model
#'   specification.
#' @param model_folder (Relative) path to the folder containing the stan model.
#' @param profile Should profiling be run during model fitting? Default is
#'   `TRUE`. If `FALSE`, will remove all profiling statements from the model
#'   before fitting. This can decrease runtime in some cases.
#' @param threads Should multihreading be enabled? Default is `FALSE`, as
#'   `EpiSewer` currently does not support within-chain parallelism.
#' @param force_recompile Should recompilation be forced before model fitting?
#'   If `FALSE` (default), the model is only recompiled when changes to the
#'   model code are detected. However, as the change detection is not fully
#'   reliable, it is sometimes necessary to force recompilation after having
#'   made changes to the stan code.
#' @param package Name of the package in which to search for stan files
#'   (defaults to EpiSewer). If NULL, will search in the normal working
#'   directory, allowing users to provide their own customized model files.
#'
#' @return
#' @export
model_stan_opts <- function(model_filename = NULL, model_folder = "stan",
                            profile = TRUE, threads = FALSE,
                            force_recompile = FALSE, package = "EpiSewer") {
  opts <- as.list(environment())
  return(opts)
}

update_compiled_stanmodel <- function(model_stan, force_recompile = FALSE) {
  n_models <- length(model_stan$model_filename)

  model_path_list <- lapply(1:n_models, function(i) {
    if (is.null(model_stan$package)) {
      file.path(
        model_stan$model_folder,
        model_stan$model_filename[[i]]
      )
    } else {
      system.file(
        model_stan$model_folder,
        model_stan$model_filename[[i]],
        package = model_stan$package
      )
    }

  })
  include_paths_list <- lapply(1:n_models, function(i) {
    if (is.null(model_stan$package)) {
      file.path(
        model_stan$model_folder
      )
    } else {
      system.file(
        model_stan$model_folder,
        package = model_stan$package
      )
    }
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
    (cmdstanr::cmdstan_model(model_path_list[[i]],
      include_paths = include_paths_list[[i]],
      dir = dirname(model_path_list[[i]]),
      cpp_options = cpp_options,
      force_recompile = force_recompile
    ))
  })

  model_stan[["load_model"]] <- lapply(1:n_models, function(i) {
    function() {
      return(stanmodel_list[[i]])
    }
  })

  return(model_stan)
}

#' Compile EpiSewer models
#'
#' @description The stan models used by EpiSewer need to be compiled for your
#'   device. This is only necessary once, after installing or updating the
#'   package. This function compiles all models from the package.
#'
#' @param force_recompile If TRUE, then models will be recompiled even if they
#'   have already been successfully compiled.
#' @param verbose If TRUE, warnings and detailed errors from the compilation are
#'   printed. This can help to diagnose compilation issues.
#'
#' @details If one or several models are not successfully compiled, please
#'   ensure that `cmdstan` is properly set up and try updating it to a newer
#'   version using [cmdstanr::install_cmdstan()]. If the problem persists,
#'   please run [EpiSewer::sewer_compile(verbose = TRUE)] and post the output in
#'   a new issue on GitHub, along with your [cmdstanr::cmdstan_version()].
#'
#' @export
sewer_compile <- function(force_recompile = FALSE, verbose = FALSE) {
  all_models <- c(
    "wastewater_Re.stan",
    "wastewater_Re_splines.stan",
    "wastewater_Re_fast.stan"
    )
  comp_success <- NULL

  for (i in 1:length(all_models)) {
    model_name <- all_models[i]
    cat("\r                                                      \r", sep = "");
    cat(sprintf("\r| Compiling model %d/%d ", i, length(all_models)), sep = "")
    #flush.console()
    success <- tryCatch(
      {
        if (verbose == FALSE) {
          suppressMessages(get_stan_model(
            model_filename = model_name,
            force_recompile = force_recompile
          ))
        } else {
          get_stan_model(
            model_filename = model_name,
            force_recompile = force_recompile
          )
        }
        TRUE
      },
      error = function(e) { FALSE }
    )
    cat(paste(sprintf(
      "\r| Compiling model %d/%d",
      i, length(all_models)), ifelse(success, "(success) ", "(failed) ")
      ), sep = "")
    Sys.sleep(0.5)
    comp_success <- c(comp_success, success)
  }
  cat("\r                                                        \r", sep = "");
  if(!all(comp_success)) {
    cli::cli_warn(
      paste(
        "The following models could not be compiled:",
        paste(all_models[!comp_success], collapse = ", ")
      )
    )
  } else {

    cli::cli_alert("All models compiled successfully.")
  }
  return(invisible(NULL))
}

#' Computes a checksum of a stanmodel, i.e. of the code from its main
#' stan file and all included stan files
get_checksum_model <- function(stanmodel, only_functions = FALSE) {
  if ("try-error" %in% class(stanmodel)) {
    cli::cli_abort("There was an error compiling the stan model.")
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
  final_checksum <- digest::digest(
    paste0(all_digests, collapse = ""),
    algo = "md5",
    serialize = FALSE
  )
  return(final_checksum)
}


