#' Load stan model
#'
#' @description Loads a specific stan model for model fitting in `EpiSewer`.
#'
#' @param model_filename File name of a specific stan model to load. If NULL,
#'   the default model will be loaded.
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
#' @param use_docker Use the model compiled in a docker container? Note that if
#'   TRUE, all other arguments are ignored and the precompiled model in the
#'   container will be used. Default is FALSE.
#'
#' @return A `list` containing the model definition and a link to the compiled
#'   stan model.
#' @export
#' @keywords internal
get_stan_model <- function(
    model_filename = NULL,
    model_folder = "stan",
    profile = TRUE,
    threads = FALSE,
    force_recompile = FALSE,
    package = "EpiSewer",
    use_docker = FALSE) {

  model_stan <- list()

  if (is.null(model_filename)) {
    model_stan[["model_filename"]] <- "EpiSewer_main.stan"
  } else {
    model_stan[["model_filename"]] <- model_filename
  }

  if (use_docker) {
    model_stan[["model_folder"]] <- "/opt/models"
    model_stan[["package"]] <- "EpiSewer-docker"
    model_filepath <- file.path(
      model_stan[["model_folder"]], model_stan[["model_filename"]]
    )
    model_stan[["load_model"]] <- list(function() {cmdstan_model(
      exe_file = tools::file_path_sans_ext(model_filepath),
      compile = FALSE
      )})
  } else {
    model_stan[["model_folder"]] <- model_folder
    model_stan[["force_recompile"]] <- force_recompile
    model_stan[["threads"]] <- threads
    model_stan[["profile"]] <- profile
    model_stan[["package"]] <- package
    model_stan <- update_compiled_stanmodel(model_stan, force_recompile)
  }

  return(model_stan)
}

#' Configure the model fitting
#'
#' @description Sets model fitting options in `EpiSewer`.
#'
#' @param sampler Which type of sampler should be used for model fitting?
#'   Currently supported are
#'    - [sampler_stan_mcmc()]: default, recommended for inference
#'    - [sampler_stan_pathfinder()]: potentially inaccurate, recommended for
#'    preview only
#' @param model Details about the model file to be used, see
#'   [model_stan_opts()]. This also makes it possible for users to supply
#'   customized models.
#'
#' @return A `list` with the model fitting options.
#' @export
#' @examples
#' ww_data <- ww_data_SARS_CoV_2_Zurich
#' ww_assumptions <- ww_assumptions_SARS_CoV_2_Zurich
#'
#' ww_fit_opts <- set_fit_opts(
#'   sampler = sampler_stan_mcmc(
#'     chains = 2,
#'     parallel_chains = 2,
#'     iter_warmup = 400,
#'     iter_sampling = 400,
#'     seed = 42 # ensures reproducibility
#'   )
#' )
#'
#' \dontrun{ww_result <- EpiSewer(
#'   data = ww_data,
#'   assumptions = ww_assumptions,
#'   fit_opts = ww_fit_opts
#' )}
set_fit_opts <- function(sampler = sampler_stan_mcmc(),
                         model = model_stan_opts()) {
  opts <- as.list(environment())
  return(opts)
}

#' Use the stan MCMC sampler
#'
#' @description This option uses stan's NUTS sampler via [cmdstanr] for Markov
#'   Chain Monte Carlo (MCMC) sampling of the `EpiSewer` model.
#'
#' @param init_pathfinder Should stan's pathfinder algorihtm be used to
#'   initialize the MCMC chains?
#' @param init_pathfinder_max_lbfgs_iters The maximum number of iterations for
#'   LBFGS during pathfinder initialization.
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
    init_pathfinder = TRUE,
    init_pathfinder_max_lbfgs_iters = NULL,
    ...) {
  opts <- c(as.list(environment()), list(...))
  if (opts$threads_per_chain == 1) {
    opts$threads_per_chain <- NULL
  }
  class(opts) <- "mcmc"
  return(opts)
}

#' Use stan's pathfinder variational inference algorithm
#'
#' @description This option uses stan's pathfinder algorithm via [cmdstanr]
#'   for variational inference. Sampling is fast but yields potentially
#'   inaccurate results with unreliable uncertainty information. Currently, this
#'   sampler is recommended only for preview purposes, not for real-world
#'   inference.
#'
#' @param ... Further arguments to pass to [cmdstanr].
#' @inheritParams cmdstanr::pathfinder
#'
#' @return A `list` with settings for the pathfinder sampler.
#' @export
sampler_stan_pathfinder <- function(
    draws = 4000,
    num_paths = 4,
    single_path_draws = NULL,
    max_lbfgs_iters = NULL,
    num_elbo_draws = NULL,
    psis_resample = FALSE,
    seed = 0,
    refresh = 1000,
    show_messages = TRUE,
    show_exceptions = FALSE,
    ...
  ) {
  opts <- c(as.list(environment()), list(...))
  class(opts) <- "pathfinder"
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
#' @return A list with details of the stan model
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
  })

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
    function() {return(stanmodel_list[[i]])}
  })

  return(model_stan)
}

#' Compile EpiSewer models
#'
#' @description The stan models used by EpiSewer need to be compiled for your
#'   device. This is only necessary once, after installing or updating the
#'   package. This function compiles all models from the package.
#'
#' @param model A character vector with a specific model to compile. If NULL
#'   (default), all models are compiled.
#' @param force_recompile If TRUE, then models will be recompiled even if they
#'   have already been successfully compiled.
#' @param verbose If TRUE, warnings and detailed errors from the compilation are
#'   printed. This can help to diagnose compilation issues.
#'
#' @details If one or several models are not successfully compiled, please
#'   ensure that `cmdstan` is properly set up and try updating it to a newer
#'   version using [cmdstanr::install_cmdstan()]. If the problem persists,
#'   please run `sewer_compile(verbose = TRUE)` and post the output in
#'   a new issue on GitHub, along with your [cmdstanr::cmdstan_version()].
#'
#' @export
sewer_compile <- function(model = NULL, force_recompile = FALSE, verbose = FALSE) {
  if (!is.null(model)) {
    if (model %in% c("main", "EpiSewer_main")) {
      model <- "EpiSewer_main.stan"
    } else {
      cli::cli_abort(paste("The model", model, "is not available."))
    }
    models_to_compile <- model
  } else {
    models_to_compile <- c("EpiSewer_main.stan")
  }
  comp_success <- NULL

  for (i in 1:length(models_to_compile)) {
    model_name <- models_to_compile[i]
    cat("\r                                                      \r", sep = "");
    cat(sprintf("\r| Compiling model %d/%d ", i, length(models_to_compile)), sep = "")
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
      i, length(models_to_compile)), ifelse(success, "(success) ", "(failed) ")
      ), sep = "")
    Sys.sleep(0.5)
    comp_success <- c(comp_success, success)
  }
  cat("\r                                                        \r", sep = "");
  if(!all(comp_success)) {
    cli::cli_warn(
      paste(
        "The following models could not be compiled:",
        paste(models_to_compile[!comp_success], collapse = ", ")
      )
    )
  } else {

    cli::cli_alert("All models compiled successfully.")
  }
  return(invisible(NULL))
}

#' Computes a checksum of a stanmodel, i.e. of the code from its main
#' stan file and all included stan files
#' @keywords internal
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

#' Run stan's pathfinder algorithm to obtain improved inits.
#'
#' @description Helper function that runs the pathfinder algorithm with default
#'   inits and returns the result to be used as inits in mcmc sampling.
#'
#' @param stanmodel_instance A CmdStanModel object of the EpiSewer model
#'   constructed using `cmdstanr::cmdstan_model()`.
#' @param job An EpiSewer job object with data and default inits.
#' @param max_iters The maximum number of iterations for LBFGS during pathfinder
#'  initialization.
#'
#' @return Fitted model, a CmdStanPathfinder object. Can be used as arguments
#'   for the inits of the `CmdStanModel$sample()` method.
#'
#' @keywords internal
get_pathfinder_inits <- function(stanmodel_instance, job, max_iters = NULL) {
  arguments <- c(
    list(data = job$data),
    init = function() job$init,
    list(
      num_paths = 1, draws = 100, seed = 0, max_lbfgs_iters = max_iters,
      psis_resample = FALSE, sig_figs = 14,
      show_messages = FALSE, show_exceptions = FALSE
      )
  )
  return(fit_stan(
    stanmodel_instance, arguments, fit_method = "pathfinder", silent = TRUE
  ))
}

#' Fit model via stan
#'
#' @description Helper function to fit an EpiSewer model using one of stan's
#'   algorithms.
#'
#' @param stanmodel_instance A CmdStanModel object of the EpiSewer model
#'   constructed using `cmdstanr::cmdstan_model()`.
#' @param arguments Arguments to be passed to stan, including the data and the
#'   sampler options.
#' @param fit_method The algorithm used for sampling, currently supported are
#'   "mcmc" and "pathfinder".
#' @param silent Should the sampling be completely silent or should status
#'   updates and warnings be shown?
#'
#' @return Either a fitted stan model with draws, or a list with error messages
#'   and sampler outputs.
#'
#' @keywords internal
fit_stan <- function(stanmodel_instance, arguments, fit_method, silent = FALSE) {
  if (silent) {
    arguments$show_messages <- FALSE
    arguments$show_exceptions <- FALSE
    sink(tempfile(), type = "out")
    on.exit(sink())
  }

  return(tryCatch(
    {
      fit_res <- withWarnings(suppress_messages_warnings(
        {
          if (fit_method == "mcmc") {
            do.call(stanmodel_instance$sample, arguments)
          } else if (fit_method == "pathfinder") {
            do.call(stanmodel_instance$pathfinder, arguments)
          } else {
            cli::cli_abort("Unknown sampler type")
          }
        },
        c(
          "Registered S3 method overwritten by 'data.table'",
          "Cannot parse stat file, cannot read file: No such file or directory",
          "cannot open file '/proc/stat': No such file or directory",
          "Mean does not exist"
        )
      ))
      if (length(fit_res$warnings) == 0) {
        fit_res <- fit_res$value
        # ensure that data is read in
        fit_res$draws("R")
      } else {
        if (!silent) {
          cat("\n")
          cli::cli_warn(
            paste(
              "There was an error while fitting the model.",
              "Only the model input is returned."
            )
          )
        }
        fit_res <- list(
          errors = unlist(lapply(
            fit_res$warnings, function(x) stringr::str_remove(x$message, "\n")
          )),
          sampler_output = invisible(force(fit_res$value$output()))
        )
      }
      return(fit_res)
    },
    error = function(err) {
      if (!silent) {
        cat("\n")
        cli::cli_warn(c(
          paste(
            "There was an error while fitting the model.",
            "Only the model input is returned."
          ),
          err$message
        ))
      }
      return(list(errors = err, sampler_output = NULL))
    }
  ))
}

#' Convert draws from a previous stan model fit to inits for MCMC sampling
#'
#' @description Helper function that converts the draws obtained (e.g. from
#'   pathfinder) to a list of inits that can be used for MCMC sampling.
#'
#' @param fit A fitted model object (e.g. obtained from pathfinder).
#' @param model The model specification, used to identify scalar parameters.
#' @param num_chains The number of chains for which to generate inits.
#'
#' @return A list of inits for MCMC sampling, with one entry per chain.
#'
#' @keywords internal
draws_to_init <- function(fit, model, num_chains = 4) {
  rvars <- posterior::as_draws_rvars(fit$draws())
  # Get parameter names (exclude diagnostics ending in __)
  param_names <- names(rvars)[!grepl("__$", names(rvars))]
  scalar_params <- read_scalar_params(model)
  lapply(seq_len(num_chains), function(i) {
    draw_idx <- sample(posterior::ndraws(rvars), 1)
    lapply(setNames(param_names, param_names), function(p) {
      val <- posterior::draws_of(
        posterior::subset_draws(rvars[p], draw = draw_idx)[[1]]
      )
      # Remove only the first (draw) dimension
      dim_without_draw <- dim(val)[-1]
      val <- array(val, dim = dim_without_draw)
      if (length(val) == 1 && (p %in% scalar_params)) {
        val <- as.numeric(val)
      }
      return(val)
    })
  })
}

#' Fit the EpiSewer model
#'
#' @description Helper function that fits the EpiSewer model using the specified
#'   sampler and model instance, with the data and inits from the job object.
#'
#' @param job An EpiSewer job object containing the data, inits, and fitting
#'   options.
#' @param model The model specification, used to identify parameters.
#' @param model_instance A compiled model object.
#' @param run_silent Should the fitting be run silently? Default is FALSE.
#'
#' @return A fitted model object if successful, or a list with error messages
#'   and sampler output if there was an error during fitting.
#'
#' @keywords internal
fit_model <- function(job, model, model_instance, run_silent = FALSE) {
  # check if data, init and fit_opts are available in the job object
  if (!"data" %in% names(job)) {
    cli::cli_abort("The job object must contain a 'data' element.")
  }
  if (!"init" %in% names(job)) {
    cli::cli_abort("The job object must contain an 'init' element.")
  }
  if (!"fit_opts" %in% names(job)) {
    cli::cli_abort("The job object must contain a 'fit_opts' element.")
  }

  scalar_vars <- read_scalar_vars(model)
  for (var in names(job[["data"]])) {
    val <- job[["data"]][[var]]
    if (length(val) == 1 && is.null(dim(val)) && !(var %in% scalar_vars)) {
      job[["data"]][[var]] <- as.array(val)
    }
  }

  scalar_params <- read_scalar_params(model)
  for (param in names(job[["init"]])) {
    val <- job[["init"]][[param]]
    if (length(val) == 1 && is.null(dim(val)) && !(param %in% scalar_params)) {
      job[["init"]][[param]] <- as.array(val)
    }
  }

  arguments <- c(
    list(data = job$data),
    init = function() job$init,
    job$fit_opts$sampler
  )

  if (run_silent) {
    sink(tempfile(), type = "out")
    on.exit(sink())
  }

  # pathfinder initialization for mcmc
  use_mcmc <- class(job$fit_opts$sampler) == "mcmc"
  init_pathfinder <- job$fit_opts$sampler$init_pathfinder
  fitting_method <- class(job$fit_opts$sampler)
  if (use_mcmc && init_pathfinder) {
    tryCatch(
      {
        cat("Initializing chains via pathfinder...\n")
        pathfind_init <- get_pathfinder_inits(
          model_instance, job,
          max_iters = job$fit_opts$sampler$init_pathfinder_max_lbfgs_iters
        )
        stopifnot(!"errors" %in% names(pathfind_init))
        options(cmdstanr_warn_inits = FALSE)
        arguments$init <- draws_to_init(
          pathfind_init, model,
          num_chains = job$fit_opts$sampler$chains
          )
      },
      error = function(e) {
        cat(paste(
          "\nPathfinder initialization failed.",
          "\nFalling back to default initialization.\n\n"
        ))
      }
    )
  } else {
    options(cmdstanr_warn_inits = TRUE)
  }

  arguments[["init_pathfinder"]] <- NULL
  arguments[["init_pathfinder_max_lbfgs_iters"]] <- NULL

  fit_res <- fit_stan(model_instance, arguments, fitting_method)

  if (!"errors" %in% names(fit_res)) {
    if (job$results_opts$fitted) {
      try(fit_res$draws(), silent = TRUE)
      try(fit_res$init(), silent = TRUE)
      try(fit_res$profiles(), silent = TRUE)
    }
    if (class(job$fit_opts$sampler) == "mcmc") {
      try(suppressMessages(fit_res$diagnostic_summary()), silent = TRUE)
      try(fit_res$sampler_diagnostics(), silent = TRUE)
    }
    try(fit_res$time())
  }

  return(fit_res)
}

#' Identify scalar variables in a stan model and write them to a text file
#'
#' @param modelfolder The folder containing the stan model file.
#' @param modelname The name of the stan model file.
#'
#' @keywords internal
write_scalar_vars <- function(modelfolder = "inst/stan", modelname = "EpiSewer_main.stan") {
  mod <- cmdstanr::cmdstan_model(file.path(modelfolder, modelname))
  vars <- mod$variables()$data
  scalar_vars <- names(vars)[sapply(vars, function(v) v$dimensions == 0)]
  writeLines(scalar_vars, file.path(modelfolder, "scalar_data_vars.txt"))
}

#' Read scalar variables from a text file
#'
#' @param model The model specification, used to identify the correct text file.
#'
#' @return A character vector with the names of the scalar variables.
#'
#' @keywords internal
read_scalar_vars <- function(model) {
  if (model$package == "EpiSewer-docker") {
    filepath <- file.path(
      model$model_folder,
      "scalar_data_vars.txt"
    )
  } else {
    filepath <- system.file(
      model$model_folder,
      "scalar_data_vars.txt",
      package = model$package
    )
  }

  scalar_vars <- readLines(filepath)
  return(scalar_vars)
}

#' Identify scalar parameters in a stan model and write them to a text file
#'
#' @param modelfolder The folder containing the stan model file.
#' @param modelname The name of the stan model file.
#'
#' @keywords internal
write_scalar_params <- function(modelfolder = "inst/stan", modelname = "EpiSewer_main.stan") {
  mod <- cmdstanr::cmdstan_model(file.path(modelfolder, modelname))
  vars <- mod$variables()$parameters
  scalar_params <- names(vars)[sapply(vars, function(v) v$dimensions == 0)]
  writeLines(scalar_params, file.path(modelfolder, "scalar_param_vars.txt"))
}

#' Read scalar parameters from a text file
#'
#' @param model The model specification, used to identify the correct text file.
#'
#' @return A character vector with the names of the scalar parameters.
#'
#' @keywords internal
read_scalar_params <- function(model) {
  if (model$package == "EpiSewer-docker") {
    filepath <- file.path(
      model$model_folder,
      "scalar_param_vars.txt"
    )
  } else {
    filepath <- system.file(
      model$model_folder,
      "scalar_param_vars.txt",
      package = model$package
    )
  }

  scalar_params <- readLines(filepath)
  return(scalar_params)
}
