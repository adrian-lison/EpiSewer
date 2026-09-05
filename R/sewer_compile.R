# The modeldata compile pipeline.
#
# Model components (see model_component()) are inert specification objects.
# modeldata_compile() executes them against a single modeldata environment in
# a fixed, package-owned slot order, after publishing every component's
# arguments (and the user's central data/assumptions) to the spec tree
# md$.spec. Cross-cutting derived variables are resolved on demand through
# the derivation registry (sewer_derivations.R).

#' Create a fresh compile-time modeldata environment
#'
#' @description The environment has `emptyenv()` as parent so that any
#'   accidental use of `with(modeldata, ...)`-style scoping fails immediately
#'   instead of resolving symbols in the wrong place.
#' @keywords internal
md_new <- function() {
  md <- new.env(parent = emptyenv())
  md$.init <- list()
  md$.metainfo <- list()
  md$.checks <- list()
  md$.spec <- list()
  md$.staging <- list()
  class(md) <- c("modeldata_env", "modeldata")
  md
}

#' The package-owned execution order of component slots
#'
#' @description The only stated order in the system. Each entry is justified
#'   by the upstream state its component bodies read; derived variables are
#'   pulled on demand via [md_need()] and do not appear here.
#' @keywords internal
component_slot_order <- function() {
  c(
    "horizon",         # h, .metainfo$forecast_horizon; read by flows,
                       # sample_effects, R
    "concentrations",  # T, dates, n_measured/n_samples/n_averaged,
                       # measured concentrations; stashes the raw
                       # measurements and partition counts for downstream use
    "noise",           # reads observation_type (spec); decides emission of
                       # dPCR_total_partitions from the stash
    "LOD",             # LOD_estimate_dPCR reads cv_type/obs_dist from noise
    "flows",           # needs T_start/T_end/measured_dates (concentrations)
                       # + forecast_horizon (horizon)
    "residence_dist",  # independent (D, residence_dist)
    "shedding_dist",   # independent (S, shedding_dist, shedding_reference)
    "incubation_dist", # branches on shedding_reference => after shedding_dist
    "generation_dist", # independent (G); needed by seeding-length derivations
    "load_per_case",   # calibrate: needs stashed measurements/flows (cases
                       # branch) or the crude load curve (min_cases branch)
    "load_variation",  # estimate: needs crude infection/case curves
                       # => after load_per_case
    "outliers",        # needs T
    "sample_effects",  # needs T, T_start/T_end dates, forecast_horizon
    "seeding",         # needs initial_cases_crude, length_seeding, se
    "R",               # needs length_R_modeled, length_seeding,
                       # partial_window, partial_generation, forecast_horizon
    "infection_noise", # inits need infection_curve_crude, length_I
    "damping"          # trivial (forecast_damping)
  )
}

#' Map component slots to their modules (for the model structure summary)
#' @keywords internal
slot_module_map <- function() {
  c(
    concentrations = "measurements", noise = "measurements",
    LOD = "measurements",
    outliers = "sampling", sample_effects = "sampling",
    flows = "sewage", residence_dist = "sewage",
    incubation_dist = "shedding", shedding_dist = "shedding",
    load_per_case = "shedding", load_variation = "shedding",
    generation_dist = "infections", R = "infections",
    seeding = "infections", infection_noise = "infections",
    horizon = "forecast", damping = "forecast"
  )
}

#' Flatten modules/components into one slot-keyed component list
#' @keywords internal
flatten_components <- function(x) {
  components <- list()
  add_component <- function(comp) {
    if (!inherits(comp, "model_component")) {
      cli::cli_abort(
        "All model specifications must be model_component objects.",
        .internal = TRUE
      )
    }
    if (!is.null(components[[comp$slot]])) {
      cli::cli_abort(paste0(
        "The `", comp$slot, "` component was specified twice (via `",
        components[[comp$slot]]$name, "()` and `", comp$name, "()`)."
      ))
    }
    components[[comp$slot]] <<- comp
  }
  for (element in x) {
    if (inherits(element, "model_component")) {
      add_component(element)
    } else if (inherits(element, "EpiSewer_module") || is.list(element)) {
      for (comp in element) {
        add_component(comp)
      }
    } else {
      cli::cli_abort(
        "All model specifications must be model_component objects.",
        .internal = TRUE
      )
    }
  }
  components
}

#' Compile model components into a modeldata object
#'
#' @description Publishes all component arguments and central inputs to the
#'   spec tree, checks required user inputs (aggregated over all components),
#'   executes the component bodies in the package-owned slot order, resolves
#'   remaining derivations, inserts defaults, and runs registered checks.
#'
#' @param components A list of `model_component` objects and/or module
#'   bundles. Slots may be missing (partial compiles are supported, e.g. for
#'   tests); each slot may only be specified once.
#' @param data Central observation data, as provided by [sewer_data()].
#' @param assumptions Central assumptions, as provided by
#'   [sewer_assumptions()].
#'
#' @return A `modeldata` object (plain list) with stan data variables at the
#'   top level plus `.init`, `.metainfo`, `.checks`, `.str` and a trimmed
#'   `.spec`.
#' @keywords internal
modeldata_compile <- function(components, data = list(), assumptions = list()) {
  components <- flatten_components(components)
  md <- md_new()

  # seed the spec tree: central inputs + every component's captured arguments
  spec <- list(data = data, assumptions = assumptions)
  for (comp in components) {
    spec[[comp$slot]] <- c(list(.helper = comp$name), comp$args)
  }
  md$.spec <- spec

  preflight_inputs(components, md)

  for (slot in component_slot_order()) {
    if (!is.null(components[[slot]])) {
      apply_component(components[[slot]], md)
    }
  }

  settle_derivations(md)

  # defaults for optional stan variables (zero-length dimension hack)
  defaults <- modeldata_defaults()
  for (i in seq_along(defaults)) {
    path <- names(defaults)[[i]]
    if (!md_has_path(md, path)) {
      md_set_path(md, path, defaults[[i]])
    }
  }

  for (check in md$.checks) {
    check(md, data, assumptions)
  }

  md$.str <- spec_to_structure(md$.spec)
  md$.spec <- trim_spec(md$.spec)

  as_modeldata_list(md)
}

#' Execute one component body against the modeldata environment
#' @keywords internal
apply_component <- function(comp, md) {
  md$.current <- comp$name
  result <- tryCatch(
    comp$component(modeldata = md),
    error = function(e) {
      rlang::abort(
        paste0("Error in '", comp$name, "': ", conditionMessage(e)),
        call = NULL
      )
    }
  )
  if (!inherits(result, "modeldata")) {
    cli::cli_abort(
      paste0(
        "The component `", comp$name, "()` did not return the modeldata ",
        "object. Component expressions must end with (or explicitly ",
        "return) the modeldata."
      ),
      .internal = TRUE
    )
  }
  md$.current <- NULL
  invisible(md)
}

#' Check all components' required user inputs, aggregated
#'
#' @description For every component, its declared `requires_data` /
#'   `requires_assumptions` are checked against the component's own arguments
#'   and the central data/assumptions pools. All missing inputs are reported
#'   together in a single error.
#' @keywords internal
preflight_inputs <- function(components, md) {
  spec <- md$.spec
  provided <- function(slot, item) {
    options <- stringr::str_split(item, "\\|")[[1]]
    for (o in options) {
      if (!is.null(spec[[slot]][[o]]) ||
          !is.null(spec$data[[o]]) ||
          !is.null(spec$assumptions[[o]])) {
        return(TRUE)
      }
    }
    FALSE
  }
  resolve_reqs <- function(reqs) {
    if (is.function(reqs)) {
      reqs <- reqs(spec)
    }
    reqs
  }

  error_msg <- c()
  for (comp in components) {
    missing_data <- c()
    missing_assumptions <- c()
    for (item in resolve_reqs(comp$requires_data)) {
      if (!provided(comp$slot, item)) {
        missing_data <- c(missing_data, item)
      }
    }
    for (item in resolve_reqs(comp$requires_assumptions)) {
      if (!provided(comp$slot, item)) {
        missing_assumptions <- c(missing_assumptions, item)
      }
    }
    if (length(missing_data) + length(missing_assumptions) > 0) {
      message_data <- c()
      message_assumptions <- c()
      if (length(missing_data) > 0) {
        message_data <- paste(
          "Data (also via sewer_data()):",
          paste(
            "`",
            stringr::str_replace_all(missing_data, "\\|", " or "),
            "`",
            sep = "", collapse = ", "
          )
        )
      }
      if (length(missing_assumptions) > 0) {
        message_assumptions <- paste(
          "Assumptions (also via sewer_assumptions()):",
          paste(
            "`",
            stringr::str_replace_all(missing_assumptions, "\\|", " or "),
            "`",
            sep = "", collapse = ", "
          )
        )
      }
      error_msg <- c(
        error_msg,
        paste0(
          "Please provide the following information to ", comp$name, "():"
        ),
        message_data,
        message_assumptions
      )
    }
  }
  if (length(error_msg) > 0) {
    cli::cli_abort(error_msg, call = NULL)
  }
  invisible(md)
}

#' Resolve all remaining satisfiable derivations
#'
#' @description Runs after all component bodies: any registered derivation
#'   whose requirements are (recursively) satisfiable is resolved, so that
#'   derived metainfo is complete even when no component body pulled it.
#'   Unsatisfiable units are skipped silently (partial compiles).
#' @keywords internal
settle_derivations <- function(md) {
  units <- md_registry(md)
  for (id in names(units)) {
    unit <- units[[id]]
    all_provided <- all(vapply(
      unit$provides, function(p) md_has_path(md, p), logical(1)
    ))
    if (all_provided) {
      next
    }
    satisfiable <- all(vapply(
      unit$requires, function(r) md_can(md, r), logical(1)
    ))
    if (satisfiable) {
      md_need(md, unit$provides[[1]])
    }
  }
  invisible(md)
}

#' Generate the model structure summary from the spec tree
#'
#' @description Builds the `modelstructure` object (used by `job$model` and
#'   the print methods) from the helper names and non-default scalar
#'   arguments recorded in the spec tree.
#' @keywords internal
spec_to_structure <- function(spec) {
  modules <- slot_module_map()
  str <- list(
    measurements = list(), sampling = list(), sewage = list(),
    shedding = list(), infections = list(), forecast = list()
  )
  for (slot in names(modules)) {
    slot_spec <- spec[[slot]]
    if (is.null(slot_spec)) {
      next
    }
    helper <- slot_spec$.helper
    details <- component_str_details(helper, slot_spec)
    # setNames(list(...)) keeps the element even when details is NULL/empty
    str[[modules[[slot]]]][[slot]] <- setNames(list(details), helper)
  }
  class(str) <- "modelstructure"
  str
}

#' Determine display details for a component (non-default scalar arguments)
#' @keywords internal
component_str_details <- function(helper, slot_spec) {
  f <- get0(helper, envir = rlang::ns_env("EpiSewer"))
  if (is.null(f)) {
    return(c())
  }
  defaults <- formals(f)
  details <- c()
  for (nm in intersect(names(slot_spec), names(defaults))) {
    val <- slot_spec[[nm]]
    if (is.null(val) || !is.atomic(val) || length(val) != 1) {
      next
    }
    default_value <- tryCatch(
      eval(defaults[[nm]], envir = rlang::ns_env("EpiSewer")),
      error = function(e) NULL
    )
    same <- !is.null(default_value) && is.atomic(default_value) &&
      length(default_value) == 1 && isTRUE(all.equal(val, default_value))
    if (!same) {
      details[[nm]] <- val
    }
  }
  details
}

#' Trim large objects from the spec tree after compilation
#' @keywords internal
trim_spec <- function(spec) {
  trim_entry <- function(entry) {
    keep <- vapply(
      entry,
      function(a) {
        is.null(a) || (is.atomic(a) && length(a) <= 8) || is.character(a)
      },
      logical(1)
    )
    entry[keep]
  }
  for (nm in names(spec)) {
    spec[[nm]] <- trim_entry(spec[[nm]])
  }
  spec
}

#' Convert the compile-time environment to a plain modeldata list
#' @keywords internal
as_modeldata_list <- function(md) {
  internal <- c(".staging", ".current", ".deriving", ".derivation_index")
  nms <- setdiff(ls(md, all.names = TRUE), internal)
  out <- mget(nms, envir = md, ifnotfound = list(NULL))
  out <- out[!vapply(out, is.null, logical(1))]
  out <- out[sort(names(out), method = "radix")]
  class(out) <- "modeldata"
  out
}
