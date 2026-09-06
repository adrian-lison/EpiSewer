#' Create a model component object
#'
#' @description A `model_component` bundles the specification of one model
#'   component: the public helper name, the component slot it fills, the
#'   (normalized) arguments captured from the calling helper function, the
#'   user inputs it requires, and a function that applies the component to a
#'   modeldata object. Components are inert until executed by
#'   [modeldata_compile()], which runs them in a fixed slot order.
#'
#' @param name Name of the public helper function this component was created
#'   by (e.g. `"noise_estimate_dPCR"`). Used for error attribution and the
#'   model structure summary.
#' @param slot The component slot this component fills. Must be one of
#'   [all_components()].
#' @param .e An expression that defines the component through changes to the
#'   `modeldata` object. Must evaluate to (or explicitly return) the updated
#'   `modeldata`. The expression can use all captured arguments by name.
#' @param args Arguments captured for the component. By default, all variables
#'   in the calling function's frame (i.e. its arguments and any local
#'   variables computed before the `model_component()` call). These are
#'   published to the spec tree (`modeldata$.spec[[slot]]`) before any
#'   component body runs.
#' @param requires_data Character vector of data inputs the component needs
#'   (checked in the pre-flight pass). Alternatives can be expressed with
#'   `"a|b"`. Can also be a `function(spec)` returning such a vector, for
#'   requirements that depend on the model specification.
#' @param requires_assumptions Like `requires_data`, for assumptions.
#'
#' @return An object of class `model_component`.
#' @export
#' @keywords internal
model_component <- function(name, slot, .e,
                            args = as.list(parent.frame()),
                            requires_data = character(0),
                            requires_assumptions = character(0)) {
  if (!is.character(name) || length(name) != 1) {
    cli::cli_abort("`name` must be a single character string", .internal = TRUE)
  }
  if (!is.character(slot) || length(slot) != 1 ||
      !slot %in% all_components()) {
    cli::cli_abort(
      paste0("`slot` must be one of: ", paste(all_components(), collapse = ", ")),
      .internal = TRUE
    )
  }
  if (!is.list(args)) {
    cli::cli_abort("`args` must be a list", .internal = TRUE)
  }
  args[["modeldata"]] <- NULL # defensive: never capture a modeldata object

  func <- function() {}
  formals(func) <- c(args, list(modeldata = NULL))
  body(func) <- substitute(.e)

  component <- do.call(purrr::partial, c(list(func), args))

  structure(
    list(
      name = name,
      slot = slot,
      args = args,
      requires_data = requires_data,
      requires_assumptions = requires_assumptions,
      component = component
    ),
    class = "model_component"
  )
}

#' The model_component S4 class registration
#' @exportClass model_component
#' @keywords internal
setClass("model_component")

#' Print a model_component
#' @export
#' @keywords internal
print.model_component <- function(x, ...) {
  cat("<EpiSewer model component>", "\n")
  cat("Helper:", x$name, "\n")
  cat("Slot:", x$slot, "\n")
  if (length(x$args) > 0) {
    show_args <- sapply(x$args, format_arg_value)
    cat(
      "Arguments:\n",
      paste(" ", names(x$args), "=", show_args, collapse = "\n"), "\n"
    )
  }
  invisible(x)
}

#' Compact display of an argument value
#' @keywords internal
format_arg_value <- function(a) {
  if (is.data.frame(a)) {
    paste0("<data.frame [", nrow(a), " x ", ncol(a), "]>")
  } else if (is.matrix(a)) {
    paste0("<matrix [", nrow(a), " x ", ncol(a), "]>")
  } else if (is.function(a)) {
    "<function>"
  } else if (is.null(a)) {
    "NULL"
  } else if (is.atomic(a) && length(a) <= 8) {
    paste(format(a, digits = 4), collapse = ", ")
  } else {
    paste0("<", class(a)[1], " [", length(a), "]>")
  }
}

#' Verify that an argument is a model_component for the right slot
#'
#' @description Used by the module functions to give a helpful error when
#'   something else than a suitable component helper result is supplied.
#'
#' @param component The object to check.
#' @param arg_name Name of the module function argument (equal to the slot).
#' @keywords internal
verify_is_component <- function(component, arg_name) {
  if (!inherits(component, "model_component")) {
    error_msg <- paste0(
      "The argument `", arg_name, "` expects an EpiSewer model component as ",
      "input. Please use a function of the form `", arg_name, "_`", " to ",
      "specify this argument."
    )
    all_fs <- names(rlang::ns_env("EpiSewer"))
    arg_fs <- all_fs[stringr::str_detect(
      # this excludes functions ending with _
      all_fs, paste0("^", arg_name, "_", ".*(?<!_)$")
    )]
    if (length(arg_fs) > 0) {
      functions <- paste0(arg_fs, "()")
      functions_cli <- paste0(
        "{.help [", functions, "](EpiSewer::", functions, "}"
      )
      error_msg <- paste(error_msg, "Available functions:")
      error_msg <- c(error_msg, functions_cli)
      names(error_msg) <- c("!", rep("*", length(error_msg) - 1))
    }
    cli::cli_abort(error_msg, call = rlang::caller_env())
  }
  if (component$slot != arg_name) {
    cli::cli_abort(
      paste0(
        "The argument `", arg_name, "` received the component `",
        component$name, "()`, which specifies the `", component$slot,
        "` component."
      ),
      call = rlang::caller_env()
    )
  }
  invisible(component)
}

#' Create a module bundle of model components
#'
#' @param module Name of the module (e.g. "measurements").
#' @param components Named list of `model_component` objects, keyed by slot.
#' @return The component list, classed as `EpiSewer_module`, with the module
#'   name attached as an attribute.
#' @keywords internal
new_module <- function(module, components) {
  structure(
    components,
    module = module,
    class = "EpiSewer_module"
  )
}

#' Print an EpiSewer module
#' @export
#' @keywords internal
print.EpiSewer_module <- function(x, ...) {
  cat(attr(x, "module"), "\n")
  for (slot in names(x)) {
    cat(" |-", slot, "=", paste0(x[[slot]]$name, "()"), "\n")
  }
  invisible(x)
}
