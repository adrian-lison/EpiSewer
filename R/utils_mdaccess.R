# Accessors for the compile-time modeldata environment.
#
# During compilation (see modeldata_compile()), modeldata is an environment,
# which gives md_need() memoization by reference: a derived value is computed
# at most once and lands in the modeldata itself. After compilation the
# modeldata is converted to a plain list.

#' Split a `$`-path into its levels
#' @keywords internal
md_path_levels <- function(path) {
  stringr::str_split(path, "\\$")[[1]]
}

#' Get a (possibly nested) modeldata variable by `$`-path
#'
#' @param md Modeldata (environment during compile, or plain list).
#' @param path A `$`-separated path, e.g. `"T"` or `".metainfo$length_I"`.
#' @return The value, or `NULL` if any level is absent.
#' @keywords internal
md_get_path <- function(md, path) {
  levels <- md_path_levels(path)
  value <- if (is.environment(md)) {
    get0(levels[1], envir = md, inherits = FALSE)
  } else {
    md[[levels[1]]]
  }
  for (l in levels[-1]) {
    if (is.null(value)) {
      return(NULL)
    }
    value <- value[[l]]
  }
  value
}

#' Set a (possibly nested) modeldata variable by `$`-path
#' @keywords internal
md_set_path <- function(md, path, value) {
  levels <- md_path_levels(path)
  if (length(levels) == 1) {
    if (is.environment(md)) {
      assign(levels[1], value, envir = md)
    } else {
      md[[levels[1]]] <- value
    }
  } else {
    top <- if (is.environment(md)) {
      get0(levels[1], envir = md, inherits = FALSE)
    } else {
      md[[levels[1]]]
    }
    if (is.null(top)) {
      top <- list()
    }
    # only two levels occur in practice (.metainfo$x, .init$x)
    if (length(levels) > 2) {
      cli::cli_abort("md_set_path supports at most 2 levels", .internal = TRUE)
    }
    top[[levels[2]]] <- value
    if (is.environment(md)) {
      assign(levels[1], top, envir = md)
    } else {
      md[[levels[1]]] <- top
    }
  }
  invisible(md)
}

#' Check presence of a modeldata variable by `$`-path
#' @keywords internal
md_has_path <- function(md, path) {
  !is.null(md_get_path(md, path))
}

#' Read a component argument from the spec tree
#'
#' @description The spec tree (`md$.spec`) holds the captured arguments of
#'   every model component, published before any component body runs. This is
#'   the supported way for one component (or derivation) to read facts rooted
#'   in another component's arguments.
#'
#' @param md Modeldata environment.
#' @param slot The component slot to read from (e.g. `"noise"`).
#' @param name The argument name. If missing, the whole slot spec (including
#'   `$.helper`, the helper function name) is returned.
#' @param default Value returned when the slot or argument is absent (also
#'   when the chosen helper for that slot does not have the argument). If
#'   missing, absence is an error.
#' @keywords internal
md_spec <- function(md, slot, name, default) {
  spec <- md_get_path(md, ".spec")[[slot]]
  if (missing(name)) {
    return(spec)
  }
  value <- spec[[name]]
  if (is.null(value)) {
    if (missing(default)) {
      helper <- spec[[".helper"]]
      cli::cli_abort(
        paste0(
          "The model specification does not provide `", name, "` for the `",
          slot, "` component",
          if (!is.null(helper)) paste0(" (specified via ", helper, "())"),
          "."
        ),
        .internal = TRUE
      )
    }
    return(default)
  }
  value
}

#' Retrieve a user input for a component
#'
#' @description Implements the input precedence rule: an input passed directly
#'   to the component helper beats the same input supplied centrally via
#'   [sewer_data()] / [sewer_assumptions()]. If both are supplied with
#'   different values, an error is thrown.
#'
#' @param md Modeldata environment.
#' @param slot The component slot requesting the input.
#' @param arg Name of the helper argument holding the direct input.
#' @param central Name of the input in the central data/assumptions list
#'   (defaults to `arg`).
#' @param type One of `"data"` or `"assumptions"`.
#' @return The input value (may be `NULL` if not supplied anywhere).
#' @keywords internal
take_input <- function(md, slot, arg, central = arg,
                       type = c("data", "assumptions")) {
  type <- rlang::arg_match(type)
  direct <- md_spec(md, slot, arg, default = NULL)
  pool <- md_get_path(md, ".spec")[[type]]
  central_value <- pool[[central]]

  if (!is.null(direct) && !is.null(central_value)) {
    # Assumptions with defaults in sewer_assumptions() are only conflict-
    # checked when they were explicitly supplied by the user.
    check <- TRUE
    if (type == "assumptions") {
      explicit <- attr(pool, "explicit")
      if (is.null(explicit)) {
        # plain list without explicitness info: fall back to checking all
        # names except those with non-NULL defaults in sewer_assumptions()
        defaults <- sewer_assumptions()
        explicit <- setdiff(
          names(pool), names(defaults)[!sapply(defaults, is.null)]
        )
      }
      check <- central %in% explicit
    }
    if (check && !isTRUE(all.equal(direct, central_value))) {
      helper <- md_spec(md, slot, ".helper", default = slot)
      if (type == "data") {
        cli::cli_abort(paste0(
          "You provided different data for `", central,
          "` in sewer_data() and ", helper, "()."
        ))
      } else {
        cli::cli_abort(paste0(
          "You provided different assumptions for `", central,
          "` in sewer_assumptions() and ", helper, "()."
        ))
      }
    }
  }

  if (!is.null(direct)) direct else central_value
}

#' Get a modeldata variable, deriving it on demand if necessary
#'
#' @description `md_need()` is the single resolution mechanism for
#'   cross-cutting derived variables: if the variable is already present in
#'   the modeldata it is returned; otherwise its registered derivation (see
#'   [modeldata_derivations()]) is resolved recursively, the results are
#'   memoized into the modeldata, and the value is returned. Cycles are
#'   detected and reported with the full derivation chain.
#'
#' @param md Modeldata environment.
#' @param path A `$`-separated variable path, e.g. `"se"` or
#'   `".metainfo$length_I"`.
#' @keywords internal
md_need <- function(md, path) {
  if (md_has_path(md, path)) {
    return(md_get_path(md, path))
  }
  if (!is.environment(md)) {
    cli::cli_abort(
      paste0("`", path, "` is not available in the compiled modeldata."),
      .internal = TRUE
    )
  }

  index <- md_derivation_index(md)
  unit_id <- index[[path]]
  if (is.null(unit_id)) {
    descriptions <- modeldata_descriptions()
    description <- descriptions[[path]]
    current <- md_get_path(md, ".current")
    cli::cli_abort(paste0(
      "`", path, "`",
      if (!is.null(description)) paste0(" (", description, ")"),
      " is required",
      if (!is.null(current)) paste0(" by ", current, "()"),
      ", but is not present in the model and has no registered derivation. ",
      "A model component that provides it may be missing."
    ))
  }

  deriving <- md_get_path(md, ".deriving")
  if (path %in% deriving) {
    chain <- paste(c(deriving, path), collapse = " <- ")
    cli::cli_abort(
      paste0("Circular dependency between derived variables: ", chain),
      .internal = TRUE
    )
  }
  md_set_path(md, ".deriving", c(deriving, path))

  unit <- md_registry(md)[[unit_id]]
  for (req in unit$requires) {
    md_need(md, req)
  }
  result <- unit$fn(md)
  for (p in unit$provides) {
    md_set_path(md, p, result[[p]])
  }

  md_set_path(md, ".deriving", deriving)
  md_get_path(md, path)
}

#' Check whether a modeldata variable is present or derivable
#' @keywords internal
md_can <- function(md, path, visiting = character(0)) {
  if (md_has_path(md, path)) {
    return(TRUE)
  }
  if (!is.environment(md)) {
    return(FALSE)
  }
  if (path %in% visiting) {
    return(FALSE) # cycle: not derivable this way
  }
  index <- md_derivation_index(md)
  unit_id <- index[[path]]
  if (is.null(unit_id)) {
    return(FALSE)
  }
  unit <- md_registry(md)[[unit_id]]
  all(vapply(
    unit$requires,
    function(req) md_can(md, req, visiting = c(visiting, path)),
    logical(1)
  ))
}

#' Get the derivation registry for a modeldata environment
#'
#' @description Returns the package registry by default; a modeldata
#'   environment can carry an override in `.derivation_registry` (used by
#'   unit tests).
#' @keywords internal
md_registry <- function(md) {
  registry <- md_get_path(md, ".derivation_registry")
  if (is.null(registry)) {
    registry <- modeldata_derivations()
  }
  registry
}

#' Get (and lazily build) the provides-path -> derivation-unit index
#' @keywords internal
md_derivation_index <- function(md) {
  index <- md_get_path(md, ".derivation_index")
  if (is.null(index)) {
    units <- md_registry(md)
    index <- list()
    for (i in seq_along(units)) {
      for (p in units[[i]]$provides) {
        if (!is.null(index[[p]])) {
          cli::cli_abort(
            paste0("Derivation for `", p, "` registered twice."),
            .internal = TRUE
          )
        }
        index[[p]] <- i
      }
    }
    md_set_path(md, ".derivation_index", index)
  }
  index
}
