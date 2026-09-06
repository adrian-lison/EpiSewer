# Unit tests for the model_component / modeldata_compile core machinery.

make_test_component <- function(name = "widget_observe",
                                slot = "concentrations",
                                value = 1,
                                requires_data = character(0),
                                requires_assumptions = character(0)) {
  # emulates a component helper: local args are captured automatically
  helper <- function(value_arg = value, extra = "abc") {
    model_component(name, slot,
      {
        modeldata$test_value <- value_arg
        modeldata$.metainfo$test_extra <- extra
        return(modeldata)
      },
      requires_data = requires_data,
      requires_assumptions = requires_assumptions
    )
  }
  helper()
}

test_that("model_component captures arguments and executes against modeldata", {
  comp <- make_test_component(value = 42)
  expect_s3_class(comp, "model_component")
  expect_equal(comp$args$value_arg, 42)
  expect_equal(comp$args$extra, "abc")

  md <- md_new()
  result <- comp$component(modeldata = md)
  expect_true(inherits(result, "modeldata"))
  expect_equal(md$test_value, 42)
  expect_equal(md$.metainfo$test_extra, "abc")
})

test_that("model_component validates name and slot", {
  expect_error(
    model_component(42, "noise", {
      modeldata
    }),
    "single character"
  )
  expect_error(
    model_component("x", "not_a_slot", {
      modeldata
    }),
    "must be one of"
  )
})

test_that("modeldata_compile publishes args to the spec tree before bodies run", {
  # the concentrations body reads an argument of the noise component
  comp_a <- local({
    x <- 5
    model_component("comp_a", "concentrations", {
      modeldata$seen_noise_arg <- md_spec(modeldata, "noise", "my_arg")
      return(modeldata)
    })
  })
  comp_b <- local({
    my_arg <- "from_noise"
    model_component("comp_b", "noise", {
      return(modeldata)
    })
  })
  md <- modeldata_compile(list(comp_a, comp_b))
  expect_equal(md$seen_noise_arg, "from_noise")
})

test_that("md_spec errors without default and falls back with default", {
  md <- md_new()
  md$.spec <- list(noise = list(.helper = "noise_estimate", replicates = TRUE))
  expect_true(md_spec(md, "noise", "replicates"))
  expect_equal(md_spec(md, "noise", "absent", default = 99), 99)
  expect_equal(md_spec(md, "unspecified_slot", "x", default = 7), 7)
  expect_error(md_spec(md, "noise", "absent"), "does not provide")
})

test_that("take_input prefers direct args and detects conflicts", {
  md <- md_new()
  md$.spec <- list(
    data = list(measurements = data.frame(a = 1)),
    assumptions = list(),
    concentrations = list(
      .helper = "concentrations_observe",
      measurements = data.frame(a = 2)
    )
  )
  expect_error(
    take_input(md, "concentrations", "measurements", type = "data"),
    "different data for `measurements` in sewer_data\\(\\) and concentrations_observe\\(\\)"
  )
  # equal values: no conflict, direct wins
  md$.spec$data$measurements <- data.frame(a = 2)
  expect_equal(
    take_input(md, "concentrations", "measurements", type = "data"),
    data.frame(a = 2)
  )
  # only central
  md$.spec$concentrations$measurements <- NULL
  expect_equal(
    take_input(md, "concentrations", "measurements", type = "data"),
    data.frame(a = 2)
  )
})

test_that("take_input exempts defaulted assumptions from conflict checks", {
  md <- md_new()
  # plain list without explicitness info (like the packaged assumption sets):
  # min_cases has a non-NULL default in sewer_assumptions() -> exempt
  md$.spec <- list(
    data = list(),
    assumptions = list(min_cases = 10, generation_dist = c(0.5, 0.5)),
    load_per_case = list(.helper = "load_per_case_calibrate", min_cases = 20),
    generation_dist = list(
      .helper = "generation_dist_assume", generation_dist = c(0.4, 0.6)
    )
  )
  expect_equal(
    take_input(md, "load_per_case", "min_cases", type = "assumptions"),
    20
  )
  expect_error(
    take_input(
      md, "generation_dist", "generation_dist", type = "assumptions"
    ),
    "different assumptions for `generation_dist`"
  )
  # with explicitness info, an explicitly supplied min_cases is checked
  attr(md$.spec$assumptions, "explicit") <- c("min_cases", "generation_dist")
  expect_error(
    take_input(md, "load_per_case", "min_cases", type = "assumptions"),
    "different assumptions for `min_cases`"
  )
})

test_that("preflight aggregates missing inputs over all components", {
  comp_a <- make_test_component(
    name = "widget_observe", slot = "concentrations",
    requires_data = "measurements"
  )
  comp_b <- local({
    limit_of_detection <- NULL
    model_component("lod_thing", "LOD",
      {
        modeldata
      },
      requires_assumptions = "limit_of_detection"
    )
  })
  err <- expect_error(
    modeldata_compile(list(comp_a, comp_b)),
    "Please provide the following information to widget_observe"
  )
  expect_match(conditionMessage(err), "lod_thing")
  expect_match(conditionMessage(err), "limit_of_detection")
  expect_match(conditionMessage(err), "measurements")

  # provided centrally: no error
  md <- modeldata_compile(
    list(comp_a, comp_b),
    data = list(measurements = data.frame(a = 1)),
    assumptions = list(limit_of_detection = 2)
  )
  expect_true(inherits(md, "modeldata"))
})

test_that("preflight supports alternatives and function-form requirements", {
  comp <- make_test_component(
    name = "calib", slot = "load_per_case",
    requires_data = "min_cases|cases"
  )
  expect_error(
    modeldata_compile(list(comp)),
    "min_cases or cases"
  )
  expect_no_error(
    modeldata_compile(list(comp), data = list(cases = data.frame(a = 1)))
  )

  comp_fn <- local({
    model_component("conditional", "incubation_dist",
      {
        modeldata
      },
      requires_assumptions = function(spec) {
        if (identical(spec$assumptions$mode, "on")) "special_input"
        else character(0)
      }
    )
  })
  expect_no_error(modeldata_compile(list(comp_fn)))
  expect_error(
    modeldata_compile(
      list(comp_fn), assumptions = list(mode = "on")
    ),
    "special_input"
  )
})

test_that("md_need resolves chained derivations with memoization", {
  registry <- list(
    a = md_derive("a", requires = "base", fn = function(md) {
      counter <- md_get_path(md, "counter_a")
      md$counter_a <- (if (is.null(counter)) 0 else counter) + 1
      list(a = md$base + 1)
    }),
    b = md_derive("b", requires = "a", fn = function(md) {
      list(b = md$a * 10)
    }),
    multi = md_derive(
      c("c1", ".metainfo$c2"),
      requires = "b",
      fn = function(md) {
        list(c1 = md$b + 1, ".metainfo$c2" = md$b + 2)
      }
    )
  )
  md <- md_new()
  md$.derivation_registry <- registry
  md$base <- 1

  expect_equal(md_need(md, "c1"), 21)
  expect_equal(md$a, 2)
  expect_equal(md$b, 20)
  expect_equal(md$.metainfo$c2, 22)
  # memoized: repeated pulls do not recompute
  expect_equal(md_need(md, "a"), 2)
  expect_equal(md$counter_a, 1)
})

test_that("md_need reports cycles with the derivation chain", {
  registry <- list(
    x = md_derive("x", requires = "y", fn = function(md) list(x = md$y)),
    y = md_derive("y", requires = "x", fn = function(md) list(y = md$x))
  )
  md <- md_new()
  md$.derivation_registry <- registry
  err <- expect_error(md_need(md, "x"), "Circular dependency")
  expect_match(conditionMessage(err), "x <- y <- x")
})

test_that("md_need errors helpfully for underivable variables", {
  md <- md_new()
  md$.derivation_registry <- list()
  md$.current <- "some_component"
  err <- expect_error(md_need(md, "nonexistent"), "no registered derivation")
  expect_match(conditionMessage(err), "some_component")
})

test_that("md_can checks recursive derivability", {
  registry <- list(
    a = md_derive("a", requires = "base", fn = function(md) list(a = 1)),
    b = md_derive("b", requires = "a", fn = function(md) list(b = 2))
  )
  md <- md_new()
  md$.derivation_registry <- registry
  expect_false(md_can(md, "b"))
  md$base <- 0
  expect_true(md_can(md, "b"))
  expect_false(md_can(md, "unknown"))
})

test_that("settle pass resolves satisfiable derivations", {
  registry <- list(
    a = md_derive("a", requires = "base", fn = function(md) list(a = md$base + 1)),
    unsat = md_derive("u", requires = "missing_thing", fn = function(md) list(u = 1))
  )
  comp <- local({
    model_component("basecomp", "concentrations", {
      modeldata$base <- 10
      return(modeldata)
    })
  })
  # inject registry through a wrapper component executed first
  md <- md_new()
  md$.derivation_registry <- registry
  md$.spec <- list(
    data = list(), assumptions = list(),
    concentrations = c(list(.helper = "basecomp"), comp$args)
  )
  apply_component(comp, md)
  settle_derivations(md)
  expect_equal(md$a, 11)
  expect_false(md_has_path(md, "u"))
})

test_that("apply_component attributes errors and enforces the return contract", {
  bad <- local({
    model_component("boomer", "noise", {
      cli::cli_abort("something exploded")
    })
  })
  md <- md_new()
  err <- expect_error(apply_component(bad, md), "Error in 'boomer'")
  expect_match(conditionMessage(err), "something exploded")

  noreturn <- local({
    model_component("forgetful", "noise", {
      modeldata$x <- 1
      "not the modeldata"
    })
  })
  md <- md_new()
  expect_error(apply_component(noreturn, md), "did not return the modeldata")
})

test_that("as_modeldata_list sorts, drops NULLs and internal slots", {
  md <- md_new()
  md$zeta <- 1
  md$alpha <- 2
  md$dropme <- NULL
  md$.staging <- list(big = 1:100)
  md$.deriving <- character(0)
  out <- as_modeldata_list(md)
  expect_s3_class(out, "modeldata")
  expect_false(any(c(".staging", ".deriving", "dropme") %in% names(out)))
  expect_true(which(names(out) == "alpha") < which(names(out) == "zeta"))
})

test_that("spec_to_structure builds a modelstructure from the spec tree", {
  spec <- list(
    data = list(), assumptions = list(),
    noise = list(.helper = "noise_estimate", replicates = TRUE),
    horizon = list(.helper = "horizon_assume", horizon = 7)
  )
  str <- spec_to_structure(spec)
  expect_s3_class(str, "modelstructure")
  expect_equal(names(str$measurements$noise), "noise_estimate")
  expect_equal(str$measurements$noise$noise_estimate[["replicates"]], TRUE)
  expect_equal(names(str$forecast$horizon), "horizon_assume")
  expect_output(print(str), "noise_estimate")
})

test_that("duplicate slots are rejected and modules are flattened", {
  comp1 <- make_test_component(name = "c1", slot = "noise")
  comp2 <- make_test_component(name = "c2", slot = "noise")
  expect_error(
    modeldata_compile(list(comp1, comp2)),
    "specified twice"
  )
  module <- new_module("measurements", list(noise = comp1))
  md <- modeldata_compile(list(module))
  expect_equal(md$test_value, 1)
})

test_that("verify_is_component gives helpful errors", {
  expect_error(
    verify_is_component(list(), "noise"),
    "expects an EpiSewer model component"
  )
  comp <- make_test_component(slot = "LOD", name = "LOD_none")
  expect_error(
    verify_is_component(comp, "noise"),
    "specifies the `LOD` component"
  )
  expect_no_error(verify_is_component(comp, "LOD"))
})
