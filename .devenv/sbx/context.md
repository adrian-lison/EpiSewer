## EpiSewer development environment (sbx)

R, CmdStan and all package dependencies live in the conda env `/opt/episewer-env`, which is already on `PATH`. `CMDSTAN` is set. Do not use `conda activate`, and do not create a new environment.

**Running code and tests**
- Load the package: `Rscript -e 'devtools::load_all()'`.
- Run all tests: `Rscript -e 'devtools::test()'`. For one file, use `devtools::test(filter = "fit_default")`.
- `tests/testthat/setup_fit.R` runs `sewer_compile(model = "main")`. The first run compiles the Stan model, which takes about 1–2 minutes, and caches the executable at `inst/stan/EpiSewer_main`.
- `test-docker-backend.R` skips itself, because there is no `episewer` docker image here.
- Package check: `Rscript -e 'devtools::check(document = FALSE, vignettes = FALSE)'`.
- Documentation: `Rscript -e 'devtools::document()'`, after editing roxygen comments.

**Stan models (`inst/stan/`)**
- Fast syntax and type check, with no C++ compile: `$CMDSTAN/bin/stanc --include-paths=inst/stan inst/stan/EpiSewer_main.stan --o=/tmp/check.hpp`.
- After changing `.stan` files:
  - Recompile with `Rscript -e 'devtools::load_all(); sewer_compile(force_recompile = TRUE)'`.
  - Refresh the metadata with `Rscript docker/prepare_docker.R`, which updates `stan_digest.txt` and the `scalar_*_vars.txt` files.
- Never commit compiled executables, such as `inst/stan/EpiSewer_main`.

**Linting and style**
- After editing an R file, run `Rscript -e 'lintr::lint("R/<file>.R")'` and `Rscript -e 'styler::style_file("R/<file>.R")'`.
- Fix the lints you introduced. Do not restyle files you did not touch.
- Repo-wide: `lintr::lint_package()`. It uses the repo `.lintr` config.

**Dependencies**
- `episewer-check-deps` checks that every DESCRIPTION dependency is installed. Add `--install` to fetch missing ones from CRAN.
- If you add a dependency, also add it to `.devenv/env.yml` as `r-<name>`, so the dev environments stay in sync.
