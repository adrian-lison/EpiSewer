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

**pkgdown site (README and vignettes in the browser)**
- To build the full site and serve it, run `episewer-site`. It takes about 40 s. Output goes to `docs/`, which is gitignored, and it is served on sandbox port 8000.
- If only the README or home page changed, run `Rscript -e 'pkgdown::build_home()'` and then `episewer-site --no-build`. For vignette changes, rerun `episewer-site` or `episewer-site --article <name>`. Plain `pkgdown::build_article()` fails, because the vignettes need the package to be installed.
- To show the user only one vignette, run `episewer-site --article <name>`, where `<name>` is the `.Rmd` file name without its extension, e.g. `detailed-example`.
  - It takes about 10–15 s.
  - It builds just that vignette into a separate folder and serves only it, and `/` redirects to it. Other pages are not exposed.
  - Give the user the URL it prints.
- Stop the server with `episewer-site --stop`.

**Temporary ("local") vignettes vs. official vignettes**
- If the user asks for a *temporary* or *local* vignette, they want it only to inspect your work visually.
  - Create it as `vignettes/local/<name>.Rmd`, following the YAML header and setup chunk of the existing vignettes.
  - Build and expose it on its own with `episewer-site --article <name>`, and give the user the printed URL.
- Temporary vignettes are **never committed by default.**
  - `vignettes/local/` is gitignored and `.Rbuildignore`d. Never `git add -f` it.
  - Do not link to it from official docs.
  - Do not list it in `_pkgdown.yml`. The `internal` section with `starts_with("local/")` already covers it and must stay.
  - When the user asks you to commit, commit everything else as usual and mention that the temporary vignette was left uncommitted.
- Only when the user explicitly says to make a temporary vignette *official* or *public*:
  1. Move it with `mv vignettes/local/<name>.Rmd vignettes/<name>.Rmd`, along with any figures or files it uses. Fix relative paths if needed.
  2. Add `<name>` to a suitable section under `articles:` in `_pkgdown.yml`. If no section obviously fits, ask which one.
  3. Rebuild the full site with `episewer-site` to check that it renders and appears in the navbar.
  4. Commit the vignette and the `_pkgdown.yml` change.
- Site builds write plot files into `vignettes/figures/`, because the vignettes set `fig.path = "figures/..."`. Do not commit new or changed files there unless the user asked to re-render the vignettes; restore or delete them before committing.
- If `README.Rmd` changed, render `README.md` first: `Rscript -e 'rmarkdown::render("README.Rmd", output_format = "github_document"); unlink("README.html")'`.
- You cannot publish ports yourself. If the user cannot reach the site, tell them to run `sbx port publish <sandbox> 8000:8000` on the host and open http://localhost:8000.

**Dependencies**
- `episewer-check-deps` checks that every DESCRIPTION dependency is installed. Add `--install` to fetch missing ones from CRAN.
- If you add a dependency, also add it to `.devenv/env.yml` as `r-<name>`, so the dev environments stay in sync.
