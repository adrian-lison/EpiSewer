# Post-install bootstrap for the sbx sandbox env (counterpart of ../env-post.R).
# Runs at image build time. Use it for packages that are excluded from the
# conda solve (see env-sbx-exclude.txt) but still needed, and for non-conda deps.
options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing <- function(pkgs, ...) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    message("Installing: ", paste(missing, collapse = ", "))
    install.packages(missing, ...)
  }
}

# CRAN packages not taken from conda-forge
install_if_missing(c("extraDistr"))

# GitHub packages
if (!requireNamespace("emo", quietly = TRUE)) {
  remotes::install_github("hadley/emo", upgrade = "never")
}

# vscDebugger is deliberately omitted (no VS Code inside the sandbox).

message("sbx post-installation completed.")
