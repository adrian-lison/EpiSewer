# post-deploy bootstrap for dev environment
options(
  repos = c(CRAN = "https://cloud.r-project.org")
)

cran_pkgs <- c("extraDistr", "draw")
to_install <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]

if (length(to_install) > 0) {
  message(paste("Installing missing CRAN packages:", paste(to_install, collapse = ", ")))
  install.packages(to_install)
}

if (!requireNamespace("emo", quietly = TRUE)) {
  message("Installing missing package: hadley/emo from github")
  remotes::install_github("hadley/emo")
}

if (!requireNamespace("vscDebugger", quietly = TRUE)) {
  message("Installing missing package: vscDebugger from R-universe")
  install.packages("vscDebugger", repos = "https://manuelhentschel.r-universe.dev")
}

message("Post-installation completed.")