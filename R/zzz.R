.onAttach <- function(...) {
  inform_update_available()
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
}

is_loading_for_tests <- function() {
  !interactive() && identical(Sys.getenv("DEVTOOLS_LOAD"), "EpiSewer")
}

inform_update_available <- function() {
  pkg_name <- "EpiSewer"
  installed_version <- utils::packageVersion("EpiSewer")

  if (curl::has_internet()) {
    try({
      # Get DESCRIPTION file on main from Github
      github_url <- "https://github.com/adrian-lison/EpiSewer/raw/main/DESCRIPTION"
      description_file <- readr::read_file(github_url)

      # Extract latest version from DESCRIPTION
      latest_version <- desc::desc(text = description_file)$get_version()

      if (installed_version < latest_version) {
        packageStartupMessage(paste0(
          'A newer version of EpiSewer is available (v',
          latest_version,
          ').\nPlease update by running ',
          '\'remotes::install_github("adrian-lison/EpiSewer", dependencies = TRUE)\'',
          ' and check for potential breaking changes.'
          )
        )
      }
    }, silent = TRUE)
  }
  return()
}
