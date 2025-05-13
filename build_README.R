allcontributors::add_contributors(format = "text")

rmarkdown::render("README.Rmd", output_format = "github_document")
unlink("README.html")
