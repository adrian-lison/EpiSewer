allcontributors::add_contributors(format = "text")

# Build markdown file
rmarkdown::render("README.Rmd", output_format = "github_document")
unlink("README.html")

# Build jupyter notebook (requires jupytext)
# You can install it via conda using
# conda create -n jupytext-env -c conda-forge jupytext
system2(
  "/opt/anaconda3/bin/conda",
  args = c(
    "run", "-n", "jupytext-env",
    "jupytext", "--to", "notebook", "README.Rmd"
    )
)
