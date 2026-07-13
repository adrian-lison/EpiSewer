#!/usr/bin/env bash
set -eo pipefail

ENV_PATH=".conda/env"

# Create or update the conda environment
mamba env create --prefix "$ENV_PATH" -f .devenv/env.yml || \
mamba env update --prefix "$ENV_PATH" -f .devenv/env.yml --prune

# Run the post-installation R script to install further R packages
echo "Running post-installation R script..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_PATH"
Rscript .devenv/env-post.R