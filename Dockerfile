# Docker image for fitting EpiSewer using cmdstanr.
#
# Optimized for image size and efficiency using a multi-stage build.
# Note: Cmdstan is not included in final image, only the compiled stan model is.

ARG CMDSTAN_VERSION=2.37.0
ARG CMDSTANR_VERSION=0.9.0

# ── Stage 1: builder ──────────────────────────────────────────────────────────
FROM rocker/r-ver:4.4.2 AS builder

ARG CMDSTAN_VERSION
ENV CMDSTAN_VERSION=${CMDSTAN_VERSION}
ENV CMDSTAN=/opt/cmdstan-${CMDSTAN_VERSION}

ARG CMDSTANR_VERSION
ENV CMDSTANR_VERSION=${CMDSTANR_VERSION}

# The following part (installation of cmdstan and cmdstanr configuration) is
# adapted from https://github.com/storopoli/cmdstanr-docker by Jose Storopoli.
# See https://github.com/storopoli/cmdstanr-docker/blob/main/LICENSE for the
# license (MIT License).

RUN echo 'options(Ncpus = max(1L, parallel::detectCores() - 1L), mc.cores = max(1L, parallel::detectCores() - 1L))' >> "${R_HOME}/etc/Rprofile.site"

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential curl \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/

RUN curl -OL https://github.com/stan-dev/cmdstan/releases/download/v$CMDSTAN_VERSION/cmdstan-$CMDSTAN_VERSION.tar.gz \
  && tar xzf cmdstan-$CMDSTAN_VERSION.tar.gz \
  && rm -rf cmdstan-$CMDSTAN_VERSION.tar.gz

COPY docker/make/local $CMDSTAN/make/local

# Build cmdstan
RUN cd cmdstan-$CMDSTAN_VERSION && make -j2 build

# Install cmdstanr
RUN Rscript -e "install.packages('remotes', repos = getOption('repos'))" \
    && Rscript -e "remotes::install_github(paste0('stan-dev/cmdstanr@v', Sys.getenv('CMDSTANR_VERSION')))"

# Compile EpiSewer model — produces executable at /opt/models/EpiSewer_main
COPY inst/stan/EpiSewer_main.stan /opt/models/EpiSewer_main.stan
COPY inst/stan/functions/ /opt/models/functions/
RUN Rscript -e "cmdstanr::set_cmdstan_path('$CMDSTAN'); cmdstanr::cmdstan_model('/opt/models/EpiSewer_main.stan')"


# ── Stage 2: runtime ──────────────────────────────────────────────────────────
FROM rocker/r-ver:4.4.2 AS runtime

ARG CMDSTAN_VERSION
ENV CMDSTAN_VERSION=${CMDSTAN_VERSION}

ARG CMDSTANR_VERSION
ENV CMDSTANR_VERSION=${CMDSTANR_VERSION}

RUN echo 'options(Ncpus = max(1L, parallel::detectCores() - 1L), mc.cores = max(1L, parallel::detectCores() - 1L))' >> "${R_HOME}/etc/Rprofile.site"

# Install cmdstanr + other R packages (no compiler needed at runtime)
RUN Rscript -e "install.packages('remotes', repos = getOption('repos'))" \
    && Rscript -e "remotes::install_github(paste0('stan-dev/cmdstanr@v', Sys.getenv('CMDSTANR_VERSION')))" \
    && Rscript -e "install.packages('stringr', dependencies=TRUE)"

# Copy only compiled Stan executable from builder
COPY --from=builder /opt/models/EpiSewer_main /opt/models/EpiSewer_main

# Copy required cmdstan files from builder to use as cmdstan dummy
COPY --from=builder /opt/cmdstan-${CMDSTAN_VERSION}/makefile /opt/cmdstan-stub/makefile
COPY --from=builder /opt/cmdstan-${CMDSTAN_VERSION}/stan/lib/stan_math/lib/tbb/libtbb.so.2 \
     /opt/cmdstan-${CMDSTAN_VERSION}/stan/lib/stan_math/lib/tbb/libtbb.so.2
ENV CMDSTAN=/opt/cmdstan-stub

# Copy Stan model metadata files
COPY inst/stan/scalar_data_vars.txt /opt/models/scalar_data_vars.txt
COPY inst/stan/scalar_param_vars.txt /opt/models/scalar_param_vars.txt

# Copy R scripts
COPY R/utils_warnings.R /opt/utils_warnings.R
COPY R/utils_stan.R /opt/utils_stan.R
COPY docker/fit_EpiSewer.R /opt/fit_EpiSewer.R

WORKDIR /root
ENTRYPOINT ["Rscript"]
CMD ["--help"]
