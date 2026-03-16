FROM rocker/r-ver:4.4.2

# nvironment variables
ENV CMDSTAN_VERSION=2.37.0
ENV CMDSTAN=/opt/cmdstan-$CMDSTAN_VERSION
ENV CMDSTANR_VERSION=0.9.0

# set Rprofile
RUN echo 'options(Ncpus = max(1L, parallel::detectCores() - 1L), mc.cores = max(1L, parallel::detectCores() - 1L))' >> "${R_HOME}/etc/Rprofile.site"

# install required packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential curl \
    && rm -rf /var/lib/apt/lists/*

# set workdir for /opt/cmdstan-CMDSTAN_VERSION
WORKDIR /opt/

# download and extract cmdstan based on CMDSTAN_VERSION from github
RUN curl -OL https://github.com/stan-dev/cmdstan/releases/download/v$CMDSTAN_VERSION/cmdstan-$CMDSTAN_VERSION.tar.gz \
  && tar xzf cmdstan-$CMDSTAN_VERSION.tar.gz \
  && rm -rf cmdstan-$CMDSTAN_VERSION.tar.gz

# copy the make/local to CMDSTAN dir
COPY docker/make/local $CMDSTAN/make/local

# build cmdstan using 2 threads
RUN cd cmdstan-$CMDSTAN_VERSION && make -j2 build

# install cmdstanr
RUN Rscript -e "install.packages('remotes', repos = getOption('repos'))" \
    && Rscript -e "remotes::install_github(paste0('stan-dev/cmdstanr@v', Sys.getenv('CMDSTANR_VERSION')))"

# install further R packages
RUN Rscript -e "install.packages('stringr', dependencies=TRUE)"

# compile EpiSewer stan model
COPY inst/stan/EpiSewer_main.stan /opt/models/EpiSewer_main.stan
COPY inst/stan/functions/ /opt/models/functions/
COPY inst/stan/scalar_data_vars.txt /opt/models/scalar_data_vars.txt
COPY inst/stan/scalar_param_vars.txt /opt/models/scalar_param_vars.txt
RUN Rscript -e "cmdstanr::set_cmdstan_path('$CMDSTAN'); cmdstanr::cmdstan_model('/opt/models/EpiSewer_main.stan')"

# copy additional R scripts to /opt
COPY R/utils_warnings.R /opt/utils_warnings.R
COPY R/utils_stan.R /opt/utils_stan.R
COPY docker/fit_EpiSewer.R /opt/fit_EpiSewer.R

# Final preparations
WORKDIR /root
ENTRYPOINT ["Rscript"]
CMD ["--help"]
