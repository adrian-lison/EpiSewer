# Recipe for the episewer-sbx mixin (see episewer-sbx.yaml and sbx/README.md).
#
# Builds the EpiSewer R/CmdStan conda env from env.yml (minus the exclusions in
# sbx/env-sbx-exclude.txt, plus sbx/env-sbx-extra.yml) in a build stage, then
# ships only /opt/episewer-env as an overlay onto the Claude sandbox image.
# The conda env is self-contained (R, compilers, CmdStan), so no apt packages
# are needed on the hardened base image.

ARG MICROMAMBA_IMAGE=mambaorg/micromamba:2.0.5

# ── build stage: create the conda env ─────────────────────────────────────────
FROM ${MICROMAMBA_IMAGE} AS build
USER root
ENV PREFIX=/opt/episewer-env
WORKDIR /tmp/devenv

COPY env.yml sbx/env-sbx-exclude.txt sbx/env-sbx-extra.yml ./

# env.yml minus excluded packages -> env.filtered.yml
RUN set -eu; \
    pkgs="$(sed -e 's/#.*//' -e 's/[[:space:]]//g' env-sbx-exclude.txt | grep -v '^$' | paste -sd'|' -)"; \
    if [ -n "$pkgs" ]; then \
      grep -Ev "^[[:space:]]*-[[:space:]]*(${pkgs})([<>=!~ ].*)?[[:space:]]*$" env.yml > env.filtered.yml; \
    else cp env.yml env.filtered.yml; fi; \
    echo "Excluded: ${pkgs}"

RUN micromamba create -y -p "$PREFIX" -f env.filtered.yml -f env-sbx-extra.yml \
    && micromamba clean -a -y

# Bake the conda activation environment (compilers etc.) into fixed config,
# because conda activation scripts do not run inside the sandbox.
COPY sbx/env-post-sbx.R sbx/episewer-check-deps ./
RUN set -e; \
    eval "$(micromamba shell hook -s posix)"; micromamba activate "$PREFIX"; \
    echo "CMDSTAN=$CMDSTAN  CC=${CC:-}  CXX=${CXX:-}"; \
    # CmdStan: pin the conda toolchain (conda's make/local already sets the TBB
    # options; CC/CXX normally come from activation scripts, which won't run)
    printf 'CC=%s\nCXX=%s\n' "$(command -v "$CC")" "$(command -v "$CXX")" >> "$CMDSTAN/make/local"; \
    cat "$CMDSTAN/make/local"; \
    # R site profile: CRAN mirror and parallelism
    printf '%s\n' \
      'local({' \
      '  options(repos = c(CRAN = "https://cloud.r-project.org"))' \
      '  n <- max(1L, parallel::detectCores() - 1L)' \
      '  options(Ncpus = n, mc.cores = n)' \
      '})' >> "$(R RHOME)/etc/Rprofile.site"; \
    Rscript env-post-sbx.R; \
    install -m 0755 episewer-check-deps "$PREFIX/bin/episewer-check-deps"; \
    # smoke test
    Rscript -e 'cat("CmdStan", as.character(cmdstanr::cmdstan_version()), "\n")'; \
    chown -R 1000:1000 "$PREFIX"

# ── overlay: only the env lands on the Claude sandbox image ───────────────────
FROM scratch
COPY --from=build /opt/episewer-env /opt/episewer-env
ENV PATH=/opt/episewer-env/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    CMDSTAN=/opt/episewer-env/bin/cmdstan \
    TZ=UTC \
    TZDIR=/opt/episewer-env/share/zoneinfo
# No ENTRYPOINT/CMD/USER: the claude workload's launch command stays.
