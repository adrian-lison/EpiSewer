# Running Claude Code in a Docker Sandbox (sbx)

This folder and the kit files in `.devenv/` let Claude Code develop EpiSewer inside an isolated [Docker Sandbox](https://docs.docker.com/ai/sandboxes/). Inside the sandbox it can run R code and tests, and compile and check the Stan models.

`.devenv/` is a v3 **mixin kit**. Pass it to sbx as `.devenv/sbx/episewer-sbx`, a symlink to `..`. sbx derives the kit name from the directory name and rejects names that start with a dot, and the build context must be `.devenv/` so that `env.yml` is reachable without copying it. It layers an R + CmdStan conda env onto the built-in `claude` sandbox. It requires sbx ≥ v0.54 and a running Docker Desktop, which builds the kit image.

| File | Purpose |
|---|---|
| `.devenv/episewer-sbx.yaml` | Kit descriptor: network allowlist, dependency check at creation, agent instructions |
| `.devenv/episewer-sbx.dockerfile` | Recipe: builds `/opt/episewer-env` from `env.yml` and ships it as an overlay |
| `sbx/env-sbx-exclude.txt` | Packages from `env.yml` that are *not* installed in the sandbox |
| `sbx/env-sbx-extra.yml` | Conda packages installed *only* in the sandbox, such as lintr and styler |
| `sbx/env-post-sbx.R` | Sandbox counterpart of `env-post.R` (no vscDebugger) |
| `sbx/episewer-check-deps` | Drift guard: are all DESCRIPTION dependencies installed? |
| `sbx/context.md` | Instructions given to Claude inside the sandbox |

`env.yml` stays the single definition of the environment, and the sandbox env is **env.yml − exclusions + extras**. If a package fails to install in the sandbox (for example, no linux-aarch64 build on conda-forge):

1. Add it to `sbx/env-sbx-exclude.txt` with a reason.
2. If EpiSewer needs it, install it from CRAN or GitHub in `sbx/env-post-sbx.R`.

`episewer-check-deps` runs in the background on every sandbox start (log: `/tmp/episewer-check-deps.log`). It reports, and installs from CRAN, anything in `DESCRIPTION` that the env is missing.

## Start a sandbox

```bash
sbx login                      # once
sbx run claude --clone --name episewer \
  --kit ./.devenv/sbx/episewer-sbx \
  --kit ~/Documents/dev/sbx-kit-global-v3 \
  --cpus 6 -m 12g .
```

- The first run builds the kit image. The conda solve, the R packages and `emo` take about 10 minutes. Later runs use the cache.
- `--clone` runs Claude on an in-container clone of the repo. Your macOS build artifacts (`.conda/`, `inst/stan/EpiSewer_main`) never enter the sandbox, and the sandbox's Linux binaries never overwrite yours.
  - **Commit or stash local changes first.** Only committed state is cloned.
- Re-attach later with `sbx run --name episewer`.
- Stop or remove the sandbox with `sbx stop episewer` or `sbx rm episewer`.

## Get Claude's work back

Claude's commits are available on the host as the git remote `sandbox-episewer`:

```bash
git fetch sandbox-episewer
git log --oneline ..sandbox-episewer/<branch>
git merge sandbox-episewer/<branch>      # or cherry-pick / open a PR
```

## After changing the environment

When you edit `env.yml` or anything in `sbx/`, the next `sbx run` with a *new* sandbox rebuilds the kit. Existing sandboxes keep their image, so recreate them: `sbx rm episewer`, then run the command above.

To debug blocked downloads, run `sbx policy log episewer`.

To build the env stage locally without sbx:

```bash
docker build --target build -f .devenv/episewer-sbx.dockerfile .devenv
```
