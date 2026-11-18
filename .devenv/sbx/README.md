# Running Claude Code in a Docker Sandbox (sbx)

This folder and the kit files in `.devenv/` let Claude Code develop EpiSewer inside an isolated [Docker Sandbox](https://docs.docker.com/ai/sandboxes/). Inside the sandbox it can run R code and tests, and compile and check the Stan models.

`.devenv/` is a v3 **mixin kit**. sbx derives the kit name from the directory name and rejects names that start with a dot, and the kit directory is the build context, so it must contain `env.yml`. `sbx/sbx-run` therefore mirrors `.devenv/` into `~/.cache/sbx-kits/episewer-sbx/` and passes that directory to sbx. There is deliberately no symlink in the repo: a symlink pointing at its parent creates a cycle, which breaks `pak`/`R CMD build` (e.g. the pkgdown GitHub Action), even though `.devenv` is `.Rbuildignore`d. It layers an R + CmdStan conda env onto the built-in `claude` sandbox. It requires sbx ≥ v0.54 and a running Docker Desktop, which builds the kit image.

| File | Purpose |
|---|---|
| `.devenv/episewer-sbx.yaml` | Kit descriptor: network allowlist, dependency check at creation, agent instructions |
| `.devenv/episewer-sbx.dockerfile` | Recipe: builds `/opt/episewer-env` from `env.yml` and ships it as an overlay |
| `sbx/env-sbx-exclude.txt` | Packages from `env.yml` that are *not* installed in the sandbox |
| `sbx/env-sbx-extra.yml` | Conda packages installed *only* in the sandbox, such as lintr and styler |
| `sbx/env-post-sbx.R` | Sandbox counterpart of `env-post.R` (no vscDebugger) |
| `sbx/episewer-check-deps` | Drift guard: are all DESCRIPTION dependencies installed? |
| `sbx/episewer-site` | Builds the pkgdown site (or a single vignette with `--article NAME`) and serves it on sandbox port 8000 |
| `sbx/context.md` | Instructions given to Claude inside the sandbox |
| `sbx/sbx-run` | Host-side launcher: stages the kit under a valid name and runs `sbx run` |

`env.yml` stays the single definition of the environment, and the sandbox env is **env.yml − exclusions + extras**. If a package fails to install in the sandbox (for example, no linux-aarch64 build on conda-forge):

1. Add it to `sbx/env-sbx-exclude.txt` with a reason.
2. If EpiSewer needs it, install it from CRAN or GitHub in `sbx/env-post-sbx.R`.

`episewer-check-deps` runs in the background on every sandbox start (log: `/tmp/episewer-check-deps.log`). It reports, and installs from CRAN, anything in `DESCRIPTION` that the env is missing.

## Start a sandbox

```bash
sbx login                      # once
.devenv/sbx/sbx-run
```

This runs `sbx run claude --clone --name episewer --kit ~/.cache/sbx-kits/episewer-sbx --kit ~/Documents/dev/sbx-kit-global-v3 --cpus 6 -m 12g -p 8000:8000 <repo>`. Extra `sbx run` options are passed through. Set `EPISEWER_SBX_NAME` to use another sandbox name, and `EPISEWER_SBX_GLOBAL_KIT` for another global kit. For a custom `sbx` command, `.devenv/sbx/sbx-run --stage-only` refreshes the staged kit and prints its path.

- The first run builds the kit image. The conda solve, the R packages and `emo` take about 10 minutes. Later runs use the cache.
- `--clone` runs Claude on an in-container clone of the repo. Your macOS build artifacts (`.conda/`, `inst/stan/EpiSewer_main`) never enter the sandbox, and the sandbox's Linux binaries never overwrite yours.
  - **Commit or stash local changes first.** Only committed state is cloned.
- Re-attach later with `sbx run --name episewer`.
- Stop or remove the sandbox with `sbx stop episewer` or `sbx rm episewer`.

## View the pkgdown site (README, vignettes) in your browser

Ask Claude to build and serve the site. It runs `episewer-site`, which builds `docs/` and serves it on sandbox port 8000. Then open <http://localhost:8000>.

To look at a single vignette, ask Claude to expose just that one. For example, "serve only the detailed-example vignette" runs `episewer-site --article detailed-example`. Only that vignette is built (about 15 s instead of 40 s) and served, and <http://localhost:8000> redirects to it.

`-p 8000:8000` in the command above publishes the port when the sandbox is created. For a sandbox created without it, publish the port later. This works while the sandbox is running:

```bash
sbx port publish episewer 8000:8000      # list ports: sbx port ls episewer
```

### Temporary vignettes

You can ask Claude for a **temporary** (or **local**) vignette, for example "make a temporary vignette showing the new plot options and serve it". Claude writes it to `vignettes/local/<name>.Rmd` and serves only that vignette.

- `vignettes/local/` is gitignored and `.Rbuildignore`d, so temporary vignettes are never committed or shipped, even when you ask Claude to commit its changes.
- In `_pkgdown.yml`, the `internal` articles section (`starts_with("local/")`) lets full site builds succeed while a temporary vignette exists, without listing it on the site.
- To keep a temporary vignette, tell Claude to make it **official** or **public**. Claude moves it to `vignettes/`, adds it to the `articles:` index in `_pkgdown.yml`, checks the full site build and commits both.
- In clone mode, temporary vignettes exist only inside the sandbox. They disappear when you remove the sandbox (`sbx rm`).

## Get Claude's work back

Claude's commits are available on the host as the git remote `sandbox-episewer`:

```bash
git fetch sandbox-episewer
git log --oneline ..sandbox-episewer/<branch>
git merge sandbox-episewer/<branch>      # or cherry-pick / open a PR
```

## After changing the environment

When you edit `env.yml` or anything in `sbx/`, the next `sbx-run` for a *new* sandbox restages and rebuilds the kit. Existing sandboxes keep their image, so recreate them: `sbx rm episewer`, then run the command above.

To debug blocked downloads, run `sbx policy log episewer`.

To build the env stage locally without sbx:

```bash
docker build --target build -f .devenv/episewer-sbx.dockerfile .devenv
```
