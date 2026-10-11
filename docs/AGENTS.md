# Rules for coding agents

Read README.md, then PLAN.md, then this file before starting any task. Work on one PLAN.md task ID at a time and name it in commit messages (e.g. `P2.2: standard DEB model`).

## Scientific integrity (most important)

1. **Never invent parameter values.** Every model parameter lives in `params/*.toml` with `value`, `units` and `source`. If you cannot find a published value, stop and leave a `TODO(source)` for the human. Do not guess a "reasonable" number.
2. **Never invent data.** Synthetic data is allowed only in tests and only when clearly labeled as synthetic in the file or function name.
3. **Cite equations.** Every function implementing a model equation has a docstring with the equation in plain text or LaTeX and the paper it comes from.
4. **Units in names or docs.** Use explicit unit suffixes for variables where ambiguity is possible (`temp_c`, `shell_height_mm`, `dt_days`). Temperature in °C at boundaries; convert to Kelvin only inside the Arrhenius function.

## Language boundaries

- `data-pipeline/` is Python. It downloads, cleans and writes contract files. It contains no model logic.
- The repository root is the Julia package `OysterFarming` (`src/`, `test/`). It reads contract files and runs the model. It does no web requests.
- The only interface between them is the data contract in PLAN.md. If a task needs a schema change, update PLAN.md, the writer, the reader and both validators in the same change.

## Julia conventions

- Plain functions and structs. Avoid macros and metaprogramming unless a package requires them.
- Make structs concretely typed or parametric; avoid `Any` fields (performance and clarity).
- Model functions are pure: inputs in, outputs out, no global state. Randomness takes an explicit `rng` argument.
- Pin package versions in `[compat]`. If an API you remember fails, check the installed version's documentation rather than guessing.
- Tests in `test/` using the `Test` stdlib, one file per source file.

## Python conventions

- Managed with uv. Type hints on public functions.
- Fetchers save raw responses to `data/raw/<source>/` and never silently overwrite.
- Tests use saved sample responses in `data-pipeline/tests/fixtures/`. Tests never call live servers.

## Testing requirements

- Every new function gets at least one test. Every model component gets the property tests listed for it in PLAN.md.
- Run the full test suite for the language you touched before finishing a task.
- A task is done only when its "Done when" line in PLAN.md is satisfied. Say explicitly if it is not.

## Data and files

- `data/raw/` is gitignored. `data/processed/` holds contract files; commit them only if small (under ~10 MB), otherwise document how to regenerate them.
- Scripts in `scripts/` are named after their task ID (`p2_4_gigas_check.jl`) and produce figures into `outputs/<task_id>/`.

## When to stop and ask

- A needed parameter or dataset has no clear source.
- A task's "Done when" cannot be met as written.
- A design choice would change the data contract or the thesis scope.
