# OysterFarming

A farm-scale simulation and decision model for eastern oyster (*Crassostrea virginica*) aquaculture in Maine.

- Project overview and research question: [docs/README.md](docs/README.md)
- Work plan and data contract: [docs/PLAN.md](docs/PLAN.md)
- Tech stack: [docs/TECH_STACK.md](docs/TECH_STACK.md)
- Rules for coding agents: [docs/AGENTS.md](docs/AGENTS.md)

## Quick start

Julia model (this repository root is the Julia package):

```
julia --project -e "using Pkg; Pkg.instantiate(); Pkg.test()"
```

Python data pipeline:

```
cd data-pipeline
uv run pytest
```
