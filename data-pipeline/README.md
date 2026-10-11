# oyster-data

Python data pipeline for OysterFarming. It downloads and cleans public data and writes the Arrow files described in `docs/PLAN.md` ("Data contract"). It contains no model logic.

```
uv run pytest                          # run tests (no network access)
uv run python scripts/make_fixtures.py # regenerate synthetic contract fixtures in data/fixtures/
```
