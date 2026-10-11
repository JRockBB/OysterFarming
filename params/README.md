# Parameter files

Every model parameter lives in a TOML file here. Each top-level table is one parameter:

```toml
title = "What this parameter set is"   # optional metadata (not a table, ignored by the loader)

[p_Am]
value = 1.0
units = "J cm^-2 d^-1"
source = "Full citation, including table or equation number"
notes = "optional: conversions, assumptions"
```

Rules (from docs/AGENTS.md):

- `value`, `units` and `source` are required. `OysterFarming.load_params` refuses a file with any parameter missing them.
- Never guess a value. If no published value is found, write `source = "TODO(source)"` and leave a note. The loader refuses `TODO(source)`, so an incomplete file cannot be used by accident.
- Synthetic parameter sets are allowed only in `test/fixtures/` and their file names start with `synthetic_`.

## Files

| file | status |
|---|---|
| `deb_c_gigas_pouvreau2006.toml` | Pacific oyster DEB set (Pouvreau et al. 2006), for method testing (P2.4) |
| `deb_c_virginica.toml` | Eastern oyster DEB set; incomplete, entries marked `TODO(source)` |
