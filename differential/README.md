# Differential tests

Two tiers, independent of the package's unit tests (`test/`). One-time
setup, from the repository root (all commands below run from the root):

```
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

- **`hecke-differential/`** — runs locally, no external CAS:

  ```
  julia --project=. differential/hecke-differential/run.jl
  ```

- **`magma-oracle/`** — the cross-CAS tier: a generated corpus with a
  ≥ 100-cases-per-exported-function quota, shipped as chunked Magma
  scripts plus a Julia comparator. See `magma-oracle/README.md` for run
  instructions.

The package's finite-place Hilbert symbol delegates to Hecke, so for that
primitive the Magma tier is the only fully external authority.
