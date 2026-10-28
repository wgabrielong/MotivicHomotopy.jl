# Hecke-differential suite for the number-field layer

Runs locally, no external CAS:

```
julia --project=. differential/hecke-differential/run.jl
```

Every exported function whose number-field support the `number-fields`
branch added or changed accumulates at least 100 differential cases; the
runner enforces the quota on full runs and exits nonzero on any divergence
or shortfall. The corpus is deterministic (fixed-seed LCG) and a full run
takes about a minute.

Oracles: Hecke's `is_isometric` (global and local), local `is_isotropic`,
`hasse_invariant`, and `discriminant`, each pinned on hand-known values in
`gates.jl` before any suite runs — a gate failure aborts the run. Hecke's
global `is_isotropic` is broken over number fields (definite forms are
reported isotropic) and is never called; global isotropy is assembled from
local verdicts (`oracle_isotropic_global` in `oracles.jl`), with binary
forms decided by `is_square`. The suites also check identities the
implementation does not use in that form: Hilbert reciprocity, the Hasse
sum formula, the Pfister dichotomy, certified-numeric signatures, and a
recomputation of the divisorial-sum scalar. `local_hilbert_symbol`
delegates to Hecke at finite places, so it is checked through those
identities rather than by direct comparison. `transfer_gw` computes with
Hecke's trace, so it is checked against routes that do not: Hecke's
restriction of scalars of the quadratic space over each field of the
matrix, and the étale algebra presenting the field, whose transfer is built
from multiplication matrices — over each field, a relative quadratic
extension of it, and a non-simple field.

## Layout

| file | contents |
|---|---|
| `run.jl` | runner: gates → suites → per-function case table |
| `harness.jl` | case counting, divergence reporting, the fixed-seed LCG |
| `gates.jl` | Hecke convention gates |
| `oracles.jl` | Hecke wrappers and Hasse–Minkowski global isotropy |
| `generators.jl` | field matrix with cost tiers; element/form/dense pools |
| `suite_*.jl` | one file per function group |

`run.jl symbols decomposition` runs a subset; `DIFF_MAX_TIER=1` restricts
to the cheaper fields. Quotas apply to full runs only.

Every divergence prints one line with the exported function, the failed
identity, and the reproducing data (field, diagonal, place, expected vs
got).
