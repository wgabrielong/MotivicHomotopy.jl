# Magma differential oracle

Generated corpus giving ≥ 100 Magma-checked cases for each of 27 exported
functions over the 13-field matrix.

## Running it (someone with Magma)

Run each `oracle_ext_<i>.m` in Magma — any order, any subset, parallel
sessions are fine, every chunk is self-contained — and concatenate the
output into one text file. Each printed line is one datum,
`D|record|key|value`. The scripts use only core intrinsics
(`NumberField`, `Decomposition`, `HilbertSymbol`, `IsSquare`,
`RealPlaces` / `Evaluate` / `Sign`, `Matrix` / `Determinant`);
`OrthogonalizeGram` is wrapped in try/catch and degrades to a `skip` line.
Expected total runtime is well under two hours; the online calculator's
~2-minute limit means running chunks one at a time, possibly trimmed.

Then, on our side (from the repository root, after
`julia --project=. -e 'using Pkg; Pkg.instantiate()'` once per clone):

```
julia --project=. differential/magma-oracle/compare.jl <results.txt>
```

Partial results files are fine: absent records are skipped and counted,
and the quota is enforced only when the file is complete.

## Reading the output

- The package's own outputs (anisotropic parts, discriminant
  representatives, divisorial-sum scalars, diagonalizations) appear in the
  scripts as literals; Magma prints square-class and symbol verdicts about
  them, and `compare_core.jl` draws the conclusions (each theorem cited at
  its point of use).
- A key may repeat within a record — once per prime ideal above the stated
  rational prime, or once per real place. The comparator matches repeats
  as multisets; booleans that pair two forms at the same place are
  computed inside Magma and printed as a single line.
- Records whose data would involve primes above 10⁹ or element literals
  beyond 300 characters are excluded at generation time (the count is
  printed); those shapes are covered by the local Hecke tier.
- `results-hecke-simulated.txt` is the same schema computed with Hecke
  primitives. It validates the parse/compare pipeline and is **not** Magma
  output; for the `hs@*` / `hasse@*` keys it is only a plumbing check,
  since the package's finite-place symbol delegates to Hecke.

## Regenerating

```
julia --project=. differential/magma-oracle/generate_all.jl
```

rebuilds the corpus, the chunk scripts, and the simulated results
deterministically, then runs the comparison against the simulation and
enforces the quota.

## Findings (resolved)

Corpus generation surfaced a form on which `anisotropic_part` failed in
[KR21] Algorithm 8:

    K = Q(sqrt 5),  u = 1000003,  v = 1 - 2 sqrt 5
    q = <1, -u, v, -u*v>

The F₂ system stayed unsolvable through eight rounds of the blind
smallest-norm enlargement of S (|S| = 35). Diagnosis: with K(√d)/K
ramified at eleven primes, the symbol patterns reachable from E_S form a
small subgroup (genus theory), so the blind strategy needs on the order
of 2^(#ramified) candidate primes — beyond any fixed retry bound. The
E_S basis and the dimension chain were verified sound (independent
Hecke-`hasse_invariant` referee).

Fixed by the closing-prime enlargement in step (5g) of
`_kr_algorithm_8_loop` (`src/nf_anisotropic_part.jl`, a flagged deviation
from the printed step): a principal prime whose generator's sign/symbol
column provably closes the system is searched for with symbols alone
before the blind batch is used as fallback. The reproducer is pinned in
`test/test_number_fields_aniso.jl` and its record now generates
normally; the fix lives on `number-fields` and is carried here.

## Covered functions and exclusions

Covered (≥ 100 each, enforced on complete runs): the 27 functions in
`EXT_QUOTA_FUNCTIONS` in `compare_core.jl`. Not Magma-checkable, covered
by the Hecke tier instead: `GWClass` validation, `form_rank`,
`stable_part`, `sum_decomposition_string`.
