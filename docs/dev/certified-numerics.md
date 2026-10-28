# Certified numerical computation — design narrative

This note records what the certified-computation update adds, the choices made
while building it, and how it was validated. The user-facing reference is the
*Certified computation over ℝ and ℂ* section of the manual's front page and
the docstrings of the affected functions; this file explains *why* the update
looks the way it does.

## What the update adds

Every function that makes a decision from floating-point input now accepts an
opt-in keyword `certify::Bool = false`:

- **The ℂ unstable degrees** (`global_unstable_A1_degree`,
  `local_unstable_A1_degree` on `HomotopyContinuation` expressions). The
  default path finds roots numerically via a companion matrix and makes six
  `linear_tolerance` decisions that determine the rank and multiplicity of the
  output form.
- **The ℝ/ℂ float-Gram form layer** (`GWClass` / `GWuClass` construction,
  `form_rank`, `form_signature`, `anisotropic_dimension`, `gw_witt_index`,
  `is_anisotropic_form` / `is_isotropic_form`, `anisotropic_part`,
  `is_isomorphic_form`). The default paths make bare float `sign` / `iszero` /
  SVD-rank decisions, including the determinant-sign well-definedness gate of
  `GWuClass` over ℝ.

With `certify = true` the result is certified **relative to the input read
exactly**: every `Float64` / `ComplexF64` entry or coefficient is an exact
dyadic (Gaussian) rational, every discrete output — rank, multiplicity,
signature counts, determinant sign, pointedness, isometry verdict — is exactly
correct for that exact input, and returned scalars are the correctly rounded
values of the true scalars. Nothing is certified about the real-world quantity
the floats approximate: values that are close but not bit-equal are distinct.
On exact base fields (ℚ, 𝔽_q, number fields) the keyword is accepted and is a
no-op, because those paths are always exact.

## Choices made, and why

**Exact rational arithmetic, not interval arithmetic.** The original request
suggested certifying determinant signs with interval arithmetic. The
implementation uses exact ℚ / ℚ(i) arithmetic via Oscar instead, because it
strictly dominates here: every decision point takes input that is exactly
rational (a float *is* a rational number), a ball straddling zero can never
decide `det == 0` and would need an exact fallback anyway, and the Gram
matrices and polynomials in scope are small enough that exact arithmetic is
immediate (the certified signature of the 12×12 Hilbert matrix takes
milliseconds). This also added **no dependencies** — Oscar was already a direct
dependency. The one place the package does use ball arithmetic,
`src/nf_arithmetic.jl`, is the case where exactness is impossible (real
embedding images are irrational); that contrast is deliberate.

**No root-finding in the certified ℂ path.** The float path's scalar is a
product over numerically-found roots. Read exactly, that product *is* a
resultant: the tolerance cancellation becomes the exact `gcd` in ℚ(i)[x], rank
and multiplicity become exact degrees and deflation counts, and the scalar
becomes an exact resultant / evaluation formula. Oscar's monic gcd keeps the
leading coefficients of numerator and denominator intact, which is exactly what
makes the exact formula match the float path, which pairs the monic numerator
with the denominator divided by the original leading coefficient of the
numerator after cancellation.

**Opt-in keyword, defaults unchanged.** The float paths of this package
follow the Macaulay2 package's behavior, including its exact-`== 0` pivot tests
on floats, apart from the corrections listed for v1.0.0 in the manual's version
history. Certification therefore had
to be additive: `certify = false` (the default) is unchanged everywhere, and
every certified branch carries a comment marking it as an extension beyond M2.

**Error loudly, never guess.** Certification refuses in three situations
rather than returning something unfaithful: non-finite (`NaN` / `Inf`) input;
an exact scalar that would round to zero or infinity in the return type; and —
for the ℂ degrees — a *divergence* between the exact reading and the default
tolerance clustering on any discrete decision (rank, is-root, multiplicity).
The divergence error shows both readings; passing a smaller `linear_tolerance`
(to separate nearby roots) or a larger one (to merge a numerical root cluster
— a fivefold root's numerical cluster is ~10⁻³ wide, far beyond the default
10⁻⁶) accepts one reading explicitly. Pointedness needs no separate guard:
cancellation drops numerator and denominator degrees together in both
readings, so once the ranks agree the pointedness verdicts agree.

**One deliberate non-certification.** In `is_isomorphic_form` on unstable
classes, the ``k^×``-scalar comparison over ℝ/ℂ keeps its `linear_tolerance`
even under `certify`. Two float scalars each carry rounding drift from however
they were produced (e.g. divisorial sums), so exact equality would declare
mathematically equal classes distinct. `certify` governs the stable-part
comparison only; the docstring records the caveat.

**Versioning.** The version number is unchanged: this branch builds on the
unreleased number-fields work, so the minor bump and its version-history entry
belong to whichever release ships first.

## Validation

- **Unit tests** (`test/test_certified.jl`): the Macaulay2-captured ℂ cases
  rerun with `certify = true` and asserted with `==` against exact values
  (e.g. the scalar `-1.9999999999999987…` becomes exactly `-2.0`); divergence
  guard tests in both directions with the `linear_tolerance` escape hatch;
  overflow / underflow refusal; certified constructors on determinants below
  the float underflow threshold; and the default behavior locked alongside as
  bug-compat. The 12×12 Hilbert matrix is a genuine divergence witness: the
  default SVD rank is 11, the certified exact rank is 12 (Hilbert matrices are
  provably nonsingular).
- **Differential tests against Macaulay2** (the cross-system suite from the
  separate differential-tests repository, oracle pattern: inputs and expected
  outputs are generated once with the Macaulay2 package `A1BrouwerDegrees`,
  then the Julia package is run on the identical inputs and compared). Both
  runners were rerun against the pre-generated Macaulay2 oracles with this
  package substituted for the original port, after the update:
  - *Exact-field runner* (`run_diff.jl`): **8383 assertions over 8300 oracle
    cases in 83 case families, all passing with zero M2 divergences.** The
    families cover, over ℚ: the scalar invariants (`form_rank`,
    `form_signature`, `integral_discriminant`, `anisotropic_dimension`,
    `gw_witt_index`, isotropy predicates, `hasse_witt_invariant`,
    `anisotropic_dimension_qqp`, Hilbert symbols, `padic_valuation`) and the
    class-valued operations (`diagonal_class`, `diagonalize_via_congruence`,
    `diagonal_entries`, direct sum, tensor product, `is_isomorphic_form`, the
    form builders, `anisotropic_part`, `sum_decomposition` and its string);
    the same operation set again over 𝔽_q (odd q) and over ℝ and ℂ as float
    Gram matrices — i.e. exactly the float form layer the `certify` keyword
    was added to, exercised through its **default** paths; étale-algebra
    functions (`multiplication_matrix`, `algebra_trace`, `algebra_norm`,
    `transfer_gw`); the stable degrees (`global_A1_degree`,
    `local_A1_degree`, `local_algebra_basis`); and the unstable layer
    (`global_unstable_A1_degree`, `local_unstable_A1_degree`,
    `divisorial_sum`, `GWuClass` construction, the unstable form builders).
    Form-valued outputs are compared up to isomorphism via a canonical
    invariant computed independently in each system (rank; signature;
    discriminant; the set of primes with Hasse–Witt invariant −1), never by
    matrix equality, and known upstream M2 bugs are tracked as `BUG:` cases.
  - *Numerical ℂ runner* (`run_num.jl`): **400 assertions over 200
    rational-function cases**, comparing the numerical unstable degrees
    (`num_global_unstable`, `num_local_unstable`) against the Macaulay2 CC
    oracle — ranks exactly, scalars to relative tolerance 10⁻⁵ (both systems
    are inexact on this path).

  All of it passes, confirming the default (`certify = false`) behavior of
  every touched function is unchanged.
- The full in-repo suite passes (2850 tests) and the Documenter build is
  clean.

These runs were made at a74e014, before the v1.0.0 corrections. Since then
`transfer_gw` and the normalization of the ℂ unstable degrees by the leading
coefficient of the numerator differ from the Macaulay2 package by design, as
the manual's version history records, and the certified ℂ formulas follow the
corrected normalization.

## Viewing the documentation locally

From the repository root:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

then open `docs/build/index.html` in a browser (on macOS:
`open docs/build/index.html`). The certified-computation section is on the
front page (*Home*), and the `certify` keyword is documented per-function on
the *Working with forms* and *Computing A¹-Brouwer degrees* pages. The first
command is only needed once per checkout; rebuild after edits with the second
command alone.
