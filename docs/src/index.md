```@meta
CurrentModule = MotivicHomotopy
```

# MotivicHomotopy.jl

```@docs
MotivicHomotopy
```

## Installation

```julia
using Pkg
Pkg.add("MotivicHomotopy")
```

The development version, from the `dev` branch, can be installed with

```julia
Pkg.add(url="https://github.com/JuliaMotivicHomotopy/MotivicHomotopy.jl", rev="dev")
```

The documentation is published for both: [stable](https://JuliaMotivicHomotopy.github.io/MotivicHomotopy.jl/stable/)
is built from `main`, which holds the released code, and [dev](https://JuliaMotivicHomotopy.github.io/MotivicHomotopy.jl/dev/)
from the `dev` branch.

## Quick start

```julia
using Oscar, MotivicHomotopy

# The A¹-degree of z ↦ z² : rank 2 (complex degree), signature 0 (real degree)
S, (x,) = polynomial_ring(QQ, ["x"])
beta = global_A1_degree([x^2 + 1])
form_rank(beta), form_signature(beta)        # (2, 0)

# Local degrees sum to the global degree
S, (x, y) = polynomial_ring(QQ, ["x", "y"])
f = [x^3 - x^2 - y, y]
d1 = local_A1_degree(f, ideal(S, [x - 1, y]))
d2 = local_A1_degree(f, ideal(S, [x, y]))
is_isomorphic_form(global_A1_degree(f), gw_direct_sum(d1, d2))   # true

# Witt decomposition of a form
sum_decomposition_string(GWClass(QQ[1 2 3; 2 4 5; 3 5 6]))       # "H + <1>"

# Forms over a number field: invariants at places and the anisotropic part
# (Koprowski–Rothkegel; see the "Working with forms" page)
K, s = quadratic_field(-7)
q = diagonal_form(K, (1, 1, 1, 1))          # anisotropic at the dyadic places
anisotropic_dimension(q)                     # 4
is_isomorphic_form(anisotropic_part(q), q)   # true

P = prime_decomposition(maximal_order(K), 2)[1][1]
local_anisotropic_dimension(q, P)            # 4   (the completion at P is Q₂)
local_hilbert_symbol(-1, -1, P)              # -1
sum_decomposition_string(gw_direct_sum(q, hyperbolic_form(K)))
                                             # "H + <1> + <2> + <1> + <2>"
```

The *unstable* degrees compute numerically over ``\mathbb{C}`` when given
`HomotopyContinuation` expressions (the stable degrees are exact everywhere —
over ``\mathbb{C}`` they are just the identity form of the algebra dimension):

```julia
import HomotopyContinuation; HomotopyContinuation.@var x
global_unstable_A1_degree((x-1)*(x-2)*(x-3), (x-1)*(x-4))
```

Over ``\mathbb{R}`` the degree functions ask you to compute over ``\mathbb{Q}`` and base-change the
result. Forms over ``\mathbb{R}`` / ``\mathbb{C}`` are represented by `Matrix{Float64}` /
`Matrix{ComplexF64}` Gram matrices; pass `Float64` / `ComplexF64` where a
base field is expected.

### Certified computation over ``\mathbb{R}`` and ``\mathbb{C}``

Every function that makes decisions from float input — the unstable
``\mathbb{C}`` degrees, and the form-layer functions on float Gram matrices
(`GWClass` / `GWuClass` construction, [`form_rank`](@ref),
[`form_signature`](@ref), [`anisotropic_dimension`](@ref),
[`anisotropic_part`](@ref), [`is_isomorphic_form`](@ref), …) — accepts an
opt-in `certify = true` keyword. It certifies the result **relative to the
input read exactly**: every `Float64` / `ComplexF64` entry or coefficient is
an exact dyadic (Gaussian) rational, and with `certify = true` each discrete
output — rank, multiplicity, signature counts, determinant sign, isometry
verdict — is exactly correct for that exact input, computed over
``\mathbb{Q}`` or ``\mathbb{Q}(i)`` with no float arithmetic and no
tolerances, while returned scalars are the correctly rounded values of the
true scalars. Nothing is certified about the real-world quantity the floats
approximate: values that are close but not bit-equal are distinct.

When certification cannot return a faithful answer it errors instead of
guessing: on non-finite (`NaN` / `Inf`) input, when an exact scalar would
round to zero or infinity, and — for the ``\mathbb{C}`` degrees — when the
exact reading and the default `linear_tolerance` root clustering disagree on
a discrete decision (the error shows both readings; adjust
`linear_tolerance` until they agree to accept one explicitly). One deliberate
exception is documented at [`is_isomorphic_form`](@ref): the unstable
``k^×``-scalar comparison over ``\mathbb{R}`` / ``\mathbb{C}`` keeps its
tolerance even under `certify`, because float scalars carry rounding drift
from however they were produced. With `certify = false` (the default)
every computation takes the default path, unchanged by the keyword; on exact
base fields the keyword is accepted and is a no-op.

The full reference is split by task: [Working with forms](@ref), [Arithmetic methods](@ref), and [Computing A¹-Brouwer degrees](@ref).

## Conventions

- **`==` is literal Gram-matrix equality**; mathematical equality in
  ``\text{GW}(k)`` is [`is_isomorphic_form`](@ref).
- Accessors and operations whose natural names would collide with Oscar
  exports carry a `gw_` prefix ([`gw_matrix`](@ref), not `gram_matrix`;
  [`gw_witt_index`](@ref), not `witt_index`).
- Degree Gram matrices are expressed in a standard-monomial basis; an
  equivalent form in another basis is the same class, so compare with
  [`is_isomorphic_form`](@ref) rather than entrywise.

## Version history

MotivicHomotopy.jl began as a Julia port of the Macaulay2 package
A1BrouwerDegrees. Its initial release consolidates the functionality of that
package's versions 1.1 and 2.0 into a single starting point; subsequent
releases listed below build on it.

- **v0.1.0** — based on versions 1.1 and 2.0 of the Macaulay2 package
  A1BrouwerDegrees. Implements the computation of local and global 𝔸¹-Brouwer
  degrees and of Grothendieck–Witt classes and their invariants; the
  computation of unstable local and global 𝔸¹-Brouwer degrees and manipulation
  of the unstable Grothendieck–Witt group; and the extension of
  Grothendieck–Witt class manipulations over fields to finite étale algebras
  over fields, including transfers along finite étale extensions. Developed by
  S. Atherton, N. Borisov, T. Brazelton, S. Dutta, F. Espino, T. Hagedorn,
  Z. Han, J. Lopez Garcia, J. Louwsma, Y. Luo, G. Ong, R. Sagayaraj, and
  A. Tawfeek.
- **v1.0.0** — first stable release. Raises the minimum supported versions to
  Oscar 1.8 and HomotopyContinuation 2.22, and reorganizes the manual into
  task-based pages. The API is that of v0.1.0, with the following corrections,
  the first of them breaking.
  - For a finite étale algebra ``L/k``, [`transfer_gw`](@ref) now returns the
    transfer along the trace ``\text{Tr}_{L/k}``, which sends a form of rank
    ``n`` over ``L`` to a form of rank ``n[L : k]``; v0.1.0 returned a form of
    rank ``n``, the traces of the entries of a diagonalization.
  - Over ``\mathbb{Q}``, [`anisotropic_part`](@ref), and with it
    [`sum_decomposition`](@ref) and [`sum_decomposition_string`](@ref), no
    longer returns a form of the wrong rank for some forms of anisotropic
    dimension at least 3.
  - The unstable degrees over ``\mathbb{C}`` ([`global_unstable_A1_degree`](@ref)
    and [`local_unstable_A1_degree`](@ref) with HomotopyContinuation input)
    now normalize by the leading coefficient of the numerator, as the exact
    computation does; for a numerator that is not monic the scalar was off by
    a power of that coefficient.
  - Over an étale algebra that is not a field, [`diagonalize_via_congruence`](@ref),
    and with it [`diagonal_class`](@ref), no longer scales by a pivot that is
    a zero divisor, which made the result degenerate on a factor of the
    algebra.
  - Over a quotient ring built with a monomial ordering other than the
    default, [`multiplication_matrix`](@ref), [`algebra_trace`](@ref) and
    [`algebra_norm`](@ref) no longer return wrong values.
  - [`sum_decomposition_string`](@ref) works over finite fields that are not
    prime fields. [`divisorial_sum`](@ref) accepts classes over distinct
    finite fields of the same order, as [`gw_direct_sum`](@ref) does.
    [`hyperbolic_form`](@ref) and [`pfister_form`](@ref) accept finite étale
    algebras, and [`diagonal_form`](@ref) over `Float64` accepts rationals.
    [`local_algebra_basis`](@ref) rejects an ideal that is not prime, and
    [`padic_valuation`](@ref) a modulus that is not prime.
