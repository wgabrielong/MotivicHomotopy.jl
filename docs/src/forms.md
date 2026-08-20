```@meta
CurrentModule = MotivicHomotopy
```

# Working with forms

Everything for building and computing with Grothendieck–Witt classes:
constructors, the ring operations, simplified representatives, the Witt
decomposition and anisotropic part, invariants, and the isomorphism test.

Supported base fields: ``\mathbb{Q}``, finite fields of odd characteristic,
any absolute simple number field (`AbsSimpleNumField` — the fields produced
by `quadratic_field`, `number_field(f)`, and
`Hecke.rationals_as_number_field()`), and ``\mathbb{R}`` / ``\mathbb{C}``
via `Matrix{Float64}` / `Matrix{ComplexF64}` Gram matrices. Decisions on
float Gram matrices can be made exactly from the entries read as exact
rationals with the opt-in `certify = true` keyword (see *Certified
computation* on the manual's front page). Non-simple
fields (`number_field([f, g])`) are rejected with an error naming the fix:
convert with `simple_extension(K)` first.

Over a number field the anisotropic part runs the Koprowski–Rothkegel
algorithms [KR23] on top of the Koprowski–Czogała local theory [KC18]
(the reference implementation of the same algorithms is the CQF Magma
package [Kop20]). Two things to know when working over ``K``:

- **Anisotropic parts are GRH-conditional** with the default settings:
  the pipeline uses Hecke's `class_group` and `sunit_group_fac_elem`,
  whose defaults assume GRH. Invariants that avoid the class group —
  rank, signatures, discriminants, Hilbert symbols, Hasse invariants,
  local and global anisotropic dimensions, isometry testing — are
  unconditional.
- **Representatives are size-reduced, never normalised.** Over
  ``\mathbb{Q}`` every square class has a canonical squarefree integer
  representative; over ``K`` none exists (the class group obstructs
  exponent reduction and the unit group leaves the associate open), so
  the output entries of [`anisotropic_part`](@ref),
  [`sum_decomposition`](@ref), and [`diagonal_class`](@ref) are merely
  kept small within their square classes. Two isomorphic forms may still
  print differently — compare with [`is_isomorphic_form`](@ref), never by
  looking at entries. (Internally the algorithms keep factored `FacElem`
  arithmetic, and the real embeddings travel coupled to certified
  isolating intervals, so none of this affects correctness — only what
  you see printed.)

## Classes and ring operations

Stable and unstable classes, accessors, orthogonal sum and tensor product
(with the integer multiple `n * beta` and tensor power `beta^n`), base
change ``\text{GW}(\mathbb{Q}) → \text{GW}(K)``, and the divisorial sum.
GW classes form a **semiring**: sums and products exist, formal
differences do not, and ``⟨-1⟩ ⊗ β`` is the opposite form rather than an
additive inverse.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["gw_classes.jl", "gwu_classes.jl"]
```

## Building and decomposing forms

Standard forms, diagonalization by congruence, simplified representatives,
and the Witt (sum) decomposition and anisotropic part.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["building_forms.jl", "matrix_methods.jl", "simplified_representatives.jl", "decomposition.jl"]
```

## Invariants and classification

Rank, signatures, discriminants, Hasse–Witt invariants, anisotropic
dimension and Witt index, the isotropy predicates, and the isomorphism
test. Conventions: [`signed_discriminant`](@ref) is
``(-1)^{n(n-1)/2} a_1 ⋯ a_n`` ([KR23] §2), compared over ``K`` only up to
squares; [`hasse_witt_invariant`](@ref) computes the **Hasse** invariant
``∏_{i<j}(a_i, a_j)_v`` (Lam's convention, [Lam05] V.3 — the name is
Macaulay2 heritage, and it is not the Witt/Cassels invariant); places of a
number field are always passed as prime ideals of the maximal order or as
infinite places, never as rational primes; and `==` is literal
Gram-matrix equality everywhere, with mathematical equivalence tested by
[`is_isomorphic_form`](@ref).

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["invariants.jl", "anisotropic_dimension.jl", "isotropy.jl", "isomorphism.jl"]
```

```@docs
local_anisotropic_dimension
```

## Étale algebras and transfer

Multiplication matrices, trace and norm of an étale algebra, and the
transfer of a Grothendieck–Witt class.

```@autodocs
Modules = [MotivicHomotopy]
Pages   = ["trace_norm.jl", "transfer.jl"]
```

## References

- [KR23] P. Koprowski and B. Rothkegel, *The anisotropic part of a quadratic form over a number field*, Journal of Symbolic Computation, 2023 (arXiv:2109.04172).
- [KC18] P. Koprowski and A. Czogała, *Computing with quadratic forms over number fields*, J. Symbolic Comput. 89 (2018), 129–145.
- [Kop20] P. Koprowski, *CQF Magma package*, ACM Commun. Comput. Algebra 54 (2020), 53–56.
- [DKR21] M. K. Darkey-Mensah, P. Koprowski and B. Rothkegel, *The anisotropic part of a quadratic form over a global function field*, ISSAC '21, 115–122.
- [LW92] D. B. Leep and A. R. Wadsworth, *The Hasse norm theorem mod squares*, J. Number Theory 42 (1992), 337–348.
- [Lam05] T. Y. Lam, *Introduction to quadratic forms over fields*, GSM 67, AMS, 2005.
- [Voi13] J. Voight, *Identifying the matrix ring*, in Quadratic and higher degree forms, Dev. Math. 31, Springer, 2013, 255–298.
