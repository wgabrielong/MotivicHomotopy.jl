"""
    MotivicHomotopy

Computing local and global ``\\mathbb{A}^1``-Brouwer degrees and studying the symmetric
bilinear forms they produce, over ``\\mathbb{Q}``, finite fields of odd characteristic, ``\\mathbb{R}``,
and ``\\mathbb{C}``.

The ``\\mathbb{A}^1``-Brouwer degree of an endomorphism of affine space with isolated
zeros is a class in the Grothendieck–Witt ring ``\\text{GW}(k)`` of symmetric bilinear
forms; for pointed rational functions ``\\mathbb{P}^1_k → \\mathbb{P}^1_k`` the degree lives in
the unstable group ``\\text{GW}^u(k)``. The package provides:

- **Classes**: [`GWClass`](@ref), [`GWuClass`](@ref) with accessors
  ([`gw_matrix`](@ref), [`gw_scalar`](@ref), [`gw_algebra`](@ref),
  [`gw_base_field`](@ref), [`stable_part`](@ref)) and operations
  ([`gw_direct_sum`](@ref), [`gw_tensor_product`](@ref), the integer
  multiple `n * beta` and tensor power `beta^n`,
  [`gw_base_change`](@ref), [`divisorial_sum`](@ref)).

  GW classes form a **semiring** here, as in the Macaulay2 package: sums
  (orthogonal sum) and products (tensor product) exist, formal differences
  do not, so there is deliberately no negation, subtraction, or unary
  minus. ``⟨-1⟩ ⊗ β`` is the *opposite form* (all signs flipped), not an
  additive inverse — ``β ⊕ (⟨-1⟩ ⊗ β)`` is hyperbolic, not zero. Genuine
  additive inverses live in the Grothendieck group of formal differences,
  which this package does not model.
- **Form constructors**: [`diagonal_form`](@ref), [`hyperbolic_form`](@ref),
  [`pfister_form`](@ref) and unstable variants.
- **Simplification**: [`diagonalize_via_congruence`](@ref),
  [`diagonal_class`](@ref), [`diagonal_entries`](@ref),
  [`sum_decomposition`](@ref), [`sum_decomposition_string`](@ref),
  [`anisotropic_part`](@ref).
- **Invariants**: [`form_rank`](@ref), [`form_signature`](@ref),
  [`form_signatures`](@ref),
  [`integral_discriminant`](@ref), [`signed_discriminant`](@ref),
  [`hilbert_symbol_padic`](@ref),
  [`hilbert_symbol_real`](@ref), [`local_hilbert_symbol`](@ref),
  [`hasse_witt_invariant`](@ref),
  [`anisotropic_dimension`](@ref), [`anisotropic_dimension_qqp`](@ref),
  [`local_anisotropic_dimension`](@ref), [`gw_witt_index`](@ref),
  [`is_isotropic_form`](@ref), [`is_anisotropic_form`](@ref),
  [`is_isomorphic_form`](@ref).
- **Degrees**: [`global_A1_degree`](@ref), [`local_A1_degree`](@ref),
  [`global_unstable_A1_degree`](@ref), [`local_unstable_A1_degree`](@ref).
- **Étale algebras**: [`multiplication_matrix`](@ref),
  [`algebra_trace`](@ref), [`algebra_norm`](@ref), [`transfer_gw`](@ref).

Which arithmetic runs is determined by the base field. Over ``\\mathbb{Q}`` and finite
fields of odd characteristic everything is computed exactly. **Number
fields** (any `AbsSimpleNumField`) are supported for the form layer:
construction, the ring operations, invariants at places
([`local_hilbert_symbol`](@ref), [`local_anisotropic_dimension`](@ref),
[`hasse_witt_invariant`](@ref) at prime ideals, [`form_signatures`](@ref)),
isometry testing, Witt decomposition, and the anisotropic part via the
Koprowski–Rothkegel algorithms — see the manual page *Working with forms*
for conventions and the GRH caveat. Over ``\\mathbb{C}`` a form is
determined by its rank, so the stable degrees are the identity form of the
algebra dimension (a discrete invariant from the exact computation); only the
*unstable* degrees, which carry a ``k^×``-scalar, use numerical
(HomotopyContinuation) computation over ``\\mathbb{C}``. Over ``\\mathbb{R}`` the degree functions direct
the user to compute over ``\\mathbb{Q}`` and base-change the result.
Functions that decide from float input accept an opt-in `certify = true`
keyword: the floats are read as the exact (Gaussian) rationals they are and
every discrete decision is computed exactly over ``\\mathbb{Q}`` or
``\\mathbb{Q}(i)``, erroring rather than guessing when a faithful answer
cannot be returned — see the manual's *Certified computation* section.

This package is a Julia port of the Macaulay2 package `A1BrouwerDegrees`
(v2.0); v1.1 is published in the Journal of Software for Algebra and Geometry
(14, 2024). It was ported by Claude Fable 5 (Anthropic), and its behavior was
verified against the Macaulay2 package.
"""
module MotivicHomotopy

using Oscar
import LinearAlgebra
import HomotopyContinuation as HC

export GWClass, GWuClass,
    gw_matrix, gw_scalar, gw_algebra, gw_base_field, stable_part,
    gw_direct_sum, gw_tensor_product, gw_base_change, divisorial_sum,
    multiplication_matrix, algebra_trace, algebra_norm,
    diagonalize_via_congruence, diagonal_class, diagonal_entries,
    padic_valuation,
    diagonal_form, hyperbolic_form, pfister_form,
    diagonal_unstable_form, hyperbolic_unstable_form,
    form_rank, form_signature, form_signatures, integral_discriminant,
    signed_discriminant, relevant_primes,
    hilbert_symbol_padic, hilbert_symbol_real, local_hilbert_symbol,
    hasse_witt_invariant,
    anisotropic_dimension, anisotropic_dimension_qqp,
    local_anisotropic_dimension, gw_witt_index,
    is_anisotropic_form, is_isotropic_form, is_isomorphic_form,
    anisotropic_part, sum_decomposition, sum_decomposition_string,
    global_A1_degree, local_A1_degree, local_algebra_basis,
    transfer_gw, global_unstable_A1_degree, local_unstable_A1_degree

# Loaded first: the GW well-definedness check needs the trace form of an étale
# algebra.
include("trace_norm.jl")
include("arithmetic.jl")
include("field_kind.jl")
include("nf_arithmetic.jl")
include("nf_singular_elements.jl")
include("nf_signs.jl")
include("gw_classes.jl")
include("gwu_classes.jl")
include("matrix_methods.jl")
# After gw_classes.jl (FloatGram) and matrix_methods.jl (the exact congruence
# elimination it reuses): the certified (certify = true) float decisions.
include("certified.jl")
include("simplified_representatives.jl")
include("building_forms.jl")
include("hilbert_symbols.jl")
include("invariants.jl")
include("nf_invariants.jl")
include("anisotropic_dimension.jl")
include("nf_anisotropic_dimension.jl")
include("nf_anisotropic_part.jl")
include("isotropy.jl")
include("isomorphism.jl")
include("decomposition.jl")
include("degrees.jl")
include("transfer.jl")
include("unstable_degrees.jl")

# The ℂ numerical backend (HomotopyContinuation.jl) adds methods to the degree
# functions for floating-point input; HC names are qualified `HC.` so they do
# not collide with Oscar's exports.
include("numerical.jl")

end # module
