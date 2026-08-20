# Port of M2 Code/AnisotropicDimension.m2: isHyperbolicQQp,
# getAnisotropicDimensionQQp, getAnisotropicDimensionQQ, getAnisotropicDimension,
# getWittIndex. Koprowski-Czogala Algorithms 6, 8, 9.
#
# Note the SIGNED discriminants: Koprowski-Czogala use (-1)^(r(r±1)/2) * disc —
# port the exponents exactly.

const _UNSUPPORTED_FIELD_ERR = "Base field not supported; only implemented over " *
    "QQ, RR, CC, and finite fields of characteristic not 2"

# M2: isHyperbolicQQp — is beta totally hyperbolic over QQ_p? (Algorithm 6)
function _is_hyperbolic_qqp(beta::GWClass, p)
    field_kind(gw_base_field(beta)) isa RationalKind ||
        error("GrothendieckWittClass is not over QQ")
    is_prime(ZZ(p)) || error("second argument must be a prime number")
    r = form_rank(beta)
    isodd(r) && return false
    d = signed_discriminant(beta)
    _is_padic_square(d, p) || return false
    m = div(r, 2)
    hw_hyperbolic = hilbert_symbol_padic(-1, -1, p)^div(m * (m - 1), 2)
    hw_hyperbolic == hasse_witt_invariant(beta, p)
end

# M2: getAnisotropicDimensionQQp — always 0..4 (Algorithm 8)
"""
    anisotropic_dimension_qqp(beta::GWClass, p)

The anisotropic dimension of a rational form over the ``p``-adic completion
``\\mathbb{Q}_p``: the rank of the anisotropic part of `beta` base-changed to ``\\mathbb{Q}_p``.
Every form of rank ≥ 5 over ``\\mathbb{Q}_p`` is isotropic, so the result is always
0, 1, 2, 3, or 4. This implements [KC18, Algorithm 8].

# Examples
```julia-repl
julia> anisotropic_dimension_qqp(diagonal_form(QQ, (1, -1, 2)), 2)
1
```

# References
- [KC18] P. Koprowski and A. Czogała, *Computing with quadratic forms over number fields*, Journal of Symbolic Computation, 2018.

See also [`anisotropic_dimension`](@ref).
"""
function anisotropic_dimension_qqp(beta::GWClass, p)
    field_kind(gw_base_field(beta)) isa RationalKind ||
        error("GrothendieckWittClass is not over QQ")
    is_prime(ZZ(p)) || error("second argument must be a prime number")
    r = form_rank(beta)
    if iseven(r)
        _is_hyperbolic_qqp(beta, p) && return 0
        d = signed_discriminant(beta)
        return _is_padic_square(d, p) ? 4 : 2
    else
        c = (-1)^div(r * (r + 1), 2) * integral_discriminant(beta)
        gamma = gw_direct_sum(beta, diagonal_form(QQ, (c,)))
        return _is_hyperbolic_qqp(gamma, p) ? 1 : 3
    end
end

# M2: getAnisotropicDimensionQQ — max over RR (|signature|), QQ_2, and the
# relevant primes (Algorithm 9)
function _anisotropic_dimension_qq(beta::GWClass)
    field_kind(gw_base_field(beta)) isa RationalKind ||
        error("GrothendieckWittClass is not over QQ")
    dims = Int[abs(form_signature(beta)), anisotropic_dimension_qqp(beta, 2)]
    for p in relevant_primes(beta)
        push!(dims, anisotropic_dimension_qqp(beta, p))
    end
    maximum(dims)
end

# M2: getAnisotropicDimension, per field
"""
    anisotropic_dimension(beta::GWClass; certify = false)
    anisotropic_dimension(A; certify = false)

The anisotropic dimension of a symmetric bilinear form over ``\\mathbb{Q}``, a number
field, ``\\mathbb{R}``, ``\\mathbb{C}``, or a
finite field of odd characteristic. With `certify = true` a float Gram matrix
is read exactly and the underlying sign counts (over ``\\mathbb{R}``) or rank
(over ``\\mathbb{C}``) are computed exactly over ``\\mathbb{Q}`` or
``\\mathbb{Q}(i)``; on exact base fields the keyword is a no-op. By the Witt decomposition theorem any
nondegenerate form decomposes uniquely as ``β ≅ n·\\mathbb{H} ⊕ β_a`` with ``β_a``
anisotropic; the anisotropic dimension is the rank of ``β_a``.

Over ``\\mathbb{Q}`` it is the maximum of the anisotropic dimensions over all completions:
``\\lvert \\text{signature} \\rvert`` at the real place and
[`anisotropic_dimension_qqp`](@ref) at 2 and the [`relevant_primes`](@ref)
([KC18, Algorithm 9]). Over a number field the same maximum runs over the
real places and the finite places of ``𝔓(q) ∪ \\{\\text{dyadic}\\}``
([`local_anisotropic_dimension`](@ref)), **plus one global correction** the
bare maximum misses: a form of even rank whose signed discriminant is not a
square in ``K`` has anisotropic dimension at least 2, witnessed at a place
outside that finite set (the concrete witness is
``⟨1, -(35+6\\sqrt{34})⟩`` over ``\\mathbb{Q}(\\sqrt{34})``, anisotropic with
every local dimension on the finite set zero).

# Examples
```julia-repl
julia> anisotropic_dimension(diagonal_form(QQ, (1, -1, 2)))
1
```

# References
- [KC18] P. Koprowski and A. Czogała, *Computing with quadratic forms over number fields*, Journal of Symbolic Computation, 2018.

See also [`gw_witt_index`](@ref), [`anisotropic_part`](@ref),
[`is_anisotropic_form`](@ref).
"""
function anisotropic_dimension(A::MatElem; certify::Bool = false)
    R = base_ring(A)
    _require_form_supported(R, _UNSUPPORTED_FIELD_ERR)
    _is_symmetric_square(A) || error("Matrix is not symmetric")
    _anisotropic_dimension(field_kind(R), A)
end

_anisotropic_dimension(::RationalKind, A::MatElem) =
    _anisotropic_dimension_qq(GWClass(_nondegenerate_part_diagonal(A)))

function _anisotropic_dimension(::FiniteKind, A::MatElem)
    R = base_ring(A)
    diagA = diagonalize_via_congruence(A)
    r = form_rank(diagA)   # number of nonzero diagonal entries
    isodd(r) && return 1
    if is_square(det(_nondegenerate_part_diagonal(diagA))) ==
       is_square(R((-1)^div(r, 2)))
        return 0
    else
        return 2
    end
end

_anisotropic_dimension(::FieldKind, ::MatElem) = error(_UNSUPPORTED_FIELD_ERR)

function anisotropic_dimension(A::Matrix{Float64}; certify::Bool = false)
    _is_symmetric_square(A) || error("Matrix is not symmetric")
    # certify=true: exact counts over QQ from the input itself, beyond M2; the
    # default keeps M2's float pre-diagonalization, so certification never
    # inherits a float pivot decision.
    pos, neg = certify ? _certified_pos_neg_counts(A) :
        _pos_neg_diag_counts(diagonalize_via_congruence(A))
    abs(pos - neg)
end

function anisotropic_dimension(A::Matrix{ComplexF64}; certify::Bool = false)
    _is_symmetric_square(A) || error("Matrix is not symmetric")
    form_rank(A; certify = certify) % 2
end

anisotropic_dimension(alpha::GWClass; certify::Bool = false) =
    anisotropic_dimension(gw_matrix(alpha); certify = certify)

# M2: getWittIndex — (rank - anisotropic dimension)/2
"""
    gw_witt_index(beta::GWClass; certify = false)

The Witt index of a form over ``\\mathbb{Q}``, ``\\mathbb{R}``, ``\\mathbb{C}``, or a finite field of odd
characteristic: the number ``n`` of hyperbolic summands in the Witt
decomposition ``β ≅ n·\\mathbb{H} ⊕ β_a`` ([L05, I.4.3]), computed as
(rank − anisotropic dimension)/2. `certify` is forwarded to
[`anisotropic_dimension`](@ref).

The name carries the `gw_` prefix because Oscar exports `witt_index`.

# Examples
```julia-repl
julia> gw_witt_index(diagonal_form(QQ, (1, -1, 2)))
1
```

# References
- [L05] T. Y. Lam, *Introduction to quadratic forms over fields*, American Mathematical Society, 2005.

See also [`anisotropic_dimension`](@ref), [`sum_decomposition`](@ref).
"""
gw_witt_index(alpha::GWClass; certify::Bool = false) =
    div(form_rank(alpha) - anisotropic_dimension(alpha; certify = certify), 2)
