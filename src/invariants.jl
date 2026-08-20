# Port of M2 Code/GWInvariants.m2: getRank, countPos/NegDiagEntries (internal),
# getSignature, getIntegralDiscriminant, getRelevantPrimes, getHasseWittInvariant.

# --- M2: getRank ------------------------------------------------------------------
# On a class: the matrix size. On a raw matrix: the actual rank (drops degeneracy).

"""
    form_rank(beta::GWClass; certify = false)
    form_rank(M; certify = false)

The rank of a symmetric bilinear form. On a [`GWClass`](@ref) (which is
nondegenerate by construction) this is the size of the Gram matrix; on a raw
matrix it is the matrix rank, so degenerate directions are not counted.

With `certify = true` a float Gram matrix is read exactly (each entry as the
dyadic rational it is) and the rank is computed exactly over ``\\mathbb{Q}``
or ``\\mathbb{Q}(i)`` instead of by the default floating-point SVD; on exact
base fields the keyword is accepted and is a no-op. See the manual for the
certification semantics.

# Examples
```julia-repl
julia> form_rank(diagonal_form(QQ, (3, 5, 7, 11)))
4
```

See also [`form_signature`](@ref), [`anisotropic_dimension`](@ref).
"""
form_rank(beta::GWClass; certify::Bool = false) = _n_of(beta.matrix)
form_rank(M::MatElem; certify::Bool = false) = nrows(M) == 0 ? 0 : rank(M)
# certify=true: exact rank over QQ / QQ(i), beyond M2; default branch unchanged
form_rank(M::FloatGram; certify::Bool = false) =
    certify ? _certified_rank(M) :
    (size(M, 1) == 0 ? 0 : LinearAlgebra.rank(M))

_n_of(M::MatElem) = nrows(M)
_n_of(M::FloatGram) = size(M, 1)

# --- M2: countPosDiagEntries / countNegDiagEntries (internal; QQ and RR only) ------

function _pos_neg_diag_counts(M::Union{MatElem, Matrix{Float64}};
                              certify::Bool = false)
    M isa MatElem && !(field_kind(base_ring(M)) isa RationalKind) &&
        error("Only implemented over QQ and RR")
    # certify=true: exact counts over QQ, beyond M2; a MatElem input is already
    # exact, so only the float branch changes. The certified helper lifts and
    # rediagonalizes M itself, so no float diagonalization enters the counts.
    certify && M isa Matrix{Float64} && return _certified_pos_neg_counts(M)
    _is_symmetric_square(M) || error("Matrix is not symmetric")
    A = _is_diag(M) ? M : diagonalize_via_congruence(M)
    pos = 0
    neg = 0
    for i in 1:_n_of(A)
        A[i, i] > 0 && (pos += 1)
        A[i, i] < 0 && (neg += 1)
    end
    (pos, neg)
end
_pos_neg_diag_counts(::Any; certify::Bool = false) =
    error("Only implemented over QQ and RR")

# --- M2: getSignature (QQ and RR only) ----------------------------------------------
# Over QQ, reuses Hecke's signature_tuple on the quadratic space (verified identical
# to M2's diagonal-count signature on a 200-form sweep). The RR (float) path keeps
# M2's diagonal counts.

"""
    form_signature(beta::GWClass; certify = false)

The signature of a symmetric bilinear form over ``\\mathbb{Q}`` or ``\\mathbb{R}``: after diagonalizing,
the number of positive diagonal entries minus the number of negative ones.
Together with the rank it classifies forms over ``\\mathbb{R}``; over ``\\mathbb{Q}`` it is one of the
invariants entering [`is_isomorphic_form`](@ref).

With `certify = true` a float Gram matrix is read exactly and the diagonal
sign counts are computed over ``\\mathbb{Q}`` instead of in floating point;
over ``\\mathbb{Q}`` itself the keyword is accepted and is a no-op (the
computation is already exact).

# Examples
```julia-repl
julia> form_signature(GWClass([0.0 0.0 1.0; 0.0 1.0 0.0; 1.0 0.0 0.0]))
1

julia> form_signature(diagonal_form(QQ, (1, -1, 1)))
1
```

See also [`form_rank`](@ref), [`integral_discriminant`](@ref),
[`is_isomorphic_form`](@ref).
"""
function form_signature(beta::GWClass; certify::Bool = false)
    M = gw_matrix(beta)
    if M isa MatElem && field_kind(base_ring(M)) isa RationalKind
        st = signature_tuple(quadratic_space(QQ, M))   # (pos, zero, neg)
        return st[1] - st[3]
    end
    pos, neg = _pos_neg_diag_counts(M; certify = certify)
    pos - neg
end

# --- signatures over a number field ------------------------------------------------
# Over K a form has one signature per real place, so the single-signature
# method above stays QQ/RR-only and the number-field API is per-place.

"""
    form_signature(beta::GWClass, sigma)

The signature of a form over a number field ``K`` at one real place (or real
embedding) `sigma` of ``K``: the number of diagonal entries of a
diagonalization that are positive under `sigma` minus the number that are
negative. `sigma` may be an `InfPlc` or a `NumFieldEmb`; it must be real and
must belong to the base field of `beta`.

See also [`form_signatures`](@ref).
"""
function form_signature(beta::GWClass, sigma::Union{InfPlc, Hecke.NumFieldEmb})
    M = gw_matrix(beta)
    M isa MatElem && _check_nonsimple(base_ring(M))
    (M isa MatElem && field_kind(base_ring(M)) isa NumberFieldKind) ||
        error("form_signature at a place requires a class over a number field")
    number_field(sigma) === base_ring(M) ||
        error("the place does not belong to the base field of the class")
    is_real(sigma) || error("the place is not real")
    sum(Int(sign(d, sigma)) for d in diagonal_entries(beta); init = 0)
end

"""
    form_signatures(beta::GWClass)

The signatures of a form over a number field ``K`` at all real places of
``K``, as a `Vector{Int}` indexed like `real_places(K)`. Empty when ``K``
has no real embeddings.

See also [`form_signature`](@ref).
"""
function form_signatures(beta::GWClass)
    M = gw_matrix(beta)
    M isa MatElem && _check_nonsimple(base_ring(M))
    (M isa MatElem && field_kind(base_ring(M)) isa NumberFieldKind) ||
        error("form_signatures requires a class over a number field")
    D = diagonal_entries(beta)
    [sum(Int(sign(d, sigma)) for d in D; init = 0)
     for sigma in real_places(base_ring(M))]
end

# --- M2: getIntegralDiscriminant (QQ only) ------------------------------------------

"""
    integral_discriminant(beta::GWClass)

A squarefree integral representative of the discriminant of a form over ``\\mathbb{Q}``:
the square class of the determinant of any Gram matrix representing `beta`,
normalized to a squarefree integer. The discriminant is one of the invariants
classifying rational forms (see [`is_isomorphic_form`](@ref)).

# Examples
```julia-repl
julia> beta = GWClass(QQ[1 4 7; 4 3 -1; 7 -1 5]);

julia> integral_discriminant(beta)
-269
```

See also [`form_signature`](@ref), [`hasse_witt_invariant`](@ref).
"""
function integral_discriminant(beta::GWClass)
    field_kind(gw_base_field(beta)) isa RationalKind ||
        error("GrothendieckWittClass is not over QQ")
    _squarefree_part(prod(diagonal_entries(beta); init = one(QQ)))
end

# --- signed discriminant ([KR21] §2 convention) -------------------------------------

"""
    signed_discriminant(beta::GWClass)

The signed discriminant ``disc(q) = (-1)^{n(n-1)/2} a_1 ⋯ a_n`` of a diagonal
representative ``⟨a_1, …, a_n⟩`` — the convention of [KR21] §2 (Szymiczek,
Def. 15.2.1). This is the discriminant the number-field algorithms are stated
in; [`integral_discriminant`](@ref) is the *unsigned* determinant square
class and stays as the M2-faithful ℚ API.

Over ℚ the result is the squarefree integer representative, i.e. literally
``(-1)^{n(n-1)/2}`` times [`integral_discriminant`](@ref). Over a finite
field it is a field element. Over a number field it is an ``O_K``-integral
representative, well-defined only up to nonzero squares (there is no
canonical squarefree representative over ``K`` — see
`_squarefree_content_rep`).

# References
- [KR21] P. Koprowski and B. Rothkegel, *The anisotropic part of a quadratic form over a number field*, arXiv:2109.04172.

See also [`integral_discriminant`](@ref), [`hasse_witt_invariant`](@ref).
"""
function signed_discriminant(beta::GWClass)
    M = gw_matrix(beta)
    kind = M isa MatElem ? field_kind(base_ring(M)) : UnsupportedKind()
    n = form_rank(beta)
    sgn = iseven(div(n * (n - 1), 2)) ? 1 : -1
    if kind isa RationalKind
        return sgn * integral_discriminant(beta)
    elseif kind isa FiniteKind
        R = base_ring(M)
        return R(sgn) * prod(diagonal_entries(beta); init = one(R))
    elseif kind isa NumberFieldKind
        return _signed_discriminant_nf(beta)
    end
    M isa MatElem && _check_nonsimple(base_ring(M))
    error("signed_discriminant is only implemented over QQ, finite fields, " *
          "and number fields")
end

# --- M2: getRelevantPrimes (QQ only) ------------------------------------------------
# Primes dividing the entries of the squarefree diagonal representative, in
# first-seen order with each entry's factors ascending (M2's order — Test 24
# expects [23, 5, 47]).

"""
    relevant_primes(beta::GWClass)

A finite list of primes containing every prime at which the Hasse–Witt
invariant of the rational form `beta` can be nontrivial. The Hasse–Witt
invariants of a form equal 1 at all but finitely many primes ([S73, IV
§3.3]); since they are products of Hilbert symbols of the diagonal entries, it
suffices to take the primes dividing the entries of a squarefree diagonal
representative.

# Examples
```julia-repl
julia> relevant_primes(diagonal_form(QQ, (6, 7, 22)))
4-element Vector{ZZRingElem}:
 2
 3
 7
 11
```

# References
- [S73] J. P. Serre, *A course in arithmetic*, Springer-Verlag, 1973.

See also [`hasse_witt_invariant`](@ref).
"""
function relevant_primes(beta::GWClass)
    kind = field_kind(gw_base_field(beta))
    # over a number field: 𝔓(q) of the FIXED diagonalisation held by the
    # class's _KRForm (memo §6), as prime ideals — never recomputed per call
    kind isa NumberFieldKind && return _relevant_primes_nf(beta)
    kind isa RationalKind || error("GrothendieckWittClass is not over QQ")
    D = diagonal_entries(diagonal_class(beta))
    L = ZZRingElem[]
    for x in D
        for p in _prime_factors(ZZ(x))
            p in L || push!(L, p)
        end
    end
    L
end

# --- M2: getHasseWittInvariant ------------------------------------------------------
# Product of pairwise Hilbert symbols over i < j (M2's convention; some sources
# use i <= j — do not change this).

"""
    hasse_witt_invariant(beta::GWClass, p)
    hasse_witt_invariant(L::AbstractVector, p)

The Hasse–Witt invariant at the prime `p` of a form over ``\\mathbb{Q}``: for a
diagonalization ``⟨a_1, …, a_n⟩``, the product ``∏_{i<j} (a_i, a_j)_p`` of
pairwise Hilbert symbols (see [`hilbert_symbol_padic`](@ref)). The list
variant takes the diagonal entries directly.

The invariant equals 1 for all but finitely many primes — for `p` not
dividing any entry of a squarefree diagonal representative it is
automatically 1 — so only the [`relevant_primes`](@ref) need checking.

# Examples
```julia-repl
julia> beta = GWClass(QQ[1 4 7; 4 3 -1; 7 -1 5]);

julia> hasse_witt_invariant(beta, 7)
1

julia> hasse_witt_invariant([6, 7, 22], 2)
-1
```

See also [`hilbert_symbol_padic`](@ref), [`relevant_primes`](@ref),
[`is_isomorphic_form`](@ref).
"""
function hasse_witt_invariant(L::AbstractVector, p)
    is_prime(ZZ(p)) || error("second argument must be a prime number")
    f = [_squarefree_part(x) for x in L]
    a = 1
    for i in 1:(length(f) - 1), j in (i + 1):length(f)
        a *= hilbert_symbol_padic(f[i], f[j], p)
    end
    a
end

# The class-level variant reuses Hecke's hasse_invariant (∏_{i<j} convention,
# verified identical to M2's on a 200-form sweep at 2/3/all relevant primes).
# The list variant above stays as M2's own product so that arbitrary
# user-supplied diagonal lists behave exactly as in M2.
function hasse_witt_invariant(beta::GWClass, p)
    field_kind(gw_base_field(beta)) isa RationalKind ||
        error("method is only implemented over the rational numbers")
    is_prime(ZZ(p)) || error("second argument must be a prime number")
    hasse_invariant(quadratic_space(QQ, gw_matrix(beta)), ZZ(p))
end
