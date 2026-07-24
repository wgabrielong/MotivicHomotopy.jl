# Invariants of forms over a number field K: the fixed-diagonalisation value
# object _KRForm (memo §6), Hilbert symbols and Hasse invariants at places,
# the signed discriminant, and the S-normalised discriminant representative
# of memo §8.
#
# [KR21] P. Koprowski and B. Rothkegel, *The anisotropic part of a quadratic
# form over a number field*, arXiv:2109.04172.

# The primes of the maximal order above 2. Deliberately no caching layer:
# Hecke memoises the maximal order and the decomposition is cheap.
_dyadic_primes(K::AbsSimpleNumField) =
    [P for (P, _) in prime_decomposition(maximal_order(K), 2)]

# --- the fixed diagonalisation (memo §6) --------------------------------------------

# One fixed diagonalisation of a form over K, threaded through the whole KR
# pipeline. 𝔓(q) depends on the CHOSEN diagonal entries, not the isometry
# class, so diag, odd_primes = 𝔓(diag), and S = 𝔓 ∪ {dyadic primes} are
# computed once here and never re-derived from a Gram matrix; Algorithm 8's
# retry loop appends to S and never rebuilds it. The entries of diag are
# O_K-integral (each raw entry scaled by a denominator square — square-class
# neutral), so 𝔓 comes from integral-ideal factorisations.
struct _KRForm
    K::AbsSimpleNumField
    diag::Vector{AbsSimpleNumFieldElem}
    odd_primes::Vector{AbsNumFieldOrderIdeal{AbsSimpleNumField, AbsSimpleNumFieldElem}}
    S::Vector{AbsNumFieldOrderIdeal{AbsSimpleNumField, AbsSimpleNumFieldElem}}
    # the (embedding, isolating interval) pairs and the η's of [KR21]
    # Algorithm 2. They depend only on K and are deterministic (Hecke caches
    # the embeddings; the cuts are certified), so every class in a reduction
    # chain carries identical copies — the driver threads one η context.
    pairs::Vector
    etas::Vector{AbsSimpleNumFieldElem}
end

function _KRForm(beta::GWClass)
    M = gw_matrix(beta)
    (M isa MatElem && field_kind(base_ring(M)) isa NumberFieldKind) ||
        error("_KRForm requires a class over a number field")
    K = base_ring(M)
    OK = maximal_order(K)
    # deliberately NOT diagonal_class: its cache slot is overwritten by
    # sum_decomposition, so it is not stable within a session
    raw = _diag_of(diagonalize_via_congruence(M))
    diag = [a * K(denominator(a, OK))^2 for a in raw]
    odd = AbsNumFieldOrderIdeal{AbsSimpleNumField, AbsSimpleNumFieldElem}[]
    for a in diag
        for (P, e) in factor(ideal(OK, OK(a)))
            isodd(e) && !any(Q -> Q == P, odd) && push!(odd, P)
        end
    end
    S = copy(odd)
    for P in _dyadic_primes(K)
        any(Q -> Q == P, S) || push!(S, P)
    end
    pairs = _real_embeddings_with_intervals(K)
    _KRForm(K, diag, odd, S, pairs, _kr_etas(K, pairs))
end

# One _KRForm per class instance, held in the class's own cache so every
# invariant reads the same diagonalisation.
_kr_form(beta::GWClass) = get!(() -> _KRForm(beta), beta.cache, :kr_form)::_KRForm

# --- Hilbert symbols at places of K -------------------------------------------------

# The Hilbert symbol (a, b)_v of two nonzero elements of K at a place v —
# one implementation for every place kind. At a prime ideal it delegates to
# Hecke's hilbert_symbol (quadratic-defect algorithm, valid at ramified
# dyadic primes; convention pinned in the tests); at a real place the symbol
# is -1 exactly when σ(a) < 0 AND σ(b) < 0; at a complex place it is 1.
function _hilbert_symbol(a::AbsSimpleNumFieldElem, b::AbsSimpleNumFieldElem,
                         v::Union{AbsNumFieldOrderIdeal, InfPlc})
    (iszero(a) || iszero(b)) && error("Hilbert symbol arguments must be nonzero")
    v isa AbsNumFieldOrderIdeal && return hilbert_symbol(a, b, v)
    is_real(v) || return 1
    (sign(a, v) < 0 && sign(b, v) < 0) ? -1 : 1
end

# --- the exported K-facing symbol -----------------------------------------------

"""
    local_hilbert_symbol(a, b, v) -> Int

The Hilbert symbol ``(a, b)_v ∈ \\{1, -1\\}`` of two nonzero elements of a
number field ``K`` at a place `v` of ``K`` — the number-field counterpart of
[`hilbert_symbol_padic`](@ref) (which is ``\\mathbb{Q}_p``-only) and
[`hilbert_symbol_real`](@ref) (``\\mathbb{R}``-only). `v` is either a prime
ideal of the maximal order of ``K`` (including ramified dyadic primes) or an
infinite place of ``K`` (`InfPlc`); at a complex place the symbol is
identically 1. `a` and `b` must be elements of ``K`` or rationals/integers,
which are coerced.

``(a, b)_v = 1`` iff ``z^2 = a x^2 + b y^2`` has a nonzero solution over the
completion ``K_v``. The symbol is bimultiplicative and satisfies Hilbert
reciprocity: ``∏_v (a, b)_v = 1`` over all places of ``K``.

# Examples
```julia-repl
julia> K, s = quadratic_field(-7);

julia> P = prime_decomposition(maximal_order(K), 2)[1][1];  # 2 splits: K_P = Q_2

julia> local_hilbert_symbol(-1, -1, P)
-1

julia> K2, t = quadratic_field(2);

julia> local_hilbert_symbol(-1, t, real_places(K2)[1]) *
       local_hilbert_symbol(-1, t, real_places(K2)[2])
-1
```

# References
- [Lam05] T. Y. Lam, *Introduction to quadratic forms over fields*, GSM 67, AMS, 2005.
- [Voi13] J. Voight, *Identifying the matrix ring*, in Quadratic and higher degree forms, Dev. Math. 31, Springer, 2013, 255–298 — the quadratic-defect algorithm behind the symbol at (ramified) dyadic primes.

See also [`hilbert_symbol_padic`](@ref), [`hasse_witt_invariant`](@ref),
[`local_anisotropic_dimension`](@ref).
"""
function local_hilbert_symbol(a, b, v::Union{AbsNumFieldOrderIdeal, InfPlc})
    K = v isa AbsNumFieldOrderIdeal ? number_field(order(v)) : number_field(v)
    _check_nonsimple(K)
    aa = _coerce_into(K, a)
    bb = _coerce_into(K, b)
    (aa === nothing || bb === nothing) &&
        error("the first two arguments must lie in the base field of the place")
    v isa AbsNumFieldOrderIdeal && !is_prime(v) &&
        error("the place must be a prime ideal or an infinite place")
    _hilbert_symbol(aa, bb, v)
end

# Hasse invariant s_v = prod_{i<j} (a_i, a_j)_v of an explicit diagonal
# vector; the KR driver applies it to extended diagonals that exist only as
# lists.
function _hasse_witt_from_diagonal(D::Vector, v)
    s = 1
    for i in 1:(length(D) - 1), j in (i + 1):length(D)
        s *= _hilbert_symbol(D[i], D[j], v)
    end
    s
end

# --- Hasse-Witt invariant at a place ------------------------------------------------
# Delegates to Hecke's hasse_invariant, which is Lam's prod_{i<j} convention
# (verified in source, Spaces.jl:215, and pinned by hand-recomputation in the
# tests — memo §5). Never witt_invariant: that is the Cassels/Kirschmer
# invariant.

function _check_nf_class_and_field(beta::GWClass, L)
    M = gw_matrix(beta)
    M isa MatElem && _check_nonsimple(base_ring(M))
    (M isa MatElem && field_kind(base_ring(M)) isa NumberFieldKind) ||
        error("method requires a class over a number field")
    L === base_ring(M) || error("the place does not belong to the base field of the class")
    nothing
end

function hasse_witt_invariant(beta::GWClass, P::AbsNumFieldOrderIdeal)
    _check_nf_class_and_field(beta, number_field(order(P)))
    is_prime(P) || error("second argument must be a prime ideal")
    hasse_invariant(quadratic_space(gw_base_field(beta), gw_matrix(beta)), P)
end

function hasse_witt_invariant(beta::GWClass, sigma::InfPlc)
    _check_nf_class_and_field(beta, number_field(sigma))
    hasse_invariant(quadratic_space(gw_base_field(beta), gw_matrix(beta)), sigma)
end

# --- relevant primes and signed discriminant over K ---------------------------------

function _relevant_primes_nf(beta::GWClass)
    copy(_kr_form(beta).odd_primes)
end

function _signed_discriminant_nf(beta::GWClass)
    kr = _kr_form(beta)
    n = length(kr.diag)
    sgn = iseven(div(n * (n - 1), 2)) ? 1 : -1
    kr.K(sgn) * prod(kr.diag; init = one(kr.K))
end

# --- isometry over K (Hasse-Minkowski) ----------------------------------------------
# Equal rank, equal signature at every real place, equal discriminant square
# class, equal Hasse invariant at every finite place. Off 𝔓 ∪ {dyadics} both
# Hasse invariants are +1 (the fixed diagonals have even valuation there and
# unit symbols at odd primes are trivial), so the union of the two S sets
# covers every place that can distinguish the forms.
function _is_isomorphic_form_nf(alpha::GWClass, beta::GWClass)
    form_rank(alpha) == form_rank(beta) || return false
    form_signatures(alpha) == form_signatures(beta) || return false
    is_square(signed_discriminant(alpha) * signed_discriminant(beta)) || return false
    places = copy(_kr_form(alpha).S)
    for P in _kr_form(beta).S
        any(Q -> Q == P, places) || push!(places, P)
    end
    for P in places
        hasse_witt_invariant(alpha, P) == hasse_witt_invariant(beta, P) || return false
    end
    true
end

# --- S-normalised discriminant representative (memo §8, corrected invariant) --------

# An element with valuation 1 at P and valuation 0 at every other prime of
# S. uniformizer(P) alone does NOT control the valuation at the other
# members of S, so the CRT correction is required, not optional:
# π ≡ uniformizer(P) (mod P²) forces valuation exactly 1 at P, and π ≡ 1
# (mod Q) forces valuation 0 at each other Q ∈ S.
function _local_uniformizer(P::AbsNumFieldOrderIdeal, S::Vector)
    OK = order(P)
    K = number_field(OK)
    others = [Q for Q in S if Q != P]
    isempty(others) && return K(uniformizer(P))
    moduli = vcat([P^2], [Q^1 for Q in others])
    residues = vcat([uniformizer(P)], [one(OK) for _ in others])
    K(crt(residues, moduli))
end

# A representative d' of the square class of d, normalised for the
# congruences of [KR21] Algorithms 6/7. The guaranteed (and asserted)
# invariant: d' ∈ O_K, and ord_P(d') ≡ ord_P(d) (mod 2) for every P ∈ S.
# The valuations at S are NOT forced into {0, 1} — the final
# denominator-clearing square need not be a unit on S. Harmless: the
# algorithms read only the parity, and their odd branch needs only
# ord_P(d') > 0, which odd parity plus integrality gives.
function _disc_rep_normalized(d::AbsSimpleNumFieldElem, S::Vector)
    iszero(d) && error("discriminant representative of 0 requested")
    K = parent(d)
    OK = maximal_order(K)
    dp = d
    for P in S
        e = fld(valuation(d, P), 2)          # fld, not div: v may be negative
        iszero(e) && continue
        dp *= _local_uniformizer(P, S)^(-2 * e)
    end
    dp *= K(denominator(dp, OK))^2

    isone(denominator(dp, OK)) ||
        error("_disc_rep_normalized postcondition failed: not integral")
    for P in S
        iseven(valuation(dp, P) - valuation(d, P)) ||
            error("_disc_rep_normalized postcondition failed: parity at $P")
    end
    dp
end
