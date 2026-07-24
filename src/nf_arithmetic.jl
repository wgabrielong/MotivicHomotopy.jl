# Number-field arithmetic helpers: real embeddings coupled with isolating
# intervals, the η elements of [KR21] Algorithm 2, and square-class
# representative normalisation.
#
# There is deliberately no per-field context object or module-level cache here:
# Hecke already
# memoises `maximal_order` on the field (MaxOrd.jl:41) and `class_group` on
# the order (ClassGrpCtx attribute, Clgp.jl:63), and `real_embeddings` is
# deterministic because `complex_embeddings` is cached on the field
# (NfAbs.jl:27). Callers use those directly.

# The real embeddings of K, each paired with a rational open interval (a, b)
# isolating its image σ(θ) of the primitive element: a < σ(θ) < b, no other
# embedding's image in [a, b]. The interval is derived FROM the embedding
# (its arb enclosure, refined until all enclosures are pairwise disjoint) and
# returned in the same tuple, so alignment holds by construction. The
# coupling is load-bearing: the η's of [KR21] Algorithm 2 are indexed by
# embedding, and a misalignment would make matrix A of Algorithm 8 sign-wrong
# — silently changing the output form. Do not split the pairs into parallel
# containers.
function _real_embeddings_with_intervals(K::AbsSimpleNumField)
    embs = real_embeddings(K)
    pairs = Tuple{eltype(embs), Tuple{QQFieldElem, QQFieldElem}}[]
    isempty(embs) && return pairs
    theta = gen(K)
    prec = 64
    while true
        balls = [real(sigma(theta, prec)) for sigma in embs]
        intervals = _isolating_intervals(balls)
        if intervals !== nothing
            return [(embs[i], intervals[i]) for i in eachindex(embs)]
        end
        prec *= 2
        prec > 2^16 &&
            error("failed to isolate the real embeddings of $K at precision $(prec ÷ 2)")
    end
end

# Rational endpoints separating pairwise-disjoint arb balls, in the input
# order; `nothing` when the balls overlap or a candidate endpoint cannot be
# certified at this precision (caller refines and retries). Every candidate
# is certified by arb comparisons, which return true only when provable, so
# a successful return is exact.
#
# The endpoints are chosen with SMALL numerator and denominator (coarsest
# certifiable continued-fraction approximation), not as raw binary
# midpoints: the η elements of [KR21] Algorithm 2 are built from these
# endpoints, and 2^prec-denominator endpoints propagate — through the t²
# integralisation and the corrected Algorithm 4 modulus
# k_p = 1 + ord_p(4) + ord_p(α₁) — into congruence moduli and CRT solutions
# so large that factoring the extended form's entries becomes infeasible.
function _isolating_intervals(balls)
    n = length(balls)
    for i in 1:n, j in (i + 1):n
        overlaps(balls[i], balls[j]) && return nothing
    end
    perm = sortperm(balls)
    sorted = balls[perm]

    # coarsest rational near f that provably lies strictly between lo and hi
    # (either bound may be `nothing`)
    function certified_cut(f, lo, hi)
        for tol in (0.25, 1e-3, 1e-6, 1e-12, 1e-24, 1e-48)
            c = QQ(rationalize(BigInt, f; tol = tol))
            (lo === nothing || lo < c) || continue
            (hi === nothing || c < hi) || continue
            return c
        end
        nothing
    end

    mid(ball) = BigFloat(midpoint(ball))
    cuts = Vector{Union{Nothing, QQFieldElem}}(undef, n + 1)
    cuts[1] = certified_cut(mid(sorted[1]) - 1, nothing, sorted[1])
    cuts[n + 1] = certified_cut(mid(sorted[n]) + 1, sorted[n], nothing)
    for j in 1:(n - 1)
        cuts[j + 1] = certified_cut((mid(sorted[j]) + mid(sorted[j + 1])) / 2,
                                    sorted[j], sorted[j + 1])
    end
    any(isnothing, cuts) && return nothing

    intervals = Vector{Tuple{QQFieldElem, QQFieldElem}}(undef, n)
    for j in 1:n
        intervals[perm[j]] = (cuts[j], cuts[j + 1])
    end
    intervals
end

# The elements η_i = (θ - a_i)(θ - b_i) of [KR21] Algorithm 2, one per
# (embedding, interval) pair: η_i is negative under the i-th embedding and
# positive under every other real one, so products of η's realise any
# prescribed sign pattern. Deliberately not cached globally — the η's belong
# in the per-run _KRForm alongside diag and S.
function _kr_etas(K::AbsSimpleNumField, pairs)
    theta = gen(K)
    [(theta - a) * (theta - b) for (_, (a, b)) in pairs]
end

# --- square-class representatives over K --------------------------------------------

# A representative b of the square class of a, normalised as far as the
# arithmetic of K allows: b is integral, and when the returned flag is true
# the content ideal b·O_K is squarefree. Two obstructions make this weaker
# than the QQ squarefree part: writing b₀·O_K = C·I² with C squarefree, a
# representative with content C exists iff I is principal (flag false
# otherwise, integral representative returned unreduced); and a successful
# reduction fixes b only up to units-times-squares — over a formally real K
# that residual ambiguity is a sign pattern at the real places, whose control
# is [KR21] Algorithm 4's job, not this function's. Callers must not assume
# anything about which unit multiple they receive.
function _squarefree_content_rep(a::AbsSimpleNumFieldElem)
    iszero(a) && error("square-class representative of 0 requested")
    K = parent(a)
    OK = maximal_order(K)
    b = a * K(denominator(a, OK))^2          # integral, same square class
    fac = factor(ideal(OK, OK(b)))           # integral ideal: exponents all >= 0
    I = prod((P^div(e, 2) for (P, e) in fac if e >= 2); init = ideal(OK, 1))
    is_one(I) && return (b, true)
    fl, g = is_principal_with_data(I)
    fl || return (b, false)
    (divexact(b, K(g)^2), true)
end

# A SIZE-REDUCED representative of the square class of a — not a canonical
# one, because none exists over K (exponent reduction needs the leftover
# ideal class trivial; the unit group leaves the associate open), so
# isomorphic forms may still print differently. Takes the shortest, as
# printed, of three class-checked candidates: the input (never worse than
# before), _squarefree_content_rep (bites when the square content is
# principal — typical for the KR driver's outputs), and Hecke's
# reduce_mod_powers (bites on non-principal content and unit bloat; can
# return large-square-content inputs unchanged, which is why the content
# route is tried too). Applied ONLY at the user-facing boundary
# (anisotropic_part, sum_decomposition, diagonal_class over K) — never
# inside the KR algorithms, whose FacElem arithmetic exists precisely so
# nothing is evaluated.
function _size_reduced_square_class_rep(a::AbsSimpleNumFieldElem)
    iszero(a) && error("square-class representative of 0 requested")
    content_rep, _ = _squarefree_content_rep(a)
    hecke_rep = evaluate(Hecke.reduce_mod_powers(a, 2))
    best = a
    for b in (content_rep, hecke_rep)
        is_square(a * b) ||
            error("_size_reduced_square_class_rep postcondition failed: class changed")
        length(string(b)) < length(string(best)) && (best = b)
    end
    best
end
