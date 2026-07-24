# [KR21] Algorithms 2-4: elements with prescribed signs at the real
# embeddings, optionally with congruence and local-square conditions at a
# finite set of primes. Dependencies are explicit in every signature; the
# (embedding, interval) pairs and the eta vector come from
# _real_embeddings_with_intervals / _kr_etas (nf_arithmetic.jl).
#
# [KR21] P. Koprowski and B. Rothkegel, *The anisotropic part of a quadratic
# form over a number field*, arXiv:2109.04172.

"""
    _kr_algorithm_2(K::AbsSimpleNumField, pairs, etas, I)
        -> AbsSimpleNumFieldElem

[KR21] Algorithm 2: an element ``ρ`` with ``σ_i(ρ) < 0`` for `i ∈ I` and
``σ_i(ρ) > 0`` otherwise, where the real embeddings ``σ_i`` are the first
components of `pairs`. `pairs` and `etas` must come from the same
[`_real_embeddings_with_intervals`](@ref) / [`_kr_etas`](@ref) call; `I` is a
collection of pair indices.

It works because at ``ξ_i = σ_i(θ)`` the two linear factors of ``η_i`` have
opposite signs, while at ``ξ_j``, ``j ≠ i``, they have the same sign — which
needs the isolating intervals pairwise disjoint, guaranteed upstream but
re-asserted here because a violation silently flips signs in matrix A of
[KR21] Algorithm 8.
"""
function _kr_algorithm_2(K::AbsSimpleNumField, pairs::Vector, etas::Vector, I)
    length(pairs) == length(etas) ||
        error("pairs and etas must come from the same construction")
    Iset = Set{Int}(I)
    all(i -> 1 <= i <= length(pairs), Iset) || error("index set out of range")
    intervals = sort([iv for (_, iv) in pairs]; by = first)
    for j in 1:(length(intervals) - 1)
        intervals[j][1] < intervals[j][2] || error("empty isolating interval")
        intervals[j][2] <= intervals[j + 1][1] || error("isolating intervals overlap")
    end

    # Steps (1)-(2), the real roots xi_i and their isolating intervals, and
    # step (3), the eta_i, were computed upstream (see the docstring).
    # Step (4): rho = prod_{i in I} eta_i.
    rho = isempty(Iset) ? one(K) : prod(etas[i] for i in Iset)

    for (i, (sigma, _)) in enumerate(pairs)
        sign(rho, sigma) == (i in Iset ? -1 : 1) ||
            error("Algorithm 2 postcondition failed: wrong sign at real embedding $i")
    end
    rho
end

"""
    _kr_algorithm_3(K::AbsSimpleNumField, primes, exponents, lambdas)
        -> AbsSimpleNumFieldElem

[KR21] Algorithm 3: a totally positive ``α ∈ K`` with
``α ≡ λ_i (mod 𝔭_i^{k_i})`` for every `i`. `primes` are distinct prime
ideals of the maximal order, `exponents` the ``k_i ≥ 1``, and `lambdas` the
targets, which must be algebraic integers (callers clear denominators by a
square first — that is the integrality repair of memo §8).
"""
function _kr_algorithm_3(K::AbsSimpleNumField, primes::Vector, exponents::Vector{Int},
                         lambdas::Vector)
    length(primes) == length(exponents) == length(lambdas) ||
        error("primes, exponents, and lambdas must have equal length")
    isempty(primes) && return one(K)
    OK = maximal_order(K)
    all(P -> order(P) === OK && is_prime(P), primes) ||
        error("congruences must be at prime ideals of the maximal order of K")
    for i in eachindex(primes), j in (i + 1):length(primes)
        primes[i] == primes[j] && error("repeated prime in the congruence list")
    end
    all(>=(1), exponents) || error("congruence exponents must be positive")
    targets = [K(lam) for lam in lambdas]
    all(lam -> isone(denominator(lam, OK)), targets) ||
        error("congruence targets must be algebraic integers")

    # Step (1): CRT for beta = lambda_i mod p_i^(k_i). Hecke's crt returns an
    # UNREDUCED representative whose coefficients blow up with the number of
    # congruences; reduce modulo the product ideal (contained in every
    # p_i^(k_i), so all congruences survive) or the element becomes too large
    # to factor downstream.
    moduli = [primes[i]^exponents[i] for i in eachindex(primes)]
    residues = [OK(lam) for lam in targets]
    beta = length(residues) == 1 ? residues[1] : crt(residues, moduli)
    beta = mod(beta, prod(moduli))

    # Steps (2)-(3): the rational primes below the p_i and, for each,
    # m(p) = max{k_i : p_i lies over p}.
    m = Dict{ZZRingElem, Int}()
    for (P, k) in zip(primes, exponents)
        p = minimum(P)
        m[p] = max(get(m, p, 0), k)
    end

    # Step (4): s = prod p^m(p). Since ord_{p_i}(s) = e(p_i|p) m(p) >= k_i,
    # adding any integer multiple of s preserves every congruence.
    s = prod(p^k for (p, k) in m)

    # Steps (5)-(6): a positive integer t with t*s > sigma_j(-beta) for every
    # real embedding, then alpha = beta + t*s. s is a rational integer, so
    # sigma_j(t*s) = t*s at EVERY real embedding — the shift is the same in
    # every real place, which is the whole reason growing t alone makes alpha
    # totally positive.
    t = ZZRingElem(1)
    alpha = K(beta) + t * s
    while !is_totally_positive(alpha)
        t *= 2
        alpha = K(beta) + t * s
    end

    is_totally_positive(alpha) ||
        error("Algorithm 3 postcondition failed: not totally positive")
    for (P, k, lam) in zip(primes, exponents, targets)
        diff = alpha - lam
        (iszero(diff) || valuation(diff, P) >= k) ||
            error("Algorithm 3 postcondition failed: congruence at $P")
    end
    alpha
end

"""
    _kr_algorithm_4(K::AbsSimpleNumField, pairs, etas, I, S::Vector)
        -> AbsSimpleNumFieldElem

[KR21] Algorithm 4: an element ``ρ`` with ``σ_i(ρ) < 0`` exactly for `i ∈ I`
that is additionally a local square at every prime in `S`. `pairs`/`etas` as
in [`_kr_algorithm_2`](@ref); `S` a vector of prime ideals of the maximal
order.

Implements a **corrected** congruence modulus (the paper's exponent is too
small whenever ``ord_p(α₁) > 0``); see the comment at the exponent
computation below.
"""
function _kr_algorithm_4(K::AbsSimpleNumField, pairs::Vector, etas::Vector, I, S::Vector)
    OK = maximal_order(K)

    # Step (1): alpha_1 with the requested sign pattern (Algorithm 2).
    alpha1 = _kr_algorithm_2(K, pairs, etas, I)
    # Not in KR (integrality repair, memo §8): the eta's have rational
    # interval endpoints, so alpha_1 need not be integral, but Algorithm 3
    # requires integral congruence targets. Scaling by t^2 > 0 changes
    # neither the signs nor the square class.
    alpha1 *= K(denominator(alpha1, OK))^2

    # CORRECTION to [KR21] (memo §8, Correction A). KR use the modulus
    # exponent 1 + ord_p(4) and conclude rho = alpha_1 * alpha_2
    # = alpha_1^2 * (alpha_2/alpha_1) is a local square at p. The Local
    # Square Theorem needs
    #     ord_p(alpha_2/alpha_1 - 1) = ord_p(alpha_2 - alpha_1) - ord_p(alpha_1)
    #                                >= 1 + ord_p(4),
    # so the exponent must be k_p = 1 + ord_p(4) + ord_p(alpha_1); KR's works
    # only when alpha_1 is a unit at p, which nothing establishes — its norm
    # is free to be divisible by primes of S, the dyadic ones included.
    # ord_p(alpha_1) is read AFTER the t^2 scaling above: the scaling raises
    # it by 2 e(p|2) per factor of 2 in t at exactly the dyadic primes, so
    # reading it earlier would reintroduce the gap the correction closes.
    ks = Int[1 + valuation(K(4), P) + valuation(alpha1, P) for P in S]

    # Step (2): a totally positive alpha_2 = alpha_1 mod p^(k_p) (Algorithm 3),
    # hence alpha_2/alpha_1 is a local square at every p in S.
    alpha2 = _kr_algorithm_3(K, S, ks, [alpha1 for _ in S])

    # Step (3): rho = alpha_1 * alpha_2.
    rho = alpha1 * alpha2

    # Postconditions, separately: the signs are alpha_1's contribution and the
    # local squareness alpha_2's — a half-correct rho must not escape (the
    # published Macaulay2 defect is exactly the first failing while the second
    # holds).
    Iset = Set{Int}(I)
    for (i, (sigma, _)) in enumerate(pairs)
        sign(rho, sigma) == (i in Iset ? -1 : 1) ||
            error("Algorithm 4 postcondition failed: wrong sign at real embedding $i")
    end
    for P in S
        is_local_square(rho, P) ||
            error("Algorithm 4 postcondition failed: not a local square at $P")
    end
    rho
end
