# [KR21] Algorithms 5-8, Observation 5.3, and the anisotropic-part driver
# over a number field. Algorithms 5-7 are the reduction steps: each takes a
# form q with dim_a(q) = d and returns alpha in K^* with
# dim_a(q ⊥ <-alpha>) = d - 1. Algorithm 8 constructs the binary base case;
# Algorithms 7 and 8 are sequential siblings in the reduction chain —
# neither calls the other.
#
# [KR21] P. Koprowski and B. Rothkegel, *The anisotropic part of a quadratic
# form over a number field*, arXiv:2109.04172.

"""
    _kr_extend(beta::GWClass, alpha) -> GWClass

`beta ⊥ ⟨-alpha⟩`, with a consistency check at every extension: the signed
discriminant is multiplied by ``(-1)^n · (-alpha)`` where ``n`` is the rank
before extension, because the sign exponent moves from ``n(n-1)/2`` to
``(n+1)n/2`` — a difference of ``n``. Our representatives are well-defined
only up to squares (integralisation scales entries by denominator squares),
so the relation is asserted as a square-class identity. It costs nothing
and catches sign-convention drift in exactly the place drift would enter.
"""
function _kr_extend(beta::GWClass, alpha)
    K = gw_base_field(beta)
    a = K(alpha)
    iszero(a) && error("cannot extend by <0>")
    n = form_rank(beta)
    gamma = gw_direct_sum(beta, diagonal_form(K, -a))
    expected = K(iseven(n) ? 1 : -1) * (-a) * signed_discriminant(beta)
    is_square(expected * signed_discriminant(gamma)) ||
        error("signed-discriminant extension relation failed: sign-convention drift")
    gamma
end

"""
    _kr_algorithm_5(beta::GWClass) -> AbsSimpleNumFieldElem

[KR21] Algorithm 5: for ``dim_a(q) = d ≥ 4``, an ``α`` with
``dim_a(q ⊥ ⟨-α⟩) = d - 1``.
"""
function _kr_algorithm_5(beta::GWClass)
    K = gw_base_field(beta)
    d = anisotropic_dimension(beta)
    d >= 4 || error("anisotropic dimension of form is not >= 4")
    # Step (1): non-real K — over a non-real global field every form of
    # dimension >= 5 is isotropic, so d = 4, and parity gives
    # dim_a(q ⊥ <-1>) = 3.
    _is_formally_real(K) || return one(K)
    # Steps (2)-(3): the real embeddings and the extreme-signature index sets
    # (the η context threaded through the class's _KRForm).
    kr = _kr_form(beta)
    pairs = kr.pairs
    etas = kr.etas
    sigs = [form_signature(beta, sigma) for (sigma, _) in pairs]
    Iminus = [i for i in eachindex(sigs) if sigs[i] == -d]
    # Step (4): alpha positive on I+, negative on I-. KR leave the signs off
    # I+ ∪ I- unconstrained; Algorithm 2's eta-product happens to fix them
    # (positive everywhere off I-), a harmless over-specification nothing
    # downstream relies on.
    # Step (5): output.
    _kr_algorithm_2(K, pairs, etas, Iminus)
end

# The congruence targets shared by Algorithms 6 and 7: for every p in S,
# force ord_p(alpha) to the OPPOSITE parity of ord_p(disc q) — residue
# disc - 1 (mod p) makes alpha a p-unit where the valuation is odd, residue
# pi_p (mod p^2) makes ord_p(alpha) = 1 where it is even. Uses the
# S-normalised discriminant representative, so disc - 1 is an order element
# at the odd-valuation primes (memo §8). Returns (primes, exponents,
# targets, dn).
function _dim3_congruences(beta::GWClass)
    kr = _kr_form(beta)
    K = kr.K
    primes = copy(kr.S)
    dn = _disc_rep_normalized(signed_discriminant(beta), primes)
    exponents = Int[]
    targets = AbsSimpleNumFieldElem[]
    for P in primes
        if isodd(valuation(dn, P))
            push!(exponents, 1)
            push!(targets, dn - 1)
        else
            push!(exponents, 2)
            push!(targets, K(uniformizer(P)))
        end
    end
    (primes, exponents, targets, dn)
end

"""
    _kr_algorithm_6(beta::GWClass) -> AbsSimpleNumFieldElem

[KR21] Algorithm 6: over a NON-REAL number field with ``dim_a(q) = 3``, an
``α`` with ``dim_a(q ⊥ ⟨-α⟩) = 2``.

The two congruence branches do one thing: they force ``ord_p(α)`` to have
the **opposite parity** to ``ord_p(disc q)`` at every ``p ∈ S``, which is
what obstructs ``α ≡ disc(q)`` modulo squares locally — and that local
square-class identity is exactly what an anisotropic
``(q_a ⊥ ⟨-α⟩) ⊗ K_p`` would force, since the unique anisotropic quaternary
over ``K_p`` has trivial discriminant.
"""
function _kr_algorithm_6(beta::GWClass)
    K = gw_base_field(beta)
    _is_formally_real(K) &&
        error("Algorithm 6 requires a non-real field; use Algorithm 7")
    anisotropic_dimension(beta) == 3 || error("anisotropic dimension of form is not 3")
    OK = maximal_order(K)
    # Step (1): S = P(q) ∪ {dyadic primes} — the class's fixed _KRForm.
    primes, exponents, targets, dn = _dim3_congruences(beta)
    # Step (2): CRT for the congruences above, reduced modulo the product
    # ideal — Hecke's crt returns an unreduced representative whose ideal
    # would be infeasible to factor when the extension is re-diagonalised.
    # The reduction is safe for the parity property because the product of
    # the p^(k_p) is contained in every single p^(k_p), so the residue class
    # modulo each modulus survives, and that class pins ord_p(alpha):
    # k_p = 2 at the even-valuation primes, so the forced ord_p(alpha) = 1
    # survives the reduction; k_p = 1 at the odd-valuation ones, where
    # alpha ≡ dn - 1 ≡ -1 stays a p-unit.
    moduli = [primes[i]^exponents[i] for i in eachindex(primes)]
    residues = [OK(t) for t in targets]
    beta = length(residues) == 1 ? residues[1] : crt(residues, moduli)
    alpha = K(mod(beta, prod(moduli)))
    # Step (3): output, after checking the parity property the congruences
    # exist to produce.
    for P in primes
        isodd(valuation(alpha, P) + valuation(dn, P)) ||
            error("Algorithm 6 postcondition failed: parity at $P")
    end
    alpha
end

"""
    _kr_algorithm_7(beta::GWClass) -> AbsSimpleNumFieldElem

[KR21] Algorithm 7: over a FORMALLY REAL number field with ``dim_a(q) = 3``,
an ``α = α_1 · α_2`` with ``dim_a(q ⊥ ⟨-α⟩) = 2``.

The two factors are load-bearing at **disjoint** sets of places. ``α_1``
(Algorithm 4) carries ``sgn σ_i(α_1) = sgn(sgn σ_i(q))`` at the real places
and is a local square at every ``p ∈ S`` — so it contributes nothing at the
finite places. ``α_2`` (Algorithm 3) is totally positive — invisible at the
real places — and carries the Algorithm 6 congruences at ``S``. Dropping
``α_1`` therefore gives code that is correct at every finite place and
wrong at the real places, and it fails silently: every non-real test field
and every form of trivial signature still passes. That exact omission is
the published bug in the Macaulay2 package's dimension-3 branch, which is
why the two factors come from separate functions with separately asserted
postconditions and are combined in one visible line.

Algorithm 7 and Algorithm 8 are sequential siblings in the reduction chain
(7 drops ``dim_a`` from 3 to 2; 8 then constructs the binary anisotropic
part); neither calls the other.
"""
function _kr_algorithm_7(beta::GWClass)
    K = gw_base_field(beta)
    _is_formally_real(K) ||
        error("Algorithm 7 requires a formally real field; use Algorithm 6")
    anisotropic_dimension(beta) == 3 || error("anisotropic dimension of form is not 3")
    kr = _kr_form(beta)
    # Step (1): S = P(q) ∪ {dyadic primes} = kr.S.
    # Step (2): alpha_1 — the signature signs, local square at every p in S
    # (Algorithm 4 asserts both halves itself). dim_a = 3 forces odd rank,
    # so no real place has signature 0 and the index set is well-defined.
    pairs = kr.pairs
    etas = kr.etas
    I = [i for (i, (sigma, _)) in enumerate(pairs) if form_signature(beta, sigma) < 0]
    alpha1 = _kr_algorithm_4(K, pairs, etas, I, kr.S)
    # Step (3): alpha_2 — totally positive with the Algorithm 6 congruences
    # (Algorithm 3 asserts positivity and every congruence itself).
    primes, exponents, targets, dn = _dim3_congruences(beta)
    alpha2 = _kr_algorithm_3(K, primes, exponents, targets)
    # Step (4): the product, in one visible line.
    alpha = alpha1 * alpha2
    # Postconditions of the product, one per factor's contribution: the sign
    # pattern (alpha_2 > 0 cannot disturb it) ...
    for (i, (sigma, _)) in enumerate(pairs)
        sign(alpha, sigma) == (i in I ? -1 : 1) ||
            error("Algorithm 7 postcondition failed: wrong sign at real embedding $i")
    end
    # ... and the finite-place parity (alpha_1 has even valuation at S,
    # being a local square there).
    for P in primes
        isodd(valuation(alpha, P) + valuation(dn, P)) ||
            error("Algorithm 7 postcondition failed: parity at $P")
    end
    alpha
end
# --- Algorithm 8: the binary base case ----------------------------------------------

# (prod b^e, d)_P computed factor-wise from the FacElem — the Hilbert symbol
# is bimultiplicative, so only the odd-exponent bases contribute. Memoised
# per (base element, prime index): base elements recur across the columns of
# matrix B, and on a high-degree field evaluating an S-unit is the
# difference between seconds and never finishing.
function _hilbert_symbol_facelem(fe::FacElem, d, P, Pidx::Int,
                                 cache::Dict{Tuple{AbsSimpleNumFieldElem, Int}, Int})
    s = 1
    for (b, e) in fe.fac
        isodd(e) || continue
        s *= get!(() -> _hilbert_symbol(b, d, P), cache, (b, Pidx))
    end
    s
end

# The next `k` primes to append to S in step (5g): deterministic order —
# increasing norm (ties by rational prime, then position in the
# decomposition) — skipping the members of S.
function _next_primes_batch(K::AbsSimpleNumField, S::Vector, k::Int)
    OK = maximal_order(K)
    limit = ZZRingElem(20)
    while true
        cands = Tuple{ZZRingElem, ZZRingElem, Int, Any}[]
        p = ZZRingElem(2)
        while p <= limit
            for (idx, (P, _)) in enumerate(prime_decomposition(OK, p))
                any(Q -> Q == P, S) || push!(cands, (norm(P), p, idx, P))
            end
            p = next_prime(p)
        end
        sort!(cands; by = t -> (t[1], t[2], t[3]))
        length(cands) >= k && return [t[4] for t in cands[1:k]]
        limit *= 2
    end
end

# alpha = prod_{eps_j = 1} beta_j, evaluated ONCE after merging the factor
# dictionaries and reducing every exponent mod 2 — only the square class of
# alpha matters, and the mod-2 reduction is what keeps the output small.
function _evaluate_mod_squares(K::AbsSimpleNumField, basis, sol)
    bases = AbsSimpleNumFieldElem[]
    exps = ZZRingElem[]
    for j in eachindex(basis)
        isone(sol[j, 1]) || continue
        for (b, e) in basis[j].fac
            _push_factor!(bases, exps, b, e)
        end
    end
    alpha = one(K)
    for (b, e) in zip(bases, exps)
        isodd(e) && (alpha *= b)
    end
    alpha
end

# Steps (1)-(4) of [KR21] Algorithm 8: the padding, the discriminant, S, and
# the real embeddings that impose conditions.
function _kr_algorithm_8_setup(beta::GWClass)
    K = gw_base_field(beta)
    anisotropic_dimension(beta) == 2 || error("anisotropic dimension of form is not 2")
    # Step (1): pad the Witt index to 0 mod 4. NOT cosmetic, do not drop: KR's
    # correctness proof uses the divisibility in exactly one line at the very
    # end — it is what makes the Hasse invariants of q and of the candidate
    # <alpha, -alpha d> agree at the DYADIC primes, where the hyperbolic
    # correction factor need not be trivial. Dropping it corrupts only dyadic
    # Hasse data, silently (see the padding regression test).
    w = gw_witt_index(beta)
    q = beta
    if w % 4 != 0
        q = gw_direct_sum(q, hyperbolic_form(K, 2 * (4 - w % 4)))
    end
    kr = _kr_form(q)
    # Step (2): the signed discriminant of the (padded) form.
    d = signed_discriminant(q)
    # Step (3): S = P(q) ∪ {dyadic primes} — the fixed _KRForm data. The loop
    # grows a local copy; no diagonalisation is ever re-derived.
    S = copy(kr.S)
    # Step (4): only the real embeddings with sigma(d) < 0 impose conditions
    # ([KR21] Lemma 5.2: over a real closed field <a,b> ⊥ w×H is hyperbolic
    # exactly when its discriminant is positive, so sigma(d) > 0 places are
    # automatically similar to the candidate).
    negd = [sigma for (sigma, _) in kr.pairs if sign(d, sigma) < 0]
    for sigma in negd
        abs(form_signature(q, sigma)) == 2 ||
            error("Algorithm 8: signature is not ±2 at an embedding with negative discriminant")
    end
    (q = q, d = d, S = S, negd = negd)
end

# Step (5) of [KR21] Algorithm 8: the retry loop. Matrices A and B are
# filled FACTOR-WISE from the FacElem representation — sign(FacElem, σ) is
# Hecke-native odd-exponent sign multiplication, and the Hilbert symbol
# column uses _hilbert_symbol_facelem — no S-unit is ever evaluated.
function _kr_algorithm_8_loop(q::GWClass, d, S::Vector, negd::Vector;
                              max_rounds::Int = 8)
    K = gw_base_field(q)
    symcache = Dict{Tuple{AbsSimpleNumFieldElem, Int}, Int}()
    for _ in 1:max_rounds
        # Step (5a): F2-basis of E_S (Algorithm 1, FacElem end to end).
        basis = _kr_algorithm_1(K, S)
        m = length(basis)
        # Steps (5b)-(5c): the right-hand side from the signatures and the
        # Hasse invariants of the padded form.
        v = Int[form_signature(q, sigma) == -2 ? 1 : 0 for sigma in negd]
        wv = Int[hasse_witt_invariant(q, P) == -1 ? 1 : 0 for P in S]
        # Steps (5d)-(5e): the sign block A and the symbol block B.
        rows = length(negd) + length(S)
        AB = zero_matrix(GF(2), rows, m)
        for j in 1:m
            fe = basis[j]
            for (i, sigma) in enumerate(negd)
                sign(fe, sigma) < 0 && (AB[i, j] = 1)
            end
            for (i, P) in enumerate(S)
                _hilbert_symbol_facelem(fe, d, P, i, symcache) == -1 &&
                    (AB[length(negd) + i, j] = 1)
            end
        end
        # Step (5f): solve over F2 and exit on success.
        ok, sol = can_solve_with_solution(AB, matrix(GF(2), rows, 1, vcat(v, wv));
                                          side = :right)
        if ok
            alpha = _evaluate_mod_squares(K, basis, sol)
            # Step (6): the binary anisotropic part.
            qa = diagonal_form(K, (alpha, -alpha * d))
            anisotropic_dimension(qa) == 2 ||
                error("Algorithm 8 postcondition failed: output is not anisotropic")
            return qa
        end
        # Step (5g): enlarge S in a BATCH. Each retry recomputes
        # sunit_group_fac_elem(S) from scratch inside Algorithm 1 — the
        # dominant cost — so the number of retries matters, not the number of
        # primes appended per retry; KR themselves append three at once in §6.
        append!(S, _next_primes_batch(K, S, 3))
    end
    error("Algorithm 8: the F2 system stayed unsolvable after $max_rounds rounds " *
          "(|S| = $(length(S))). Lemma 5.1 guarantees termination, so raise " *
          "max_rounds or investigate the E_S basis")
end

"""
    _kr_algorithm_8(beta::GWClass; max_rounds = 8) -> GWClass

[KR21] Algorithm 8: for ``dim_a(q) = 2``, the binary anisotropic part
``⟨α, -α·d⟩`` directly, split at the paper's own seam — setup (steps
(1)-(4)) and the retry loop (step (5), with the step-(6) output).

Termination is [KR21] Lemma 5.1: some anisotropic part of ``q`` has
coefficients that are ``(S ∪ {𝔮})``-singular for a single additional place
``𝔮``, so after finitely many enlargements of S the F₂ system acquires a
solution; `max_rounds` merely converts a bug into an informative error
instead of a spin.
"""
_kr_algorithm_8(beta::GWClass; max_rounds::Int = 8) =
    (s = _kr_algorithm_8_setup(beta);
     _kr_algorithm_8_loop(s.q, s.d, s.S, s.negd; max_rounds = max_rounds))

"""
    _kr_observation_5_3(beta::GWClass) -> GWClass

[KR21] Observation 5.3: if ``dim_a(q) = 1``, then ``⟨disc(q)⟩`` is an
anisotropic part of ``q``.
"""
function _kr_observation_5_3(beta::GWClass)
    K = gw_base_field(beta)
    diagonal_form(K, (K(signed_discriminant(beta)),))
end

# --- the driver ----------------------------------------------------------------------

"""
    _anisotropic_part_nf(beta::GWClass) -> GWClass

The anisotropic part of a form over a number field, assembling [KR21]
Algorithms 5-8 and Observation 5.3: Algorithm 5 reduces while
``dim_a ≥ 4``; Algorithm 6 (non-real) or 7 (formally real) at
``dim_a = 3``; Algorithm 8 constructs the binary part at ``dim_a = 2``;
Observation 5.3 handles ``dim_a = 1``; the empty form is ``dim_a = 0``.
Every ``⟨-α⟩`` dropped by a reduction step contributes ``⟨α⟩`` to the
anisotropic part and is reassembled at the end, exactly as [KR21] §6 ends
with ``q_a = ⟨1406⟩ ⊥ q'_a``.
"""
function _anisotropic_part_nf(beta::GWClass)
    K = gw_base_field(beta)
    coeffs = AbsSimpleNumFieldElem[]
    q = beta
    while anisotropic_dimension(q) >= 4
        alpha = _kr_algorithm_5(q)
        push!(coeffs, K(alpha))
        q = _kr_extend(q, alpha)
    end
    if anisotropic_dimension(q) == 3
        alpha = _is_formally_real(K) ? _kr_algorithm_7(q) : _kr_algorithm_6(q)
        push!(coeffs, K(alpha))
        q = _kr_extend(q, alpha)
    end
    da = anisotropic_dimension(q)
    core = if da == 2
        _kr_algorithm_8(q)
    elseif da == 1
        _kr_observation_5_3(q)
    else
        diagonal_form(K, ())
    end
    isempty(coeffs) && return core
    gw_direct_sum(diagonal_form(K, Tuple(coeffs)), core)
end
