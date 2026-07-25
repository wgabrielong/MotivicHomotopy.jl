# Anisotropic dimensions over number fields (KC18 Algorithm 9 shape with the
# global discriminant correction) and, below, the KR reduction algorithms.
# Includes the memo §10 boundary verification of the max-over-places formula.

using Test
using Oscar
using MotivicHomotopy

const MH = MotivicHomotopy

# The shared number-field test matrix (see nf_test_fields.jl).
include("nf_test_fields.jl")

# Independent isotropy through Hecke's LOCAL theory. Hecke 0.39.21's GLOBAL
# is_isotropic is wrong over number fields (it reports the definite <1,1,1>
# over Q(sqrt 2) isotropic; the QQ specialisation is correct), so the memo
# §10 boundary check is built from the local predicates: for rank >= 3 a
# form is isotropic iff isotropic at every real place and every P in S (off
# S it is a unit form at an odd prime, isotropic in rank >= 3); rank 2 is
# the exact global square test; rank <= 1 is never isotropic. Only the
# place list is shared with our S-machinery — the local computations are
# Hecke's own.
function _independent_isotropic(beta)
    M = gw_matrix(beta)
    K = base_ring(M)
    r = form_rank(beta)
    r <= 1 && return false
    r == 2 && return is_square(-det(M))
    V = quadratic_space(K, M)
    all(is_isotropic(V, v) for v in real_places(K)) || return false
    all(is_isotropic(V, P) for P in MH._kr_form(beta).S)
end

@testset "local anisotropic dimensions" begin
    OK7 = maximal_order(K7)
    P7 = prime_decomposition(OK7, 7)[1][1]
    q4 = diagonal_form(K7, (1, -3, -s7, 3 * s7))    # <1,-u,-pi,u*pi> at P7
    @test MH._local_anisotropic_dimension(q4, P7) == 4
    @test is_square(signed_discriminant(q4))         # the dim_a = 4 case has trivial disc
    @test anisotropic_dimension(q4) == 4
    @test is_anisotropic_form(q4)

    h = hyperbolic_form(K7, 4)
    @test MH._local_anisotropic_dimension(h, P7) == 0
    @test anisotropic_dimension(h) == 0
    @test gw_witt_index(h) == 2

    beta = diagonal_form(K2, (1, 1, 1))
    for sigma in real_places(K2)
        @test MH._local_anisotropic_dimension(beta, sigma) == 3
    end
    cplx = [P for P in infinite_places(K7) if !is_real(P)]
    @test MH._local_anisotropic_dimension(diagonal_form(K7, (1, 1, 1)), cplx[1]) == 1
    @test MH._local_anisotropic_dimension(diagonal_form(K7, (1, 1)), cplx[1]) == 0
end

@testset "global dimension: the narrow-class blind spot" begin
    # memo §10 finding (i): eps is a totally positive non-square unit that IS
    # a local square at the unique dyadic prime, so <1,-eps> is locally
    # isotropic at every place of real ∪ S yet globally anisotropic; the
    # global discriminant correction detects it, the max over real ∪ S alone
    # would not
    eps = 35 + 6 * s34
    OK34 = maximal_order(K34)
    P2_34 = prime_decomposition(OK34, 2)[1][1]
    @test is_totally_positive(eps)
    @test is_unit(OK34(eps))
    @test !is_square(eps)
    @test is_local_square(eps, P2_34)
    blind = diagonal_form(K34, (1, -eps))
    @test MH._local_anisotropic_dimension(blind, P2_34) == 0
    @test all(sigma -> MH._local_anisotropic_dimension(blind, sigma) == 0,
              real_places(K34))
    @test anisotropic_dimension(blind) == 2
    @test is_anisotropic_form(blind)
    @test !_independent_isotropic(blind)
    # and the square counterpart really is hyperbolic
    @test anisotropic_dimension(diagonal_form(K34, (1, -eps^2))) == 0
end

@testset "memo §10 boundary check: dima < rank vs Hecke local isotropy" begin
    fixtures = GWClass[
        diagonal_form(K2, (1, 1, 1)),
        diagonal_form(K2, (1, -1, 1)),
        diagonal_form(K2, (1, s2, 2, 3 - s2)),
        diagonal_form(K2, (s2, -s2)),
        diagonal_form(K2, (1, -2)),
        diagonal_form(K5, (1, 1, 1)),
        diagonal_form(K5, (s5, 3, 1, -2)),
        diagonal_form(K5, (1, 5)),
        diagonal_form(K7, (1, -3, -s7, 3 * s7)),
        diagonal_form(K7, (s7, 3, 1, -2)),
        diagonal_form(K7, (1, 1, 1, 1, 1)),
        diagonal_form(K7,
            (-3 - 9s7, -1, -2 - 6s7, 1 - s7, -6 + 4s7, -3 + 2s7, 4 - 4s7)),
        hyperbolic_form(K7, 6),
        diagonal_form(Kc, (1, ac, ac - 1)),
        diagonal_form(Kc, (1, 1, 1, 1)),
        diagonal_form(Kc, (ac, -ac, 2, 3)),
        diagonal_form(K34, (1, -(35 + 6 * s34))),
        diagonal_form(K34, (1, 1, 1)),
        diagonal_form(K34, (s34, 2, -1)),
        GWClass(map_entries(x -> KQ(x), QQ[1 4 7; 4 3 -1; 7 -1 5])),
        GWClass(map_entries(x -> KQ(x), QQ[2 1; 1 3])),
        # the rest of the field matrix:
        diagonal_form(Kr5, (1, sr5, 2)),
        diagonal_form(K21, (1, s21)),
        diagonal_form(K32, (1, a32, a32^2 - 3)),
        diagonal_form(Kq, (aq, 1 - aq)),
    ]
    for beta in fixtures
        r = form_rank(beta)
        d = anisotropic_dimension(beta)
        @test 0 <= d <= r
        @test iseven(r - d)
        @test gw_witt_index(beta) == div(r - d, 2)
        @test _independent_isotropic(beta) == (d < r)
        @test is_isotropic_form(beta) == (d < r)
    end
end

# --- KR Algorithms 5-7: the reduction chain ------------------------------------------


# one reduction step, dispatched by dim_a and formal reality as the driver does
function _reduce_once(beta)
    d = anisotropic_dimension(beta)
    K = gw_base_field(beta)
    alpha = d >= 4 ? MH._kr_algorithm_5(beta) :
            (MH._is_formally_real(K) ? MH._kr_algorithm_7(beta) :
                                       MH._kr_algorithm_6(beta))
    MH._kr_extend(beta, alpha)
end

@testset "_kr_extend: signed-discriminant extension relation" begin
    beta = diagonal_form(K2, (1, s2, 3))
    n = form_rank(beta)
    for alpha in (K2(2), s2, -1 - s2)
        gamma = MH._kr_extend(beta, alpha)      # the relation is asserted inside
        @test form_rank(gamma) == n + 1
        @test is_square(signed_discriminant(gamma) *
                        (K2((-1)^n) * (-alpha) * signed_discriminant(beta)))
    end
    @test_throws ErrorException MH._kr_extend(beta, K2(0))
end

@testset "KR Algorithm 5 (dim_a >= 4)" begin
    cases = GWClass[
        diagonal_form(K7, (1, -3, -s7, 3 * s7)),      # non-real, dima 4
        diagonal_form(K2, (1, 1, 1, 1)),              # definite
        diagonal_form(K2, (-1, -1, -1, -1)),          # negative definite
        diagonal_form(K2, (1, 1, 1, 1, 1)),           # dima 5
        diagonal_form(Kc, (1, 1, 1, 1)),              # totally real cubic
        diagonal_form(Kr5, (1, 1, 1, 1)),
        diagonal_form(K32, (1, 1, 1, 1)),             # mixed signature
    ]
    for q in cases
        d = anisotropic_dimension(q)
        @test d >= 4
        alpha = MH._kr_algorithm_5(q)
        @test anisotropic_dimension(MH._kr_extend(q, alpha)) == d - 1
    end
    # non-real: output is 1 (d = 4 is forced; parity does the rest)
    @test MH._kr_algorithm_5(diagonal_form(K7, (1, -3, -s7, 3 * s7))) == one(K7)
    @test_throws ErrorException MH._kr_algorithm_5(diagonal_form(K2, (1, 1, 1)))
end

@testset "KR Algorithm 6 (dim_a = 3, non-real): KR §6 intermediate" begin
    q6 = diagonal_form(K7,
        (-3 - 9s7, -1, -2 - 6s7, 1 - s7, -6 + 4s7, -3 + 2s7, 4 - 4s7))
    @test anisotropic_dimension(q6) == 3
    # the paper's prime data: P(q) ∪ dyadics = both dyadic primes, (3), and
    # one prime above 37 — four primes, exactly as their §6 states
    @test sort([minimum(P) for P in MH._kr_form(q6).S]) == ZZRingElem[2, 2, 3, 37]
    alpha = MH._kr_algorithm_6(q6)
    @test anisotropic_dimension(MH._kr_extend(q6, alpha)) == 2
    # KR's own alpha is 1406; ours is a different CRT representative (6290 at
    # the time of writing — not canonical, so not pinned). Both must work:
    @test anisotropic_dimension(MH._kr_extend(q6, K7(1406))) == 2

    # further non-real inputs: Q(sqrt-7) unit form, Q(sqrt-5) (h = 2, ramified
    # dyadic prime)
    for q in (diagonal_form(K7, (1, 1, 1)),
              diagonal_form(K5, (1, 2, s5)),
              diagonal_form(K5, (1, s5, 2)))
        @test anisotropic_dimension(q) == 3
        a = MH._kr_algorithm_6(q)
        @test anisotropic_dimension(MH._kr_extend(q, a)) == 2
    end

    @test_throws ErrorException MH._kr_algorithm_6(diagonal_form(K2, (1, 1, 1)))
end

@testset "KR Algorithm 7 (dim_a = 3, real) and the alpha_1 regression" begin
    cases = GWClass[
        diagonal_form(K2, (1, 1, 1)),
        diagonal_form(K2, (-1, -1, -1)),
        diagonal_form(Kr5, (1, 1, 1)),
        diagonal_form(Kc, (-1, -1, -1)),
        diagonal_form(K32, (1, 1, 1)),                # mixed signature field
        diagonal_form(K32, (-1, -1, -1)),
    ]
    for q in cases
        @test anisotropic_dimension(q) == 3
        alpha = MH._kr_algorithm_7(q)
        @test anisotropic_dimension(MH._kr_extend(q, alpha)) == 2
    end

    # THE regression, explicitly: with sgn σ(q) = -3 somewhere, alpha_2 alone
    # carries every finite-place condition but the wrong real-place signs, and
    # the reduction fails LOUDLY on this input (dim_a stays 4, not 2). This is
    # the published Macaulay2 defect shape; it is invisible on non-real fields
    # and on forms of trivial signature.
    q = diagonal_form(K2, (-1, -1, -1))
    primes, exponents, targets, _ = MH._dim3_congruences(q)
    alpha2 = MH._kr_algorithm_3(K2, primes, exponents, targets)
    @test is_totally_positive(alpha2)                  # hence sign-wrong for this q
    @test anisotropic_dimension(MH._kr_extend(q, alpha2)) == 4
    @test anisotropic_dimension(MH._kr_extend(q, alpha2)) != 2

    @test_throws ErrorException MH._kr_algorithm_7(diagonal_form(K7, (1, 1, 1)))
end

@testset "reduction chain: dim_a = n down to 2, one step at a time" begin
    starts = GWClass[
        diagonal_form(K2, (1, 1, 1, 1, 1, 1)),        # 6 -> 2 through Alg 5 then 7
        diagonal_form(Kc, (-1, -1, -1, -1, -1)),      # 5 -> 2
        diagonal_form(K7, (1, -3, -s7, 3 * s7)),      # 4 -> 2 through Alg 5 then 6
        diagonal_form(K7,
            (-3 - 9s7, -1, -2 - 6s7, 1 - s7, -6 + 4s7, -3 + 2s7, 4 - 4s7)),
        diagonal_form(K32, (1, 1, 1, 1)),             # mixed signature, 4 -> 2
    ]
    for q0 in starts
        q = q0
        d = anisotropic_dimension(q)
        while d > 2
            q = _reduce_once(q)
            dnew = anisotropic_dimension(q)
            @test dnew == d - 1
            @test dnew >= 2
            d = dnew
        end
        @test d == 2
    end
end

# --- KR Algorithm 8, Observation 5.3, and the driver ---------------------------------

@testset "KR §6 fixture, end to end" begin
    # Every §6 assertion in one place (consolidated: the fixture values moved
    # here from test_number_fields.jl; the former "paper oracle" and
    # "AMENDMENT 1" testsets merged in). Comparisons use is_isomorphic_form,
    # never == — the anisotropic part is unique only up to isometry, and
    # Algorithm 8 depends on a choice of E_S basis and of solution vector.
    q6 = diagonal_form(K7,
        (-3 - 9s7, -1, -2 - 6s7, 1 - s7, -6 + 4s7, -3 + 2s7, 4 - 4s7))

    # the paper's printed fixture values
    @test signed_discriminant(q6) == -61056 - 342912 * s7    # disc q, literally
    @test sort([minimum(P) for P in MH._kr_form(q6).S]) == ZZRingElem[2, 2, 3, 37]
    # non-real K: #{p : s_p = -1} is even; for the Algorithm 8 input (padded
    # to Witt index 0 mod 4) it is exactly two. WHICH two primes was verified
    # three independent ways (Hecke; hand product; the Q2 oracle through the
    # split completion): {one dyadic prime, one prime over 37} — the paper's
    # (♣) row enumeration is NOT the listing order d1,d2,p1..p7, so only the
    # canonical count and the verified prime set are asserted.
    negq = [P for P in MH._kr_form(q6).S if hasse_witt_invariant(q6, P) == -1]
    @test iseven(length(negq))
    q6pad = gw_direct_sum(gw_direct_sum(q6, diagonal_form(K7, -1406)),
                          hyperbolic_form(K7, 2))
    negpad = [P for P in MH._kr_form(q6pad).S if hasse_witt_invariant(q6pad, P) == -1]
    @test length(negpad) == 2
    @test sort([minimum(P) for P in negpad]) == ZZRingElem[2, 37]

    # the dim_a chain: paper's dim_a(q) = 3, dropped to 2 by <-1406>
    @test anisotropic_dimension(q6) == 3
    @test gw_witt_index(q6) == 2
    @test is_isotropic_form(q6)
    qp = gw_direct_sum(q6, diagonal_form(K7, -1406))
    @test anisotropic_dimension(qp) == 2
    @test gw_witt_index(qp) == 3

    # the alpha_KR chain of the paper's section 6, in order:
    aKR = (-27 - 19s7) // 2
    # (a) rank 7 has sign exponent 21 (odd), rank 8 has 28 (even), so the
    # extension flips the sign and disc(q') = -(-1406) * disc(q) — exactly.
    @test signed_discriminant(qp) == 1406 * signed_discriminant(q6)
    # (b) the paper's printed third coefficient is LITERAL arithmetic:
    # -alpha_KR * disc(q') equals it exactly (verified — no square factor).
    @test -aKR * signed_discriminant(qp) == 30903025152 - 7324337664 * s7
    # (c) the PAPER's answer is a correct anisotropic part, by our machinery
    # and independent of any choice we make.
    qaKR = diagonal_form(K7, (aKR, -aKR * signed_discriminant(qp)))
    @test anisotropic_dimension(qaKR) == 2
    @test is_anisotropic_form(qaKR)
    @test is_isomorphic_form(qp, gw_direct_sum(qaKR, hyperbolic_form(K7, 6)))
    # (d) the acceptance test for Algorithm 8 on §6: OUR alpha yields a form
    # isometric to THEIRS — independent of the E_S basis, the solution
    # vector, and any row order. (If (d) fails while (c) passes, the fault is
    # localised to the alpha construction.)
    qa8 = MH._kr_algorithm_8(qp)
    @test is_isomorphic_form(qa8, qaKR)
    # (e) the assembled driver ends exactly as the paper does:
    # q_a = <1406> ⊥ q'_a, of dimension 3.
    t6 = @elapsed qa_full = anisotropic_part(q6)
    @info "KR §6 driver wall-clock" seconds = round(t6; digits = 3)
    @test form_rank(qa_full) == 3
    @test is_isomorphic_form(qa_full, gw_direct_sum(diagonal_form(K7, 1406), qaKR))

    # the paper oracle for the full published part (hecke-trust.md §2(b)):
    # q_a = <1406, (-27-19√-7)/2, 30903025152-7324337664√-7>; the last line
    # also checks Hecke's is_isometric, an independent code path from the
    # broken global is_isotropic (invariant comparison, no represents()).
    qa_paper = diagonal_form(K7,
        (1406, (-27 - 19s7) // 2, 30903025152 - 7324337664s7))
    @test is_anisotropic_form(qa_paper)
    @test is_isomorphic_form(qa_full, qa_paper)
    decomp = gw_direct_sum(qa_paper, hyperbolic_form(K7, 4))
    @test is_isomorphic_form(q6, decomp)
    @test is_isometric(quadratic_space(K7, gw_matrix(q6)),
                       quadratic_space(K7, gw_matrix(decomp)))
end

@testset "AMENDMENT 2: the padding step has its own regression input" begin
    # Sweep result (recorded): the §6 form does NOT exercise step (1) — its
    # dyadic Hasse set is identical padded and unpadded — but the corpus
    # sweep found witnesses at every residue w mod 4 in {1,2,3}. First one:
    # <1, sqrt2> ⊥ H over Q(sqrt2), where the dyadic Hasse invariant is -1
    # unpadded and +1 padded, so a driver without step (1) targets the wrong
    # dyadic invariant.
    OK2 = maximal_order(K2)
    P2 = prime_decomposition(OK2, 2)[1][1]
    q0 = diagonal_form(K2, (1, s2))
    for w in (1, 3)                                     # w mod 4 = 1 and 3
        q = gw_direct_sum(q0, hyperbolic_form(K2, 2w))
        @test anisotropic_dimension(q) == 2 && gw_witt_index(q) % 4 == w % 4
        qpad = gw_direct_sum(q, hyperbolic_form(K2, 2 * (4 - w)))
        @test hasse_witt_invariant(q, P2) == -1
        @test hasse_witt_invariant(qpad, P2) == 1       # padding changes it
        qa = anisotropic_part(q)
        # the true anisotropic part carries the PADDED form's dyadic Hasse
        # invariant, not the unpadded one — step (1) is load-bearing
        @test hasse_witt_invariant(qa, P2) == hasse_witt_invariant(qpad, P2)
        @test hasse_witt_invariant(qa, P2) != hasse_witt_invariant(q, P2)
        @test is_isomorphic_form(q, gw_direct_sum(qa, hyperbolic_form(K2, 2 * (w + 0))))
    end
    # w mod 4 = 2, over the split-dyadic field
    OK7 = maximal_order(K7)
    q7 = gw_direct_sum(diagonal_form(K7, (1, s7)), hyperbolic_form(K7, 4))
    @test gw_witt_index(q7) % 4 == 2
    q7pad = gw_direct_sum(q7, hyperbolic_form(K7, 4))
    @test any(hasse_witt_invariant(q7, P) != hasse_witt_invariant(q7pad, P)
              for (P, _) in prime_decomposition(OK7, 2))
    qa7 = anisotropic_part(q7)
    @test is_isomorphic_form(q7, gw_direct_sum(qa7, hyperbolic_form(K7, 4)))
end

@testset "anisotropic_part over K: the full specification, and idempotence" begin
    Keps = 35 + 6 * s34
    battery = GWClass[
        gw_direct_sum(diagonal_form(K2, (1, s2)), hyperbolic_form(K2, 2)),
        diagonal_form(K2, (1, 1, 1)),
        diagonal_form(K2, (s2, -s2)),                 # hyperbolic: empty part
        diagonal_form(K2, (1, 1, 1, 1, 1, 1)),        # chain from dim_a 6
        diagonal_form(K7, (1, 1, 1)),                 # Alg 6 then 8
        diagonal_form(K7, (1, -3, -s7, 3 * s7)),      # Alg 5 -> 6 -> 8
        hyperbolic_form(K7, 6),
        diagonal_form(K7, (1, -1, 3)),                # dim_a 1: Observation 5.3
        diagonal_form(K7, 3),                         # rank 1
        diagonal_form(K5, (1, 2, s5)),                # h = 2, Alg 6 then 8
        diagonal_form(K34, (1, -Keps)),                # blind spot, Alg 8 direct
        diagonal_form(Kc, (-1, -1, -1)),               # Alg 7 then 8, cubic
        # the rest of the field matrix:
        diagonal_form(Kr5, (1, 1, 1)),                 # inert dyadic, Alg 7 -> 8
        diagonal_form(K21, (1, 3, s21)),               # C = Z/2 x Z/2: lambda_b in Alg 8
        diagonal_form(K32, (1, 1, 1)),                 # mixed signature, Alg 7 -> 8
        diagonal_form(Kq, (1, -aq, 2)),                # totally real quartic
        diagonal_form(Ks10, (-1, -1, -1)),             # ramified dyadic + h = 2, Alg 7
        diagonal_form(Ks10, (s10, 2, -1)),             # mixed signature over Q(sqrt10)
        diagonal_form(K162, (-1, -1, -1)),             # WILD dyadic (e = 4), Alg 7
        diagonal_form(K162, (w16, 1, 1)),              # mixed signs, four real places
        diagonal_form(K2i, (1, 1, 1)),                 # ramified dyadic, non-real
        diagonal_form(K2i, (1, s2i)),                  # rank 2 over the control field
        diagonal_form(K34, (1, -s34)),                 # lambda_b WITH a real-place row
    ]
    for q in battery
        K = gw_base_field(q)
        qa = anisotropic_part(q)
        # the full specification — everything else is a probe:
        @test is_anisotropic_form(qa)
        @test form_rank(qa) == anisotropic_dimension(q)
        w = gw_witt_index(q)
        @test is_isomorphic_form(q, gw_direct_sum(qa, hyperbolic_form(K, 2w)))
    end
    # idempotence on a cross-section (definite, Alg-8 output, Obs 5.3 output)
    for q in (battery[2], battery[11], battery[8])
        qa = anisotropic_part(q)
        @test is_isomorphic_form(anisotropic_part(qa), qa)
    end
end

@testset "Algorithm 8 details: A-block, shape relation, batching" begin
    # a formally real field with sigma(d) < 0 somewhere: the A block is
    # non-empty (over the imaginary quadratic fields above it is empty and
    # that half of the solver is otherwise untested)
    q = diagonal_form(K2, (1, s2))
    setup = MH._kr_algorithm_8_setup(q)
    @test length(setup.negd) == 1                      # d = -sqrt2 < 0 at one embedding
    qa = MH._kr_algorithm_8(q)
    @test is_isomorphic_form(qa, q)                    # q was already anisotropic

    # shape relation on trivial-class-group fields: dim E_S = r1 + r2 + |S|,
    # preserved as our BATCHED enlargement appends primes (KR append their
    # own three in §6; ours differ, the relation must not)
    for (K, form) in ((K7, diagonal_form(K7, (1, s7))),
                      (K2, diagonal_form(K2, (1, s2))))
        S = copy(MH._kr_form(form).S)
        for _ in 1:2
            @test length(MH._kr_algorithm_1(K, S)) == sum(signature(K)) + length(S)
            batch = MH._next_primes_batch(K, S, 3)
            @test length(batch) == 3
            @test all(!any(Q -> Q == P, S) for P in batch)
            append!(S, batch)
        end
    end
end

@testset "sum_decomposition over K (M2-mirrored cache overwrite)" begin
    q6 = diagonal_form(K7,
        (-3 - 9s7, -1, -2 - 6s7, 1 - s7, -6 + 4s7, -3 + 2s7, 4 - 4s7))
    kr_before = MH._kr_form(q6)
    r = sum_decomposition(q6)
    @test is_isomorphic_form(r, q6)
    @test form_rank(r) == 7
    s = sum_decomposition_string(q6)
    @test startswith(s, "2H + <")
    # per the recorded decision: the number-field path mirrors M2/QQ and
    # OVERWRITES the cached diagonal class with the decomposition rep ...
    @test diagonal_class(q6) === r
    # ... and _KRForm is immune: it lives under its own cache key and reads
    # the Gram matrix directly, never diagonal_class
    @test MH._kr_form(q6) === kr_before
end

@testset "dyadic stress: ramified, inert, split" begin
    # Forms whose isotropy pattern is decided at a dyadic place, and forms
    # whose coefficients all have even dyadic valuation — over the ramified
    # (Q(sqrt2)), inert (Q(sqrt5)), and split (Q(sqrt-7)) fields. Premises
    # are asserted before the property.
    OK2 = maximal_order(K2)
    P2r = prime_decomposition(OK2, 2)[1][1]
    OKr5 = maximal_order(Kr5)
    P2i = prime_decomposition(OKr5, 2)[1][1]
    dys7 = MH._dyadic_primes(K7)

    # inert: <4, 4√5, 12> is isotropic AT the dyadic place (local dim 1)
    # while anisotropic of dimension 3 at an odd prime and at a real place —
    # and every coefficient has even dyadic valuation
    q = diagonal_form(Kr5, (4, 4 * sr5, 12))
    kr = MH._kr_form(q)
    @test all(iseven(valuation(a, P2i)) for a in kr.diag)
    @test MH._local_anisotropic_dimension(q, P2i) == 1
    @test maximum([MH._local_anisotropic_dimension(q, P) for P in kr.odd_primes]) == 3
    @test anisotropic_dimension(q) == 3
    qa = anisotropic_part(q)
    @test is_anisotropic_form(qa) && form_rank(qa) == 3
    @test is_isomorphic_form(q, qa)

    # ramified: <1,1,1,1> is isotropic at the dyadic place yet globally
    # anisotropic (definite at both real places)
    q = diagonal_form(K2, (1, 1, 1, 1))
    @test MH._local_anisotropic_dimension(q, P2r) == 0
    @test anisotropic_dimension(q) == 4
    @test is_anisotropic_form(q)

    # split: <1,1,1> over Q(sqrt-7) is anisotropic ONLY at the dyadic places
    # (P(q) is empty, there are no real places) — the dyadic computation
    # alone carries the answer
    q = diagonal_form(K7, (1, 1, 1))
    @test isempty(MH._kr_form(q).odd_primes)
    @test [MH._local_anisotropic_dimension(q, P) for P in dys7] == [3, 3]
    @test anisotropic_dimension(q) == 3
    qa = anisotropic_part(q)
    @test is_anisotropic_form(qa) && form_rank(qa) == 3
    @test is_isomorphic_form(q, qa)

    # even dyadic valuation throughout with the answer decided dyadically,
    # in all three splitting behaviours
    for (K, entries, dys) in ((K2, (4, -4 * (1 + 2 * s2)), [P2r]),
                              (Kr5, (4, -4 * (2 + sr5)), [P2i]),
                              (K7, (4, -4 * (2 + s7)), dys7))
        q = diagonal_form(K, entries)
        kr = MH._kr_form(q)
        for P in dys
            @test all(iseven(valuation(a, P)) for a in kr.diag)
            @test MH._local_anisotropic_dimension(q, P) == 2
        end
        @test anisotropic_dimension(q) == 2
        qa = anisotropic_part(q)
        @test is_anisotropic_form(qa) && form_rank(qa) == 2
        @test is_isomorphic_form(q, qa)
    end
end

@testset "scale probe: dimension 10 over a degree-10 totally real field (KR §7)" begin
    # KR's Magma timings for their dimension-10 forms over a degree-10
    # totally real field: 0.24-3.42 s, mean 2.82 s (their §7). Correctness
    # is asserted; the wall clock is reported via @info, not asserted.
    Zx, zx = polynomial_ring(ZZ, "zx")
    K10, c10 = number_field(map_coefficients(QQ, cos_minpoly(25, zx)), "c10")
    @test degree(K10) == 10
    @test signature(K10) == (10, 0)
    q10 = diagonal_form(K10,
        (1, c10, c10^2 - 2, 3, -c10 - 1, 2c10, c10^3 - 3c10, 5, c10 - 2, -7))
    t10 = @elapsed qa10 = anisotropic_part(q10)
    @info "dimension-10 / degree-10 totally real wall-clock" seconds = round(t10; digits = 2)
    @test is_anisotropic_form(qa10)
    @test form_rank(qa10) == anisotropic_dimension(q10)
    @test is_isomorphic_form(q10,
        gw_direct_sum(qa10, hyperbolic_form(K10, 10 - form_rank(qa10))))
end

@testset "boundary size reduction over K (square-class preserving, NOT canonical)" begin
    # the exported boundary output vs the raw KR driver output: same class,
    # entries reduced — across the whole field matrix
    for (K, s) in ((K7, s7), (K2, s2), (Kr5, sr5), (K5, s5), (K21, s21),
                   (K34, s34), (Kc, ac), (K32, a32), (Kq, aq),
                   (Ks10, s10), (K162, w16), (K2i, s2i))
        for q in (gw_direct_sum(diagonal_form(K, (3, -s)), hyperbolic_form(K)),
                  diagonal_form(K, (1, 2, s)))
            raw = MH._anisotropic_part_nf(
                GWClass(MH._nondegenerate_part_diagonal(gw_matrix(q))))
            qa = anisotropic_part(q)
            @test is_isomorphic_form(qa, raw)
            @test is_isomorphic_form(q,
                gw_direct_sum(qa, hyperbolic_form(K, 2 * gw_witt_index(q))))
            for a in diagonal_entries(qa)
                @test is_square(a * MH._size_reduced_square_class_rep(a))
            end
        end
    end

    # the helper itself, on a hostile sweep: a nonsquare unit, a perfect
    # square, a negative, non-integral elements, big products, and the wild
    # dyadic field — class preserved and printed size never worse
    hostile = [
        (K34, 35 + 6 * s34), (K34, (6 + s34)^2), (K34, -(35 + 6 * s34)),
        (K5, (2 + s5)^4 * 3), (K5, s5 // 4),
        (Ks10, (3 + s10)^6 * (1 + s10)),
        (K162, (1 + w16)^7), (K2, -(1 + s2)^3),
    ]
    for (K, a) in hostile
        b = MH._size_reduced_square_class_rep(a)
        @test is_square(a * b)
        @test length(string(b)) <= length(string(a))
    end
    # a nonsquare unit must stay a nonsquare — reduction may NOT collapse
    # unit classes to 1
    @test !is_square(MH._size_reduced_square_class_rep(35 + 6 * s34))
end

@testset "regression: the audit's Q(sqrt-5) rank-2 boundary blow-up" begin
    # The API coverage audit observed a 20-digit coefficient in
    # sum_decomposition_string(<3, -sqrt-5> + H) and 29-character
    # anisotropic_part entries over Q(sqrt-5). The boundary reduction exists
    # because of this input, so it is pinned here by name. Size improvement
    # is reported via @info, deliberately not asserted (representative size
    # is not a contract; the square class is).
    q = gw_direct_sum(diagonal_form(K5, (3, -s5)), hyperbolic_form(K5))
    raw = MH._anisotropic_part_nf(
        GWClass(MH._nondegenerate_part_diagonal(gw_matrix(q))))
    qa = anisotropic_part(q)
    @test is_isomorphic_form(qa, raw)
    before = maximum([length(string(a)) for a in diagonal_entries(raw)]; init = 0)
    after = maximum([length(string(a)) for a in diagonal_entries(qa)]; init = 0)
    @info "Q(sqrt-5) audit regression: max anisotropic_part entry (chars)" before after
    str = sum_decomposition_string(q)
    @info "Q(sqrt-5) audit regression: sum_decomposition_string" str
    @test occursin("H", str)
    @test is_isomorphic_form(sum_decomposition(q), q)

    # Q(sqrt10), h = 2: where the representative question bites hardest —
    # exponent reduction can strand a non-principal ideal class
    for q10 in (gw_direct_sum(diagonal_form(Ks10, (3, -s10)), hyperbolic_form(Ks10)),
                diagonal_form(Ks10, ((3 + s10)^2 * s10, -1)),
                diagonal_form(Ks10, (1 + s10, 2, -3)))
        raw10 = MH._anisotropic_part_nf(
            GWClass(MH._nondegenerate_part_diagonal(gw_matrix(q10))))
        qa10 = anisotropic_part(q10)
        @test is_isomorphic_form(qa10, raw10)
        b = maximum([length(string(a)) for a in diagonal_entries(raw10)]; init = 0)
        a_ = maximum([length(string(a)) for a in diagonal_entries(qa10)]; init = 0)
        @info "Q(sqrt10) boundary sizes: max entry (chars)" before = b after = a_
    end
end

@testset "reduction and the diagonal_class cache (M2-mirrored overwrite intact)" begin
    # both writers of the :diagonal_class slot now produce boundary-reduced
    # entries, so the cache can never hold a mix — checked in both orders
    mk() = gw_direct_sum(diagonal_form(K5, (3, -s5)), hyperbolic_form(K5))
    q1 = mk()
    d1 = diagonal_class(q1)
    r1 = sum_decomposition(q1)
    @test diagonal_class(q1) === r1          # the M2-mirrored overwrite, kept
    @test is_isomorphic_form(r1, d1)
    q2 = mk()
    r2 = sum_decomposition(q2)
    @test diagonal_class(q2) === r2
    # whichever writer ran, every cached entry represents its own reduced
    # square class
    for r in (r1, r2), a in diagonal_entries(r)
        @test is_square(a * MH._size_reduced_square_class_rep(a))
    end
end

@testset "simplification surface across the field matrix (rank 0, 1, and beyond)" begin
    for (K, s) in ((K7, s7), (K2, s2), (Kr5, sr5), (K5, s5), (K21, s21),
                   (K34, s34), (Kc, ac), (K32, a32), (Kq, aq),
                   (Ks10, s10), (K162, w16), (K2i, s2i), (KQ, KQ(5)))
        # rank 0
        E = diagonal_form(K, ())
        @test nrows(diagonalize_via_congruence(gw_matrix(E))) == 0
        @test form_rank(diagonal_class(E)) == 0
        @test sum_decomposition(E) == E
        @test sum_decomposition_string(E) == "empty form"
        # rank 1
        b1 = diagonal_form(K, (2 + s,))
        @test diagonalize_via_congruence(gw_matrix(b1)) == gw_matrix(b1)
        d1 = diagonal_class(b1)
        @test form_rank(d1) == 1
        @test is_isomorphic_form(d1, b1)
        r1 = sum_decomposition(b1)
        @test form_rank(r1) == 1
        @test is_isomorphic_form(r1, b1)
        s1 = sum_decomposition_string(b1)
        @test startswith(s1, "<") && !occursin("H", s1)
        # rank 2, nondiagonal input
        M = matrix(K, 2, 2, [K(2), K(1), K(1), 2 + s])
        D = diagonalize_via_congruence(M)
        @test is_diagonal(D)
        @test is_isomorphic_form(GWClass(D), GWClass(M))
        beta = GWClass(M)
        dc = diagonal_class(beta)
        @test is_diagonal(gw_matrix(dc))
        @test is_isomorphic_form(dc, beta)
        @test diagonal_class(beta) === dc    # cached: repeated calls are free
        # a hyperbolic summand shows up in the string, and the decomposition
        # round-trips
        q = gw_direct_sum(diagonal_form(K, (2 + s,)), hyperbolic_form(K))
        @test gw_witt_index(q) == 1
        @test occursin("H", sum_decomposition_string(q))
        @test is_isomorphic_form(sum_decomposition(q), q)
    end
end

@testset "Algorithm 8 closing-prime enlargement (Q(sqrt5) large-prime regression)" begin
    # <1, -u, v, -u*v> with u = 1000003, v = 1 - 2 sqrt5: K(sqrt d)/K is
    # ramified at eleven primes, so the blind step-(5g) enlargement can
    # never close the F2 system inside the retry bound; the closing-prime
    # search must. The form has anisotropic dimension 4 with signatures
    # [0, 0], so the whole reduction chain (Algorithms 5, 7, 8) runs.
    u = Kr5(1000003)
    v = 1 - 2 * sr5
    beta = diagonal_form(Kr5, (Kr5(1), -u, v, -u * v))
    @test anisotropic_dimension(beta) == 4
    @test form_signatures(beta) == [0, 0]
    part = anisotropic_part(beta)
    @test form_rank(part) == 4
    @test gw_witt_index(beta) == 0
    @test is_isomorphic_form(part, beta)
    @test is_anisotropic_form(part)
end
