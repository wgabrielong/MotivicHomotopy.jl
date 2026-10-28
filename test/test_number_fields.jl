# Number-field groundwork: the FieldKind trait, the arithmetic context,
# GW classes over K, per-place signatures, and the square-class
# representative normalisation. No KR21 mathematics yet — those algorithms
# arrive with later work; here K-classes construct, print, diagonalize, and
# the invariant layers reject them exactly as before.

using Test
using Oscar
using MotivicHomotopy

const MH = MotivicHomotopy

# The shared number-field test matrix (see nf_test_fields.jl).
include("nf_test_fields.jl")

@testset "field matrix properties (asserted, not trusted)" begin
    # (field, (r1, r2), dyadic (e, f) pairs, class number)
    for (K, sig, dyadic, h) in (
            (K7, (0, 1), [(1, 1), (1, 1)], 1),        # 2 split
            (K2, (2, 0), [(2, 1)], 1),                # 2 ramified
            (Kr5, (2, 0), [(1, 2)], 1),               # 2 inert
            (K5, (0, 1), [(2, 1)], 2),
            (K21, (0, 1), [(2, 1)], 4),               # Cl = Z/2 x Z/2
            (K34, (2, 0), [(2, 1)], 2),
            (Kc, (3, 0), [(1, 3)], 1),
            (K32, (1, 1), [(3, 1)], 1),               # 2 totally ramified, e = 3
            (Kq, (4, 0), [(4, 1)], 1),        # 2 is totally ramified here TOO
                                              # (found by this assertion — the
                                              # suite had never pinned Kq's
                                              # dyadic structure)
            (Ks10, (2, 0), [(2, 1)], 2),
            (K162, (4, 0), [(4, 1)], 1),              # 2 TOTALLY ramified, e = 4
            (K2i, (0, 1), [(2, 1)], 1))
        @test signature(K) == sig
        OK = maximal_order(K)
        dec = prime_decomposition(OK, 2)
        @test sort([(e, degree(P)) for (P, e) in dec]) == sort(dyadic)
        C, mC = class_group(OK)
        @test order(C) == h
    end
    # K21's class group is elementary abelian of rank 2, not Z/4
    C21, _ = class_group(maximal_order(K21))
    @test elementary_divisors(C21) == [2, 2]
    # K162 really is the maximal real subfield of Q(zeta_16): its defining
    # polynomial is the minimal polynomial of 2 cos(2 pi / 16)
    let
        Zx, zx = polynomial_ring(ZZ, "zx")
        g = cos_minpoly(16, zx)
        f = defining_polynomial(K162)
        @test degree(f) == degree(g) &&
              all(coeff(f, i) == coeff(g, i) for i in 0:degree(g))
    end
    # the corrected Algorithm 4 modulus over K162: ord_P(4) = 2 e(P|2) = 8,
    # so k_P = 1 + ord_P(4) + ord_P(alpha_1) >= 9 — the "wild" claim is real
    P16 = MH._dyadic_primes(K162)[1]
    @test valuation(K162(4), P16) == 8

    # FINDING (verified, not assumed): Q(sqrt10) can NEVER exercise the
    # lambda_b branch of Algorithm 1. Its dyadic prime is non-principal and
    # its class generates Cl = Z/2, so C_S = Cl / <[P] : P in S> is trivial
    # for every admissible S (S always contains the dyadic prime) — the same
    # mechanism as the recorded K5 divergence. The field of the matrix where
    # lambda_b coexists with real places is K34, whose dyadic prime is
    # PRINCIPAL (6 + sqrt34 has norm 2).
    for (K, dyadic_principal, lambda_possible) in
            ((Ks10, false, false), (K34, true, true), (K5, false, false))
        OK = maximal_order(K)
        C, mC = class_group(OK)
        dys = MH._dyadic_primes(K)
        @test all(is_principal(P) for P in dys) == dyadic_principal
        CS, _ = quo(C, [preimage(mC, P) for P in dys])
        CSs, _ = snf(CS)
        has2 = any(iseven(order(CSs[i])) for i in 1:ngens(CSs))
        @test has2 == lambda_possible
    end
    @test (6 + s34) * maximal_order(K34) == MH._dyadic_primes(K34)[1]
end

@testset "FieldKind trait and real-ness predicates" begin
    @test MH.field_kind(QQ) isa MH.RationalKind
    @test MH.field_kind(GF(5)) isa MH.FiniteKind
    @test MH.field_kind(GF(2)) isa MH.UnsupportedKind
    @test MH.field_kind(K7) isa MH.NumberFieldKind
    @test MH.field_kind(Kc) isa MH.NumberFieldKind
    @test MH.field_kind("not a field") isa MH.UnsupportedKind
    # deliberately no Type-level arms for the float backends
    @test MH.field_kind(Float64) isa MH.UnsupportedKind
    @test MH.field_kind(ComplexF64) isa MH.UnsupportedKind

    @test MH._exact_supported(QQ)
    @test MH._exact_supported(GF(7))
    @test !MH._exact_supported(GF(2))
    @test !MH._exact_supported(K7)

    @test MH._is_formally_real(QQ)
    @test !MH._is_formally_real(GF(3))
    @test !MH._is_formally_real(K7)      # imaginary quadratic
    @test MH._is_formally_real(K2)
    @test !MH._is_formally_real(K5)
    @test MH._is_formally_real(Kc)

    @test MH._is_totally_real(QQ)
    @test !MH._is_totally_real(GF(3))
    @test MH._is_totally_real(K2)
    @test MH._is_totally_real(Kc)
    @test !MH._is_totally_real(K7)
    @test !MH._is_totally_real(K5)
end

@testset "Real embeddings with isolating intervals, and the KR etas" begin
    # non-real field: no pairs
    @test isempty(MH._real_embeddings_with_intervals(K7))

    for K in (K2, Kc)
        pairs = MH._real_embeddings_with_intervals(K)
        @test length(pairs) == signature(K)[1]
        @test first.(pairs) == real_embeddings(K)   # coupled in embedding order
        theta = gen(K)
        for (i, (sigma, (a, b))) in enumerate(pairs)
            @test a < b
            # the interval isolates THIS embedding's image of theta ...
            @test sign(theta - a, sigma) == 1
            @test sign(theta - b, sigma) == -1
            # ... and excludes every other embedding's image
            for (j, (tau, _)) in enumerate(pairs)
                j == i && continue
                @test sign(theta - a, tau) == sign(theta - b, tau)
            end
        end

        # the property the etas exist for: eta_i negative exactly at sigma_i
        etas = MH._kr_etas(K, pairs)
        @test length(etas) == length(pairs)
        for (i, eta) in enumerate(etas)
            for (j, (tau, _)) in enumerate(pairs)
                @test sign(eta, tau) == (i == j ? -1 : 1)
            end
        end
    end
end

@testset "Non-simple number fields get a targeted error" begin
    Kns, _ = number_field([x^2 - 2, x^2 - 3], "g")
    beta = GWClass(matrix(Kns, 2, 2, [Kns(1), Kns(0), Kns(0), Kns(2)]))
    @test_throws "simple_extension" anisotropic_dimension(beta)
    @test_throws "simple_extension" anisotropic_part(beta)
    @test_throws "simple_extension" is_isotropic_form(beta)
    @test_throws "simple_extension" form_signatures(beta)
    @test_throws "simple_extension" is_isomorphic_form(beta, beta)
end

@testset "GW classes over K: construction, accessors, printing" begin
    for (K, s) in ((K7, s7), (K2, s2), (K5, s5), (Kc, ac))
        # from a Gram matrix over K
        M = matrix(K, 2, 2, [K(2), s, s, K(3)])
        beta = GWClass(M)
        @test gw_matrix(beta) == M
        @test gw_algebra(beta) === K
        @test gw_base_field(beta) === K
        @test form_rank(beta) == 2

        # from a diagonal list (K elements, integers, rationals mixed)
        delta = diagonal_form(K, (1, s, -2 // 3))
        @test gw_base_field(delta) === K
        @test form_rank(delta) == 3
        @test diagonal_entries(delta) == [K(1), s, K(-2 // 3)]

        # from an integer
        gamma = diagonal_form(K, 5)
        @test gw_matrix(gamma) == matrix(K, 1, 1, [K(5)])

        # hyperbolic constructor
        @test diagonal_entries(hyperbolic_form(K, 4)) ==
              [K(1), K(-1), K(1), K(-1)]

        # elements of a different number field are rejected
        other = K === K7 ? s2 : s7
        @test_throws ErrorException diagonal_form(K, (other,))

        # well-definedness over K: symmetry and non-degeneracy
        @test_throws ErrorException GWClass(matrix(K, 2, 2, [K(1), s, K(0), K(1)]))
        @test_throws ErrorException GWClass(matrix(K, 2, 2, [K(1), s, s, s^2]))

        # equality is literal Gram equality; printing appends the field
        @test beta == GWClass(M)
        @test beta != delta
        str = sprint(show, MIME("text/plain"), beta)
        @test occursin("\nover ", str)
        @test occursin(string(K), str)
        @test occursin("over " * string(K), sprint(show, beta))
    end

    # printing over the existing backends is unchanged: matrix only
    q = GWClass(QQ[2 1; 1 3])
    @test sprint(show, MIME("text/plain"), q) ==
          sprint(show, MIME("text/plain"), gw_matrix(q))
    @test sprint(show, q) == sprint(show, gw_matrix(q))
    f = GWClass(GF(7)[1 0; 0 3])
    @test sprint(show, MIME("text/plain"), f) ==
          sprint(show, MIME("text/plain"), gw_matrix(f))
    r = GWClass([2.0 0.0; 0.0 -1.0])
    @test sprint(show, MIME("text/plain"), r) ==
          sprint(show, MIME("text/plain"), gw_matrix(r))
end

@testset "Diagonalisation over K" begin
    for (K, s) in ((K7, s7), (K2, s2), (K5, s5), (Kc, ac))
        # zero-pivot branch: antidiagonal Gram matrix
        M0 = matrix(K, 2, 2, [K(0), K(1), K(1), K(0)])
        D0 = diagonalize_via_congruence(M0)
        @test is_diagonal(D0)
        @test all(!iszero(D0[i, i]) for i in 1:2)
        @test is_square(det(M0) * det(D0))       # congruent: det differs by a square

        # dense matrix with generator entries
        M = matrix(K, 3, 3, [K(2), s, K(1), s, K(3), K(0), K(1), K(0), s + 2])
        @test !iszero(det(M))
        D = diagonalize_via_congruence(M)
        @test is_diagonal(D)
        @test form_rank(D) == 3
        @test is_square(det(M) * det(D))

        # diagonal_class over K SIZE-reduces each entry within its square
        # class (not a canonical form — none exists over K), so the diagonal
        # representative is pinned up to isometry plus the entrywise
        # square-class relation to the raw diagonalization, never literally.
        beta = GWClass(M)
        dc = diagonal_class(beta)
        @test is_diagonal(gw_matrix(dc))
        @test is_isomorphic_form(dc, beta)
        @test all(is_square(gw_matrix(dc)[i, i] * D[i, i]) for i in 1:nrows(D))

        # is_isomorphic_form supports K since the invariants layer landed
        @test is_isomorphic_form(GWClass(M), GWClass(M))
        @test is_isomorphic_form(GWClass(M), GWClass(D))   # congruent to its diagonalisation
    end
end

@testset "form_signature / form_signatures per real place" begin
    # totally real cubic: independent check by numeric evaluation of the
    # diagonal entries at the (hardcoded) real roots of x^3 - 3x - 1
    roots = [-1.5320888862379560, -0.3472963553338607, 1.8793852415718168]
    beta = diagonal_form(Kc, (1, ac, ac - 1, ac + 2))
    D = diagonal_entries(beta)
    approx(d, r) = sum(Float64(coeff(d, i)) * r^i for i in 0:2)
    expected = [sum(sign(approx(d, r)) for d in D) for r in roots]
    sigs = form_signatures(beta)
    @test length(sigs) == 3
    @test sort(sigs) == sort(Int.(expected))
    @test sigs == [form_signature(beta, P) for P in real_places(Kc)]

    # real quadratic: hand-computable signs. sqrt(2) is positive under one
    # embedding and negative under the other.
    gamma = diagonal_form(K2, (1, s2))
    sig2 = form_signatures(gamma)
    @test sort(sig2) == [0, 2]
    @test sig2 == [form_signature(gamma, P) for P in real_places(K2)]
    # embeddings are accepted as well as places
    @test sort([form_signature(gamma, e) for e in real_embeddings(K2)]) == [0, 2]

    # non-real field: no real places, empty vector
    @test form_signatures(diagonal_form(K7, (1, s7))) == Int[]

    # guardrails
    @test_throws ErrorException form_signatures(GWClass(QQ[1 0; 0 1]))
    @test_throws ErrorException form_signature(GWClass(QQ[1 0; 0 1]),
                                               real_places(K2)[1])
    @test_throws ErrorException form_signature(gamma, real_places(Kc)[1])
    # the 1-arg method stays QQ/RR-only
    @test_throws ErrorException form_signature(beta)
end

@testset "_squarefree_content_rep" begin
    OK5 = maximal_order(K5)
    OK7 = maximal_order(K7)

    # class number 2: (2) = P^2 with P non-principal, so <2> cannot be reduced
    b, fl = MH._squarefree_content_rep(K5(2))
    @test !fl
    @test b == K5(2)

    # (4) = P^4 and P^2 = (2) is principal, so <4> reduces to a unit
    b, fl = MH._squarefree_content_rep(K5(4))
    @test fl
    @test is_one(b * OK5)                        # content is the unit ideal
    @test is_square(b * K5(4))                   # same square class

    # class number 1: (s7^3) = p^3 with p = (s7); reduces to content p
    b, fl = MH._squarefree_content_rep(s7^3)
    @test fl
    @test isone(denominator(b, OK7))             # integral
    @test is_square(b * s7^3)                    # same square class
    @test all(e <= 1 for (_, e) in factor(b * OK7))

    # denominators are cleared by squares
    b, fl = MH._squarefree_content_rep(K7(9 // 4))
    @test fl
    @test isone(denominator(b, OK7))
    @test is_square(b * K7(9 // 4))
    @test is_one(b * OK7) || all(e <= 1 for (_, e) in factor(b * OK7))

    @test_throws ErrorException MH._squarefree_content_rep(K7(0))
end

# --- KR21 Algorithms 1-4 --------------------------------------------------------------

# Support primes of a FacElem lying outside S, deduplicated and compared by ==
# (ideal hashing is not relied on).
function _test_support_outside(fe, OK, S)
    out = []
    for (b, _) in fe.fac
        J = b * OK
        is_one(J) && continue
        for (P, _) in factor(J)
            any(Q -> Q == P, S) && continue
            any(Q -> Q == P, out) || push!(out, P)
        end
    end
    out
end

@testset "KR Algorithm 2: all sign patterns, complete alignment" begin
    @test signature(Kq) == (4, 0)
    for K in (Kc, Kq)
        pairs = MH._real_embeddings_with_intervals(K)
        etas = MH._kr_etas(K, pairs)
        r = length(pairs)
        @test r == degree(K)
        # complete alignment check: eta_i negative at sigma_i, positive elsewhere
        for i in 1:r, j in 1:r
            @test sign(etas[i], pairs[j][1]) == (i == j ? -1 : 1)
        end
        # every one of the 2^r sign patterns is realised
        for mask in 0:(2^r - 1)
            I = [i for i in 1:r if (mask >> (i - 1)) & 1 == 1]
            rho = MH._kr_algorithm_2(K, pairs, etas, I)
            for (i, (sigma, _)) in enumerate(pairs)
                @test sign(rho, sigma) == (i in I ? -1 : 1)
            end
        end
    end

    # overlapping intervals are rejected before any product is formed
    embs2 = real_embeddings(K2)
    badpairs = [(embs2[1], (QQ(-2), QQ(1))), (embs2[2], (QQ(0), QQ(2)))]
    @test_throws "overlap" MH._kr_algorithm_2(K2, badpairs, [K2(1), K2(1)], Int[])
end

@testset "KR Algorithm 3: totally positive with congruences" begin
    for (K, s) in ((K2, s2), (Kr5, sr5), (K7, s7))
        OK = maximal_order(K)
        dy = [P for (P, _) in prime_decomposition(OK, 2)]
        podd = prime_decomposition(OK, 11)[1][1]
        primes = vcat(dy, [podd])
        ks = vcat([3 + valuation(K(4), P) for P in dy], [2])
        lambdas = [K(i + 1) + s for i in eachindex(primes)]
        alpha = MH._kr_algorithm_3(K, primes, ks, lambdas)
        @test is_totally_positive(alpha)
        for (P, k, lam) in zip(primes, ks, lambdas)
            d = alpha - lam
            @test iszero(d) || valuation(d, P) >= k
        end
    end

    # non-integral targets are refused, not silently accepted
    P = prime_decomposition(maximal_order(K2), 3)[1][1]
    @test_throws "integers" MH._kr_algorithm_3(K2, [P], [1], [K2(1 // 3)])
end

@testset "KR Algorithm 4: signs + local squares; the modulus correction is load-bearing" begin
    OK2 = maximal_order(K2)
    P2 = prime_decomposition(OK2, 2)[1][1]       # ramified dyadic: e = 2
    @test valuation(K2(4), P2) == 4
    embs = real_embeddings(K2)
    eneg = embs[findfirst(e -> is_negative(s2, e), embs)]
    epos = embs[findfirst(e -> is_positive(s2, e), embs)]

    # Hand-built (valid) isolating intervals whose eta has POSITIVE
    # P2-valuation: eta_neg = (theta + 2)(theta - 0) = 2 + 2 theta, and
    # 1 + theta is a unit (norm -1), so ord_P2(eta_neg) = ord_P2(2) = 2.
    pairs = [(eneg, (QQ(-2), QQ(0))), (epos, (QQ(0), QQ(2)))]
    etas = MH._kr_etas(K2, pairs)
    @test etas == [2 + 2 * s2, 2 - 2 * s2]

    alpha1 = MH._kr_algorithm_2(K2, pairs, etas, [1])
    @test alpha1 == 2 + 2 * s2
    @test isone(denominator(alpha1, OK2))        # integralisation is a no-op here
    @test valuation(alpha1, P2) == 2             # the premise: v_P > 0 after integralisation

    S = [P2]
    rho = MH._kr_algorithm_4(K2, pairs, etas, [1], S)
    @test is_local_square(rho, P2)
    @test sign(rho, eneg) == -1
    @test sign(rho, epos) == 1

    # The correction is load-bearing on this input: KR's modulus
    # P^(1 + ord_P 4) = P^5 admits a COMPLIANT alpha_2 whose product with
    # alpha_1 is not a local square. With alpha2bad = alpha1 + theta^5 the
    # difference has valuation exactly 5, so u = alpha2bad/alpha1 has
    # ord(u - 1) = 5 - 2 = 3; in Q_2(sqrt 2) every square unit w^2 = (1+z)^2
    # has ord(w^2 - 1) = ord(2z + z^2) in {0, 2} or >= 4 — never 3 — so u is
    # not a square and neither is rho_bad. (Re-running the full pipeline with
    # the uncorrected exponent can pass by accident: CRT and the +t*s step are
    # free to overshoot the modulus. The witness form is deterministic.)
    alpha2bad = alpha1 + s2^5
    @test valuation(alpha2bad - alpha1, P2) == 5 == 1 + valuation(K2(4), P2)
    @test !is_local_square(alpha1 * alpha2bad, P2)

    # Negative control: alpha_2 alone carries the congruences but the WRONG
    # signs whenever the requested pattern is non-trivial (it is totally
    # positive by construction; the pattern asks for -1 at eneg).
    v = valuation(alpha1, P2)
    ks = [1 + valuation(K2(4), P2) + v]
    alpha2 = MH._kr_algorithm_3(K2, S, ks, [alpha1])
    @test is_totally_positive(alpha2)
    @test sign(alpha2, eneg) == 1
    d = alpha2 - alpha1
    @test iszero(d) || valuation(d, P2) >= ks[1]
end

@testset "KR Algorithm 1 (Czogala): F2-basis of E_S" begin
    # Q(sqrt-7), h = 1, with the nine finite primes of the [KR21 S6] example:
    # both dyadic primes, the inert (3), (19), (5), the ramified (sqrt-7),
    # one prime over 11, and both primes over the split 37.
    OK7 = maximal_order(K7)
    dy7 = [P for (P, _) in prime_decomposition(OK7, 2)]
    @test length(dy7) == 2
    split37 = [P for (P, _) in prime_decomposition(OK7, 37)]
    @test length(split37) == 2
    S9 = vcat(dy7, [prime_decomposition(OK7, p)[1][1] for p in (3, 19, 5, 7, 11)], split37)
    @test length(S9) == 9
    basis7 = MH._kr_algorithm_1(K7, S9)
    # trivial class group: dim = torsion (1) + S-unit rank (0 + 1 - 1 + 9) = 10
    @test length(basis7) == 10
    # S-singularity: even valuation at every support prime outside S
    for fe in basis7
        for P in _test_support_outside(fe, OK7, S9)
            @test iseven(valuation(fe, P))
        end
    end

    # Q(sqrt-5), h = 2: the class group is GENERATED by the class of the
    # unique dyadic prime P2 (Minkowski bound < 3), and Algorithm 1 requires
    # P2 in S, so C_S is trivial for every admissible S and the lambda branch
    # cannot fire — dim E_S equals the unit part 1 + |S|, NOT more.
    # [Divergence from the review prompt's expectation; reported.]
    OK5 = maximal_order(K5)
    P2_5 = prime_decomposition(OK5, 2)[1][1]
    @test !is_principal(P2_5)
    @test length(MH._kr_algorithm_1(K5, [P2_5])) == 2
    P3_5 = prime_decomposition(OK5, 3)[1][1]
    @test length(MH._kr_algorithm_1(K5, [P2_5, P3_5])) == 3

    # Q(sqrt-21), C = Z/2 x Z/2: modding out [P2] leaves Z/2, so exactly one
    # lambda_b fires and the basis is strictly larger than the unit part.
    OK21 = maximal_order(K21)
    P2_21 = prime_decomposition(OK21, 2)[1][1]
    U21, _ = sunit_group_fac_elem([P2_21])
    basis21 = MH._kr_algorithm_1(K21, [P2_21])
    @test length(basis21) == ngens(U21) + 1
    @test length(basis21) == 3
    for fe in basis21
        for P in _test_support_outside(fe, OK21, [P2_21])
            @test iseven(valuation(fe, P))
        end
    end

    # the dyadic precondition is enforced
    @test_throws "dyadic" MH._kr_algorithm_1(K7, [prime_decomposition(OK7, 3)[1][1]])
end

@testset "signed_discriminant: refactor is behaviour-neutral" begin
    # over QQ: literally (-1)^(n(n-1)/2) * integral_discriminant, on a sweep
    # covering both parities of n(n-1)/2 and both discriminant signs
    sweep = [
        GWClass(QQ[2 1; 1 3]),
        GWClass(QQ[1 4 7; 4 3 -1; 7 -1 5]),
        GWClass(QQ[9 1 7 4; 1 10 3 2; 7 3 6 7; 4 2 7 5]),
        diagonal_form(QQ, (3, -3, 2, 5, 1, -9)),
        diagonal_form(QQ, (1, -1)),
        diagonal_form(QQ, (6, 7, 22, -5, 11)),
        diagonal_form(QQ, 7),
    ]
    for beta in sweep
        n = form_rank(beta)
        @test signed_discriminant(beta) ==
              (-1)^div(n * (n - 1), 2) * integral_discriminant(beta)
    end

    # over a finite field: the (-1)-power times the diagonal product
    delta = GWClass(GF(13)[9 1 7 4; 1 10 3 2; 7 3 6 7; 4 2 7 5])
    n = form_rank(delta)
    @test signed_discriminant(delta) ==
          GF(13)((-1)^div(n * (n - 1), 2)) * prod(diagonal_entries(delta))
end

# --- NF invariants layer --------------------------------------------------------------

@testset "NF invariants: relevant primes, dyadic primes, KR form" begin
    q = diagonal_form(K7, (s7, 3, 1))
    kr = MH._kr_form(q)
    @test MH._kr_form(q) === kr                       # one fixed diagonalisation per class
    rp = relevant_primes(q)
    @test sort([minimum(P) for P in rp]) == ZZRingElem[3, 7]
    push!(rp, rp[1])                                  # returned vector is a copy
    @test length(relevant_primes(q)) == 2
    @test sort([minimum(P) for P in kr.S]) == ZZRingElem[2, 2, 3, 7]
    @test all(isone(denominator(a, maximal_order(K7))) for a in kr.diag)

    @test length(MH._dyadic_primes(K7)) == 2          # 2 split
    @test length(MH._dyadic_primes(K2)) == 1          # 2 ramified
    @test length(MH._dyadic_primes(Kc)) == 1          # x^3-3x-1 irreducible mod 2
end

@testset "Hilbert symbol and Hasse pinning (incl. ramified dyadic)" begin
    # real-place symbol: -1 iff both arguments negative — pinned against
    # Hecke's InfPlc method and against the sign definition
    for K in (K2, Kc)
        th = gen(K)
        for sigma in real_places(K),
            a in (K(2), th, -th, th - 1, K(-3)),
            b in (K(1), -th, th + 2, K(-2))
            ours = MH._hilbert_symbol(a, b, sigma)
            @test ours == hilbert_symbol(a, b, sigma)
            @test ours == ((sign(a, sigma) < 0 && sign(b, sigma) < 0) ? -1 : 1)
        end
    end

    # Hasse invariant: hand recomputation of prod_{i<j} (a_i, a_j)_v against
    # the Hecke hasse_invariant delegation — at an odd prime, at the RAMIFIED
    # dyadic prime of Q(sqrt 2), and at a real place
    OK2 = maximal_order(K2)
    P2 = prime_decomposition(OK2, 2)[1][1]
    P7odd = prime_decomposition(OK2, 7)[1][1]
    forms2 = [diagonal_form(K2, (1, s2, 2, 3 - s2)),
              diagonal_form(K2, (s2, -1, 2 + s2)),
              GWClass(matrix(K2, 2, 2, [K2(2), s2, s2, K2(3)]))]
    for beta in forms2
        D = diagonal_entries(beta)
        for v in vcat(Any[P2, P7odd], real_places(K2))
            byhand = 1
            for i in 1:(length(D) - 1), j in (i + 1):length(D)
                byhand *= MH._hilbert_symbol(D[i], D[j], v)
            end
            @test hasse_witt_invariant(beta, v) == byhand
            @test MH._hasse_witt_from_diagonal(D, v) == byhand
        end
    end

    # all splitting behaviours of 2 are pinned: ramified above (Q(sqrt 2)),
    # SPLIT here (Q(sqrt -7), both dyadic primes), INERT (Q(sqrt 5), residue
    # degree 2), and WILDLY ramified (Q(zeta_16)^+, e = 4)
    for (K, forms) in ((K7, [diagonal_form(K7, (s7, 3, 1, -2)),
                             diagonal_form(K7, (1 - s7, -1, 2 + s7))]),
                       (Kr5, [diagonal_form(Kr5, (1, sr5, 2, 3 - sr5)),
                               diagonal_form(Kr5, (sr5, -1, 2 + sr5))]),
                       (K162, [diagonal_form(K162, (w16, 1 - w16, 2)),
                               diagonal_form(K162, (2 + w16, -1, w16^2 - 2))]))
        dys = MH._dyadic_primes(K)
        @test length(dys) == (K === K7 ? 2 : 1)
        for beta in forms, P in dys
            D = diagonal_entries(beta)
            byhand = 1
            for i in 1:(length(D) - 1), j in (i + 1):length(D)
                byhand *= MH._hilbert_symbol(D[i], D[j], P)
            end
            @test hasse_witt_invariant(beta, P) == byhand
            @test MH._hasse_witt_from_diagonal(D, P) == byhand
        end
    end
end

@testset "Hilbert reciprocity (dyadic-nontrivial inputs)" begin
    for (K, s) in ((K2, s2), (K7, s7), (K5, s5),
                   (Ks10, s10), (K162, w16), (K2i, s2i))
        OK = maximal_order(K)
        cases = [(s, 1 + s), (K(2), 1 + s), (2 + s, s), (K(-1), s),
                 (3 - s, 2 * s), (s * (1 + s), K(2)), (K(-2), s - 1)]
        for (a, b) in cases
            places = Any[]
            for x in (a, b), (P, _) in factor(x * OK)
                any(Q -> Q == P, places) || push!(places, P)
            end
            for P in MH._dyadic_primes(K)
                any(Q -> Q == P, places) || push!(places, P)
            end
            pr = prod(MH._hilbert_symbol(a, b, P) for P in places; init = 1)
            for sigma in real_places(K)
                pr *= MH._hilbert_symbol(a, b, sigma)
            end
            @test pr == 1
        end
    end
end

@testset "exported local invariants: local_hilbert_symbol / local_anisotropic_dimension" begin
    # the exported names agree with the internals they wrap, on a sweep that
    # includes ramified (K2), split (K7), and wildly ramified (K162) dyadics
    for (K, s) in ((K2, s2), (K7, s7), (K162, w16))
        cases = [(K(-1), K(-1)), (K(2), s), (1 + s, K(-2)), (s, 3 - s)]
        for P in MH._dyadic_primes(K), (a, b) in cases
            @test local_hilbert_symbol(a, b, P) == MH._hilbert_symbol(a, b, P)
        end
        for sigma in real_places(K), (a, b) in cases
            @test local_hilbert_symbol(a, b, sigma) == MH._hilbert_symbol(a, b, sigma)
        end
    end

    # absolute pins at a split dyadic prime of K7: the completion is Q_2
    # itself, so the exported K symbol must equal the QQ-only export
    P2s = first(MH._dyadic_primes(K7))
    for (a, b) in ((-1, -1), (2, 5), (3, 7), (-2, 3), (2, 3))
        @test local_hilbert_symbol(a, b, P2s) == hilbert_symbol_padic(a, b, 2)
    end
    # real places: (-1,-1) is -1 at every real place; (-1, sqrt2) is -1 at
    # exactly one of the two real places of Q(sqrt2)
    sigmas2 = real_places(K2)
    @test all(local_hilbert_symbol(-1, -1, sigma) == -1 for sigma in sigmas2)
    @test prod(local_hilbert_symbol(-1, s2, sigma) for sigma in sigmas2) == -1
    # complex places carry the trivial symbol
    cplc = [v for v in infinite_places(K7) if !is_real(v)]
    @test local_hilbert_symbol(-1, -1, only(cplc)) == 1

    # Hilbert reciprocity through the exported name, over the wild field
    let a = w16, b = 2 + w16
        OK162 = maximal_order(K162)
        places = Any[]
        for x in (a, b), (P, _) in factor(x * OK162)
            any(Q -> Q == P, places) || push!(places, P)
        end
        for P in MH._dyadic_primes(K162)
            any(Q -> Q == P, places) || push!(places, P)
        end
        pr = prod(local_hilbert_symbol(a, b, P) for P in places; init = 1)
        for sigma in real_places(K162)
            pr *= local_hilbert_symbol(a, b, sigma)
        end
        @test pr == 1
    end

    # coercion and rejection
    @test local_hilbert_symbol(2, 3, P2s) isa Int
    @test_throws ErrorException local_hilbert_symbol(0, 3, P2s)
    @test_throws ErrorException local_hilbert_symbol(s2, 3, P2s)   # wrong field
    @test_throws ErrorException local_hilbert_symbol(2, 3, MH._dyadic_primes(K7)[1]^2)

    # local_anisotropic_dimension: agreement with the internal on a sweep ...
    for (K, form) in ((K2, diagonal_form(K2, (1, s2, 2, 3 - s2))),
                      (K7, diagonal_form(K7, (s7, 3, 1, -2))),
                      (K162, diagonal_form(K162, (w16, 1 - w16, 2))))
        for P in MH._dyadic_primes(K)
            @test local_anisotropic_dimension(form, P) ==
                  MH._local_anisotropic_dimension(form, P)
        end
        for sigma in real_places(K)
            @test local_anisotropic_dimension(form, sigma) ==
                  MH._local_anisotropic_dimension(form, sigma)
        end
        # ... and the global dimension dominates every local one
        for P in MH._kr_form(form).S
            @test anisotropic_dimension(form) >= local_anisotropic_dimension(form, P)
        end
    end
    # absolute pins: <1,1,1,1> over K7 is the quaternion norm form at the
    # split dyadics (K_P = Q_2), hyperbolic forms vanish locally everywhere,
    # a real place reads |signature|, a complex place the rank mod 2
    @test local_anisotropic_dimension(diagonal_form(K7, (1, 1, 1, 1)), P2s) == 4
    @test local_anisotropic_dimension(hyperbolic_form(K7, 4), P2s) == 0
    @test local_anisotropic_dimension(diagonal_form(K2, (1, 1)), sigmas2[1]) == 2
    cplc32 = [v for v in infinite_places(K32) if !is_real(v)]
    @test local_anisotropic_dimension(diagonal_form(K32, (1, 1, 1)), only(cplc32)) == 1
    # rejection: non-prime ideal, place of a different field
    @test_throws ErrorException local_anisotropic_dimension(
        diagonal_form(K7, (1, 1)), MH._dyadic_primes(K7)[1]^2)
    @test_throws ErrorException local_anisotropic_dimension(
        diagonal_form(K7, (1, 1)), MH._dyadic_primes(K2)[1])
end

@testset "Hasse product formula (product over all places is 1)" begin
    for (K, form) in ((K2, diagonal_form(K2, (1, s2, 2, 3 - s2))),
                      (K7, diagonal_form(K7, (s7, 3, 1, -2))),
                      (Kc, diagonal_form(Kc, (ac, ac - 1, 2))),
                      (Ks10, diagonal_form(Ks10, (s10, 2, -1))),
                      (K162, diagonal_form(K162, (w16, 1 - w16, 2))),
                      (K2i, diagonal_form(K2i, (s2i, 3, 1 + s2i))))
        kr = MH._kr_form(form)
        pr = prod(hasse_witt_invariant(form, P) for P in kr.S; init = 1)
        for sigma in real_places(K)
            pr *= hasse_witt_invariant(form, sigma)
        end
        @test pr == 1
    end

end

@testset "Agreement with Hecke: is_isometric and per-place signature" begin
    for (K, s) in ((K2, s2), (K7, s7))
        f1 = diagonal_form(K, (1, s, -2))
        M1 = gw_matrix(f1)
        T = matrix(K, 3, 3, [K(1), K(2), K(0), K(0), K(1), s, K(0), K(0), K(1)])
        f2 = GWClass(transpose(T) * M1 * T)           # same class, different Gram
        f3 = diagonal_form(K, (1, s, 2))
        for (a, b) in ((f1, f2), (f1, f3), (f2, f3), (f1, f1))
            @test is_isomorphic_form(a, b) ==
                  is_isometric(quadratic_space(K, gw_matrix(a)),
                               quadratic_space(K, gw_matrix(b)))
        end
        @test is_isomorphic_form(f1, f2)
    end
    beta = diagonal_form(K2, (1, s2, -2, 3 - s2))
    V = quadratic_space(K2, gw_matrix(beta))
    for sigma in real_places(K2)
        st = signature_tuple(V, sigma)                # (pos, zero, neg)
        @test form_signature(beta, sigma) == st[1] - st[3]
    end
end

@testset "_disc_rep_normalized and _local_uniformizer" begin
    OK2 = maximal_order(K2)
    P2 = prime_decomposition(OK2, 2)[1][1]
    P7o = prime_decomposition(OK2, 7)[1][1]
    S = [P2, P7o]
    # uniformizer(P) alone does not control the valuation at the other
    # members of S; the CRT-corrected one does
    pi2 = MH._local_uniformizer(P2, S)
    @test valuation(pi2, P2) == 1 && valuation(pi2, P7o) == 0
    pi7 = MH._local_uniformizer(P7o, S)
    @test valuation(pi7, P7o) == 1 && valuation(pi7, P2) == 0

    # includes the ill-formed case of memo §8: odd NEGATIVE valuation at the
    # ramified dyadic prime
    for dd in ((3 - s2) // s2^3, K2(2) * 49, s2^5 * (3 - s2), K2(9 // 4))
        d = K2(dd)
        dp = MH._disc_rep_normalized(d, S)
        @test isone(denominator(dp, OK2))             # corrected invariant: integral,
        for P in S                                    # parity preserved on S — NOT {0,1}
            @test iseven(valuation(dp, P) - valuation(d, P))
        end
        @test is_square(dp * d)                       # same square class
    end
end

@testset "transfer_gw over number fields" begin
    # A number field L over its base field k is a finite étale k-algebra, and
    # transfer_gw computes the transfer along Tr_{L/k} with Hecke's trace on
    # Hecke's basis of L. Each value is checked against a route sharing no code
    # with that: the étale algebra k[t]/(f) presenting L, whose transfer is
    # built from multiplication matrices, Hecke's restriction of scalars, or a
    # computation by hand.
    Qx, x = QQ["x"]
    K, a = number_field(x^2 - 2, "a")

    # By hand, on the basis 1, sqrt2: Tr(1) = 2 and Tr(sqrt2) = 0, so tr<1> is
    # diag(2, 4), and tr<3 + 5 sqrt2> is diag(2, 4) [3 10; 5 3] = [6 20; 20 12].
    @test gw_matrix(transfer_gw(GWClass(matrix(K, 1, 1, [K(1)])))) == QQ[2 0; 0 4]
    @test gw_matrix(transfer_gw(GWClass(matrix(K, 1, 1, [3 + 5 * a])))) ==
          QQ[6 20; 20 12]

    for f in (x^2 - 2, x^2 + 1, x^3 - 2, x^4 - x^3 - x^2 + x + 1)
        L, c = number_field(f, "c")
        S, (t,) = polynomial_ring(QQ, ["t"])
        A, _ = quo(S, ideal(S, [f(t)]))
        # c goes to t, so an element goes to its power-basis coordinates in t
        into_A(e) = (coords = coordinates(e);
                     A(sum(coords[i] * t^(i - 1) for i in eachindex(coords))))
        for B in (matrix(L, 1, 1, [2 + c]),
                  matrix(L, 2, 2, [L(3), L(1), L(1), 1 + c]),
                  matrix(L, 3, 3, [L(1), c, L(0), c, L(2), L(1), L(0), L(1), 5 - c]))
            G = gw_matrix(transfer_gw(GWClass(B)))
            @test base_ring(G) === QQ
            @test nrows(G) == nrows(B) * degree(L)
            @test G == gw_matrix(transfer_gw(GWClass(map(into_A, B))))
            @test G == map_entries(QQ, gram_matrix(
                restrict_scalars(quadratic_space(L, B), QQ)[1]))
        end
    end

    # A relative field transfers to its base field: for L = K(sqrt3) over
    # K = QQ(sqrt2), by hand tr<1> = diag(2, 6) over K on the basis 1, sqrt3.
    Kz, z = K["z"]
    L, r = number_field(z^2 - 3, "r")
    @test gw_matrix(transfer_gw(GWClass(matrix(L, 1, 1, [L(1)])))) ==
          matrix(K, 2, 2, [2, 0, 0, 6])
    Ky, (y,) = polynomial_ring(K, ["y"])
    E, _ = quo(Ky, ideal(Ky, [y^2 - 3]))
    into_E(e) = (coords = coordinates(e);
                 E(sum(coords[i] * y^(i - 1) for i in eachindex(coords))))
    for B in (matrix(L, 1, 1, [1 + r]), matrix(L, 2, 2, [L(3), r, r, 1 + a * r]))
        G = gw_matrix(transfer_gw(GWClass(B)))
        @test base_ring(G) === K
        @test G == gw_matrix(transfer_gw(GWClass(map(into_E, B))))
    end

    # A non-simple field, which the invariant layers refuse: on Hecke's basis
    # 1, u, v, uv of QQ(sqrt2, sqrt3), tr<1> is the trace form diag(4, 8, 12, 24),
    # as for the étale algebra QQ[u, v]/(u^2 - 2, v^2 - 3) on its standard
    # monomials 1, u, v, uv.
    N, _ = number_field([x^2 - 2, x^2 - 3], ["u", "v"])
    S2, (u, v) = polynomial_ring(QQ, ["u", "v"])
    biquadratic, _ = quo(S2, ideal(S2, [u^2 - 2, v^2 - 3]))
    trace_form = diagonal_matrix(QQ, [4, 8, 12, 24])
    @test gw_matrix(transfer_gw(GWClass(matrix(N, 1, 1, [N(1)])))) == trace_form
    @test gw_matrix(transfer_gw(diagonal_form(biquadratic, biquadratic(1)))) ==
          trace_form
end
