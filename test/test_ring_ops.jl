# The GW semiring surface: identities, integer multiple, tensor power, and
# base change — generic across base fields — plus the number-field pinning of
# the operations the API coverage audit
# found working but untested over K: gw_tensor_product and its * alias, the
# additive and multiplicative identities, and pfister_form.

using Test
using Oscar
using MotivicHomotopy

const MH = MotivicHomotopy

include("nf_test_fields.jl")

# the full number-field matrix, with a nonrational generator per field
# (KQ's generator is 1, so it gets a plain nonsquare instead)
_matrix_fields() = ((K7, s7), (K2, s2), (Kr5, sr5), (K5, s5), (K21, s21),
                    (K34, s34), (Kc, ac), (K32, a32), (Kq, aq),
                    (Ks10, s10), (K162, w16), (K2i, s2i), (KQ, KQ(5)))

@testset "additive and multiplicative identities over the full field matrix" begin
    for (K, s) in _matrix_fields()
        beta = diagonal_form(K, (1, 2, s))
        empty = diagonal_form(K, ())
        @test form_rank(empty) == 0
        @test gw_direct_sum(empty, beta) == beta
        @test gw_direct_sum(beta, empty) == beta
        one_form = diagonal_form(K, 1)
        @test gw_tensor_product(one_form, beta) == beta
        @test gw_tensor_product(beta, one_form) == beta
    end
end

@testset "tensor product over the full field matrix (properties)" begin
    for (K, s) in _matrix_fields()
        beta = diagonal_form(K, (1, 2, s))     # rank 3
        gamma = diagonal_form(K, (3, -s))      # rank 2
        t = gw_tensor_product(beta, gamma)
        @test t == beta * gamma                               # the * alias
        # rank multiplies
        @test form_rank(t) == form_rank(beta) * form_rank(gamma)
        # det(A ⊗ B) = det(A)^m det(B)^n for A n x n, B m x m — exactly, so
        # the discriminant square class behaves as it must
        @test det(gw_matrix(t)) ==
              det(gw_matrix(beta))^form_rank(gamma) *
              det(gw_matrix(gamma))^form_rank(beta)
        # tensoring ISOMORPHIC forms gives isomorphic results — the property
        # a bad diagonalisation would break and a rank check cannot see:
        # beta2 is beta in a different basis, gamma2 rescales gamma's
        # entries by nonzero squares
        M = gw_matrix(beta)
        T = matrix(K, 3, 3, [K(1), K(1), K(0),
                             K(0), K(1), s,
                             K(0), K(0), K(1)])
        beta2 = GWClass(transpose(T) * M * T)
        gamma2 = diagonal_form(K, (3 * (1 + s)^2, -s * K(4)))
        @test is_isomorphic_form(beta, beta2)
        @test is_isomorphic_form(gamma, gamma2)
        @test is_isomorphic_form(gw_tensor_product(beta, gamma),
                                 gw_tensor_product(beta2, gamma2))
    end
end

@testset "pfister_form over the full field matrix" begin
    for (K, s) in _matrix_fields()
        # <<s, 2>> = <1,-s> ⊗ <1,-2>, rank 4, square discriminant
        p2 = pfister_form(K, (s, 2))
        @test form_rank(p2) == 4
        @test is_square(det(gw_matrix(p2)))
        @test p2 == gw_tensor_product(diagonal_form(K, (1, -s)),
                                      diagonal_form(K, (1, -2)))
        # one-fold: <<a>> = <1, -a>; <<square>> is hyperbolic
        @test pfister_form(K, (s,)) == diagonal_form(K, (1, -s))
        @test is_isomorphic_form(pfister_form(K, (s^2,)), hyperbolic_form(K))
        # a Pfister form is hyperbolic or anisotropic, never in between
        w = gw_witt_index(p2)
        @test w == 0 || w == 2
    end
end

@testset "integer multiple n * beta (all base fields)" begin
    K7beta = diagonal_form(K7, (1, s7))
    backends = Any[
        diagonal_form(QQ, (1, 3)),
        diagonal_form(GF(7), (1, 3)),
        K7beta,
        diagonal_form(K162, (1, w16)),
        GWClass([2.0 0.0; 0.0 3.0]),
        GWClass(ComplexF64[2 0; 0 3]),
    ]
    for beta in backends
        z = 0 * beta
        @test form_rank(z) == 0
        @test 1 * beta == beta
        @test 2 * beta == gw_direct_sum(beta, beta)
        @test beta * 2 == 2 * beta
        @test form_rank(3 * beta) == 3 * form_rank(beta)
        @test_throws ErrorException (-1) * beta
    end
    # over K the multiple interacts correctly with the Witt machinery:
    # 2 * H is totally hyperbolic, and n * beta scales the Witt index
    @test gw_witt_index(2 * hyperbolic_form(K7)) == 2
    @test gw_witt_index(3 * hyperbolic_form(K2)) == 3
end

@testset "tensor power beta^n (all base fields)" begin
    backends = Any[
        diagonal_form(QQ, (1, 3)),
        diagonal_form(GF(7), (1, 3)),
        diagonal_form(K7, (1, s7)),
        diagonal_form(K162, (1, w16)),
        GWClass([2.0 0.0; 0.0 3.0]),
        GWClass(ComplexF64[2 0; 0 3]),
    ]
    for beta in backends
        u = beta^0
        @test form_rank(u) == 1
        @test beta^1 == beta
        @test beta^2 == gw_tensor_product(beta, beta)
        @test form_rank(beta^3) == form_rank(beta)^3
        @test_throws ErrorException beta^(-1)
    end
    # beta^0 is the multiplicative identity <1>
    beta = diagonal_form(K2, (1, s2))
    @test beta^0 == diagonal_form(K2, 1)
    @test gw_tensor_product(beta^0, beta) == beta
end

@testset "gw_base_change GW(Q) -> GW(K)" begin
    beta = diagonal_form(QQ, (3, -5, 7))
    gamma = diagonal_form(QQ, (1, -1))
    for (K, s) in _matrix_fields()
        bc = gw_base_change(K, beta)
        @test gw_base_field(bc) === K
        # commutes with rank
        @test form_rank(bc) == form_rank(beta)
        # discriminant up to squares: the image's signed discriminant is the
        # rational one, read in K
        @test is_square(signed_discriminant(bc) * K(signed_discriminant(beta)))
        # commutes with orthogonal sum, literally
        @test gw_base_change(K, gw_direct_sum(beta, gamma)) ==
              gw_direct_sum(gw_base_change(K, beta), gw_base_change(K, gamma))
        # and with the tensor product
        @test gw_base_change(K, gw_tensor_product(beta, gamma)) ==
              gw_tensor_product(gw_base_change(K, beta), gw_base_change(K, gamma))
        # hyperbolic goes to hyperbolic
        @test is_isomorphic_form(gw_base_change(K, gamma), hyperbolic_form(K))
    end
    # base change can merge square classes: <2> = <1> over Q(sqrt2), not
    # over Q(sqrt-7)
    two = diagonal_form(QQ, (2,))
    @test is_isomorphic_form(gw_base_change(K2, two), diagonal_form(K2, 1))
    @test !is_isomorphic_form(gw_base_change(K7, two), diagonal_form(K7, 1))
    # only GW(Q) -> GW(K) is defined
    @test_throws ErrorException gw_base_change(K2, diagonal_form(K2, 1))
end

@testset "GW^u over K: validation, comparison, and the layer around them" begin
    # the two previously broken feet of the GW^u layer over K: construction
    # must REJECT a scalar in the wrong square class (previously accepted
    # with a warning via the étale-algebra fallback), and is_isomorphic_form
    # must work (previously errored: no K arm in _scalar_field_supported)
    for (K, s) in _matrix_fields()
        M = gw_matrix(diagonal_form(K, (2, 3 + s)))
        d = det(M)
        # scalar off by the SQUARE (1+s)^2: same class, accepted
        u = GWuClass(M, d * (1 + s)^2)
        @test gw_scalar(u) == d * (1 + s)^2
        @test stable_part(u) isa GWClass && form_rank(stable_part(u)) == 2
        # scalar off by the NONSQUARE s: wrong class, rejected — not warned
        @test !is_square(K(s))
        @test_throws ErrorException GWuClass(M, d * s)

        # constructors and accessors
        du = diagonal_unstable_form(K, (2, 3 + s))
        @test gw_scalar(du) == d
        hu = hyperbolic_unstable_form(K)
        @test gw_scalar(hu) == K(-1)
        @test is_isomorphic_form(stable_part(hu), hyperbolic_form(K))

        # direct sum: matrices block-sum, scalars multiply
        two = gw_direct_sum(du, hu)
        @test gw_scalar(two) == d * K(-1)
        @test gw_matrix(two) ==
              gw_matrix(gw_direct_sum(stable_part(du), stable_part(hu)))

        # comparison: same stable class in another basis, equal scalar -> true
        T = matrix(K, 2, 2, [K(1), s, K(0), K(1)])
        u2 = GWuClass(transpose(T) * M * T, d)
        u1 = GWuClass(M, d)
        @test is_isomorphic_form(u1, u2)
        # equal stable class, scalars differing by a square -> FALSE: the
        # k^x-coordinate compares exactly, as over Q (M2's rule)
        u4 = GWuClass(M, d * K(4))
        @test !is_isomorphic_form(u1, u4)

        # diagonal_class / sum_decomposition carry the scalar through, and
        # the M2-mirrored cache overwrite holds for GWu over K too
        dc = diagonal_class(u)
        @test gw_scalar(dc) == gw_scalar(u)
        @test is_diagonal(gw_matrix(dc))
        @test is_isomorphic_form(dc, u)
        r = sum_decomposition(u)
        @test gw_scalar(r) == gw_scalar(u)
        @test is_isomorphic_form(r, u)
        @test diagonal_class(u) === r
    end
    # different base fields never compare
    @test_throws ErrorException is_isomorphic_form(
        diagonal_unstable_form(K7, (1, s7)), diagonal_unstable_form(K5, (1, s5)))
end

@testset "divisorial_sum over K: validated, not merely run" begin
    # Until the two GWu feet were repaired this output could not be checked
    # at all. Two independent validations:
    # (1) the [I+24] definition, recomputed explicitly in the test: block-sum
    #     stable part, scalar = prod(scalars) * prod_{i<j} (r_i - r_j)^(2 m_i m_j)
    for (K, s) in ((K7, s7), (K2, s2), (Ks10, s10), (K162, w16))
        a1 = diagonal_unstable_form(K, (2 + s,))          # rank 1
        a2 = diagonal_unstable_form(K, (3, -s))           # rank 2
        a3 = diagonal_unstable_form(K, (K(5),))           # rank 1
        roots = [K(0), K(1), s]
        ds = divisorial_sum([a1, a2, a3], roots)
        expected_scalar = gw_scalar(a1) * gw_scalar(a2) * gw_scalar(a3) *
            (roots[1] - roots[2])^(2 * 1 * 2) *
            (roots[1] - roots[3])^(2 * 1 * 1) *
            (roots[2] - roots[3])^(2 * 2 * 1)
        @test gw_scalar(ds) == expected_scalar
        @test gw_matrix(ds) == block_diagonal_matrix(
            [gw_matrix(a1), gw_matrix(a2), gw_matrix(a3)])
        # the result is a well-defined unstable class: reconstructing it
        # re-runs the (now correct) validation
        @test GWuClass(gw_matrix(ds), gw_scalar(ds)) == ds
        @test is_isomorphic_form(ds, ds)
    end
    # (2) the M2-validated Q path as an oracle: the same rational data
    #     through QQ and through KQ = rationals-as-number-field must agree
    #     under coercion (the docstring example: local degrees of
    #     (x^2+x-2)/(3x+5) at -2 and 1, divisorial sum = global degree)
    aQ = GWuClass(matrix(QQ, 1, 1, [QQ(1, 3)]))
    bQ = GWuClass(matrix(QQ, 1, 1, [QQ(8, 3)]))
    dsQ = divisorial_sum([aQ, bQ], [-2, 1])
    aK = GWuClass(matrix(KQ, 1, 1, [KQ(1) // 3]))
    bK = GWuClass(matrix(KQ, 1, 1, [KQ(8) // 3]))
    dsK = divisorial_sum([aK, bK], [KQ(-2), KQ(1)])
    @test gw_scalar(dsK) == KQ(gw_scalar(dsQ))
    @test gw_matrix(dsK) == map_entries(KQ, gw_matrix(dsQ))
    @test gw_scalar(dsQ) == QQ(8)                       # the docstring value
end
