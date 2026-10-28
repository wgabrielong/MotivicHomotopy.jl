# The certify = true paths (certified.jl and the certified section of
# numerical.jl). Everything asserted here is exact: the certified claim is
# relative to the float input read exactly (each Float64 / ComplexF64 as the
# dyadic (Gaussian) rational it is), so scalars are tested with ==, never with
# a tolerance. Default-path values asserted alongside are bug-compat locks:
# they pin the M2-ported behavior the certified branch deliberately does not
# change.

using Test
using MotivicHomotopy
using Oscar
import HomotopyContinuation as HC
import LinearAlgebra

HC.@var x

@testset "certified CC unstable degrees: M2 captures become exact" begin
    f = (x - 1) * (x - 2) * (x - 3)
    g = (x - 1) * (x - 4)
    G = global_unstable_A1_degree(f, g; certify = true)
    @test gw_matrix(G) == Matrix{ComplexF64}(LinearAlgebra.I, 2, 2)
    @test gw_scalar(G) == -2.0 + 0.0im            # default path: -1.9999999999999987…
    @test gw_scalar(global_unstable_A1_degree(f, g)) != -2.0 + 0.0im   # lock the drift

    fd = (x - 1)^2 * (x - 5)
    gd = x - 2
    GG = global_unstable_A1_degree(fd, gd; certify = true)
    @test size(gw_matrix(GG), 1) == 3
    @test gw_scalar(GG) == -3.0 + 0.0im           # M2 capture: -2.9999999999999982

    L1 = local_unstable_A1_degree(fd, gd, 1; certify = true)
    @test size(gw_matrix(L1), 1) == 2
    @test gw_scalar(L1) == -0.0625 + 0.0im        # -1/16 exactly
    L5 = local_unstable_A1_degree(fd, gd, 5; certify = true)
    @test size(gw_matrix(L5), 1) == 1
    @test gw_scalar(L5) == 0.1875 + 0.0im         # 3/16 exactly
    L2 = local_unstable_A1_degree(f, g, 2; certify = true)
    @test gw_scalar(L2) == 2.0 + 0.0im
end

@testset "certified CC unstable degrees: high multiplicity" begin
    # (x-2)^5 / (x-3): rank 5, scalar (-1)^(5*4/2) * (2-3)^5 = -1
    G5 = global_unstable_A1_degree((x - 2)^5, (x - 3); certify = true)
    @test size(gw_matrix(G5), 1) == 5
    @test gw_scalar(G5) == -1.0 + 0.0im

    # locally at 2 the numerical quintuple-root cluster is ~1e-3 wide — far
    # beyond the default 1e-6 — so the default reading counts a wrong
    # multiplicity and the divergence guard must refuse
    @test_throws "cannot certify" local_unstable_A1_degree(
        (x - 2)^5, (x - 3), 2; certify = true)
    # widening linear_tolerance past the cluster makes the readings agree
    L = local_unstable_A1_degree((x - 2)^5, (x - 3), 2;
                                 certify = true, linear_tolerance = 1e-2)
    @test size(gw_matrix(L), 1) == 5
    @test gw_scalar(L) == -1.0 + 0.0im            # (-1)^10 * ((2-3)/1)^5
end

@testset "certified CC unstable degrees: non-monic numerators" begin
    # The pair (f, g) stands for the rational function f/g, normalized by the
    # leading coefficient of f as on the exact and float paths. For
    # f/g = (x^2-3x+2)/((x-4)/2) the scalars are -Res = -(g(1) g(2)) = -3/2
    # globally, then 3/2 at 1 and -1 at 2; all are dyadic, so certify = true
    # must return them exactly.
    for (f, g) in ((2x^2 - 6x + 4, x - 4), (4x^2 - 12x + 8, 2x - 8))
        @test gw_scalar(global_unstable_A1_degree(f, g; certify = true)) == -1.5 + 0.0im
        @test gw_scalar(local_unstable_A1_degree(f, g, 1; certify = true)) == 1.5 + 0.0im
        @test gw_scalar(local_unstable_A1_degree(f, g, 2; certify = true)) == -1.0 + 0.0im
    end
    # a negative leading coefficient: f/g = (x^3 - x)/(-3/2), so the global
    # scalar is -(-3/2)^3 = 27/8 and the local one at 0 is
    # (-3/2)/((0-1)(0+1)) = 3/2
    three = HC.Expression(3)
    @test gw_scalar(global_unstable_A1_degree(-2x^3 + 2x, three; certify = true)) ==
          3.375 + 0.0im
    @test gw_scalar(local_unstable_A1_degree(-2x^3 + 2x, three, 0; certify = true)) ==
          1.5 + 0.0im
    # scaling f and g together does not change f/g, so it does not change the
    # certified scalar either
    f = (x - 1) * (x - 2) * (x - 3)
    g = (x - 1) * (x - 4)
    unscaled = gw_scalar(global_unstable_A1_degree(f, g; certify = true))
    for c in (2, -3)
        @test gw_scalar(global_unstable_A1_degree(c * f, c * g; certify = true)) == unscaled
        @test gw_scalar(local_unstable_A1_degree(c * f, c * g, 2; certify = true)) ==
              gw_scalar(local_unstable_A1_degree(f, g, 2; certify = true))
    end
end

@testset "certified CC unstable degrees: divergence guard" begin
    # near-common root: the tolerance reading cancels 1 against 1+1e-8 (rank
    # 1), the exact reading has a trivial gcd (rank 2)
    @test_throws "cannot certify" global_unstable_A1_degree(
        (x - 1) * (x - 3), x - (1 + 1e-8); certify = true)
    # escape hatch: a tolerance below the root gap accepts the exact reading
    Ge = global_unstable_A1_degree((x - 1) * (x - 3), x - (1 + 1e-8);
                                   certify = true, linear_tolerance = 1e-12)
    @test size(gw_matrix(Ge), 1) == 2

    # near-double root: exact multiplicity of 1 is 1, tolerance counts 2
    @test_throws "cannot certify" local_unstable_A1_degree(
        (x - 1) * (x - (1 + 1e-8)) * (x - 5), x - 2, 1; certify = true)

    # r within tolerance of a root but not exactly a root: the readings
    # disagree on whether r is a zero at all
    @test_throws "cannot certify" local_unstable_A1_degree(
        (x - 1) * (x - 3), x - 2, 1 + 1e-9; certify = true)
    # far from every root both readings agree it is not a zero: plain error
    @test_throws "not a zero" local_unstable_A1_degree(
        (x - 1) * (x - 3), x - 2, 1 + 1e-3; certify = true)

    # exact common root: both readings cancel it, guard passes (also covered
    # by the capture testset above; kept explicit for the error paths)
    @test_throws "not pointed" global_unstable_A1_degree(
        x - 1, (x - 1) * (x - 4); certify = true)
end

@testset "certified CC unstable degrees: unrepresentable scalars" begin
    # overflow: scalar carries lc(g)^d = (1e200)^2
    @test_throws "not representable" global_unstable_A1_degree(
        (x - 1) * (x - 2), 1e200 * (x - 4); certify = true)
    # underflow: scalar is -(1e-40)^10 = -1e-400, below the subnormal range.
    # linear_tolerance = 0.5 makes the tolerance reading agree on the
    # multiplicity (the numerical cluster of the tenfold root 0 is ~1e-2 wide)
    @test_throws "not representable" local_unstable_A1_degree(
        x^10, 0 * x + 1e-40, 0; certify = true, linear_tolerance = 0.5)
end

@testset "certified float-Gram rank" begin
    M = [1.0 0.0; 0.0 1e-20]
    @test form_rank(M) == 1                       # bug-compat lock: SVD tolerance
    @test form_rank(M; certify = true) == 2

    Mc = ComplexF64[1.0 im; im -1.0]              # exact determinant 0
    @test form_rank(Mc; certify = true) == 1
    @test anisotropic_dimension(Mc; certify = true) == 1
    # rank decides CC isometry: certified sees the degeneracy the SVD default
    # also sees here, but the near-singular case diverges
    Mc2 = ComplexF64[1.0 0.0; 0.0 1e-20]
    @test !is_isomorphic_form(Mc2, Matrix{ComplexF64}(LinearAlgebra.I, 2, 2))
    @test is_isomorphic_form(Mc2, Matrix{ComplexF64}(LinearAlgebra.I, 2, 2);
                             certify = true)
end

@testset "certified float-Gram signature (Hilbert matrices)" begin
    hilbert(n) = [1.0 / (i + j - 1) for i in 1:n, j in 1:n]
    # Hilbert matrices are positive definite, so the certified signature is
    # provably n and the certified rank provably n
    for n in (10, 12)
        @test form_signature(GWClass(hilbert(n)); certify = true) == n
    end
    @test form_rank(hilbert(12); certify = true) == 12
    # bug-compat locks of the default path at n = 12: the float congruence
    # happens to keep all signs, but the SVD rank drops to 11
    @test form_signature(GWClass(hilbert(12))) == 12
    @test form_rank(hilbert(12)) == 11
end

@testset "certified GW / GWu constructors" begin
    M = LinearAlgebra.diagm([1e-200, 1e-200])     # exact det 1e-400 > 0
    @test_throws ErrorException GWClass(M)        # bug-compat: float det underflows
    @test gw_matrix(GWClass(M; certify = true)) == M

    @test_throws ErrorException GWuClass(M, 1e-300)
    @test gw_scalar(GWuClass(M, 1e-300; certify = true)) == 1e-300
    # the omitted scalar is the exact determinant correctly rounded; here that
    # rounding underflows to zero, which is refused
    @test_throws "not representable" GWuClass(M; certify = true)

    # a genuinely sign-flipped scalar is still rejected under certify
    @test_throws ErrorException GWuClass([2.0 0.0; 0.0 3.0], -6.0; certify = true)
    @test gw_scalar(GWuClass([2.0 0.0; 0.0 3.0]; certify = true)) == 6.0
    @test gw_scalar(GWuClass(ComplexF64[2.0 0.0; 0.0 3.0]; certify = true)) ==
        6.0 + 0.0im
    @test gw_scalar(GWuClass(GWClass([2.0 0.0; 0.0 3.0]); certify = true)) == 6.0

    # non-finite input is refused outright
    @test_throws "non-finite" form_rank([NaN 0.0; 0.0 1.0]; certify = true)
    @test_throws "non-finite" GWClass([NaN 0.0; 0.0 1.0]; certify = true)
end

@testset "certify is a no-op on exact base fields" begin
    @test form_signature(diagonal_form(QQ, (1, -1, 1)); certify = true) == 1
    @test anisotropic_dimension(diagonal_form(QQ, (1, -1, 2)); certify = true) == 1
    @test is_isotropic_form(diagonal_form(QQ, (1, -1)); certify = true)
    @test gw_scalar(GWuClass(QQ[2 0; 0 3]; certify = true)) == QQ(6)
    u = GWuClass(QQ[1 0; 0 1], 4)
    @test is_isomorphic_form(u, GWuClass(QQ[1 0; 0 1], 4); certify = true)
end

@testset "certified anisotropic part and GWu isometry over floats" begin
    A = [1.0 0.0; 0.0 -1e-20]
    @test anisotropic_part(A; certify = true) == zeros(Float64, 0, 0)
    @test anisotropic_part(A) == zeros(Float64, 0, 0)   # bug-compat lock
    @test gw_witt_index(GWClass(A); certify = true) == 1

    # certify governs the stable parts; the k^x-scalar comparison stays a
    # linear_tolerance comparison even under certify (documented)
    u1 = GWuClass([2.0 0.0; 0.0 3.0], 6.0)
    u2 = GWuClass([1.0 0.0; 0.0 6.0], 6.0 + 1e-9)
    @test is_isomorphic_form(u1, u2; certify = true)
end
