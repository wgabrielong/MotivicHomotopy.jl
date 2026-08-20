# The ℂ numerical backend, on HomotopyContinuation.jl. Numerical computation is
# needed ONLY for the unstable degrees over ℂ: there the class carries a scalar
# (a k^× factor), so the actual complex roots must be found. These add methods to
# the unstable-degree functions for floating-point input; HC names are qualified
# `HC.` so they never collide with Oscar's exports.
#
# The stable degrees do NOT use this backend. Over ℂ a symmetric bilinear form is
# determined by its rank, and the stable A¹-degree is the identity form of rank
# equal to the dimension of the (global or local) coordinate algebra — a discrete
# invariant obtained from the exact computation over ℚ, with no need to solve for
# roots. See `global_A1_degree` / `local_A1_degree` in degrees.jl.
#
# The unstable methods below find univariate roots via the companion matrix,
# cancel roots of numerator and denominator that agree within `linear_tolerance`,
# and form the scalar from the root-product resultant.
#
# With `certify = true` the same two entry points instead run the certified
# path (see the "certified path" section): the coefficients are read exactly
# and everything is computed over ℚ(i), with the default pipeline's discrete
# decisions re-derived as a divergence guard.

# ---------- univariate polynomial machinery ----------------------------------------

_to_c64(c) = c isa HC.Expression ?
    ComplexF64(HC.ModelKit.to_number(c)) : ComplexF64(c)

# ascending coefficient vector of a univariate expression
function _poly_coeffs(f::HC.Expression, var)
    E, C = HC.exponents_coefficients(f, [var])
    d = maximum(E)
    cv = zeros(ComplexF64, d + 1)
    for (k, e) in enumerate(E[1, :])
        cv[e + 1] = _to_c64(C[k])
    end
    cv
end

# degree-many roots, via the companion matrix
function _poly_roots(f::HC.Expression, var)
    cv = _poly_coeffs(f, var)
    d = length(cv) - 1
    d <= 0 && return ComplexF64[]
    C = zeros(ComplexF64, d, d)
    for i in 2:d
        C[i, i - 1] = 1
    end
    for i in 1:d
        C[i, d] = -cv[i] / cv[d + 1]
    end
    LinearAlgebra.eigvals(C)
end

_lead_coefficient(f::HC.Expression, var) = _poly_coeffs(f, var)[end]

function _the_variable(f::HC.Expression, g::HC.Expression)
    vars = unique!(vcat(HC.variables(f), HC.variables(g)))
    length(vars) == 1 ||
        error("both input polynomials must be defined over the same univariate polynomial ring")
    vars[1]
end

# cancel roots of L1 against L2 within tol, pairwise
function _remove_common_approx(L1, L2, tol)
    L1new = ComplexF64[]
    L2new = copy(L2)
    for r in L1
        i = findfirst(s -> abs(r - s) < tol, L2new)
        if i === nothing
            push!(L1new, r)
        else
            deleteat!(L2new, i)
        end
    end
    (L1new, L2new)
end

# ---------- certified path (certify = true) -----------------------------------------
# Deliberate extension beyond M2, which has no certification. The coefficients
# are read through the same `_to_c64` conversion as the default path — so
# certification refers to exactly the polynomial the default path sees — and
# then lifted exactly into ℚ(i) (certified.jl). From there NO root is ever
# found: the tolerance cancellation becomes the exact gcd, rank and
# multiplicity become exact degrees, and the scalar becomes a resultant /
# evaluation formula, so every discrete output is exact and the returned
# scalar is correctly rounded. Before returning, the default pipeline's
# discrete decisions are recomputed and compared; a mismatch means the input
# has roots the two readings cluster differently, and certification cannot
# know which reading the caller meant — so it errors instead of guessing.

# f as an exact element of K[x], K = ℚ(i), from its ComplexF64 coefficients
function _exact_qi_poly(Kx::PolyRing, f::HC.Expression, var)
    ascending = _poly_coeffs(f, var)
    Kx([_exact_qi(c) for c in ascending])
end

function _certified_divergence_error(what::String, exact_value, tolerance_value,
                                     linear_tolerance::Real)
    error("cannot certify: reading the coefficients as exact gives $what = " *
          "$exact_value, but clustering roots within linear_tolerance = " *
          "$linear_tolerance gives $what = $tolerance_value. The input has " *
          "roots closer than linear_tolerance that are not exactly equal, or " *
          "a repeated root whose numerical root cluster is wider than " *
          "linear_tolerance. Adjust linear_tolerance until the two readings " *
          "agree (smaller to separate nearby roots, larger to merge a root " *
          "cluster), or compute over an exact field.")
end

# Shared head of both certified degrees: exact lift, gcd cancellation, and the
# rank guard against the default pipeline. Returns (f_reduced, g_reduced, d).
# Pointedness needs no separate guard: cancellation drops numerator and
# denominator degrees by the same count in both readings, so once the ranks
# agree the pointedness verdicts agree too.
function _certified_reduced_pair(f::HC.Expression, g::HC.Expression, var,
                                 linear_tolerance::Real)
    K, _ = _gaussian_rationals()
    Kx, _ = polynomial_ring(K, :x)
    fK = _exact_qi_poly(Kx, f, var)
    gK = _exact_qi_poly(Kx, g, var)
    common = gcd(fK, gK)
    # Oscar's polynomial gcd over a field is monic, so cancelling it leaves
    # lc(g_reduced) = lc(g) — the exact scalar formulas below then match the
    # float path, which multiplies by the ORIGINAL leading coefficient of g
    # after cancellation (see `lcg` in the default bodies).
    isone(leading_coefficient(common)) ||
        error("internal error: expected a monic gcd")
    f_reduced = divexact(fK, common)
    g_reduced = divexact(gK, common)
    d = degree(f_reduced)
    r1_tol, _ = _remove_common_approx(_poly_roots(f, var), _poly_roots(g, var),
                                      linear_tolerance)
    d == length(r1_tol) ||
        _certified_divergence_error("the rank (degree of the reduced numerator)",
                                    d, length(r1_tol), linear_tolerance)
    (f_reduced, g_reduced, d)
end

function _certified_global_unstable(f::HC.Expression, g::HC.Expression,
                                    linear_tolerance::Real)
    var = _the_variable(f, g)
    K, _ = _gaussian_rationals()
    f_reduced, g_reduced, d = _certified_reduced_pair(f, g, var, linear_tolerance)
    d > degree(g_reduced) ||
        error("the rational function is not pointed after reduction")
    monic_f = divexact(f_reduced, leading_coefficient(f_reduced))
    # the float path's scalar (-1)^(d(d-1)/2) ∏_a lc(g) ∏_b (a-b), with a over
    # the roots of f_reduced and b over the roots of g_reduced, IS the
    # resultant of the monic numerator with g_reduced — computed here with no
    # root ever found
    exact_scalar = K((-1)^div(d^2 - d, 2)) * resultant(monic_f, g_reduced)
    scalar = _round_to_c64(exact_scalar)
    (iszero(scalar) || !isfinite(scalar)) && error(_CERT_UNREPRESENTABLE_ERR)
    GWuClass(Matrix{ComplexF64}(LinearAlgebra.I, d, d), scalar)
end

function _certified_local_unstable(f::HC.Expression, g::HC.Expression, r::Number,
                                   linear_tolerance::Real)
    var = _the_variable(f, g)
    K, _ = _gaussian_rationals()
    f_reduced, g_reduced, d = _certified_reduced_pair(f, g, var, linear_tolerance)
    d > degree(g_reduced) ||
        error("the rational function is not pointed after reduction")

    # r read through the same ComplexF64 conversion as the default path
    rc = ComplexF64(r)
    rK = _exact_qi(rc)
    monic_f = divexact(f_reduced, leading_coefficient(f_reduced))

    # is-root guard: the default reading tests |∏(rc - a)| over the numerical
    # roots against the tolerance (replicated verbatim), the exact reading
    # evaluates the reduced numerator at r
    is_zero_exact = iszero(evaluate(monic_f, rK))
    r1_tol, _ = _remove_common_approx(_poly_roots(f, var), _poly_roots(g, var),
                                      linear_tolerance)
    frr = prod((rc - a for a in r1_tol); init = one(ComplexF64))
    is_zero_tol = !(abs(frr) > linear_tolerance)
    is_zero_exact == is_zero_tol ||
        _certified_divergence_error("whether r is a zero after reduction",
                                    is_zero_exact, is_zero_tol, linear_tolerance)
    is_zero_exact ||
        error("the field element is not a zero of the function after reduction")

    # k = exact multiplicity of r in the reduced numerator; deflation stops at
    # the first quotient h with h(r) != 0, i.e. h = monic_f / (x - r)^k
    x = gen(parent(monic_f))
    k = 0
    h = monic_f
    while iszero(evaluate(h, rK))
        h = divexact(h, x - rK)
        k += 1
    end
    k_tol = count(a -> abs(a - rc) < linear_tolerance, r1_tol)
    k == k_tol ||
        _certified_divergence_error("the multiplicity of r", k, k_tol,
                                    linear_tolerance)

    # the float path's scalar (-1)^(k(k-1)/2) (grr / ld_denom)^k: grr is
    # g_reduced(r) (since lc(g_reduced) = lc(g)), and ld_denom — the product
    # over the roots of the numerator other than r — is h(r)
    g_at_r = evaluate(g_reduced, rK)
    h_at_r = evaluate(h, rK)
    exact_scalar = K((-1)^div(k^2 - k, 2)) * divexact(g_at_r, h_at_r)^k
    scalar = _round_to_c64(exact_scalar)
    (iszero(scalar) || !isfinite(scalar)) && error(_CERT_UNREPRESENTABLE_ERR)
    GWuClass(Matrix{ComplexF64}(LinearAlgebra.I, k, k), scalar)
end

# ---------- unstable degrees over CC ------------------------------------------------

function global_unstable_A1_degree(f::HC.Expression, g::HC.Expression;
                                   linear_tolerance::Real = 1e-6,
                                   certify::Bool = false)
    linear_tolerance > 0 || error("linearTolerance must be a positive number")
    certify && return _certified_global_unstable(f, g, linear_tolerance)
    var = _the_variable(f, g)
    lcg = _lead_coefficient(g, var)
    r1, r2 = _remove_common_approx(_poly_roots(f, var), _poly_roots(g, var),
                                   linear_tolerance)
    d = length(r1)
    d > length(r2) || error("the rational function is not pointed after reduction")
    # (-1)^((d^2-d)/2) * resultant(fr, gr), with fr monic from the roots
    scalar = ComplexF64((-1)^div(d^2 - d, 2))
    for a in r1
        scalar *= lcg * prod((a - b for b in r2); init = one(ComplexF64))
    end
    GWuClass(Matrix{ComplexF64}(LinearAlgebra.I, d, d), scalar)
end

function local_unstable_A1_degree(f::HC.Expression, g::HC.Expression, r::Number;
                                  linear_tolerance::Real = 1e-6,
                                  certify::Bool = false)
    linear_tolerance > 0 || error("linearTolerance must be a positive number")
    certify && return _certified_local_unstable(f, g, r, linear_tolerance)
    var = _the_variable(f, g)
    lcg = _lead_coefficient(g, var)
    r1, r2 = _remove_common_approx(_poly_roots(f, var), _poly_roots(g, var),
                                   linear_tolerance)
    d = length(r1)
    d > length(r2) || error("the rational function is not pointed after reduction")
    rc = ComplexF64(r)
    frr = prod((rc - a for a in r1); init = one(ComplexF64))
    abs(frr) > linear_tolerance &&
        error("the field element is not a zero of the function after reduction")
    k = count(a -> abs(a - rc) < linear_tolerance, r1)
    roots_not_r = filter(a -> abs(a - rc) >= linear_tolerance, r1)
    ld_denom = prod((rc - a for a in roots_not_r); init = one(ComplexF64))
    grr = lcg * prod((rc - b for b in r2); init = one(ComplexF64))
    scalar = ComplexF64((-1)^div(k^2 - k, 2)) * (grr / ld_denom)^k
    GWuClass(Matrix{ComplexF64}(LinearAlgebra.I, k, k), scalar)
end
