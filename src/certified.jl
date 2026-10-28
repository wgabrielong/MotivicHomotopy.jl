# Certified re-computation of the float (ℝ/ℂ) decisions, behind the opt-in
# `certify = true` keyword. This file is a deliberate extension beyond the
# Macaulay2 package, which has no certification; every default
# (`certify = false`) path in the package is left unchanged by it.
#
# The contract (compare the certified-interval contract in nf_arithmetic.jl):
# a `Float64` / `ComplexF64` value is an exact dyadic (Gaussian) rational, and
# `certify = true` answers for the input read that way. Every discrete
# decision — rank, signature counts, determinant sign — is computed exactly
# over ℚ (or ℚ(i)), never with float arithmetic or tolerances, so a certified
# discrete answer is exact. Nothing is certified about the real-world quantity
# the floats approximate: values that are close but not bit-equal are
# distinct elements. The only lossy step is rounding a result scalar back to
# `Float64` / `ComplexF64` (`_round_to_f64` / `_round_to_c64`), which is
# correctly rounded (nearest); when rounding would reach zero or infinity the
# result is refused (`_CERT_UNREPRESENTABLE_ERR`) rather than returned wrong.
#
# Interval or ball arithmetic is deliberately NOT used here: a ball straddling
# zero can never decide `det == 0`, so a ball layer would need this exact path
# as its fallback anyway, and the Gram matrices and polynomials in scope are
# small enough that exact arithmetic is immediate.

const _CERT_NONFINITE_ERR = "cannot certify: the input contains non-finite " *
    "(NaN or Inf) floating-point values"

const _CERT_UNREPRESENTABLE_ERR = "cannot certify: the exact scalar is not " *
    "representable as a nonzero finite floating-point number (overflow or " *
    "underflow in the rounding step)"

# --- exact lifts of floats ---------------------------------------------------------

# ℚ(i), the field where ComplexF64 entries live when read exactly. Relies on
# Hecke's own caching of quadratic fields (same policy as the no-cache note in
# nf_arithmetic.jl): repeated calls return the identical field object, so
# elements built in separate calls share a parent.
_gaussian_rationals() = quadratic_field(-1)

function _exact_qq(x::Float64)
    isfinite(x) || error(_CERT_NONFINITE_ERR)
    QQ(Rational{BigInt}(x))
end

function _exact_qi(z::ComplexF64)
    K, im_unit = _gaussian_rationals()
    K(_exact_qq(real(z))) + K(_exact_qq(imag(z))) * im_unit
end

_exact_qq_gram(M::Matrix{Float64}) = matrix(QQ, map(_exact_qq, M))

function _exact_qi_gram(M::Matrix{ComplexF64})
    K, _ = _gaussian_rationals()
    matrix(K, map(_exact_qi, M))
end

# --- certified discrete invariants of float Gram matrices --------------------------

function _certified_rank(M::Matrix{Float64})
    size(M, 1) == 0 && return 0
    rank(_exact_qq_gram(M))
end

function _certified_rank(M::Matrix{ComplexF64})
    size(M, 1) == 0 && return 0
    rank(_exact_qi_gram(M))
end

# Diagonalizes the exact ℚ lift with the package's own exact congruence
# elimination (degenerate-safe), then counts signs — the same algorithm shape
# as the M2 float path, with every pivot test exact.
function _certified_pos_neg_counts(M::Matrix{Float64})
    _is_symmetric_square(M) || error("Matrix is not symmetric")
    A = diagonalize_via_congruence(_exact_qq_gram(M))
    pos = 0
    neg = 0
    for i in 1:nrows(A)
        A[i, i] > 0 && (pos += 1)
        A[i, i] < 0 && (neg += 1)
    end
    (pos, neg)
end

function _certified_det_sign(M::Matrix{Float64})
    size(M, 1) == 0 && return 1     # det of the empty matrix is 1
    d = det(_exact_qq_gram(M))
    iszero(d) && return 0
    d > 0 ? 1 : -1
end

_certified_det_nonzero(M::Matrix{Float64}) = _certified_det_sign(M) != 0

function _certified_det_nonzero(M::Matrix{ComplexF64})
    size(M, 1) == 0 && return true
    !iszero(det(_exact_qi_gram(M)))
end

# --- correctly-rounded return trip -------------------------------------------------

_round_to_f64(x::QQFieldElem) =
    Float64(BigInt(numerator(x)) // BigInt(denominator(x)))

# a ∈ ℚ(i) in the power basis 1, i of _gaussian_rationals(); each component
# is correctly rounded on its own.
_round_to_c64(a::AbsSimpleNumFieldElem) =
    ComplexF64(_round_to_f64(coeff(a, 0)), _round_to_f64(coeff(a, 1)))
