# Port of M2 Code/GrothendieckWittClasses.m2.
#
# Representation: M2's GrothendieckWittClass is a HashTable {matrix, cache}.
# Here: a struct holding the Gram matrix plus a Dict cache (M2 caches only the
# diagonal class; the Dict reproduces that mutable cache on an immutable class).
#
# M2's RR_53 / CC_53 are represented by Matrix{Float64} / Matrix{ComplexF64};
# exact fields and étale algebras use Oscar matrices (MatElem).

const FloatGram = Union{Matrix{Float64}, Matrix{ComplexF64}}

"""
    GWClass(M)

The isomorphism class of the nondegenerate symmetric bilinear form with Gram
matrix `M`, as an element of the Grothendieck–Witt ring ``\\text{GW}(k)`` of a field
(or finite étale algebra over a field) of characteristic not 2.

Given a basis ``e_1, …, e_n`` of a ``k``-vector space ``V``, a symmetric
bilinear form ``β : V × V → k`` is encoded by its Gram matrix
``(β(e_i, e_j))_{i,j}``; a change of basis replaces the Gram matrix by a
congruent one, so a symmetric matrix determines the form up to congruence.
The constructor checks that `M` is symmetric, nondegenerate, and defined over
a supported coefficient ring, and errors otherwise.

Supported Gram matrices:
- an Oscar matrix (`MatElem`) over `QQ`, a finite field of odd characteristic,
  or a finite étale algebra (a zero-dimensional `MPolyQuoRing` with
  nondegenerate trace form);
- a plain Julia matrix of real or complex numbers, stored as `Matrix{Float64}`
  or `Matrix{ComplexF64}` — these play the role of forms over ``\\mathbb{R}`` and ``\\mathbb{C}``.

Equality `==` compares Gram matrices literally (same base ring, same entries).
Mathematical equality in ``\\text{GW}(k)`` is tested with
[`is_isomorphic_form`](@ref).

For float Gram matrices, `GWClass(M; certify = true)` decides nondegeneracy
by the exact determinant of `M` read entrywise as exact rationals, instead of
the default floating-point determinant; on exact base fields the keyword is a
no-op.

The Gram matrix, coefficient algebra, and base field are recovered with
[`gw_matrix`](@ref), [`gw_algebra`](@ref), and [`gw_base_field`](@ref); a
diagonal representative with [`diagonal_class`](@ref). Further invariants:
[`form_rank`](@ref), [`form_signature`](@ref),
[`integral_discriminant`](@ref), [`hasse_witt_invariant`](@ref),
[`anisotropic_dimension`](@ref), [`anisotropic_part`](@ref),
[`sum_decomposition`](@ref), [`is_isotropic_form`](@ref),
[`is_anisotropic_form`](@ref).

# Examples
```julia-repl
julia> beta = GWClass(QQ[2 1; 1 3])
[2   1]
[1   3]

julia> gw_matrix(beta)
[2   1]
[1   3]

julia> gw_base_field(beta)
Rational field

julia> GWClass([0.0 1.0; 1.0 0.0])    # a form over the real numbers
2×2 Matrix{Float64}:
 0.0  1.0
 1.0  0.0
```
"""
struct GWClass{T}
    matrix::T
    cache::Dict{Symbol, Any}
end

# --- well-definedness (M2: isWellDefinedGW) ------------------------------------

_ring_char(R) = characteristic(R)
_ring_char(R::MPolyQuoRing) = characteristic(coefficient_ring(base_ring(R)))

function _is_well_defined_gw(M::MatElem; certify::Bool = false)
    _is_symmetric_square(M) || return false
    iszero(det(M)) && return false
    R = base_ring(M)
    (R isa Field || is_finite_etale_algebra(R)) || return false
    _ring_char(R) != 2 || return false
    true
end

function _is_well_defined_gw(M::FloatGram; certify::Bool = false)
    # certify=true: nondegeneracy decided by the exact determinant over
    # QQ / QQ(i), beyond M2; default branch unchanged. Non-finite entries are
    # refused loudly before any other test — NaN breaks even the symmetry
    # comparison (NaN != NaN), and a silent `false` would misreport the reason.
    if certify
        all(isfinite, M) || error(_CERT_NONFINITE_ERR)
        _is_symmetric_square(M) || return false
        return _certified_det_nonzero(M)
    end
    _is_symmetric_square(M) || return false
    !iszero(LinearAlgebra.det(M))
end

const _GW_ERR = "makeGWClass called on a matrix that does not represent a " *
    "nondegenerate symmetric bilinear form over a field of characteristic not 2"

# --- constructors (M2: makeGWClass) --------------------------------------------

function GWClass(M::MatElem; certify::Bool = false)
    _is_well_defined_gw(M) || error(_GW_ERR)
    GWClass{typeof(M)}(M, Dict{Symbol, Any}())
end

function GWClass(M::AbstractMatrix{<:Real}; certify::Bool = false)
    Mf = Matrix{Float64}(M)
    _is_well_defined_gw(Mf; certify = certify) || error(_GW_ERR)
    GWClass{Matrix{Float64}}(Mf, Dict{Symbol, Any}())
end

function GWClass(M::AbstractMatrix{<:Complex}; certify::Bool = false)
    Mf = Matrix{ComplexF64}(M)
    _is_well_defined_gw(Mf; certify = certify) || error(_GW_ERR)
    GWClass{Matrix{ComplexF64}}(Mf, Dict{Symbol, Any}())
end

# --- accessors (M2: getMatrix / getAlgebra / getBaseField) ----------------------

"""
    gw_matrix(beta)

The Gram matrix of a [`GWClass`](@ref) or [`GWuClass`](@ref): a symmetric
matrix over the coefficient algebra of the class (`MatElem` for exact fields
and étale algebras, `Matrix{Float64}` / `Matrix{ComplexF64}` over ``\\mathbb{R}`` / ``\\mathbb{C}``).

# Examples
```julia-repl
julia> beta = GWClass(QQ[2 1; 1 3]);

julia> gw_matrix(beta)
[2   1]
[1   3]
```

See also [`gw_scalar`](@ref), [`gw_algebra`](@ref), [`gw_base_field`](@ref).
"""
gw_matrix(beta::GWClass) = beta.matrix

_algebra_of(M::MatElem) = base_ring(M)
_algebra_of(::Matrix{Float64}) = Float64
_algebra_of(::Matrix{ComplexF64}) = ComplexF64

function _base_field_of(M::MatElem)
    R = base_ring(M)
    (R isa QQField || R isa FinField) && return R
    R isa Field && return R
    if R isa MPolyQuoRing
        # M2: requires the zero ideal to be prime, then applies toField.
        # Julia has no toField; we return the (domain) quotient ring itself.
        is_prime(ideal(R, [zero(R)])) ||
            error("the Grothendieck-Witt class is not defined over a field")
        return R
    end
    error("the Grothendieck-Witt class is not defined over a field")
end
_base_field_of(::Matrix{Float64}) = Float64
_base_field_of(::Matrix{ComplexF64}) = ComplexF64

"""
    gw_algebra(beta)

The coefficient algebra over which a [`GWClass`](@ref) or [`GWuClass`](@ref)
is defined — a field or a finite étale algebra over a field. Over ``\\mathbb{R}`` / ``\\mathbb{C}``
(float-backed classes) this returns the type `Float64` / `ComplexF64`.

# Examples
```julia-repl
julia> gw_algebra(GWClass(QQ[2 1; 1 3]))
Rational field
```

See also [`gw_base_field`](@ref), [`gw_matrix`](@ref).
"""
gw_algebra(beta::GWClass) = _algebra_of(beta.matrix)

"""
    gw_base_field(beta)

The base field of a [`GWClass`](@ref) or [`GWuClass`](@ref), when the class is
defined over a field. For a class over an étale algebra this checks that the
zero ideal is prime (i.e. the algebra is a field) and errors otherwise; for a
class that already lives over `QQ` or a finite field it returns that field.
Over ``\\mathbb{R}`` / ``\\mathbb{C}`` (float-backed classes) it returns the type `Float64` /
`ComplexF64`.

# Examples
```julia-repl
julia> gw_base_field(GWClass(QQ[2 1; 1 3]))
Rational field
```

See also [`gw_algebra`](@ref).
"""
gw_base_field(beta::GWClass) = _base_field_of(beta.matrix)

# --- equality (M2 `===`: literal Gram matrix, cache ignored) --------------------

_gram_equal(A::MatElem, B::MatElem) =
    base_ring(A) === base_ring(B) && size(A) == size(B) && A == B
_gram_equal(A::Matrix{T}, B::Matrix{T}) where {T <: Union{Float64, ComplexF64}} =
    size(A) == size(B) && A == B
_gram_equal(::Any, ::Any) = false

Base.:(==)(a::GWClass, b::GWClass) = _gram_equal(a.matrix, b.matrix)
Base.hash(a::GWClass, h::UInt) = hash(a.matrix, hash(:GWClass, h))

# M2 prints the underlying matrix (net GrothendieckWittClass := net getMatrix).
# Over a number field the matrix entries do not identify the field (the same
# polynomial expression means different things in different fields), so the
# field is appended; over every other backend the output is the matrix alone,
# exactly as before.
function Base.show(io::IO, beta::GWClass)
    show(io, beta.matrix)
    _show_field_suffix(io, beta.matrix)
end
function Base.show(io::IO, mime::MIME"text/plain", beta::GWClass)
    show(io, mime, beta.matrix)
    _show_field_suffix(io, beta.matrix)
end

_show_field_suffix(::IO, ::Any) = nothing
function _show_field_suffix(io::IO, M::MatElem)
    R = base_ring(M)
    R isa AbsSimpleNumField && print(io, "\nover ", R)
    nothing
end

# --- operations (M2: addGW / multiplyGW) ----------------------------------------

const _FIELD_MISMATCH_ERR = "these classes have different underlying fields"

# M2 allows two Galois fields of the same order; map the second into the first.
function _align_second(Mb::MatElem, Mg::MatElem)
    Kb, Kg = base_ring(Mb), base_ring(Mg)
    if Kb isa FinField && Kg isa FinField
        order(Kb) == order(Kg) || error(_FIELD_MISMATCH_ERR)
        Kb === Kg && return Mg
        # prime fields transfer through ZZ (M2: substitute(-, Kb)); a same-order
        # transfer between distinct non-prime fields would need an embedding.
        return map_entries(x -> Kb(lift(ZZ, x)), Mg)
    end
    Kb === Kg || error(_FIELD_MISMATCH_ERR)
    Mg
end

"""
    gw_direct_sum(beta, gamma)

The direct (block) sum of two classes over the same base field — addition in
the Grothendieck–Witt ring. For two [`GWClass`](@ref) inputs the result has
the block-diagonal Gram matrix; for two [`GWuClass`](@ref) inputs the Gram
matrices are block-summed and the scalars multiplied. `beta + gamma` is an
alias.

# Examples
```julia-repl
julia> beta1 = GWClass(QQ[1 2; 2 3]);

julia> beta2 = GWClass(QQ[3 4; 4 5]);

julia> gw_direct_sum(beta1, beta2)
[1   2   0   0]
[2   3   0   0]
[0   0   3   4]
[0   0   4   5]

julia> gw_direct_sum(GWuClass(QQ[2 1; 1 2]), GWuClass(QQ[1 2; 2 6]))
([2 1 0 0; 1 2 0 0; 0 0 1 2; 0 0 2 6], 6)
```

See also [`gw_tensor_product`](@ref), [`divisorial_sum`](@ref).
"""
function gw_direct_sum(beta::GWClass{<:MatElem}, gamma::GWClass{<:MatElem})
    Mg = _align_second(beta.matrix, gamma.matrix)
    GWClass(block_diagonal_matrix([beta.matrix, Mg]))
end

"""
    gw_tensor_product(beta, gamma)

The tensor product of two [`GWClass`](@ref)es over the same base field —
multiplication in the Grothendieck–Witt ring. The resulting Gram matrix is
the Kronecker product of the two Gram matrices. `beta * gamma` is an alias.

# Examples
```julia-repl
julia> beta1 = GWClass(QQ[1 2; 2 3]);

julia> beta2 = GWClass(QQ[3 4; 4 5]);

julia> gw_tensor_product(beta1, beta2)
[3    4    6    8]
[4    5    8   10]
[6    8    9   12]
[8   10   12   15]
```

See also [`gw_direct_sum`](@ref).
"""
function gw_tensor_product(beta::GWClass{<:MatElem}, gamma::GWClass{<:MatElem})
    Mg = _align_second(beta.matrix, gamma.matrix)
    GWClass(kronecker_product(beta.matrix, Mg))
end

function _float_block_diag(A::Matrix{T}, B::Matrix{T}) where {T}
    n1, n2 = size(A, 1), size(B, 1)
    C = zeros(T, n1 + n2, n1 + n2)
    C[1:n1, 1:n1] = A
    C[(n1 + 1):end, (n1 + 1):end] = B
    C
end

gw_direct_sum(beta::GWClass{Matrix{T}}, gamma::GWClass{Matrix{T}}) where {T} =
    GWClass(_float_block_diag(beta.matrix, gamma.matrix))
gw_tensor_product(beta::GWClass{Matrix{T}}, gamma::GWClass{Matrix{T}}) where {T} =
    GWClass(LinearAlgebra.kron(beta.matrix, gamma.matrix))

gw_direct_sum(::GWClass, ::GWClass) = error(_FIELD_MISMATCH_ERR)
gw_tensor_product(::GWClass, ::GWClass) = error(_FIELD_MISMATCH_ERR)

# Optional operator aliases; the named functions are canonical.
Base.:+(beta::GWClass, gamma::GWClass) = gw_direct_sum(beta, gamma)
Base.:*(beta::GWClass, gamma::GWClass) = gw_tensor_product(beta, gamma)

# --- semiring scalar multiple and power ----------------------------------------
# GW classes here form a SEMIRING: sums and products exist, formal
# differences do not (see the module documentation). Hence n >= 0 only.

# the rank-zero form (additive identity) and <1> (multiplicative identity)
# over the class's coefficient algebra
_gw_zero_class(alg) = GWClass(_diag_gram(alg, Any[]))
_gw_one_class(alg::Ring) = GWClass(_diag_gram(alg, [one(alg)]))
_gw_one_class(::Type{T}) where {T <: Union{Float64, ComplexF64}} =
    GWClass(_diag_gram(T, [1]))

const _SEMIRING_ERR = "GW classes form a semiring: no negation or formal " *
    "differences exist, so the integer factor must be nonnegative (see the " *
    "MotivicHomotopy module documentation)"

"""
    n * beta
    beta * n

The `n`-fold orthogonal sum ``β ⊕ ⋯ ⊕ β`` of a [`GWClass`](@ref), for an
integer ``n ≥ 0``; `0 * beta` is the rank-zero form. Negative `n` is an
error: GW classes form a semiring here — there is no negation (see the
module documentation).

# Examples
```julia-repl
julia> 2 * diagonal_form(QQ, (1, 3))
[1   0   0   0]
[0   3   0   0]
[0   0   1   0]
[0   0   0   3]

julia> form_rank(0 * diagonal_form(QQ, (1, 3)))
0
```

# References
- [Lam05] T. Y. Lam, *Introduction to quadratic forms over fields*, GSM 67, AMS, 2005 — Chapter II for the semiring of isometry classes inside the Witt–Grothendieck ring.

See also [`gw_direct_sum`](@ref) (and `beta^n` for the tensor power).
"""
function Base.:*(n::Integer, beta::GWClass)
    n >= 0 || error(_SEMIRING_ERR)
    out = _gw_zero_class(gw_algebra(beta))
    for _ in 1:n
        out = gw_direct_sum(out, beta)
    end
    out
end
Base.:*(beta::GWClass, n::Integer) = n * beta

"""
    beta ^ n

The `n`-fold tensor product ``β ⊗ ⋯ ⊗ β`` of a [`GWClass`](@ref), for an
integer ``n ≥ 0``; `beta^0` is the rank-one form ``⟨1⟩``. Negative `n` is an
error: GW classes form a semiring here and are not invertible (see the
module documentation).

# Examples
```julia-repl
julia> beta = diagonal_form(QQ, (1, 3));

julia> form_rank(beta^3)
8

julia> beta^0
[1]
```

# References
- [Lam05] T. Y. Lam, *Introduction to quadratic forms over fields*, GSM 67, AMS, 2005 — Chapter II for the tensor product of bilinear forms.

See also [`gw_tensor_product`](@ref), [`pfister_form`](@ref).
"""
function Base.:^(beta::GWClass, n::Integer)
    n >= 0 || error(_SEMIRING_ERR)
    out = _gw_one_class(gw_algebra(beta))
    for _ in 1:n
        out = gw_tensor_product(out, beta)
    end
    out
end

# Julia lowers a LITERAL beta^-1 to inv(beta) via literal_pow; route every
# literal exponent through the method above so negative literals raise the
# semiring error instead of a MethodError for inv.
Base.literal_pow(::typeof(^), beta::GWClass, ::Val{p}) where {p} = beta^p

# --- base change GW(Q) -> GW(K) --------------------------------------------------

"""
    gw_base_change(K::AbsSimpleNumField, beta::GWClass) -> GWClass

The image of a class over ``\\mathbb{Q}`` under the base-change map
``\\text{GW}(\\mathbb{Q}) → \\text{GW}(K)`` to a number field: the same Gram
matrix, reinterpreted entrywise over ``K``. Base change is a ring
homomorphism — it commutes with [`gw_direct_sum`](@ref) and
[`gw_tensor_product`](@ref) — and preserves the rank and the square class of
the discriminant (a rational square is a square in ``K``; the converse can
fail, so distinct classes over ``\\mathbb{Q}`` may collide in ``\\text{GW}(K)``:
``⟨1⟩ ≠ ⟨2⟩`` over ``\\mathbb{Q}`` but ``⟨1⟩ ≅ ⟨2⟩`` over ``\\mathbb{Q}(\\sqrt 2)``).

Only ``\\mathbb{Q} → K`` is provided: between two number fields the map
depends on a choice of embedding, which the Gram matrix does not carry.

# Examples
```julia-repl
julia> K, s = quadratic_field(2);

julia> gw_base_change(K, diagonal_form(QQ, (1, -2)))
[1    0]
[0   -2]
over Real quadratic field defined by x^2 - 2

julia> is_isomorphic_form(gw_base_change(K, diagonal_form(QQ, (1, -2))),
                          hyperbolic_form(K))
true
```

# References
- [Lam05] T. Y. Lam, *Introduction to quadratic forms over fields*, GSM 67, AMS, 2005 — VII.2 for scalar extension of quadratic forms.

See also [`GWClass`](@ref), [`gw_direct_sum`](@ref).
"""
function gw_base_change(K::AbsSimpleNumField, beta::GWClass)
    M = gw_matrix(beta)
    (M isa MatElem && base_ring(M) isa QQField) ||
        error("gw_base_change is the map GW(Q) -> GW(K): the class must be over QQ")
    GWClass(map_entries(K, M))
end
