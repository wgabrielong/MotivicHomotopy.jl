# Local and global anisotropic dimensions over a number field, following
# Koprowski-Czogala, *Computing with quadratic forms over number fields*,
# J. Symbolic Comput. 89 (2018) 129-145 [KC18], Algorithm 9 — the algorithm
# KR cite as [20]. Mirrors the QQ port in anisotropic_dimension.jl
# (KC18 Algorithms 6 and 8 specialised there to QQ_p).

# Is the diagonal form D totally hyperbolic over K_P? (KC18 Algorithm 6 shape.)
function _is_hyperbolic_at(K::AbsSimpleNumField, D::Vector, P)
    r = length(D)
    isodd(r) && return false
    sgn = iseven(div(r * (r - 1), 2)) ? 1 : -1
    d = K(sgn) * prod(D; init = one(K))          # signed discriminant of D
    is_local_square(d, P) || return false
    m = div(r, 2)
    _hasse_witt_from_diagonal(D, P) ==
        _hilbert_symbol(K(-1), K(-1), P)^div(m * (m - 1), 2)
end

# Anisotropic dimension of the diagonal form D over the completion K_P
# (KC18 Algorithm 8 shape). Always 0..4; the dim_a = 4 case is the unique
# anisotropic quaternary <1, -u, -pi, u*pi> over K_P, whose discriminant is
# trivial — hence "even rank, not hyperbolic, square discriminant -> 4".
function _local_aniso_dim_finite(K::AbsSimpleNumField, D::Vector, P)
    r = length(D)
    if iseven(r)
        _is_hyperbolic_at(K, D, P) && return 0
        sgn = iseven(div(r * (r - 1), 2)) ? 1 : -1
        d = K(sgn) * prod(D; init = one(K))
        return is_local_square(d, P) ? 4 : 2
    else
        # rank odd: adjoin the complement c making the extension a
        # discriminant-trivial even form; dim_a is 1 or 3 by parity
        sgn = iseven(div(r * (r + 1), 2)) ? 1 : -1
        c = K(sgn) * prod(D; init = one(K))
        return _is_hyperbolic_at(K, vcat(D, [c]), P) ? 1 : 3
    end
end

# The anisotropic dimension of beta over the completion at the place v: at
# most 1 at a complex place (rank parity), the absolute signature at a real
# place, and the [KC18] rank/discriminant/Hasse computation at a finite
# place, always <= 4. Finite places read the fixed diagonalisation held by
# the class's _KRForm.
function _local_anisotropic_dimension(beta::GWClass, sigma::InfPlc)
    _check_nf_class_and_field(beta, number_field(sigma))
    is_real(sigma) ? abs(form_signature(beta, sigma)) : form_rank(beta) % 2
end

function _local_anisotropic_dimension(beta::GWClass, P::AbsNumFieldOrderIdeal)
    _check_nf_class_and_field(beta, number_field(order(P)))
    is_prime(P) || error("second argument must be a prime ideal")
    kr = _kr_form(beta)
    _local_aniso_dim_finite(kr.K, kr.diag, P)
end

"""
    local_anisotropic_dimension(beta::GWClass, v) -> Int

The anisotropic dimension of a form over a number field ``K`` after
completion at the place `v` of ``K`` — the number-field counterpart of
[`anisotropic_dimension_qqp`](@ref) (which is ``\\mathbb{Q}_p``-only). `v` is
either a prime ideal of the maximal order of ``K`` or an infinite place
(`InfPlc`): at a finite place the result is at most 4 (the unique
anisotropic quaternary over a local field), at a real place it is the
absolute value of the signature, and at a complex place the rank mod 2.
The global [`anisotropic_dimension`](@ref) is realised at one of these
places, up to the even-rank/non-square-discriminant correction documented
there.

# Examples
```julia-repl
julia> K, s = quadratic_field(-7);

julia> P = prime_decomposition(maximal_order(K), 2)[1][1];  # 2 splits: K_P = Q_2

julia> local_anisotropic_dimension(diagonal_form(K, (1, 1, 1, 1)), P)
4

julia> local_anisotropic_dimension(hyperbolic_form(K, 4), P)
0
```

# References
- [KC18] P. Koprowski and A. Czogała, *Computing with quadratic forms over number fields*, J. Symbolic Comput. 89 (2018), 129–145 — Algorithms 6 and 8, the local computation used here.
- [Lam05] T. Y. Lam, *Introduction to quadratic forms over fields*, GSM 67, AMS, 2005.

See also [`anisotropic_dimension`](@ref), [`anisotropic_dimension_qqp`](@ref),
[`local_hilbert_symbol`](@ref).
"""
local_anisotropic_dimension(beta::GWClass, P::AbsNumFieldOrderIdeal) =
    _local_anisotropic_dimension(beta, P)
local_anisotropic_dimension(beta::GWClass, sigma::InfPlc) =
    _local_anisotropic_dimension(beta, sigma)

# Global anisotropic dimension over K ([KC18] Algorithm 9 shape): the
# maximum of the local anisotropic dimensions. Off real ∪ S the fixed
# diagonal is a unit form at an odd prime, where dim_a is 0 for even rank
# with locally-square discriminant, 1 for odd rank, and 2 for even rank
# with locally-non-square discriminant — so one GLOBAL correction is needed
# beyond the max over real ∪ S: for even rank whose signed discriminant is
# not a square in K, some place off S sees dim_a = 2 (Chebotarev), and that
# place need not lie in S. The blind spot is real, not theoretical: over
# K = Q(sqrt 34), eps = 35 + 6 sqrt 34 is a totally positive non-square
# unit that IS a local square at the unique dyadic prime, so <1, -eps> is
# locally isotropic at every place of real ∪ S yet globally anisotropic.
# (Memo §10 flagged the naive max-over-S reading as an assumption to check;
# this is the result of checking it.) Complex places (dim_a = rank mod 2)
# never exceed a finite place of the same parity and are not enumerated.
function _anisotropic_dimension_nf(beta::GWClass)
    kr = _kr_form(beta)
    dims = Int[form_rank(beta) % 2]
    if iseven(form_rank(beta)) && !is_square(signed_discriminant(beta))
        push!(dims, 2)
    end
    for sigma in real_places(kr.K)
        push!(dims, abs(form_signature(beta, sigma)))
    end
    for P in kr.S
        push!(dims, _local_aniso_dim_finite(kr.K, kr.diag, P))
    end
    maximum(dims)
end

_anisotropic_dimension(::NumberFieldKind, A::MatElem) =
    _anisotropic_dimension_nf(GWClass(_nondegenerate_part_diagonal(A)))
