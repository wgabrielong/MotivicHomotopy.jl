# The Hecke oracle layer. Every intrinsic used here is convention-gated in
# gates.jl on hand-known values before any suite runs:
#
#   is_isometric(V, W)        global isometry over a number field
#   is_isometric(V, W, v)     local isometry at a place
#   hasse_invariant(V, P)     Lam's prod_{i<j} (a_i, a_j)_P convention
#   is_isotropic(V, v)        LOCAL isotropy at a finite or real place
#   discriminant(V)           the signed discriminant (-1)^(n(n-1)/2) det
#
# Hecke's GLOBAL is_isotropic(V) is broken over number fields (definite
# forms are reported isotropic) and is never called; oracle_isotropic_global
# assembles the answer from local verdicts (Hasse-Minkowski).

"The diagonal Gram matrix of `diag` over `K`."
function oracle_gram(K::AbsSimpleNumField, diag::Vector)
    n = length(diag)
    M = zero_matrix(K, n, n)
    for i in 1:n
        M[i, i] = K(diag[i])
    end
    M
end

oracle_space(K::AbsSimpleNumField, diag::Vector) = quadratic_space(K, oracle_gram(K, diag))

"Hecke's verdict on global isometry of two diagonal forms (rank-0 handled here)."
function oracle_isometric(K::AbsSimpleNumField, d1::Vector, d2::Vector)::Bool
    length(d1) == length(d2) || return false
    isempty(d1) && return true
    is_isometric(oracle_space(K, d1), oracle_space(K, d2))
end

"Hecke's verdict on global isometry of two Gram matrices."
function oracle_isometric_gram(M1::MatElem, M2::MatElem)::Bool
    nrows(M1) == nrows(M2) || return false
    nrows(M1) == 0 && return true
    is_isometric(quadratic_space(base_ring(M1), M1), quadratic_space(base_ring(M2), M2))
end

"Hecke's Hasse invariant (Lam convention) of a diagonal form at a prime ideal."
function oracle_hasse(K::AbsSimpleNumField, diag::Vector, P)::Int
    isempty(diag) && return 1
    hasse_invariant(oracle_space(K, diag), P)
end

"Hecke's LOCAL isotropy of a diagonal form at a finite or infinite place."
function oracle_isotropic_at(K::AbsSimpleNumField, diag::Vector, v)::Bool
    length(diag) <= 1 && return false
    is_isotropic(oracle_space(K, diag), v)
end

"Hecke's signed discriminant of a diagonal form (1 for rank 0)."
function oracle_signed_disc(K::AbsSimpleNumField, diag::Vector)
    isempty(diag) && return one(K)
    discriminant(oracle_space(K, diag))
end

# The places where a diagonal form can be locally nontrivial: the real
# places, the dyadic primes, and every prime where some entry has odd
# valuation. Off this set the form is, up to local squares, a unit form at an
# odd prime. (An entry of even valuation is a unit times a local square, so
# even-valuation support primes genuinely impose nothing for rank >= 3.)
function oracle_support_places(K::AbsSimpleNumField, diag::Vector)
    OK = maximal_order(K)
    finite = [P for (P, _) in prime_decomposition(OK, 2)]
    for a in diag
        for (P, e) in factor(K(a) * OK)
            isodd(e) && !any(Q -> Q == P, finite) && push!(finite, P)
        end
    end
    finite, real_places(K)
end

"""
    oracle_isotropic_global(K, diag) -> Bool

Global isotropy WITHOUT Hecke's broken global intrinsic:

- rank <= 1: anisotropic;
- rank 2: <a, b> is isotropic iff -ab is a square in K;
- rank >= 3: Hasse-Minkowski over the support places. A unit form of rank
  >= 3 at an odd prime is isotropic (Chevalley-Warning plus Hensel), so
  places off the support list cannot obstruct.
"""
function oracle_isotropic_global(K::AbsSimpleNumField, diag::Vector)::Bool
    length(diag) <= 1 && return false
    length(diag) == 2 && return is_square(K(-diag[1] * diag[2]))
    finite, reals = oracle_support_places(K, diag)
    all(v -> oracle_isotropic_at(K, diag, v), reals) &&
        all(P -> oracle_isotropic_at(K, diag, P), finite)
end
