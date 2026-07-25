# Deterministic case generators. The field matrix is the package's own
# (test/nf_test_fields.jl, included by run.jl before this file); every
# element and form pool is a fixed function of that matrix and the DetRNG
# seed.

struct DiffField
    name::String
    K::AbsSimpleNumField
    gen::AbsSimpleNumFieldElem
    # cost tier 1 = cheap (class number 1, degree <= 2, full case load),
    # 2 = moderate (h = 2 quadratics), 3 = expensive (higher degree)
    tier::Int
end

"""
    diff_fields() -> Vector{DiffField}

The differential field matrix: every field of test/nf_test_fields.jl except
none — the full matrix, with cost tiers steering how many cases the expensive
suites draw per field.
"""
function diff_fields()
    [
        DiffField("Q", KQ, KQ(1), 1),
        DiffField("Q(sqrt-7)", K7, s7, 1),
        DiffField("Q(sqrt2)", K2, s2, 1),
        DiffField("Q(sqrt5)", Kr5, sr5, 1),
        DiffField("Q(sqrt-2)", K2i, s2i, 1),
        DiffField("Q(sqrt-5)", K5, s5, 2),
        DiffField("Q(sqrt-21)", K21, s21, 2),
        DiffField("Q(sqrt34)", K34, s34, 2),
        DiffField("Q(sqrt10)", Ks10, s10, 2),
        DiffField("x^3-3x-1", Kc, ac, 3),
        DiffField("x^3-2", K32, a32, 3),
        DiffField("x^4-4x^2+1", Kq, aq, 3),
        DiffField("Q(zeta16)+", K162, w16, 3),
    ]
end

# Field-specific elements, keyed by the DiffField name: the Q(sqrt34)
# blind-spot witness eps = 35 + 6 sqrt34 and the norm-2 element 6 + sqrt34;
# the [KR21] section 6 diagonal over Q(sqrt-7); w^2 - 2 = ±sqrt2 in the
# wildly ramified dyadic corner of Q(zeta16)+.
function special_elements(f::DiffField)
    g = f.gen
    K = f.K
    f.name == "Q(sqrt34)" && return [K(35 + 6 * g), -K(35 + 6 * g), K(6 + g)]
    f.name == "Q(sqrt-7)" &&
        return [-3 - 9 * g, -2 - 6 * g, 1 - g, -6 + 4 * g, -3 + 2 * g, 4 - 4 * g]
    f.name == "Q(zeta16)+" && return [g^2 - 2, 2 - g^2, g^3 - 3 * g]
    f.name == "Q(sqrt2)" && return [K(2 + g), K(1 + g)]
    f.name == "Q(sqrt10)" && return [K(3 + g), K(-3 + g)]
    AbsSimpleNumFieldElem[]
end

"""
    element_pool(f::DiffField) -> Vector{AbsSimpleNumFieldElem}

Nonzero elements of `f.K`: small rationals, generator combinations, large
square content, denominators, even dyadic valuation, plus the
field-specific witnesses of [`special_elements`](@ref). Deterministic.
"""
function element_pool(f::DiffField)
    K, g = f.K, f.gen
    raw = AbsSimpleNumFieldElem[
        K(1), K(-1), K(2), K(-2), K(3), K(-3), K(5), K(6), K(7), K(-7),
        K(10), K(15), K(-30),
        # large square content
        K(4), K(-4), K(9), K(16), K(1024), K(-1024),
        # denominators, square and nonsquare
        K(1 // 2), K(-1 // 4), K(3 // 4), K(2 // 9),
        # million-scale primes
        K(1000003), K(-999983),
        # generator combinations
        g, -g, g + 1, g - 1, g + 2, 2 * g - 1, 3 * g + 2,
        # generator with square content and denominators
        4 * g, -16 * (g + 1), g * K(1 // 2), (g + 2) * K(1 // 9),
    ]
    append!(raw, special_elements(f))
    pool = AbsSimpleNumFieldElem[]
    for a in raw
        iszero(a) && continue
        any(b -> b == a, pool) || push!(pool, a)
    end
    pool
end

# How many deterministic-random diagonals each suite draws per field, by tier.
random_forms_per_field(tier::Int) = tier == 1 ? 30 : tier == 2 ? 20 : 12

"""
    form_pool(f::DiffField, rng::DetRNG) -> Vector{Vector{AbsSimpleNumFieldElem}}

Diagonals over `f.K`: rank 0 and 1, hyperbolic and scaled-hyperbolic planes,
definite blocks, even-dyadic-valuation blocks, the field's special witnesses
as binaries `<1, -w>`, the section 6 fixture on Q(sqrt-7), and
deterministic-random diagonals of rank 2..6 from [`element_pool`](@ref).
"""
function form_pool(f::DiffField, rng::DetRNG)
    K, g = f.K, f.gen
    E = element_pool(f)
    forms = Vector{AbsSimpleNumFieldElem}[
        AbsSimpleNumFieldElem[],                       # rank 0
        [K(1)], [K(-3)], [4 * g == 0 ? K(4) : 4 * g],  # rank 1
        [K(1), K(-1)],                                 # H
        [K(3), K(-3)], [g, -g],                        # scaled H
        [K(1), K(1)], [K(1), K(1), K(1)],              # definite over real fields
        [K(1), K(1), K(1), K(1)],                      # the Q_2 quaternary
        [K(4), K(8)], [K(4), K(-16), K(64)],           # even dyadic valuation
        [K(1 // 2), K(-2), K(3 // 4)],                 # denominators
        [K(1), K(-1), K(1), K(-1), K(2)],              # stacked hyperbolics
    ]
    for w in special_elements(f)
        push!(forms, [K(1), -w])
    end
    if f.name == "Q(sqrt-7)"
        push!(forms, [-3 - 9 * g, K(-1), -2 - 6 * g, 1 - g, -6 + 4 * g, -3 + 2 * g, 4 - 4 * g])
    end
    for _ in 1:random_forms_per_field(f.tier)
        r = det_int!(rng, 7) + 1                       # rank 2..8
        push!(forms, [det_pick!(rng, E) for _ in 1:r])
    end
    forms
end

"""
    dense_pool(f::DiffField, rng::DetRNG) -> Vector{MatElem}

Non-diagonal symmetric Gram matrices: T' D T for diagonal D from the form
pool and deterministic T with small entries and nonzero determinant. These
feed the diagonalization suite, where the input must NOT already be diagonal.
"""
function dense_pool(f::DiffField, rng::DetRNG)
    K = f.K
    out = MatElem[]
    for D in form_pool(f, DetRNG(rng.state))     # forked stream: same forms
        n = length(D)
        (2 <= n <= 5) || continue
        A = oracle_gram(K, D)
        T = zero_matrix(K, n, n)
        while iszero(det(T))
            for i in 1:n, j in 1:n
                T[i, j] = K(det_int!(rng, 5) - 3)  # entries in -2..2
            end
        end
        push!(out, transpose(T) * A * T)
    end
    out
end
