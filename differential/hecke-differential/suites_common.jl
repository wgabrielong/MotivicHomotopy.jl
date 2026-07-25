# Helpers shared by the suites.

"The GWClass of a diagonal, built through the package's public constructor."
to_class(K::AbsSimpleNumField, D::Vector) = GWClass(oracle_gram(K, D))

"Reproduction string for a failure message."
case_str(fname::String, D::Vector) = "field=$(fname) diag=$(D)"
case_str(fname::String, D::Vector, v) = "field=$(fname) diag=$(D) place=$(v)"

"The first prime ideal of the maximal order of `K` above the rational prime `p`."
odd_prime_over(K::AbsSimpleNumField, p::Int) =
    prime_decomposition(maximal_order(K), p)[1][1]

"All dyadic primes of `K`."
dyadic_primes_of(K::AbsSimpleNumField) =
    [P for (P, _) in prime_decomposition(maximal_order(K), 2)]

"""
    finite_test_places(K) -> Vector

The finite places every place-indexed suite visits: all dyadic primes (the
structural hazard), one prime over 3 and one over 5 (odd control primes).
"""
finite_test_places(K::AbsSimpleNumField) =
    vcat(dyadic_primes_of(K), [odd_prime_over(K, 3), odd_prime_over(K, 5)])

"""
    numeric_signature(e, D; prec = 128) -> Int

The signature of the diagonal `D` under the real embedding `e`, decided
numerically with certified arb balls — the independent path against which
`form_signatures` (which uses exact sign computations) is compared. Raises
precision until every ball excludes zero, which terminates because the
entries are nonzero algebraic numbers.
"""
function numeric_signature(e, D::Vector; prec::Int = 128)
    total = 0
    for a in D
        p = prec
        while true
            x = real(e(a, p))
            if is_positive(x)
                total += 1
                break
            elseif is_negative(x)
                total -= 1
                break
            end
            p *= 2
            p > 2^16 && error("numeric_signature: ball never separated from zero for $a")
        end
    end
    total
end

"Multiset equality of two integer vectors."
same_multiset(a::Vector{Int}, b::Vector{Int}) = sort(a) == sort(b)

"""
    forms_for_tier(f::DiffField, forms::Vector, cap1, cap2, cap3) -> Vector

The tier-dependent prefix of a form pool: expensive suites draw fewer cases
from expensive fields. The pool order puts the curated pathological shapes
first, so truncation drops only deterministic-random cases.
"""
function forms_for_tier(f::DiffField, forms::Vector, cap1::Int, cap2::Int, cap3::Int)
    cap = f.tier == 1 ? cap1 : f.tier == 2 ? cap2 : cap3
    forms[1:min(cap, length(forms))]
end
