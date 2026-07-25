# is_isomorphic_form against Hecke's is_isometric: permutations, square
# scalings, and dense congruence transforms (isometric by construction), and
# one-entry twists (ours and Hecke's must agree either way).

function suite_isomorphism!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K
        E = element_pool(f)
        bases = [D for D in form_pool(f, rng) if 1 <= length(D) <= 5]
        for D in bases
            n = length(D)
            beta = to_class(K, D)

            # Permutation: reverse is enough to move every entry for n >= 2.
            perm = reverse(D)
            check!(res, "is_isomorphic_form",
                   is_isomorphic_form(beta, to_class(K, perm)) &&
                       oracle_isometric(K, D, perm),
                   "permutation invariance (ours AND Hecke): " * case_str(f.name, D))

            # Square scaling of every entry.
            c = det_pick!(rng, E)
            scaled = [c^2 * a for a in D]
            check!(res, "is_isomorphic_form",
                   is_isomorphic_form(beta, to_class(K, scaled)) &&
                       oracle_isometric(K, D, scaled),
                   "square-scaling invariance: " * case_str(f.name, D) * " c=$c")

            # Dense congruence transform T' A T, det T != 0: same class, but
            # the Gram matrices are wildly different.
            A = oracle_gram(K, D)
            T = zero_matrix(K, n, n)
            while iszero(det(T))
                for i in 1:n, j in 1:n
                    T[i, j] = K(det_int!(rng, 5) - 3)
                end
            end
            B = transpose(T) * A * T
            check!(res, "is_isomorphic_form",
                   is_isomorphic_form(beta, GWClass(B)) && oracle_isometric_gram(A, B),
                   "congruence transform: " * case_str(f.name, D))

            # Twist one entry by a pool element: verdict decided by Hecke,
            # ours must agree either way.
            u = det_pick!(rng, E)
            twisted = copy(D)
            twisted[det_int!(rng, n)] *= u
            check!(res, "is_isomorphic_form",
                   is_isomorphic_form(beta, to_class(K, twisted)) ==
                       oracle_isometric(K, D, twisted),
                   "twist agreement with Hecke: " * case_str(f.name, D) * " u=$u twisted=$twisted")

            # Rank mismatch is never isomorphic.
            check!(res, "is_isomorphic_form",
                   !is_isomorphic_form(beta, to_class(K, vcat(D, [K(1)]))),
                   "rank mismatch: " * case_str(f.name, D))
        end
    end
end
