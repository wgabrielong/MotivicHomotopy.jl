# The semiring layer: gw_direct_sum, gw_tensor_product, integer multiples,
# tensor powers, diagonal_form / hyperbolic_form / pfister_form, form_rank.
# Cross-function identities: the Hasse sum formula
# c(q1 ⟂ q2) = c(q1) c(q2) (det q1, det q2)_P ([Lam05] V.3.18) and the
# Pfister dichotomy (isotropic implies hyperbolic, [Lam05] X.1.7).

function suite_ring_ops!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K
        E = element_pool(f)
        pairs_pool = [D for D in form_pool(f, rng) if 1 <= length(D) <= 4]

        for _ in 1:16
            D1 = det_pick!(rng, pairs_pool)
            D2 = det_pick!(rng, pairs_pool)
            q1, q2 = to_class(K, D1), to_class(K, D2)

            # Direct sum: literal block structure, rank additivity, and the
            # Hasse sum formula at the dyadic and odd control places.
            s = gw_direct_sum(q1, q2)
            # gw_matrix, not diagonal_entries: the latter size-reduces its
            # representatives over number fields, so only the Gram matrix
            # itself is a literal witness of the block structure.
            check!(res, "gw_direct_sum",
                   gw_matrix(s) == oracle_gram(K, vcat(D1, D2)),
                   "block-diagonal literal: " * case_str(f.name, vcat(D1, D2)))
            check!(res, "form_rank", form_rank(s) == length(D1) + length(D2),
                   "rank additivity: " * case_str(f.name, vcat(D1, D2)))
            det1 = prod(D1; init = one(K))
            det2 = prod(D2; init = one(K))
            for P in finite_test_places(K)
                check!(res, "gw_direct_sum",
                       hasse_witt_invariant(s, P) ==
                           hasse_witt_invariant(q1, P) * hasse_witt_invariant(q2, P) *
                           local_hilbert_symbol(det1, det2, P),
                       "Hasse sum formula: " * case_str(f.name, vcat(D1, D2), P))
            end

            # Tensor product: literal Kronecker diagonal, rank
            # multiplicativity, per-place signature multiplicativity.
            t = gw_tensor_product(q1, q2)
            kron = [a * b for a in D1 for b in D2]
            check!(res, "gw_tensor_product", gw_matrix(t) == oracle_gram(K, kron),
                   "Kronecker literal: " * case_str(f.name, vcat(D1, D2)))
            check!(res, "form_rank", form_rank(t) == length(D1) * length(D2),
                   "rank multiplicativity: " * case_str(f.name, vcat(D1, D2)))
            for sigma in real_places(K)
                check!(res, "gw_tensor_product",
                       form_signature(t, sigma) ==
                           form_signature(q1, sigma) * form_signature(q2, sigma),
                       "signature multiplicativity: " * case_str(f.name, vcat(D1, D2), sigma))
            end

            # Integer multiple and tensor power against their definitions.
            for n in (2, 3)
                check!(res, "* (integer multiple)",
                       gw_matrix(n * q1) ==
                           gw_matrix(reduce(gw_direct_sum, fill(q1, n))),
                       "n*q = q ⟂ ... ⟂ q (n=$n): " * case_str(f.name, D1))
            end
            if length(D1) <= 2
                for n in (2, 3)
                    check!(res, "^ (tensor power)",
                           gw_matrix(q1^n) ==
                               gw_matrix(reduce(gw_tensor_product, fill(q1, n))),
                           "q^n = q ⊗ ... ⊗ q (n=$n): " * case_str(f.name, D1))
                end
            end
        end

        # diagonal_form: the public constructor produces exactly the diagonal
        # Gram matrix, and permuting the entries lands in the same class.
        for _ in 1:8
            D = det_pick!(rng, pairs_pool)
            q = diagonal_form(K, Tuple(D))
            check!(res, "diagonal_form", gw_matrix(q) == oracle_gram(K, D),
                   "literal Gram: " * case_str(f.name, D))
            check!(res, "diagonal_form",
                   oracle_isometric(K, D, reverse(D)) &&
                       is_isomorphic_form(q, diagonal_form(K, Tuple(reverse(D)))),
                   "permutation lands in the same class: " * case_str(f.name, D))
        end

        # hyperbolic_form: H^k is totally isotropic — Witt index k by the
        # package AND isotropic by the oracle; and padding any form by H^k
        # raises its Witt index by exactly k.
        for k in (1, 2, 3)
            h = hyperbolic_form(K, 2 * k)
            hd = diagonal_entries(h)
            check!(res, "hyperbolic_form",
                   gw_witt_index(h) == k && anisotropic_dimension(h) == 0 &&
                       oracle_isotropic_global(K, [K(a) for a in hd]),
                   "H^$k is split: field=$(f.name)")
            check!(res, "hyperbolic_form",
                   sum_decomposition_string(hyperbolic_form(K, 2 * k)) ==
                       (k == 1 ? "H" : "$(k)H"),
                   "H^$k decomposition string: field=$(f.name)")
        end
        for _ in 1:6
            D = det_pick!(rng, pairs_pool)
            k = det_int!(rng, 2)
            padded = gw_direct_sum(to_class(K, D), hyperbolic_form(K, 2 * k))
            check!(res, "hyperbolic_form",
                   gw_witt_index(padded) == gw_witt_index(to_class(K, D)) + k,
                   "padding by H^$k shifts the Witt index: " * case_str(f.name, D))
        end

        # pfister_form: <<a, b>> = <1,-a> ⊗ <1,-b> literally, and the Pfister
        # dichotomy isotropic => hyperbolic, decided by the oracle.
        for _ in 1:6
            a, b = det_pick!(rng, E), det_pick!(rng, E)
            p1 = pfister_form(K, (a,))
            check!(res, "pfister_form",
                   gw_matrix(p1) == gw_matrix(diagonal_form(K, (K(1), -a))),
                   "<<a>> = <1,-a>: field=$(f.name) a=$a")
            p2 = pfister_form(K, (a, b))
            manual = gw_tensor_product(diagonal_form(K, (K(1), -a)),
                                       diagonal_form(K, (K(1), -b)))
            check!(res, "pfister_form", gw_matrix(p2) == gw_matrix(manual),
                   "<<a,b>> literal tensor: field=$(f.name) a=$a b=$b")
            pd = [K(x) for x in diagonal_entries(p2)]
            if oracle_isotropic_global(K, pd)
                check!(res, "pfister_form", anisotropic_dimension(p2) == 0,
                       "Pfister dichotomy (isotropic => hyperbolic): field=$(f.name) a=$a b=$b")
            else
                check!(res, "pfister_form", anisotropic_dimension(p2) == 4,
                       "Pfister dichotomy (anisotropic stays rank 4): field=$(f.name) a=$a b=$b")
            end
        end
    end
end
