# The unstable layer: GWuClass validation (is_square(det * c); scalars
# compare exactly under isomorphism), stable_part, and divisorial_sum, whose
# scalar formula prod_i c_i * prod_{i<j} (r_i - r_j)^(2 m_i m_j) is
# recomputed here from scratch.

function suite_gwu!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K
        E = element_pool(f)
        u = first(a for a in E if !is_square(a))
        small = [D for D in form_pool(f, rng) if 1 <= length(D) <= 3]
        bases = small[1:min(9, length(small))]

        for D in bases
            M = oracle_gram(K, D)
            d = det(M)

            # Validation: det itself and det times squares construct; det
            # times a certified nonsquare, and a zero scalar, must not.
            check!(res, "GWuClass", _gwu_constructs(M, d),
                   "accepts c = det: " * case_str(f.name, D))
            check!(res, "GWuClass", _gwu_constructs(M, d * K(4)) &&
                       _gwu_constructs(M, d * det_pick!(rng, E)^2),
                   "accepts c = det times squares: " * case_str(f.name, D))
            check!(res, "GWuClass", !_gwu_constructs(M, d * u),
                   "rejects c = det times nonsquare u=$u: " * case_str(f.name, D))
            check!(res, "GWuClass", !_gwu_constructs(M, zero(K)),
                   "rejects c = 0: " * case_str(f.name, D))

            # Exact-scalar isomorphism semantics: same class under a
            # permutation with the same scalar; NOT isomorphic when the
            # scalar moves inside its square class.
            alpha = GWuClass(M, d)
            check!(res, "is_isomorphic_form",
                   is_isomorphic_form(alpha, GWuClass(oracle_gram(K, reverse(D)), d)),
                   "GWu: permutation with equal scalar: " * case_str(f.name, D))
            check!(res, "is_isomorphic_form",
                   !is_isomorphic_form(alpha, GWuClass(M, d * K(4))),
                   "GWu: scalar moved inside its square class is NOT equal: " *
                   case_str(f.name, D))

            # stable_part forgets the scalar and keeps the class.
            sp = stable_part(alpha)
            check!(res, "stable_part", gw_matrix(sp) == M,
                   "matrix kept: " * case_str(f.name, D))
            check!(res, "stable_part",
                   oracle_isometric(K, [K(a) for a in D],
                                    [K(a) for a in diagonal_entries(sp)]),
                   "class kept (Hecke isometry): " * case_str(f.name, D))
        end

        # divisorial_sum against a from-scratch recomputation of the
        # local-to-global formula, over pairs and triples with distinct roots.
        roots_pool = AbsSimpleNumFieldElem[K(0), K(1), K(-2), f.gen, f.gen + 3]
        for _ in 1:10
            k = 2 + det_int!(rng, 2) - 1              # 2 or 3 summands
            Ds = [det_pick!(rng, bases) for _ in 1:k]
            roots = _distinct_prefix(roots_pool, k)
            length(roots) == k || continue
            classes = [GWuClass(oracle_gram(K, D)) for D in Ds]

            ds = divisorial_sum(classes, roots)

            expected_scalar = prod(gw_scalar(c) for c in classes)
            for i in 1:k, j in (i + 1):k
                expected_scalar *=
                    (roots[i] - roots[j])^(2 * length(Ds[i]) * length(Ds[j]))
            end
            check!(res, "divisorial_sum", gw_scalar(ds) == expected_scalar,
                   "scalar formula recomputed: field=$(f.name) diags=$(Ds) roots=$(roots)")
            check!(res, "divisorial_sum",
                   gw_matrix(ds) == oracle_gram(K, reduce(vcat, Ds)),
                   "matrix is the block sum: field=$(f.name) diags=$(Ds)")
            check!(res, "stable_part",
                   form_rank(stable_part(ds)) == sum(length(D) for D in Ds),
                   "stable part of the sum has the total rank: field=$(f.name)")
        end
    end
end

"Does GWuClass accept the pair without throwing?"
function _gwu_constructs(M::MatElem, c)::Bool
    try
        GWuClass(M, c)
        true
    catch
        false
    end
end

"The first `k` pairwise distinct elements of `pool` (fewer if the pool repeats)."
function _distinct_prefix(pool::Vector, k::Int)
    out = eltype(pool)[]
    for r in pool
        any(x -> x == r, out) || push!(out, r)
        length(out) == k && break
    end
    out
end
