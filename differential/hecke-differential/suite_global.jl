# anisotropic_dimension, is_isotropic_form, is_anisotropic_form, against
# oracle_isotropic_global plus parity, signature, stability, and scaling
# constraints.

function suite_global!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K
        E = element_pool(f)
        forms = form_pool(f, rng)
        for D in forms
            r = length(D)
            beta = to_class(K, D)
            d = anisotropic_dimension(beta)

            check!(res, "anisotropic_dimension", 0 <= d <= r && iseven(r - d),
                   "range/parity: " * case_str(f.name, D) * " d=$d")

            # The oracle decides isotropy; d == r must match "anisotropic".
            iso = oracle_isotropic_global(K, D)
            check!(res, "anisotropic_dimension", (d == r) == !iso,
                   "d = rank iff anisotropic (oracle): " * case_str(f.name, D) * " d=$d oracle_isotropic=$iso")
            check!(res, "is_isotropic_form", is_isotropic_form(beta) == iso,
                   "vs oracle: " * case_str(f.name, D))
            check!(res, "is_anisotropic_form", is_anisotropic_form(beta) == !iso,
                   "vs oracle: " * case_str(f.name, D))

            # Exact binary value: <a, b> has d = 0 iff -ab is a square, else 2.
            r == 2 && check!(res, "anisotropic_dimension",
                   d == (is_square(K(-D[1] * D[2])) ? 0 : 2),
                   "binary exact value: " * case_str(f.name, D) * " d=$d")

            # d dominates every real signature.
            isempty(real_places(K)) || check!(res, "anisotropic_dimension",
                   d >= maximum(abs.(form_signatures(beta))),
                   "d >= |signature|: " * case_str(f.name, D))

            # Hyperbolic stability and square-scaling invariance.
            check!(res, "anisotropic_dimension",
                   anisotropic_dimension(to_class(K, vcat(D, [K(1), K(-1)]))) == d,
                   "hyperbolic stability: " * case_str(f.name, D))
            if !isempty(D)
                c = det_pick!(rng, E)
                check!(res, "anisotropic_dimension",
                       anisotropic_dimension(to_class(K, [c * a for a in D])) == d,
                       "scaling invariance by c=$c: " * case_str(f.name, D))
            end
        end
    end
end
