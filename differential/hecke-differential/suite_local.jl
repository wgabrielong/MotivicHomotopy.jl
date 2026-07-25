# local_anisotropic_dimension against Hecke's local isotropy (isotropic iff
# d < rank), plus range, parity, and hyperbolic-stability constraints.

function suite_local!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K
        forms = form_pool(f, rng)
        for D in forms
            isempty(D) && continue
            r = length(D)
            beta = to_class(K, D)
            beta_plus_H = to_class(K, vcat(D, [K(1), K(-1)]))

            for P in finite_test_places(K)
                d = local_anisotropic_dimension(beta, P)
                check!(res, "local_anisotropic_dimension",
                       0 <= d <= min(r, 4) && iseven(r - d),
                       "range/parity at a finite place: " * case_str(f.name, D, P) * " d=$d")
                r >= 2 && check!(res, "local_anisotropic_dimension",
                       (d < r) == oracle_isotropic_at(K, D, P),
                       "isotropy vs Hecke: " * case_str(f.name, D, P) * " d=$d")
                check!(res, "local_anisotropic_dimension",
                       local_anisotropic_dimension(beta_plus_H, P) == d,
                       "hyperbolic stability: " * case_str(f.name, D, P))
            end

            for sigma in real_places(K)
                d = local_anisotropic_dimension(beta, sigma)
                check!(res, "local_anisotropic_dimension",
                       d == abs(form_signature(beta, sigma)),
                       "real place = |signature|: " * case_str(f.name, D, sigma))
                r >= 2 && check!(res, "local_anisotropic_dimension",
                       (d < r) == oracle_isotropic_at(K, D, sigma),
                       "real-place isotropy vs Hecke: " * case_str(f.name, D, sigma))
            end
        end
    end
end
