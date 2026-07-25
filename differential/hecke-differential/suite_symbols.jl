# local_hilbert_symbol. The finite-place implementation delegates to Hecke's
# hilbert_symbol, so it is checked through identities rather than direct
# comparison: Hilbert reciprocity prod_v (a,b)_v = 1, the symbol axioms
# (values in {±1}, symmetry, bimultiplicativity, (a,-a) = 1), and agreement
# with hasse_invariant on binary forms.

function suite_symbols!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K
        E = element_pool(f)
        picks = [det_pick!(rng, E) for _ in 1:6]
        places = vcat(finite_test_places(K), real_places(K))

        for v in places, i in 1:length(picks), j in i:length(picks)
            a, b = picks[i], picks[j]
            s = local_hilbert_symbol(a, b, v)
            check!(res, "local_hilbert_symbol",
                   s in (-1, 1) && local_hilbert_symbol(b, a, v) == s,
                   "symbol axioms (±1, symmetry): " * case_str(f.name, [a, b], v))
            c = det_pick!(rng, picks)
            check!(res, "local_hilbert_symbol",
                   local_hilbert_symbol(a, b * c, v) ==
                       s * local_hilbert_symbol(a, c, v),
                   "bimultiplicativity: " * case_str(f.name, [a, b, c], v))
            check!(res, "local_hilbert_symbol",
                   local_hilbert_symbol(a, -a, v) == 1,
                   "(a, -a) = 1: " * case_str(f.name, [a], v))
        end

        # Hilbert reciprocity: the product over ALL places is 1, and only the
        # support places (oracle_support_places) can contribute -1.
        for _ in 1:25
            a, b = det_pick!(rng, E), det_pick!(rng, E)
            finite, reals = oracle_support_places(K, [a, b])
            total = prod(Int[local_hilbert_symbol(a, b, P) for P in finite]; init = 1) *
                    prod(Int[local_hilbert_symbol(a, b, v) for v in reals]; init = 1)
            check!(res, "local_hilbert_symbol", total == 1,
                   "Hilbert reciprocity: " * case_str(f.name, [a, b]) * " product=$total")
        end

        # Binary Hasse cross-check: c(<a,b>) = (a,b)_P by definition, and
        # Hecke's hasse_invariant reaches the value through its own quadratic
        # space machinery.
        for _ in 1:12
            a, b = det_pick!(rng, E), det_pick!(rng, E)
            for P in finite_test_places(K)
                check!(res, "local_hilbert_symbol",
                       local_hilbert_symbol(a, b, P) == oracle_hasse(K, [a, b], P),
                       "binary Hasse agreement: " * case_str(f.name, [a, b], P))
            end
        end
    end
end
