# hasse_witt_invariant, signed_discriminant, form_signature, form_signatures,
# relevant_primes — each against a source Hecke computes independently of the
# package (or, for signatures, against certified numerics).

function suite_invariants!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K
        forms = form_pool(f, rng)

        # hasse_witt_invariant vs Hecke's hasse_invariant, at every dyadic
        # prime and the odd control primes.
        for D in forms, P in finite_test_places(K)
            isempty(D) && continue
            beta = to_class(K, D)
            check!(res, "hasse_witt_invariant",
                   hasse_witt_invariant(beta, P) == oracle_hasse(K, D, P),
                   "vs Hecke hasse_invariant: " * case_str(f.name, D, P))
        end

        # signed_discriminant vs Hecke's discriminant (signed convention,
        # gate-verified), compared as square classes.
        for D in forms
            beta = to_class(K, D)
            ours = signed_discriminant(beta)
            hecke = oracle_signed_disc(K, D)
            check!(res, "signed_discriminant", is_square(K(ours) * hecke),
                   "square class vs Hecke discriminant: " * case_str(f.name, D) *
                   " ours=$ours hecke=$hecke")
        end
        # Multiplicativity: disc(q1 ⟂ q2) = disc(q1) disc(q2) (-1)^(r1 r2)
        # as square classes ([Lam05] II.5.11 for the signed discriminant).
        for _ in 1:5
            D1, D2 = det_pick!(rng, forms), det_pick!(rng, forms)
            (isempty(D1) || isempty(D2)) && continue
            d12 = signed_discriminant(to_class(K, vcat(D1, D2)))
            d1 = signed_discriminant(to_class(K, D1))
            d2 = signed_discriminant(to_class(K, D2))
            expected = K(d1) * K(d2) * K(-1)^(length(D1) * length(D2))
            check!(res, "signed_discriminant", is_square(K(d12) * expected),
                   "multiplicativity: " * case_str(f.name, vcat(D1, D2)))
        end

        # form_signatures vs certified-numeric signatures, as multisets (the
        # numeric embedding order need not match real_places order).
        embs = [e for e in complex_embeddings(K) if is_real(e)]
        for D in forms
            beta = to_class(K, D)
            sigs = form_signatures(beta)
            numeric = Int[numeric_signature(e, [K(a) for a in D]) for e in embs]
            check!(res, "form_signatures", same_multiset(sigs, numeric),
                   "vs certified numerics: " * case_str(f.name, D) *
                   " ours=$sigs numeric=$numeric")
            # form_signature at each real place is the matching entry of
            # form_signatures, and is additive in the first argument.
            for (i, sigma) in enumerate(real_places(K))
                check!(res, "form_signature",
                       form_signature(beta, sigma) == sigs[i],
                       "per-place consistency: " * case_str(f.name, D, sigma))
            end
        end
        for _ in 1:4
            D1, D2 = det_pick!(rng, forms), det_pick!(rng, forms)
            q12 = to_class(K, vcat(D1, D2))
            for sigma in real_places(K)
                check!(res, "form_signature",
                       form_signature(q12, sigma) ==
                           form_signature(to_class(K, D1), sigma) +
                           form_signature(to_class(K, D2), sigma),
                       "additivity: " * case_str(f.name, vcat(D1, D2), sigma))
            end
        end

        # relevant_primes: at any prime NOT in the returned set, every entry
        # has even valuation and Hecke's Hasse invariant is trivial — i.e.
        # the form is locally a unit form (times squares) with no invariant
        # content. Probes the primes over 3, 5, 7 that the set omits.
        for D in forms
            isempty(D) && continue
            beta = to_class(K, D)
            R = relevant_primes(beta)
            for p in (3, 5, 7)
                for (P, _) in prime_decomposition(maximal_order(K), p)
                    any(Q -> Q == P, R) && continue
                    even_vals = all(a -> iseven(valuation(K(a), P)), D)
                    check!(res, "relevant_primes",
                           even_vals && oracle_hasse(K, D, P) == 1,
                           "omitted prime carries no invariant content: " *
                           case_str(f.name, D, P))
                    break   # one omitted prime per rational p suffices
                end
            end
        end
    end
end
