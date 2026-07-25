# anisotropic_part, gw_witt_index, sum_decomposition, sum_decomposition_string,
# diagonal_class. For every form q the joint statement
#
#     q  ~Hecke~  anisotropic_part(q) ⟂ H^w      (global isometry)
#     anisotropic_part(q) is anisotropic          (oracle_isotropic_global)
#
# determines the Witt index and the isometry class of the anisotropic part
# uniquely.

function suite_decomposition!(res::DiffResults, fields::Vector{DiffField}, rng::DetRNG)
    for f in fields
        K = f.K
        forms = forms_for_tier(f, form_pool(f, rng), 36, 22, 12)
        for D in forms
            r = length(D)
            beta = to_class(K, D)

            part = anisotropic_part(beta)
            part_diag = diagonal_entries(part)
            w = gw_witt_index(beta)

            check!(res, "anisotropic_part",
                   length(part_diag) == anisotropic_dimension(beta),
                   "rank of part = anisotropic dimension: " * case_str(f.name, D))
            hyperbolic_pad = repeat([K(1), K(-1)], w)
            check!(res, "anisotropic_part",
                   oracle_isometric(K, [K(a) for a in D],
                                    vcat([K(a) for a in part_diag], hyperbolic_pad)),
                   "q ~ part ⟂ H^w (Hecke isometry): " * case_str(f.name, D) *
                   " part=$(part_diag) w=$w")
            check!(res, "anisotropic_part",
                   !oracle_isotropic_global(K, [K(a) for a in part_diag]),
                   "part is anisotropic (oracle): " * case_str(f.name, D) *
                   " part=$(part_diag)")

            check!(res, "gw_witt_index", w == div(r - length(part_diag), 2),
                   "w = (rank - dim_a)/2: " * case_str(f.name, D))
            # Adding hyperbolic planes raises the Witt index by exactly one each.
            check!(res, "gw_witt_index",
                   gw_witt_index(to_class(K, vcat(D, [K(1), K(-1)]))) == w + 1,
                   "w(q ⟂ H) = w(q) + 1: " * case_str(f.name, D))

            # sum_decomposition reassembles to the isometry class of q; the
            # string mirrors it: a "wH" (or "H") head and rank - 2w entries.
            dec = sum_decomposition(beta)
            dec_diag = diagonal_entries(dec)
            check!(res, "sum_decomposition",
                   length(dec_diag) == r &&
                       oracle_isometric(K, [K(a) for a in D], [K(a) for a in dec_diag]),
                   "reassembles to q (Hecke isometry): " * case_str(f.name, D))
            s = sum_decomposition_string(beta)
            expected_head = w == 0 ? "" : w == 1 ? "H" : "$(w)H"
            head_ok = w == 0 ? !occursin("H", s) : startswith(s, expected_head)
            pieces = count("<", s)
            check!(res, "sum_decomposition_string",
                   r == 0 ? s == "empty form" : (head_ok && pieces == r - 2 * w),
                   "string matches decomposition: " * case_str(f.name, D) * " got \"$s\" w=$w")

            # diagonal_class: isometric to q, and literally stable on a fresh
            # class (cache coherence: two calls, one object, same matrix).
            fresh = to_class(K, D)
            dc1 = diagonal_class(fresh)
            check!(res, "diagonal_class",
                   oracle_isometric(K, [K(a) for a in D],
                                    [K(a) for a in diagonal_entries(dc1)]),
                   "isometric to q (Hecke): " * case_str(f.name, D))
            check!(res, "diagonal_class", gw_matrix(diagonal_class(fresh)) == gw_matrix(dc1),
                   "cache-coherent repeat call: " * case_str(f.name, D))
        end
    end
end
