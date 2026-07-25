# Emits a results file in the exact line schema oracle_ext_*.m would print,
# computed with Hecke primitives (hilbert_symbol, is_square, sign at real
# places). It validates the parse/compare pipeline and is NOT Magma output
# (the file it writes says so in its header). The package's finite-place
# Hilbert symbol delegates to Hecke, so the hs@*/hasse@* keys here are a
# plumbing check only.

# Hecke-primitive counterparts of the Magma prelude helpers.
_sim_hasse(K, D, P) = prod(Int[Hecke.hilbert_symbol(K(D[i]), K(D[j]), P)
                               for i in 1:length(D) for j in (i + 1):length(D)];
                           init = 1)
_sim_sig(K, D, v) = sum(Int[sign(K(a), v) for a in D]; init = 0)
_sim_hsr(K, a, b, v) = (sign(K(a), v) < 0 && sign(K(b), v) < 0) ? -1 : 1
_sim_sd(K, D) = K((-1)^div(length(D) * (length(D) - 1), 2)) *
                prod(K.(D); init = one(K))

function _sim_lines(r::ExtRecord, K::AbsSimpleNumField)::Vector{Tuple{String, String}}
    OK = maximal_order(K)
    decomp(p) = [P for (P, _) in prime_decomposition(OK, p)]
    rp = real_places(K)
    out = Tuple{String, String}[]
    put(key, val) = push!(out, (key, string(val)))
    pay = r.payload

    if r.kind == :sym
        a, b = pay[:a], pay[:b]
        for p in pay[:ps], P in decomp(p)
            put("hs@$p", Hecke.hilbert_symbol(a, b, P))
        end
        pay[:use_real] && for v in rp
            put("hs@real", _sim_hsr(K, a, b, v))
        end

    elseif r.kind == :form
        D = pay[:D]
        put("sdisc_sq", is_square(_sim_sd(K, D)))
        put("sd_match", is_square(_sim_sd(K, D) * pay[:sd]))
        for v in rp
            put("sig", _sim_sig(K, D, v))
        end
        for p in pay[:support], P in decomp(p)
            put("hasse@$p", _sim_hasse(K, D, P))
        end
        length(D) == 2 && put("nps", is_square(K(-D[1] * D[2])))
        if haskey(pay, :ap)
            ap, dp = pay[:ap], pay[:padded]
            put("ap_disc_eq", is_square(_sim_sd(K, D) * _sim_sd(K, dp)))
            put("ap_sig_eq", all(_sim_sig(K, D, v) == _sim_sig(K, dp, v) for v in rp))
            for p in pay[:support]
                put("ap_hasse_eq@$p",
                    all(_sim_hasse(K, D, P) == _sim_hasse(K, dp, P) for P in decomp(p)))
            end
            length(ap) == 2 && put("ap_nps", is_square(K(-ap[1] * ap[2])))
            for v in rp
                put("ap_sig", _sim_sig(K, ap, v))
            end
            if haskey(pay, :ap_comp)
                apx = vcat(ap, [pay[:ap_comp]])
                for p in pay[:support], P in decomp(p)
                    put("ap_dim3@$p", _sim_hasse(K, apx, P) !=
                        Hecke.hilbert_symbol(K(-1), K(-1), P)^pay[:ap_e])
                end
            end
        end
        ex = pay[:extra]
        if haskey(ex, :sum)
            D1, D2 = ex[:sum]
            put("sum_disc_ok", is_square(prod(D; init = one(K)) *
                prod(D1; init = one(K)) * prod(D2; init = one(K))))
            for p in pay[:support]
                put("sum_hasse_ok@$p", all(
                    _sim_hasse(K, D, P) == _sim_hasse(K, D1, P) * _sim_hasse(K, D2, P) *
                        Hecke.hilbert_symbol(prod(D1; init = one(K)),
                                             prod(D2; init = one(K)), P)
                    for P in decomp(p)))
            end
        end
        if haskey(ex, :tensor)
            D1, D2 = ex[:tensor]
            put("tensor_disc_ok", is_square(prod(D; init = one(K)) *
                prod(D1; init = one(K))^length(D2) * prod(D2; init = one(K))^length(D1)))
            put("tensor_sig_ok", all(_sim_sig(K, D, v) ==
                _sim_sig(K, D1, v) * _sim_sig(K, D2, v) for v in rp))
        end
        if haskey(ex, :mult)
            D1, n = ex[:mult]
            put("mult_disc_ok", is_square(prod(D; init = one(K)) *
                prod(D1; init = one(K))^n))
            put("mult_sig_ok", all(_sim_sig(K, D, v) == n * _sim_sig(K, D1, v) for v in rp))
        end
        if haskey(ex, :pow)
            D1, n = ex[:pow]
            put("pow_sig_ok", all(_sim_sig(K, D, v) == _sim_sig(K, D1, v)^n for v in rp))
        end
        haskey(ex, :bc_sigq) &&
            put("bc_sig_ok", all(_sim_sig(K, D, v) == ex[:bc_sigq] for v in rp))
        if haskey(pay, :rp_omit)
            p = pay[:rp_omit]
            put("rp_clean@$p", all(_sim_hasse(K, D, P) == 1 for P in decomp(p)))
        end

    elseif r.kind == :ldim
        dx = vcat(pay[:D], [pay[:comp]])
        for p in pay[:support], P in decomp(p)
            put("ldim3@$p", _sim_hasse(K, dx, P) !=
                Hecke.hilbert_symbol(K(-1), K(-1), P)^pay[:e])
        end

    elseif r.kind == :pair
        D1, D2 = pay[:D1], pay[:D2]
        put("disc_ratio_sq", is_square(prod(D1; init = one(K)) * prod(D2; init = one(K))))
        put("sig_eq", all(_sim_sig(K, D1, v) == _sim_sig(K, D2, v) for v in rp))
        for p in pay[:support]
            put("hasse_eq@$p",
                all(_sim_hasse(K, D1, P) == _sim_hasse(K, D2, P) for P in decomp(p)))
        end

    elseif r.kind == :dense
        A = pay[:A]
        put("det_ratio_sq", is_square(det(A) * prod(pay[:ourD]; init = one(K))))
        put("ents_ratio_sq", is_square(det(A) * prod(pay[:ents]; init = one(K))))
        put("dcls_ratio_sq", is_square(det(A) * prod(pay[:dcls]; init = one(K))))
        # a real Magma run may print true/false here
        put("og_class_eq", "skip")

    elseif r.kind == :gwu
        put("scale_sq", is_square(prod(pay[:D]; init = one(K)) * pay[:c]))

    elseif r.kind == :dsum
        Ds, roots = pay[:Ds], pay[:roots]
        expected = prod(prod(Dq; init = one(K)) for Dq in Ds)
        for a in 1:length(Ds), b in (a + 1):length(Ds)
            expected *= (roots[a] - roots[b])^(2 * length(Ds[a]) * length(Ds[b]))
        end
        put("ds_scalar_ok", pay[:scalar] == expected)
    else
        error("unknown record kind $(r.kind)")
    end
    out
end

"""
    simulate_results(recs, fields, path)

Write the Hecke-simulated results file for the whole corpus to `path`.
"""
function simulate_results(recs::Vector{ExtRecord}, fields::Vector{DiffField},
                          path::String)
    kmap = Dict(f.name => f.K for f in fields)
    open(path, "w") do io
        println(io, "# HECKE-SIMULATED results — NOT Magma output. Generated by")
        println(io, "# differential/magma-oracle/generate_all.jl to validate the")
        println(io, "# pipeline; replace with the print-out of a real Magma run of the")
        println(io, "# oracle_ext_*.m chunks as soon as one is available.")
        for r in recs
            for (key, val) in _sim_lines(r, kmap[r.fname])
                println(io, "D|$(r.id)|$key|$val")
            end
        end
        println(io, "D|META|simulated|true")
    end
    path
end
