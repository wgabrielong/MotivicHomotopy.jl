# Comparison of a results file (Magma's print-out, or the Hecke simulation)
# against the package, record by record. Sources for the conclusions drawn:
#
#   - forms over a number field are classified by rank, signed discriminant,
#     Hasse invariants at all finite places, and signatures at all real
#     places ([Lam05] VI.3 over Q; [OM73] 66:4 over number fields) — the
#     :pair verdict and the Witt-equivalence certificate;
#   - a binary form <a, b> is isotropic iff -ab is a square — the rank-2
#     verdicts;
#   - an odd-rank form over a local field has anisotropic dimension 3 iff
#     its discriminant-trivialized extension is not totally hyperbolic
#     ([KC18] Algorithms 6 and 8) — the ldim3 and ap_dim3 keys;
#   - a globally anisotropic form has a locally anisotropic place
#     (Hasse-Minkowski), local anisotropy is impossible above rank 4, and a
#     real witness means |signature| = rank — the anisotropy certificates.
#
# Per-prime data is compared as multisets over the prime ideals above each
# stated rational prime.

"Parse a results file into id -> key -> values (order preserved per key)."
function parse_results(path::String)
    out = Dict{String, Dict{String, Vector{String}}}()
    for ln in eachline(path)
        startswith(ln, "D|") || continue
        parts = split(ln, '|')
        length(parts) == 4 || continue
        id, key, val = String(parts[2]), String(parts[3]), String(parts[4])
        push!(get!(get!(out, id, Dict{String, Vector{String}}()), key, String[]), val)
    end
    out
end

_istrue(v::Vector{String}) = !isempty(v) && all(==("true"), v)

# One comparison context: the shared tally plus this record's parsed lines.
struct CmpCtx
    res::DiffResults
    lines::Dict{String, Vector{String}}
    rec::ExtRecord
    K::AbsSimpleNumField
    missing_keys::Ref{Int}
end

lines_of(c::CmpCtx, key::String) = get(c.lines, key, String[])

"check! against the harness tally, with the record id in the context string."
function ccheck!(c::CmpCtx, func::String, ok::Bool, what::String)
    check!(c.res, func, ok, "$(c.rec.id) ($(c.rec.fname)): $what")
end

"Compare a multiset of ours against the parsed lines; skip (and count) if absent."
function multiset_check!(c::CmpCtx, func::String, key::String, ours::Vector{String},
                         what::String)
    L = lines_of(c, key)
    if isempty(L) && !isempty(ours)
        c.missing_keys[] += 1
        return
    end
    if sort(L) == sort(ours)
        for _ in L
            ccheck!(c, func, true, what)
        end
    else
        ccheck!(c, func, false, "$what: ours=$(sort(ours)) oracle=$(sort(L))")
    end
end

function _cmp_sym(c::CmpCtx)
    K, pay = c.K, c.rec.payload
    a, b = pay[:a], pay[:b]
    OK = maximal_order(K)
    for p in pay[:ps]
        ours = [string(local_hilbert_symbol(a, b, P))
                for (P, _) in prime_decomposition(OK, p)]
        multiset_check!(c, "local_hilbert_symbol", "hs@$p", ours, "symbols over $p")
    end
    if pay[:use_real]
        ours = [string(local_hilbert_symbol(a, b, v)) for v in real_places(K)]
        multiset_check!(c, "local_hilbert_symbol", "hs@real", ours, "symbols at real places")
    end
end

# The anisotropy certificate for the stored anisotropic part, from oracle
# lines alone. Mirrors _aniso_certifiable.
function _ap_certified(c::CmpCtx)::Bool
    pay = c.rec.payload
    ap = pay[:ap]
    r = length(ap)
    r <= 1 && return true
    if r == 2
        v = lines_of(c, "ap_nps")
        return !isempty(v) && all(==("false"), v)
    end
    sigs = [abs(parse(Int, s)) for s in lines_of(c, "ap_sig")]
    any(==(r), sigs) && return true
    if r == 3
        for p in pay[:support]
            any(==("true"), lines_of(c, "ap_dim3@$p")) && return true
        end
    end
    false
end

function _cmp_form(c::CmpCtx)
    K, pay = c.K, c.rec.payload
    D = pay[:D]
    beta = GWClass(oracle_gram(K, D))
    rp = real_places(K)

    # signed discriminant: Magma's square-class verdicts
    sq = lines_of(c, "sdisc_sq")
    isempty(sq) || ccheck!(c, "signed_discriminant",
        (sq[1] == "true") == is_square(K(signed_discriminant(beta))),
        "IsSquare(signed disc) agreement")
    m = lines_of(c, "sd_match")
    isempty(m) || ccheck!(c, "signed_discriminant", _istrue(m),
        "our representative lies in Magma's square class")

    # signatures, three consumers of the same oracle lines
    sig_lines = lines_of(c, "sig")
    ours_sigs = form_signatures(beta)
    sig_ok = false
    if !isempty(rp)
        if isempty(sig_lines)
            c.missing_keys[] += 1
        else
            sig_ok = sort(sig_lines) == sort(string.(ours_sigs))
            ccheck!(c, "form_signatures", sig_ok,
                "signature multiset: ours=$(sort(ours_sigs)) oracle=$(sort(sig_lines))")
            for (j, v) in enumerate(rp)
                ccheck!(c, "form_signature", sig_ok &&
                    form_signature(beta, v) == ours_sigs[j], "per-place signature")
                ccheck!(c, "local_anisotropic_dimension",
                    sig_ok && local_anisotropic_dimension(beta, v) == abs(ours_sigs[j]),
                    "real-place anisotropic dimension = |signature|")
            end
        end
    end

    # Hasse invariants as multisets per rational prime
    OK = maximal_order(K)
    hasse_ok = true
    for p in pay[:support]
        ours = [string(hasse_witt_invariant(beta, P))
                for (P, _) in prime_decomposition(OK, p)]
        L = lines_of(c, "hasse@$p")
        if isempty(L)
            c.missing_keys[] += 1
            hasse_ok = false
            continue
        end
        agree = sort(L) == sort(ours)
        hasse_ok &= agree
        if agree
            for _ in L
                ccheck!(c, "hasse_witt_invariant", true, "Hasse at $p")
            end
        else
            ccheck!(c, "hasse_witt_invariant", false,
                "Hasse at $p: ours=$(sort(ours)) oracle=$(sort(L))")
        end
    end
    std_ok = _istrue(m) && hasse_ok && (isempty(rp) || sig_ok)

    # rank 2: exact isotropy from IsSquare(-ab)
    nps = lines_of(c, "nps")
    if length(D) == 2 && !isempty(nps)
        iso = nps[1] == "true"
        ccheck!(c, "is_isotropic_form", is_isotropic_form(beta) == iso,
            "binary isotropy vs IsSquare(-ab)")
        ccheck!(c, "is_anisotropic_form", is_anisotropic_form(beta) == !iso,
            "binary anisotropy vs IsSquare(-ab)")
        ccheck!(c, "anisotropic_dimension",
            anisotropic_dimension(beta) == (iso ? 0 : 2), "binary exact dimension")
    end

    # the decomposition certificate
    if haskey(pay, :ap)
        w = pay[:w]
        wq_ok = _istrue(lines_of(c, "ap_disc_eq")) &&
                _istrue(lines_of(c, "ap_sig_eq")) &&
                all(_istrue(lines_of(c, "ap_hasse_eq@$p")) for p in pay[:support])
        if pay[:aniso_cert]
            full = wq_ok && _ap_certified(c)
            what = "q ~ part ⟂ H^$w with the part certified anisotropic"
            ccheck!(c, "anisotropic_part", full, what)
            ccheck!(c, "gw_witt_index", full, what)
            ccheck!(c, "sum_decomposition", full, what)
            ccheck!(c, "anisotropic_dimension",
                full && anisotropic_dimension(beta) == length(pay[:ap]), what)
            if length(D) != 2   # rank 2 already credited above via nps
                ccheck!(c, "is_isotropic_form",
                    full && is_isotropic_form(beta) == (w >= 1), what)
                ccheck!(c, "is_anisotropic_form",
                    full && is_anisotropic_form(beta) == (w == 0), what)
            end
        end
    end

    # composite relations and constructor tags
    ex = pay[:extra]
    if haskey(ex, :sum)
        ccheck!(c, "gw_direct_sum", _istrue(lines_of(c, "sum_disc_ok")) && std_ok,
            "sum discriminant relation")
        for p in pay[:support]
            ccheck!(c, "gw_direct_sum", _istrue(lines_of(c, "sum_hasse_ok@$p")),
                "Hasse sum formula at $p")
        end
    end
    if haskey(ex, :tensor)
        ccheck!(c, "gw_tensor_product",
            _istrue(lines_of(c, "tensor_disc_ok")) && std_ok, "tensor discriminant relation")
        ccheck!(c, "gw_tensor_product", _istrue(lines_of(c, "tensor_sig_ok")),
            "tensor signature relation")
    end
    if haskey(ex, :mult)
        ccheck!(c, "* (integer multiple)",
            _istrue(lines_of(c, "mult_disc_ok")) && std_ok, "multiple discriminant relation")
        ccheck!(c, "* (integer multiple)", _istrue(lines_of(c, "mult_sig_ok")),
            "multiple signature relation")
    end
    haskey(ex, :pow) && ccheck!(c, "^ (tensor power)",
        _istrue(lines_of(c, "pow_sig_ok")) && std_ok, "power signature relation")
    haskey(ex, :bc_sigq) && ccheck!(c, "gw_base_change",
        _istrue(lines_of(c, "bc_sig_ok")) && std_ok,
        "all K-signatures equal the Q-signature")
    for tag in pay[:tags]
        tag in ("diagonal_form", "pfister_form", "hyperbolic_form") &&
            ccheck!(c, tag, std_ok, "standard invariant block on the constructor output")
    end
    if haskey(pay, :rp_omit)
        p = pay[:rp_omit]
        ccheck!(c, "relevant_primes", _istrue(lines_of(c, "rp_clean@$p")),
            "omitted prime $p carries no Hasse content")
    end
end

function _cmp_ldim(c::CmpCtx)
    K, pay = c.K, c.rec.payload
    beta = GWClass(oracle_gram(K, pay[:D]))
    OK = maximal_order(K)
    for p in pay[:support]
        ours = [string(local_anisotropic_dimension(beta, P) == 3)
                for (P, _) in prime_decomposition(OK, p)]
        multiset_check!(c, "local_anisotropic_dimension", "ldim3@$p", ours,
            "odd-rank local dimension test at $p")
    end
end

function _cmp_pair(c::CmpCtx)
    K, pay = c.K, c.rec.payload
    keys_present = !isempty(lines_of(c, "disc_ratio_sq")) &&
                   !isempty(lines_of(c, "sig_eq")) &&
                   all(!isempty(lines_of(c, "hasse_eq@$p")) for p in pay[:support])
    if !keys_present
        c.missing_keys[] += 1
        return
    end
    verdict = _istrue(lines_of(c, "disc_ratio_sq")) &&
              _istrue(lines_of(c, "sig_eq")) &&
              all(_istrue(lines_of(c, "hasse_eq@$p")) for p in pay[:support])
    ours = is_isomorphic_form(GWClass(oracle_gram(K, pay[:D1])),
                              GWClass(oracle_gram(K, pay[:D2])))
    ccheck!(c, "is_isomorphic_form", ours == verdict,
        "vs the rank/disc/Hasse/signature classification: ours=$ours oracle=$verdict")
end

function _cmp_dense(c::CmpCtx)
    for (key, func) in (("det_ratio_sq", "diagonalize_via_congruence"),
                        ("ents_ratio_sq", "diagonal_entries"),
                        ("dcls_ratio_sq", "diagonal_class"))
        L = lines_of(c, key)
        if isempty(L)
            c.missing_keys[] += 1
        else
            ccheck!(c, func, _istrue(L), "determinant square-class preservation ($key)")
        end
    end
    og = lines_of(c, "og_class_eq")
    (!isempty(og) && og[1] != "skip") &&
        ccheck!(c, "diagonalize_via_congruence", og[1] == "true",
            "full class agreement via OrthogonalizeGram")
end

function _cmp_gwu(c::CmpCtx)
    L = lines_of(c, "scale_sq")
    if isempty(L)
        c.missing_keys[] += 1
        return
    end
    ccheck!(c, "GWuClass", c.rec.payload[:accept] == (L[1] == "true"),
        "constructor verdict vs IsSquare(det * c)")
end

function _cmp_dsum(c::CmpCtx)
    L = lines_of(c, "ds_scalar_ok")
    if isempty(L)
        c.missing_keys[] += 1
        return
    end
    ccheck!(c, "divisorial_sum", _istrue(L), "scalar formula recomputed by the oracle")
end

const _COMPARERS = Dict{Symbol, Function}(
    :sym => _cmp_sym, :form => _cmp_form, :ldim => _cmp_ldim,
    :pair => _cmp_pair, :dense => _cmp_dense, :gwu => _cmp_gwu,
    :dsum => _cmp_dsum)

# The exported functions the Magma tier covers, each requiring >= 100 cases
# in a complete run. Covered by the Hecke tier instead: GWClass validation,
# form_rank, stable_part, sum_decomposition_string.
const EXT_QUOTA = 100
const EXT_QUOTA_FUNCTIONS = [
    "GWuClass", "* (integer multiple)", "^ (tensor power)",
    "anisotropic_dimension", "anisotropic_part", "diagonal_class",
    "diagonal_entries", "diagonal_form", "diagonalize_via_congruence",
    "divisorial_sum", "form_signature", "form_signatures", "gw_base_change",
    "gw_direct_sum", "gw_tensor_product", "gw_witt_index",
    "hasse_witt_invariant", "hyperbolic_form", "is_anisotropic_form",
    "is_isomorphic_form", "is_isotropic_form", "local_anisotropic_dimension",
    "local_hilbert_symbol", "pfister_form", "relevant_primes",
    "signed_discriminant", "sum_decomposition",
]

"""
    compare_results(recs, fields, results_path; enforce_quota = true) -> Bool

Compare every corpus record against the parsed results file; print the
per-function table; return overall success. Records absent from the results
file are counted, not failed, so partial Magma runs remain useful — but the
quota is then only enforced if `enforce_quota`.
"""
function compare_results(recs::Vector{ExtRecord}, fields::Vector{DiffField},
                         results_path::String; enforce_quota::Bool = true)
    parsed = parse_results(results_path)
    kmap = Dict(f.name => f.K for f in fields)
    res = DiffResults()
    missing_recs = 0
    missing_keys = Ref(0)
    for r in recs
        lines = get(parsed, r.id, nothing)
        if lines === nothing
            missing_recs += 1
            continue
        end
        _COMPARERS[r.kind](CmpCtx(res, lines, r, kmap[r.fname], missing_keys))
    end
    println("records: $(length(recs)) total, $missing_recs absent from the results file")
    missing_keys[] > 0 &&
        println("keys expected but absent (guarded or partial output): $(missing_keys[])")
    report(res, enforce_quota ? EXT_QUOTA_FUNCTIONS : String[], EXT_QUOTA)
end
