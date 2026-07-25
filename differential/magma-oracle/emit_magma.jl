# Renders the corpus as chunked Magma scripts. Output schema, one
# datum per line:   D|<record id>|<key>|<value>
# A key may repeat within a record (once per prime ideal above the stated
# rational prime, once per real place); the comparator matches repeats as
# multisets. OrthogonalizeGram sits inside try/catch and degrades to a
# "skip" line; everything else uses core intrinsics only.

"An element of K as a Magma expression in the chunk's generator variable."
function mag_elem(gvar::String, a::AbsSimpleNumFieldElem)
    d = degree(parent(a))
    terms = String[]
    for i in 0:(d - 1)
        c = coeff(a, i)
        iszero(c) && continue
        cs = "($(numerator(c))/$(denominator(c)))"
        push!(terms, i == 0 ? cs : i == 1 ? "$cs*$gvar" : "$cs*$gvar^$i")
    end
    isempty(terms) ? "0" : join(terms, " + ")
end

"A diagonal as a Magma sequence with universe K."
mag_diag(Kv::String, gv::String, D::Vector) =
    "[" * Kv * " | " * join([mag_elem(gv, a) for a in D], ", ") * "]"

# Per-field Magma variable names, by position in diff_fields().
_kv(i) = "K$(i)"; _gv(i) = "w$(i)"; _ov(i) = "O$(i)"; _rv(i) = "RP$(i)"
_plv(i, p) = "PL$(i)_$(p)"

# The helper functions every chunk carries. hasse is Lam's prod_{i<j}; hsr is
# the real-place symbol (-1 iff both arguments negative); sd is the signed
# discriminant; pr the plain product.
const MAGMA_PRELUDE = """
SetColumns(0);
_<x> := PolynomialRing(Rationals());
line := procedure(id, key, val)
    printf "D|%o|%o|%o\\n", id, key, val;
end procedure;
hasse := function(d, P)
    h := 1;
    for i in [1..#d-1] do for j in [i+1..#d] do
        h *:= HilbertSymbol(d[i], d[j], P);
    end for; end for;
    return h;
end function;
hsr := function(a, b, v)
    if Evaluate(a, v) lt 0 and Evaluate(b, v) lt 0 then return -1; end if;
    return 1;
end function;
sig := function(d, v)
    s := 0; for a in d do s +:= Sign(Evaluate(a, v)); end for; return s;
end function;
pr := function(K, d)
    p := K!1; for a in d do p *:= a; end for; return p;
end function;
sd := function(K, d)
    n := #d;
    return K!((-1)^((n*(n-1)) div 2)) * pr(K, d);
end function;
"""

# ---------------------------------------------------------------------------
# per-kind emitters: each returns the Magma statements for one record and
# registers the (field, prime) decompositions it needs in `places_used`.

# Primes are ZZRingElem end to end; names and keys print them in decimal.
function _need!(places_used::Set{Tuple{Int, ZZRingElem}}, i::Int, ps)
    for p in ps
        push!(places_used, (i, ZZ(p)))
    end
end

function _emit_sym(r::ExtRecord, i::Int, pu)::Vector{String}
    gv, Kv = _gv(i), _kv(i)
    a, b = mag_elem(gv, r.payload[:a]), mag_elem(gv, r.payload[:b])
    _need!(pu, i, r.payload[:ps])
    out = ["a := $(Kv)!($a); b := $(Kv)!($b);"]
    for p in r.payload[:ps]
        push!(out, "for t in $(_plv(i, p)) do line(\"$(r.id)\", \"hs@$p\", HilbertSymbol(a, b, t[1])); end for;")
    end
    r.payload[:use_real] &&
        push!(out, "for v in $(_rv(i)) do line(\"$(r.id)\", \"hs@real\", hsr(a, b, v)); end for;")
    out
end

# The all-places-agree single boolean: sig at every real place, or hasse at
# every prime above p, compared inside Magma so pairing is trivially right.
_allreal_eq(i, d1, d2, id, key) = """
ok := true;
for v in $(_rv(i)) do if sig($d1, v) ne sig($d2, v) then ok := false; end if; end for;
line("$id", "$key", ok);"""

_allp_hasse_eq(i, p, d1, d2, id, key) = """
ok := true;
for t in $(_plv(i, p)) do if hasse($d1, t[1]) ne hasse($d2, t[1]) then ok := false; end if; end for;
line("$id", "$key", ok);"""

function _emit_form(r::ExtRecord, i::Int, pu)::Vector{String}
    gv, Kv = _gv(i), _kv(i)
    pay = r.payload
    D = pay[:D]
    support = pay[:support]
    _need!(pu, i, support)
    out = String["d := $(mag_diag(Kv, gv, D));"]
    push!(out, "line(\"$(r.id)\", \"sdisc_sq\", IsSquare(sd($Kv, d)));")
    push!(out, "line(\"$(r.id)\", \"sd_match\", IsSquare(sd($Kv, d) * ($(Kv)!($(mag_elem(gv, pay[:sd]))))));")
    push!(out, "for v in $(_rv(i)) do line(\"$(r.id)\", \"sig\", sig(d, v)); end for;")
    for p in support
        push!(out, "for t in $(_plv(i, p)) do line(\"$(r.id)\", \"hasse@$p\", hasse(d, t[1])); end for;")
    end
    length(D) == 2 &&
        push!(out, "line(\"$(r.id)\", \"nps\", IsSquare(-d[1]*d[2]));")

    if haskey(pay, :ap)
        ap, padded = pay[:ap], pay[:padded]
        push!(out, "ap := $(mag_diag(Kv, gv, ap));")
        push!(out, "dp := $(mag_diag(Kv, gv, padded));")
        push!(out, "line(\"$(r.id)\", \"ap_disc_eq\", IsSquare(sd($Kv, d) * sd($Kv, dp)));")
        push!(out, _allreal_eq(i, "d", "dp", r.id, "ap_sig_eq"))
        for p in support
            push!(out, _allp_hasse_eq(i, p, "d", "dp", r.id, "ap_hasse_eq@$p"))
        end
        length(ap) == 2 &&
            push!(out, "line(\"$(r.id)\", \"ap_nps\", IsSquare(-ap[1]*ap[2]));")
        push!(out, "for v in $(_rv(i)) do line(\"$(r.id)\", \"ap_sig\", sig(ap, v)); end for;")
        if haskey(pay, :ap_comp)
            push!(out, "apx := ap cat [$(Kv)!($(mag_elem(gv, pay[:ap_comp])))];")
            for p in support
                push!(out, "for t in $(_plv(i, p)) do line(\"$(r.id)\", \"ap_dim3@$p\", hasse(apx, t[1]) ne HilbertSymbol($(Kv)!(-1), $(Kv)!(-1), t[1])^$(pay[:ap_e])); end for;")
            end
        end
    end

    ex = pay[:extra]
    if haskey(ex, :sum)
        D1, D2 = ex[:sum]
        push!(out, "d1 := $(mag_diag(Kv, gv, D1)); d2 := $(mag_diag(Kv, gv, D2));")
        push!(out, "line(\"$(r.id)\", \"sum_disc_ok\", IsSquare(pr($Kv, d) * pr($Kv, d1) * pr($Kv, d2)));")
        for p in support
            push!(out, """
ok := true;
for t in $(_plv(i, p)) do
    if hasse(d, t[1]) ne hasse(d1, t[1]) * hasse(d2, t[1]) * HilbertSymbol(pr($Kv, d1), pr($Kv, d2), t[1]) then ok := false; end if;
end for;
line("$(r.id)", "sum_hasse_ok@$p", ok);""")
        end
    end
    if haskey(ex, :tensor)
        D1, D2 = ex[:tensor]
        n1, n2 = length(D1), length(D2)
        push!(out, "d1 := $(mag_diag(Kv, gv, D1)); d2 := $(mag_diag(Kv, gv, D2));")
        push!(out, "line(\"$(r.id)\", \"tensor_disc_ok\", IsSquare(pr($Kv, d) * pr($Kv, d1)^$(n2) * pr($Kv, d2)^$(n1)));")
        push!(out, """
ok := true;
for v in $(_rv(i)) do if sig(d, v) ne sig(d1, v) * sig(d2, v) then ok := false; end if; end for;
line("$(r.id)", "tensor_sig_ok", ok);""")
    end
    if haskey(ex, :mult)
        D1, n = ex[:mult]
        push!(out, "d1 := $(mag_diag(Kv, gv, D1));")
        push!(out, "line(\"$(r.id)\", \"mult_disc_ok\", IsSquare(pr($Kv, d) * pr($Kv, d1)^$(n)));")
        push!(out, """
ok := true;
for v in $(_rv(i)) do if sig(d, v) ne $(n) * sig(d1, v) then ok := false; end if; end for;
line("$(r.id)", "mult_sig_ok", ok);""")
    end
    if haskey(ex, :pow)
        D1, n = ex[:pow]
        push!(out, "d1 := $(mag_diag(Kv, gv, D1));")
        push!(out, """
ok := true;
for v in $(_rv(i)) do if sig(d, v) ne sig(d1, v)^$(n) then ok := false; end if; end for;
line("$(r.id)", "pow_sig_ok", ok);""")
    end
    if haskey(ex, :bc_sigq)
        push!(out, """
ok := true;
for v in $(_rv(i)) do if sig(d, v) ne $(ex[:bc_sigq]) then ok := false; end if; end for;
line("$(r.id)", "bc_sig_ok", ok);""")
    end
    if haskey(pay, :rp_omit)
        p = pay[:rp_omit]
        _need!(pu, i, (p,))
        push!(out, """
ok := true;
for t in $(_plv(i, p)) do if hasse(d, t[1]) ne 1 then ok := false; end if; end for;
line("$(r.id)", "rp_clean@$p", ok);""")
    end
    out
end

function _emit_ldim(r::ExtRecord, i::Int, pu)::Vector{String}
    gv, Kv = _gv(i), _kv(i)
    _need!(pu, i, r.payload[:support])
    out = String[
        "dx := $(mag_diag(Kv, gv, vcat(r.payload[:D], [r.payload[:comp]])));",
    ]
    for p in r.payload[:support]
        push!(out, "for t in $(_plv(i, p)) do line(\"$(r.id)\", \"ldim3@$p\", hasse(dx, t[1]) ne HilbertSymbol($(Kv)!(-1), $(Kv)!(-1), t[1])^$(r.payload[:e])); end for;")
    end
    out
end

function _emit_pair(r::ExtRecord, i::Int, pu)::Vector{String}
    gv, Kv = _gv(i), _kv(i)
    D1, D2 = r.payload[:D1], r.payload[:D2]
    _need!(pu, i, r.payload[:support])
    out = String[
        "d1 := $(mag_diag(Kv, gv, D1)); d2 := $(mag_diag(Kv, gv, D2));",
        "line(\"$(r.id)\", \"disc_ratio_sq\", IsSquare(pr($Kv, d1) * pr($Kv, d2)));",
        _allreal_eq(i, "d1", "d2", r.id, "sig_eq"),
    ]
    for p in r.payload[:support]
        push!(out, _allp_hasse_eq(i, p, "d1", "d2", r.id, "hasse_eq@$p"))
    end
    out
end

function _emit_dense(r::ExtRecord, i::Int, pu)::Vector{String}
    gv, Kv = _gv(i), _kv(i)
    A = r.payload[:A]
    n = nrows(A)
    entries = join([mag_elem(gv, A[a, b]) for a in 1:n for b in 1:n], ", ")
    _need!(pu, i, r.payload[:support])
    out = String[
        "A := Matrix($Kv, $n, $n, [$Kv | $entries]);",
        "d := $(mag_diag(Kv, gv, r.payload[:ourD]));",
        "e := $(mag_diag(Kv, gv, r.payload[:ents]));",
        "c := $(mag_diag(Kv, gv, r.payload[:dcls]));",
        "line(\"$(r.id)\", \"det_ratio_sq\", IsSquare(Determinant(A) * pr($Kv, d)));",
        "line(\"$(r.id)\", \"ents_ratio_sq\", IsSquare(Determinant(A) * pr($Kv, e)));",
        "line(\"$(r.id)\", \"dcls_ratio_sq\", IsSquare(Determinant(A) * pr($Kv, c)));",
    ]
    hasse_checks = join(["""
    for t in $(_plv(i, p)) do
        if hasse(gd, t[1]) ne hasse(d, t[1]) then ok := false; end if;
    end for;""" for p in r.payload[:support]], "\n")
    push!(out, """
try
    G := OrthogonalizeGram(A);
    gd := [$Kv | G[k][k] : k in [1..$n]];
    ok := true;
    for v in $(_rv(i)) do if sig(gd, v) ne sig(d, v) then ok := false; end if; end for;
$hasse_checks
    line("$(r.id)", "og_class_eq", ok);
catch err
    line("$(r.id)", "og_class_eq", "skip");
end try;""")
    out
end

_emit_gwu(r::ExtRecord, i::Int, pu)::Vector{String} = String[
    "d := $(mag_diag(_kv(i), _gv(i), r.payload[:D]));",
    "line(\"$(r.id)\", \"scale_sq\", IsSquare(pr($(_kv(i)), d) * ($(_kv(i))!($(mag_elem(_gv(i), r.payload[:c]))))));",
]

function _emit_dsum(r::ExtRecord, i::Int, pu)::Vector{String}
    gv, Kv = _gv(i), _kv(i)
    Ds, roots = r.payload[:Ds], r.payload[:roots]
    k = length(Ds)
    out = String[]
    for (j, D) in enumerate(Ds)
        push!(out, "q$j := $(mag_diag(Kv, gv, D));")
    end
    push!(out, "rr := $(mag_diag(Kv, gv, roots));")
    push!(out, "ex := " * join(["pr($Kv, q$j)" for j in 1:k], " * ") * ";")
    for a in 1:k, b in (a + 1):k
        push!(out, "ex *:= (rr[$a] - rr[$b])^$(2 * length(Ds[a]) * length(Ds[b]));")
    end
    push!(out, "line(\"$(r.id)\", \"ds_scalar_ok\", ($(Kv)!($(mag_elem(gv, r.payload[:scalar])))) eq ex);")
    out
end

const _EMITTERS = Dict{Symbol, Function}(
    :sym => _emit_sym, :form => _emit_form, :ldim => _emit_ldim,
    :pair => _emit_pair, :dense => _emit_dense, :gwu => _emit_gwu,
    :dsum => _emit_dsum)

"""
    emit_magma_chunks(recs, fields, dir; nchunks = 8) -> Vector{String}

Write `oracle_ext_<i>.m` files into `dir`, records in corpus order split into
`nchunks` sequential slices. Each chunk is fully self-contained (field
definitions, helper functions, and exactly the prime decompositions its
records need), so chunks can run in any order or in parallel Magma sessions.
Returns the file paths.
"""
function emit_magma_chunks(recs::Vector{ExtRecord}, fields::Vector{DiffField},
                           dir::String; nchunks::Int = 8)
    fidx = Dict(f.name => i for (i, f) in enumerate(fields))
    per = ceil(Int, length(recs) / nchunks)
    paths = String[]
    for c in 1:nchunks
        slice = recs[((c - 1) * per + 1):min(c * per, length(recs))]
        isempty(slice) && continue
        # first pass: render record bodies, collecting field/place needs
        pu = Set{Tuple{Int, ZZRingElem}}()
        used_fields = Set{Int}()
        bodies = String[]
        for r in slice
            i = fidx[r.fname]
            push!(used_fields, i)
            push!(bodies, "// record $(r.id) ($(r.kind), $(r.fname))")
            append!(bodies, _EMITTERS[r.kind](r, i, pu))
        end
        # second pass: the header those needs imply
        head = String["// oracle_ext chunk $c of $nchunks — GENERATED by",
            "// differential/magma-oracle/generate_all.jl; do not edit.",
            "// Prints one datum per line: D|record|key|value.",
            MAGMA_PRELUDE]
        for i in sort!(collect(used_fields))
            f = fields[i]
            if f.name == "Q"
                push!(head, "$(_kv(i)) := RationalsAsNumberField(); $(_gv(i)) := $(_kv(i))!1;")
            else
                pol = defining_polynomial(f.K)
                pstr = join(["($(numerator(coeff(pol, j)))/$(denominator(coeff(pol, j))))*x^$j"
                             for j in 0:degree(pol)], " + ")
                push!(head, "$(_kv(i))<$(_gv(i))> := NumberField($pstr);")
            end
            push!(head, "$(_ov(i)) := MaximalOrder($(_kv(i))); $(_rv(i)) := RealPlaces($(_kv(i)));")
        end
        for (i, p) in sort!(collect(pu))
            push!(head, "$(_plv(i, p)) := Decomposition($(_ov(i)), $p);")
        end
        path = joinpath(dir, "oracle_ext_$c.m")
        open(path, "w") do io
            for s in vcat(head, bodies)
                println(io, s)
            end
            println(io, "printf \"D|META|chunk_done|$c\\n\";")
        end
        push!(paths, path)
    end
    paths
end
