# The Magma differential corpus: the deterministic record set from
# which both the Magma script chunks and the comparison are derived. Loaded
# by generate_all.jl and compare.jl after Oscar, MotivicHomotopy,
# test/nf_test_fields.jl, and the hecke-differential harness, oracles, and
# generators. Each record carries the package outputs computed at build
# time; the Magma scripts receive them as literals and print Magma-computed
# verdicts about them (see README.md for the schema).

struct ExtRecord
    id::String
    kind::Symbol      # :sym, :form, :ldim, :pair, :dense, :gwu, :dsum, :transfer
    fname::String     # DiffField name
    payload::Dict{Symbol, Any}
end

# Magma runtime budget: no record may involve a prime above 10^9 or an
# element literal beyond 300 characters; drops are counted and printed.
# Symbol records are exempt (their places are always over 2, 3, 5).
const MAGMA_MAX_SUPPORT_PRIME = ZZ(10)^9
const MAGMA_MAX_ELEM_CHARS = 300
const _MAGMA_BUDGET_DROPPED = Ref(0)

_magma_sized(a) = length(string(a)) <= MAGMA_MAX_ELEM_CHARS

function _within_magma_budget(support::Vector{ZZRingElem}, elem_groups::Vector...)
    ok = all(p -> p <= MAGMA_MAX_SUPPORT_PRIME, support) &&
         all(_magma_sized(a) for D in elem_groups for a in D)
    ok || (_MAGMA_BUDGET_DROPPED[] += 1)
    ok
end

"Does every entry of the diagonal keep its odd-valuation support below the budget?"
_diag_support_ok(K::AbsSimpleNumField, D::Vector) =
    all(p -> p <= MAGMA_MAX_SUPPORT_PRIME, support_rational_primes(K, [K(a) for a in D]))

# The rational primes where any entry of any given diagonal has odd
# valuation, plus 2. ZZRingElem, NOT Int: anisotropic-part representatives
# out of [KR21] Algorithm 8 can carry prime support far beyond 2^63 (that is
# the large-square-content phenomenon the size-reduction boundary exists
# for), and both Hecke and Magma decompose such primes without complaint.
function support_rational_primes(K::AbsSimpleNumField, diags::Vector...)
    OK = maximal_order(K)
    ps = ZZRingElem[ZZ(2)]
    for D in diags, a in D
        iszero(K(a)) && continue
        for (P, e) in factor(K(a) * OK)
            p = minimum(P)
            isodd(e) && !(p in ps) && push!(ps, p)
        end
    end
    sort!(ps)
end

# Anisotropy of a diagonal is certifiable from unguarded Magma data when
# rank <= 1, rank == 2 (IsSquare(-ab)), rank odd <= 3 (the parity-forced
# symbol test), or a real place witnesses definiteness.
function _aniso_certifiable(K::AbsSimpleNumField, D::Vector)::Bool
    r = length(D)
    r <= 3 && return true                    # rank 2 via nps, odd via dim3
    beta = GWClass(oracle_gram(K, D))
    any(s -> abs(s) == r, form_signatures(beta))
end

# The standard :form payload. `tags` names the extra exported functions this
# record credits; `extra` carries their relation data; `with_ap = false`
# skips the anisotropic-part block.
function _form_record!(recs::Vector{ExtRecord}, id::String, f, D::Vector;
                       tags::Vector{String} = String[],
                       extra::Dict{Symbol, Any} = Dict{Symbol, Any}(),
                       with_ap::Bool = true)
    K = f.K
    beta = GWClass(oracle_gram(K, D))
    payload = Dict{Symbol, Any}(
        :D => [K(a) for a in D],
        :sd => K(signed_discriminant(beta)),
        :tags => tags,
        :extra => extra,
    )
    if with_ap
        # A failure here is a package finding: the record is dropped so the
        # kit still ships, and the input is reported for investigation.
        ap_class = try
            anisotropic_part(beta)
        catch err
            @error "anisotropic_part FAILED on a corpus form — record dropped, investigate" field = f.name diagonal = join(string.(D), " ; ") exception = err
            return nothing
        end
        ap = [K(a) for a in diagonal_entries(ap_class)]
        w = gw_witt_index(beta)
        w == div(length(D) - length(ap), 2) ||
            @error "witt index inconsistent with the anisotropic part" field = f.name diagonal = D
        payload[:ap] = ap
        payload[:w] = w
        payload[:dim] = anisotropic_dimension(beta)
        payload[:padded] = vcat(ap, repeat([K(1), K(-1)], w))
        payload[:aniso_cert] = _aniso_certifiable(K, ap)
        payload[:support] = support_rational_primes(K, payload[:D], ap)
        if isodd(length(ap))
            # the parity-forced local test for the part: extend by the
            # complement that makes the signed discriminant a square, then
            # compare Hasse-Witt against the totally hyperbolic target
            # ([KC18] Algorithm 6 shape; only Magma's symbols are consumed)
            r = length(ap)
            sgn = iseven(div(r * (r + 1), 2)) ? 1 : -1
            payload[:ap_comp] = K(sgn) * prod(ap; init = one(K))
            m = div(r + 1, 2)
            payload[:ap_e] = div(m * (m - 1), 2)
        end
    else
        payload[:support] = support_rational_primes(K, payload[:D],
            get(extra, :other_diags, Vector{AbsSimpleNumFieldElem}[])...)
    end
    # the Magma budget, over everything this record makes Magma touch
    groups = Vector[payload[:D], [payload[:sd]]]
    haskey(payload, :ap) && push!(groups, payload[:ap], payload[:padded])
    haskey(payload, :ap_comp) && push!(groups, [payload[:ap_comp]])
    for dd in get(extra, :other_diags, Vector{AbsSimpleNumFieldElem}[])
        push!(groups, dd)
    end
    _within_magma_budget(payload[:support], groups...) || return nothing

    # relevant_primes probe: the smallest small odd prime the package
    # declares irrelevant and at which every entry has even valuation.
    if !isempty(D)
        R = relevant_primes(beta)
        OK = maximal_order(K)
        for p in (3, 5, 7, 11, 13)
            decomp = prime_decomposition(OK, p)
            declared = any(any(Q -> Q == P, R) for (P, _) in decomp)
            evenvals = all(iseven(valuation(K(a), P)) for (P, _) in decomp, a in D)
            if !declared && evenvals
                payload[:rp_omit] = p
                break
            end
        end
    end
    push!(recs, ExtRecord(id, :form, f.name, payload))
    nothing
end

"""
    build_ext_corpus() -> Vector{ExtRecord}

The full corpus, deterministic (fixed DetRNG seed, fixed pools).
Record volumes are chosen so every covered exported function reaches at
least 100 comparisons; `compare_results` enforces this on complete runs.
"""
function build_ext_corpus()
    rng = DetRNG(0x4b52323345585431)          # fixed seed
    fields = diff_fields()
    recs = ExtRecord[]
    _MAGMA_BUDGET_DROPPED[] = 0
    n_id = Ref(0)
    nid(prefix) = string(prefix, lpad(string(n_id[] += 1), 4, '0'))

    for f in fields
        K = f.K
        E = element_pool(f)
        # form pools pre-filtered to the Magma budget; symbol records keep
        # the full element pool
        pool = [D for D in form_pool(f, rng) if _diag_support_ok(K, D)]
        small = [D for D in pool if 1 <= length(D) <= 4]

        # --- :sym — Hilbert symbols at the primes over 2, 3, 5 and the real places
        for _ in 1:14
            a, b = det_pick!(rng, E), det_pick!(rng, E)
            push!(recs, ExtRecord(nid("s"), :sym, f.name, Dict{Symbol, Any}(
                :a => a, :b => b, :ps => [2, 3, 5],
                :use_real => !isempty(real_places(K)))))
        end

        # --- plain :form records, the decomposition workhorses
        for D in pool[1:min(20, length(pool))]
            length(D) <= 8 || continue
            _form_record!(recs, nid("f"), f, D; tags = String["diagonal_form"])
        end

        # --- :ldim — odd-rank local-dimension records (the parity-forced
        # symbol test needs no local-square data)
        odd_forms = [D for D in pool if isodd(length(D)) && length(D) <= 7]
        for D in odd_forms[1:min(10, length(odd_forms))]
            r = length(D)
            sgn = iseven(div(r * (r + 1), 2)) ? 1 : -1
            comp = K(sgn) * prod(K.(D); init = one(K))
            m = div(r + 1, 2)
            push!(recs, ExtRecord(nid("l"), :ldim, f.name, Dict{Symbol, Any}(
                :D => [K(a) for a in D], :comp => comp,
                :e => div(m * (m - 1), 2),
                :support => support_rational_primes(K, [K(a) for a in D]))))
        end

        # --- :pair — is_isomorphic_form vs the Hasse-Minkowski classification
        for _ in 1:10
            D1 = det_pick!(rng, small)
            variant = det_int!(rng, 3)
            D2 = variant == 1 ? reverse(D1) :
                 variant == 2 ? [det_pick!(rng, E)^2 * a for a in D1] :
                 begin
                     t = copy(D1)
                     t[det_int!(rng, length(t))] *= det_pick!(rng, E)
                     t
                 end
            length(D1) == length(D2) || continue
            push!(recs, ExtRecord(nid("p"), :pair, f.name, Dict{Symbol, Any}(
                :D1 => [K(a) for a in D1], :D2 => [K(a) for a in D2],
                :support => support_rational_primes(K, [K(a) for a in D1],
                                                    [K(a) for a in D2]))))
        end

        # --- :dense — diagonalization of non-diagonal Gram matrices
        denses = dense_pool(f, rng)
        n_dense = 0
        for A in denses
            n_dense >= 9 && break
            beta = GWClass(A)
            ourD = [K(A2) for A2 in _diag_of(diagonalize_via_congruence(A))]
            ents = [K(a) for a in diagonal_entries(beta)]
            dcls = [K(a) for a in diagonal_entries(diagonal_class(beta))]
            support = support_rational_primes(K, ourD, ents)
            # re-check the budget: diagonalization mints fresh primes
            _within_magma_budget(support, ourD, ents, dcls) || continue
            n_dense += 1
            push!(recs, ExtRecord(nid("d"), :dense, f.name, Dict{Symbol, Any}(
                :A => A, :ourD => ourD, :ents => ents, :dcls => dcls,
                :support => support)))
        end

        # --- :gwu — the scalar/determinant square-class rule
        for i in 1:9
            D = det_pick!(rng, small)
            d = prod(K.(D); init = one(K))
            u = first(a for a in E if !is_square(a))
            c = i % 3 == 0 ? d * u : i % 3 == 1 ? d : d * det_pick!(rng, E)^2
            accepted = try
                GWuClass(oracle_gram(K, D), c)
                true
            catch
                false
            end
            push!(recs, ExtRecord(nid("u"), :gwu, f.name, Dict{Symbol, Any}(
                :D => [K(a) for a in D], :c => c, :accept => accepted)))
        end

        # --- :dsum — the divisorial-sum scalar formula, recomputed by Magma
        roots_pool = AbsSimpleNumFieldElem[K(0), K(1), K(-2), f.gen, f.gen + 3]
        for _ in 1:8
            k = det_int!(rng, 2) + 1
            Ds = [det_pick!(rng, small) for _ in 1:k]
            roots = AbsSimpleNumFieldElem[]
            for r in roots_pool
                any(x -> x == r, roots) || push!(roots, r)
                length(roots) == k && break
            end
            length(roots) == k || continue
            classes = [GWuClass(oracle_gram(K, [K(a) for a in D])) for D in Ds]
            ds = divisorial_sum(classes, roots)
            _magma_sized(K(gw_scalar(ds))) || (_MAGMA_BUDGET_DROPPED[] += 1; continue)
            push!(recs, ExtRecord(nid("v"), :dsum, f.name, Dict{Symbol, Any}(
                :Ds => [[K(a) for a in D] for D in Ds], :roots => roots,
                :scalar => K(gw_scalar(ds)))))
        end

        # --- composite outputs, registered as :form records with relation data
        for _ in 1:8   # direct sum, with the Hasse sum formula relation
            D1, D2 = det_pick!(rng, small), det_pick!(rng, small)
            S = [K(a) for a in _diag_of(gw_matrix(gw_direct_sum(
                    GWClass(oracle_gram(K, D1)), GWClass(oracle_gram(K, D2)))))]
            _form_record!(recs, nid("f"), f, S; tags = String["gw_direct_sum"],
                extra = Dict{Symbol, Any}(:sum => ([K(a) for a in D1], [K(a) for a in D2]),
                    :other_diags => [[K(a) for a in D1], [K(a) for a in D2]]),
                with_ap = false)
        end
        for _ in 1:8   # tensor product: disc and per-place signature relations
            D1, D2 = det_pick!(rng, small), det_pick!(rng, small)
            T = [K(a) for a in _diag_of(gw_matrix(gw_tensor_product(
                    GWClass(oracle_gram(K, D1)), GWClass(oracle_gram(K, D2)))))]
            length(T) <= 12 || continue
            _form_record!(recs, nid("f"), f, T; tags = String["gw_tensor_product"],
                extra = Dict{Symbol, Any}(:tensor => ([K(a) for a in D1], [K(a) for a in D2]),
                    :other_diags => [[K(a) for a in D1], [K(a) for a in D2]]),
                with_ap = false)
        end
        for i in 1:8   # integer multiple and tensor power
            D = det_pick!(rng, small)
            n = 2 + i % 2
            M = [K(a) for a in _diag_of(gw_matrix(n * GWClass(oracle_gram(K, D))))]
            _form_record!(recs, nid("f"), f, M; tags = String["* (integer multiple)"],
                extra = Dict{Symbol, Any}(:mult => ([K(a) for a in D], n),
                    :other_diags => [[K(a) for a in D]]),
                with_ap = false)
            # tensor powers from the rank <= 2 sub-pool (rank^n <= 8)
            tiny = [T for T in small if length(T) <= 2]
            Dp = det_pick!(rng, tiny)
            P = [K(a) for a in _diag_of(gw_matrix(GWClass(oracle_gram(K, Dp))^n))]
            _form_record!(recs, nid("f"), f, P; tags = String["^ (tensor power)"],
                extra = Dict{Symbol, Any}(:pow => ([K(a) for a in Dp], n),
                    :other_diags => [[K(a) for a in Dp]]),
                with_ap = false)
        end
        for _ in 1:10   # Pfister forms, with the decomposition block
            a, b = det_pick!(rng, E), det_pick!(rng, E)
            Pf = [K(x) for x in _diag_of(gw_matrix(pfister_form(K, (a, b))))]
            _form_record!(recs, nid("f"), f, Pf; tags = String["pfister_form"])
        end
        for i in 1:8   # hyperbolic forms, with the decomposition block
            k = i % 3 + 1
            H = [K(x) for x in _diag_of(gw_matrix(hyperbolic_form(K, 2 * k)))]
            _form_record!(recs, nid("f"), f, H; tags = String["hyperbolic_form"])
        end
        for Dq in BASE_CHANGE_QQ_DIAGONALS[1:8]   # base change from Q
            q = diagonal_form(QQ, Tuple(Dq))
            bc = gw_base_change(K, q)
            DK = [K(a) for a in _diag_of(gw_matrix(bc))]
            _form_record!(recs, nid("f"), f, DK; tags = String["gw_base_change"],
                extra = Dict{Symbol, Any}(:bc_sigq => form_signature(q)),
                with_ap = false)
        end
    end
    # --- :transfer — transfer_gw along Tr_{K/QQ}, and along Tr_{L/K} for a
    # relative quadratic extension L = K(sqrt m), recomputed by Magma from the
    # input form alone. These records come after all the others and draw from a
    # stream of their own, so the records above, and the chunks they fill, are
    # the same as without them.
    trng = DetRNG(0x5452414e53464552)         # its own fixed seed
    for f in fields
        K = f.K
        for n in (1, 1, 2, 2, 3, 3)
            B = _transfer_gram(K, [one(K), f.gen], n, trng)
            push!(recs, ExtRecord(nid("t"), :transfer, f.name,
                                  Dict{Symbol, Any}(:B => B, :m => nothing)))
        end
        m = first(m for m in (2, 3, 5, 6, 7, -1, -3) if !is_square(K(m)))
        Kz, z = K["z"]
        L, r = number_field(z^2 - m, "r")
        for n in (1, 2, 2)
            B = _transfer_gram(L, [one(L), r, L(f.gen)], n, trng)
            push!(recs, ExtRecord(nid("t"), :transfer, f.name,
                                  Dict{Symbol, Any}(:B => B, :m => m)))
        end
    end

    _MAGMA_BUDGET_DROPPED[] > 0 &&
        println("Magma budget: $(_MAGMA_BUDGET_DROPPED[]) records dropped " *
                "(support prime > 10^9 or oversized literal); the shapes remain " *
                "covered by the local hecke-differential tier")
    recs
end

"The diagonal of a square matrix as a Vector (helper shared by the emitters)."
_diag_of(M::MatElem) = [M[i, i] for i in 1:nrows(M)]

# A transferred Gram matrix as the rationals the Magma chunk prints for it: its
# entries row by row, each entry over a number field written by its coordinates
# in the power basis.
function transfer_coordinates(G::MatElem)
    entries = [G[i, j] for i in 1:nrows(G) for j in 1:ncols(G)]
    base_ring(G) isa QQField && return entries
    return reduce(vcat, [collect(coordinates(e)) for e in entries])
end

"Rationals as Magma prints them, comma-separated: p/q, or p when q = 1."
magma_rationals(v::Vector) =
    join([isone(denominator(q)) ? string(numerator(q)) :
          "$(numerator(q))/$(denominator(q))" for q in v], ",")

"The inverse of [`magma_rationals`](@ref)."
parse_magma_rationals(s::AbstractString) =
    [occursin('/', t) ? QQ(parse(BigInt, split(t, '/')[1]), parse(BigInt, split(t, '/')[2])) :
                        QQ(parse(BigInt, t)) for t in split(s, ',')]
