# Coverage measurement for the differential corpus: runs every record through
# the instrumented number-field layer (src/nf_coverage.jl, off by default,
# switched on here around each probe) and reports, per internal branch, which
# records exercise it. generate_oracle.jl includes this file to embed the
# measured matrix in oracle.m and differential/CORPUS.md; it can also be run
# standalone:
#     julia --project=. differential/magma-oracle/coverage.jl
#
# Two probe windows per record, so that driver-internal recomputations of the
# anisotropic dimension cannot contaminate the input-level records:
#   window 1: anisotropic_dimension(input) alone — reads the dimension
#             correction flag for the INPUT form;
#   window 2: anisotropic_part(input) — everything else (reduction chain,
#             Algorithm 4/8 internals, lambda_b, retries, padding).

using Oscar
using MotivicHomotopy

# input-coefficient shape flags, computed from the entries as Magma will see
# them (NOT from the package's normalised diagonal)
function _shape_flags(K, entries)
    OK = maximal_order(K)
    all_units = true
    odd_high = false
    even_dyadic = false
    odd_val = false
    dyadic_val = false
    for a in entries
        aK = K(a)
        isone(denominator(aK, OK)) || error("corpus entries must be integral: $a")
        is_unit(OK(aK)) || (all_units = false)
        for (P, e) in factor(ideal(OK, OK(aK)))
            if minimum(P) == 2
                dyadic_val = true
                e >= 2 && iseven(e) && (even_dyadic = true)
            else
                odd_val = true
                e >= 2 && (odd_high = true)
            end
        end
    end
    (all_units = all_units, odd_high = odd_high, even_dyadic = even_dyadic,
     mixed = odd_val && dyadic_val)
end

# one record -> Dict of measured branch data (empty for the :Q calibration
# record, which never enters the number-field layer)
function measure_record(K, entries)
    K === :Q && return Dict{Symbol, Any}(:calibration => true)
    beta = diagonal_form(K, Tuple(entries))

    MotivicHomotopy._KR_COVERAGE[] = Dict{Symbol, Any}()
    dim_a = anisotropic_dimension(beta)
    w1 = MotivicHomotopy._KR_COVERAGE[]

    MotivicHomotopy._KR_COVERAGE[] = Dict{Symbol, Any}()
    anisotropic_part(beta)
    w2 = MotivicHomotopy._KR_COVERAGE[]
    MotivicHomotopy._KR_COVERAGE[] = nothing

    out = Dict{Symbol, Any}(w2)
    out[:input_dim_a] = dim_a
    out[:dim_correction_decisive] = get(w1, :dim_correction_decisive, false)
    sh = _shape_flags(K, entries)
    out[:shape_all_units] = sh.all_units
    out[:shape_odd_high] = sh.odd_high
    out[:shape_even_dyadic] = sh.even_dyadic
    out[:shape_mixed] = sh.mixed
    out
end

# The matrix rows: (name, predicate on a record's measurement Dict). One row
# per branch the corpus is required to cover; the required count is >= 3
# everywhere except the explicitly-argued exceptions in CORPUS.md.
const MATRIX_ROWS = [
    ("dim_a = 0",                        m -> get(m, :input_dim_a, -1) == 0),
    ("dim_a = 1",                        m -> get(m, :input_dim_a, -1) == 1),
    ("dim_a = 2",                        m -> get(m, :input_dim_a, -1) == 2),
    ("dim_a = 3",                        m -> get(m, :input_dim_a, -1) == 3),
    ("dim_a = 4",                        m -> get(m, :input_dim_a, -1) == 4),
    ("dim_a = 5",                        m -> get(m, :input_dim_a, -1) == 5),
    ("dim_a = 6",                        m -> get(m, :input_dim_a, -1) == 6),
    ("Algorithm 5, real branch",         m -> get(m, :alg5_real_runs, 0) > 0),
    ("Algorithm 5, non-real branch",     m -> get(m, :alg5_nonreal_runs, 0) > 0),
    ("Algorithm 5 iterated (>= 2 runs)", m -> get(m, :alg5_runs, 0) >= 2),
    ("Algorithm 6 (dim 3, non-real)",    m -> get(m, :alg6_runs, 0) > 0),
    ("Algorithm 7 (dim 3, real)",        m -> get(m, :alg7_runs, 0) > 0),
    ("Algorithm 4 called",               m -> get(m, :alg4_runs, 0) > 0),
    ("Algorithm 4 correction load-bearing at a dyadic prime",
                                         m -> get(m, :alg4_alpha1_nonunit_dyadic, false)),
    ("Algorithm 4 correction load-bearing somewhere in S",
                                         m -> get(m, :alg4_alpha1_nonunit_S, false)),
    ("Algorithm 8 reached",              m -> get(m, :alg8_runs, 0) > 0),
    ("Algorithm 8: lambda_b in the E_S basis",
                                         m -> get(m, :alg1_lambda_count, 0) > 0),
    ("Algorithm 8: lambda_b AND a nonempty A block",
                                         m -> get(m, :alg1_lambda_count, 0) > 0 &&
                                              get(m, :alg8_A_rows, 0) > 0),
    ("Algorithm 8: A block empty (real half inert)",
                                         m -> get(m, :alg8_runs, 0) > 0 &&
                                              get(m, :alg8_A_rows, -1) == 0),
    ("Algorithm 8: A block nonempty",    m -> get(m, :alg8_A_rows, 0) > 0),
    ("Algorithm 8: retry loop enlarged S",
                                         m -> get(m, :alg8_retries, 0) > 0),
    ("Algorithm 8 entered with w = 0 mod 4 (no padding)",
                                         m -> get(m, :alg8_w_mod4, -1) == 0),
    ("Algorithm 8 entered with w = 1 mod 4",
                                         m -> get(m, :alg8_w_mod4, -1) == 1),
    ("Algorithm 8 entered with w = 2 mod 4",
                                         m -> get(m, :alg8_w_mod4, -1) == 2),
    ("Algorithm 8 entered with w = 3 mod 4",
                                         m -> get(m, :alg8_w_mod4, -1) == 3),
    ("Algorithm 8 padding changed a dyadic Hasse invariant",
                                         m -> get(m, :alg8_padding_changed_dyadic_hasse, false)),
    ("dimension correction decisive (blind-spot witness)",
                                         m -> get(m, :dim_correction_decisive, false)),
    ("shape: all coefficients units",    m -> get(m, :shape_all_units, false)),
    ("shape: valuation >= 2 at an odd prime",
                                         m -> get(m, :shape_odd_high, false)),
    ("shape: even 2-adic valuation >= 2",
                                         m -> get(m, :shape_even_dyadic, false)),
    ("shape: mixed odd/dyadic valuations",
                                         m -> get(m, :shape_mixed, false)),
]

# measure the whole corpus; returns (per-record Dicts, matrix as
# Vector{(rowname, count, record indices)})
function measure_corpus(corpus; verbose::Bool = true)
    measurements = Vector{Dict{Symbol, Any}}()
    for (i, (label, K, entries, _)) in enumerate(corpus)
        t = @elapsed m = measure_record(K, entries)
        push!(measurements, m)
        verbose && println("  [$i/$(length(corpus))] $(round(t; digits = 2))s  $label")
    end
    matrix = [(name, count(p, measurements),
               [i for (i, m) in enumerate(measurements) if p(m)])
              for (name, p) in MATRIX_ROWS]
    (measurements, matrix)
end

function print_matrix(matrix)
    w = maximum(length(name) for (name, _, _) in matrix)
    for (name, c, idxs) in matrix
        marker = c == 0 ? "  <-- UNCOVERED" : (c == 1 ? "  <-- single" : "")
        println(rpad(name, w), "  ", lpad(c, 3), "  records ", idxs, marker)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    include(joinpath(@__DIR__, "corpus.jl"))
    _, matrix = measure_corpus(corpus)
    println()
    print_matrix(matrix)
end
