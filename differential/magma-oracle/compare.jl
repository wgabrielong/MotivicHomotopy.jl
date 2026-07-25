# Compares the package against a results file printed by the oracle_ext_*.m
# chunks (concatenate multiple chunk outputs into one file in any order):
#
#     julia --project=. differential/magma-oracle/compare.jl <results.txt>
#
# With no argument it uses results-hecke-simulated.txt — the pipeline
# self-check, NOT a Magma verdict. Partial results files are fine: absent
# records are counted and skipped, and the >= 100-per-function quota is
# enforced only when nothing is missing.

include(joinpath(@__DIR__, "loader.jl"))

function main(args::Vector{String})
    path = isempty(args) ? joinpath(@__DIR__, "results-hecke-simulated.txt") : args[1]
    isfile(path) || error("no results file at $path")
    println("building the corpus (deterministic; recomputes the package outputs) ...")
    recs = build_ext_corpus()
    fields = diff_fields()
    parsed = parse_results(path)
    complete = all(haskey(parsed, r.id) for r in recs)
    ok = compare_results(recs, fields, path; enforce_quota = complete)
    complete || println("NOTE: results file is partial - quota not enforced on this run")
    ok || exit(1)
end

main(ARGS)
