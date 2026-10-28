# Regenerates the Magma oracle kit:
#
#     julia --project=. differential/magma-oracle/generate_all.jl
#
# 1. builds the deterministic corpus (computing the package outputs that get
#    injected into the Magma script);
# 2. writes the chunked oracle_ext_<i>.m scripts;
# 3. writes results-hecke-simulated.txt in the same line schema;
# 4. compares the corpus against the simulation, which both validates the
#    pipeline and prints the per-function case table with the >= 100 quota
#    enforced.
#
# Exit is nonzero if the simulated comparison diverges or a quota is unmet.

include(joinpath(@__DIR__, "loader.jl"))

function main()
    fields = diff_fields()
    println("building the corpus (computes anisotropic parts; takes minutes) ...")
    t = @elapsed recs = build_ext_corpus()
    println("corpus: $(length(recs)) records in $(round(t; digits = 1))s")

    # the transfer records get chunks of their own after the eight of the
    # others, so those keep their boundaries; two, so that each fits the online
    # calculator's input limit
    others = [r for r in recs if r.kind != :transfer]
    transfers = [r for r in recs if r.kind == :transfer]
    paths = vcat(emit_magma_chunks(others, fields, @__DIR__),
                 emit_magma_chunks(transfers, fields, @__DIR__; nchunks = 2, first = 9))
    total_kb = round(Int, sum(filesize, paths) / 1024)
    println("Magma chunks: $(length(paths)) files, $(total_kb) KiB")

    sim = joinpath(@__DIR__, "results-hecke-simulated.txt")
    t = @elapsed simulate_results(recs, fields, sim)
    println("Hecke simulation written in $(round(t; digits = 1))s: $sim")

    println("\ncomparing package vs the simulated results ...")
    t = @elapsed ok = compare_results(recs, fields, sim)
    println("comparison time: $(round(t; digits = 1))s")
    ok || exit(1)
end

main()
