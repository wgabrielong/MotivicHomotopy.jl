# Runner for the Hecke-differential suite over number fields.
#
#     julia --project=. differential/hecke-differential/run.jl [suite ...]
#
# With no arguments every suite runs. Naming
# suites (e.g. `run.jl symbols decomposition`) runs a subset — a development
# convenience only; the quota table is meaningful only for full runs.
#
# The run aborts before any suite if a Hecke convention gate fails, prints a
# per-exported-function case table at the end, and exits nonzero on any
# divergence or an unmet quota. Everything is deterministic: same corpus,
# same order, every run on every machine.

using Oscar
using MotivicHomotopy

const HERE = @__DIR__
include(joinpath(HERE, "..", "..", "test", "nf_test_fields.jl"))
include(joinpath(HERE, "harness.jl"))
include(joinpath(HERE, "oracles.jl"))
include(joinpath(HERE, "generators.jl"))
include(joinpath(HERE, "suites_common.jl"))
include(joinpath(HERE, "gates.jl"))
include(joinpath(HERE, "suite_symbols.jl"))
include(joinpath(HERE, "suite_invariants.jl"))
include(joinpath(HERE, "suite_local.jl"))
include(joinpath(HERE, "suite_global.jl"))
include(joinpath(HERE, "suite_decomposition.jl"))
include(joinpath(HERE, "suite_isomorphism.jl"))
include(joinpath(HERE, "suite_ring_ops.jl"))
include(joinpath(HERE, "suite_base_change.jl"))
include(joinpath(HERE, "suite_diagonalization.jl"))
include(joinpath(HERE, "suite_gwu.jl"))
include(joinpath(HERE, "suite_transfer.jl"))

# The exported functions whose number-field support this branch added or
# changed; each must accumulate at least QUOTA cases in a full run. The two
# Base-operator entries are the semiring methods `n * beta` and `beta ^ n`,
# which have no exported name of their own.
const QUOTA = 100
const QUOTA_FUNCTIONS = [
    "GWClass", "GWuClass",
    "anisotropic_dimension", "anisotropic_part",
    "diagonal_class", "diagonal_entries", "diagonal_form",
    "diagonalize_via_congruence", "divisorial_sum",
    "form_rank", "form_signature", "form_signatures",
    "gw_base_change", "gw_direct_sum", "gw_tensor_product", "gw_witt_index",
    "hasse_witt_invariant", "hyperbolic_form",
    "is_anisotropic_form", "is_isomorphic_form", "is_isotropic_form",
    "local_anisotropic_dimension", "local_hilbert_symbol",
    "pfister_form", "relevant_primes",
    "signed_discriminant", "stable_part",
    "sum_decomposition", "sum_decomposition_string",
    "transfer_gw",
    "* (integer multiple)", "^ (tensor power)",
]

const SUITES = [
    ("symbols", suite_symbols!),
    ("invariants", suite_invariants!),
    ("local", suite_local!),
    ("global", suite_global!),
    ("isomorphism", suite_isomorphism!),
    ("ring_ops", suite_ring_ops!),
    ("base_change", suite_base_change!),
    ("diagonalization", suite_diagonalization!),
    ("gwu", suite_gwu!),
    ("decomposition", suite_decomposition!),   # heaviest
    # each suite's seed comes from its position, so a new suite goes last and
    # the suites above keep drawing the same cases
    ("transfer", suite_transfer!),
]

function main(args::Vector{String})
    selected = isempty(args) ? first.(SUITES) : args
    unknown = setdiff(selected, first.(SUITES))
    isempty(unknown) || error("unknown suite(s): $(join(unknown, ", "))")

    println("running Hecke convention gates ...")
    t_gate = @elapsed run_gates()
    println("gates: $(round(t_gate; digits = 1))s")

    # DIFF_MAX_TIER=1 or 2 restricts to the cheaper fields — a development
    # convenience like suite selection; quotas apply to full runs only.
    max_tier = parse(Int, get(ENV, "DIFF_MAX_TIER", "3"))
    fields = [f for f in diff_fields() if f.tier <= max_tier]
    res = DiffResults()
    total = @elapsed for (idx, (name, suite!)) in enumerate(SUITES)
        name in selected || continue
        # per-suite deterministic seed from the suite's fixed position;
        # a subset run draws the same cases as a full run
        rng = DetRNG(0x2109041720260725 + 0x100003 * UInt64(idx))
        t = @elapsed suite!(res, fields, rng)
        println("suite $(rpad(name, 16)) $(round(t; digits = 1))s   " *
                "(cases so far: $(sum(values(res.counts))))")
    end
    println("total suite time: $(round(total; digits = 1))s")

    full_run = sort(selected) == sort(first.(SUITES)) && max_tier == 3
    ok = report(res, full_run ? QUOTA_FUNCTIONS : String[], QUOTA)
    ok || exit(1)
end

main(ARGS)
