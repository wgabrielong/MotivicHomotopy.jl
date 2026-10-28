# Shared include stack for generate_all.jl and compare.jl. Reuses the
# hecke-differential harness (DetRNG), oracle helpers (oracle_gram), and
# pathological pools so the two differential tiers draw from one corpus
# discipline.

using Oscar
using MotivicHomotopy

const _EXT_DIR = @__DIR__
const _ROOT = joinpath(_EXT_DIR, "..", "..")
include(joinpath(_ROOT, "test", "nf_test_fields.jl"))
include(joinpath(_ROOT, "differential", "hecke-differential", "harness.jl"))
include(joinpath(_ROOT, "differential", "hecke-differential", "oracles.jl"))
include(joinpath(_ROOT, "differential", "hecke-differential", "generators.jl"))
include(joinpath(_ROOT, "differential", "hecke-differential", "suite_base_change.jl"))
include(joinpath(_ROOT, "differential", "hecke-differential", "suite_transfer.jl"))
include(joinpath(_EXT_DIR, "corpus_ext.jl"))
include(joinpath(_EXT_DIR, "emit_magma.jl"))
include(joinpath(_EXT_DIR, "simulate_hecke.jl"))
include(joinpath(_EXT_DIR, "compare_core.jl"))
