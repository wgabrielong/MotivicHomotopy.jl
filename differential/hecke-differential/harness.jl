# Per-exported-function case counting and failure reporting. Every check is
# credited to exactly one exported function. A failure never throws: it is
# recorded and printed with reproduction context, and the runner exits
# nonzero at the end.

mutable struct DiffResults
    counts::Dict{String, Int}
    failures::Vector{String}
end

DiffResults() = DiffResults(Dict{String, Int}(), String[])

"""
    check!(res::DiffResults, exported_function::String, ok::Bool, context::String)

Record one differential case for `exported_function`. `context` must let a
reader reproduce the case: field, diagonal (or Gram matrix), place, expected
vs. got.
"""
function check!(res::DiffResults, exported_function::String, ok::Bool, context::String)
    res.counts[exported_function] = get(res.counts, exported_function, 0) + 1
    if !ok
        msg = "[$exported_function] DIVERGENCE: $context"
        push!(res.failures, msg)
        println(stderr, msg)
    end
    ok
end

"""
    report(res::DiffResults, quota_functions::Vector{String}, quota::Int) -> Bool

Print the per-function case table and every recorded divergence; return
whether the run is clean AND every function in `quota_functions` reached
`quota` cases.
"""
function report(res::DiffResults, quota_functions::Vector{String}, quota::Int)
    println("\n=== differential cases per exported function ===")
    width = maximum(length.(collect(keys(res.counts))); init = 20)
    for name in sort!(collect(keys(res.counts)))
        short = res.counts[name] < quota && name in quota_functions
        println(rpad(name, width + 2), lpad(res.counts[name], 6), short ? "   << BELOW QUOTA" : "")
    end
    missing_fns = [f for f in quota_functions if get(res.counts, f, 0) < quota]
    println("total cases: ", sum(values(res.counts)))
    println("divergences: ", length(res.failures))
    isempty(missing_fns) ||
        println("functions below the $quota-case quota: ", join(missing_fns, ", "))
    isempty(res.failures) && isempty(missing_fns)
end

# Deterministic pseudo-randomness: a fixed linear congruential generator
# (Knuth's MMIX constants), identical on every machine and Julia version.
mutable struct DetRNG
    state::UInt64
end

DetRNG() = DetRNG(0x2109041720260725)   # fixed seed

function det_next!(rng::DetRNG)::UInt64
    rng.state = rng.state * 0x5851f42d4c957f2d + 0x14057b7ef767814f
    rng.state
end

"A deterministic integer in 1..n."
det_int!(rng::DetRNG, n::Int) = Int(det_next!(rng) >> 33) % n + 1

"A deterministic element of a vector."
det_pick!(rng::DetRNG, v::Vector) = v[det_int!(rng, length(v))]
