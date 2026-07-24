# Branch-coverage instrumentation for the number-field layer. NOT part of the
# public API and OFF by default: every hook is a no-op unless a caller has
# installed a Dict in `_KR_COVERAGE`. It exists for exactly one consumer —
# differential/magma-oracle/generate_oracle.jl builds its corpus
# coverage-driven, and needs to observe which internal branches an input
# exercised (which reduction algorithm ran, whether the Algorithm 4 modulus
# correction was load-bearing, whether Algorithm 8 retried, ...). Nothing in
# src/ or test/ reads these records; when the flag is off the cost of a hook
# is one Ref load.

"""
    _KR_COVERAGE

Off switch for the number-field coverage hooks. `nothing` (the default)
disables all recording. The differential-corpus generator sets it to a
`Dict{Symbol, Any}` around each probe, reads the records off, and resets it
to `nothing`.
"""
const _KR_COVERAGE = Ref{Union{Nothing, Dict{Symbol, Any}}}(nothing)

_kr_coverage_enabled() = _KR_COVERAGE[] !== nothing

# increment a counter
function _kr_cov_count!(key::Symbol)
    cov = _KR_COVERAGE[]
    cov === nothing && return nothing
    cov[key] = get(cov, key, 0) + 1
    nothing
end

# accumulate a boolean with OR (true once seen true)
function _kr_cov_or!(key::Symbol, val::Bool)
    cov = _KR_COVERAGE[]
    cov === nothing && return nothing
    cov[key] = get(cov, key, false) | val
    nothing
end

# accumulate an integer with max
function _kr_cov_max!(key::Symbol, val::Integer)
    cov = _KR_COVERAGE[]
    cov === nothing && return nothing
    cov[key] = max(get(cov, key, typemin(Int)), Int(val))
    nothing
end

# record the value of the FIRST occurrence (later calls do not overwrite)
function _kr_cov_first!(key::Symbol, val)
    cov = _KR_COVERAGE[]
    cov === nothing && return nothing
    haskey(cov, key) || (cov[key] = val)
    nothing
end
