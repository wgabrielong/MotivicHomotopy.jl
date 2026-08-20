using Test

@testset "MotivicHomotopy.jl" begin
    include("test_core_types.jl")
    include("test_matrix_layer.jl")
    include("test_invariants.jl")
    include("test_decomposition.jl")
    include("test_degrees.jl")
    include("test_v2_features.jl")
    include("test_number_fields.jl")
    include("test_number_fields_aniso.jl")
    include("test_number_fields_oracle.jl")
    include("test_ring_ops.jl")
    include("test_numerical.jl")
    include("test_certified.jl")
end
