using OysterFarming
using Test

@testset "OysterFarming.jl" begin
    include("test_contract.jl")
    include("test_params.jl")
    include("test_deb.jl")
    include("test_aqua.jl")
end
