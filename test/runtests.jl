using OysterFarming
using Test

@testset "OysterFarming.jl" begin
    include("test_contract.jl")
    include("test_aqua.jl")
end
