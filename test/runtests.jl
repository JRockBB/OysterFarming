using OysterFarming
using Test

@testset "OysterFarming.jl" begin
    include("test_contract.jl")
    include("test_params.jl")
    include("test_deb.jl")
    include("test_forcing.jl")
    include("test_calibration.jl")
    include("test_farm.jl")
    include("test_nitrogen.jl")
    include("test_dermo.jl")
    include("test_stochastic.jl")
    include("test_vibrio.jl")
    include("test_economics.jl")
    include("test_aqua.jl")
end
