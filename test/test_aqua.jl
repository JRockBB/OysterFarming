using Aqua

@testset "Aqua quality checks" begin
    Aqua.test_all(OysterFarming)
end
