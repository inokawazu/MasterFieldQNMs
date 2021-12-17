using JET, Test, MasterFieldQNMs

@testset "Jet Testing" begin
    shear = MasterFieldQNMs.BlackBrane5D.Shear()

    T = ComplexF64

    w0 = T(1.4436414 − 1.0692250im)
    q0 = T(sqrt(1.4436414 − 1.0692250im))
    @testset "Performance Testing" begin
        @test_opt target_modules=(@__MODULE__,) horexp(shear, w0, q0)
        @test_opt target_modules=(@__MODULE__,) hubeny_horowitz_criticalpoint(shear, w0, q0)
    end
    @testset "Error Testing" begin
        @test_call target_modules=(@__MODULE__,) horexp(shear, w0, q0)
        @test_call target_modules=(@__MODULE__,) hubeny_horowitz_criticalpoint(shear, w0, q0)
    end

end
