using LevelSetSublimation
using UnPack
using Test


const LSS = LevelSetSublimation

# Define my own convenience function, since have lots of arrays to check
approxzero(x) = isapprox(x, 0, atol=100eps(typeof(x)))

# Run test suite
println("Starting tests\n")
ti = time()

include("tests_setup.jl")

include("tests_reinit.jl")

include("tests_levelsetgeometry.jl")
dom1 = Domain(25, 24, 1.3, 2.0)
dom2 = Domain(45, 35, 1.2, 2.1)
@testset "state variable array construction" for dom in [dom1, dom2]
# This might be redundant after using ComponentArrays
    init_prof = :circ
    Tf0 = 233u"K"
    Tvw0 = 245u"K"
    u = LevelSetSublimation.make_u0_ndim(init_prof, Tf0, Tvw0, dom)
    ϕ1 = LevelSetSublimation.make_ϕ0(init_prof, dom)
    @testset "make_u0_ndim, first segment" for i in eachindex(ϕ1)
        @test u.ϕ[i] == ϕ1[i]
    end
    @testset "second segment" for i in eachindex(u.Tf)
        @test u.Tf[i] == ustrip(u"K", Tf0)
    end
    @test u.Tvw == ustrip(u"K", Tvw0)

    # Test with views
    ϕr, Tfr = u.ϕ, u.Tf
    ϕr .= π
    Tfr .= 250
    u.Tvw = 260
    @testset "make_u0_ndim, first segment" for i in eachindex(ϕ1)
        @test ϕr[i] ≈ π
    end
    @testset "second segment" for i in eachindex(u.Tf)
        @test u.Tf[i] == 250
    end
    @test length(u.Tvw) == 1
    @test u.Tvw == 260
end

include("tests_lyopronto.jl")

@info "Getting into the empty test sets now"

@testset "T solution: nothing here" begin
    # Compare to analytical:
    # 1D r direction, no ice 
    # 1D z direction, no ice 
    # 2D no ice 
    # 1D r direction, ghost θ > θ_THRESH 
    # 1D r direction, ghost θ < θ_THRESH 
    # 1D z direction, ghost θ > θ_THRESH 
    # 1D z direction, ghost θ < θ_THRESH 
end

@testset "T derivatives for velocity: no tests yet" begin
    # 1D r direction, ghost θ > θ_THRESH 
    # 1D r direction, ghost θ < θ_THRESH 
    # 1D z direction, ghost θ > θ_THRESH 
    # 1D z direction, ghost θ < θ_THRESH 
    # 
end

# Tests for the old simulation version that did not treat mass transfer rigorously
# These are broken, but also may not be useful anyway,
# so unless something changes I am not maintaining them
# include("tests_1DT_motion.jl")

@testset "p solution: nothing here" begin
end

@testset "p derivatives for velocity: no tests yet" begin
end

@testset "velocity extrapolation: no tests yet" begin
end

@testset "radial T solution: no tests yet" begin
end

ti = time() - ti
println("\nTesting took:")
println(round(ti/60, digits = 3), " minutes")
