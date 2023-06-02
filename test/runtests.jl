using DrWatson, Test
@quickactivate :LevelSetSublimation


const LSS = LevelSetSublimation

# Define my own convenience function, since have lots of arrays to check
approxzero(x) = isapprox(x, 0, atol=100eps(typeof(x)))
# Directory function
testdir(args...) = projectdir("test", args...)


# Run test suite
println("Starting tests\n")
ti = time()

include(testdir("tests_setup.jl"))

include(testdir("tests_reinit.jl"))

@testset "T solution: nothing here" begin
    # Compare to analytical:
    # 1D r direction, no ice 
    # 1D z direction, no ice 
    # 2D no ice 
    # 1D r direction, ghost θ > θ_thresh 
    # 1D r direction, ghost θ < θ_thresh 
    # 1D z direction, ghost θ > θ_thresh 
    # 1D z direction, ghost θ < θ_thresh 
end

@testset "T derivatives for velocity: no tests yet" begin
    # 1D r direction, ghost θ > θ_thresh 
    # 1D r direction, ghost θ < θ_thresh 
    # 1D z direction, ghost θ > θ_thresh 
    # 1D z direction, ghost θ < θ_thresh 
    # 
end

@testset "p solution: nothing here" begin
end

@testset "p derivatives for velocity: no tests yet" begin
end

@testset "velocity extrapolation" begin
end

@testset "comparison to LyoPronto, default parameters" begin
    # Read in precomputed results from LyoPronto
    # These are generated in scripts/default_lyopronto_config.jl, then saved
    lyopronto_res = load(testdir("lyopronto_default_results.jld2"))
    @unpack t_lp, T_lp, m_lp, dryfrac_lp = lyopronto_res
    # Read in a config object
    include(testdir("lyopronto_default_setup.jl"))
    res = sim_from_dict(config)
    tsol, Tsol, msol, fsol = compare_lyopronto_res(t_lp, res, config)

    @test sum(t_lp .== tsol) == 100
    @test sum(isapprox.(Tsol, T_lp, atol=0.5u"K")) == 100
    @test sum(isapprox.(msol, m_lp, atol=0.05u"kg/hr/m^2")) >= 98
    @test sum(isapprox.(fsol, dryfrac_lp, atol=0.02)) == 100
    @test res["sol"].t[end]*u"s" ≈ t_lp[end] rtol=0.02
end

include(testdir("tests_radialT.jl"))

ti = time() - ti
println("\nTesting took:")
println(round(ti/60, digits = 3), " minutes")
