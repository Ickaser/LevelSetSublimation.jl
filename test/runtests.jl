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

include(testdir("tests_1DT_motion.jl"))

@testset "p solution: nothing here" begin
end

@testset "p derivatives for velocity: no tests yet" begin
end

@testset "velocity extrapolation" begin
end

include(testdir("tests_lyopronto.jl"))

include(testdir("tests_radialT.jl"))

ti = time() - ti
println("\nTesting took:")
println(round(ti/60, digits = 3), " minutes")
