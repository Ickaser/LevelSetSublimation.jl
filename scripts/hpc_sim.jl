using Pkg
Pkg.activate("..")
Pkg.instantiate()
using DrWatson
@quickactivate :LevelSetSublimation

const LSS = LevelSetSublimation

if abspath(PROGRAM_FILE) == @__FILE__ # Run as individual script
    inputfile = projectdir(ARGS[1])
    @info "Reading input from: $(projectdir(inputfile))"
else
    @info "hpc_sim.jl run not as standalone, using `example_input.jl` as input"
    # inputfile = scriptsdir("nosync_siminputs", "input_highcurvature.jl")
    inputfile = scriptsdir("example_input.jl")
end

include(projectdir(inputfile))
    
pol_kwargs = (filename=hash, prefix="sim_hpc", verbose=false, tag=true)
# dom = Domain(51, 51, 1.0, 1.0)
@time produce_or_load(sim_and_postprocess, config, datadir("sims"); pol_kwargs...)


# resultsanim(res, config, "hpc_test", seconds_length=10)
# resultsanim(res, config, "hpc_test", seconds_length=10, heatvar=:p)
