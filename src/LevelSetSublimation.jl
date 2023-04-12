module LevelSetSublimation

using DrWatson, Reexport
@reexport using SparseArrays
@reexport using DifferentialEquations
@reexport using Contour  
@reexport using Plots
@reexport using Unitful
using CSV
using ProgressMeter

export contour
if !isdefined(LevelSetSublimation, :contour)
    const contour = Contour.contour # I want access to all of Plots, because I'm lazy, so I have to specify this one.
end
const CI = CartesianIndex

# Export statements belong in each source file below.

include(srcdir("structs.jl"))
include(srcdir("get_vial_dims.jl"))
include(srcdir("levelset_plots.jl"))
include(srcdir("levelset_reinit.jl"))
include(srcdir("levelset_advect.jl"))
include(srcdir("solve_T.jl"))
include(srcdir("physical_data.jl"))
include(srcdir("solve_p.jl"))
include(srcdir("heat_motion.jl"))
include(srcdir("sim_setup.jl"))
include(srcdir("sim_from_dict.jl"))
include(srcdir("plot_sim_results.jl"))

end # module LevelSetSublimation