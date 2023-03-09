module LevelSetSublimation

using DrWatson, Reexport
@reexport using SparseArrays
@reexport using DifferentialEquations
@reexport using Contour  
@reexport using Plots

export contour
if !isdefined(LevelSetSublimation, :contour)
    const contour = Contour.contour # This used to be necessary...
end
const CI = CartesianIndex

# Export statements belong in each source file below.

include(srcdir("structs.jl"))
include(srcdir("levelset_plots.jl"))
include(srcdir("levelset_reinit.jl"))
include(srcdir("levelset_advect.jl"))
include(srcdir("solve_T.jl"))
include(srcdir("heat_motion.jl"))
include(srcdir("sim_setup.jl"))
include(srcdir("sim_from_dict.jl"))
include(srcdir("plot_sim_results.jl"))

end # module LevelSetSublimation