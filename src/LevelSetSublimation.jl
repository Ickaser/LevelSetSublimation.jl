module LevelSetSublimation

using DrWatson, Reexport
@reexport using SparseArrays
@reexport using DifferentialEquations
@reexport using Contour  
@reexport using Plots
@reexport using Unitful
using CSV
using ProgressMeter
using LinearAlgebra: norm
using NLsolve
using NonlinearSolve
using LinearSolve
using NaNMath

export contour
if !isdefined(LevelSetSublimation, :contour)
    const contour = Contour.contour # I want access to all of Plots, because I'm lazy, so I have to specify this one.
end
const CI = CartesianIndex

# Export statements belong in each source file below.

include(srcdir("structs.jl"))
include(srcdir("get_vial_dims.jl"))
include(srcdir("plotting_tools.jl"))
include(srcdir("levelset_reinit.jl"))
include(srcdir("vel_extrap.jl"))
include(srcdir("solve_T.jl"))
include(srcdir("physical_data.jl"))
include(srcdir("solve_p.jl"))
include(srcdir("front_motion.jl"))
include(srcdir("sim_setup.jl"))
include(srcdir("time_derivatives.jl"))
include(srcdir("sim_from_dict.jl"))
include(srcdir("access_sim_results.jl"))
include(srcdir("levelset_geometry.jl"))

end # module LevelSetSublimation