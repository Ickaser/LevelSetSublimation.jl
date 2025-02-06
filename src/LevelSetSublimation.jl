module LevelSetSublimation

using DrWatson, Reexport
@reexport using SparseArrays
# @reexport using DiffEqCallbacks # Already have this coming from LyoPronto
# @reexport using OrdinaryDiffEqRosenbrock # Already have this coming from LyoPronto
@reexport using OrdinaryDiffEqSSPRK
@reexport using OrdinaryDiffEqBDF
using SteadyStateDiffEq
@reexport using Plots
@reexport using Unitful
@reexport using LaTeXStrings
@reexport using CSV
@reexport using TypedTables
@reexport using PrettyTables
@reexport using StatsPlots: @df
using DataInterpolations: CubicSpline, LinearInterpolation, ExtrapolationType
using Contour: contour, lines, coordinates
using LinearAlgebra: norm, Diagonal
using UnicodePlots: spy
using NonlinearSolve
using LinearSolve
using NaNMath
using PrecompileTools
using DocStringExtensions
using Parameters

@reexport using LyoPronto
import LyoPronto: calc_psub
export calc_psub

const CI = CartesianIndex

# Export statements belong in each source file below.

# A constant which shows up in multiple places, defined globally here:
const θ_THRESH = 0.01

include(srcdir("structs.jl"))
include(srcdir("plotting_tools.jl"))
include(srcdir("levelset_reinit.jl"))
include(srcdir("vel_extrap.jl"))
include(srcdir("heat_transfer.jl"))
include(srcdir("physical_data.jl"))
include(srcdir("mass_transfer.jl"))
include(srcdir("front_motion.jl"))
include(srcdir("sim_setup.jl"))
include(srcdir("time_derivatives.jl"))
include(srcdir("sim_from_dict.jl"))
include(srcdir("access_sim_results.jl"))
include(srcdir("levelset_geometry.jl"))

include(srcdir("without_mass_transfer.jl"))

# include(srcdir("precompile.jl"))

end # module LevelSetSublimation