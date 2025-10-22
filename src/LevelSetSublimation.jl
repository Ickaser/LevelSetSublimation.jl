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
using TypedTables
@reexport using StatsPlots: @df
using DataInterpolations: CubicSpline, LinearInterpolation, ExtrapolationType, ConstantInterpolation
using ComponentArrays
using Contour: contour, lines, coordinates
using LinearAlgebra: norm, Diagonal
using UnicodePlots: spy
using NonlinearSolve
using LinearSolve
using Sparspak
using NaNMath
using PrecompileTools
using DocStringExtensions
using Parameters
using ADTypes

@reexport using LyoPronto
import LyoPronto: calc_psub
export calc_psub

const CI = CartesianIndex

# Export statements belong in each source file below.

# A constant which shows up in multiple places, defined globally here:
const θ_THRESH = 0.01

include("structs.jl")
include("plotting_tools.jl")
include("levelset_reinit.jl")
include("vel_extrap.jl")
include("heat_transfer.jl")
include("physical_data.jl")
include("mass_transfer.jl")
include("front_motion.jl")
include("sim_setup.jl")
include("time_derivatives.jl")
include("sim_from_dict.jl")
include("access_sim_results.jl")
include("levelset_geometry.jl")

include("without_mass_transfer.jl")

# include("precompile.jl")

end # module LevelSetSublimation