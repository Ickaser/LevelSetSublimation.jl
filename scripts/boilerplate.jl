using DrWatson
@quickactivate
using SparseArrays

# using StaticArrays
# using DomainSets
using Contour
# using Interpolations
using Plots
using DifferentialEquations


include(srcdir("structs.jl"))
include(srcdir("levelset_plots.jl"))
include(srcdir("levelset_reinit.jl"))
include(srcdir("levelset_advect.jl"))
include(srcdir("solve_T.jl"))
include(srcdir("coupled_motion.jl"))