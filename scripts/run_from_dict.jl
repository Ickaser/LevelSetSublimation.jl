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

include(srcdir("sim_from_dict.jl"))

# ------------------------------- Parameters to control

# casename = "box_all"

sim_dt = 1.0

# T_params = make_artificial_params()

Q_gl = 1.0
Q_sh = 1.0
Q_ic = 1.0
Q_ck = 0.0
k = 1.0
Tf = 250.0
ΔH = 1.0 # 
ρf = 100.0 
# " Takes T as Kelvin, returns P in Pa"
# function calc_psub(T)
#     ai = [-0.212144006e2,  0.273203819e2,  -0.610598130e1]
#     bi = [0.333333333e-2,  0.120666667e1,  0.170333333e1]
#     θ = T/275.16
#     lnπ = sum(ai .* θ .^bi) / θ
#     exp(lnπ)*611.657
# end
# Rw = 8.3145 / .018 # J/molK * mol/kg
# calc_ρvap(T) = calc_psub(T)/Rw/T
T_params = Dict{Symbol, Any}()
@pack! T_params = Q_gl, Q_sh, Q_ic, Q_ck, k, Tf, ΔH, ρf

println("Parameters loaded")

# ----------------- Domain setup


nr = 50
nz = 45
# rmin = 0.0
# zmin = 0.0
rmax = 1.0
zmax = 1.0

# Default "bandwidth" size: 20% of domain size, extending on both sides of front
# This default is built into Domain constructor
# bwr = ceil(Int, 0.2*(rmax-rmin)/dr)
# bwz = ceil(Int, 0.2*(zmax-zmin)/dz)

dom = Domain(nr, nz, rmax, zmax)
# dom = Domain(nr, nz, rmin, rmax, zmin, zmax)
# dom = Domain(nr, nz, rmin, rmax, zmin, zmax, bwr, bwz)
# dom = Domain(nr, nz, rmax, zmax, bwr, bwz)


# --------------------- Set up initial contour ϕ0

# Flat on top
ϕ0 = [ z - 0.999  for r in dom.rgrid, z in dom.zgrid]
# initshape = "flat"

# Cylinder
# ϕ0 = [r - 0.999 for r in dom.rgrid, z in dom.zgrid]

# Box
# ϕ0 = [max(r,z)-0.999 for r in dom.rgrid, z in dom.zgrid]

# Ellipse
# ϕ0 = [r^2 + 2z^2 - 1.5 for r in dom.rgrid, z in dom.zgrid]
# initshape = "ellipse"

# Separated ellipse
# ϕ0 = [1.5r^2 + 6(z-0.5)^2 - 1.0 for r in dom.rgrid, z in dom.zgrid]
# initshape = "separated"

# Circle
# ϕ0 = [1.1r^2 + 1.1z^2 - 1.0 for r in dom.rgrid, z in dom.zgrid]
# initshape = "circle"

# Line
# ϕ0 = [0.5r + z - 0.9 + 0.001 for r in dom.rgrid, z in dom.zgrid]
# ϕ0 = [r + z - 1.5 + 0.001 for r in dom.rgrid, z in dom.zgrid]
# initshape = "cone"

# ----------------------- COmbine into a dict

# one_config = Dict(
#     "T_params" => T_params,
#     "sim_dt" => sim_dt,
#     "dom" => dom,
#     "ϕ0" => ϕ0,
# )

# To do multiple combinations, include a vector here
all_configs = Dict(
    "T_params" => T_params,
    # "sim_dt" => sim_dt,
    "sim_dt" => [sim_dt*0.5, sim_dt],
    "dom" => dom,
    "ϕ0" => ϕ0,
)


# -------------- Use DrWatson to run all simulations everything

list_configs = dict_list(all_configs)

pol_kwargs = (filename=hash, prefix="simdat", verbose=false, tag=true)

for conf in list_configs
    produce_or_load(sim_from_dict, conf, datadir("sims", "testing"); pol_kwargs...)
end
# @tagsave