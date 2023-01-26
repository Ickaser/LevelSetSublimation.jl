"Assume that DrWatson is already imported and the LevelSetMethods project activated"
# using DrWatson
# @quickactivate :LevelSetMethods


# ------------------------------- Parameters to control

# casename = "box_all"
ϕ0type = :flat # Options: :flat, :cyl, :box, :cone, etc. See make_ϕ0

sim_dt = 1.0

# T_params = make_artificial_params()

Q_gl = 0.0
Q_sh = 1.0
Q_ic = 0.0
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


nr = 51
nz = 51
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

dom_fine = Domain(2*nr, 2*nz, rmax, zmax)

# --------------------- Set up initial contour ϕ0

# ----------------------- COmbine into a dict

# To do multiple combinations, include a vector here, per DrWatson
all_configs = Dict(
    "T_params" => T_params,
    "sim_dt" => sim_dt,
    # "sim_dt" => [sim_dt*0.5, sim_dt],
    # "dom" => dom,
    "dom" => [dom, dom_fine],
    "ϕ0type" => ϕ0type,
)


# -------------- Use DrWatson to run all simulations everything

list_configs = dict_list(all_configs)

pol_kwargs = (filename=hash, prefix="simdat", verbose=false, tag=true)

println("Starting $(length(list_configs)) simulations")
for conf in list_configs
    # println(conf["dom"].nr)
    @time simres, simdatfile = produce_or_load(sim_from_dict, conf, datadir("sims", "testing"); pol_kwargs...)

    # display(simres)
    casename = "sim_$(hash(conf))"
    sumplot = summaryplot(simres, conf)
    savefig(plotsdir("summary_"*casename*".svg"))
    resultsanim(simres, conf, casename)

end
println("Finished running (or loading) simulations")
# @tagsave

# println
