using DrWatson
@quickactivate
using SparseArrays

using StaticArrays
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

println("Packages and code loaded")
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

# --------------------- Get thermal parameters

# T_params = make_artificial_params()

Q_gl = 1.0
Q_sh = 1.0
Q_ic = 1.0
Q_ck = 0.0
k = 1.0
Tf = 250.0
ΔH = 1.0
# ΔHsub = 678.0 # u"cal/g"
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

casename = "box_all"

# --------------------- Choose ϕ0

# Ellipse
# ϕ0 = [r^2 + 2z^2 - 1.5 for r in dom.rgrid, z in dom.zgrid]
# initshape = "ellipse"

# Separated ellipse
# ϕ0 = [1.5r^2 + 6(z-0.5)^2 - 1.0 for r in dom.rgrid, z in dom.zgrid]
# initshape = "separated"

# Circle
# ϕ0 = [1.1r^2 + 1.1z^2 - 1.0 for r in dom.rgrid, z in dom.zgrid]
# initshape = "circle"

# Box
# ϕ0 = [-1.0 for r in dom.rgrid, z in dom.zgrid]
# ϕ0[end,:] .= 1.0
# ϕ0[:,end] .= 1.0
# initshape = "box"

# Line
# ϕ0 = [0.5r + z - 0.9 + 0.001 for r in dom.rgrid, z in dom.zgrid]
# ϕ0 = [r + z - 1.5 + 0.001 for r in dom.rgrid, z in dom.zgrid]
# initshape = "cone"

# Flat on top
# ϕ0 = [ z - 0.999  for r in dom.rgrid, z in dom.zgrid]
# initshape = "flat"

# Cylinder
# ϕ0 = [r - 0.999 for r in dom.rgrid, z in dom.zgrid]

# Box
ϕ0 = [max(r,z)-0.99 for r in dom.rgrid, z in dom.zgrid]

# ----------------- Plot ϕ0

p = heat(ϕ0, dom)
markfront(ϕ0, dom)
plot_contour(ϕ0,dom, c=:white)
plot_contour(reinitialize_ϕ(ϕ0, dom, 1.0), dom, c=:black)
# display(p)

println("ϕ0 set up")

# ------------- Set up for simulation

" Stuff everything in a global function. May need to rethink this "
function take_time_step(Ti, ϕi, params, dt=1.0)
    # Precompute for velocity
    Qice = compute_Qice(ϕi, dom, params)
    icesurf = compute_icesurf(ϕi, dom, params)
    Qice_surf = Qice / icesurf
    
    prop_t = 1.0
    # # if minimum(ϕip1) < -min(5dr, 5dz) # When not much ice left, repair function needs more time
    # if minimum(ϕi) < -min(5dr, 5dz) || maximum(ϕi) > max(20dr, 20dz) # When contour is far from some regions, needs more time
    #     prop_t *= 2
    # end

    Bf = identify_B(ϕ, dom)
    
    frontfunc(ir, iz) = compute_frontvel_withT(Ti, ϕi, ir, iz, dom, params, Qice_surf)
    vf = vector_extrap_from_front(ϕi, Bf, frontfunc, dom, prop_t)
    
    ϕip1 = advect_ϕ(ϕi, vf, dom, dt)
    
    # repair_t = 1.0
    # # if minimum(ϕip1) < -min(5dr, 5dz) # When not much ice left, repair function needs more time
    # if minimum(ϕip1) < -min(5dr, 5dz) || maximum(ϕip1) > max(20dr, 20dz) # When contour is far from some regions, needs more time
    #     repair_t *= 2
    # end
    
    reinitialize_ϕ!(ϕip1, dom, prop_t)
    
    Tip1 = solve_T(ϕip1, dom, params)
    return Tip1, ϕip1
end
function multistep(n, dt, T0, ϕ0)
    Ti = copy(T0)
    ϕi = copy(ϕ0)
    for i in 1:n
        Ti, ϕi = take_time_step(Ti, ϕi, T_params, dt)
    end
    return Ti, ϕi
end

reinitialize_ϕ!(ϕ0, dom, 1.0)
T0 = solve_T(ϕ0, dom, T_params)

maxsteps = 100

full_T = fill(1.0, (maxsteps+1, dom.nr, dom.nz))
full_ϕ = fill(1.0, (maxsteps+1, dom.nr, dom.nz))
plots = []
Ti = copy(T0)
ϕi = copy(ϕ0)
full_T[1,:,:] .= Ti
full_ϕ[1,:,:] .= ϕi

sim_dt = 1.0

println("Ready to start simulation")

# --------------- Run simulation

for i in 1:maxsteps
    @time global Ti, ϕi = take_time_step(Ti, ϕi, T_params, sim_dt)
    # @time Ti, ϕi = multistep(3, 10.0, Ti, ϕi)
    
    # Store solutions for later plotting
    global full_T[i+1,:,:] .= round.(Ti, sigdigits=10) # Get rid of numerical noise
    global full_ϕ[i+1,:,:] .= ϕi

    # println("Completed: $i")
    
    if minimum(ϕi) > 0
        println("Sublimation finished")
        global full_T = full_T[1:i+1, :,:]
        global full_ϕ = full_ϕ[1:i+1, :,:]
        break
    end
end

# ------------------- Plot results

nt = size(full_T, 1) 
plots = []
if nt >= 6
    frames = round.(Int, range(1, nt, length=6))
else
    frames = 1:nt
end

for f in frames
    # Default: plot heat
    local p = plot(aspect_ratio=:equal)
    plot_cylheat(full_T[f,:,:], dom)
    plot_cylcont(full_ϕ[f,:,:], dom, c=:white)
    plot!(title="timestep=$(f-1)")
    
    # Debug: plot heat and level set
    # p1 = heat(full_T[f,:,:])
    # plot_contour(full_ϕ[f,:,:], c=:white)
    # # p2 = plot(aspect_ratio=:equal, xlim=(0,1), ylim=(0,1))
    # p2 = heat(full_ϕ[f,:,:])
    # plot_contour(full_ϕ[f,:,:], c=:white)
    # p = plot(p1, p2)

    # p = plot_frontvel(full_ϕ[f,:,:], full_T[f,:,:])

    # p2 = heat(full_ϕ[f,:,:])
    # markfront(full_ϕ[f,:,:])
    # plot_contour(full_ϕ[f,:,:], c=:black)
    # p = plot(p1, p2)

    push!(plots, p)
end

# ---------------- Save plots
bigplot = plot(plots..., size=(500*2, 200*3), layout=(3,2))
savefig(plotsdir("$(casename)_evolution.svg"))
display(bigplot)

freshplot()
plot!(size=(800,500))
anim = @animate for i ∈ 1:nt
    freshplot()
    plot!(title="timestep=$(i-1)")
    plot_cylheat(full_T[i,:,:])
    plot_cylcont(full_ϕ[i,:,:])
end

fps = ceil(Int, nt/2)  # nt/x makes x-second long animation
gif(anim, plotsdir("$(casename)_evol.gif"), fps=fps)
# println(full_T[end,:,:])