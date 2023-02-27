export take_time_step, multistep
export sim_from_dict

# using DrWatson
# @quickactivate
# using SparseArrays

# # using StaticArrays
# # using DomainSets
# using Contour
# # using Interpolations
# using Plots
# using DifferentialEquations


# include(srcdir("structs.jl"))
# include(srcdir("levelset_plots.jl"))
# include(srcdir("levelset_reinit.jl"))
# include(srcdir("levelset_advect.jl"))
# include(srcdir("solve_T.jl"))
# include(srcdir("coupled_motion.jl"))

# println("Packages and code loaded")

"""
    take_time_step(Ti, ϕi, dom::Domain, params, dt=1.0)

Return ``T_{i+1}, ϕ_{i+1}``, given ``T_i, ϕ_i``.

A fair amount of logic happens inside here. For example, the CFL condition is enforced 
inside each of the separate time integrations. However, under high heating, it is
possible to advect the interface past the edge of the band where the level set function 
is maintained.
"""
function take_time_step(Ti, ϕi, dom::Domain, params, dt=1.0)
    
    prop_t = 1.0 # Time step for reinitialization and similar steps

    # --------- Extrapolation by PDE
    # vf = extrap_v_pde(ϕi, Ti, dom, params)

    # ---------- Extrapolation by fast marching (more recent)
    vf = extrap_v_fastmarch(ϕi, Ti, dom, params)

    # ----------
    
    ϕip1 = advect_ϕ(ϕi, vf, dom, dt)
    
    
    if sum(ϕip1 .<= 0) == 0 # No ice cells left
        # ϕip1[1,1] = 0 # Artifically add a tiny amount of ice
        Tip1 = solve_T(ϕip1, dom, params)
        return Tip1, ϕip1
    end



    
    reinitialize_ϕ!(ϕip1, dom, prop_t)
    
    Tip1 = solve_T(ϕip1, dom, params)
    return Tip1, ϕip1
end
"""
    multistep(n, dt, T0, ϕ0, dom::Domain, T_params)

Return `Ti` and `ϕi` after `n` timesteps of `dt`.

I anticipate this being useful mostly if the timestep for stability is small and don't need to store every time step.
"""
function multistep(n, dt, T0, ϕ0, dom::Domain, T_params)
    Ti = copy(T0)
    ϕi = copy(ϕ0)
    for i in 1:n
        Ti, ϕi = take_time_step(Ti, ϕi, dom, T_params, dt)
    end
    return Ti, ϕi
end

# ----------- Actual start of a simulation
"""
    sim_from_dict(fullconfig; maxsteps=1000)

Run a simulation from `fullconfig`, returning T and ϕ at each time step.

Implicitly, have maximum number of time steps of 1000.
`fullconfig` should have the following fields:
- `ϕ0type`, types listed for [`make_ϕ0`](@ref)
- `dom`, an instance of [`Domain`](@ref)
- `sim_dt`, simulation time step (used in advection only)
- `T_params`, which in turn has fields
    - `Tf`: ice temperature
    - `Q_gl`, `Q_sh` : heat flux from glass and shelf, respectively
    - `Q_ic`, `Q_ck` : volumetric heating in ice and cake, respectively
    - `k`: thermal conductivity of cake
"""
function sim_from_dict(fullconfig; maxsteps=1000)

    # ------------------- Get simulation parameters

    @unpack T_params, ϕ0type, dom, sim_dt = fullconfig

    # println("Inside: dom.nr = $(dom.nr)")
    ϕ0 = make_ϕ0(ϕ0type, dom)

    # ---------------------- Set up first step
    reinitialize_ϕ!(ϕ0, dom, 1.0)
    T0 = solve_T(ϕ0, dom, T_params)

    # maxsteps = 1000

    full_T = fill(1.0, (maxsteps+1, dom.nr, dom.nz))
    full_ϕ = fill(1.0, (maxsteps+1, dom.nr, dom.nz))
    Ti = copy(T0)
    ϕi = copy(ϕ0)
    full_T[1,:,:] .= Ti
    full_ϕ[1,:,:] .= ϕi


    # println("Ready to start simulation")

    # --------------- Run simulation

    for i in 1:maxsteps
        Ti, ϕi = take_time_step(Ti, ϕi, dom, T_params, sim_dt)
        # @time Ti, ϕi = multistep(3, 10.0, Ti, ϕi)
        
        # Store solutions for later plotting
        full_T[i+1,:,:] .= round.(Ti, sigdigits=10) # Get rid of numerical noise
        full_ϕ[i+1,:,:] .= ϕi

        # println("Completed: $i")
        
        if minimum(ϕi) > 0
            println("Sublimation finished after $i timesteps")
            full_T = full_T[1:i+1, :,:]
            full_ϕ = full_ϕ[1:i+1, :,:]
            break
        end
    end
    @strdict full_T full_ϕ
end