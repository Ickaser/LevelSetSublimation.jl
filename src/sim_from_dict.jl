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
    # Precompute for velocity
    Qice = compute_Qice(ϕi, dom, params)
    icesurf = compute_icesurf(ϕi, dom)
    Qice_surf = Qice / icesurf
    
    prop_t = 1.0 # Time step for reinitialization and similar steps
    # # if minimum(ϕip1) < -min(5dr, 5dz) # When not much ice left, repair function needs more time
    # if minimum(ϕi) < -min(5dr, 5dz) || maximum(ϕi) > max(20dr, 20dz) # When contour is far from some regions, needs more time
    #     prop_t *= 2
    # end

    Bf = identify_B(ϕi, dom)
    
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
"""
    multistep(n, dt, T0, ϕ0, dom::Domain)

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
function sim_from_dict(fullconfig)

    # ------------------- Get simulation parameters

    @unpack T_params, ϕ0type, dom, sim_dt = fullconfig

    # println("Inside: dom.nr = $(dom.nr)")
    ϕ0 = make_ϕ0(ϕ0type, dom)

    # ---------------------- Set up first step
    reinitialize_ϕ!(ϕ0, dom, 1.0)
    T0 = solve_T(ϕ0, dom, T_params)

    maxsteps = 1000

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