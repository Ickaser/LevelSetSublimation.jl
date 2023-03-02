export take_time_step, multistep
export sim_from_dict, sim_from_dict_ODE

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

"""
    ϕevol_RHS!(dϕ_flat, ϕ, p, t)
"""
function ϕevol_RHS!(dϕ_flat, ϕ, p, t)
    dom = p[1]
    T_params = p[2]
    dϕ = reshape(dϕ_flat, dom.nr, dom.nz)
    ϕ = reshape(ϕ, dom.nr, dom.nz)
    # debug
    if isnan(sum(ϕ))
        ϕ[isnan.(ϕ)] .= 1
        @info "Nan in ϕ"
    end

    T = solve_T(ϕ, dom, T_params)
    vf = extrap_v_fastmarch(ϕ, T, dom, T_params)
    vr = @view vf[:,:,1]
    vz = @view vf[:,:,2]
    

    for iz in 1:dom.nz, ir in 1:dom.nr
        irs = max.(1, min.(dom.nr, ir-3:ir+3)) # Pad with boundary values
        izs = max.(1, min.(dom.nz, iz-3:iz+3))
        dϕdr_l, dϕdr_r = wenodiffs_local(ϕ[irs, iz]..., dom.dr)
        dϕdz_l, dϕdz_r = wenodiffs_local(ϕ[ir, izs]..., dom.dz)
        rcomp = vr[ir,iz] > 0 ? dϕdr_l*vr[ir,iz] : dϕdr_r*vr[ir,iz]
        zcomp = vz[ir,iz] > 0 ? dϕdz_l*vz[ir,iz] : dϕdz_r*vz[ir,iz]
        dϕ[ir, iz] = -rcomp - zcomp
        if isnan(dϕ[ir, iz])
            println("Here be NaN: t=$t, ir=$ir, iz=$iz")
            @show vr[ir,iz] vz[ir,iz]
        end
    end
    return nothing
end


function reinit_wrap(integ)
    dom = integ.p[1]
    ϕ = reshape(integ.u, dom.nr, dom.nz)
    r_pre = get_subf_r(ϕ, dom)
    z_pre = get_subf_z(ϕ, dom)
    # vol_pre = sum(reshape(dom.rgrid, :, 1) .* (ϕ .<= 0) ) * dom.dr * dom.dz * 2π
    reinitialize_ϕ!(ϕ, dom, alg=BS3()) 
    r_post = get_subf_r(ϕ, dom)
    z_post = get_subf_z(ϕ, dom)
    rmove = r_post - r_pre
    zmove = z_post - z_pre
    # vol_post = sum(reshape(dom.rgrid, :, 1) .* (ϕ .<= 0) ) * dom.dr * dom.dz * 2π
    # volchange = vol_post - vol_pre
    @info "Reinit at t=$(integ.t), "
    @debug "com of interface moved:" rmove, zmove
    # @debug "vol change:" volchange
end

function next_reinit_time(integ)
    dom = integ.p[1]
    dϕ = similar(integ.u)
    ϕevol_RHS!(dϕ, integ.u, integ.p, integ.t)
    max_dϕdt = maximum(abs.(dϕ))
    mindx = min(dom.dr, dom.dz)
    domfrac = 0.15
    minlen = min(domfrac*dom.rmax , domfrac*dom.zmax)
    # Reinit next when interface should have moved across 10% of domain 
    dt = max_dϕdt / minlen
    return integ.t + dt
end

function sim_from_dict_ODE(fullconfig;tf = 100, reinit_time = 0.5)

    # ------------------- Get simulation parameters

    @unpack T_params, ϕ0type, dom, sim_dt = fullconfig

    # println("Inside: dom.nr = $(dom.nr)")
    ϕ0 = make_ϕ0(ϕ0type, dom)

    # reinitialize_ϕ!(ϕ0, dom, 1.0) # Don't reinit because callback handles this

    ϕ0_flat = reshape(ϕ0, :)

    # ---- Set up for ODEProblem
    prob_p = (dom, T_params)
    tspan = (0, tf)
    prob = ODEProblem(ϕevol_RHS!, ϕ0_flat, tspan, prob_p)

    # --- Set up reinitialization callback

    # cb1 = PeriodicCallback(reinit_wrap, reinit_time, initial_affect=true)
    # cb1 = PeriodicCallback(x->println("called back"), reinit_time, initial_affect=true)
    cb1 = IterativeCallback(next_reinit_time, reinit_wrap,  initial_affect = true)

    # --- Set up simulation end callback

    cond(u, t, integ) = minimum(u)
    cb2 = ContinuousCallback(cond, terminate!)

    cbs = CallbackSet(cb1, cb2)

    # --- Solve

    # sol = solve(prob, SSPRK43(), callback=cbs; )
    sol = solve(prob, Tsit5(), callback=cbs; )
    # sol = solve(prob, Tsit5(), callback=cb2; ) # No reinit
end



