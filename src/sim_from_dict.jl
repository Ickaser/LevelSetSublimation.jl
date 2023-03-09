export take_time_step, multistep
export sim_from_dict_old, sim_from_dict

export ϕevol_RHS


# ---------- Fully adaptive time stepping functions

"""
    ϕevol_RHS!(dϕ_flat, ϕ_flat, p, t)
"""
function ϕevol_RHS!(dϕ_flat, ϕ_flat, p, t)
    dom = p[1]
    T_params = p[2]
    dϕ = reshape(dϕ_flat, dom.nr, dom.nz)
    ϕ = reshape(ϕ_flat, dom.nr, dom.nz)
    # debug
    # if isnan(sum(ϕ))
    #     ϕ[isnan.(ϕ)] .= 1
    #     @info "NaN in ϕ"
    # end

    T = solve_T(ϕ, dom, T_params)
    vf = extrap_v_fastmarch(ϕ, T, dom, T_params)
    vr = @view vf[:,:,1]
    vz = @view vf[:,:,2]
    

    indmin = CI(1, 1)
    indmax = CI(dom.nr, dom.nz)
    rshift = [CI(i, 0) for i in -3:3]
    zshift = [CI(0, i) for i in -3:3]
    for ind in CartesianIndices(ϕ)
    # for iz in 1:dom.nz, ir in 1:dom.nr
        # irs = max.(1, min.(dom.nr, ir-3:ir+3)) # Pad with boundary values
        # izs = max.(1, min.(dom.nz, iz-3:iz+3))
        rst = max.(min.([ind].+rshift, [indmax]), [indmin]) # Pad beyond boundary with boundary
        zst = max.(min.([ind].+zshift, [indmax]), [indmin]) # Pad beyond boundary with boundary
        dϕdr_w, dϕdr_e = wenodiffs_local(ϕ[rst]..., dom.dr)
        dϕdz_s, dϕdz_n = wenodiffs_local(ϕ[zst]..., dom.dz)
        # Boundary cases: use internal derivative
        if rst[5] == ind # Right boundary
            dϕdr = dϕdr_w
        elseif rst[3] == ind # Left boundary
            dϕdr = dϕdr_e
        else
            dϕdr = (vr[ind] > 0 ? dϕdr_w : dϕdr_e)
        end
        # Boundary cases: use internal derivative
        if zst[5] == ind # Top boundary
            dϕdz = dϕdz_s
        elseif zst[3] == ind # Bottom boundary
            dϕdz = dϕdz_n
        else
            dϕdz = (vz[ind] > 0 ? dϕdz_s : dϕdz_n)
        end

        rcomp = dϕdr*vr[ind]
        zcomp = dϕdz*vz[ind]
        dϕ[ind] = -rcomp - zcomp
        # if isnan(dϕ[ir, iz])
        #     println("Here be NaN: t=$t, ir=$ir, iz=$iz")
        #     @show vr[ir,iz] vz[ir,iz]
        # end
        # @info "check" ind rcomp zcomp vz[ind] dϕdz_s dϕdz_n
    end
    return nothing
end


"""
    ϕevol_RHS(ϕ, dom, T_params)
    
Wraps a call on `ϕevol_RHS!`, for convenience in debugging
"""
function ϕevol_RHS(ϕ, dom, T_params)
    p = (dom, T_params)
    dϕ = zeros(dom.nr, dom.nz)
    dϕ_flat = reshape(dϕ, :)
    ϕ_flat = reshape(ϕ, :)
    ϕevol_RHS!(dϕ_flat, ϕ_flat, p, 0.0)
    return dϕ
end

function reinit_wrap(integ)
    dom = integ.p[1]
    ϕ = reshape(integ.u, dom.nr, dom.nz)
    # r_pre = get_subf_r(ϕ, dom)
    # z_pre = get_subf_z(ϕ, dom)
    # vol_pre = sum(reshape(dom.rgrid, :, 1) .* (ϕ .<= 0) ) * dom.dr * dom.dz * 2π
    reinitialize_ϕ!(ϕ, dom, alg=BS3()) 
    # r_post = get_subf_r(ϕ, dom)
    # z_post = get_subf_z(ϕ, dom)
    # rmove = r_post - r_pre
    # zmove = z_post - z_pre
    # vol_post = sum(reshape(dom.rgrid, :, 1) .* (ϕ .<= 0) ) * dom.dr * dom.dz * 2π
    # volchange = vol_post - vol_pre
    # @info "Reinit at t=$(integ.t), "
    # @debug "com of interface moved:" rmove, zmove
    # @debug "vol change:" volchange
end

function next_reinit_time(integ)
    dom = integ.p[1]
    dϕ_flat = similar(integ.u)
    ϕevol_RHS!(dϕ_flat, integ.u, integ.p, integ.t)

    B = identify_B(reshape(integ.u, dom.nr, dom.nz), dom)
    dϕ = reshape(dϕ_flat, dom.nr, dom.nz)
    max_dϕdt = maximum(abs.(dϕ[B]))

    # max_dϕdt = maximum(abs.(dϕ_flat))

    # Reinit next when interface should have moved across half of band around interface 
    domfrac = min(0.6 * dom.bwfrac, 0.1) # Minimum of 0.6 of the band, or 0.25 of domain size.
    minlen = min(domfrac*dom.rmax , domfrac*dom.zmax,) 
    dt = minlen / max_dϕdt 
    @info "Reinit at t=$(integ.t)"#, next at t=$(integ.t+dt)" 
    return integ.t + dt
end

function sim_from_dict(fullconfig;tf = 100)

    # ------------------- Get simulation parameters

    @unpack T_params, ϕ0type, dom, sim_dt = fullconfig

    # println("Inside: dom.nr = $(dom.nr)")
    ϕ0 = make_ϕ0(ϕ0type, dom)

    reinitialize_ϕ!(ϕ0, dom, 100.0) # Don't reinit because callback handles this

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
    return Dict("ϕsol"=>sol)
end



# ------------ Semi-manual time stepping
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
function sim_from_dict_old(fullconfig; maxsteps=1000)

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
