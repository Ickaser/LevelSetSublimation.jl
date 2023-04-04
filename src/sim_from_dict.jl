export take_time_step, multistep
export sim_from_dict_old, sim_from_dict

export ϕevol_RHS


# ---------- Fully adaptive time stepping functions

"""
    ϕevol_RHS!(du, u, p, t)

Evaluate local time rate of change for `u` and put results in `du`.

`u` and `du` are both structured as follows:
First `dom.ntot` values are `ϕ`, reshaped; `dom.ntot+1` index is frozen temperature `Tf`

Parameters `p` assumed to be `(dom::Domain, params)`
This is a right-hand-side for ∂ₜϕ = -v⋅∇ϕ, where v = `(vr, vz)` is evaluated by computing and extrapolating front velocity
using `compute_frontvel_mass`.
"""
function ϕevol_RHS!(du, u, p, t)
    dom = p[1]
    params = p[2]
    ntot = dom.ntot
    dϕ = reshape(du[1:ntot], dom.nr, dom.nz)
    ϕ = reshape(u[1:ntot], dom.nr, dom.nz)

    Tf = u[ntot+1]
    p_sub = calc_psub(Tf) # Returned as Pa
    @pack! params = Tf, p_sub # Update current temperature and sublimation pressure 
    # debug
    # if isnan(sum(ϕ))
    #     ϕ[isnan.(ϕ)] .= 1
    #     @info "NaN in ϕ"
    # end


    # Plug current value of Tf into params, where it then gets passed around


    T = solve_T(ϕ, dom, params)
    p = solve_p(ϕ, T, dom, params)
    vf = extrap_v_fastmarch(ϕ, T, p, dom, params)
    vr = @view vf[:,:,1]
    vz = @view vf[:,:,2]
    
    Qice = compute_Qice(ϕ, T, p, dom, params)
    @unpack ρf, Cpf = params
    dTdt = Qice / ρf / Cpf / compute_icevol(ϕ, dom)
    du[ntot+1] = dTdt

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
    vtot = hypot.(vr, vz)
    # display(heat(dϕ, dom))
    @info "eval derivatives" Tf p_sub-params[:p_ch] extrema(dϕ) extrema(vtot)
    return nothing
end


"""
    ϕevol_RHS(ϕ, dom::Domain, params)
    ϕevol_RHS(ϕ, config)
    
Compute the time derivative of `ϕ` with given parameters.

Wraps a call on `ϕevol_RHS!`, for convenience in debugging and elsewhere that efficiency is less important
"""
function ϕevol_RHS(ϕ, dom::Domain, params)
    p = (dom, params)
    du = zeros(dom.ntot+1)
    # dϕ = zeros(dom.nr, dom.nz)

    # dϕ_flat = reshape(dϕ, :)
    u = similar(ϕ, dom.ntot+1)
    u[1:dom.ntot] .= reshape(ϕ, :)
    u[dom.ntot+1] = params[:Tf]
    ϕevol_RHS!(du, u, p, 0.0)
    return du
end
function ϕevol_RHS(ϕ, config)
    ϕevol_RHS(ϕ, config[:dom], config[:params])
end

"""
    reinit_wrap(integ)

Thin wrapper to reinitialize the state of the level set function.

Calls `reinitialize_ϕ!(ϕ, dom)`, so uses the default reinitialization setup.
Used internally in an `IterativeCallback`, as implemented in `DiffEqCallbacks`.
"""
function reinit_wrap(integ; verbose=false)
    if verbose
        @info "Reinit at t=$(integ.t)"
    end
dom = integ.p[1]
ϕ = reshape((@view integ.u[1:dom.ntot]), dom.nr, dom.nz)
# reinitialize_ϕ!(ϕ, dom) 
reinitialize_ϕ_HCR!(ϕ, dom, tol=1e-6) 
end

"""
next_reinit_time(integ)

Compute the next time reinitialization should be necessary, given the current integrator state.
Used internally in an `IterativeCallback`, as implemented in `DiffEqCallbacks`.
"""
function next_reinit_time(integ)
    dom = integ.p[1]
    du = similar(integ.u)
    du_0 = copy(du)
    ϕevol_RHS!(du, integ.u, integ.p, integ.t)
    
    @info "problems?" extrema(du_0 .- du)

    # B = identify_B(reshape(integ.u, dom.nr, dom.nz), dom)
    # dϕ = reshape(dϕ_flat, dom.nr, dom.nz)
    # max_dϕdt = maximum(abs.(dϕ[B]))

    # The main region of concern is the frozen region near interface
    # Find the largest value of dϕdt in that region
    ϕ = reshape(integ.u[1:dom.ntot], dom.nr, dom.nz)
    dϕ = reshape(du[1:dom.ntot], dom.nr, dom.nz)
    B = identify_B(ϕ, dom)
    B⁻ = B .& (ϕ .<= 0)
    max_dϕdt = maximum(abs.(dϕ[B⁻]))

    # max_dϕdt = maximum(abs.(dϕ_flat))

    # Reinit next when interface should have moved across half of band around interface 
    domfrac = min(0.5 * dom.bwfrac, 0.1) # Minimum of 0.5 of the band, or 0.1 of domain size.
    minlen = min(domfrac*dom.rmax , domfrac*dom.zmax, integ.t*max_dϕdt + dom.dz)  # Also: at early times, do more often
    dt = minlen / max_dϕdt 
    # dt = 0.5 * minlen / max_dϕdt 
    @info "Reinit at t=$(integ.t), dt=$dt" extrema(du) extrema(dϕ[B⁻])#, next at t=$(integ.t+dt)" 
    return integ.t + dt
end

function cond_reinit(u, t, integ)
    dom = integ.p[1]
    ϕ = reshape(u, size(dom))
    err = sdf_err_L1(ϕ, dom)
    tol = 1e-5 # Roughly dx^4
    @info "error from sdf" err-tol
    return err-tol
end

"""
    sim_from_dict(fullconfig; tf=100, verbose=false)

Given a simulation configuration `fullconfig`, run a simulation

Maximum simulation time is specified by `tf`.
`verbose=true` will put out some info messages about simulation progress, i.e. at each reinitialization.

`fullconfig` should have the following fields:
- `ϕ0type`, types listed for [`make_ϕ0`](@ref)
- `dom`, an instance of [`Domain`](@ref)
- `params`, which in turn has fields
    - `Tf`: ice temperature
    - `Q_gl`, `Q_sh` : heat flux from glass and shelf, respectively
    - `Q_ic`, `Q_ck` : volumetric heating in ice and cake, respectively
    - `k`: thermal conductivity of cake
    - `ϵ` : porosity of porous medium
    - `l` : dusty gas model constant: characteristic length for Knudsen diffusion
    - `κ` : dusty gas model constant: length^2 corresponding loosely to Darcy's Law permeability
    - `R` : universal gas constant, with appropriate units 
    - `Mw`: molecular weight of species (water), with appropriate units
    - `μ` : dynamic viscosity of species (water), with appropriate units

If you are getting a warning about instability, it can often be fixed by tinkering with the reinitialization behavior.
That shouldn't be true but it seems like it is.


"""
function sim_from_dict(fullconfig; tf=100, verbose=false)

    # ------------------- Get simulation parameters

    @unpack params, ϕ0type, dom = fullconfig

    params = deepcopy(params)

    ϕ0 = make_ϕ0(ϕ0type, dom)
    reinitialize_ϕ_HCR!(ϕ0, dom, maxsteps=1000) # Don't reinit if using IterativeCallback
    ϕ0_flat = reshape(ϕ0, :)
    
    # Full array of starting state variables
    u0 = similar(ϕ0_flat, dom.ntot+1) # Add 1 to length: Tf
    u0[1:dom.ntot] .= ϕ0_flat
    u0[dom.ntot+1] = params[:Tf]
    @info "Initial Tf:" u0[dom.ntot+1]

    # ---- Set up ODEProblem
    prob_p = (dom, params)
    tspan = (0, tf)
    prob = ODEProblem(ϕevol_RHS!, u0, tspan, prob_p)

    # --- Set up reinitialization callback

    # cb1 = PeriodicCallback(reinit_wrap, reinit_time, initial_affect=true)
    if verbose
        cb1 = IterativeCallback(next_reinit_time, x->reinit_wrap(x, verbose=true),  initial_affect = true)
    else
        cb1 = IterativeCallback(next_reinit_time, reinit_wrap,  initial_affect = true)
    end

    # cb1 = ContinuousCallback(cond_reinit, reinit_wrap)

    # --- Set up simulation end callback

    # When the minimum value of ϕ is 0, front has disappeared
    cond_end(u, t, integ) = minimum(u) 
    # ContinuousCallback gets thrown when `cond` evaluates to 0
    # `terminate!` ends the solve there
    cb2 = ContinuousCallback(cond_end, terminate!)

    # --- Put callbacks together
    cbs = CallbackSet(cb1, cb2)

    # --- Solve
    sol = solve(prob, SSPRK43(), callback=cbs; )
    # sol = solve(prob, Tsit5(), callback=cbs; )
    # sol = solve(prob, Tsit5(), callback=cb2; ) # No reinit
    return Dict("ϕsol"=>sol)
end



# ------------ Semi-manual time stepping. Obsolete, I think  -------------------------
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
    multistep(n, dt, T0, ϕ0, dom::Domain, params)

Return `Ti` and `ϕi` after `n` timesteps of `dt`.

I anticipate this being useful mostly if the timestep for stability is small and don't need to store every time step.
"""
function multistep(n, dt, T0, ϕ0, dom::Domain, params)
    Ti = copy(T0)
    ϕi = copy(ϕ0)
    for i in 1:n
        Ti, ϕi = take_time_step(Ti, ϕi, dom, params, dt)
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
- `params`, which in turn has fields
    - `Tf`: ice temperature
    - `Q_gl`, `Q_sh` : heat flux from glass and shelf, respectively
    - `Q_ic`, `Q_ck` : volumetric heating in ice and cake, respectively
    - `k`: thermal conductivity of cake
"""
function sim_from_dict_old(fullconfig; maxsteps=1000)

    # ------------------- Get simulation parameters

    @unpack params, ϕ0type, dom, sim_dt = fullconfig

    # println("Inside: dom.nr = $(dom.nr)")
    ϕ0 = make_ϕ0(ϕ0type, dom)

    # ---------------------- Set up first step
    reinitialize_ϕ!(ϕ0, dom, 1.0)
    T0 = solve_T(ϕ0, dom, params)

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
        Ti, ϕi = take_time_step(Ti, ϕi, dom, params, sim_dt)
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
