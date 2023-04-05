export take_time_step, multistep
export sim_from_dict_old, sim_from_dict

export ϕevol_RHS

# --------- Convenience functions that need a home

function current_state(u, dom)
    ϕ = reshape(u[1:dom.ntot], size(dom))
    Tf = u[dom.ntot+1]
    return ϕ, Tf
end

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
    dϕ = reshape((@view du[1:ntot]), dom.nr, dom.nz)
    ϕ, Tf = current_state(u, dom)
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
    # @info "after vf" params[:p_sub]
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
    # p1 = heat(vr, dom)
    # p2 = heat(vz, dom)
    # display(plot(p1, p2))
    # display(heat(dϕ, dom))
    # if minimum(dϕ) < 0
    #     @info "negative sublimation" findall(dϕ .<0)
    # end
    # @info "eval derivatives" Tf p_sub-params[:p_ch] extrema(dϕ) extrema(T) extrema(p)
    return nothing
end


"""
    ϕevol_RHS(u, dom::Domain, params)
    ϕevol_RHS(u, config)
    
Compute the time derivative of `u` with given parameters.

`u` has `dom.ntot` entries for `ϕ` and one for `Tf`.

Wraps a call on `ϕevol_RHS!`, for convenience in debugging and elsewhere that efficiency is less important
"""
function ϕevol_RHS(u, dom::Domain, params)
    p = (dom, params)
    du = similar(u)
    # dϕ = zeros(dom.nr, dom.nz)

    # dϕ_flat = reshape(dϕ, :)
    # u = similar(ϕ, dom.ntot+1)
    # u[1:dom.ntot] .= reshape(ϕ, :)
    # u[dom.ntot+1] = params[:Tf]
    ϕevol_RHS!(du, u, p, 0.0)
    return du
end
function ϕevol_RHS(u, config)
    ϕevol_RHS(u, config[:dom], config[:params])
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
    ϕevol_RHS!(du, integ.u, integ.p, integ.t)

    # The main region of concern is the frozen region near interface
    # Find the largest value of dϕdt in that region
    ϕ = reshape(integ.u[1:dom.ntot], dom.nr, dom.nz)
    dϕ = reshape(du[1:dom.ntot], dom.nr, dom.nz)
    B = identify_B(ϕ, dom)
    B⁻ = B .& (ϕ .<= 0)
    max_dϕdt = maximum(abs.(dϕ[B⁻]))

    # Reinit next when interface should have moved across half of band around interface 
    domfrac = min(0.5 * dom.bwfrac, 0.1) # Minimum of 0.5 of the band, or 0.1 of domain size.
    minlen = min(domfrac*dom.rmax , domfrac*dom.zmax, integ.t*max_dϕdt + dom.dz)  # Also: at early times, do more often
    dt = minlen / max_dϕdt 
    # dt = 0.5 * minlen / max_dϕdt 
    @info "Reinit at t=$(integ.t), dt=$dt" minlen extrema(dϕ[B⁻])#, next at t=$(integ.t+dt)" 
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
- `Tf0`, an initial ice temperature
- `cparams`, which in turn has fields
    - `Q_gl`, `Q_sh` : heat flux from glass and shelf, respectively
    - `Q_ic`, `Q_ck` : volumetric heating in ice and cake, respectively
    - `k`: thermal conductivity of cake
    - `ρf`: density of ice
    - `Cpf`: heat capacity of ice
    - `ΔH` : heat of sublimation of ice
    - `p_ch` : pressure at top of cake
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

    @unpack cparams, ϕ0type, dom = fullconfig

    params = deepcopy(cparams)

    ϕ0 = make_ϕ0(ϕ0type, dom)
    reinitialize_ϕ_HCR!(ϕ0, dom, maxsteps=1000) # Don't reinit if using IterativeCallback
    ϕ0_flat = reshape(ϕ0, :)
    
    # Full array of starting state variables
    u0 = similar(ϕ0_flat, dom.ntot+1) # Add 1 to length: Tf
    u0[1:dom.ntot] .= ϕ0_flat
    u0[dom.ntot+1] = Tf = params[:Tf]
    p_sub = calc_psub(Tf)
    params[:p_sub] = p_sub
    @info "Initial Tf:" Tf p_sub

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
    sol = solve(prob, SSPRK43(), callback=cbs; ) # Adaptive timestepping
    # sol = solve(prob, SSPRK33(), dt=1e-4, callback=cbs; )
    # sol = solve(prob, Tsit5(), callback=cbs; )
    # sol = solve(prob, Tsit5(), callback=cb2; ) # No reinit
    return Dict("ϕsol"=>sol)
end



