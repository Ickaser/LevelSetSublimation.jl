export calc_uTfTp_res, get_t_Tf, get_t_Tf_subflux, compare_lyopronto_res
export get_subf_z, get_subf_r, get_ϕ, get_SA
export get_eff_Rp
export virtual_thermocouple
"$(SIGNATURES)"
function get_t_Tf(sim)
    @unpack sol, dom, config = sim
    t = sol.t * u"s"
    Tf = map(sol.t) do ti
        Tfr = calc_Tf_res(ti, sim)
        return Tfr[1]*u"K"
    end
    return t, Tf
end

"$(SIGNATURES)"
function get_t_Tf_subflux(sim)
    @unpack sol, dom = sim
    t, Tf = get_t_Tf(sim)
    mdt = map(sol.t) do ti
        params = calc_params_at_t(ti, sim.config)
        uTfTp = calc_uTfTp_res(ti, sim)
        Tfi = uTfTp[2]
        md = compute_topmassflux(uTfTp..., dom, params) * u"kg/s"
        if sign(md) == -1
            @info "md=$md" ti Tfi calc_psub(Tfi)
        end
        md = max(zero(md), md)
        md
    end
    return t, Tf, mdt
end


"$(SIGNATURES)"
calc_fillvol(dom) = π*dom.rmax^2*dom.zmax*u"m^3"

"$(SIGNATURES)"
function get_eff_Rp(sim)
    @unpack sol, dom = sim
    t, Tf = get_t_Tf(sim)
    hf0 = dom.zmax*u"m"
    Ap = π*dom.rmax^2*u"m^2"
    fillvol = Ap*hf0
    h_md_A_t = map(sol.t) do ti
        params = calc_params_at_t(ti, sim.config)
        uTfTp = calc_uTfTp_res(ti, sim)
        ϕ = sim.sol(ti).ϕ
        Vf = compute_icevol_H(ϕ, dom)*u"m^3"
        Asub = compute_icesurf_δ(ϕ, dom)*u"m^2"
        hd = (1-Vf/fillvol)*hf0
        md = compute_topmassflux(uTfTp..., dom, params) * u"kg/s"
        if sign(md) == -1
            Tfi = uTfTp[2]
            @info "md=$md" ti Tfi calc_psub.(Tfi)
        end
        md = max(zero(md), md)
        [hd, md, Asub]
    end
    hd = [u[1] for u in h_md_A_t]
    md = [u[2] for u in h_md_A_t]
    Asub = [u[3] for u in h_md_A_t]
    psub = calc_psub.(Tf)
    pch = sim.config[:paramsd][3].pch.(sol.t*u"s")
    Rp = @. (psub - pch)*Ap/md
    return t, hd, Rp, Asub
end

"$(SIGNATURES)"
function compare_lyopronto_res(sim)
    t = uconvert.(u"hr", sim.sol.t .* u"s")
    return compare_lyopronto_res(t, sim)
end

"$(SIGNATURES)"
function compare_lyopronto_res(ts, sim)
    @unpack sol, dom, config = sim
    ts_ndim = ustrip.(u"s", ts)
    cyc = ts_ndim .< sol.t[end]
    if ~all(cyc)
        @warn "Some requested times are beyond the end of the simulation. Trimming to fit."
    end
    Tf = similar(ts_ndim[cyc])*u"K"
    md = similar(ts_ndim[cyc]).*u"kg/s"
    for (i, ti) in enumerate(ts_ndim[cyc])
        params = calc_params_at_t(ti, sim.config)
        uTfTp = calc_uTfTp_res(ti, sim)
        Tf[i] = uTfTp[3][1,1]*u"K"
        mdi = compute_topmassflux(uTfTp..., dom, params) * u"kg/s"
        if sign(mdi) == -1
            @info "md=$mdi" ti Tfi calc_psub(Tfi)
        end
        mdi = max(zero(mdi), mdi)
        md[i] = mdi
    end
    totvol = π*dom.rmax^2 * dom.zmax
    dryfrac = map(ts_ndim[cyc]) do ti
        1 - compute_icevol_H(sol(ti).ϕ, dom) / totvol
    end

    return ts[cyc], Tf, md, dryfrac
end

function LyoPronto.obj_expT(sim::SimResults, pdfit::PrimaryDryFit;
    tweight=1.0, verbose=false, Tvw_weight=1.0)
    (;sol, dom, config) = sim
    # if sol.retcode !== ReturnCode.Terminated || length(sol.u) <= 1
    #     verbose && @warn "ODE solve did not reach end of drying. Either parameters are bad, or tspan is not large enough." sol.retcode sol.prob.p.hf0 sol[end]
    #     return NaN
    # end
    solt = sol.t .* u"s"
    tmd = solt[end]
    nt = length(sol.t) - 1
    i_solstart = searchsortedfirst(pdfit.t, solt[begin]) 
    # Identify if the solution is pre-interpolated to the time points in pdfit.t

    # Compute temperature objective for all frozen temperatures
    ftrim = solt[begin] .< pdfit.t .< tmd
    tf_trim = pdfit.t[ftrim]
    
    Tfmd = map(tf_trim) do t
        calc_Tf_res(ustrip(u"s", t), sim)[1]*u"K"
    end
    # Sometimes the interpolation procedure of the solution produces wild temperatures, as in below absolute zero.
    # This bit replaces any subzero values with the previous positive temperature, and notifies that it happened.
    if any(Tfmd .< 0u"K")
        subzero = findall(Vector(Tfmd .< 0u"K"))
        Tfmd[subzero] .= Tfmd[subzero[1] - 1] 
        verbose && @info "bad interpolation" subzero Tfmd[subzero]
    end
    Tfobj = 0.0u"K^2"
    for (Tf, iend) in zip(pdfit.Tfs, pdfit.Tf_iend)
        trim = min(iend, length(Tfmd))
        Tfobj += sum(abs2, (Tf[i_solstart:trim] .- Tfmd[begin:trim-i_solstart+1]))/(trim-i_solstart+1)
    end
    if ismissing(pdfit.Tvws) # No vial wall temperatures
        Tvwobj = 0.0u"K^2"
    elseif ismissing(pdfit.Tvw_iend) # Only an endpoint temperature provided
        Tvwend = pdfit.Tvws
        Tvwobj = (sol[3, end]*u"K" - uconvert(u"K", Tvwend))^2
    else # Regular case of fitting to at least one full temperature series
        vwtrim = sol.t[begin]*u"hr" .< pdfit.t .< tmd
        tvw_trim = pdfit.t[vwtrim]
        Tvwmd = map(tvw_trim) do t# .- 273.15
            sol(ustrip(u"s", t)).Tvw*u"K"
        end
        # Compute temperature objective for all vial wall temperatures
        Tvwobj = 0.0u"K^2"
        for (Tvw, iend) in zip(pdfit.Tvws, pdfit.Tvw_iend)
            trim = min(iend, length(Tvwmd))
            Tvwobj += sum(abs2, (Tvw[i_solstart:trim] .- Tvwmd[begin:trim-i_solstart+1]))/(trim-i_solstart+1)
        end
    end
    tobj = if ismissing(pdfit.t_end) 
        # No drying time provided
        0.0u"hr^2"
    elseif pdfit.t_end isa Tuple 
        # See if is inside window and scale appropriately
        mid_t = (pdfit.t_end[1] + pdfit.t_end[2]) / 2.0
        if tmd < pdfit.t_end[1]
            (mid_t - tmd)^2
        elseif tmd > pdfit.t_end[2]
            (mid_t - tmd)^2
        else # Inside window, so no error
            0.0u"hr^2"
        end
    else 
        # Compare to a single drying time
        (pdfit.t_end - tmd)^2
    end
    verbose && @info "loss call" tmd tobj Tfobj Tvwobj 
    return ustrip(u"K^2", Tfobj + Tvw_weight*Tvwobj) + tweight*ustrip(u"hr^2", tobj)

end

"$(SIGNATURES)"
function gen_sumplot(config, var=:T, casename="test")
    pol_kwargs = (filename=hash, prefix="simdat", verbose=false, tag=true)
    @time simres, simdatfile = produce_or_load(sim_from_dict, config,
            datadir("sims", casename); pol_kwargs...)
    summaryplot(simres, config, heatvar=var, layout=(4,1))
end
"$(SIGNATURES)"
function gen_anim(config, var=:T, casename="test")
    pol_kwargs = (filename=hash, prefix="simdat", verbose=false, tag=true)
    @time simres, simdatfile = produce_or_load(sim_from_dict, config,
            datadir("sims", casename); pol_kwargs...)
    resultsanim(simres, config, casename, heatvar=var)
    return simres
end

"$(SIGNATURES)"
function calc_params_at_t(t, config)
    @unpack paramsd = config
    params = params_nondim_setup(paramsd)
    return (params[1], params[2], params[3](t))
end

"""
    $(SIGNATURES)

Interpolate a set of SavedValues of Tf.
Useful so that I use the same interpolation everywhere it shows up.
"""
function interp_saved_Tf(saved_Tf) 
    if length(saved_Tf.t) <= 2  
        return ConstantInterpolation(saved_Tf.saveval, saved_Tf.t, extrapolation=ExtrapolationType.Constant)
    elseif any(diff(saved_Tf.t) == 0)
        return LinearInterpolation(saved_Tf.saveval, saved_Tf.t)
    else
        return CubicSpline(saved_Tf.saveval, saved_Tf.t)
    end
end

"""
    $(SIGNATURES)

Get Tf at a given time point, regardless of simulation's time integration. 
"""
function calc_Tf_res(t, sim)
    @unpack sol, dom = sim
    if sol isa CombinedSolution
        if t <= sol.tsplit
            return sol.sol1(t).Tf
        else
            sol.Tf2(t)
        end
    elseif !isnothing(sim.Tf)
        return sim.Tf(t)
    else # Used a DAE or implicit solve
        return sol(t).Tf
    end
end

"""
    $(SIGNATURES)
"""
function calc_uTfTp_res(t, sim)
    @unpack sol, dom, config = sim
    params = calc_params_at_t(t, config)
    u = sol(t)
    Tf = calc_Tf_res(t, sim)
    identify_dry(dom, u.ϕ)
    T = solve_T(u, Tf, dom, params)
    p = solve_p(u, Tf, T, dom, params)
    return u, Tf, T, p
end

"""
    $(SIGNATURES)
"""
function virtual_thermocouple(sim::SimResults) 
    virtual_thermocouple([(0, 0)], sim)
end
function virtual_thermocouple(locs, sim::SimResults)
    # Avoid the very last time--tends to be poorly-behaved
    # evalt = range(0.0, sim.sol.t[end]*0.99, length=100)
    virtual_thermocouple(locs, sim.sol.t, sim)
end
function virtual_thermocouple(locs, t, sim::SimResults)
    for loc in locs
        if length(loc) != 2
            @error "Each location should be a 2-tuple or similar."
        end
        if any(loc .< 0) || any(loc .> 1)
            @error "Location, in r and z, should be given as number between 0 and 1 inclusive." loc
        end
    end
    @unpack sol, dom, config = sim

    inds = map(locs) do (r,z)
        CI(round(Int, r*(dom.nr-1))+1, round(Int, z*(dom.nz-1))+1)
    end
    Tdat = map(t) do ti 
        params = calc_params_at_t(ti, config)
        u = sol(ti)
        identify_dry(dom, u.ϕ)
        Tf = calc_Tf_res(ti, sim)
        T = solve_T(u, Tf, dom, params)
        Tloc = T[inds]
        return Tloc
    end
    return stack(Tdat)'
end
virtual_thermocouple(sim::Dict) = virtual_thermocouple(sim["sim"])

"""
    $(SIGNATURES)

Compute surface area over time for the given simulaiton results.
    
Returns (ts, SA_t)
"""
function get_SA(res)
    ts = res["sol"].t
    get_SA(ts, res)
end
function get_SA(ts, res)
    @unpack sol, dom = res
    SA_t = map(ts) do ti
        ϕ = sol(ti).ϕ
        SA = compute_icesurf_δ(ϕ, dom)
    end
    return ts, SA_t
end

"""
    $(SIGNATURES)

Compute the average 𝑧 position of the sublimation front.
"""
function get_subf_z(ϕ, dom)
    δ = compute_discrete_δ(ϕ, dom)
    ave_z = sum(δ .* permutedims(dom.zgrid) .*dom.rgrid) / sum(δ .* dom.rgrid)
    return ave_z
end
"""
    $(SIGNATURES)

Compute the average 𝓇 position of the sublimation front.
"""
function get_subf_r(ϕ, dom)
    δ = compute_discrete_δ(ϕ, dom)
    ave_r = sum(δ .* permutedims(dom.zgrid) .*dom.rgrid) / sum(δ .* permutedims(dom.zgrid))
    return ave_r
end


