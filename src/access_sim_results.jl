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
            @info "md=$md" ti Tfi calc_psub(Tfi)
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
            Tf = sol.sol1(t).Tf
        else
            Tf = sol.Tf2(t)
        end
    elseif haskey(sim, :Tf)
        Tf = sim.Tf(t)
    else # Used a DAE or implicit solve
        Tf = sol(t).Tf
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
    T = solve_T(u, Tf, dom, params)
    p = solve_p(u, Tf, T, dom, params)
    return u, Tf, T, p
end

"""
    $(SIGNATURES)
"""
function virtual_thermocouple(sim::NamedTuple) 
    virtual_thermocouple([(0, 0)], sim)
end
function virtual_thermocouple(locs, sim::NamedTuple)
    # Avoid the very last time--tends to be poorly-behaved
    # evalt = range(0.0, sim.sol.t[end]*0.99, length=100)
    virtual_thermocouple(locs, sim.sol.t, sim)
end
function virtual_thermocouple(locs, t, sim::NamedTuple)
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


