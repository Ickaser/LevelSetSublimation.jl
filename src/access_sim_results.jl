export calc_uϕTp_res, calc_uTp_res, get_t_Tf, get_t_Tf_subflux, compare_lyopronto_res
export get_subf_z, get_subf_r, get_ϕ

function get_t_Tf(simresults::Dict)
    @unpack sol, dom = simresults
    t = sol.t .* u"s"
    Tf = sol[dom.ntot+1,:] .* u"K"
    return t, Tf
end

function get_t_Tf_subflux(simresults::Dict, simconfig::Dict)
    @unpack sol, dom = simresults
    t = sol.t .* u"s"
    Tf = sol[dom.ntot+1,:] .* u"K"
    md = map(sol.t) do ti
        params = calc_params_at_t(ti, simconfig)
        uTp = calc_uTp_res(ti, simresults, simconfig)
        Tfi = ϕ_T_from_u(uTp[1], dom)[2]
        md = compute_topmassflux(uTp..., dom, params) * u"kg/s"
        if sign(md) == -1
            @info "md=$md" ti Tfi calc_psub(Tfi)
        end
        md = max(zero(md), md)
        md
    end
    return t, Tf, md
end

function compare_lyopronto_res(simresults::Dict, simconfig::Dict)
    @unpack sol, dom = simresults
    t = uconvert.(u"hr", sol.t .* u"s")
    return compare_lyopronto_res(t, simresults, simconfig)
end

function compare_lyopronto_res(ts, simresults::Dict, simconfig::Dict)
    @unpack sol, dom = simresults
    # t = uconvert.(u"hr", sol.t .* u"s")
    ts_ndim = ustrip.(u"s", ts)
    Tf = sol(ts_ndim, idxs=dom.ntot+1).u .* u"K"
    md = map(ts_ndim) do ti
        params = calc_params_at_t(ti, simconfig)
        uTp = calc_uTp_res(ti, simresults, simconfig)
        Tfi = ϕ_T_from_u(uTp[1], dom)[2]
        mdi = compute_topmassflux(uTp..., dom, params) * u"kg/s"
        if sign(mdi) == -1
            @info "md=$mdi" ti Tfi calc_psub(Tfi)
        end
        mdi = max(zero(mdi), mdi)
        mdi
    end
    mfd = uconvert.(u"kg/hr", md) / (π*(dom.rmax*u"m")^2)
    totvol = π*dom.rmax^2 * dom.zmax
    dryfrac = map(ts_ndim) do ti
        ϕ = ϕ_T_from_u(sol(ti), dom)[1]
        1 - compute_icevol(ϕ, dom) / totvol
    end
    return ts, Tf, mfd, dryfrac
end

function gen_sumplot(config, var=:T, casename="test")
    pol_kwargs = (filename=hash, prefix="simdat", verbose=false, tag=true)
    @time simres, simdatfile = produce_or_load(sim_from_dict, config,
            datadir("sims", casename); pol_kwargs...)
    summaryplot(simres, config, heatvar=var, layout=(4,1))
end
function gen_anim(config, var=:T, casename="test")
    pol_kwargs = (filename=hash, prefix="simdat", verbose=false, tag=true)
    @time simres, simdatfile = produce_or_load(sim_from_dict, config,
            datadir("sims", casename); pol_kwargs...)
    resultsanim(simres, config, casename, heatvar=var)
    return simres
end

function calc_params_at_t(t::Float64, simconfig::Dict)
    @unpack cparams, controls = simconfig
    
    params, meas_keys, ncontrols = params_nondim_setup(cparams, controls)

    t_samp = get(ncontrols, :t_samp, 0.0)
    if meas_keys !== nothing
        # Interpolation here
        if t > t_samp[end] # Past end of sampling interval
            for ki in meas_keys
                params[ki] = ncontrols[ki][end]
            end
        else
            tim = findlast(t_samp .<= t)
            tim = clamp(tim, 1, length(t_samp)-1)
            tip = tim + 1
            tfrac = clamp((t - t_samp[tim]) / (t_samp[tip] - t_samp[tim]), 0, 1)
            for ki in meas_keys
                params[ki] = (ncontrols[ki][tip] - ncontrols[ki][tim])*tfrac  + ncontrols[ki][tim]
            end
        end
    end
    return params
end

function calc_uϕTp_res(t::Float64, simresults::Dict, simconfig::Dict; p0=nothing)
    @unpack sol, dom = simresults
    u, T, p = calc_uTp_res(t, simresults, simconfig; p0=p0)    
    ϕ = ϕ_T_from_u(u, dom)[1]
    return u, ϕ, T, p
end


function calc_uTp_res(t::Float64, simresults::Dict, simconfig::Dict; p0=nothing)
    @unpack sol, dom = simresults
    params = calc_params_at_t(t, simconfig)
    
    u = sol(t)
    # p_sub = calc_psub(Tf)
    T = solve_T(u, dom, params)
    if haskey(simconfig, :dudt_func) && simconfig[:dudt_func] == dudt_heatonly!
        return u, ϕ, T, zeros(size(ϕ))
    end
    if isnothing(p0)
        p = solve_p(u, T, dom, params)
    else
        p = solve_p(u, T, dom, params, p0)
    end
    return u, T, p
end


"""
    get_subf_z(ϕ, dom)

Compute the average 𝑧 position of the sublimation front.
"""
function get_subf_z(ϕ, dom)
    cl = contour(dom.rgrid, dom.zgrid, ϕ, 0.0)
    ls = lines(cl)
    if length(ls) == 0 # No sublimation front: average z is 0
        zbar = 0
    elseif length(ls) > 1
        @warn "Interface has more than one contiguous component"
        zbar = 0
        for line in ls
            rs, zs = coordinates(line)
            zbar += sum(zs) / length(zs)
        end
    else
        rs, zs = coordinates(ls[1])
        zbar = sum(zs) / length(zs)
    end
    zbar
end
"""
    get_subf_r(ϕ, dom)

Compute the average 𝓇 position of the sublimation front.
"""
function get_subf_r(ϕ, dom)
    cl = contour(dom.rgrid, dom.zgrid, ϕ, 0.0)
    ls = lines(cl)
    if length(ls) == 0 # No sublimation front: average z is 0
        rbar = 0
    elseif length(ls) > 1
        @warn "Interface has more than one contiguous component"
        rbar = 0
        for line in ls
            rs, zs = coordinates(line)
            rbar += sum(rs) / length(rs)
        end
    else
        line = ls[1]
        rs, zs = coordinates(line)
        rbar = sum(rs) / length(rs)
    end
    rbar
end

