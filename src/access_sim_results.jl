export calc_uϕTp_res, calc_uTfTp_res, get_t_Tf, get_t_Tf_subflux, compare_lyopronto_res
export get_subf_z, get_subf_r, get_ϕ, get_SA
export virtual_thermocouple

function get_t_Tf(simresults::Dict)
    @unpack sol, dom = simresults
    t = sol.t .* u"s"
    Tf = sol[dom.ntot+1,:] .* u"K"
    return t, Tf
end

function get_t_Tf_subflux(simresults::Dict, simconfig::Dict)
    @unpack sol, dom = simresults
    t = sol.t .* u"s"
    Tf_g = fill(245.0, dom.nr)
    Tf = sol[dom.ntot+1,:] .* u"K"
    md = map(sol.t) do ti
        params = calc_params_at_t(ti, simconfig)
        uTfTp = calc_uTfTp_res(ti, simresults, simconfig, Tf0=Tf_g)
        Tf_g = uTfTp[2]
        Tfi = uTfTp[1][iTf(dom)]
        md = compute_topmassflux(uTfTp..., dom, params) * u"kg/s"
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
    # Tf = sol(ts_ndim, idxs=dom.ntot+1).u .* u"K"
    Tf_g = fill(245.0, dom.nr)
    Tf = similar(ts_ndim)
    md = similar(ts_ndim).*u"kg/s"
    for (i, ti) in enumerate(ts_ndim)
        params = calc_params_at_t(ti, simconfig)
        uTfTp = calc_uTfTp_res(ti, simresults, simconfig, Tf0=Tf_g)
        Tf_g = uTfTp[2]
        # Tfi = uTfTp[1][iTf(dom)]
        Tf[i] = uTfTp[3][1,1]
        mdi = compute_topmassflux(uTfTp..., dom, params) * u"kg/s"
        if sign(mdi) == -1
            @info "md=$mdi" ti Tfi calc_psub(Tfi)
        end
        mdi = max(zero(mdi), mdi)
        md[i] = mdi
    end
    mfd = uconvert.(u"kg/hr", md) / (π*(dom.rmax*u"m")^2)
    totvol = π*dom.rmax^2 * dom.zmax
    dryfrac = map(ts_ndim) do ti
        ϕ = ϕ_T_from_u(sol(ti), dom)[1]
        1 - compute_icevol_H(ϕ, dom) / totvol
    end

    completed = dryfrac .== 1
    cyc = .~completed
    return ts[cyc], Tf[cyc].*u"K", mfd[cyc], dryfrac[cyc]
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
    @unpack paramsd = simconfig
    
    params = params_nondim_setup(paramsd)

    return (params[1], params[2], params[3](t))
end

function calc_uϕTp_res(t::Float64, simresults::Dict, simconfig::Dict; Tf0=nothing)
    @unpack sol, dom = simresults
    u, Tf, T, p = calc_uTfTp_res(t, simresults, simconfig; Tf0=Tf0)    
    ϕ = ϕ_T_from_u(u, dom)[1]
    return u, ϕ, T, p
end


function calc_uTfTp_res(t::Float64, simresults::Dict, simconfig::Dict; Tf0=nothing)
    @unpack sol, dom = simresults
    params = calc_params_at_t(t, simconfig)
    
    u = sol(t)
    # if haskey(simconfig, :dudt_func) && simconfig[:dudt_func] == dudt_heatmass_dae!
    #     Tf = sol(t, idxs=(dom.nr*dom.nz+1):(dom.nr*(dom.nz+1)))
    # elseif haskey(simconfig, :dudt_func) && simconfig[:dudt_func] == dudt_heatmass_implicit!
    #     Tf = sol(t, idxs=(dom.nr*dom.nz+1):(dom.nr*(dom.nz+1)))
    # else
    #     Tf = pseudosteady_Tf(u, dom, params, Tf0)
    # end
    Tf = sol(t, idxs=(dom.nr*dom.nz+1):(dom.nr*(dom.nz+1)))

    # p_sub = calc_psub(Tf)
    T = solve_T(u, Tf, dom, params)
    if haskey(simconfig, :dudt_func) && simconfig[:dudt_func] == dudt_heatonly!
        return u, Tf, T, zeros(size(dom))
    end
    # if isnothing(p0)
    #     p = solve_p(u, Tf, T, dom, params)
    # else
    #     p = solve_p(u, Tf, T, dom, params, p0)
    # end
    p = solve_p(u, Tf, T, dom, params)
    return u, Tf, T, p
end

function virtual_thermocouple(simresults::Dict, simconfig::Dict) 
    simt = simresults["sol"].t
    # Avoid the very last time--tends to be poorly-behaved
    evalt = range(0.0, simt[end]*0.99, length=100)
    virtual_thermocouple(0, 0, evalt, simresults, simconfig)
end
function virtual_thermocouple(t::TT, simresults::Dict, simconfig::Dict) where TT <: AbstractArray
    virtual_thermocouple(0, 0, t, simresults, simconfig)
end
function virtual_thermocouple(rpos, zpos, simresults::Dict, simconfig::Dict)
    simt = simresults["sol"].t
    evalt = range(0.0, simt[end]*0.99, length=100)
    virtual_thermocouple(rpos, zpos, evalt, simresults, simconfig)
end
function virtual_thermocouple(rpos, zpos, t::TT, simresults::Dict, simconfig::Dict) where TT <: AbstractArray
    if any(rpos .< 0) || any(rpos .> 1)
        @error "Radial position of thermocouple should be given as number between 0 and 1 inclusive." rpos
    elseif any(zpos .< 0) || any(zpos .> 1)
        @error "Axial (vertical) position of thermocouple should be given as number between 0 and 1 inclusive." zpos
    end
    if length(rpos) != length(zpos)
        @error "Number of radial positions should match number of vertical positions." rpos zpos
    end
    @unpack sol, dom = simresults
    # Tf = fill(245.0, dom.nr)
    ri = @. round(Int, rpos*(dom.nr-1)) + 1
    zi = @. round(Int, zpos*(dom.nz-1)) + 1
    inds = [CI(ir, iz) for (ir, iz) in zip(ri, zi)]
    Tdat = map(t) do ti 
        params = calc_params_at_t(ti, simconfig)
        u = sol(ti)
        # if haskey(simconfig, :dudt_func) && simconfig[:dudt_func] == dudt_heatmass_dae!
        #     Tf = sol(ti, idxs= (dom.nr*dom.nz+1):(dom.nr*(dom.nz+1)))
        # elseif haskey(simconfig, :dudt_func) && simconfig[:dudt_func] == dudt_heatmass_implicit!
        #     Tf = sol(ti, idxs= (dom.nr*dom.nz+1):(dom.nr*(dom.nz+1)))
        # else
        #     Tf = pseudosteady_Tf(u, dom, params, Tf)
        # end
        Tf = sol(ti, idxs= (dom.nr*dom.nz+1):(dom.nr*(dom.nz+1)))
        # @info "check" ti Tf[1]
        T = solve_T(u, Tf, dom, params)
        Tloc = T[inds]
    end
    return stack(Tdat)'
end

"""
    get_SA(res::Dict)
    get_SA(ts, res::Dict)

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
        u = sol(ti)
        ϕ = ϕ_T_from_u(u, dom)[1]
        SA = compute_icesurf_δ(ϕ, dom)
    end
    return ts, SA_t
end

"""
    get_subf_z(ϕ, dom)

Compute the average 𝑧 position of the sublimation front.
"""
function get_subf_z(ϕ, dom)
    δ = compute_discrete_δ(ϕ, dom)
    ave_z = sum(δ .* permutedims(dom.zgrid) .*dom.rgrid) / sum(δ .* dom.rgrid)
    return ave_z
end
"""
    get_subf_r(ϕ, dom)

Compute the average 𝓇 position of the sublimation front.
"""
function get_subf_r(ϕ, dom)
    δ = compute_discrete_δ(ϕ, dom)
    ave_r = sum(δ .* permutedims(dom.zgrid) .*dom.rgrid) / sum(δ .* permutedims(dom.zgrid))
    return ave_r
end


