export summaryplot, resultsanim,  gen_sumplot, gen_anim
export plotframe, finalframe 
export calc_ϕTp_res, get_t_Tf, get_t_Tf_subflux, compare_lyopronto_res
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
        uTp = calc_ϕuTp_res(ti, simresults, simconfig)[2:4]
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
        uTp = calc_ϕuTp_res(ti, simresults, simconfig)[2:4]
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

    params = calc_params_at_t(t, simconfig)
    
    u = sol(t)
    ϕ = ϕ_T_from_u(u, dom)[1]
    # p_sub = calc_psub(Tf)
    T = solve_T(u, dom, params)
    if simconfig[:dudt_func] == dudt_heatonly!
        return u, ϕ, T, zeros(size(ϕ))
    end
    if isnothing(p0)
        p = solve_p(u, T, dom, params)
    else
        p = solve_p(u, T, dom, params, p0)
    end
    return u, ϕ, T, p
end

"""
    plotframe(t::Float64, simresults::Dict, simconfig::Dict; maxT=nothing, heatvar=:T)

Unpack simulation results and plot the state at time `t`.

`heatvar = :T` or `=:ϕ` or `=:p` decides whether temperature, level set function, or pressure is plotted as colors.
If given, `maxT` sets an upper limit for the associated colorbar.
"""
function plotframe(t::Float64, simresults::Dict, simconfig::Dict; maxT=nothing, heatvar=:T, p0=nothing)
    @unpack sol, dom = simresults

    u, ϕ, T, p = calc_uϕTp_res(t, simresults, simconfig; p0=p0)
    if heatvar == :ϕ 
        heatvar_vals = ϕ
        clab = "ϕ, m"
        cmap = :algae
        cont_c = :black
    elseif heatvar == :T 
        # T = solve_T(u, dom, params)
        heatvar_vals = T .- 273.15
        clab = " \nT, °C"
        cmap = :thermal
        if maximum(ϕ_T_from_u(u, dom)[2]) > maximum(T) # Tf > T
            cont_c = :black
        else
            cont_c = :white
        end
    elseif heatvar == :p
        # T = solve_T(u, dom, params)
        # p = solve_p(u, T, dom, params, p0)
        heatvar_vals = ustrip.(u"mTorr", p.*u"Pa")
        clab = "p, mTorr"
        cmap = :ice
        cont_c = :black
        
    else
        @warn "Invalid value of heatvar passed to `plotframe`. Should be :ϕ, :T, or :p." heatvar
    end

    clims = extrema(heatvar_vals)
    if clims[2] - clims[1] < 1e-4 && heatvar != :ϕ
        clims = clims .+ (0, 0.10)
    end

    tr = round(t/3600, digits=2)
    local pl = plot(aspect_ratio=:equal, xlim=(-dom.rmax,dom.rmax), ylim=(dom.zmin,dom.zmax))
    # plot_cylheat(heatvar_vals, dom; clims=clims)
    heatmap!(dom.rgrid, dom.zgrid, heatvar_vals', c=cmap, clims=clims)
    heatmap!(dom.rgrid .- dom.rmax, dom.zgrid, heatvar_vals[end:-1:begin, :]', c=cmap) # plot reflected
    plot_cylcont(ϕ, dom, c=cont_c)
    if heatvar == :ϕ
        Plots.contour!(dom.rgrid, dom.zgrid, ϕ', color=:black)
    end
    plot!(xlabel="t = $tr hr")
    plot!(colorbar_title=clab)
    plot!(x_ticks = ([-dom.rmax, 0, dom.rmax], ["-R", "0", "R"]), )
    plot!(y_ticks = ([0, dom.zmax], ["0", "L"]),  )
    # plot!(x_ticks=[-dom.rmax, 0, dom.rmax], xlabel="radius")
    return pl, heatvar_vals
end

function finalframe(simresults, simconfig; kwargs...)
    t = simresults["sol"].t[end]
    plotframe(t, simresults, simconfig; kwargs...)[1]
end

"""
    summaryplot(simresults::Dict, simconfig; layout=(3,2), heatvar=:T)

Return a 2x3 plot of simulation results from start to finish.

`simresults` should have a field `"sol"` , which is passed to `get_ϕ(sol, t, dom::Domain)` .  
`heatvar` determines what is plotted as a heatmap in the results (`:T` or `:ϕ`, currently.)
"""
function summaryplot(simresults::Dict, simconfig; layout=(3,2), heatvar=:T)
    @unpack sol, dom = simresults

    tf = sol.t[end]

    plots = []
    nplots = prod(layout)
    frames = range(0.0, tf*0.99, length=nplots)

    max_heat = 250.0
    min_heat = 250.0

    # T_nm1 = solve_T(sol(frames[end-1]), dom, cparams)
    # maxT = maximum(T_nm1)
    heatvals = fill(0.0, size(simresults["dom"]))

    for f in frames
        # p = plotframe(f, simresults, simconfig, maxT=maxT, heatvar=heatvar)
        pl, heatvals = plotframe(f, simresults, simconfig, heatvar=heatvar, p0= heatvals)
        ext_heat = extrema(heatvals)
        if f == 0
            min_heat, max_heat = ext_heat
        else
            min_heat = min(min_heat, ext_heat[1])
            max_heat = max(max_heat, ext_heat[2])
        end
        push!(plots, pl)
    end

    # for p_i in plots
    #     plot!(p_i, clims=(min_heat, max_heat))
    # end
    plsize = (1.25*(2*dom.rmax/dom.zmax)*200 * layout[2], 200*layout[1])
    bigplot = plot(plots..., size=plsize, layout=layout)
end

"""
    resultsanim(simresults, simconfig, casename; seconds_length=5, heatvar=:T)

Generate a .gif of the given simresults, with filename `casename_heatvar_evol.gif`.

Pass either `:p` or `:T` as `heatvar`. Passing `ϕ` will probably cause filename problems

"""
function resultsanim(simresults, simconfig, casename; seconds_length=5, heatvar=:T)
    @unpack sol, dom = simresults
    @unpack cparams= simconfig

    tf = sol.t[end]

    fps = 30
    frames = range(0, tf, length=seconds_length*fps)
    # if heatvar == :p
    heatvals = fill(0.0, size(dom))
    # end
    anim = @animate for ti ∈ frames
        pl, heatvals = plotframe(ti, simresults, simconfig, heatvar=heatvar, p0=heatvals)
        # heat_p_min = heat_ex[1] - 0.1*max(1e-3, heat_ex[2]-heat_ex[1])
        # heat_p_max = heat_ex[2] + 0.1*max(1e-3, heat_ex[2]-heat_ex[1])
        # plot!(p, clims=(heat_p_min, heat_p_max))
    end

    # fps = ceil(Int, nt/2)  # nt/x makes x-second long animation
    fname = "$(casename)_$(heatvar)_$(hash(simconfig))_evol.gif"
    gif(anim, plotsdir(fname), fps=fps)
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

