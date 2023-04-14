export summaryplot, resultsanim, plotframe,  gen_sumplot, gen_anim
export get_subf_z, get_subf_r, get_ϕ


function gen_sumplot(config, var=:T, casename="test")
    pol_kwargs = (filename=hash, prefix="simdat", verbose=false, tag=true)
    @time simres, simdatfile = produce_or_load(sim_from_dict, config,
            datadir("sims", casename); pol_kwargs...)
    summaryplot(simres, config, heatvar=var)
end
function gen_anim(config, var=:T, casename="test")
    pol_kwargs = (filename=hash, prefix="simdat", verbose=false, tag=true)
    @time simres, simdatfile = produce_or_load(sim_from_dict, config,
            datadir("sims", casename); pol_kwargs...)
    resultsanim(simres, config, casename, heatvar=var)
    return simres
end


"""
    plotframe(t::Float64, simresults::Dict, simconfig::Dict; maxT=nothing, heatvar=:T)

Unpack simulation results and plot the state at time `t`.

`heatvar = :T` or `=:ϕ` or `=:p` decides whether temperature, level set function, or pressure is plotted as colors.
If given, `maxT` sets an upper limit for the associated colorbar.
"""
function plotframe(t::Float64, simresults::Dict, simconfig::Dict; maxT=nothing, heatvar=:T, p0=nothing)
    @unpack sol, dom = simresults
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
            tip = findfirst(t_samp .> t)
            tip = clamp(tip, 2, length(t_samp))
            tim = tip - 1
            for ki in meas_keys
                params[ki] = (ncontrols[ki][tip] - ncontrols[ki][tim]) / (t_samp[tip] - t_samp[tim]) * (t - t_samp[tim]) + ncontrols[ki][tim]
            end
        end
    end
    
    u = sol(t)
    ϕ = ϕ_T_from_u(u, dom)[1]
    # p_sub = calc_psub(Tf)
    T = solve_T(u, dom, params)
    p = solve_p(u, T, dom, params, p0=p0)
    if heatvar == :T 
        heatvar_vals = T
        clab = "T, Kelvin"
    elseif heatvar == :ϕ 
        heatvar_vals = ϕ
        clab = "ϕ, m"
    elseif heatvar == :p
        heatvar_vals = ustrip.(u"mTorr", p.*u"Pa") # Either ϕ, 
        clab = "p, mTorr"
    else
        @warn "Invalid value of heatvar passed to `plotframe`. Should be :ϕ, :T, or :p." heatvar
    end

    tr = round(t/3600, digits=2)
    cont_c = (argmin(heatvar_vals) > argmax(heatvar_vals)) ? :black : :white
    local pl = plot(aspect_ratio=:equal)
    plot_cylheat(heatvar_vals, dom; maxT=maxT)
    plot_cylcont(ϕ, dom, c=cont_c)
    if heatvar == :ϕ
        Plots.contour!(dom.rgrid, dom.zgrid, ϕ', color=:black)
    end
    plot!(title="t = $tr hr")
    plot!(colorbar_title=clab)
    return pl, heatvar_vals
end

"""
    summaryplot(simresults::Dict, simconfig; layout=(3,2), heatvar=:T)

Return a 2x3 plot of simulation results from start to finish.

`simresults` should have a field `"sol"` , which is passed to `get_ϕ(sol, t, dom::Domain)` .  
`heatvar` determines what is plotted as a heatmap in the results (`:T` or `:ϕ`, currently.)
"""
function summaryplot(simresults::Dict, simconfig; layout=(3,2), heatvar=:T)
    @unpack sol = simresults

    tf = sol.t[end]

    plots = []
    nplots = prod(layout)
    frames = range(0.0, tf, length=nplots)

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

    bigplot = plot(plots..., size=(500*layout[2], 200*layout[1]), layout=layout)
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

