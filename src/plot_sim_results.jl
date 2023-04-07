export summaryplot, resultsanim, plotframe, get_ϕ_Tf
export get_subf_z, get_subf_r, get_ϕ



"""
    plotframe(t::Float64, simresults::Dict, simconfig::Dict; maxT=nothing, heatvar=:T)

Unpack simulation results and plot the state at time `t`.

`heatvar = :T` or `=:ϕ` or `=:p` decides whether temperature, level set function, or pressure is plotted as colors.
If given, `maxT` sets an upper limit for the associated colorbar.
"""
function plotframe(t::Float64, simresults::Dict, simconfig::Dict; maxT=nothing, heatvar=:T)
    @unpack ϕsol = simresults
    @unpack dom, cparams = simconfig
    params = deepcopy(cparams)
    u = ϕsol(t)
    ϕ = ϕ_T_from_u(u, dom)[1]
    # p_sub = calc_psub(Tf)
    T = solve_T(u, dom, params)
    p = solve_p(u, T, dom, params)
    if heatvar == :T 
        heatvar_vals = T
    elseif heatvar == :ϕ 
        heatvar_vals = ϕ
    elseif heatvar == :p
        heatvar_vals = p # Either ϕ, 
    else
        @warn "Invalid value of heatvar passed to `plotframe`. Should be :ϕ, :T, or :p." heatvar
    end

    tr = round(t, sigdigits=3)
    local pl = plot(aspect_ratio=:equal)
    plot_cylheat(heatvar_vals, dom; maxT=maxT)
    plot_cylcont(ϕ, dom, c=:white)
    # Plots.contour!(dom.rgrid, dom.zgrid, ϕ')
    plot!(title="timestep=$tr")
    return pl
end

"""
    summaryplot(simresults::Dict, simconfig; layout=(3,2), heatvar=:T)

Return a 2x3 plot of simulation results from start to finish.

`simresults` should have a field `"ϕsol"` , which is passed to `get_ϕ(ϕsol, t, dom::Domain)` .  
`heatvar` determines what is plotted as a heatmap in the results (`:T` or `:ϕ`, currently.)
"""
function summaryplot(simresults::Dict, simconfig; layout=(3,2), heatvar=:T)
    @unpack ϕsol = simresults
    @unpack dom, cparams= simconfig

    tf = ϕsol.t[end]

    plots = []
    nplots = prod(layout)
    frames = range(0.0, tf, length=nplots)

    # T_nm1 = solve_T(ϕsol(frames[end-1]), dom, cparams)
    # maxT = maximum(T_nm1)

    for f in frames
        # p = plotframe(f, simresults, simconfig, maxT=maxT, heatvar=heatvar)
        p = plotframe(f, simresults, simconfig, heatvar=heatvar)
        push!(plots, p)
    end

    bigplot = plot(plots..., size=(500*layout[2], 200*layout[1]), layout=layout)
end

"""
    resultsanim(simresults, simconfig, casename; seconds_length=5, heatvar=:T)

Generate a .gif of the given simresults, with filename `casename_heatvar_evol.gif`.

Pass either `:p` or `:T` as `heatvar`. Passing `ϕ` will probably cause filename problems

TODO: generate names in the style of `produce_or_load`.
"""
function resultsanim(simresults, simconfig, casename; seconds_length=5, heatvar=:T)
    @unpack ϕsol = simresults
    @unpack dom, cparams= simconfig

    tf = ϕsol.t[end]

    # Maximum T for plotting
    # ϕ = reshape(ϕsol(0.9*tf), dom.nr, dom.nz)
    # T = solve_T(ϕ, dom, params)
    # maxT = maximum(T)

    fps = 30
    frames = range(0, tf, length=seconds_length*fps)
    # freshplot()
    # plot!(size=(800,500))
    anim = @animate for ti ∈ frames
        # freshplot()
        p = plotframe(ti, simresults, simconfig, heatvar=heatvar)
        # plot(aspect_ratio=:equal, size=(800,500))
        # ϕ = reshape(ϕsol(ti), dom.nr, dom.nz)
        # T = solve_T(ϕ, dom, params)
        # tr = round(ti, sigdigits=3)
        # plot!(title="timestep=$tr")
        # # plot_cylheat(T, dom, maxT=maxT)
        # plot_cylheat(T, dom)
        # plot_cylcont(ϕ, dom, c=:white)
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

