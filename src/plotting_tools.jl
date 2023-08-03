export heat, plot_contour, freshplot, markfront, markcells
export plot_cylheat, plot_cylcont, arrows
export summaryplot, resultsanim,  gen_sumplot, gen_anim
export plotframe, finalframe 

"""
    heat(field, dom::Domain)

Return a heatmap of `field`, scaled to fit on `dom`.

This function generates a *new* plot.
"""
function heat(field, dom::Domain)
    return heatmap(dom.rgrid, dom.zgrid, field'; 
                    xlim=(dom.rmin,dom.rmax), ylim=(dom.zmin, dom.zmax),
                     aspect_ratio=:equal)
end
"""
    plot_contour(phi, dom::Domain; lab="", c=:black)

In a mutating fashion, add 0-level contour of `phi` to current plot.
"""
function plot_contour(phi, dom::Domain; lab="", c=:black)
    # cl = levels(contours(dom.rgrid,dom.zgrid,phi, [0.0]))[1]
    cl = contour(dom.rgrid,dom.zgrid,phi, 0.0)
    p = nothing
    for seg in lines(cl)
        xs, ys = coordinates(seg) # coordinates of this line segment
        p = plot!(xs, ys, color=c, label=lab) 
    end
    p
end
"""
    freshplot(dom::Domain)  

Generate an empty plot, scaled to show only `dom`.
"""
freshplot(dom::Domain) = plot(aspect_ratio=:equal, xlim=(dom.rmin,dom.rmax), ylim=(dom.zmin,dom.zmax))

"""
    plot_cylheat(T, dom::Domain)

In a mutating fashion, add a "cylindrical" heatmap of `T` to the current plot.
("Cylindrical" meaning reflected across x=0 axis.)
"""
function plot_cylheat(T, dom::Domain; maxT=nothing, clims=nothing)
    if maxT === nothing maxT = maximum(T)
    end
    minT = minimum(T)
    if isnothing(clims)
        clims = (minT, maxT)
    end

    heatmap!(dom.rgrid, dom.zgrid, T', c=:thermal, clims=clims)
    heatmap!(dom.rgrid .- dom.rmax, dom.zgrid, T[end:-1:begin, :]', c=:thermal) # plot reflected
    plot!(xlim=(-dom.rmax,dom.rmax), ylim=(dom.zmin,dom.zmax))
end
"""
    plot_cylcont(ϕ, dom::Domain; lab="", c=:black)

In a mutating fashion, add a "cylindrical" 0-level contour of `ϕ` to the current plot.
("Cylindrical" meaning reflected across x=0 axis.)
"""
function plot_cylcont(phi, dom::Domain; lab="", c=:black)
    # cl = levels(contours(dom.rgrid,dom.zgrid,phi, [0.0]))[1]
    cl = contour(dom.rgrid,dom.zgrid,phi, 0.0)
    for seg in lines(cl)
        xs, ys = coordinates(seg) # coordinates of this line segment
        plot!(xs, ys, color=c, label=lab) 
        plot!(-xs, ys, color=c, label="") 
    end
end
"""
    markfront(phi, dom::Domain; lab="", c=:white)

Add a star marker to Γ (front) cells based on `ϕ` for current plot, in mutating fashion.

Internally calls `Γ_cells`.
"""
function markfront(phi, dom::Domain; lab="", c=:white)
    frontcells = Γ_cells(phi, dom)
    rs = map(x-> dom.rgrid[Tuple(x)[1]] , frontcells)
    zs = map(x-> dom.zgrid[Tuple(x)[2]] , frontcells)
    scatter!(rs, zs, color=c, label=lab)
end

"""
    markcells(cells, dom::Domain; lab="", c=:white)

Add a star marker for each CartesianIndex in  `cells` to the current plot.
"""
function markcells(cells, dom::Domain; lab = "", c=:white, kwargs...)
    rs = map(x-> dom.rgrid[Tuple(x)[1]], cells)
    zs = map(x-> dom.zgrid[Tuple(x)[2]], cells)
    scatter!(rs, zs, color=c, label=lab; kwargs...)
end


"""
    arrows(Vf, dom::Domain)

Add "quiver" to current plot, with velocity field `Vf` of shape (nx, ny, 2).
"""
function arrows(Vf, dom::Domain)
    rmsh = reshape([r for r in dom.rgrid, z in dom.zgrid], :)
    zmsh = reshape([z for r in dom.rgrid, z in dom.zgrid], :)
    Vrmsh = reshape(Vf[:,:,1], :)
    Vzmsh = reshape(Vf[:,:,2], :)

    # freshplot()
    return quiver!(rmsh, zmsh, velocity=(Vrmsh, Vzmsh))
end

"""
    function plot_frontvel(ϕ, T, dom::Domain)

Calculate, then plot the front velocity given `ϕ` and `T`.

Meant for debugging, mostly. Scales all velocity arrows to have length 0.5.
Generates a freshplot().
"""
function plot_frontvel(ϕ, T, dom::Domain, params)
    front_cells = findall(identify_Γ(ϕ, dom) .& (ϕ .> 0))
    xs = []
    ys = []
    vrs = []
    vzs = []
    for cell in front_cells
        push!(xs, dom.rgrid[Tuple(cell)[1]])
        push!(ys, dom.zgrid[Tuple(cell)[2]])
        vr, vz = compute_frontvel_withT(ϕ, T, Tuple(cell)..., dom, params)
        push!(vrs, vr)
        push!(vzs, vz)
        # push!(vrs, get_front_vr(T, ϕ, Tuple(cell)..., params) )
        # push!(vzs, get_front_vz(T, ϕ, Tuple(cell)..., params) )
    end
    maxv = max(maximum(abs.(vrs)), maximum(abs.(vzs)))
    println("Maximum front velocity: $maxv")
    vrs ./= maxv * 2
    vzs ./= maxv * 2

    freshplot(dom)
    quiver!(xs, ys, quiver=(vrs, vzs))
end

"""
    plotframe(t::Float64, simresults::Dict, simconfig::Dict; maxT=nothing, heatvar=:T)

Unpack simulation results and plot the state at time `t`.

`heatvar = :T` or `=:ϕ` or `=:p` decides whether temperature, level set function, or pressure is plotted as colors.
If given, `maxT` sets an upper limit for the associated colorbar.
"""
function plotframe(t::Float64, simresults::Dict, simconfig::Dict; maxT=nothing, heatvar=:T, p0=nothing)
    @unpack sol, dom = simresults

    if heatvar == :ϕ 
        ϕ = ϕ_T_from_u(sol(t), dom)[1]
        heatvar_vals = ϕ
        clab = "ϕ, m"
        cmap = :algae
        cont_c = :black
    elseif heatvar == :T 
        u, ϕ, T, p = calc_uϕTp_res(t, simresults, simconfig; p0=p0)
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
        u, ϕ, T, p = calc_uϕTp_res(t, simresults, simconfig; p0=p0)
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