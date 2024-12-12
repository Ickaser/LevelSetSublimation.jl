export heat, freshplot, markfront, markcells
export plot_cylheat, arrows
export LevelSet
export summaryplot, resultsanim,  gen_sumplot, gen_anim
export summaryT
export plotframe, finalframe 
export blankplothrC

"""
    heat(field, dom::Domain)

Return a heatmap of `field`, scaled to fit on `dom`.

This function generates a *new* plot.
"""
function heat(field, dom::Domain; kwargs...)
    return heatmap(dom.rgrid, dom.zgrid, field'; 
                    xlim=(dom.rmin,dom.rmax), ylim=(dom.zmin, dom.zmax),
                     aspect_ratio=:equal, kwargs...)
end

"""
    LevelSet(phi, dom::Domain)

A small struct, used in dispatch for plotting.

The plot recipe can be variously called as:

    plot(LevelSet(phi, dom), reflect=true)
    plot(LevelSet(phi, dom), reflect=false)
    plot!(LevelSet(phi, dom), reflect=true)
    plot!(LevelSet(phi, dom), reflect=false)
"""
struct LevelSet
    phi
    dom::Domain
end
@recipe function f(ls::LevelSet; reflect=true)
    cl = contour(ls.dom.rgrid,ls.dom.zgrid,ls.phi, 0.0)
    for seg in lines(cl)
        xs, ys = coordinates(seg)
        @series begin
            seriestype := :path
            color --> :black
            label --> ""
            (xs, ys)
        end
        if reflect
            @series begin
                primary := false # no legend entry
                seriestype := :path
                color --> :black
                (-xs, ys)
            end
        end
    end
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
    plotframe(t::Float64, sim; maxT=nothing, heatvar=:T)

Unpack simulation results and plot the state at time `t`.

`heatvar = :T` or `=:ϕ` or `=:p` decides whether temperature, level set function, or pressure is plotted as colors.
If given, `maxT` sets an upper limit for the associated colorbar.
"""
function plotframe(t, sim; heatvar=:T, clims=nothing)
    @unpack sol, dom = sim

    if heatvar == :ϕ 
        ϕ = reshape(sol(t, idxs=iϕ(dom)), size(dom))
        heatvar_vals = ϕ
        clab = "ϕ, m"
        # cmap = :algae
        cmap = :linear_green_5_95_c69_n256
        cont_c = :black
    elseif heatvar == :T 
        u, Tf, T, p = calc_uTfTp_res(t, sim)
        ϕ = reshape(u[iϕ(dom)], size(dom))
        T .-= 273.15
        Tvw = u[iTvw(dom)] - 273.15
        Tsh = ustrip(u"°C", sim.config[:paramsd][3].Tsh(t*u"s"))
        pl = plotframe_T(t, ϕ, T, Tvw, Tsh, dom; clims=clims)
        return pl, T, Tf
    elseif heatvar == :p
        u, Tf, T, p = calc_uTfTp_res(t, sim)
        ϕ = reshape(u[iϕ(dom)], size(dom))
        # T = solve_T(u, dom, params)
        # p = solve_p(u, T, dom, params, p0)
        heatvar_vals = ustrip.(u"mTorr", p.*u"Pa")
        clab = "p, mTorr"
        cmap = :ice
        cont_c = :black
        
    else
        @warn "Invalid value of heatvar passed to `plotframe`. Should be :ϕ, :T, or :p." heatvar
    end

    if isnothing(clims)
        clims = extrema(heatvar_vals)
        if clims[2] - clims[1] < 1e-4 && heatvar != :ϕ
            clims = clims .+ (0, 0.10)
        end
    end

    tr = round(t/3600, digits=2)
    local pl = plot(aspect_ratio=:equal, xlim=(-dom.rmax,dom.rmax), ylim=(dom.zmin,dom.zmax))
    # plot_cylheat(heatvar_vals, dom; clims=clims)
    heatmap!(dom.rgrid, dom.zgrid, heatvar_vals', c=cmap, clims=clims)
    heatmap!(dom.rgrid .- dom.rmax, dom.zgrid, heatvar_vals[end:-1:begin, :]', c=cmap, clims=clims) # plot reflected
    plot!(LevelSet(ϕ, dom), c=cont_c)
    if heatvar == :ϕ
        Plots.contour!(dom.rgrid, dom.zgrid, ϕ', color=:black)
    end
    plot!(xlabel="t = $tr hr")
    plot!(colorbar_title=clab)
    plot!(x_ticks = ([-dom.rmax, 0, dom.rmax], [L"-R", "0", L"R"]), )
    plot!(y_ticks = ([0, dom.zmax], ["0", L"h_{f0}"]),  )
    # plot!(x_ticks=[-dom.rmax, 0, dom.rmax], xlabel="radius")
    return pl, heatvar_vals, Tf
end

function plotframe_T(t, ϕ, T, Tvw, Tsh, dom; clims=nothing)
    cmap = :linear_bmy_10_95_c78_n256
    cg = cgrad(cmap)
    if isnothing(clims)
        clims = extrema(T)
        if clims[2] - clims[1] < 1e-4
            clims = clims .+ (0, 0.10)
        end
    end
    cpick = x->cg[(x-clims[1])/(clims[2]-clims[1])]

    hextend = 1.05
    vextend = 1.05
    lwall = Plots.Shape([(-dom.rmax*hextend, dom.zmax), 
             (-dom.rmax*hextend, 0),
             (-dom.rmax, 0),
             (-dom.rmax, dom.zmax)])
    rwall = Plots.Shape([( dom.rmax*hextend, dom.zmax), 
             ( dom.rmax*hextend, 0),
             ( dom.rmax, 0),
             ( dom.rmax, dom.zmax)])
    shelf = Plots.Shape([( dom.rmax*hextend, 0), 
             (-dom.rmax*hextend, 0),
             (-dom.rmax*hextend, -(vextend-1)*dom.zmax),
             ( dom.rmax*hextend, -(vextend-1)*dom.zmax)])

    # plots = map(zip(frames, ϕs, Ts, Tvws, Tshs)) do (t, ϕ, T, Tvw, Tsh)
    if (T[1] - clims[1])/(clims[2]-clims[1]) > 0.65 # Tf relatively high
        cont_c = :black
    else
        cont_c = :white
    end
    tr = round(t/3600, digits=2)

    # Set up plot
    pl = plot(aspect_ratio = :equal, showaxis=false,
            xlim=(-1, 1) .* (dom.rmax*hextend) , ylim=(dom.zmin-(vextend-1)*dom.zmax, dom.zmax))
    # plot temperatures
    heatmap!(dom.rgrid, dom.zgrid, T', c=cmap, clims=clims, )
    heatmap!(dom.rgrid .- dom.rmax, dom.zgrid, T[end:-1:begin, :]', c=cmap, clims=clims) # plot reflected
    # plot ice surface
    plot!(LevelSet(ϕ, dom), c=cont_c)
    # Plot shelf and wall Ts
    plot!(lwall, c=cpick(Tvw) , lw=0, linealpha=0, label="")
    plot!(rwall, c=cpick(Tvw) , lw=0, linealpha=0, label="")
    plot!(shelf, c=cpick(Tsh), lw=0, linealpha=0, label="")
    # Vial outline
    plot!([-1, -1, 1, 1] .* (dom.rmax), [dom.zmax, 0.0, 0.0, dom.zmax], c=:black, label=:none)
    plot!(xlabel="t = $tr hr")
#     plot!(colorbar_title=clab)
    plot!(x_ticks = ([-dom.rmax, 0, dom.rmax], [L"-R", "0", L"R"]), )
    plot!(y_ticks = ([0, dom.zmax], ["0", L"h_{f0}"]),  )
    # end
    return pl
end

function summaryT(sim; layout=(3,2), clims=nothing, tstart=0.01, tend=0.99)
    @unpack sol, dom = sim

    tf = sol.t[end]
    nplots = prod(layout)
    frames = range(tf*tstart, tf*tend, length=nplots)

    Ts = map(frames) do f
        u, Tf, T, p = calc_uTfTp_res(f, sim)
        T .- 273.15
    end
    ϕs = map(frames) do f
        reshape(sim.sol(f, idxs = iϕ(dom)), dom.nr, dom.nz)
    end
    Tvws = map(frames) do f
        sim.sol(f, idxs = dom.nr*(dom.nz+1)+1) - 273.15
    end
    Tshs = map(frames) do f
        Tsh = ustrip(u"°C", sim.config[:paramsd][3].Tsh(f*u"s"))
    end

    if clims === nothing
        cmax = max(maximum(maximum.(Ts)), maximum(Tvws), maximum(Tshs))
        cmin = min(minimum(minimum.(Ts)), minimum(Tvws), minimum(Tshs))
        clims = (cmin, cmax)
    end

    # cmap = :plasma
    cmap = :linear_bmy_10_95_c78_n256
    cg = cgrad(cmap)
    cpick = x->cg[(x-clims[1])/(clims[2]-clims[1])]

    lwall = Plots.Shape([(-dom.rmax-2.5dom.dr, dom.zmax), 
             (-dom.rmax-2.5dom.dr, 0),
             (-dom.rmax, 0),
             (-dom.rmax, dom.zmax)])
    rwall = Plots.Shape([( dom.rmax+2.5dom.dr, dom.zmax), 
             ( dom.rmax+2.5dom.dr, 0),
             ( dom.rmax, 0),
             ( dom.rmax, dom.zmax)])
    shelf = Plots.Shape([( dom.rmax+2.5dom.dr, 0), 
             (-dom.rmax-2.5dom.dr, 0),
             (-dom.rmax-2.5dom.dr, -2.5dom.dz),
             ( dom.rmax+2.5dom.dr, -2.5dom.dz)])

    plots = map(zip(frames, ϕs, Ts, Tvws, Tshs)) do (t, ϕ, T, Tvw, Tsh)
        if (T[1] - clims[1])/(clims[2]-clims[1]) > 0.65 # Tf relatively high
            cont_c = :black
        else
            cont_c = :white
        end
        tr = round(t/3600, digits=2)

        # Set up plot
        pl = plot(aspect_ratio = :equal, showaxis=false,
             xlim=(-1, 1) .* (dom.rmax+2.5dom.dr) , ylim=(dom.zmin-1.5dom.dz, dom.zmax))
        # plot temperatures
        heatmap!(dom.rgrid, dom.zgrid, T', c=cmap, clims=clims, cbar=nothing)
        heatmap!(dom.rgrid .- dom.rmax, dom.zgrid, T[end:-1:begin, :]', c=cmap, clims=clims, cbar=nothing) # plot reflected
        # plot ice surface
        plot!(LevelSet(ϕ, dom), c=cont_c)
        # Plot shelf and wall Ts
        plot!(lwall, c=cpick(Tvw) , lw=0, linealpha=0, label="")
        plot!(rwall, c=cpick(Tvw) , lw=0, linealpha=0, label="")
        plot!(shelf, c=cpick(Tsh), lw=0, linealpha=0, label="")
        # Vial outline
        plot!([-1, -1, 1, 1] .* (dom.rmax), [dom.zmax, 0.0, 0.0, dom.zmax], c=:black, label=:none)
        plot!(xlabel="t = $tr hr")
    #     plot!(colorbar_title=clab)
        plot!(x_ticks = ([-dom.rmax, 0, dom.rmax], [L"-R", "0", L"R"]), )
        plot!(y_ticks = ([0, dom.zmax], ["0", L"h_{f0}"]),  )
    end
    l = @layout [grid(layout...) a{0.03w}]

    c2 = plot([Inf], [Inf], zcolor=[Inf], clims=clims,
                 xlims=(1,1.1), label="", c=cmap, colorbar_title="Temperature [°C]", framestyle=:none)

    p_all = plot(plots..., c2, layout=l, link=:all)

    plsize = ((2*dom.rmax/dom.zmax)*200 * layout[2], 200*layout[1])
    plot!(size=plsize)
end

function finalframe(sim; kwargs...)
    t = sim.sol.t[end]
    plotframe(t, sim; kwargs...)[1]
end

"""
    summaryplot(sim; layout=(3,2), heatvar=:T)

Return a 2x3 plot of simulation results from start to finish.

`sim` should have a field `"sol"` , which is passed to `get_ϕ(sol, t, dom::Domain)` .  
`heatvar` determines what is plotted as a heatmap in the results (`:T` or `:ϕ`, currently.)
"""
function summaryplot(sim; layout=(3,2), tstart=0, tend=0.95, heatvar=:T)
    @unpack sol, dom = sim

    tf = sol.t[end]

    plots = []
    nplots = prod(layout)
    frames = range(tf*tstart, tf*tend, length=nplots)

    max_heat = 250.0
    min_heat = 250.0

    # T_nm1 = solve_T(sol(frames[end-1]), dom, cparams)
    # maxT = maximum(T_nm1)
    heatvals = fill(0.0, size(sim.dom))

    for f in frames
        # p = plotframe(f, sim, maxT=maxT, heatvar=heatvar)
        pl, heatvals = plotframe(f, sim, heatvar=heatvar)
        ext_heat = extrema(heatvals)
        if f == 0
            min_heat, max_heat = ext_heat
        else
            min_heat = min(min_heat, ext_heat[1])
            max_heat = max(max_heat, ext_heat[2])
        end
        push!(plots, pl)
    end

    plsize = (1.25*(2*dom.rmax/dom.zmax)*200 * layout[2], 200*layout[1])
    bigplot = plot(plots..., size=plsize, layout=layout)
end

"""
    resultsanim(sim, casename; seconds_length=5, heatvar=:T)

Generate a .gif of the given sim, with filename `casename_heatvar_evol.gif`.

Pass either `:p` or `:T` as `heatvar`. Passing `ϕ` will probably cause filename problems

"""
function resultsanim(sim, casename; seconds_length=5, heatvar=:T, clims=nothing)
    @unpack sol, dom = sim

    tf = sol.t[end]

    fps = 30
    frames = range(0, tf, length=seconds_length*fps)
    # if heatvar == :p
    heatvals = fill(0.0, size(dom))
    Tf_g = fill(245.0, dom.nr)
    # end
    anim = @animate for ti ∈ frames
        pl, heatvals, Tf_g = plotframe(ti, sim, heatvar=heatvar, Tf0=Tf_g, clims=clims)
        # heat_p_min = heat_ex[1] - 0.1*max(1e-3, heat_ex[2]-heat_ex[1])
        # heat_p_max = heat_ex[2] + 0.1*max(1e-3, heat_ex[2]-heat_ex[1])
        # plot!(p, clims=(heat_p_min, heat_p_max))
    end

    # fps = ceil(Int, nt/2)  # nt/x makes x-second long animation
    fname = "$(casename)_$(heatvar)_$(hash(sim.config))_evol.mp4"
    mp4(anim, plotsdir(fname), fps=fps)
end

blankplothrC() = plot(u"hr", u"°C", ylabel = "Temperature", xlabel = "Time", unitformat=:square)

@userplot TPlotModel
@recipe function f(tpm::TPlotModel)
    time, Ts = tpm.args
    step = size(time, 1) ÷ 10
    n = size(Ts, 2)
    pal = palette(:Oranges_4, rev=true) 
    
    for i in axes(Ts, 2)
        T = Ts[:,i]
        @series begin
            seriestype := :samplemarkers
            step := step
            offset := step÷n *(i-1) + 1
            markersize --> 7
            seriescolor --> pal[i]
            return time, T
        end
    end
end

@userplot VT_Plot
@recipe function f(vtp::VT_Plot)
    sim, locs = vtp.args
    time = sim.sol.t*u"s"
    Tfs = virtual_thermocouple(locs, sim)*u"K"
    Tvw = sim.sol.(sim.sol.t, idxs=iTvw(sim.dom))*u"K"
    @series begin
        TPlotModel((time, Tfs))
    end
    @series begin
        TPlotModVW((time, Tvw))
    end
end

@userplot PlaceThermocouples
@recipe function f(pt::PlaceThermocouples)
    dom, locs = pt.args
    locr = [loc[1] for loc in locs]
    locz = [loc[2] for loc in locs]
    @series begin
        seriestype := :scatter
        label --> ""
        markersize --> 7
        dom.rmax.*locr, dom.zmax.*locz
    end
end

@userplot LabelThermocouples
@recipe function f(lt::LabelThermocouples)
    dom, locs, labels = lt.args
    if isnothing(labels)
        @info "check" labels
        texts = [text(L"$T_text{f%$i}$", size=10) for i in 1:length(locs)]
    else
        texts = labels
    end
    locr = [loc[1] for loc in locs]
    locz = [loc[2] for loc in locs]
    @series begin
        @info "check" texts
        seriestype := :scatter
        markeralpha := 0
        primary := false
        text
        series_annotations := texts
        return dom.rmax.*locr, dom.zmax.*locz
    end
end

@userplot TPlotModVW
@recipe function f(tpmv::TPlotModVW)
    time, T = tpmv.args
    step = size(time, 1) ÷ 10
    color = palette(:Oranges_4)[end]
    
    @series begin
        seriestype := :samplemarkers
        step := step
        markershape --> :dtriangle
        markersize --> 7
        seriescolor --> color
        linestyle := :dash
        time, T
    end
end