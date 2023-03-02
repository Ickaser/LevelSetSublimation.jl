export summaryplot, resultsanim, plotframe
export get_subf_z, get_subf_r

function plotframe(f::Int, simresults::Dict, simconfig::Dict; maxT=nothing)
    @unpack full_ϕ, full_T = simresults
    @unpack dom = simconfig
    local p = plot(aspect_ratio=:equal)
    plot_cylheat(full_T[f,:,:], dom; maxT=maxT)
    plot_cylcont(full_ϕ[f,:,:], dom, c=:white)
    plot!(title="timestep=$(f-1)")
    return p
end

"""
    summaryplot(simresults::Dict, simconfig; layout=(3,2))

Return a 2x3 plot of simulation results from start to finish.
"""
function summaryplot(simresults::Dict, simconfig; layout=(3,2))
    @unpack full_ϕ, full_T = simresults
    @unpack dom = simconfig

    nt = size(full_T, 1) 
    plots = []
    nplots = prod(layout)
    if nt >= nplots
        frames = round.(Int, range(1, nt-1, length=nplots))
    else
        frames = 1:nt
    end

    maxT = maximum(full_T[1:end-1,:,:])

    for f in frames
        # Default: plot heat
        p = plotframe(f, simresults, simconfig, maxT=maxT)
        
        # Debug: plot heat and level set
        # p1 = heat(full_T[f,:,:])
        # plot_contour(full_ϕ[f,:,:], c=:white)
        # # p2 = plot(aspect_ratio=:equal, xlim=(0,1), ylim=(0,1))
        # p2 = heat(full_ϕ[f,:,:])
        # plot_contour(full_ϕ[f,:,:], c=:white)
        # p = plot(p1, p2)

        # p = plot_frontvel(full_ϕ[f,:,:], full_T[f,:,:])

        # p2 = heat(full_ϕ[f,:,:])
        # markfront(full_ϕ[f,:,:])
        # plot_contour(full_ϕ[f,:,:], c=:black)
        # p = plot(p1, p2)

        push!(plots, p)
    end

    bigplot = plot(plots..., size=(500*layout[2], 200*layout[1]), layout=layout)
end

"""
    summaryplot(sol::ODESolution, simconfig; layout=(3,2))

Return a plot of simulation results (number from `layout`) from start to finish.
"""
function summaryplot(sol::ODESolution, simconfig; layout=(3,2))
    
    tf=sol.t[end]

    @unpack dom, T_params = simconfig

    plots = []
    nplots = prod(layout)
    frames = range(0.0, tf, length=nplots)

    maxT = maximum(solve_T(reshape(sol(frames[end-1]), dom.nr, dom.nz), dom, T_params))
    minT = T_params[:Tf]
    for ti in frames
        # Default: plot heat
        ϕ = reshape(sol(ti), dom.nr, dom.nz)
        T = solve_T(ϕ, dom, T_params)

        p = freshplot(dom)
        plot_cylheat(T, dom,maxT= maxT)
        # plot_cylheat(ϕ, dom)
        plot_cylcont(ϕ, dom, c=:white)
        plot!(title="t=$ti")

        push!(plots, p)
    end

    bigplot = plot(plots..., size=(500*layout[2], 200*layout[1]), layout=layout)
end
"""
    resultsanim(simresults, casename)

Generate a .gif of the given simresults, with filename casename_evol.gif.
TODO: generate names in the style of produce_or_load.
"""
function resultsanim(simresults, simconfig, casename)
    @unpack full_ϕ, full_T = simresults
    @unpack dom = simconfig
    nt = size(full_T, 1) 
    # freshplot()
    # plot!(size=(800,500))
    anim = @animate for i ∈ 1:nt
        # freshplot()
        plot(aspect_ratio=:equal, size=(800,500))
        plot!(title="timestep=$(i-1)")
        plot_cylheat(full_T[i,:,:], dom)
        plot_cylcont(full_ϕ[i,:,:], dom)
    end

    # fps = ceil(Int, nt/2)  # nt/x makes x-second long animation
    seconds_length = 3
    fname = "$(casename)_$(hash(simconfig))_evol.gif"
    gif(anim, plotsdir(fname), fps=ceil(Int, nt/seconds_length))
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

