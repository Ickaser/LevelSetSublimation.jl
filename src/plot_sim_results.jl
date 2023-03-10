export summaryplot, resultsanim, plotframe
export get_subf_z, get_subf_r, get_ϕ

# function plotframe(f::Int, simresults::Dict, simconfig::Dict; maxT=nothing)
#     @unpack full_ϕ, full_T = simresults
#     @unpack dom = simconfig
#     local p = plot(aspect_ratio=:equal)
#     plot_cylheat(full_T[f,:,:], dom; maxT=maxT)
#     plot_cylcont(full_ϕ[f,:,:], dom, c=:white)
#     plot!(title="timestep=$(f-1)")
#     return p
# end
function get_ϕ(sol::ODESolution, t, dom::Domain)
    reshape(sol(t), dom.nr, dom.nz)
end


function plotframe(t::Float64, simresults::Dict, simconfig::Dict; maxT=nothing, heatvar=:T)
    @unpack ϕsol = simresults
    @unpack dom, T_params = simconfig
    ϕ = reshape(ϕsol(t), dom.nr, dom.nz)

    tr = round(t, sigdigits=3)
    if heatvar == :T
        T = solve_T(ϕ, dom, T_params)
        local p = plot(aspect_ratio=:equal)
        plot_cylheat(T, dom; maxT=maxT)
        plot_cylcont(ϕ, dom, c=:white)
        plot!(title="timestep=$tr")
        return p
    elseif heatvar == :ϕ
        local p = plot(aspect_ratio=:equal)
        plot_cylheat(ϕ, dom;)
        plot_cylcont(ϕ, dom, c=:white)
        plot!(title="timestep=$tr")
        return p
    else
        @warn "Invalid variable to be plotted as heatmap" heatvar
        local p = plot(aspect_ratio=:equal)
        plot_cylcont(ϕ, dom, c=:white)
        plot!(title="timestep=$tr")
        return p
    end
end

"""
    summaryplot(simresults::Dict, simconfig; layout=(3,2), heatvar=:T)

Return a 2x3 plot of simulation results from start to finish.

`simresults` should have a field `"ϕsol"` , which is passed to `get_ϕ(ϕsol, t, dom::Domain)` .  
`heatvar` determines what is plotted as a heatmap in the results (`:T` or `:ϕ`, currently.)
"""
function summaryplot(simresults::Dict, simconfig; layout=(3,2), heatvar=:T)
    @unpack ϕsol = simresults
    @unpack dom, T_params= simconfig

    tf = ϕsol.t[end]

    plots = []
    nplots = prod(layout)
    frames = range(0.0, tf, length=nplots)

    T_nm1 = solve_T(get_ϕ(ϕsol, frames[end-1], dom), dom, T_params)
    maxT = maximum(T_nm1)

    for f in frames
        p = plotframe(f, simresults, simconfig, maxT=maxT, heatvar=heatvar)
        push!(plots, p)
    end

    bigplot = plot(plots..., size=(500*layout[2], 200*layout[1]), layout=layout)
end

"""
    resultsanim(simresults, simconfig, casename; seconds_length=3)

Generate a .gif of the given simresults, with filename casename_evol.gif.
TODO: generate names in the style of produce_or_load.
"""
function resultsanim(simresults, simconfig, casename; seconds_length=5)
    @unpack ϕsol = simresults
    @unpack dom, T_params= simconfig

    tf = ϕsol.t[end]

    # Maximum T for plotting
    ϕ = reshape(ϕsol(0.9*tf), dom.nr, dom.nz)
    T = solve_T(ϕ, dom, T_params)
    maxT = maximum(T)

    fps = 30
    frames = range(0, tf, length=seconds_length*fps)
    # freshplot()
    # plot!(size=(800,500))
    anim = @animate for ti ∈ frames
        # freshplot()
        plot(aspect_ratio=:equal, size=(800,500))
        ϕ = reshape(ϕsol(ti), dom.nr, dom.nz)
        T = solve_T(ϕ, dom, T_params)
        tr = round(ti, sigdigits=3)
        plot!(title="timestep=$tr")
        plot_cylheat(T, dom, maxT=maxT)
        plot_cylcont(ϕ, dom, c=:white)
    end

    # fps = ceil(Int, nt/2)  # nt/x makes x-second long animation
    fname = "$(casename)_$(hash(simconfig))_evol.gif"
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

