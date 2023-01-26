export summaryplot, resultsanim

"""
    summaryplot(simresults::Dict)

Return a 2x3 plot of simulation results from start to finish.
"""
function summaryplot(simresults::Dict, simconfig)
    @unpack full_ϕ, full_T = simresults
    @unpack dom = simconfig

    nt = size(full_T, 1) 
    plots = []
    if nt >= 6
        frames = round.(Int, range(1, nt, length=6))
    else
        frames = 1:nt
    end

    for f in frames
        # Default: plot heat
        local p = plot(aspect_ratio=:equal)
        plot_cylheat(full_T[f,:,:], dom)
        plot_cylcont(full_ϕ[f,:,:], dom, c=:white)
        plot!(title="timestep=$(f-1)")
        
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

    bigplot = plot(plots..., size=(500*2, 200*3), layout=(3,2))
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
    gif(anim, plotsdir("$(casename)_evol.gif"), fps=ceil(Int, nt/seconds_length))
end

function get_subf_z(ϕ, dom)
end