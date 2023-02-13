export heat, plot_contour, freshplot, markfront
export plot_cylheat, plot_cylcont, arrows

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
    for seg in lines(cl)
        xs, ys = coordinates(seg) # coordinates of this line segment
        plot!(xs, ys, color=c, label=lab) 
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
function plot_cylheat(T, dom::Domain)
    heatmap!(dom.rgrid, dom.zgrid, T', c=:thermal)
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

Internally calls [Γ_cells](@ref).
"""
function markfront(phi, dom::Domain; lab="", c=:white)
    frontcells = Γ_cells(phi, dom)
    rs = map(x-> dom.rgrid[Tuple(x)[1]] , frontcells)
    zs = map(x-> dom.zgrid[Tuple(x)[2]] , frontcells)
    scatter!(rs, zs, color=c, label=lab)
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