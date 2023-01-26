export heat, plot_contour, freshplot, markfront
export plot_cylheat, plot_cylcont, arrows

function heat(field, dom::Domain)
    return heatmap(dom.rgrid, dom.zgrid, field'; 
                    xlim=(dom.rmin,dom.rmax), ylim=(dom.zmin, dom.zmax),
                     aspect_ratio=:equal)
end
function plot_contour(phi, dom::Domain; lab="", c=:black)
    cl = levels(contour(dom.rgrid,dom.zgrid,phi, 0.0))
    for seg in lines(cl)
        xs, ys = coordinates(seg) # coordinates of this line segment
        plot!(xs, ys, color=c, label=lab) 
    end
end
freshplot(dom::Domain) = plot(aspect_ratio=:equal, 
                    xlim=(dom.rmin,dom.rmax), ylim=(dom.zmin,dom.zmax))

function plot_cylheat(T, dom::Domain)
    heatmap!(dom.rgrid, dom.zgrid, T', c=:thermal)
    heatmap!(dom.rgrid .- dom.rmax, dom.zgrid, T[end:-1:begin, :]', c=:thermal) # plot reflected
    plot!(xlim=(-dom.rmax,dom.rmax), ylim=(dom.zmin,dom.zmax))
end
function plot_cylcont(phi, dom::Domain; lab="", c=:black)
    cl = levels(contour(dom.rgrid,dom.zgrid,phi, 0.0))
    for seg in lines(cl)
        xs, ys = coordinates(seg) # coordinates of this line segment
        plot!(xs, ys, color=c, label=lab) 
        plot!(-xs, ys, color=c, label="") 
    end
end
function markfront(phi, dom::Domain; lab="", c=:white)
    frontcells = Γ_cells(phi, dom)
    rs = map(x-> dom.rgrid[Tuple(x)[1]] , frontcells)
    zs = map(x-> dom.zgrid[Tuple(x)[2]] , frontcells)
    scatter!(rs, zs, color=c, label=lab)
end

function arrows(Vf, dom::Domain)
    rmsh = reshape([r for r in dom.rgrid, z in dom.zgrid], :)
    zmsh = reshape([z for r in dom.rgrid, z in dom.zgrid], :)
    Vrmsh = reshape(Vf[:,:,1], :)
    Vzmsh = reshape(Vf[:,:,2], :)

    # freshplot()
    return quiver!(rmsh, zmsh, velocity=(Vrmsh, Vzmsh))
end