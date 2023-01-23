
function heat(field)
    return heatmap(rgrid, zgrid, field'; xlim=(0.0, 1.0), ylim=(0.0, 1.0), aspect_ratio=:equal)
end
function plot_contour(phi; lab="", c=:black)
    cl = levels(contours(rgrid,zgrid,phi, [0]))[1]
    for seg in lines(cl)
        xs, ys = coordinates(seg) # coordinates of this line segment
        plot!(xs, ys, color=c, label=lab) 
    end
end
freshplot() = plot(aspect_ratio=:equal, xlim=(rmin,rmax), ylim=(zmin,zmax))

function plot_cylheat(T)
    heatmap!(rgrid, zgrid, T', c=:thermal)
    heatmap!(rgrid .- rmax, zgrid, T[end:-1:begin, :]', c=:thermal) # plot reflected
    plot!(xlim=(-rmax,rmax), ylim=(zmin,zmax))
end
function plot_cylcont(phi; lab="", c=:black)
    cl = levels(contours(rgrid,zgrid,phi, [0]))[1]
    for seg in lines(cl)
        xs, ys = coordinates(seg) # coordinates of this line segment
        plot!(xs, ys, color=c, label=lab) 
        plot!(-xs, ys, color=c, label="") 
    end
end
function markfront(phi; lab="", c=:white)
    frontcells = findall(identify_Γ(ϕ0))
    rs = map(x-> rgrid[Tuple(x)[1]] , frontcells)
    zs = map(x-> zgrid[Tuple(x)[2]] , frontcells)
    scatter!(rs, zs, color=c, label=lab)
end

function arrows(Vf)
    rmsh = reshape([r for r in rgrid, z in zgrid], :)
    zmsh = reshape([z for r in rgrid, z in zgrid], :)
    Vrmsh = reshape(Vf[:,:,1], :)
    Vzmsh = reshape(Vf[:,:,2], :)

    # freshplot()
    return quiver!(rmsh, zmsh, velocity=(Vrmsh, Vzmsh))
end