
struct Domain{I,F}
    nr::I
    nz::I
    ntot::I
    rmin::F
    rmax::F
    zmin::F
    zmax::F
    rgrid::AbstractVector{F}
    zgrid::AbstractVector{F}
    dr::F
    dz::F
    dr1::F
    dz1::F
    dr2::F
    dz2::F
    bwr::I
    bwz::I
end

function Domain(nr::I, nz::I, rmax::F, zmax::F) where {I,F}
    rmin = 0.0
    zmin = 0.0
    return Domain(nr, nz, rmin, rmax, zmin, zmax)
end

function Domain(nr::I, nz::I, rmin::F, rmax::F, 
    zmin::F, zmax::F) where {I,F}

    ntot = nr*nz

    rgrid = range(rmin, rmax, length=nr)
    zgrid = range(zmin, zmax, length=nz)

    dr = rgrid[2] - rgrid[1]
    dz = zgrid[2] - zgrid[1]
    dr1 = 1/dr
    dz1 = 1/dz
    dr2 = 1/dr/dr
    dz2 = 1/dz/dz

    bwr = ceil(Int, 0.2*(rmax-rmin)/dr)
    bwz = ceil(Int, 0.2*(rmax-rmin)/dr)
    return Domain(
    nr::I, nz::I, ntot::I, rmin::F, rmax::F, zmin::F, zmax::F,
    rgrid::AbstractVector{F}, zgrid::AbstractVector{F},
    dr::F, dz::F, dr1::F, dz1::F, dr2::F, dz2::F,
    bwr::I, bwz::I
    )
end

function Domain(nr::I, nz::I, rmax::F, zmax::F,
    bwr::I, bwz::I) where {I,F}

    ntot = nr*nz

    rmin = 0.0
    zmin = 0.0

    rgrid = range(rmin, rmax, length=nr)
    zgrid = range(zmin, zmax, length=nz)

    dr = rgrid[2] - rgrid[1]
    dz = zgrid[2] - zgrid[1]
    dr1 = 1/dr
    dz1 = 1/dz
    dr2 = 1/dr/dr
    dz2 = 1/dz/dz

    return Domain(
    nr::I, nz::I, ntot::I, rmin::F, rmax::F, zmin::F, zmax::F,
    rgrid::AbstractVector{F}, zgrid::AbstractVector{F},
    dr::F, dz::F, dr1::F, dz1::F, dr2::F, dz2::F,
    bwr::I, bwz::I
    )
end

function Domain(nr::I, nz::I, rmin::F, rmax::F, zmin::F, zmax::F,
    bwr::I, bwz::I) where {I,F}
    ntot = nr*nz

    rgrid = range(rmin, rmax, length=nr)
    zgrid = range(zmin, zmax, length=nz)

    dr = rgrid[2] - rgrid[1]
    dz = zgrid[2] - zgrid[1]
    dr1 = 1/dr
    dz1 = 1/dz
    dr2 = 1/dr/dr
    dz2 = 1/dz/dz

    return Domain(
    nr::I, nz::I, ntot::I, rmin::F, rmax::F, zmin::F, zmax::F,
    rgrid::AbstractVector{F}, zgrid::AbstractVector{F},
    dr::F, dz::F, dr1::F, dz1::F, dr2::F, dz2::F,
    bwr::I, bwz::I
    )
end


