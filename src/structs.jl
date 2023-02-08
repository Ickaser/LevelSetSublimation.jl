export Domain, LevelSet2DFields

"""
Stores information on the grid and domain size for simulation.

**Constructors:**  
```julia
Domain(nr, nz, rmax, zmax) 
Domain(nr, nz, rmax, zmax, bwr, bwz)
Domain(nr, nz, rmin, rmax, zmin, zmax)
Domain(nr, nz, rmin, rmax, zmin, zmax, bwr, bwz)
```

If not given, `rmin` and `zmin` default to 0.0, while `bwr` and `bwz` default to 20% of domain size.
Note that `nr`, `nz`, `bwr`, and `bwz` are assumed to be integers,
while `rmin`, `rmax`, `zmin`, and `zmax` are assumed to be floats

**Fields:** (given for 𝑟, likewise named for 𝑧)
- `nr` - number of grid points in 𝑟
- `ntot` - total number of grid points = `nr*nz`
- `rmin`, `rmax` - range of domain in 𝑟
- `rgrid` - vector of locations of grid points in 𝑟
- `dr` - grid spacing in 𝑟
- `dr1` - `1/dr`
- `dr2` - `1/dr^2`
- `bwr` - integer width of band around interface in  which level set is treated
"""
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

    dr = step(rgrid)
    dz = step(zgrid)
    dr1 = 1/dr
    dz1 = 1/dz
    dr2 = dr1*dr1
    dz2 = dz1*dz1

    bwr = ceil(Int, 0.2*nr)
    bwz = ceil(Int, 0.2*nz)
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


struct LevelSet2DFields{F}
    ϕ::Matrix{F}
    outside_B::F
    Γ::Vector{CartesianIndex{2}}
    B::Vector{CartesianIndex{2}}
    BnΓ::Vector{CartesianIndex{2}}
    ΩnB::Vector{CartesianIndex{2}}
    R::Vector{CartesianIndex{2}}
    C::Vector{CartesianIndex{2}}
end
