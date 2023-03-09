export Domain, LevelSet2DFields

"""
Stores information on the grid and domain size for simulation.

**Constructors:**  
```julia
Domain(nr, nz, rmax, zmax) 
Domain(nr, nz, rmax, zmax, bwfrac)
Domain(nr, nz, rmin, rmax, zmin, zmax)
Domain(nr, nz, rmin, rmax, zmin, zmax, bwfrac)
```

If not given, `rmin` and `zmin` default to 0.0, while `bwfrac` defaults to 0.2.
Note that `nr`, `nz`, `bwr`, and `bwz` are assumed to be integers,
while `rmin`, `rmax`, `zmin`, and `zmax` are assumed to be floats

**All Fields:** (given for 𝑟, likewise named for 𝑧)
- `nr` - number of grid points in 𝑟
- `rmin`, `rmax` - range of domain in 𝑟
- `bwfrac` - fraction of domain included in band around interface
- `rgrid` - vector of locations of grid points in 𝑟
- `dr` - grid spacing in 𝑟
- `dr1` - `1/dr`
- `dr2` - `1/dr^2`
- `bwr` - `=ceil(Int, bwfrac*nr)` integer width of band around interface in  which level set is treated
- `ntot` - total number of grid points = `nr*nz`
"""
struct Domain{I,F}
    nr::I
    nz::I
    rmin::F
    rmax::F
    zmin::F
    zmax::F
    bwfrac::F
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
    ntot::I
end

# Pretty printing
function Base.show(io::IO, d::Domain) 
    if d.rmin == 0 && d.zmin == 0 && d.bwfrac == 0.2
        return print(io, "Domain($(d.nr), $(d.nz), $(d.rmax), $(d.zmax))")
    elseif d.rmin == 0 && d.zmin == 0
        return print(io, "Domain($(d.nr), $(d.nz), $(d.rmax), $(d.zmax), $(d.bwfrac))")
    elseif d.bwfrac == 0.2
        return print(io, "Domain($(d.nr), $(d.nz), $(d.rmin), $(d.rmax), $(d.zmin), $(d.zmax))")
    else
        return print(io, "Domain($(d.nr), $(d.nz), $(d.rmin), $(d.rmax), $(d.zmin), $(d.zmax), $(d.bwfrac))")
    end
end

function Domain(nr::I, nz::I, rmax::F, zmax::F) where {I,F}

    rmin = 0.0
    zmin = 0.0
    bwfrac = 0.2

    return Domain(nr, nz, rmin, rmax, zmin, zmax, bwfrac)
end

function Domain(nr::I, nz::I, rmin::F, rmax::F, 
    zmin::F, zmax::F) where {I,F}

    bwfrac = 0.2

    return Domain(nr::I, nz::I, rmin::F, rmax::F, 
                zmin::F, zmax::F, bwfrac::F)
end

function Domain(nr::I, nz::I, rmax::F, zmax::F,
    bwfrac::F) where {I,F}

    rmin = 0.0
    zmin = 0.0

    return Domain(nr::I, nz::I, rmin::F, rmax::F, 
                zmin::F, zmax::F, bwfrac::F)
end

function Domain(nr::I, nz::I, rmin::F, rmax::F, zmin::F, zmax::F,
    bwfrac::F) where {I,F}

    ntot = nr*nz

    rgrid = range(rmin, rmax, length=nr)
    zgrid = range(zmin, zmax, length=nz)

    dr = rgrid[2] - rgrid[1]
    dz = zgrid[2] - zgrid[1]
    dr1 = 1/dr
    dz1 = 1/dz
    dr2 = 1/dr/dr
    dz2 = 1/dz/dz

    bwr = ceil(Int, bwfrac*nr)
    bwz = ceil(Int, bwfrac*nz)

    return Domain(
    nr::I, nz::I, rmin::F, rmax::F, zmin::F, zmax::F, bwfrac::F,
    rgrid::AbstractVector{F}, zgrid::AbstractVector{F},
    dr::F, dz::F, dr1::F, dz1::F, dr2::F, dz2::F,
    bwr::I, bwz::I, ntot::I
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
