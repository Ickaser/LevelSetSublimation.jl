export Domain

"""
Stores information on the grid and domain size for simulation.

**Constructors:**  
```julia
Domain(simconfig::Dict) 
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

# For convenience, extend `size` for Domains
Base.size(d::Domain) = (d.nr, d.nz)

function Domain(simconfig::Dict)
    @unpack vialsize, fillvol = simconfig
    simgridsize = get(simconfig, :simgridsize, (51,51))

    r_vial = get_vial_radii(vialsize)[1]
    Ap = π * r_vial^2
    z_fill = fillvol/Ap *ρ_wat/ρ_ice 
    rmax = ustrip(u"m", r_vial)
    zmax = ustrip(u"m", z_fill)

    dom = Domain(simgridsize..., rmax, zmax)
    return dom
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


struct CombinedSolution
    sol1
    sol2
    tsplit
    isplit
    t
    u
    Tf2
end

CombinedSolution(sol1, sol2) = CombinedSolution(sol1, sol2, sol1.t[end], length(sol1.t), vcat(sol1.t, sol2.t), vcat(sol1.u, sol2.u), nothing)
CombinedSolution(sol1, sol2, Tf2) = CombinedSolution(sol1, sol2, sol1.t[end], length(sol1.t), vcat(sol1.t, sol2.t), vcat(sol1.u, sol2.u), Tf2)
function (cs::CombinedSolution)(t) 
    if t <= cs.tsplit
        return cs.sol1(t)
    else
        return cs.sol2(t)
    end
end
function Base.getindex(cs::CombinedSolution, i)
    if i <= isplit
        return cs.sol1[i]
    elseif i <= isplit + length(cs.sol2.t)
        return cs.sol2[i]
    else
        throw(ArgumentError("Index must be less than combined solution length"))
    end
end
function Base.show(io::IO, cs::CombinedSolution)
    print(io, "CombinedSolution: ")
    # print(io, "sol1: $(cs.sol1), ")
    # print(io, "sol2: $(cs.sol2), ")
    print(io, "tsplit: $(cs.tsplit), ")
    print(io, "isplit: $(cs.isplit)")
end