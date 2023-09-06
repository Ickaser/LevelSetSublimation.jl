dom = Domain(25, 21, 2.0, 1.3)
ϕ1 = make_ϕ0(:flat, dom)
ϕ2 = make_ϕ0(:rad, dom)
ϕ4 = make_ϕ0(:circ, dom)

dϕdx_all_1 = LSS.dϕdx_all_WENO(ϕ1, dom)
dϕdx_all_2 = LSS.dϕdx_all_WENO(ϕ2, dom)
# Skip the box shape for now: doesn't have as nice derivatives
dϕdx_all_4 = LSS.dϕdx_all_WENO(ϕ4, dom)
dϕdx_all_4_alt = LSS.dϕdx_all_WENO_loc(ϕ4, dom)

rbroad = reshape(dom.rgrid, dom.nr, 1)
zbroad = reshape(dom.zgrid, 1, dom.nz)
dϕdr_4_anl = 1.1sqrt(dom.rmax*dom.zmax) * @. rbroad / dom.rmax^2 / [sqrt(r^2/dom.rmax^2 + z^2/dom.zmax^2) for r in dom.rgrid, z in dom.zgrid]
dϕdr_4_anl[1,1] = dϕdr_4_anl[2,1] # Avoid singularity
dϕdz_4_anl = 1.1sqrt(dom.rmax*dom.zmax) * @. zbroad / dom.zmax^2 / [sqrt(r^2/dom.rmax^2 + z^2/dom.zmax^2) for r in dom.rgrid, z in dom.zgrid]
dϕdz_4_anl[1,1] = dϕdz_4_anl[1,2] # Avoid singularity


@testset "WENO derivative tests: tested on ellipse front" begin
    @test sum(dϕdx_all_4[1] .== dϕdx_all_4_alt[1]) == dom.ntot  # New and old implementations match
    @test sum(dϕdx_all_4[2] .== dϕdx_all_4_alt[2]) == dom.ntot  # New and old implementations match
    @test sum(dϕdx_all_4[3] .== dϕdx_all_4_alt[3]) == dom.ntot  # New and old implementations match
    @test sum(dϕdx_all_4[4] .== dϕdx_all_4_alt[4]) == dom.ntot  # New and old implementations match

    @test sum(approxzero.(dϕdx_all_1[1])) == dom.ntot  # East derivatives: 0
    @test sum(approxzero.(dϕdx_all_1[2])) == dom.ntot  # West derivatives: 0
    @test sum(dϕdx_all_1[3] .≈ 1) == dom.ntot-dom.nr # North derivatives: 1 
    @test sum(dϕdx_all_1[4] .≈ 1) == dom.ntot-dom.nr # South derivatives: 1 

    @test sum(dϕdx_all_2[1] .≈ 1) == dom.ntot-dom.nz  # East derivatives: 1
    @test sum(dϕdx_all_2[2] .≈ 1) == dom.ntot-dom.nz  # West derivatives: 1
    @test sum(approxzero.(dϕdx_all_2[3])) == dom.ntot  # North derivatives: 0
    @test sum(approxzero.(dϕdx_all_2[4])) == dom.ntot  # South derivatives: 0

    # This should maybe not be tested on this domain with a discontinuous slope?
    @test sum(isapprox.(dϕdx_all_4[1], dϕdr_4_anl, atol=dom.dr)) >= dom.ntot - dom.nz # North derivatives: allow for breakage at north boundary
    @test sum(isapprox.(dϕdx_all_4[2], dϕdr_4_anl, atol=dom.dr)) >= dom.ntot - dom.nz - 5 # South derivatives: break at south edge, center discontinuity
    @test sum(isapprox.(dϕdx_all_4[3], dϕdz_4_anl, atol=dom.dz)) >= dom.ntot - dom.nr # East derivatives : allow for breakage at east boundary
    @test sum(isapprox.(dϕdx_all_4[4], dϕdz_4_anl, atol=dom.dz)) >= dom.ntot - dom.nr - 5 # West derivatives : allow for breakage at west boundary, center discontinuity

    @test_broken sum(isapprox.(dϕdx_all_4[1], dϕdx_all_4[2], atol=dom.dr)) == dom.ntot # East and west are within grid size
    @test_broken sum(isapprox.(dϕdx_all_4[3], dϕdx_all_4[4], atol=dom.dz)) == dom.ntot # North and south are within grid size
end

ϕ_pre = make_ϕ0(:circ, dom)
ϕ_post = reinitialize_ϕ_HCR(ϕ_pre, dom, maxsteps=500, tol = 0.01)
const sdf_err_L∞ = LSS.sdf_err_L∞

vol_pre = LSS.compute_icevol(ϕ_pre, dom)
vol_post = LSS.compute_icevol(ϕ_post, dom)
surf_pre = LSS.compute_icesurf(ϕ_pre, dom)
surf_post = LSS.compute_icesurf(ϕ_post, dom)

@testset "Reinitialization function: tested only on ellipse" begin
    @test sdf_err_L∞(ϕ_pre, dom) > 0.01
    @test sdf_err_L∞(ϕ_post, dom) < 0.01
    @test vol_pre ≈ vol_post rtol = 1/max(dom.nr, dom.nz)^2 # Second order mass conservation
    @test surf_pre ≈ surf_post rtol = 1/max(dom.nr, dom.nz)^2 # Second order area conservation
end
