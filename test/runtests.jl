using DrWatson, Test
@quickactivate :LevelSetSublimation

# Here you include files using `srcdir`
# include(srcdir("file.jl"))

const LSS = LevelSetSublimation

# Define my own convenience function, since have lots of arrays to check
approxzero(x) = isapprox(x, 0, atol=10eps(typeof(x)))


# Run test suite
println("Starting tests")
ti = time()

dom = Domain(25, 21, 2.0, 1.3)
ϕ1 = make_ϕ0(:flat, dom)
ϕ2 = make_ϕ0(:rad, dom)
ϕ3 = make_ϕ0(:box, dom)
ϕ4 = make_ϕ0(:circ, dom)
circ_area = π/4*(dom.rmax/1.1)*(dom.zmax/1.1)
dom_area = dom.rmax*dom.zmax

@testset "ϕ0 setup tests" begin
    @test maximum(ϕ1) ≈ dom.zmax*1e-4
    @test minimum(ϕ1) ≈ -dom.zmax*(1-1e-4)
    @test sum(ϕ1 .> 0) == dom.nr
    @test maximum(ϕ2) ≈ dom.rmax*1e-4
    @test minimum(ϕ2) ≈ -dom.rmax*(1-1e-4)
    @test sum(ϕ2 .> 0) == dom.nz
    @test maximum(ϕ3) ≈ max(dom.rmax, dom.zmax)*1e-4
    @test minimum(ϕ3) ≈ -min(dom.rmax, dom.zmax)*(1-1e-4)
    @test sum(ϕ3 .> 0) == dom.nz + dom.nr - 1
    @test maximum(ϕ4) ≈ (1.1*sqrt(2) - 1.0).*sqrt(dom.rmax*dom.zmax)
    @test minimum(ϕ4) ≈ -sqrt(dom.rmax*dom.zmax)
    @test sum(ϕ4 .< 0)/dom.ntot ≈ circ_area/dom_area rtol=1/min(dom.nr,dom.nz)
end

dϕdx_all_1 = LSS.dϕdx_all_WENO(ϕ1, dom)
dϕdx_all_2 = LSS.dϕdx_all_WENO(ϕ2, dom)
# Skip the box shape for now: doesn't have as nice derivatives
dϕdx_all_4 = LSS.dϕdx_all_WENO(ϕ4, dom)

rbroad = reshape(dom.rgrid, dom.nr, 1)
zbroad = reshape(dom.zgrid, 1, dom.nz)
dϕdr_4_anl = 1.1sqrt(dom.rmax*dom.zmax) * @. rbroad / dom.rmax^2 / [sqrt(r^2/dom.rmax^2 + z^2/dom.zmax^2) for r in dom.rgrid, z in dom.zgrid]
dϕdr_4_anl[1,1] = dϕdr_4_anl[2,1] # Avoid singularity
dϕdz_4_anl = 1.1sqrt(dom.rmax*dom.zmax) * @. zbroad / dom.zmax^2 / [sqrt(r^2/dom.rmax^2 + z^2/dom.zmax^2) for r in dom.rgrid, z in dom.zgrid]
dϕdz_4_anl[1,1] = dϕdz_4_anl[1,2] # Avoid singularity


@testset "WENO derivative tests" begin
    @test sum(approxzero.(dϕdx_all_1[1])) == dom.ntot  # East derivatives: 0
    @test sum(approxzero.(dϕdx_all_1[2])) == dom.ntot  # West derivatives: 0
    @test sum(dϕdx_all_1[3][:,begin:end-1] .≈ 1) == dom.ntot - dom.nr # North derivatives: 1 except on north boundary
    @test sum(approxzero.(dϕdx_all_1[3][:,end])) == dom.nr
    @test sum(dϕdx_all_1[4][:,begin+1:end] .≈ 1) == dom.ntot - dom.nr # South derivatives: 1 except on south boundary
    @test sum(approxzero.(dϕdx_all_1[4][:,begin])) == dom.nr

    @test sum(dϕdx_all_2[1][begin:end-1,:] .≈ 1) == dom.ntot - dom.nz # East derivatives: 1 except on east boundary
    @test sum(approxzero.(dϕdx_all_2[1][end,:])) == dom.nz
    @test sum(dϕdx_all_2[2][begin+1:end,:] .≈ 1) == dom.ntot - dom.nz # West derivatives: 1 except on west boundary
    @test sum(approxzero.(dϕdx_all_2[2][begin,:])) == dom.nz
    @test sum(approxzero.(dϕdx_all_2[3])) == dom.ntot  # North derivatives: 0
    @test sum(approxzero.(dϕdx_all_2[4])) == dom.ntot  # South derivatives: 0

    @test_broken sum(isapprox.(dϕdx_all_4[1][begin:end-1,:], dϕdr_4_anl[begin:end-1,:], atol=dom.dr)) == dom.ntot - dom.nz # North derivatives: break at north boundary
    @test_broken(isapprox.(dϕdx_all_4[2][begin+1:end,:], dϕdr_4_anl[begin+1:end,:], atol=dom.dr)) == dom.ntot - dom.nz # South derivatives: break at south boundary
    @test_broken sum(isapprox.(dϕdx_all_4[3][:,begin:end-1], dϕdz_4_anl[:,begin:end-1], atol=dom.dz)) == dom.ntot - dom.nr # East derivatives: break at east boundary
    @test_broken sum(isapprox.(dϕdx_all_4[4][:,begin+1:end], dϕdz_4_anl[:,begin+1:end], atol=dom.dz)) == dom.ntot - dom.nr # West derivatives: break ast west boundary
    @test sum(isapprox.(dϕdx_all_4[1][begin:end-1,:], dϕdr_4_anl[begin:end-1,:], atol=dom.dr^2)) >= dom.ntot - 2dom.nz # Allow breakage at both boundaries
    @test sum(isapprox.(dϕdx_all_4[2][begin+1:end,:], dϕdr_4_anl[begin+1:end,:], atol=dom.dr^2)) >= dom.ntot - 2dom.nz # Allow breakage at both boundaries
    @test sum(isapprox.(dϕdx_all_4[3][:,begin:end-1], dϕdz_4_anl[:,begin:end-1], atol=dom.dz^2)) >= dom.ntot - 2dom.nr # Allow breakage at both boundaries
    @test sum(isapprox.(dϕdx_all_4[4][:,begin+1:end], dϕdz_4_anl[:,begin+1:end], atol=dom.dz^2)) >= dom.ntot - 2dom.nr # Allow breakage at both boundaries
end

ti = time() - ti
println("\nTest took total time of:")
println(round(ti/60, digits = 3), " minutes")
