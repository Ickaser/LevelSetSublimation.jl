using DrWatson, Test
@quickactivate :LevelSetSublimation


const LSS = LevelSetSublimation

# Define my own convenience function, since have lots of arrays to check
approxzero(x) = isapprox(x, 0, atol=100eps(typeof(x)))
# Directory function
testdir(args...) = projectdir("test", args...)


# Run test suite
println("Starting tests\n")
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


@testset "T solution: nothing here" begin
    # Compare to analytical:
    # 1D r direction, no ice 
    # 1D z direction, no ice 
    # 2D no ice 
    # 1D r direction, ghost θ > θ_thresh 
    # 1D r direction, ghost θ < θ_thresh 
    # 1D z direction, ghost θ > θ_thresh 
    # 1D z direction, ghost θ < θ_thresh 
end

@testset "T derivatives for velocity: no tests yet" begin
    # 1D r direction, ghost θ > θ_thresh 
    # 1D r direction, ghost θ < θ_thresh 
    # 1D z direction, ghost θ > θ_thresh 
    # 1D z direction, ghost θ < θ_thresh 
    # 
end

@testset "p solution: nothing here" begin
end

@testset "p derivatives for velocity: no tests yet" begin
end

@testset "velocity extrapolation" begin
end

@testset "comparison to LyoPronto, default parameters" begin
    # Read in precomputed results from LyoPronto
    # These are generated in scripts/default_lyopronto_config.jl, then saved
    lyopronto_res = load(testdir("lyopronto_default_results.jld2"))
    @unpack t_lp, T_lp, m_lp, dryfrac_lp = lyopronto_res
    # Read in a config object
    include(testdir("lyopronto_default_setup.jl"))
    res = sim_from_dict(config)
    tsol, Tsol, msol, fsol = compare_lyopronto_res(t_lp, res, config)

    @test sum(t_lp .== tsol) == 100
    @test sum(isapprox.(Tsol, T_lp, atol=0.5u"K")) == 100
    @test sum(isapprox.(msol, m_lp, atol=0.05u"kg/hr/m^2")) >= 98
    @test sum(isapprox.(fsol, dryfrac_lp, atol=0.02)) == 100
    @test res["sol"].t[end]*u"s" ≈ t_lp[end] rtol=0.02
end

include(testdir("tests_radialT.jl"))

ti = time() - ti
println("\nTesting took:")
println(round(ti/60, digits = 3), " minutes")
