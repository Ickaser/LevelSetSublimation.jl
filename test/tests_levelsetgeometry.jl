# using DrWatson, Test
# @quickactivate :LevelSetSublimation

# const LSS = LevelSetSublimation

rmax = 2.2
zmax = 1.0
dom1 = Domain(25, 24,   rmax, zmax)
dom2 = Domain(47, 49,   rmax, zmax)

ϕ1a = make_ϕ0(:circ, dom1)
reinitialize_ϕ_HCR!(ϕ1a, dom1)
ϕ2a = make_ϕ0(:circ, dom2) 
reinitialize_ϕ_HCR!(ϕ2a, dom2)

a = dom1.rmax/1.1
c = dom1.zmax/1.1
if c < a # Oblate spheroid
    e = sqrt(1 - c^2/a^2)
    spheroid_surf = 2π*a^2*(1 + (1-e^2)/e*atanh(e))  / 2 # Half for half a spheroid
    # spheroid_surf = 2π*a^2 +  π*c^2/e*log((1+e)/(1-e))  / 2 # Half for half a spheroid
else # Prolate spheroid
    e = sqrt(1 - a^2/c^2)
    spheroid_surf = 2π*a^2*(1 + c/a/e*asin(e))  / 2
end

spheroid_vol = 4/3*π*a*a*c/2

@testset "surface areas" begin
    @test spheroid_surf ≈ LSS.compute_icesurf_δ(ϕ1a, dom1) atol=dom1.dz
    @test spheroid_surf ≈ LSS.compute_icesurf_δ(ϕ2a, dom2) atol=dom2.dz
    @test spheroid_surf ≈ LSS.compute_icesurf(ϕ1a, dom1)  atol=dom1.dz 
    @test spheroid_surf ≈ LSS.compute_icesurf(ϕ2a, dom2)  atol=dom2.dz 
end

@testset "volumes" begin
    @test spheroid_vol ≈ LSS.compute_icevol_H(ϕ1a, dom1) atol=dom1.dz
    @test spheroid_vol ≈ LSS.compute_icevol_H(ϕ2a, dom2) atol=dom2.dz
    @test spheroid_vol ≈ LSS.compute_icevol(ϕ1a, dom1)  atol=dom1.dz 
    @test spheroid_vol ≈ LSS.compute_icevol(ϕ2a, dom2)  atol=dom2.dz 
end

