
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