using DrWatson, Test
@quickactivate :LevelSetSublimation

const LSS = LevelSetSublimation

dom1 = Domain(25, 24, 1.3, 2.0)
dom2 = Domain(45, 35, 1.2, 2.1)

function testing_u_funcs(dom::Domain)
    init_prof = :circ
    Tf0 = 233u"K"
    Tvw0 = 245u"K"
    u = LSS.make_u0_ndim(init_prof, Tf0, Tvw0, dom)
    @test length(u) == dom.ntot + dom.nr + 1
    ϕ1 = LSS.make_ϕ0(init_prof, dom)
    @testset "make_u0_ndim, first segment" for i in 1:dom.ntot
        @test u[i] == ϕ1[i]
    end
    @testset "second segment" for i in dom.ntot+1:dom.ntot+dom.nr
        @test u[i] == ustrip(u"K", Tf0)
    end
    @test u[end] == ustrip(u"K", Tvw0)

    # Test with views
    ϕr, Tfr, Tvwr = LSS.ϕ_T_from_u_view(u, dom)
    @testset "make_u0_ndim, first segment" for i in 1:dom.ntot
        @test ϕr[i] == ϕ1[i]
    end
    @testset "second segment" for i in 1:dom.nr
        @test Tfr[i] == ustrip(u"K", Tf0)
    end
    @test length(Tvwr) == 1
    @test Tvwr[1] == ustrip(u"K", Tvw0)

    ϕr .= π
    Tfr .= 250
    Tvwr .= 260
    @testset "views working, first seg" for i in 1:dom.ntot
        @test u[i] ≈ π
    end
    @testset "second segment" for i in dom.ntot+1:dom.ntot+dom.nr
        @test u[i] == 250
    end
    @test u[end] == 260

    # Reset, test without views
    u = LSS.make_u0_ndim(init_prof, Tf0, Tvw0, dom)
    ϕ2, Tf2, Tvw2 = LSS.ϕ_T_from_u(u, dom)
    @testset "ϕ_T_from_u, first segment" for i in 1:dom.ntot
        @test ϕ2[i] == ϕ1[i]
    end
    @testset "second segment" for i in 1:dom.nr
        @test Tf2[i] == ustrip(u"K", Tf0)
    end
    @test Tvw2 == ustrip(u"K", Tvw0)

    # Try writing, confirm that original array didn't change
    ϕ2 .= π
    Tf2 .= 250
    Tvw2 = 260
    @testset "non-view working, first seg" for i in 1:dom.ntot
        @test u[i] == ϕ1[i]
    end
    @testset "second segment" for i in dom.ntot+1:dom.ntot+dom.nr
        @test u[i] == ustrip(u"K", Tf0)
    end
    @test u[end] == ustrip(u"K", Tvw0)
    
end

@testset "size of ODE array, related functions" begin
    testing_u_funcs(dom1)
    testing_u_funcs(dom2)




end

@testset "radial T" begin
    @test true
end