using DrWatson, Test
@quickactivate :LevelSetSublimation

const LSS = LevelSetSublimation

dom1 = Domain(25, 24, 1.3, 2.0)
dom2 = Domain(45, 35, 1.2, 2.1)

# This is probably redundant after using COmponentArrays
function testing_u_funcs(dom::Domain)
    init_prof = :circ
    Tf0 = 233u"K"
    Tvw0 = 245u"K"
    u = LSS.make_u0_ndim(init_prof, Tf0, Tvw0, dom)
    ϕ1 = LSS.make_ϕ0(init_prof, dom)
    @testset "make_u0_ndim, first segment" for i in eachindex(ϕ1)
        @test u.ϕ[i] == ϕ1[i]
    end
    @testset "second segment" for i in eachindex(u.Tf)
        @test u.Tf[i] == ustrip(u"K", Tf0)
    end
    @test u.Tvw == ustrip(u"K", Tvw0)

    # Test with views
    ϕr, Tfr = u.ϕ, u.Tf
    ϕr .= π
    Tfr .= 250
    u.Tvw = 260
    @testset "make_u0_ndim, first segment" for i in eachindex(ϕ1)
        @test ϕr[i] ≈ π
    end
    @testset "second segment" for i in eachindex(u.Tf)
        @test u.Tf[i] == 250
    end
    @test length(u.Tvw) == 1
    @test u.Tvw == 260

end

@testset "size of ODE array, related functions" begin
    testing_u_funcs(dom1)
    testing_u_funcs(dom2)
end

@testset "radial T" begin
    @test true
end