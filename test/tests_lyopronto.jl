
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