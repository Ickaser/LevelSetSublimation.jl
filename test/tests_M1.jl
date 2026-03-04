function M1_setup()

    paramsd = make_M1_properties()
    fillvol = 5u"mL"
    vialsize = "6R"
    simgridsize = (41, 35)
    config = (; paramsd, fillvol, vialsize, simgridsize, time_integ=Val(:dae_then_exp))
    return config
end

@testset "smoke test: M1 execution" begin
    # Prep for comparison
    config = M1_setup()
    @info "Starting M1 simulation. Might take 20 minutes"
    @time sim = sim_from_dict(config)["sim"]

    @test 5u"hr" < sim.sol.t[end]*u"s" < 20u"hr" 
end