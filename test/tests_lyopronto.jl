function lyoprontolike_setup()

    base_props = LevelSetSublimation.base_props
    vialsize = "6R"
    fillvol = 2u"mL"

    # R_p values to mass transfer
    Rp0 = 1.4u"cm^2*hr*Torr/g"
    A1 = 16u"cm*hr*Torr/g"
    Tguess = 250u"K"
    Mw = 18.002u"g/mol" # Water vapor
    l_bulk = upreferred(sqrt(u"R"*Tguess/Mw) / A1)
    ϵ = 0.95 
    κ = 0.0u"m^2" 
    kd = LyoPronto.k_sucrose * (1-ϵ)
    m_v = LyoPronto.get_vial_mass(vialsize)
    A_p, A_v = π.*LyoPronto.get_vial_radii(vialsize) .^2

    Kvwf = 0.0u"W/m^2/K" # Eliminate 2D effects
    Bd = 0.0u"Ω/m^2" # No microwave
    Bvw = 0.0u"Ω/m^2" # No microwave
    Bf = 0.0u"Ω/m^2" # No microwave
    tcprops = TimeConstantProperties(ϵ, l_bulk, κ, Rp0, kd, Kvwf, m_v, A_v, Bd, Bf, Bvw)

    Tsh = RampedVariable([238.15u"K", 293.15u"K"], 1u"K/minute")
    pch = RampedVariable(150u"mTorr")
    P_per_vial = RampedVariable(0u"W")
    f_RF = RampedVariable(8.0u"GHz") # Not used, so can set to whatever

    KC = 2.75e-4*u"cal/s/K/cm^2" 
    KP = 8.93e-4*u"cal/s/K/cm^2/Torr" 
    KD = 0.46*u"1/Torr"
    Kshf = RpFormFit(KC * (A_v/A_p), KP* (A_v/A_p) , KD)# Corrected for Av/Ap
    tvprops = TimeVaryingProperties(f_RF, P_per_vial, Tsh, pch, Kshf)
    paramsd = (base_props, tcprops, tvprops)
    Tf0 = -35.857u"°C"
    simgridsize = (41, 31)
    config = (; paramsd, fillvol, vialsize, simgridsize, Tf0, time_integ=Val(:exp_newton))

    lp_params = LyoPronto.ParamObjPikal(
        RpFormFit(Rp0, A1, 0.0u"cm^-1"), fillvol/A_p, 0.05u"g/mL", 1.0u"g/mL",
        RpFormFit(KC, KP, KD), A_v, A_p,
        pch, Tsh
    )
    lp_sim = solve(ODEProblem(lp_params), LyoPronto.odealg_chunk2)

    return config, lp_sim
end

@testset "comparison to LyoPronto, default parameters" begin
    # Prep for comparison
    config, lp_sim = lyoprontolike_setup()
    t_lp = lp_sim.t .* u"hr"
    dryfrac_lp = 1 .- lp_sim[1,:] ./ lp_sim[1,1] .|> NoUnits
    T_lp = lp_sim[2,:]u"K"
    m_lp = [-LyoPronto.calc_md_Q(ui, lp_sim.prob.p, ti)[1]|>u"kg/s" for (ui, ti) in zip(lp_sim.u, lp_sim.t)]

    @info "Starting lyopronto comparison simulation. Might take 20 minutes"
    @time res = sim_from_dict(config)
    tsol, Tsol, msol, fsol = compare_lyopronto_res(t_lp, res, config)

    @test sum(t_lp .== tsol) == 100
    @test sum(isapprox.(Tsol, T_lp, atol=0.5u"K")) == 100
    @test sum(isapprox.(msol, m_lp, atol=0.05u"kg/hr")) >= 98
    @test sum(isapprox.(fsol, dryfrac_lp, atol=0.02)) == 100
    @test res["sol"].t[end]*u"s" ≈ t_lp[end] rtol=0.02
end