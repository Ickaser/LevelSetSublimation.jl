using JLD2

function lyopronto_sim()

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
    A_v = π*LyoPronto.get_vial_radii(vialsize)[2]^2

    Kvwf = 0.0u"W/m^2/K" # Eliminate 2D effects
    Bd = 0.0u"Ω/m^2" # No microwave
    Bvw = 0.0u"Ω/m^2" # No microwave
    Bf = 0.0u"Ω/m^2" # No microwave
    tcprops = TimeConstantProperties(ϵ, l_bulk, κ, Rp0, kd, Kvwf, m_v, A_v, Bd, Bf, Bvw)

    Tsh = RampedVariable([238.15u"K", 293.15u"K"], 1u"K/minute")
    pch = RampedVariable(150u"mTorr")
    P_per_vial = RampedVariable(0u"W")
    f_RF = RampedVariable(8.0u"GHz") # Not used, so can set to whatever

    KC = 2.75e-4*u"cal/s/K/cm^2" * (3.8/3.14) # Corrected for Av/Ap
    KP = 8.93e-4*u"cal/s/K/cm^2/Torr" * (3.8/3.14) # Corrected for Av/Ap
    KD = 0.46*u"1/Torr"
    Kshf = RpFormFit(KC, KP, KD)

    tvprops = TimeVaryingProperties(f_RF, P_per_vial, Tsh, pch, Kshf)

    paramsd = (base_props, tcprops, tvprops)

    Tf0 = -35.857u"°C"
    simgridsize = (41, 31)

    config = (; paramsd, fillvol, vialsize, simgridsize, Tf0, time_integ=Val(:dae_then_exp))

    return config
end


# -----------------------------
# Read in LyoPronto data

# using CSV
# lpfname = datadir("lyopronto", "output_saved_230512_1714.csv")
# lpdat = CSV.File(lpfname)

# slice = floor.(Int, range(1, length(lpdat), length=100))
# t_lp = lpdat["Time [hr]"][slice]*u"hr"
# T_lp = lpdat["Sublimation Temperature [C]"][slice]*u"°C"
# Tsh_lp = lpdat["Shelf Temperature [C]"][slice]*u"°C"
# m_lp = lpdat["Sublimation Flux [kg/hr/m^2]"][slice]*u"kg/hr/m^2"
# dryfrac_lp = lpdat["Percent Dried"][slice] / 100

@testset "comparison to LyoPronto, default parameters" begin
    # Read in precomputed results from LyoPronto
    # These are generated in scripts/default_lyopronto_config.jl, then saved
    lyopronto_res = load(joinpath((@__DIR__),"lyopronto_default_results.jld2"))
    @unpack t_lp, T_lp, m_lp, dryfrac_lp = lyopronto_res
    # Read in a config object
    config = lyopronto_sim()
    @info "Starting lyopronto comparison simulation. Might take 20 minutes"
    @time res = sim_from_dict(config)
    tsol, Tsol, msol, fsol = compare_lyopronto_res(t_lp, res, config)

    @test sum(t_lp .== tsol) == 100
    @test sum(isapprox.(Tsol, T_lp, atol=0.5u"K")) == 100
    @test sum(isapprox.(msol, m_lp, atol=0.05u"kg/hr/m^2")) >= 98
    @test sum(isapprox.(fsol, dryfrac_lp, atol=0.02)) == 100
    @test res["sol"].t[end]*u"s" ≈ t_lp[end] rtol=0.02
end