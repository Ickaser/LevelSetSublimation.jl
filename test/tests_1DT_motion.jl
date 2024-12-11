# Pulled in largely as-is from 1D_verification_heat.ipynb

function case1_analyt_compare(simres, config)
    @unpack sol, dom = simres
    @unpack Kshf, ρf, ΔH, ϵ = config[:cparams]
    @unpack Tsh = config[:controls]
    @unpack Tf0 = config
    Q_sh = Kshf*(Tsh(0)-Tf0)
    vz = ustrip(u"m/s", -Q_sh/ ρf/ΔH / ϵ)

    ts = sol.t
    zs = [get_subf_z(reshape(sol(ti, idxs=iϕ(dom)), size(dom)), dom) for ti in ts]
    zs_analyt = dom.zmax*(1-1e-4) .+ vz.*ts

    return ts, zs, zs_analyt
end

function case2_analyt_compare(simres, config)
    @unpack sol, dom = simres
    @unpack Kshf, ρf, ΔH, ϵ = config[:cparams]
    QRFf = config[:controls][:QRFf]

    ts = sol.t
    rs = [get_subf_r(reshape(sol(ti, idxs=iϕ(dom)), size(dom)),dom) for ti in ts]
    rs_analyt = dom.rmax*(1-1e-4) * exp.(ustrip.(NoUnits, (-QRFf(0)/ρf/ΔH/ϵ/2 * ts*u"s") ))

    return ts, rs, rs_analyt
end

function case3_analyt_compare(simres, config)
    @unpack sol, dom = simres
    @unpack Kvwf, kd, ϵ, ρf, ΔH = config[:cparams]
    Tvw = config[:Tvw0]
    Tf = config[:Tf0]

    ts = sol.t *u"s"
    rs = [get_subf_r(reshape(sol(ustrip(u"s", ti), idxs=iϕ(dom)), size(dom)),dom) for ti in ts] *u"m"

    # analytical easier in terms of ξ, not t
    R = dom.rmax*u"m"
    R0 = R*(1 - 1e-4)
    rkk = R*Kvwf/kd
    A = rkk*(Tvw-Tf)
    B = rkk
    rrs = rs
    B = uconvert(NoUnits, B)

    term1 =@. 1/4*rrs^2 * (B*(2log(R/rrs)+1) + 2)
    term2 =   1/4*R0^2  * (B*(2log(R/R0 )+1) + 2)
    ts_analyt = @. ϵ*ρf*ΔH/A *(term2-term1)/kd

    return ts, rs, ts_analyt
end

@testset "1D motion, governed by T only" begin
# Setup
simgrid_coarse = (51, 51)
simgrid_fine = (101, 101)
pol_kwargs = (filename=hash, prefix="simdat", verbose=false, tag=true)

# Base parameter set
params_base = make_default_params()
params_base[:QRFf] *= 0
params_base[:Q_ck] *= 0
params_base[:Kvwf] *= 0
params_base[:Kshf] *= 0

# TODO: this is broken right now
# time_integ = LSS.dudt_heatonly!
init_prof = :flat
Tf0 = 233.15u"K"
vialsize = "10R"
fillvol = 2u"mL"
simgridsize = simgrid_coarse

QRFvw = 0.0u"W" # = volumetric * relevant vial volume
Tsh = 233.15u"K"
QRFf = 0.0u"W/cm^3"
pch = 100u"mTorr"
controls = Dict{Symbol, Any}()
@pack! controls = QRFvw, Tsh, QRFf, pch

cparams = deepcopy(params_base)
config_base = Dict{Symbol, Any}()
@pack! config_base = cparams, init_prof, Tf0, controls, vialsize, fillvol, simgridsize, time_integ

@testset "Case 1: Vertical Motion, Shelf Heating" begin
    
    casename = "sim_1D_z"

    config_1a = deepcopy(config_base)
    config_1a[:init_prof] = :flat
    config_1a[:controls][:Tsh] = RampedVariable(243.15u"K")
    config_1a[:cparams][:Kshf] = 10u"W/m^2/K"

    config_1b = deepcopy(config_1a)
    config_1b[:simgridsize] = simgrid_fine

    simres1a = sim_from_dict(config_1a, tf=1e6)
    simres1b = sim_from_dict(config_1b, tf=1e6)
    # simres1a, simdatfile1a = produce_or_load(sim_from_dict, config_1a,
    #             datadir("sims", "1D"); pol_kwargs...)
    # simres1b, simdatfile1b = produce_or_load(sim_from_dict, config_1b,
    #             datadir("sims", "1D"); pol_kwargs...)

    t1a, r1a, ra1a = case1_analyt_compare(simres1a, config_1a)
    t1b, r1b, ra1b = case1_analyt_compare(simres1b, config_1b)

    err1a = maximum(abs.(r1a - ra1a))
    err1b = maximum(abs.(r1b - ra1b))
    @test err1a < dom.dz
    @test err1b < dom.dz
end

@testset "Case 2: Radial motion, volumetric ice heating" begin
    casename = "sim_1D_r_ice"

    config_2a = deepcopy(config_base)
    config_2a[:init_prof] = :cyl
    config_2a[:controls][:QRFf] = RampedVariable(.1u"W/cm^3")

    config_2b = deepcopy(config_2a)
    config_2b[:simgridsize] = simgrid_fine

    # simres2a = sim_from_dict(config_2a)
    # simres2b = sim_from_dict(config_2b)
    simres2a, simdatfile2a = produce_or_load(sim_from_dict, config_2a,
                datadir("sims", "1D"); pol_kwargs...)
    simres2b, simdatfile2b = produce_or_load(sim_from_dict, config_2b,
                datadir("sims", "1D"); pol_kwargs...)

    t2a, r2a, ra2a = case2_analyt_compare(simres2a, config_2a)
    t2b, r2b, ra2b = case2_analyt_compare(simres2b, config_2b)

    err2a = maximum(abs.(r2a - ra2a))
    err2b = maximum(abs.(r2b - ra2b))
    @test err2a < dom.dz
    @test err2b < dom.dz
end

@testset "Case 3: Radial motion, constant Tvw" begin
    
    casename = "sim_1D_r_glass"

    config_3a = deepcopy(config_base)
    config_3a[:init_prof] = :cyl
    config_3a[:Tvw0] = 243.15u"K"
    config_3a[:cparams][:Kvwf] = 50u"W/m^2/K"

    config_3b = deepcopy(config_3a)
    config_3b[:simgridsize] = simgrid_fine

    # simres3a, simdatfile3a = produce_or_load(sim_from_dict, config_3a,
    #             datadir("sims", "1D"); pol_kwargs...)
    # simres3b, simdatfile3b = produce_or_load(sim_from_dict, config_3b,
    #             datadir("sims", "1D"); pol_kwargs...)
    simres3a = sim_from_dict(config_3a, tf=1e7)
    simres3b = sim_from_dict(config_3b, tf=1e7)
    t3a, r3a, ta3a = case3_analyt_compare(simres3a, config_3a)
    t3b, r3b, ta3b = case3_analyt_compare(simres3b, config_3b)

    maxerr3a = maximum((abs.(t3a - ta3a)./ta3a))
    maxerr3b = maximum((abs.(t3b - ta3b)./ta3b))
    aveerr3a = sum((abs.(t3a - ta3a)./ta3a)[begin+1:end])/length(t3a)
    aveerr3b = sum((abs.(t3b - ta3b)./ta3b)[begin+1:end])/length(t3b)
    
    @test_broken maxerr3a < .01
    @test_broken maxerr3b < .01
    @test 2aveerr3b < aveerr3a
end

end # Wrapping testset