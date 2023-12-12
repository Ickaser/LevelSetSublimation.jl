using Optimization
using OptimizationPolyalgorithms
using SciMLSensitivity

# ------------- Read in experimental data
using CSV

dat1 = CSV.File(datadir("exp_raw", "3_8_23_Man05_10C_constPower10W_17vials_5ml_6R", "TemperatureLog_030223.csv"))
dat2 = CSV.File(datadir("exp_raw", "3_8_23_Man05_10C_constPower10W_17vials_5ml_6R", "Mannitol_5ml_Temp.csv"), comment="#", stripwhitespace=true)
start_ti = findfirst(dat2["SPLYO.VACUUM_SP.F_CV"] .== 100)
Pch = dat2["SPLYO.CHAMBER_CM.F_CV"][start_ti:end]
Tsh = dat2["SPLYO.SHELF_SP.F_CV"][start_ti:end]
Pch_pir = dat2["SPLYO.CHAMBER_PIRANI.F_CV"][start_ti:end]

tdat_raw = dat1["Elapsed [s]"]*u"s"
T1dat_raw = dat1["T1 [C]"]*u"°C"
T3dat_raw = dat1["T3 [C]"]*u"°C"
T4dat_raw = dat1["T4 [C]"]*u"°C"
sti = 92
# fti = length(tdat_raw)
fti = argmax(T4dat_raw)
tdat = tdat_raw[sti:fti]
T1dat = T1dat_raw[sti:fti]
T3dat = T3dat_raw[sti:fti]
T4dat = T4dat_raw[sti:fti]

tdat .-= tdat[1]
tdat = uconvert.(u"hr", tdat)


trim = floor.(Int64, range(1, length(tdat), length=101))

texp = tdat[trim]
Tf_exp = T4dat[trim]
Tw_exp = T3dat[trim]


# --------- Base config
cparams = make_default_params()
init_prof = :flat
vialsize = "6R"
fillvol = 5u"mL"
simgridsize = (51, 51)
Rp0 = 1.4u"cm^2*hr*Torr/g"
@pack! cparams = Rp0 # Same as above line
A1 = 16u"cm*hr*Torr/g"
Tguess = 260u"K"
l_from_A1 = sqrt(cparams[:R]*Tguess/cparams[:Mw]) / A1
cparams[:l] = l_from_A1
cparams[:κ] = 0u"m^2"
cparams[:ϵ] = 0.95
cparams[:m_cp_gl] = 7.9u"g" * 840u"J/kg/K"
cparams[:Kw] = 50.0u"W/m^2/K"
Tsh = RampedVariable([233.15u"K", 283.15u"K"], [1u"K/minute"], [10u"hr"])
p_ch = RampedVariable(100u"mTorr")
Q_ic = RampedVariable(0.12u"W/cm^3")
Q_gl_RF = RampedVariable(0.06u"W") # = volumetric * relevant vial volume
Tf0 = Tsh(0u"s") # Set to initial shelf temperature
Kc = 3.58e-4u"cal/s/K/cm^2" 
Kd = 11.6e-4u"cal/s/K/cm^2/Torr"
Kp = 0.46u"1/Torr"
cparams[:Kv] = Kc + Kd*p_ch(0)/(1+Kp*p_ch(0))
controls = Dict{Symbol, Any}()
@pack! controls = Q_gl_RF, Tsh, Q_ic, p_ch
base_config = Dict{Symbol, Any}()
@pack! base_config = cparams, init_prof, Tf0, controls, vialsize, fillvol
# base_config[:dudt_func] = dudt_heatmass_implicit!

# ------------- Set up loss function

guess_params = [.12, .06, 20]
# Loss function
function loss(calib_p)
    # calib_params = Dict{Symbol, Any}()
    if any(calib_p .< 0)
        return Inf, nothing
        # calib_p = abs.(calib_p)
    end
    # E_RFp_l = calib_p[1]
    # E_RFv_l = calib_p[2]*u"cm^3"
    # Kvp_l = calib_p[3]*u"W/K/m^2"
    # Kvp_l = calib_p[3]*u"cal/s/K/cm^2"

    Q_ic_l = RampedVariable(calib_p[1]*u"W/cm^3")
    Q_gl_RF_l = RampedVariable(calib_p[2]*u"W")
    Kvp_l = calib_p[3]*u"W/K/m^2"

    newconfig = deepcopy(base_config)
    newconfig[:controls][:Q_ic] = Q_ic_l
    newconfig[:controls][:Q_gl_RF] = Q_gl_RF_l
    newconfig[:cparams][:Kw] = Kvp_l

    @time res = sim_from_dict(newconfig)
    @unpack sol, dom = res

    max_t = min(tdat[end], sol.t[end]*u"s")
    ts_avail = findlast(tdat .< max_t)
    subtrim = floor.(Int64, range(1, ts_avail, length=101))
    tcomp = tdat[subtrim]

    @time Tf_sim = virtual_thermocouple(ustrip.(u"s", tcomp), res, newconfig)
    Tw_sim = sol(ustrip.(u"s", tcomp), idxs=dom.nr*(dom.nz+1)+1)

    # @info "sizes" length(tdat) length(sol.t) size(sol) sol.t tdat
    loss1 = sum(abs2, Tf_sim .- ustrip.(u"K", T4dat[subtrim])) /101 # Product T
    loss2 = sum(abs2, Tw_sim .- ustrip.(u"K", T3dat[subtrim])) /101 # Vial T
    loss3 = ustrip(u"hr", (max_t - tdat[end]))^2 # Total time
    if typeof(calib_p[1]) <: AbstractFloat
        @info "After solving:" loss1 loss2 loss3 calib_p
    end
    return loss1+loss2+ loss3, res
end

# @info "Running loss function once to make sure it works"
# @time loss(guess_params)[1]

# plot_callback = function (p, l, pred)
#     @info "callback" l
#     if isnothing(pred)
#         return false
#     end
#     # pl1 = plot(pred, idxs=[1,2])
#     pl1 = plot(pred.t .*u"hr", pred[1,:].*u"K", idxs=[1], c=:red)
#     plot!(tdat, T4dat)
#     pl2 = plot(pred.t .*u"hr", pred[2,:].*u"K", c=:red)
#     # pl2 = plot(tdat, T1dat)
#     plot!(tdat, T3dat)
#     display(plot(pl1, pl2, layout=(2,1)))
#     return false
# end

# adtype = Optimization.AutoForwardDiff()
# adtype = Optimization.AutoZygote()
adtype = Optimization.AutoFiniteDiff()
# adtype = Optimization.AutoEnzyme()
optf = Optimization.OptimizationFunction((x, p) -> loss(x), adtype)
optprob = Optimization.OptimizationProblem(optf, guess_params)

@info "About to start optimizing"
# @time result_params = Optimization.solve(optprob, PolyOpt(), callback=plot_callback, maxiters=100)
@time result_params = Optimization.solve(optprob, PolyOpt(), maxiters=100)
# @time result_params = Optimization.solve(optprob, PolyOpt())
    
@info "Optimized!" result_params result_params.u
result_res = loss(result_params.u)[2]
safesave(datadir("opts", "2023-11-02.jld2"), @strdict result_params result_res