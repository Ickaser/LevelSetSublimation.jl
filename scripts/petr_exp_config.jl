using CSV

# dom = Domain(51, 51, 1.0, 1.0)
cparams = p = make_default_params()

cparams[:Cpf] = 2050u"J/kg/K"
Tf0 = 233.15u"K"

vialsize = "6R"
fillvol = 5u"mL"
simgridsize = (51, 51)

# R_p values to mass transfer
Rp0 = 1.4u"cm^2*hr*Torr/g"
@pack! cparams = Rp0
A1 = 16u"cm*hr*Torr/g"
Tguess = 260u"K"
l_bulk = sqrt(cparams[:R]*Tguess/cparams[:Mw]) / A1

# r_vial = get_vial_radii(vialsize)[1]
# L = fillvol / π/r_vial^2
# dz_approx = L / simgridsize[2]
# l_surf = l_bulk * (A1 / (R0/dz_approx ))
# l_surf = l_bulk * .05 # Experimental fit attempt

# l_bulk = upreferred(l_bulk)



init_prof = :flat

# # Shortly after ramp stops, no more t_samp: no need for callback to interfere
# t_samp = range(0, 2, step=1//60) .* u"hr"
# # T_sh = [233.15
# Tsh = fill(10.0u"°C", length(t_samp))
# # Tsh[1:101] .= range(-40u"°C", 10u"°C", length=101)
# Tsh[1:101] .-= range(50u"K", 0u"K", length=101) # Ramp from -40 to 10
# # @show Tsh
Tsh = RampedVariable([233.15u"K", 283.15u"K"], [1u"K/minute"], [10u"hr"])

# ----------- Fitting parameters! 
# With all four, get: 
# Q_ic = 0.5u"W/cm^3"
# Q_gl_RF = 0.16u"W"
# cparams[:Kgl] = 29.0u"W/m^2/K"
# cparams[:m_cp_gl] = 3.3u"g" * 840u"J/kg/K"

# With fixed vial mass including whole vial...
cparams[:m_cp_gl] = 7.9u"g" * 840u"J/kg/K"
Q_ic = RampedVariable(0.16u"W/cm^3")
Q_gl_RF = RampedVariable(0.13u"W") # = volumetric * relevant vial volume
cparams[:Kgl] = 100.0u"W/m^2/K"

# With fixed vial mass, varying l_surf
# cparams[:m_cp_gl] = 7.9u"g" * 840u"J/kg/K"
# Q_ic = 0.30u"W/cm^3"
# Q_gl_RF = 0.02u"W" # = volumetric * relevant vial volume
# cparams[:Kgl] = 600.0u"W/m^2/K"
# l_surf = .05*l_bulk

# Fixed mass, varying l_surf, less aggressive time matching
# cparams[:m_cp_gl] = 7.9u"g" * 840u"J/kg/K"
# Q_ic = 2.64u"W/cm^3"
# Q_gl_RF = 0.15u"W" # = volumetric * relevant vial volume
# cparams[:Kgl] = 225.0u"W/m^2/K"
# l_surf = 1.13*l_bulk



# -------------------------

# l_arr = fill(l_bulk, simgridsize)

cparams[:l] = l_bulk

p_ch = RampedVariable(100u"mTorr")

controls = Dict{Symbol, Any}()
# @pack! controls = t_samp, Q_gl_RF, Tsh, Q_ic, p_ch
@pack! controls = Q_gl_RF, Tsh, Q_ic, p_ch



config = Dict{Symbol, Any}()
@pack! config = cparams, init_prof, Tf0, controls, vialsize, fillvol


# # Set up stuff to make debugging easier
params, meas_keys, ncontrols = params_nondim_setup(cparams, controls)

dom = Domain(config)

u0 = make_u0_ndim(init_prof, Tf0, Tf0, dom)
# ϕ0 = make_ϕ0(init_prof, dom)   
# u0 = zeros(dom.ntot+2)
# u0[1:dom.ntot] = reshape(ϕ0, :)
# u0[dom.ntot+1] = ustrip(u"K", Tf0)
# u0[dom.ntot+2] = ustrip(u"K", Tf0)
# reinitialize_ϕ_HCR!(ϕ0, dom)
# T0 = solve_T(u0, dom, params)
# p0 = solve_p(u0, T0, dom, params)



exfname = datadir("exp_raw", "3_8_23_Man05_10C_constPower10W_17vials_5ml_6R", "TemperatureLog_030223.csv")
exdat = CSV.File(exfname)
t_ex = exdat["Elapsed [s]"]*u"s"
T1_ex = exdat["T1 [C]"]*u"°C"