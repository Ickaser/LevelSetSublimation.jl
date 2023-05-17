using CSV

# dom = Domain(51, 51, 1.0, 1.0)
cparams = p = make_default_params()

cparams[:Cpf] = 2050u"J/kg/K"
Tf0 = -35.857u"°C"

vialsize = "6R"
fillvol = 2u"mL"
simgridsize = (51, 51)

ϕ0type = :flat # Not always necessary


# ---- Variables which can be controlled during a run
# Set points

t_samp = range(0, 2, step=1//60) .* u"hr"
# Ramp from -35 to 20
Tsh = fill(20.0u"°C", length(t_samp))
Tsh[1:101] .= make_ramp(-35u"°C", 20u"°C", 1u"K/minute", t_samp[1:101])
p_ch = cparams[:p_ch] = 150u"mTorr"
Q_gl_RF = 0u"W"
Q_ic = 0u"W/cm^3"

# Kv value for heat transfer
KC = 2.75e-4*u"cal/s/K/cm^2" # Original
KP = 8.93e-4*u"cal/s/K/cm^2/Torr"
KD = 0.46*u"1/Torr"
cparams[:Kv] = (KC + KP*p_ch/(1+KD*p_ch)) 
cparams[:Kv] *= 3.8/3.14 # Correct for Av/Ap factor

# R_p values to mass transfer
Rp0 = 1.4u"cm^2*hr*Torr/g"
@pack! cparams = Rp0
A1 = 16u"cm*hr*Torr/g"
Tguess = 250u"K"
l_bulk = sqrt(cparams[:R]*Tguess/cparams[:Mw]) / A1
l_bulk = upreferred(l_bulk)
cparams[:l] = l_bulk
cparams[:ϵ] = 0.95 

cparams[:κ] *= 0 # Multiply by 0, to match dimensions
cparams[:Kgl] *= 0





# ----------- Fitting parameters! 
# With all four, get: 
# Q_ic = 0.5u"W/cm^3"
# Q_gl_RF = 0.16u"W"
# cparams[:Kgl] = 29.0u"W/m^2/K"
# cparams[:m_cp_gl] = 3.3u"g" * 840u"J/kg/K"

# With fixed vial mass including whole vial...
cparams[:m_cp_gl] = 7.9u"g" * 840u"J/kg/K"

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

controls = Dict{Symbol, Any}()
@pack! controls = t_samp, Q_gl_RF, Tsh, Q_ic, p_ch



config = Dict{Symbol, Any}()
@pack! config = cparams, ϕ0type, Tf0, controls, vialsize, fillvol


# --------------------------------------
# # Set up stuff to make debugging easier
# params, meas_keys, ncontrols = params_nondim_setup(cparams, controls)

# r_vial = get_vial_radii(vialsize)[1]
# z_fill = fillvol / π / r_vial^2

# rmax = ustrip(u"m", r_vial)
# zmax = ustrip(u"m", z_fill)

# dom = Domain(51, 51, rmax, zmax)

# ϕ0 = make_ϕ0(ϕ0type, dom)   
# u0 = zeros(dom.ntot+2)
# u0[1:dom.ntot] = reshape(ϕ0, :)
# u0[dom.ntot+1] = ustrip(u"K", Tf0)
# u0[dom.ntot+2] = ustrip(u"K", Tf0)
# reinitialize_ϕ_HCR!(ϕ0, dom)
# T0 = solve_T(u0, dom, params)
# p0 = solve_p(u0, T0, dom, params)


# -----------------------------
# Read in LyoPronto data

# exfname = datadir("exp_raw", "3_8_23_Man05_10C_constPower10W_17vials_5ml_6R", "TemperatureLog_030223.csv")
lpfname = datadir("lyopronto", "output_saved_230512_1714.csv")
lpdat = CSV.File(lpfname)
t_lp = lpdat["Time [hr]"]*u"hr"
T_lp = lpdat["Sublimation Temperature [C]"]*u"°C"
Tsh_lp = lpdat["Shelf Temperature [C]"]*u"°C"
m_lp = lpdat["Sublimation Flux [kg/hr/m^2]"]*u"kg/hr/m^2"
dryfrac_lp = lpdat["Percent Dried"] / 100