
# This function generates a complete dictionary of parameters
# Many of them are physical constants that don't need tinkering,
# so I have not set their values here. But you can list those and
# change them if you want; a handful (Tsh, p_ch, Q_ic, Q_ck) 
# are controlled by the "controls" dict, for controlled variables.
cparams = make_default_params()

init_prof = :flat
vialsize = "6R"
fillvol = 5u"mL"
simgridsize = (51, 51)

# R_p values to dusty gas mass transfer coefficients
Rp0 = 1.4u"cm^2*hr*Torr/g"
# cparams[:Rp0] = Rp0
@pack! cparams = Rp0 # Same as above line
A1 = 16u"cm*hr*Torr/g"
Tguess = 260u"K"
l_from_A1 = sqrt(cparams[:R]*Tguess/cparams[:Mw]) / A1
# Empirical dusty gas parameter: like pore size
# Can be either a scalar, or an array matching the simulation grid size
cparams[:l] = l_from_A1
# l_varying = fill(l_from_A1, simgridsize)
# l_varying .*= range(0.5, 1.5, length=simgridsize[1]) # Make pores larger at outer radius
# # l_varying .*= permutedims(range(0.5, 1.5, length=simgridsize[1])) # Make pores larger at top of cake
# cparams[:l] = l_varying

# Empirical dusty gas parameter: like permeability
cparams[:κ] = 0u"m^2"

# Porosity: approximately 1 - (solids content)
cparams[:ϵ] = 0.95

# Vial mass times heat capacity
cparams[:m_cp_gl] = 7.9u"g" * 840u"J/kg/K"

# Vial wall to ice/cake heat transfer coefficenit
cparams[:Kw] = 50.0u"W/m^2/K"

# Shelf ramp: [setpoints], [ramprate between setpoints], [hold time after ramp]
# Same syntax holds for other variables; if constant, just pass the constant value 
Tsh = RampedVariable([233.15u"K", 283.15u"K"], [1u"K/minute"], [10u"hr"])
p_ch = RampedVariable(100u"mTorr")
Q_ic = RampedVariable(0.12u"W/cm^3")
Q_gl_RF = RampedVariable(0.06u"W") # = volumetric * relevant vial volume

# Initial temperature: used for wall initial and as initial guess for frozen temperature
# Tf0 = 233.15u"K" # Manual
Tf0 = Tsh(0u"s") # Set to initial shelf temperature

# Kv
Kc = 3.58e-4u"cal/s/K/cm^2" 
Kd = 11.6e-4u"cal/s/K/cm^2/Torr"
Kp = 0.46u"1/Torr"
# TODO: actually compute Kv as a function of chamber pressure internally
cparams[:Kv] = Kc + Kd*p_ch(0)/(1+Kp*p_ch(0))

controls = Dict{Symbol, Any}()
@pack! controls = Q_gl_RF, Tsh, Q_ic, p_ch

config = Dict{Symbol, Any}()
@pack! config = cparams, init_prof, Tf0, controls, vialsize, fillvol

#config[:dudt_func] = LSS.dudt_heatmass_dae!

