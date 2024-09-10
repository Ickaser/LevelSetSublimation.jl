cparams = make_default_params()

Tf0 = -35.857u"°C"
vialsize = "6R"
fillvol = 2u"mL"

simgridsize = (51, 51)
init_prof = :flat 


# ---- Variables which can be controlled during a run
# Set points


Tsh = RampedVariable([238.15u"K", 293.15u"K"], [1u"K/minute"], [10u"hr"])
p_ch = RampedVariable(150u"mTorr")
QRFvw = RampedVariable(0u"W")
QRFf = RampedVariable(0u"W/cm^3")

# Kshf value for heat transfer
KC = 2.75e-4*u"cal/s/K/cm^2" # Original
KP = 8.93e-4*u"cal/s/K/cm^2/Torr"
KD = 0.46*u"1/Torr"
cparams[:Kshf] = (KC + KP*p_ch(0u"s")/(1+KD*p_ch(0u"s"))) 
cparams[:Kshf] *= 3.8/3.14 # Correct for Av/Ap factor
cparams[:Kshf] *= .95 # Correct for ice thickness resistance


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
cparams[:Kvwf] *= 0

# Assemble parameters into objects -------------------------

controls = Dict{Symbol, Any}()
@pack! controls = QRFvw, Tsh, QRFf, p_ch

config = Dict{Symbol, Any}()
@pack! config = cparams, Tf0, init_prof, controls, vialsize, fillvol


# --------------------------------------
# # Set up stuff to make debugging easier
# params, meas_keys, ncontrols = params_nondim_setup(cparams, controls)

# r_vial = get_vial_radii(vialsize)[1]
# z_fill = fillvol / π / r_vial^2

# rmax = ustrip(u"m", r_vial)
# zmax = ustrip(u"m", z_fill)

# dom = Domain(51, 51, rmax, zmax)

# ϕ0 = make_ϕ0(init_prof, dom)   
# u0 = zeros(dom.ntot+2)
# u0[1:dom.ntot] = reshape(ϕ0, :)
# u0[dom.ntot+1] = ustrip(u"K", Tf0)
# u0[dom.ntot+2] = ustrip(u"K", Tf0)
# reinitialize_ϕ_HCR!(ϕ0, dom)
# T0 = solve_T(u0, dom, params)
# p0 = solve_p(u0, T0, dom, params)


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