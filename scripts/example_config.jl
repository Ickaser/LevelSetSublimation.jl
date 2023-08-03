# dom = Domain(51, 51, 1.0, 1.0)
cparams = make_default_params()
init_prof = :circ
Tf0 = 233.15u"K"
Q_gl_RF = RampedVariable(0.002u"W") # = volumetric * relevant vial volume
# t_samp = (0:0.1:1) .* u"hr"
# Tsh = 263.15u"K"
Tsh = RampedVariable([233.15u"K", 263.15u"K"], [1u"K/minute"], [10u"hr"])
Q_ic = RampedVariable(0.3u"W/cm^3")
p_ch = RampedVariable(100u"mTorr")

controls = Dict{Symbol, Any}()
@pack! controls = Q_gl_RF, Tsh, Q_ic, p_ch


vialsize = "10R"
fillvol = 2u"mL"

config = Dict{Symbol, Any}()
@pack! config = cparams, init_prof, Tf0, controls, vialsize, fillvol


# Set up stuff to make debugging easier
params, ncontrols = params_nondim_setup(cparams, controls)

r_vial = get_vial_radii(vialsize)[1]
z_fill = fillvol / π / r_vial^2

rmax = ustrip(u"m", r_vial)
zmax = ustrip(u"m", z_fill)

dom = Domain(51, 51, rmax, zmax)

ϕ0 = make_ϕ0(init_prof, dom)   
u0 = zeros(dom.ntot+2)
u0[1:dom.ntot] = reshape(ϕ0, :)
# u0[dom.ntot+1] = ustrip(u"K", Tf0)
# u0[dom.ntot+2] = ustrip(u"K", Tf0)
# reinitialize_ϕ_HCR!(ϕ0, dom)
# T0 = solve_T(u0, Tf0, dom, params)
# p0 = solve_p(u0, Tf0, T0, dom, params)
