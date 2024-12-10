

# -------------------------
# Some base physical properties, which only change if using
# vials which are not borosilicate or if using solvent other than water
base_props = PhysicalProperties()

# ----------------------------
# Simulation geometry inputs: these names must be exact
vialsize = "6R"
fillvol = 5u"mL"
simgridsize = (41, 31)

dom = Domain(@dict vialsize fillvol simgridsize)

# -------------------------
# Case-specific inputs, documented in TimeConstantProperties
# These names can be changed freely, as long as they are passed
# to the TimeConstantProperties constructor in the right order

Rp0 = 1.4u"cm^2*hr*Torr/g"

# Varying pore structure constant in vertical direction
A1 = 16u"cm*hr*Torr/g"
A2 = 0.5u"1/cm"
Tguess = 260u"K"
# hd = (dom.zmax .- dom.zgrid)*u"m"
# dRpdh = @. A1/(1 + A2*hd)^2 # This is a constant for A2=0
# Varying l in vertical direction sometimes makes solve_p go wrong. Why? Stencil shift I think from theta large to small
# Very difficult to understand why, but maybe the small-theta stencil at Robin boundary next to Stefan is somehow ill-conditioned or in wrong units or something
# l = [sqrt(base_props.R*Tguess/base_props.Mw) / dRp for r in dom.rgrid, dRp in dRpdh] # Varying in height
l = sqrt(base_props.R*Tguess/base_props.Mw) / A1 # Constant in height

κ = 0u"m^2"
ϵ = 0.95
kd = LevelSetSublimation.k_sucrose*(1-ϵ)
Kvwf = 20.0u"W/m^2/K"
m_v = get_vial_mass(vialsize)
A_v = get_vial_radii(vialsize)[2]^2 *π
B_d = 0u"Ω/m^2"
B_f = 5e7u"Ω/m^2"
B_vw = 2e6u"Ω/m^2"

tcp = TimeConstantProperties(ϵ, l, κ, Rp0,
                            kd, Kvwf, m_v, A_v,
                            B_d, B_f, B_vw)

# ------------------------------
# Process parameters which may change over time.
# Documented in TimeVaryingProperties; names can be chosen, but order must be matched
# RampedVariable is a convenience type from LyoPronto.jl, documented there

Tsh = RampedVariable([233.15u"K", 283.15u"K"], 1u"K/minute")
pch = RampedVariable(100u"mTorr")
f_RF = RampedVariable(12u"GHz")
P_per_vial = RampedVariable(0.5u"W")

# Kv
Kc = 3.58e-4u"cal/s/K/cm^2" 
Kd = 11.6e-4u"cal/s/K/cm^2/Torr"
Kp = 0.46u"1/Torr"
Kshf = RpFormFit(Kc, Kd, Kp)

tvp = TimeVaryingProperties(f_RF, P_per_vial, Tsh, pch, Kshf)

# ---------------------------------
# Pack everything into a dictionary for simulating or saving
# The dictionary keys need to be exact, but 

paramsd = (base_props, tcp, tvp)

# This is one way of constructing the dictionary, if your variable 
# names match the required dictionary keys
config = @dict paramsd vialsize fillvol simgridsize
config[:dudt_func] = dudt_heatmass_dae!

# alt_simgridsize = (21, 21)
# # Or we can do the following, if variables have other names
# config = Dict{Symbol, Any}(
#     :paramsd => paramsd,
#     :vialsize => vialsize,
#     :fillvol => fillvol,
#     :simgridsize => alt_simgridsize,
# )

# -----------------------------
# Now, to run the simulation:

# just run it once and see progress in stdout...
res = sim_from_dict(config, verbose=true)

# # or, check first to see if the simulation has been run and load those results if it has already been run
# # if not already run, simulation runs quietly
# # See DrWatson documentation for more on this
# res, fname = produce_or_load(sim_and_postprocess, config; 
#     filename=hash, prefix=datadir("sims", "ex"))


# -------------------------------
# Generate some plots
summaryT(res, config)