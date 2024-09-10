
# A listing of physical properties for use in simulation, with references to sources

# -----------------------------------
# Sublimation pressure, to more significant figures than necessary but hey
# Left without units becuase this is used inside functions


# Triple point, without units: from IAPWS 1995, as rereleased in 2018
const psub_Tt = 273.16 # Triple point T, kelvin
const psub_pt = 611.654_771 # Triple point p, Pa


const psub_ei = [
 20.996_966_510_789_7,
  3.724_374_782_713_62,
-13.920_548_321_552_4,
 29.698_876_501_356_6,
-40.197_239_263_594_4,
 29.788_048_105_021_5,
 -9.130_509_635_477_21,
]
"""
    calc_psub_highorder(T)

Compute pressure (in Pascals) of sublimation at temperature `T` in Kelvin

From Feistel and Wagner, 2006
"""
function calc_psub_highorder(T)
    # if clamp(T, 20.0, 273.0) != T
    #     @warn "Invalid temperature for sublimation pressure correlation. Clamped to 20 and 273 K" T
    #     T = clamp(T, 20.0, 273.0)
    # end
    η = sum(psub_ei .* (T/psub_Tt) .^ (0:6))
    lnπ = 1.5log(T/psub_Tt) + (1 - psub_Tt/T) * η
    return exp(lnπ) * psub_pt
end

# -----------------------------
# Sublimation pressure, correlation used in LyoPRONTO: 
# simple Arrhenius form is almost certainly less accurate but probably not significantly so
# --------------------



# Heat of sublimation ----------
# Approximate: coresponds to 254K or 225K
# From Feistel and Wagner, 2006
const ΔH = 2838.0 * u"kJ/kg"

# Water vapor viscosity in dilute limit ---------------
# From Hellmann and Vogel, 2015
# 250K: 8.054 μPa*s
# 260K: 8.383 μPa*s
# 270K: 8.714 μPa*s
const μ = 8.1 * u"μPa*s"


# Molecular weight
const Mw = 18.015 * u"g/mol"
# Universal gas constant
const R = 8.3145* u"J/mol/K"

# Thermal conductivity of ice
# From Slack, 1980
# 200K : .032 W/cmK
# 250K : .024 W/cmK
# 273K : .0214 W/cmK
const kf = 2.4* u"W/m/K"

# Ice density
# IAPWS for ice: use triple point value for simplicity
const ρ_ice = 0.916e3 * u"kg/m^3"

const ρ_wat = 997.0u"kg/m^3"

# Ice heat capacity
# IAPWS for ice: use triple point value for simplicity
const Cp_ice = 0.209e4 * u"J/kg/K"

# Water vapor heat capacity near triple point
const Cp_vap = 1.86u"kJ/kg/K"

# From McCarthy and Fabre, 1989 book chapter
# Thermal conductivity of sucrose, Caster grade powder
const ρ_sucrose = 892u"kg/m^3"
const k_sucrose = 0.139u"W/m/K"

# From NIST, multiple sources
const Cp_sucrose = 422.50u"J/mol/K" / 342.3u"g/mol"

# From a solar paper and nanomaterials paper.
# Not clear what temperature (or even state) this is for
const k_mannitol1 = 0.60u"W/m/K"
const k_mannitol2 = 0.74u"W/m/K"

# Glass thermal conductivity
# Taken from papers making composites
# const k_gl = 2.0u"W/m/K"
# From Bansal and Doremus, rough estimate based an a sodium borosilicate glass
const k_gl = 1.0u"W/m/K"

# Glass heat capacity
# Rough estimate from Bansal and Doremus 1986, "Handbook of Glass Properties"
const cp_gl = 8.0u"J/kg/K"

# Typical value of porosity, heuristically
const ϵ_typical = 0.92

# Guessed values of porous structure parameters
const l_guess = 1e-5 * u"m" # 10 micron pore size
const κ_guess = 1e-10 * u"m^2" # 100 Darcys

# Typical value at around p_ch = 20 Pa
const Kshf_guess = 15.0u"W/m^2/K"
const p_ch_guess = 20.0u"Pa"

