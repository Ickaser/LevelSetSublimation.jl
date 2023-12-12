export get_vial_radii
# Get dimensions from here: https://www.schott.com/en-us/products/vials/-/media/project/onex/products/v/vials/application-variants/schott-brochure-schott-vials-english-20092017.pdf

# Added to "vial_sizes.csv", read in here
const VIAL_DIMS = CSV.File(srcdir("vial_sizes.csv"))

"""
    get_vial_radii(vialsize::String)

Return inner and outer radius for passed ISO vial size.

Uses a table provided by Schott, stored internally in a CSV.
"""
function get_vial_radii(vialsize::String)
    alldims = [row for row in VIAL_DIMS if row.Size == vialsize]
    if length(alldims) != 1
        @error "bad vial size passed" vialsize
    end
    alldims = alldims[1] # Extract object corresponding to row of table
    rad_o = alldims.d1 / 2 * u"mm"
    rad_i = rad_o - alldims.s1 * u"mm"
    return rad_i, rad_o
end

"""
    get_vial_thickness(vialsize::String)

Return vial wall thickness for given ISO vial size.

Uses a table provided by Schott, stored internally in a CSV.
"""
function get_vial_thickness(vialsize::String)
    alldims = [row for row in VIAL_DIMS if row.Size == vialsize]
    if length(alldims) != 1
        @error "bad vial size passed" vialsize
    end
    alldims = alldims[1] # Extract object corresponding to row of table
    thickness = alldims.s1*u"mm"
    return thickness
end