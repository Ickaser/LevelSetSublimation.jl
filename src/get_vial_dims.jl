export get_vial_rad
# Get dimensions from here: https://www.schott.com/en-us/products/vials/-/media/project/onex/products/v/vials/application-variants/schott-brochure-schott-vials-english-20092017.pdf

# Added to "vial_sizes.csv", read in here

function get_vial_rad(vialsize::String)
    alldims = [row for row in CSV.File(srcdir("vial_sizes.csv")) if row.Size==vialsize]
    if length(alldims) != 1
        @error "bad vial size passed" vialsize
    end
    alldims = alldims[1] # Extract object corresponding to row of table
    rad_o = alldims.d1 / 2 * u"mm"
    rad_i = rad_o - alldims.s1*u"mm"
    
end