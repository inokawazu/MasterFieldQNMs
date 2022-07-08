module MasterFieldQNMs

export acoef, bcoef, horexp, hubeny_horowitz_criticalpoint

# main code
include("master_field_equation.jl")
include("transformations.jl")
include("horizon_expansion.jl")

# two use cases
include("adsblackbrane5d.jl")
include("boostedeoms.jl")

end # module
