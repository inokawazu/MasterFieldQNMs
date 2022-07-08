module MasterFieldQNMs

export # transformations
       acoef, bcoef, 
       # horizon_expansion
       horexp, hubeny_horowitz_criticalpoint, hubeny_horowitz_qnm

using TaylorSeries: Taylor1 # horizon_expansion
using NLsolve # horizon_expansion
using Roots: muller # horizon_expansion 

# main code
include("master_field_equation.jl")
include("transformations.jl")
include("horizon_expansion.jl")

# two use cases
include("adsblackbrane5d.jl")
include("boostedeoms.jl")

end # module
