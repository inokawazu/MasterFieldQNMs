module BoostedEOMs

import ..MasterFieldEquations: MasterFieldEquation, acoef, bcoef, indicialexponent

mutable struct BoostedEquation{T <: Real} <: MasterFieldEquation
    a::T
    unboosted::MasterFieldEquation
    wtrans::Function
    qtrans::Function
end

function boostedcoef(be::BoostedEquation)::Function
    w = be.wtrans
    q = be.qtrans
    a = be.a
    return fun -> ((u, nu, j) -> fun(u, w(a, nu, j), q(a, nu, j)))
end

function boostedindicial(be::BoostedEquation)::Function
    w = be.wtrans
    q  = be.qtrans
    a = be.a
    return fun -> ((nu, j) -> fun(w(a, nu, j), q(a, nu, j)))
end

acoef(be::BoostedEquation)::Function = boostedcoef(be)(acoef(be.unboosted))
bcoef(be::BoostedEquation)::Function = boostedcoef(be)(bcoef(be.unboosted))

function indicialexponent(be::BoostedEquation)::Function
    return boostedindicial(be)(indicialexponent(be.unboosted))
end

end # end module
