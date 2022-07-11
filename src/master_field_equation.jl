abstract type MasterFieldEquation end

acoef(mfe::MasterFieldEquation) = error("$(mfe) must have an 'a(u,w,q)' coefficent.")

bcoef(mfe::MasterFieldEquation) = error("$(mfe) must have an 'b(u,w,q)' coefficent.")

ccoef(_::MasterFieldEquation) = (x, _...,) -> one(x)

indicialexponent(mfe::MasterFieldEquation) = error("$(mfe) must have an indicial exponent, α(w,q)")

horizonlocation(_::MasterFieldEquation) = 1
