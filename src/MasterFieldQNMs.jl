module MasterFieldQNMs


# Begin MasterFieldEquation
abstract type MasterFieldEquation{ T } end

acoef(mfe::MasterFieldEquation) = error("$(mfe) must have an 'a(u,w,q)' coefficent.")
bcoef(mfe::MasterFieldEquation) = error("$(mfe) must have an 'b(u,w,q)' coefficent.")
indicialexponent(mfe::MasterFieldEquation) = error("$(mfe) must have an indicial exponent, α")

horizonlocation(_::MasterFieldEquation{T}) where T = one(T)
# End MasterFieldEquation

# Begin Transformation functions
function atrans(a::Function, α::Number)
    return (x...,) -> a(x...) + 2*α
end 

function btrans(b::Function, a::Function, α::Number)
    return (x...,) -> b(x...) + a(x...)*α + α*(α-1)
end 
# End Transformation functions

# Begin horizon expansion
using TaylorSeries: Taylor1
using NLsolve

function horexp(mfe::MasterFieldEquation{T}, w::Complex{T}, q::Complex{T}; hororder = 20) where T
  horloc = horizonlocation(mfe)
  α = indicialexponent(mfe) # -im*w/2  ingoing boundary condition at t  he horizon
  t = Taylor1(T, hororder)
  a = acoef(mfe)
  b = bcoef(mfe)

  ahor = atrans(a,α)(t+horloc, w, q)
  bhor = btrans(b, a, α)(t+horloc, w, q)

  lhs = T[
    ( m == n    ? (n-1)*n                                         : 0.0) +
    ((n+1) > m  ? ahor.coeffs[(n+1)-m]*m + bhor.coeffs[(n+1)-m]   : 0.0)
    for n in 1:hororder, m in 1:hororder
   ]

  rhs = T[-bhor.coeffs[n+1] for n in 1:hororder]

  return Taylor1(T[1; lhs\rhs])
end

function hubeny_horowitz_criticalpoint(
    mfe::MasterFieldEquation{T}, w0::Complex{T}, q0::Complex{T}; hororder = 20
  ) where T

  deriv_diff(f, h=sqrt(4*eps(T))) = x -> sum(
                                             co*f(x + h*(n-1)) 
                                             for (n, co) 
                                             in enumerate([-137/60, 5, -5, 10/3, -5/4, 1/5])
                                            ) / h

    function g!(F, x)
      w′,q′ = x

      function f(w) 
        he = horexp(mfe, w, q′; hororder=hororder)
        return Taylor1(he)(-1)  
      end

      F[1] = f(w′)
      F[2] = deriv_diff(f)(w′)
      return F
    end

    nlsolve(g!, [w0, q0])
end

# End horizon expansion

end # module
