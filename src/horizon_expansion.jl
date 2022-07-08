module HorizonExpansion # Begin horizon expansion

export horexp, hubeny_horowitz_criticalpoint
using ..MasterFieldEquations, ..Transformations

using TaylorSeries: Taylor1
using NLsolve

function horexp(
    mfe::MasterFieldEquation, w::Complex{T}, q::Complex{T}; hororder = 20
  ) where T <: Number
  horloc = horizonlocation(mfe) |> T
  α = indicialexponent(mfe)(w, q) |> Complex{T} # ingoing boundary condition at the horizon
  t = Taylor1(Complex{T}, hororder)
  a = acoef(mfe)
  b = bcoef(mfe)

  ahor = atrans(a,α)(t+horloc, w, q)
  bhor = btrans(b,a,α)(t+horloc, w, q)

  lhs = Complex{T}[
                   (m == n    ? T((n-1)*n)                                      : T(0.0)) +
                   ((n+1) > m  ? ahor.coeffs[(n+1)-m]*m + bhor.coeffs[(n+1)-m]   : T(0.0))
                   for n in 1:hororder, m in 1:hororder
                  ]

  rhs = Complex{T}[-bhor.coeffs[n+1] for n in 1:hororder]

  return Taylor1(Complex{T}[Complex{T}(1); lhs\rhs])
end

function hubeny_horowitz_criticalpoint(
    mfe::MasterFieldEquation, w0::Complex{T}, q0::Complex{T}; hororder = 20
  ) where T <: Number

  deriv_diff(f, h=sqrt(eps(T))) = x -> sum(
                                             co*f(x + h*(n-1)) 
                                             for (n, co) 
                                             in enumerate(T[-137/60, 5, -5, 10/3, -5/4, 1/5])
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

end # End horizon expansion
