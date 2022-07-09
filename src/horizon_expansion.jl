
function horexp(
    mfe::MasterFieldEquation, w::Number, q::Number; hororder = 20
  )
    w, q = promote(w, q)
    T = real(eltype(w))

    horloc::T = horizonlocation(mfe)
    α::Complex{T}  = indicialexponent(mfe)(w, q) # ingoing boundary condition at the horizon
    t = Taylor1(Complex{T}, hororder)
    a = acoef(mfe)
    b = bcoef(mfe)

    ahor = atrans(a,α)(t+horloc, w, q)
    bhor = btrans(b,a,α)(t+horloc, w, q)

    out_series = Array{Complex{T}}(undef, hororder)
    out_series[1] = one(Complex{T})

    for ord in 2:hororder
        m = ord - 1
        s = sum(
                (bhor.coeffs[ord - k] + (k-1)*ahor.coeffs[ord - k])*out_series[k] 
                for k in 1:ord-1
               )
        out_series[ord] = - s/(m*(m-1) + m*ahor.coeffs[1])
    end

    return out_series
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

      # TODO: remove hard coded "boundary"
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

# TODO: convergence is weak.
function hubeny_horowitz_qnm(
    mfe::MasterFieldEquation, w0::Number, q::Number; hororder = 20,
    xatol=nothing, xrtol=nothing, maxevals=1000
    ) 

    # TODO: remove hard coded "boundary"
    function f(freq) 
        he = horexp(mfe, freq, q; hororder=hororder)
        return Taylor1(he)(-1)  
    end

    muller(f, w0; xatol=xatol, xrtol=xrtol, maxevals=maxevals)
end
