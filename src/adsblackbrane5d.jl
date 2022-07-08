module BlackBrane5D

import ..MasterFieldEquation, ..acoef, ..bcoef, ..indicialexponent

indexp(w, _) = -im*w/2  

struct Shear <: MasterFieldEquation end

function shear_a(u, w, q) 
    f = 1-u^2
    f′= -2*u
    return ((w^2-q^2*f)*f-u*w^2*f′)/(u*(-(u+1))*(q^2*f-w^2))   
end

function shear_b(u, w, q)
    f = 1-u^2
    return (w^2-q^2*f)/(u*(1+u)^2)
end

acoef(_::Shear) = shear_a
bcoef(_::Shear) = shear_b

indicialexponent(_::Shear) = indexp

struct Sound <: MasterFieldEquation end

function sound_a(u, w, q) 
    nu = -(3*w^2*(1+u^2) + q^2*(2*u^2-3*u^4-3))
    de = u*(-(1+u))*(3*w^2+q^2*(u^2-3))
    return nu/de
end

function sound_b(u, w, q)
    nu = 3*w^4+q^4*(3-4*u^2+u^4) + q^2*(4*u^5-4*u^3+4*u^2*w^2-6*w^2)
    de = u*(1+u)^2*(3*w^2+q^2*(u^2-3))
    return nu/de
end

acoef(_::Sound) = sound_a
bcoef(_::Sound) = sound_b

indicialexponent(_::Sound) = indexp

end # BlackBrane5D module
