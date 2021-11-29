module MasterFieldQNMs


# Begin MasterFieldEquation
abstract type MasterFieldEquation{ T } end

acoef(mfe::MasterFieldEquation) = error("$(mfe) must have an 'a(u,w,q)' coefficent.")
bcoef(mfe::MasterFieldEquation) = error("$(mfe) must have an 'b(u,w,q)' coefficent.")

horizonlocation(_::MasterFieldEquation{T}) where T = one(T)
# End MasterFieldEquation

function atrans(a::Function, α::Number)
    return (x...,) -> a(x...) + 2*α
end 

function btrans(b::Function, a::Function, α::Number)
    return (x...,) -> b(x...) + a(x...)*α + α*(α-1)
end 

end # module
