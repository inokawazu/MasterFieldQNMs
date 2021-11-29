module MasterFieldQNMs


# Begin MasterFieldEquation
abstract type MasterFieldEquation end

acoef(mfe::MasterFieldEquation) = error("$(mfe) must have an 'a(u,w,q)' coefficent.")
bcoef(mfe::MasterFieldEquation) = error("$(mfe) must have an 'b(u,w,q)' coefficent.")
# End MasterFieldEquation

end # module
