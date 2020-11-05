import Hecke: multiplicative_jordan_decomposition

export
    issemisimple,
    isunipotent,
    multiplicative_jordan_decomposition


########################################################################
#
# Semisimple / Unipotent
#
########################################################################

function multiplicative_jordan_decomposition(x::MatrixGroupElem)
   a,b = multiplicative_jordan_decomposition(x.elm)
   return MatrixGroupElem(x.parent,a), MatrixGroupElem(x.parent,b)
end

issemisimple(x::MatrixGroupElem) = iscoprime(Int(order(x)), Int(characteristic(x.parent.ring)))
isunipotent(x::MatrixGroupElem) = isone(x) || ispower(Int(order(x)))[2]==Int(characteristic(x.parent.ring))
