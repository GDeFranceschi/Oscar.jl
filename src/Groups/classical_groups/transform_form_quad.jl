# change of basis from [1,1,0,l1] to [1,1,0,l2] (if exists)
function _matrix_for_2(l1::FieldElem, l2::FieldElem)
   R,t = PolynomialRing(parent(l1), "t")
   f = t^2*(4*l1-1)+t*(1-4*l1)+(l1-l2)
   @assert !isirreducible(f) "The forms are not congruent"
   b = roots(f)[1]
   return matrix(parent(l1),2,2,[1,0,b,1-2*b])
end

# returns z such that z*x*transpose(z) == matrix([1,1,0,something])
# assume x[2,1]=0, x[1,2]!=0
function _matrix_for_2_bis(x::MatElem{T}) where T
   (x[1,1]==0 && x[2,2]==0) && return matrix(base_ring(x),2,2,[1,x[1,2]^-1,1,0])
   x[1,1]==0 && return matrix(base_ring(x),2,2,[0,sqrt(x[2,2]^-1), sqrt(x[2,2])*x[1,2]^-1,0])
   return matrix(base_ring(x),2,2,[sqrt(x[1,1]^-1),0,0,sqrt(x[1,1])*x[1,2]^-1])   
end


# assume x1, x2 have dim=2 and are upper triangular
function _matrix_change_quad(x1::MatElem{T}, x2::MatElem{T}) where T
   z1 = _matrix_for_2_bis(x1)
   z2 = _matrix_for_2_bis(x2)

   l1 = (z1*x1*transpose(z1))[2,2]
   l2 = (z2*x2*transpose(z2))[2,2]
   y = _matrix_for_2(l1,l2)
   return z2^-1*y*z1
end

