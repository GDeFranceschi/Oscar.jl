

########################################################################
#
# MyIsConjugate
#
########################################################################

function pol_elementary_divisors(A::MatElem{T}) where T
   a,_,c = _rational_canonical_form_setup(A)
   L = refine_for_jordan(a,c,A)[2]
   V = Vector(undef, length(L))
   for i in 1:length(L)
      V[i] = (L[i][2],L[i][3])
   end

# sorting the vector in order to have all equal polynomials in a consecutive bunch
# and block size in increasing order
   for i in 1:length(L)-1
   for j in i+1:length(L)
      if V[i][1]==V[j][1]
         V[i+1],V[j] = V[j],V[i+1]
         if V[i+1][2]<V[i][2] V[i],V[i+1] = V[i+1],V[i]  end
      end
   end
   end

   return V
end

function generalized_jordan_block(f::T, n::Int) where T<:PolyElem
   d = degree(f)
   JB = diagonal_join([companion_matrix(f) for i in 1:n])
   pos = 1
   for i in 1:n-1
      insert_block!(JB, identity_matrix(base_ring(f),degree(f)),pos,pos+degree(f))
      pos += degree(f)
   end
   return JB
end

# TODO is there a way to accelerate the process? pol_elementary_divisors and generalized_jordan_block repeat parts of the same code.
function generalized_jordan_form(A::MatElem{T}; with_pol=false) where T
   V = pol_elementary_divisors(A)
   GJ = diagonal_join([generalized_jordan_block(v[1],v[2]) for v in V])
   a = rational_canonical_form(A)[2]
   gj = rational_canonical_form(GJ)[2]
   if with_pol return GJ, gj^-1*a, V
   else return GJ, gj^-1*a
   end
end


function isconjugate_gl(G::MatrixGroup, x::MatrixGroupElem, y::MatrixGroupElem)
   isdefined(G,:descr) || throw(ArgumentError("Group must be general or special linear group"))
   if G.descr==:GL || G.descr==:SL
      Jx,ax = jordan_normal_form(x.elm)
      Jy,ay = jordan_normal_form(y.elm)
      if Jx != Jy return false, nothing end
      z = inv(ax)*ay
      if G.descr==:GL return true, G(z) end
      ED = pol_elementary_divisors(x.elm)
      l = gcd([k[2] for k in ED])
      l = gcd(l, order(G.ring)-1)
      d = det(z)
      if isone(d^( div(order(G.ring)-1,l)))
         corr = _elem_given_det(x, d^-1)
         return true, G(corr*z)
      else return false, nothing
      end
   else
      return isconjugate(G,x,y)
   end
end


function partitions(n::Int, m::Int)         # partitions where the biggest block is at most m
   if n==0 return [Int[]] end
   if m==1 return [[1 for i in 1:n]] end
   L = []
   for i in 0:div(n,m)
      t = [m for i in 1:i]
      for k in partitions(n-m*i, m-1)
         push!(L, vcat(t,k))
      end
   end

   return L
end

partitions(n::Int) = partitions(n,n)

