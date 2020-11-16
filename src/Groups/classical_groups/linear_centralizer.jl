
# return an element in the centralizer of x in GL(n,F) with determinant d

# first: brute force way
function _elem_given_det(x,d)
   C,e = centralizer(GL(x.parent.deg, x.parent.ring),x)
   U,fa = unit_group(x.parent.ring)
   GA,ea = sub(U, [preimage(fa,det(g)) for g in gens(C)])
   l = preimage(ea,preimage(fa,d))
   return prod([C[i]^Int(l[i]) for i in 1:ngens(C)])
end

# returns as matrices
# TODO: uses gen(F) for a generator of the multiplicative group of F
function _gens_for_GL(n::Int, F::Ring)
   if n==1 return matrix(F,1,1,[gen(F)]) end
   if order(F)==2
      h1 = identity_matrix(F,n)
      h1[1,2] = 1
      h2 = zero_matrix(F,n,n)
      for i in 1:n-1 h2[i+1,i]=1 end
      h2[1,n] = 1
   else
      h1 = identity_matrix(F,n)
      h1[1,1] = gen(F)
      h2 = zero_matrix(F,n,n)
      for i in 1:n-1 h2[i+1,i]=-1 end
      h2[1,1] = -1
      h2[1,n] = 1
   end
   return h1,h2      
end

# V = vector of integers of the dimensions of Jordan blocks
# return the generators for the centralizers of the unipotent element
# assumes V is sorted (e.g. [1,1,1,2,3,3])
function _centr_unipotent(F::Ring, V::AbstractVector{Int}) 
   
   # L = multiset(V)
   L=[[V[1],1]]
   for i in 2:length(V)
      if V[i]==L[length(L)][1]
         L[length(L)][2] += 1
      else
         L = cat(L, [[V[i],1]]; dims=1)
      end
   end
   listgens = MatElem[]
   
   

end
