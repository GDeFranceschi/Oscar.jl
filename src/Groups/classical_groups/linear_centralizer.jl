
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
function _gens_for_GL(n::Int, F::Ring)
   if n==1 return [matrix(F,1,1,[primitive_element(F)])] end
   if order(F)==2
      h1 = identity_matrix(F,n)
      h1[1,2] = 1
      h2 = zero_matrix(F,n,n)
      for i in 1:n-1 h2[i+1,i]=1 end
      h2[1,n] = 1
   else
      h1 = identity_matrix(F,n)
      h1[1,1] = primitive_element(F)
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
   n = sum(V)
   _lambda = gen(F)  # yes, gen(F) is correct; we don't need a primitive element in this case

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

   # generators for GL + internal diagonal blocks
   pos=1
   for l in L
      for x in _gens_for_GL(l[2],F)
         z = block_matrix(l[2],l[2],[x[i,j]*identity_matrix(F,l[1]) for i in 1:l[2] for j in 1:l[2]])
         z = insert_block(identity_matrix(F,n),z,pos,pos)
         listgens = cat(listgens,[z]; dims=1)
      end
      if l[1]>1
         for i in 1:l[1]-1
         for j in 1:degree(F)
            z = identity_matrix(F,l[1])
            for k in 1:l[1]-i z[k,i+k]=_lambda^j end
            z = insert_block(identity_matrix(F,n),z,pos,pos)
            listgens = cat(listgens,[z]; dims=1)
         end
         end
      end
      pos += l[1]*l[2]
   end

   # external diagonal blocks
   pos=1
   for i in 1:length(L)-1
      pos += L[i][1]*L[i][2]
      # block above diagonal
      z = identity_matrix(F,n)
      for j in 1:L[i][1] z[pos-L[i][1]+j-1,pos+L[i+1][1]-L[i][1]+j-1]=1 end
      listgens = cat(listgens,[z]; dims=1)
      # block below diagonal
      z = identity_matrix(F,n)
      for j in 1:L[i][1] z[pos+j-1,pos-L[i][1]+j-1]=1 end
      listgens = cat(listgens,[z]; dims=1)
   end

   return listgens
end
