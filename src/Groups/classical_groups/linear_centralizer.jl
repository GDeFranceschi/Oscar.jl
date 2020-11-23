
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

# returns as matrices
# does the matrix above with F = F[x]/(f), but every entry is replaced by a diagonal join of D corresponding blocks
# ASSUMPTION: deg(f) > 1
function _gens_for_GL_matrix(f::PolyElem, n::Int, F::Ring; D=1)
   C = companion_matrix(f)
   CP = evaluate(_centralizer(f),C)            # matrix of maximal order in the centralizer of the companion matrix

   if n==1 return [diagonal_join([CP for i in 1:D])] end
   h1 = identity_matrix(F,n*degree(f)*D)
   insert_block!(h1,diagonal_join([CP for i in 1:D]),1,1)
   h2 = zero_matrix(F,n*degree(f)*D,n*degree(f)*D)
   for i in 1:(n-1)*degree(f)*D h2[i+degree(f)*D,i]=-1 end
   for i in 1:degree(f)*D
      h2[i,i]=-1
      h2[i,i+(n-1)*degree(f)*D]=1
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


# V = vector of integers of the dimensions of Jordan blocks
# return the generators for the centralizers of the unipotent element
# assumes V is sorted (e.g. [1,1,1,2,3,3])
# does the same as above, but every entry is replaced by the corresponding block matrix
function _centr_block_unipotent(f::PolyElem, F::Ring, V::AbstractVector{Int})
   d = degree(f)
   d>1 || return _centr_unipotent(F,V)
   n = sum(V)*d
   C = companion_matrix(f)

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
      for x in _gens_for_GL_matrix(f,l[2],F; D=l[1])
         z = insert_block(identity_matrix(F,n),x,pos,pos)
         listgens = cat(listgens,[z]; dims=1)
      end
      if l[1]>1
         for i in 1:l[1]-1
         c = identity_matrix(F,d)
         for j in 1:degree(F)*d
            z = identity_matrix(F,l[1]*d)
#            for k in 1:l[1]-i z[k,i+k]=_lambda^j end
            for k in 0:l[1]-i-1 insert_block!(z,c,d*k+1 ,d*(k+i)+1) end
            z = insert_block(identity_matrix(F,n),z,pos,pos)
            c *= C           # every time, the block C^(j-1) is inserted
            listgens = cat(listgens,[z]; dims=1)
         end
         end
      end
      pos += l[1]*l[2]*d
   end

   # external diagonal blocks
   pos=1
   for i in 1:length(L)-1
      pos += L[i][1]*L[i][2]*d
      # block above diagonal
      z = identity_matrix(F,n)
#      for j in 1:L[i][1] z[pos-L[i][1]+j-1,pos+L[i+1][1]-L[i][1]+j-1]=1 end
      for j in 1:L[i][1]*d z[pos-L[i][1]*d+j-1,pos+(L[i+1][1]-L[i][1])*d+j-1]=1 end
      listgens = cat(listgens,[z]; dims=1)
      # block below diagonal
      z = identity_matrix(F,n)
      for j in 1:L[i][1]*d z[pos+j-1,pos-L[i][1]*d+j-1]=1 end
      listgens = cat(listgens,[z]; dims=1)
   end

   return listgens
end


function _centralizer(x::MatElem)
   _,cbm,ED = generalized_jordan_form(x; with_pol=true)    # cbm = change basis matrix
   n=nrows(x)
   listgens = MatElem[]

   i=1
   pos=1
   f=ED[1][1]
   V=[ED[1][2]]
   while i <= length(ED)
      if i<length(ED) && ED[i+1][1]==f
         i+=1
         V = cat(V,[ED[i][2]]; dims=1)
      else
         for z in _centr_block_unipotent(f,base_ring(x),V)
            listgens = cat(listgens, [insert_block(identity_matrix(base_ring(x),n),z,pos,pos)]; dims=1)
         end
         pos += degree(f)*sum(V)
         i+=1
         if i<=length(ED)
            f = ED[i][1]
            V = [ED[i][2]]
         end
      end
   end

   return listgens, cbm
end



########################################################################
#
# User level functions
#
########################################################################

pol_elementary_divisors(x::MatrixGroupElem) = pol_elementary_divisors(x.elm)

function centralizer(G::MatrixGroup, x::MatrixGroupElem)
   if isdefined(G,:descr) && G.descr==:GL
      V,a = _centralizer(x.elm)
      am = inv(a)
      L = [G(am*v*a) for v in V]
      return MatrixGroup(G.deg, G.ring, L), Nothing          # do not return the embedding of the centralizer into G to do not compute G.X
   end
   C = GAP.Globals.Centralizer(G.X, x.X)
   return _as_subgroup(G, C)
end
