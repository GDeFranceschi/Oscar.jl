
#return true, m if f and g are congruent, and m is such that g = mfm*
function iscongruent(f::SesquilinearForm, g::SesquilinearForm)
   base_ring(f)==base_ring(g) || return false, Nothing      # TODO: do we want to return an ERROR in these cases?
   nrows(f.matrix)==nrows(g.matrix) || return false, Nothing
   f.descr==g.descr || return false, Nothing

   mf = f.mat_iso(GAP.Globals.BaseChangeToCanonical(f.X))
   mg = g.mat_iso(GAP.Globals.BaseChangeToCanonical(g.X))
   if f.descr==:hermitian
      e = div(degree(base_ring(f)),2)
      Vero = mf*f.matrix*transpose(frobenius(mf,e))==mg*g.matrix*transpose(frobenius(mg,e))
   else
      Vero = mf*f.matrix*transpose(mf)==mg*g.matrix*transpose(mg)
   end
   if Vero
      return true, mg^-1*mf
   else
      return false, Nothing
   end
end



###############################################################################################################

# Algorithm

###############################################################################################################

# TODO to be removed when the complement function is written
# computes the orthogonal complement (w.r.t the standard inner product) for W in V
function complement(V::AbstractAlgebra.Generic.FreeModule{T}, W::AbstractAlgebra.Generic.Submodule{T}) where T <: FieldElem
   @assert issubmodule(V,W) "The second argument is not a subspace of the first one"

   e = W.map
   H = zero_matrix(base_ring(V), dim(W), dim(V))
   for i in 1:dim(W)
      v = e(gen(W,i))
      for j in 1:dim(V)  H[i,j] = v[j]   end
   end
   d, K = kernel(H)
   _gens = [ sum([K[j,i]*gen(V,j) for j in 1:dim(V)]) for i in 1:d ]
   return sub(V,_gens)
end


# returns C,A,d where B*A = C, C = [C 0] and d = Rank(C)
function radical(B::MatElem{T}, F::Field, n::Int, m::Int; e=0, Symmetric=false) where T <: FieldElem

   V = VectorSpace(F,m);
   K = Symmetric ? kernel(ModuleHomomorphism(V,V,B))[1] : kernel(ModuleHomomorphism(V,V,transpose(B)))[1]
   U,emb = complement(V,K)
   d = dim(U)
   A = zero_matrix(F,m,m)
   for i in 1:d
   for j in 1:m
      A[i,j] = emb(gen(U,i)[j])
   end
   end
   for i in 1:m-d
   for j in 1:m
      A[d+i,j] = K.map(gen(K,i))[j]
   end
   end

   if Symmetric then
      return A*B*transpose(map(y -> frobenius(y,e),A)), A, d
   else
      A = transpose(A)
      return B*A, A, d
   end

end
