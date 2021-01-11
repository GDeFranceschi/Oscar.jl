import AbstractAlgebra: Field

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
# based on paper of James Wilson, Optimal algorithms of Gram-Schmidt type

###############################################################################################################

# returns x,y such that ax^2+by^2 = c
function _solve_eqn(a::T, b::T, c::T) where T <: FieldElem
   F = parent(a)
   found = false
   w = primitive_element(F)
   x = one(F)
   while found==false
      x *= w
      s = (c - a*x^2)*b^-1
      vero, y = issquare(s)
      if vero
         found = true
         return x,y
      end
   end
end



# if Symmetric, returns C,A,d where A*B*transpose(frobenius(A,e)) = C, C = block_matrix(2,2,[C,0,0,0]) and d = rank(C)
# else returns C,A,d where B*A = C, C = [C 0] and d = rank(C)
# Assumption: if Symmetric==true, then nr=nc always
function find_radical(B::MatElem{T}, F::Field, nr::Int, nc::Int; e=0, Symmetric=false) where T <: FieldElem

   V1 = VectorSpace(F,nc)
   V2 = VectorSpace(F,nr)
   K = Symmetric ? kernel(ModuleHomomorphism(V1,V2,B))[1] : kernel(ModuleHomomorphism(V1,V2,transpose(B)))[1]
   U,emb = complement(V1,K)
   d = dim(U)
   A = matrix(vcat(typeof(emb(gen(U,1)))[emb(v) for v in gens(U)], typeof(K.map(gen(K,1)))[K.map(v) for v in gens(K)] ))
#=   A = zero_matrix(F,m,m)
   for i in 1:d
   for j in 1:m
      A[i,j] = emb(gen(U,i))[j]
   end
   end
   for i in 1:m-d
   for j in 1:m
      A[d+i,j] = K.map(gen(K,i))[j]
   end
   end
=#

   if Symmetric
      return A*B*transpose(map(y -> frobenius(y,e),A)), A, d
   else
      A = transpose(A)
      return B*A, A, d
   end

end





# returns D, A such that A*B*transpose(frobenius(A)) = D and 
# D is diagonal matrix (or with blocks [0 1 s 0])
# f = dimension of the zero block in B in the isotropic case
function block_anisotropic_elim(B::MatElem{T}, _type; isotr=false, f=0)  where T <: FieldElem

   d = nrows(B)
   F = base_ring(B)
   if d in (0,1)
      return B, identity_matrix(F,d)
   end

   if _type=="orthogonal"
      degF=0
      s=1
   elseif _type=="symplectic"
      degF=0
      s=-1
   elseif _type=="unitary"
      degF=div(degree(F),2)
      s=1
   end

   # conjugate transpose
   star(X; exp=degF) = transpose(map(y -> frobenius(y,exp),X))

   if isotr
      q = characteristic(F)^degF
      g = d-f
      U = submatrix(B,1,f+1,f,g)
      V = submatrix(B,f+1,f+1,g,g)
      C,A,e = find_radical(U,F,f,g)
      # I expect C always to be of rank f
      C = submatrix(C,1,1,f,f)                     
      Vprime = star(A)*V*A
      Z = submatrix(Vprime, 1,1,f,f)
      Y = submatrix(Vprime,f+1,1,g-f,f)
      Bprime = submatrix(Vprime,f+1,f+1,g-f,g-f)
      TR = zero_matrix(F,f,f)
      D = zero_matrix(F,f,f)
      for i in 1:f
         D[i,i] = Z[i,i]
         for j in i+1:f  TR[i,j] = Z[i,j] end
      end
      Aarray = MatElem{T}[]  #TODO type?
      Barray = MatElem{T}[]
      for i in 1:f
         alpha = D[i,i]
         if alpha != 0
            push!(Barray, matrix(F,2,2,[-alpha^-1,0,0,alpha]))
            push!(Aarray, matrix(F,2,2,[1,-alpha^-1,0,1]))
         else
            push!(Barray, matrix(F,2,2,[0,1,s,0]))
            push!(Aarray, matrix(F,2,2,[1,0,0,1]))
         end
      end
      B0,A0 = block_anisotropic_elim(Bprime,_type)
      B1 = diagonal_join(Barray)
      B1 = diagonal_join(B1,B0)
      C = C^-1
      Temp = vcat(C,-TR*C)
      Temp = vcat(Temp,-Y*C)
      Temp = hcat(Temp, vcat(zero_matrix(F,f,g),star(A)))
      P = zero_matrix(F,2*f,2*f)
      for i in 1:f
         P[2*i-1,i] = 1
         P[2*i,i+f] = 1
      end
      A1 = diagonal_join(Aarray)*P
      A1 = diagonal_join(A1,A0)
      return B1, A1*Temp
   else
      c,f = Int(ceil(d/2)), Int(floor(d/2))
      B0 = submatrix(B,1,1,c,c)
      U = submatrix(B,c+1,1,f,c)
      V = submatrix(B,c+1,c+1,f,f)
      B1,A0,e = find_radical(B0,F,c,c; e=degF, Symmetric=true)
      B1 = submatrix(B1,1,1,e,e)
      U = U*star(A0)
      U1 = submatrix(U,1,1,f,e)
      U2 = submatrix(U,1,e+1,f,c-e)
      Z = V-s*U1*B1^-1*star(U1)
      D1,A1 = block_anisotropic_elim(B1,_type)
      Temp = zero_matrix(F,d-e,d-e)
      insert_block!(Temp,s*star(U2),1,c-e+1)
      insert_block!(Temp,U2,c-e+1,1)
      insert_block!(Temp,Z,c-e+1,c-e+1)
      if c-e==0
         D2,A2 = block_anisotropic_elim(Temp,_type)
      else
         D2,A2 = block_anisotropic_elim(Temp, _type; isotr=true, f=c-e)
      end
      Temp = hcat(-U1*B1^-1, zero_matrix(F,f,c-e))*A0
      Temp = vcat(A0,Temp)
      Temp = insert_block(identity_matrix(F,d),Temp,1,1)

      return diagonal_join(D1,D2), diagonal_join(A1,A2)*Temp
   end
end



# assume B is nondegenerate


# returns D, A such that A*B*transpose(frobenius(A)) = D and 
# D is diagonal matrix (or with blocks [0 1 s 0])
# f = dimension of the zero block in B in the isotropic case
function block_herm_elim(B::MatElem{T}, _type) where T <: FieldElem
   d = nrows(B)
   F = base_ring(B)

   if d==1
      return B, identity_matrix(F,1)
   end

   c = Int(ceil(d/2))
   B2 = submatrix(B,1,1,c,c)
   if B2==0
      D,A = block_anisotropic_elim(B,_type; isotr=true, f=c)
   else
      D,A = block_anisotropic_elim(B,_type)
   end

   return D,A
end
