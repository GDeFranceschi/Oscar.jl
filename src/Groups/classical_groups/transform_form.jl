import AbstractAlgebra: Field

export
    iscongruent

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

# returns x,y such that ax^2+by^2 = c
# at the moment, it simply search for x such that (c-ax^2)/b is a square
# TODO: is there a faster way to find x and y?
# TODO: it would be better if this is deterministic. This depends on gen(F) and issquare(F).
function _solve_eqn(a::T, b::T, c::T) where T <: FieldElem
   F = parent(a)  
   ch = Int(characteristic(F))
   dg = degree(F)
   w = gen(F)
   for i in collect(Iterators.product([0:ch-1 for j in 1:dg]...))
      x = sum([w^(j-1)*i[j] for j in 1:dg])
      s = (c - a*x^2)*b^-1
      vero, y = issquare(s)
      if vero return x,y end
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
      return A*B*conjugate_transpose(A), A, d
   else
      A = transpose(A)
      return B*A, A, d
   end

end





# returns D, A such that A*B*transpose(frobenius(A)) = D and 
# D is diagonal matrix (or with blocks [0 1 s 0])
# f = dimension of the zero block in B in the isotropic case

# based on paper of James Wilson, Optimal algorithms of Gram-Schmidt type
function block_anisotropic_elim(B::MatElem{T}, _type::Symbol; isotr=false, f=0)  where T <: FieldElem

   d = nrows(B)
   F = base_ring(B)
   if d in (0,1)
      return B, identity_matrix(F,d)
   end

   if _type==:symmetric
      degF=0
      s=1
   elseif _type==:alternating
      degF=0
      s=-1
   elseif _type==:hermitian
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



# returns D such that D*B*conjugatetranspose(D) is the standard basis
# it modifies the basis_change_matrix of the function block_herm_elim
# TODO: not done for orthogonal

function _to_standard_form(B::MatElem{T}, _type::Symbol)  where T <: FieldElem
   F = base_ring(B)
   n = nrows(B)
   A,D = block_herm_elim(B, _type)

   if _type==:alternating
      our_perm = vcat(1:2:n, reverse(2:2:n))
      D = permutation_matrix(F,our_perm)*D
   elseif _type==:hermitian
      w = primitive_element(F)
      q = Int(sqrt(order(F)))
      Z = identity_matrix(F,n)
      # turn the elements on the main diagonal into 1
      for i in 1:n
         if A[i,i]!=0
            lambda = _disc_log(w^(q+1),A[i,i])
            Z[i,i] = w^-lambda
            A[i,i] = 1
         end
      end
      D = Z*D
      # moving all hyperbolic lines at the end
      Z = identity_matrix(F,n)
      our_permut = Array(1:n)
      NOZ = 0      # Number Of Zeros on the diagonal before A[i,i]
      for i in 1:n
         if A[i,i]==0
            NOZ += 1
         else
            j = i
            while j>i-NOZ
               swap_cols!(A,j,j-1)
               j -= 1
            end
            j = i
            while j > i-NOZ
               swap_rows!(A,j,j-1)
               j -= 1
            end
            for j in 1:NOZ
               our_permut[i+1-j] = our_permut[i-j]
            end
            our_permut[i-NOZ] = i
         end
      end
      Z = permutation_matrix(F, our_permut)
      D = Z*D
      # turn 2x2 identities into 2x2 anti-diagonal blocks
      Z = identity_matrix(F,n)
      if isodd(q)
         b = (1+w^(div((q-1)^2,2)))^-1
         a = b*w^(div((1-q),2))
         d = 1
         c = w^(div((q-1),2))
      else
         b = (1+w^(q-1))^-1
         a = b*w^(q-1)
         d = 1
         c = 1
      end
      S = matrix(F,2,2,[a,b,c,d])
      if div(n-NOZ,2)==0
         S = zero_matrix(F,0,0)
      else
         S = diagonal_join([S for i in 1:div(n-NOZ,2)])
      end
      # turn into standard GAP form
      sec_perm = Int[]
      for i in 1:div(n,2)
         sec_perm = vcat([i,n+1-i],sec_perm)
      end
      if isodd(n)
         insert_block!(Z,S,2,2)
         sec_perm = vcat([div(n+1,2)],sec_perm)
      else
         insert_block!(Z,S,1,1)
      end
      D = transpose(permutation_matrix(F,sec_perm))*Z*D
   end

   return D
end



###############################################################################################################

# Change of basis between two matrices

###############################################################################################################

# modifies A by eliminating all hyperbolic lines and turning A into a diagonal matrix
# return the matrix Z such that Z*A*transpose(Z) is diagonal
function _elim_hyp_lines(A::MatElem{T}) where T <: FieldElem
   F = base_ring(A)
   n = nrows(A)
   b = matrix(F,2,2,[1,1,1,-1])  # change of basis from matrix([0,1,1,0]) to matrix([2,0,0,-2])
   Z = identity_matrix(F,n)

   i = 1
   while i <= n
      if A[i,i]==0
         A[i,i]=2
         A[i,i+1]=0
         A[i+1,i]=0
         A[i+1,i+1]=-2
         insert_block!(Z,b,i,i)
         i+=2
      else
         i+=1
      end
   end

   return Z
end


# return true, D such that D*B1*conjugatetranspose(D)=B2
# return false, Nothing if D does not exist
# TODO: orthogonal only in odd char, at the moment
function _change_basis(B1::MatElem{T}, B2::MatElem{T}, _type::Symbol)  where T <: FieldElem

   if _type==:alternating || _type==:hermitian
      D1 = _to_standard_form(B1,_type)
      D2 = _to_standard_form(B2,_type)
      return true, D2^-1*D1
   elseif _type==:symmetric
      F = base_ring(B1)
      n = nrows(B1)
      A1,D1 = block_herm_elim(B1, _type)
      A2,D2 = block_herm_elim(B2, _type)
      q = order(F)
      # eliminate all hyperbolic lines and turn A1,A2 into diagonal matrices
      # TODO: assure that the function _elim_hyp_lines actually modifies A1 and A2
      D1 = _elim_hyp_lines(A1)*D1
      D2 = _elim_hyp_lines(A2)*D2
      issquare( prod(diagonal(A1))*prod(diagonal(A2)) )[1] || return false, Nothing
      # move all the squares on the diagonal at the begin
      _squares = [i for i in 1:n if issquare(A1[i,i])[1]]
      our_perm = vcat(_squares, [i for i in 1:n if !(i in _squares)])      # TODO is there a more elengant way?
      P = permutation_matrix(F,our_perm)
      s1 = length(_squares)
      D1 = P*D1
      A1 = P*A1*transpose(P)
      _squares = [i for i in 1:n if issquare(A2[i,i])[1]]
      our_perm = vcat(_squares, [i for i in 1:n if !(i in _squares)])      # TODO is there a more elengant way?
      P = permutation_matrix(F,our_perm)
      s2 = length(_squares)
      D2 = P*D2
      A2 = P*A2*transpose(P)
      # get same number of squares on the two diagonals of A1 and A2 by modifying A1
      if s1!=s2
         s = min(s1,s2)+1
         w = A1[s,s]*A2[s,s]     # I'm sure this is not a square
         a,b = _solve_eqn(F(1),F(1),w)
         L = identity_matrix(F,n)
         for i in 0:div(abs(s1-s2),2)-1
            L[s+2*i,s+2*i] = a
            L[s+2*i,s+2*i+1] = b*issquare(A1[s+2*i,s+2*i]*A1[s+2*i+1,s+2*i+1]^-1)[2]
            L[s+2*i+1,s+2*i] = b
            L[s+2*i+1,s+2*i+1] = -a*issquare(A1[s+2*i,s+2*i]*A1[s+2*i+1,s+2*i+1]^-1)[2]
         end
         D1 = L*D1
         A1 = L*A1*transpose(L)
      end
      # change matrix from A1 to A2
      S = diagonal_matrix([issquare(A2[i,i]*A1[i,i]^-1)[2] for i in 1:n])
      return true, D2^-1*S*D1
   end

end

"""
    iscongruent(f::SesquilinearForm{T}, g::SesquilinearForm{T}) where T <: RingElem
Return (`true`, `C`) if there exists a matrix `C` such that `CAC* = B`, or `CAC*-B` is skew-symmetric (for quadratic forms), where `A` and `B` are the Gram matrices of `f` and `g` respectively, and `C*` is the transpose-conjugate matrix of `C`. If such `C` does not exist, then return (`false`, `Nothing`).
"""
function iscongruent(f::SesquilinearForm{T}, g::SesquilinearForm{T}) where T <: RingElem

   @assert base_ring(f)==base_ring(g) "The forms have not the same base ring"
   @assert nrows(f.matrix)==nrows(g.matrix) "The forms act on vector spaces of different dimensions"
   f.descr==g.descr || return false, Nothing
   
   if f.descr==:quadratic
      if iseven(characteristic(base_ring(f)))            # in this case we use the GAP algorithms
         GAP.Globals.IsometricCanonicalForm(f.X)==GAP.Globals.IsometricCanonicalForm(g.X) || return false, Nothing
         B1 = GAP.Globals.BaseChangeToCanonical(f.X)
         B2 = GAP.Globals.BaseChangeToCanonical(g.X)
         return true, f.mat_iso(B2)^-1*f.mat_iso(B1)
      else
         return _change_basis(f.matrix+transpose(f.matrix), g.matrix+transpose(g.matrix), :symmetric)
      end
   else
      return _change_basis(f.matrix, g.matrix, f.descr)
   end
end
