<<<<<<< HEAD

# TODO : in this file, some functions for matrices and vectors are defined just to make other files work,
# such as forms.jl, transform_form.jl, linear_conjugate.jl and linear_centralizer.jl
# TODO : functions in this file are only temporarily, and often inefficient.
# TODO: once similar working methods are defined in other files or packages (e.g. Hecke), 
# functions in this file are to be removed / moved / replaced
# TODO: when this happens, files mentioned above need to be modified too.

=======
>>>>>>> GAP: deal with matrices, vectors, finite field elements
import AbstractAlgebra: FieldElem, map, Ring
import Hecke: evaluate, multiplicative_jordan_decomposition, PolyElem, _rational_canonical_form_setup, refine_for_jordan

export
    block_matrix,
    complement,
    conjugate_transpose,
    diagonal_join,
    insert_block,
    insert_block!,
    isconjugate_gl,
    ishermitian_matrix,
    isskewsymmetric_matrix,
<<<<<<< HEAD
=======
    partitions,
>>>>>>> GAP: deal with matrices, vectors, finite field elements
    permutation_matrix,
    submatrix



<<<<<<< HEAD
=======

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

>>>>>>> GAP: deal with matrices, vectors, finite field elements
########################################################################
#
# Matrix manipulation
#
########################################################################

"""
    submatrix(A::MatElem{T}, i::Int, j::Int, m::Int, n::Int)
<<<<<<< HEAD

=======
>>>>>>> GAP: deal with matrices, vectors, finite field elements
Return the `m x n` submatrix of `A` rooted at `(i,j)`
"""
function submatrix(A::MatElem, i::Int, j::Int, nr::Int, nc::Int)
   return matrix(base_ring(A),nr,nc, [A[s,t] for s in i:nr+i-1 for t in j:nc+j-1])
end

# exists already in Hecke _copy_matrix_into_matrix
"""
    insert_block(A::MatElem, B::MatElem, i,j)
<<<<<<< HEAD

=======
>>>>>>> GAP: deal with matrices, vectors, finite field elements
Return the matrix `A` with the block `B` inserted at the position `(i,j)`.
"""
function insert_block(A::MatElem{T}, B::MatElem{T}, i::Int, j::Int) where T <: RingElem
   C = deepcopy(A)
<<<<<<< HEAD
   return insert_block!(C,B,i,j)
=======
   for s in 1:nrows(B)
   for t in 1:ncols(B)
      C[i+s-1,j+t-1] = B[s,t]
   end
   end
   return C
>>>>>>> GAP: deal with matrices, vectors, finite field elements
end

"""
    insert_block!(A::MatElem, B::MatElem, i,j)
<<<<<<< HEAD

Insert the block `B` in the matrix `A` at the position `(i,j)`.
"""
function insert_block!(A::MatElem{T}, B::MatElem{T}, i::Int, j::Int) where T <: RingElem
   for s in 1:nrows(B), t in 1:ncols(B)
      A[i+s-1,j+t-1] = B[s,t]
   end
=======
Insert the block `B` in the matrix `A` at the position `(i,j)`.
"""
function insert_block!(A::MatElem{T}, B::MatElem{T}, i::Int, j::Int) where T <: RingElem
   for s in 1:nrows(B)
   for t in 1:ncols(B)
      A[i+s-1,j+t-1] = B[s,t]
   end
   end
>>>>>>> GAP: deal with matrices, vectors, finite field elements
   return A
end

# exists already in Hecke cat(V...; dims=(1,2))
"""
    diagonal_join(V::AbstractVector{<:MatElem})
    diagonal_join(V::T...) where T <: MatElem
<<<<<<< HEAD

Return the diagonal join of the matrices in `V`.
"""
function diagonal_join(V::AbstractVector{T}) where T <: MatElem
   nr = sum(nrows, V)
   nc = sum(ncols, V)
=======
Return the diagonal join of the matrices in `V`.
"""
function diagonal_join(V::AbstractVector{T}) where T <: MatElem
   nr = sum([nrows(v) for v in V])
   nc = sum([ncols(v) for v in V])
>>>>>>> GAP: deal with matrices, vectors, finite field elements
   B = zero_matrix(base_ring(V[1]), nr,nc)
   pos_i=1
   pos_j=1
   for k in 1:length(V)
      insert_block!(B,V[k],pos_i,pos_j)
      pos_i += nrows(V[k])
      pos_j += ncols(V[k])
   end
   return B
end

diagonal_join(V::T...) where T <: MatElem = diagonal_join(collect(V))

#=
diagonal_join(V::T...) where T <: MatElem = cat(V; dims=(1,2))
=#

"""
    block_matrix(m::Int, n::Int, V::AbstractVector{T}) where T <: MatElem
<<<<<<< HEAD

Given a sequence `V` of matrices, return the `m x n` block matrix, where the `(i,j)`-block is the `((i-1)*n+j)`-th element of `V`. The sequence `V` must have length `mn` and the dimensions of the matrices of `V` must be compatible with the above construction.
=======
Return the matrix constructed from the given block matrices, which should be given as a sequence `V` of `m*n` matrices (given in row major order, in other words listed across rows). 
>>>>>>> GAP: deal with matrices, vectors, finite field elements
"""
function block_matrix(m::Int, n::Int, V::AbstractVector{T}) where T <: MatElem
   length(V)==m*n || throw(ArgumentError("Wrong number of inserted blocks"))
   n_rows=0
   for i in 1:m
      for j in 1:n
         nrows(V[n*(i-1)+j])==nrows(V[n*(i-1)+1]) || throw(ArgumentError("Invalid matrix dimension"))
         ncols(V[n*(i-1)+j])==ncols(V[j]) || throw(ArgumentError("Invalid matrix dimension"))
      end
      n_rows += nrows(V[n*(i-1)+1])
   end
<<<<<<< HEAD
   n_cols = sum(ncols(V[j]) for j in 1:n)
=======
   n_cols = sum([ncols(V[j]) for j in 1:n])
>>>>>>> GAP: deal with matrices, vectors, finite field elements
   B = zero_matrix(base_ring(V[1]), n_rows, n_cols)
   pos_i=1
   for i in 1:m
      pos_j=1
      for j in 1:n
         insert_block!(B,V[n*(i-1)+j],pos_i,pos_j)
         pos_j += ncols(V[n*(i-1)+j])
      end
      pos_i += nrows(V[n*(i-1)+1])
   end
   return B
end

"""
    matrix(A::Array{AbstractAlgebra.Generic.FreeModuleElem{T},1})
<<<<<<< HEAD

=======
>>>>>>> GAP: deal with matrices, vectors, finite field elements
Return the matrix whose rows are the vectors in `A`. Of course, vectors in `A` must have the same length and the same base ring.
"""
function matrix(A::Array{AbstractAlgebra.Generic.FreeModuleElem{T},1}) where T <: FieldElem
   c = length(A[1].v)
<<<<<<< HEAD
   @assert all(x -> length(x.v)==c, A) "Vectors must have the same length"
   X = zero_matrix(base_ring(A[1]), length(A), c)
   for i in 1:length(A), j in 1:c
      X[i,j] = A[i][j]
=======
   X = zero_matrix(base_ring(A[1]), length(A), c)
   for i in 1:length(A)
      @assert length(A[i].v)==c "Vectors must have the same length"
      for j in 1:c X[i,j] = A[i][j] end
>>>>>>> GAP: deal with matrices, vectors, finite field elements
   end

   return X
end

"""
<<<<<<< HEAD
    conjugate_transpose(x::MatElem{T}) where T <: FinFieldElem

If the base ring of `x` is `GF(q^2)`, return the matrix `transpose( map ( y -> y^q, x) )`.
 An error is signalled if the base ring does not have even degree.
"""
function conjugate_transpose(x::MatElem{T}) where T <: FinFieldElem
=======
    conjugate_transpose(x::MatElem{T}) where T <: FieldElem
If the base ring of `x` is `GF(q^2)`, return the matrix `transpose( map ( y -> y^q, x) )`. An error is returned if the base ring has no even degree.
"""
function conjugate_transpose(x::MatElem{T}) where T <: FieldElem
>>>>>>> GAP: deal with matrices, vectors, finite field elements
   iseven(degree(base_ring(x))) || throw(ArgumentError("The base ring must have even degree"))
   e = div(degree(base_ring(x)),2)
   return transpose(map(y -> frobenius(y,e),x))
end


# computes a complement for W in V (i.e. a subspace U of V such that V is direct sum of U and W)
"""
    complement(V::AbstractAlgebra.Generic.FreeModule{T}, W::AbstractAlgebra.Generic.Submodule{T})
<<<<<<< HEAD

=======
>>>>>>> GAP: deal with matrices, vectors, finite field elements
Return a complement for `W` in `V`, i.e. a subspace `U` of `V` such that `V` is direct sum of `U` and `W`.
"""
function complement(V::AbstractAlgebra.Generic.FreeModule{T}, W::AbstractAlgebra.Generic.Submodule{T}) where T <: FieldElem
   @assert issubmodule(V,W) "The second argument is not a subspace of the first one"
   if dim(W)==0 return sub(V,basis(V)) end

   e = W.map

   H = matrix( vcat([e(g) for g in gens(W)], [zero(V) for i in 1:(dim(V)-dim(W)) ]) )
<<<<<<< HEAD
=======
   d = dim(W)
>>>>>>> GAP: deal with matrices, vectors, finite field elements
   A_left = identity_matrix(base_ring(V), dim(V))
   A_right = identity_matrix(base_ring(V), dim(V))
   for rn in 1:dim(W)     # rn = row number
      cn = rn    # column number
      while H[rn,cn]==0 cn+=1 end   # bring on the left the first non-zero entry
      swap_cols!(H,rn,cn)
      swap_rows!(A_right,rn,cn)
      for j in rn+1:dim(W)
         add_row!(H,H[j,rn]*H[rn,rn]^-1,rn,j)
         add_column!(A_left,A_left[j,rn]*A_left[rn,rn]^-1,j,rn)
      end
   end
   for j in dim(W)+1:dim(V)  H[j,j]=1  end
   H = A_left*H*A_right
   _gens = [V([H[i,j] for j in 1:dim(V)]) for i in dim(W)+1:dim(V) ]

   return sub(V,_gens)
end

"""
    permutation_matrix(F::Ring, Q::AbstractVector{T}) where T <: Int
<<<<<<< HEAD
    permutation_matrix(F::Ring, p::PermGroupElem)

Return the permutation matrix over the ring `R` corresponding to the sequence `Q` or to the permutation `p`. If `Q` is a sequence, then `Q` must contain exactly once every integer from 1 to some `n`.
"""
function permutation_matrix(F::Ring, Q::AbstractVector{T}) where T <: Base.Integer
=======
Return the permutation matrix over the ring `R` corresponding to `Q`. Here, `Q` must contain exactly once every integer from 1 to some `n`.
"""
function permutation_matrix(F::Ring, Q::AbstractVector{T}) where T <: Int
>>>>>>> GAP: deal with matrices, vectors, finite field elements
   @assert Set(Q)==Set(1:length(Q)) "Invalid input"
   Z = zero_matrix(F,length(Q),length(Q))
   for i in 1:length(Q) Z[i,Q[i]] = 1 end
   return Z
end

<<<<<<< HEAD
permutation_matrix(F::Ring, p::PermGroupElem) = permutation_matrix(F, listperm(p))

"""
    evaluate(f::PolyElem, X::MatElem)

Evaluate the polynomial `f` in the matrix `X`.
"""
function evaluate(f::PolyElem, X::MatElem)
   B = zero(X)
   C = one(X)
   for i in 0:degree(f)
      B += C*base_ring(X)(coeff(f,i))
      C *= X
   end
   return B
end
=======
"""
    evaluate(f::PolyElem, X::MatElem)
Evaluate the polynomial `f` in the matrix `X`.
"""
evaluate(f::PolyElem, X::MatElem) = sum([X^i*base_ring(X)(coeff(f,i)) for i in 0:degree(f)])
>>>>>>> GAP: deal with matrices, vectors, finite field elements

########################################################################
#
# New properties
#
########################################################################

# TODO: not sure whether this definition of skew-symmetric is standard (for fields of characteristic 2)
"""
    isskewsymmetric_matrix(B::MatElem{T}) where T <: Ring
<<<<<<< HEAD

=======
>>>>>>> GAP: deal with matrices, vectors, finite field elements
Return whether the matrix `B` is skew-symmetric, i.e. `B = -transpose(B)` and `B` has zeros on the diagonal. Returns `false` if `B` is not a square matrix.
"""
function isskewsymmetric_matrix(B::MatElem{T}) where T <: RingElem
   n = nrows(B)
   n==ncols(B) || return false

   for i in 1:n
      B[i,i]==0 || return false
      for j in i+1:n
         B[i,j]==-B[j,i] || return false
      end
   end

   return true
end

"""
    ishermitian_matrix(B::MatElem{T}) where T <: Ring
<<<<<<< HEAD

Return whether the matrix `B` is hermitian, i.e. `B = conjugate_transpose(B)`. Returns `false` if `B` is not a square matrix, or the field has not even degree.
"""
function ishermitian_matrix(B::MatElem{T}) where T <: FinFieldElem
=======
Return whether the matrix `B` is hermitian, i.e. `B = conjugate_transpose(B)`. Returns `false` if `B` is not a square matrix, or the field has not even degree.
"""
function ishermitian_matrix(B::MatElem{T}) where T <: RingElem
>>>>>>> GAP: deal with matrices, vectors, finite field elements
   n = nrows(B)
   n==ncols(B) || return false
   e = degree(base_ring(B))
   iseven(e) ? e = div(e,2) : return false

   for i in 1:n
      for j in i:n
         B[i,j]==frobenius(B[j,i],e) || return false
      end
   end

   return true
end

########################################################################
#
# New operations
#
########################################################################

Base.getindex(V::AbstractAlgebra.Generic.FreeModule, i::Int) = gen(V, i)

# scalar product
Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},u::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = (v.v*transpose(u.v))[1]

Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatElem{T}) where T <: FieldElem = v.parent(v.v*x)
Base.:*(x::MatElem{T},u::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = x*transpose(u.v)

# evaluation of the form x into the vectors v and u
Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatElem{T},u::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = (v.v*x*transpose(u.v))[1]


Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatrixGroupElem{T}) where T <: FieldElem = v.parent(v.v*x.elm)
Base.:*(x::MatrixGroupElem{T},u::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = x.elm*transpose(u.v)

# evaluation of the form x into the vectors v and u
Base.:*(v::AbstractAlgebra.Generic.FreeModuleElem{T},x::MatrixGroupElem{T},u::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = (v.v*x.elm*transpose(u.v))[1]

map(f::Function, v::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: FieldElem = v.parent(map(f,v.v))