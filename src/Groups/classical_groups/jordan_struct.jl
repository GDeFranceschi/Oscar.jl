import Hecke: evaluate, multiplicative_jordan_decomposition, PolyElem, _rational_canonical_form_setup, refine_for_jordan

export
    block_matrix,
    diagonal_join,
    generalized_jordan_form,
    insert_block,
    insert_block!,
    isconjugate_gl,
    issemisimple,
    isunipotent,
    multiplicative_jordan_decomposition,
    pol_elementary_divisors,
    submatrix


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








########################################################################
#
# Matrix manipulation
#
########################################################################

"""
    submatrix(A::MatElem{T}, m::Int, n::Int, i::Int, j::Int)

Return the `m x n` submatrix of `A` rooted at `(i,j)`
"""
function submatrix(A::MatElem, nr::Int, nc::Int, i::Int, j::Int)
   return matrix(base_ring(A),nr,nc, [A[s,t] for s in i:nr+i-1 for t in j:nc+j-1])
end

# exists already in Hecke _copy_matrix_into_matrix
"""
    insert_block(A::MatElem, B::MatElem, i,j)

Return the matrix `A` with the block `B` inserted at the position `(i,j)`.
"""
function insert_block(A::MatElem{T}, B::MatElem{T}, i::Int, j::Int) where T <: RingElem
   C = deepcopy(A)
   for s in 1:nrows(B)
   for t in 1:ncols(B)
      C[i+s-1,j+t-1] = B[s,t]
   end
   end
   return C
end

"""
    insert_block!(A::MatElem, B::MatElem, i,j)

Insert the block `B` in the matrix `A` at the position `(i,j)`.
"""
function insert_block!(A::MatElem{T}, B::MatElem{T}, i::Int, j::Int) where T <: RingElem
   for s in 1:nrows(B)
   for t in 1:ncols(B)
      A[i+s-1,j+t-1] = B[s,t]
   end
   end
   return A
end

# exists already in Hecke cat(V...; dims=(1,2))
"""
    diagonal_join(V::AbstractVector{<:MatElem})
    diagonal_join(V::T...) where T <: MatElem

Return the diagonal join of the matrices in `V`.
"""
function diagonal_join(V::AbstractVector{T}) where T <: MatElem
   nr = sum([nrows(v) for v in V])
   nc = sum([ncols(v) for v in V])
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

Return the matrix constructed from the given block matrices, which should be given as a sequence `V` of `m*n` matrices (given in row major order, in other words listed across rows). 
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
   n_cols = sum([ncols(V[j]) for j in 1:n])
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



########################################################################
#
# Misc
#
########################################################################

function evaluate(f::PolyElem, x::MatElem)
   return sum([coefficients(f)[i]*x^i for i in 0:degree(f)])
end
