import AbstractAlgebra: FieldElem
import Hecke: gram_matrix

export 
    alternating_form,
    hermitian_form,
    preserved_quadratic_forms,
    preserved_sesquilinear_forms,
    SesquilinearForm,
    symmetric_form

mutable struct SesquilinearForm{T<:RingElem}
   matrix::MatElem{T}
   descr::Symbol       # quadratic, bilinear or hermitian
   pol::AbstractAlgebra.Generic.MPoly{T}     # only for quadratic forms
   X::GapObj

   SesquilinearForm(B::MatElem{T},sym) where T = new{T}(B,sym)
end


GAP.Packages.load("forms")

"""
    preserved_quadratic_forms(G::MatrixGroup)
Return a generating set for the vector space of quadratic forms preserved by `G`.

!!! warning "Note:"
    The process involves random procedures, so the function may return different outputs every time.

"""
function preserved_quadratic_forms(G::MatrixGroup)
   L = GAP.Globals.PreservedQuadraticForms(G.X)
   return [G.mat_iso(GAP.Globals.GramMatrix(L[i])) for i in 1:length(L)]
end

"""
    preserved_sesquilinear_forms(G::MatrixGroup)
Return a generating set for the vector space of sesquilinear forms preserved by `G`.

!!! warning "Note:"
    The process involves random procedures, so the function may return different outputs every time.

"""
function preserved_sesquilinear_forms(G::MatrixGroup)
   L = GAP.Globals.PreservedSesquilinearForms(G.X)
   return [G.mat_iso(GAP.Globals.GramMatrix(L[i])) for i in 1:length(L)]
end

"""
    alternating_form(B::MatElem{T}; check=true)
Return the alternating form with Gram matrix `B`. If `check` is set as `false`, it does not check whether the matrix is skew-symmetric.
"""
function alternating_form(B::MatElem{T}; check=true) where T <: FieldElem
   if check
      for i in 1:nrows(B)
      for j in i:nrows(B)
         @assert B[i,j]==-B[j,i] "The matrix is not skew-symmetric"
      end
      end
      if characteristic(base_ring(B))==2
         for i in 1:nrows(B)
            @assert B[i,i]==0 "The matrix is not skew-symmetric"
         end
      end
   end
   f = SesquilinearForm(B, :alternating)
   return f
end

"""
    symmetric_form(B::MatElem{T}; check=true)
Return the symmetric form with Gram matrix `B`. If `check` is set as `false`, it does not check whether the matrix is symmetric.
"""
function symmetric_form(B::MatElem{T}; check=true) where T <: FieldElem
   if check
      for i in 1:nrows(B)
      for j in i+1:nrows(B)
         @assert B[i,j]==B[j,i] "The matrix is not symmetric"
      end
      end
   end
   f = SesquilinearForm(B, :symmetric)
   return f
end

"""
    hermitian_form(B::MatElem{T}; check=true)
Return the hermitian form with Gram matrix `B`. If `check` is set as `false`, it does not check whether the matrix is hermitian.
"""
function hermitian_form(B::MatElem{T}; check=true) where T <: FieldElem
   if check
      @assert iseven(degree(base_ring(B))) "The matrix is not hermitian"
      e = div(degree(base_ring(B)),2)
      for i in 1:nrows(B)
      for j in i:nrows(B)
         @assert B[i,j]==frobenius(B[j,i],e) "The matrix is not hermitian"
      end
      end
   end
   f = SesquilinearForm(B, :hermitian)
   return f
end


########################################################################
#
# Show
#
########################################################################


function _assign_description(sym::Symbol)
   if sym== :alternating print("Alternating")
   elseif sym== :symmetric print("Symmetric")
   else error("unsupported description")
   end
end


function Base.show(io::IO, B::SesquilinearForm)
   _assign_description(B.descr)
   println(" form with Gram matrix ")
   show(io, B.matrix)
end



########################################################################
#
# Fields of the variable
#
########################################################################

gram_matrix(B::SesquilinearForm) = B.matrix

