import AbstractAlgebra: FieldElem, Generic.MPoly
import Hecke: base_ring, defining_polynomial, gram_matrix

export 
    alternating_form,
    corresponding_bilinear_form,
    corresponding_quadratic_form,
    hermitian_form,
    isalternating_form,
    ishermitian_form,
    isquadratic_form,
    issymmetric_form,
    preserved_quadratic_forms,
    preserved_sesquilinear_forms,
    quadratic_form,
    SesquilinearForm,
    symmetric_form



# descr is always defined
# matrix is always defined except when descr="quadratic"; in such a case, at least one of matrix and pol is defined
# NOTE: the fields ring_iso and mat_iso are always defined if the field X is

mutable struct SesquilinearForm{T<:RingElem}
   matrix::MatElem{T}
   descr::Symbol       # quadratic, symmetric, alternating or hermitian
   pol::AbstractAlgebra.Generic.MPoly{T}     # only for quadratic forms
   X::GapObj
   ring_iso::GenRingIso
   mat_iso::GenMatIso

   SesquilinearForm(B::MatElem{T},sym) where T = new{T}(B,sym)

   function SesquilinearForm(f::MPoly{T},sym) where T
      @assert sym==:quadratic "Only quadratic forms are described by polynomials"
      r = new{T}()
      r.pol = f
      r.descr = :quadratic
      return r
   end
end





########################################################################
#
# Properties
#
########################################################################


# return whether B is skew-symmetric matrix
function _is_skew_symmetric(B::MatElem)
   for i in 1:nrows(B)
   for j in i:nrows(B)
      if B[i,j] != -B[j,i] return false end
   end
   end
   if characteristic(base_ring(B))==2
      for i in 1:nrows(B)
         if B[i,i]!=0 return false end
      end
   end

   return true
end

# return whether B is symmetric matrix
function _is_symmetric(B::MatElem)
   for i in 1:nrows(B)
   for j in i+1:nrows(B)
      if B[i,j]!=B[j,i] return false end
   end
   end

   return true
end

# return whether B is hermitian matrix
function _is_hermitian(B::MatElem)
   if isodd(degree(base_ring(B))) return false end
   e = div(degree(base_ring(B)),2)
   for i in 1:nrows(B)
   for j in i:nrows(B)
      if B[i,j]!=frobenius(B[j,i],e) return false end
   end
   end

   return true
end

"""
    isalternating_form(f::SesquilinearForm)
Return whether the form `f` is an alternating form.
"""
isalternating_form(f::SesquilinearForm) = f.descr==:alternating

"""
    ishermitian_form(f::SesquilinearForm)
Return whether the form `f` is a hermitian form.
"""
ishermitian_form(f::SesquilinearForm) = f.descr==:hermitian

"""
    isquadratic_form(f::SesquilinearForm)
Return whether the form `f` is a quadratic form.
"""
isquadratic_form(f::SesquilinearForm) = f.descr==:quadratic

"""
    issymmetric_form(f::SesquilinearForm)
Return whether the form `f` is a symmetric form.
"""
issymmetric_form(f::SesquilinearForm) = f.descr==:symmetric

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




########################################################################
#
# Constructors
#
########################################################################


"""
    alternating_form(B::MatElem{T}; check=true)
Return the alternating form with Gram matrix `B`. If `check` is set as `false`, it does not check whether the matrix is skew-symmetric.
"""
function alternating_form(B::MatElem{T}; check=true) where T <: FieldElem
   if check
      @assert _is_skew_symmetric(B) "The matrix is not skew-symmetric"
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
      @assert _is_symmetric(B) "The matrix is not symmetric"
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
      @assert _is_hermitian(B) "The matrix is not hermitian"
   end
   f = SesquilinearForm(B, :hermitian)
   return f
end

# turns the matrix of a quadratic form into an upper triangular matrix of the same form
# (two matrices A,B represent the same quadratic form iff A-B is skew-symmetric)
function _upper_triangular_version(C::MatElem)
   B = deepcopy(C)
   for i in 1:nrows(B)
   for j in i+1:nrows(B)
      B[i,j]+=B[j,i]
      B[j,i]=0
   end
   end
   return B
end

"""
    quadratic_form(B::MatElem{T})
Return the quadratic form with Gram matrix `B`.
"""
function quadratic_form(B::MatElem{T}; check=true) where T <: FieldElem
   f = SesquilinearForm(_upper_triangular_version(B), :quadratic)
   return f
end

"""
    quadratic_form(f::MPoly{T})
Return the quadratic form described by the polynomial `f`. Here, `f` must be a homogeneous polynomial of degree 2.
"""
function quadratic_form(f::MPoly{T}; check=true) where T <: FieldElem
   @assert total_degree(f)==2 "The polynomial must have degree 2"
   @assert ishomogeneous(f) "The polynomial is not homogeneous"
   
   f = SesquilinearForm(f, :quadratic)
   return f
end


########################################################################
#
# Show
#
########################################################################


function _assign_description(sym::Symbol)
   if sym== :alternating print("Alternating")
   elseif sym== :hermitian print("Hermitian")
   elseif sym== :symmetric print("Symmetric")
   elseif sym== :quadratic print("Quadratic")
   else error("unsupported description")
   end
end


function Base.show(io::IO, B::SesquilinearForm)
   _assign_description(B.descr)
   println(" form with Gram matrix ")
   show(io, "text/plain", B.matrix)
end




########################################################################
#
# Basic
#
########################################################################

#TODO: checking whether two quadratic forms coincide by checking their polynomials is not possible yet.
==(B::SesquilinearForm, C::SesquilinearForm) = B.matrix==C.matrix && B.descr==C.descr

function base_ring(B::SesquilinearForm)
   if isdefined(B,:matrix) return base_ring(B.matrix)
   else return base_ring(B.pol)
   end
end

"""
    corresponding_bilinear_form(Q::SesquilinearForm)
Given a quadratic form `Q`, return the bilinear form `B` defined by `B(u,v) = Q(u+v)-Q(u)-Q(v)`.
"""
function corresponding_bilinear_form(B::SesquilinearForm)
   B.descr==:quadratic || throw(ArgumentError("The form must be a quadratic form"))
   M = B.matrix+transpose(B.matrix)
   if characteristic(base_ring(B))==2 return alternating_form(M)
   else return symmetric_form(M)
   end
end

"""
    corresponding_quadratic_form(Q::SesquilinearForm)
Given a symmetric form `f`, returns the quadratic form `Q` defined by `Q(v) = f(v,v)/2`. It is defined only in odd characteristic.
"""
function corresponding_quadratic_form(B::SesquilinearForm)
   B.descr==:symmetric || throw(ArgumentError("The form must be a symmetric form"))
   characteristic(base_ring(B))!=2 || throw(ArgumentError("Corresponding quadratic form not uniquely determined"))
   M = deepcopy(B.matrix)
   l = inv(base_ring(B)(2))
   for i in 1:nrows(M)
      for j in i+1:nrows(M)
         M[j,i]=0
      end
      M[i,i]*=l
   end
   return quadratic_form(M)
end



########################################################################
#
# Fields of the variable
#
########################################################################

"""
    gram_matrix(B::SesquilinearForm)
Return the Gram matrix of a sesquilinear or quadratic form `B`.
"""
gram_matrix(B::SesquilinearForm) = B.matrix

"""
    defining_polynomial(f::SesquilinearForm)
Return the polynomial that defines the quadratic form `f`.
"""
defining_polynomial(B::SesquilinearForm) = B.pol


function assign_from_description(f::SesquilinearForm)
   if f.descr==:quadratic f.X=GAP.Globals.QuadraticFormByMatrix(f.mat_iso(f.matrix))
   elseif f.descr==:symmetric || f.descr==:alternating f.X=GAP.Globals.BilinearFormByMatrix(f.mat_iso(f.matrix))
   elseif f.descr==:hermitian f.X=GAP.Globals.HermitianFormByMatrix(f.mat_iso(f.matrix),f.ring_iso.codomain)
   else error("unsupported description")
   end
end


function Base.getproperty(f::SesquilinearForm, sym::Symbol)

   if isdefined(f,sym) return getfield(f,sym) end

   if sym === :matrix && f.descr==:quadratic  # assume only the polynomial is defined
      d = nvars(parent(f.pol))
      B = zero_matrix( base_ring(f.pol), d, d )
      V = collect(exponent_vectors(f.pol))
      C = collect(coeffs(f.pol))
      for i in 1:length(V)
         for j in 1:d
            if V[i][j] !=0
               global x = j
               break
            end
         end
         for j in 1:d
            if V[i][d+1-j] !=0
               global y = d+1-j
               break
            end
         end
         B[x,y] = C[i]
      end
      f.matrix = B

   elseif sym === :pol
      @assert f.descr == :quadratic "Polynomial defined only for quadratic forms"
      R = PolynomialRing(base_ring(f.matrix), nrows(f.matrix) )[1]
      p = zero(R)
      for i in 1:nrows(f.matrix)
      for j in i:nrows(f.matrix)
         p += f.matrix[i,j] * R[i]*R[j]
      end
      end
      f.pol = p

   elseif sym === :mat_iso || sym === :ring_iso
      fm = gen_mat_iso(nrows(f.matrix), base_ring(f.matrix))
      f.ring_iso=fm.fr
      f.mat_iso=fm

   elseif sym == :X
      if !isdefined(f, :X)
         if !isdefined(f,:ring_iso)
            fm = gen_mat_iso(nrows(f.matrix), base_ring(f.matrix))
            f.ring_iso=fm.fr
            f.mat_iso=fm
         end
         assign_from_description(f)
      end

   end

   return getfield(f, sym)

end



########################################################################
#
# Operations
#
########################################################################

function Base.:*(f::SesquilinearForm, l::FieldElem)
   l !=0 || throw(ArgumentError("Zero is not admitted"))
   parent(l)==base_ring(f) || throw(ArgumentError("The scalar does not belong to the base ring of the form"))
   if !isdefined(f,:matrix)
      return SesquilinearForm(l*f.pol, f.descr)
   else
      g = SesquilinearForm(l*f.matrix, f.descr)
      if isdefined(f,:pol) g.pol=f.pol end
      return g
   end
end

Base.:*(l::FieldElem, f::SesquilinearForm) = f*l

function Base.:^(f::SesquilinearForm{T}, x::MatElem{T}; check=false) where T <: RingElem
   @assert base_ring(f)==base_ring(x) "Incompatible base rings"
   @assert nrows(f.matrix)==nrows(x) "Incompatible dimensions"

   if check @assert rank(x)==nrows(x) "Matrix not invertible" end
   if f.descr==:hermitian
      m = x^-1*f.matrix*conjugate_transpose(x^-1)
   else
      m = x^-1*f.matrix*transpose(x^-1)
   end
   return SesquilinearForm(m, f.descr)
end

Base.:^(f::SesquilinearForm{T}, x::MatrixGroupElem{T}; check=false) where T <: RingElem = f^x.elm

function (f::SesquilinearForm{T})(v::AbstractAlgebra.Generic.FreeModuleElem{T},w::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: RingElem
   @assert f.descr!=:quadratic "Quadratic forms requires only one argument"

   if f.descr==:hermitian
      return v*f.matrix*map( y->frobenius(y,div(degree(base_ring(w)),2)),w)
   else
      return v*f.matrix*w
   end
end

function (f::SesquilinearForm{T})(v::AbstractAlgebra.Generic.FreeModuleElem{T}) where T <: RingElem
   @assert f.descr==:quadratic "Sesquilinear forms requires two arguments"
   return v*f.matrix*v
end
