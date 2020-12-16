import AbstractAlgebra: FieldElem

Base.getindex(V::AbstractAlgebra.Generic.FreeModule, i::Int) = gens(G, i)

# this just allows to write V([1,z,0]) instead of V([F(1),z,F(0)])
function (V::AbstractAlgebra.Generic.FreeModule)(l::Array{T,1}) where T
   if T<:typeof(base_ring(V)) return V(l)
   else return V( [base_ring(V)(j) for j in l] )
   end
end

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
