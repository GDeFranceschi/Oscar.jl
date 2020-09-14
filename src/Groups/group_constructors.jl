################################################################################
#
#  Some basic constructors
#  
################################################################################

export
    abelian_group,
    alternating_group,
    cyclic_group,
    dihedral_group,
    free_abelian_group,
    free_group,
    general_linear_group,
    isabelian,
    isalternating_group,
    iscyclic,
    isdihedral_group,
    isquaternion_group,
    issymmetric_group,
    mathieu_group,
    omega_group,
    orthogonal_group,
    quaternion_group,
    special_linear_group,
    special_orthogonal_group,
    special_unitary_group,
    symmetric_group,
    symplectic_group,
    unitary_group,
    GL, GO, GU, SL, SO, Sp, SU


_gap_filter(::Type{PermGroup}) = GAP.Globals.IsPermGroup
_gap_filter(::Type{PcGroup}) = GAP.Globals.IsPcGroup
_gap_filter(::Type{FPGroup}) = GAP.Globals.IsFpGroup
_gap_filter(::Type{MatrixGroup}) = GAP.Globals.IsMatrixGroup

# TODO: matrix group handling usually is more complex: there usually
# is another extra argument then to specify the base field
#_gap_filter(::Type{MatrixGroup}) = GAP.Globals.IsMatrixGroup

"""
    symmetric_group(n::Int64)
    symmetric_group(::Type{T}, n::Int)
Return the full symmetric group over a set of `n` elements. The group is returned of type `T` for `T` in {`PermGroup`, `PcGroup`}. If `T` is not specified, then `T` is set as `PermGroup`.
"""
symmetric_group(n::Int64) = symmetric_group(PermGroup, n)

function symmetric_group(::Type{T}, n::Int) where T <: GAPGroup
  if n < 1
    throw(ArgumentError("n must be a positive integer"))
  end
  return T(GAP.Globals.SymmetricGroup(_gap_filter(T), n))
end

function issymmetric_group(G::GAPGroup)
  return GAP.Globals.IsSymmetricGroup(G.X)
end

"""
    alternating_group(n::Int64)
    alternating_group(::Type{T}, n::Int)
Return the full alternating group over a set of `n` elements. The group is returned of type `T` for `T` in {`PermGroup`, `PcGroup`}. If `T` is not specified, then `T` is set as `PermGroup`.
"""
alternating_group(n::Int64) = alternating_group(PermGroup, n)

function alternating_group(::Type{T}, n::Int) where T <: GAPGroup
  if n < 1
    throw(ArgumentError("n must be a positive integer"))
  end
  return T(GAP.Globals.AlternatingGroup(_gap_filter(T), n))
end

function isalternating_group(G::GAPGroup)
  return GAP.Globals.IsAlternatingGroup(G.X)
end

cyclic_group(n::Int) = cyclic_group(PcGroup, n)

"""
    cyclic_group(::Type{T}, n::Int)
Return the cyclic group of order `n` and type `T`. If the type is not specified, the group is returned of type `PcGroup`.
"""
function cyclic_group(::Type{T}, n::Int) where T <: GAPGroup
  return T(GAP.Globals.CyclicGroup(_gap_filter(T), n))
end

function iscyclic(G::GAPGroup)
  return GAP.Globals.IsCyclic(G.X)
end

# already defined in Hecke
#=
function abelian_group(v::Vector{Int})
  for i = 1:length(v)
    iszero(v[i]) && error("Cannot represent an infinite group as a polycyclic group")
  end
  v1 = GAP.julia_to_gap(v)
  return PcGroup(GAP.Globals.AbelianGroup(v1))
end
=#

"""
    abelian_group(::Type{T}, v::Vector{Int}) where T <: Group -> PcGroup
Return the direct product of cyclic groups of order v[1] x v[2] x ... x v[n], as group of type `T`. Here, `T` must be of type `PermGroup`, `FPGroup` or `PcGroup`.
"""
function abelian_group(::Type{T}, v::Vector{Int}) where T <: GAPGroup
  v1 = GAP.julia_to_gap(v)
  return T(GAP.Globals.AbelianGroup(_gap_filter(T), v1))
end

"""
    isabelian(G::Group)
Return whether `G` is abelian.
"""
function isabelian(G::GAPGroup)
  return GAP.Globals.IsAbelian(G.X)
end

function mathieu_group(n::Int)
  @assert n in Int[9, 10, 11, 12, 21, 22, 23, 24]
  return PermGroup(GAP.Globals.MathieuGroup(n), n)
end


################################################################################
#
# begin FpGroups
#
################################################################################

"""
    free_group
There are four ways to define a free group.
- `free_group(n::Int) -> FPGroup`; return the free group of rank `n`, with generators printed as `"f1"`,`"f2"`,`"f3"`, etc.
- `free_group(L::String...) -> FPGroup`; return the free group with length(`L`) generators, printed as `L[1]`, `L[2]`, `L[3]`, etc.
- `free_group(L::Array{String,1}) -> FPGroup`; same as above.
- `free_group(n::Int, s::String) -> FPGroup`; return the free group of rank `n`, with generators printed as `"s1"`, `"s2"`, `"s3"`, etc.

!!! warning "Note"
    In every case, it is *not* defined a variable named as the generators are printed.
"""
function free_group(n::Int)
   return FPGroup(GAP.Globals.FreeGroup(n))
end


function free_group(L::Array{String,1})
   J=GAP.julia_to_gap([GAP.julia_to_gap(x) for x in L])
   return FPGroup(GAP.Globals.FreeGroup(J))
end

function free_group(L::String...)
   J=GAP.julia_to_gap([GAP.julia_to_gap(x) for x in L])
   return FPGroup(GAP.Globals.FreeGroup(J))
end

function free_group(n::Int, s::String)
   return FPGroup(GAP.Globals.FreeGroup(n,GAP.julia_to_gap(s)))
end

# FIXME: a function `free_abelian_group` with the same signature is
# already being defined by Hecke
#function free_abelian_group(n::Int)
#  return FPGroup(GAP.Globals.FreeAbelianGroup(n))
#end

function free_abelian_group(::Type{FPGroup}, n::Int)
 return FPGroup(GAP.Globals.FreeAbelianGroup(n))
end


# for the definition of group modulo relations, see the quo function in the sub.jl section

function free_group(G::FPGroup)
   return FPGroup(GAP.Globals.FreeGroupOfFpGroup(G.X))
end

################################################################################
#
# end FpGroups
#
################################################################################

"""
    dihedral_group(n::Int)
    dihedral_group(::Type{T}, n::Int)
Return the dihedral group of order `n` of type `T`, where `T` is in {`PcGroup`,`PermGroup`,`FPGroup`}. In the first case, the type is set as `PcGroup`.
"""
dihedral_group(n::Int) = dihedral_group(PcGroup, n)

function dihedral_group(::Type{T}, n::Int) where T <: GAPGroup
  @assert iseven(n)
  return T(GAP.Globals.DihedralGroup(_gap_filter(T), n))
end

"""
    quaternion_group(n::Int)
    quaternion_group(::Type{T}, n::Int)
Return the quaternion group of order `n` of type `T`, where `T` is in {`PcGroup`,`PermGroup`,`FPGroup`,`MatrixGroup`}. In the first case, the type is set as `PcGroup`.
"""
quaternion_group(n::Int) = quaternion_group(PcGroup, n)

function quaternion_group(::Type{T}, n::Int) where T <: GAPGroup 
   @assert iszero(mod(n, 4))
  return T(GAP.Globals.QuaternionGroup(_gap_filter(T), n))
end

function isquaternion_group(G::GAPGroup)
  return GAP.Globals.IsQuaternionGroup(G.X)
end

################################################################################
#
# start isometry groups
#
################################################################################

"""
    general_linear_group(n::Int, q::Int)
    general_linear_group(T::Type, n::Int, q::Int)
    GL = general_linear_group
return the general linear group of dimension `n` over the field GF(`q`). It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
function general_linear_group(::Type{T}, n::Int, F::FqNmodFiniteField) where T <: GAPGroup
   q=Int(order(F))
   if T==MatrixGroup
      return T(GAP.Globals.GL(_gap_filter(T),n,q),F)
   else
      return T(GAP.Globals.GL(_gap_filter(T),n,q))
   end
end

general_linear_group(n::Int, F::FqNmodFiniteField) = general_linear_group(MatrixGroup, n, F)

function general_linear_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup
   a,p=ispower(q)
   isprime(p) || throw(ArgumentError("q must be a prime power"))
   F=GF(p,a)[1]
   return general_linear_group(T,n,F)
end

general_linear_group(n::Int, q::Int) = general_linear_group(MatrixGroup, n, q)

"""
    special_linear_group(n::Int, q::Int)
    special_linear_group(T::Type, n::Int, q::Int)
    SL = special_linear_group
return the special linear group of dimension `n` over the field GF(`q`). It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
function special_linear_group(::Type{T}, n::Int, F::FqNmodFiniteField) where T <: GAPGroup
   q=Int(order(F))
   if T==MatrixGroup
      return T(GAP.Globals.SL(_gap_filter(T),n,q),F)
   else
      return T(GAP.Globals.SL(_gap_filter(T),n,q))
   end
end

special_linear_group(n::Int, F::FqNmodFiniteField) = special_linear_group(MatrixGroup, n, F)

function special_linear_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup
   a,p=ispower(q)
   isprime(p) || throw(ArgumentError("q must be a prime power"))
   F=GF(p,a)[1]
   return special_linear_group(T,n,F)
end

special_linear_group(n::Int, q::Int) = special_linear_group(MatrixGroup, n, q)

"""
    symplectic_group(n::Int, q::Int)
    symplectic_group(T::Type, n::Int, q::Int)
    Sp = symplectic_group
return the special linear group of dimension `n` over the field GF(`q`), for `n` even. It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
function symplectic_group(::Type{T}, n::Int, F::FqNmodFiniteField) where T <: GAPGroup
   q=Int(order(F))
   iseven(n) || throw(ArgumentError("dimension must be even"))
   if T==MatrixGroup
      return T(GAP.Globals.Sp(_gap_filter(T),n,q),F)
   else
      return T(GAP.Globals.Sp(_gap_filter(T),n,q))
   end
end

symplectic_group(n::Int, F::FqNmodFiniteField) = symplectic_group(MatrixGroup, n, F)

function symplectic_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup
   a,p=ispower(q)
   isprime(p) || throw(ArgumentError("q must be a prime power"))
   F=GF(p,a)[1]
   return symplectic_group(T,n,F)
end

symplectic_group(n::Int, q::Int) = symplectic_group(MatrixGroup, n, q)

"""
    unitary_group(n::Int, q::Int)
    unitary_group(T::Type, n::Int, q::Int)
    GU = unitary_group
return the unitary group of dimension `n` over the field GF(`q^2`). It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
function unitary_group(::Type{T}, n::Int, F::FqNmodFiniteField) where T <: GAPGroup
   q=Int(order(F))
   issquare(fmpz(q)) || throw(ArgumentError("The field must have even degree"))
   q=Int(sqrt(q))
   if T==MatrixGroup
      return T(GAP.Globals.GU(_gap_filter(T),n,q),F)
   else
      return T(GAP.Globals.GU(_gap_filter(T),n,q))
   end
end

unitary_group(n::Int, F::FqNmodFiniteField) = unitary_group(MatrixGroup, n, F)

function unitary_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup
   a,p=ispower(q)
   isprime(p) || throw(ArgumentError("q must be a prime power"))
   F=GF(p,2*a)[1]
   return unitary_group(T,n,F)
end

unitary_group(n::Int, q::Int) = unitary_group(MatrixGroup, n, q)


"""
    special_unitary_group(n::Int, q::Int)
    special_unitary_group(T::Type, n::Int, q::Int)
    SU = special_unitary_group
return the special unitary group of dimension `n` over the field GF(`q^2`). It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
function special_unitary_group(::Type{T}, n::Int, F::FqNmodFiniteField) where T <: GAPGroup
   q=Int(order(F))
   issquare(fmpz(q)) || throw(ArgumentError("The field must have even degree"))
   q=Int(sqrt(q))
   if T==MatrixGroup
      return T(GAP.Globals.SU(_gap_filter(T),n,q),F)
   else
      return T(GAP.Globals.SU(_gap_filter(T),n,q))
   end
end

special_unitary_group(n::Int, F::FqNmodFiniteField) = special_unitary_group(MatrixGroup, n, F)

function special_unitary_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup
   a,p=ispower(q)
   isprime(p) || throw(ArgumentError("q must be a prime power"))
   F=GF(p,2*a)[1]
   return special_unitary_group(T,n,F)
end

special_unitary_group(n::Int, q::Int) = special_unitary_group(MatrixGroup, n, q)


"""
    orthogonal_group(n::Int, q::Int)
    orthogonal_group(T::Type, n::Int, q::Int)
    orthogonal_group(e::Int, n::Int, q::Int)
    orthogonal_group(::Type{T}, e::Int, n::Int, q::Int)
    GO = orthogonal_group
return the orthogonal group of dimension `n` over the field GF(`q`) of type `e`, where `e` in {`+1`,`-1`} for `n` even and `e`=`0` for `n` odd. If `n` is odd, `e` can be omitted. It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
function orthogonal_group(::Type{T}, e::Int, n::Int, F::FqNmodFiniteField) where T <: GAPGroup
   q=Int(order(F))
   if isodd(n)
      e==0 || throw(ArgumentError("Error, sign <e> must be 0 in odd dimension"))
   else
      e in (-1,+1) || throw(ArgumentError("Error, sign <e> must be -1, +1 in even dimension"))
   end
   if T==MatrixGroup
      return T(GAP.Globals.GO(_gap_filter(T),e,n,q),F)
   else
      return T(GAP.Globals.GO(_gap_filter(T),e,n,q))
   end
end

orthogonal_group(n::Int, F::FqNmodFiniteField) = orthogonal_group(MatrixGroup, 0, n, F)

orthogonal_group(::Type{T}, n::Int, F::FqNmodFiniteField) where T <: GAPGroup = orthogonal_group(T, 0, n, F)

orthogonal_group(e::Int, n::Int, F::FqNmodFiniteField) = orthogonal_group(MatrixGroup, e, n, F)

function orthogonal_group(::Type{T}, e::Int, n::Int, q::Int) where T <: GAPGroup
   a,p=ispower(q)
   isprime(p) || throw(ArgumentError("q must be a prime power"))
   F=GF(p,a)[1]
   return orthogonal_group(T,e,n,F)
end

orthogonal_group(n::Int, q::Int) = orthogonal_group(MatrixGroup, 0, n, q)

orthogonal_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup = orthogonal_group(T, 0, n, q)

orthogonal_group(e::Int, n::Int, q::Int) = orthogonal_group(MatrixGroup, e, n, q)

"""
    special_orthogonal_group(n::Int, q::Int)
    special_orthogonal_group(T::Type, n::Int, q::Int)
    special_orthogonal_group(e::Int, n::Int, q::Int)
    special_orthogonal_group(::Type{T}, e::Int, n::Int, q::Int)
    SO = special_orthogonal_group
return the special orthogonal group of dimension `n` over the field GF(`q`) of type `e`, where `e` in {`+1`,`-1`} for `n` even and `e`=`0` for `n` odd. If `n` is odd, `e` can be omitted. It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
function special_orthogonal_group(::Type{T}, e::Int, n::Int, F::FqNmodFiniteField) where T <: GAPGroup
   q=Int(order(F))
   if isodd(n)
      e==0 || throw(ArgumentError("Error, sign <e> must be 0 in odd dimension"))
   else
      e in (-1,+1) || throw(ArgumentError("Error, sign <e> must be -1, +1 in even dimension"))
   end
   if T==MatrixGroup
      return T(GAP.Globals.SO(_gap_filter(T),e,n,q),F)
   else
      return T(GAP.Globals.SO(_gap_filter(T),e,n,q))
   end
end

special_orthogonal_group(n::Int, F::FqNmodFiniteField) = special_orthogonal_group(MatrixGroup, 0, n, F)

special_orthogonal_group(::Type{T}, n::Int, F::FqNmodFiniteField) where T <: GAPGroup = special_orthogonal_group(T, 0, n, F)

special_orthogonal_group(e::Int, n::Int, F::FqNmodFiniteField) = special_orthogonal_group(MatrixGroup, e, n, F)

function special_orthogonal_group(::Type{T}, e::Int, n::Int, q::Int) where T <: GAPGroup
   a,p=ispower(q)
   isprime(p) || throw(ArgumentError("q must be a prime power"))
   F=GF(p,a)[1]
   return special_orthogonal_group(T,e,n,F)
end

special_orthogonal_group(n::Int, q::Int) = special_orthogonal_group(MatrixGroup, 0, n, q)

special_orthogonal_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup = special_orthogonal_group(T, 0, n, q)

special_orthogonal_group(e::Int, n::Int, q::Int) = special_orthogonal_group(MatrixGroup, e, n, q)

"""
    omega_group(n::Int, q::Int)
    omega_group(T::Type, n::Int, q::Int)
    omega_group(e::Int, n::Int, q::Int)
    omega_group(::Type{T}, e::Int, n::Int, q::Int)
return the Omega group of dimension `n` over the field GF(`q`) of type `e`, where `e` in {`+1`,`-1`} for `n` even and `e`=`0` for `n` odd. If `n` is odd, `e` can be omitted. It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
function omega_group(::Type{T}, e::Int, n::Int, F::FqNmodFiniteField) where T <: GAPGroup
   q=Int(order(F))
   if isodd(n)
      e==0 || throw(ArgumentError("Error, sign <e> must be 0 in odd dimension"))
   else
      e in (-1,+1) || throw(ArgumentError("Error, sign <e> must be -1, +1 in even dimension"))
   end
   if T==MatrixGroup
      return T(GAP.Globals.Omega(_gap_filter(T),e,n,q),F)
   else
      return T(GAP.Globals.Omega(_gap_filter(T),e,n,q))
   end
end

omega_group(n::Int, F::FqNmodFiniteField) = omega_group(MatrixGroup, 0, n, F)

omega_group(::Type{T}, n::Int, F::FqNmodFiniteField) where T <: GAPGroup = omega_group(T, 0, n, F)

omega_group(e::Int, n::Int, F::FqNmodFiniteField) = omega_group(MatrixGroup, e, n, F)

function omega_group(::Type{T}, e::Int, n::Int, q::Int) where T <: GAPGroup
   a,p=ispower(q)
   isprime(p) || throw(ArgumentError("q must be a prime power"))
   F=GF(p,a)[1]
   return omega_group(T,e,n,F)
end

omega_group(n::Int, q::Int) = omega_group(MatrixGroup, 0, n, q)

omega_group(::Type{T}, n::Int, q::Int) where T <: GAPGroup = omega_group(T, 0, n, q)

omega_group(e::Int, n::Int, q::Int) = omega_group(MatrixGroup, e, n, q)


"""
    general_linear_group(n::Int, q::Int)
    general_linear_group(T::Type, n::Int, q::Int)
    GL(n::Int, q::Int)
return the general linear group of dimension `n` over the field GF(`q`). It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
const GL = general_linear_group

"""
    special_linear_group(n::Int, q::Int)
    special_linear_group(T::Type, n::Int, q::Int)
    SL(n::Int, q::Int)
return the special linear group of dimension `n` over the field GF(`q`). It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
const SL = special_linear_group

"""
    symplectic_group(n::Int, q::Int)
    symplectic_group(T::Type, n::Int, q::Int)
    Sp(n::Int, q::Int)
return the special linear group of dimension `n` over the field GF(`q`), for `n` even. It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
const Sp = symplectic_group

"""
    unitary_group(n::Int, q::Int)
    unitary_group(T::Type, n::Int, q::Int)
    GU(n::Int, q::Int)
return the unitary group of dimension `n` over the field GF(`q^2`). It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
const GU = unitary_group

"""
    special_unitary_group(n::Int, q::Int)
    special_unitary_group(T::Type, n::Int, q::Int)
    SU(n::Int, q::Int)
return the special unitary group of dimension `n` over the field GF(`q^2`). It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
const SU = special_unitary_group

"""
    orthogonal_group(n::Int, q::Int)
    orthogonal_group(T::Type, n::Int, q::Int)
    orthogonal_group(e::Int, n::Int, q::Int)
    orthogonal_group(::Type{T}, e::Int, n::Int, q::Int)
    GO = orthogonal_group
return the orthogonal group of dimension `n` over the field GF(`q`) of type `e`, where `e` in {`+1`,`-1`} for `n` even and `e`=`0` for `n` odd. If `n` is odd, `e` can be omitted. It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
const GO = orthogonal_group

"""
    special_orthogonal_group(n::Int, q::Int)
    special_orthogonal_group(T::Type, n::Int, q::Int)
    special_orthogonal_group(e::Int, n::Int, q::Int)
    special_orthogonal_group(::Type{T}, e::Int, n::Int, q::Int)
    SO = special_orthogonal_group
return the special orthogonal group of dimension `n` over the field GF(`q`) of type `e`, where `e` in {`+1`,`-1`} for `n` even and `e`=`0` for `n` odd. If `n` is odd, `e` can be omitted. It is returned of type `T` for `T` = `MatrixGroup` or `PermGroup`.
"""
const SO = special_orthogonal_group


################################################################################
#
# end isometry groups
#
################################################################################
