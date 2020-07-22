export
    acting_domain,
    double_coset,
    double_cosets,
    elements,
    isbicoset,
    isleft,
    isright,
    left_acting_group,
    left_coset,
    left_cosets,
    left_transversal,
    representative,
    right_acting_group,
    right_coset,
    right_cosets,
    right_transversal

# T=type of the group, S=type of the element
"""
    GroupCoset{T<: Group, S <: GAPGroupElem}
Type of group cosets. It is displayed as `H * x` (right cosets) or `x * H` (left cosets), where `H` is a subgroup of a group `G` and `x` is an element of `G`. Two cosets are equal if, and only if, they are both left (resp. right) and they contain the same elements.
"""
struct GroupCoset{T<: GAPGroup, S <: GAPGroupElem} 
   G::T                    # big group containing the subgroup and the element
   H::T                    # subgroup
   repr::S                 # element
   side::Symbol            # says if the coset is left or right
   X::GapObj               # underlying GAP object of the subgroup H
end

Base.hash(x::GroupCoset) = 0 # FIXME

function _group_coset(G::GAPGroup, H::GAPGroup, repr::GAPGroupElem, side::Symbol, X::GapObj)
  return GroupCoset{typeof(G), typeof(repr)}(G, H, repr, side, X)
end

function ==(x::GroupCoset, y::GroupCoset)
   answer = false
   if x.side == y.side && x.H == y.H
      if x.side== :left
         if (x.repr)^-1*(y.repr) in x.H
            answer = true
         end
      else
         if (y.repr)*(x.repr)^-1 in x.H
            answer = true
         end
      end
   end

   return answer
end


"""
    right_coset(H::Group, g::GAPGroupElem)
    *(H::Group, g::GAPGroupElem)
Return the coset `Hg`.
"""
function right_coset(H::GAPGroup, g::GAPGroupElem)
   @assert elem_type(H) == typeof(g)
   if !GAP.Globals.IsSubset(parent(g).X, H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return _group_coset(parent(g), H, g, :right, H.X)
end

"""
    left_coset(H::Group, g::GAPGroupElem)
    *(g::GAPGroupElem, H::Group)
Return the coset `gH`.
"""
function left_coset(H::GAPGroup, g::GAPGroupElem)
   @assert elem_type(H) == typeof(g)
   if !GAP.Globals.IsSubset(parent(g).X, H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return _group_coset(parent(g), H, g, :left, H.X)
end

function show(io::IO, x::GroupCoset)
   a = GAP.gap_to_julia(GAP.Globals.StringViewObj(x.H.X))
   b = GAP.gap_to_julia(GAP.Globals.StringViewObj(x.repr.X))
   if x.side == :right
      print(io, "Right coset   ", a, " * ", b)
   else
      print(io, "Left coset   ", b, " * ", a)
   end
   return nothing
end

"""
    isleft(c::GroupCoset)
Return whether the coset `c` is a left coset of its acting domain.
"""
isleft(c::GroupCoset) = c.side == :left

"""
    isright(c::GroupCoset)
Return whether the coset `c` is a right coset of its acting domain.
"""
isright(c::GroupCoset) = c.side == :right

Base.:*(H::GAPGroup, g::GAPGroupElem) = right_coset(H,g)
Base.:*(g::GAPGroupElem, H::GAPGroup) = left_coset(H,g)

function Base.:*(c::GroupCoset, y::GAPGroupElem)
   @assert y in c.G "element not in the group"
   if c.side == :right
      return right_coset(c.H, c.repr*y)
   else
      return left_coset(c.H^y, c.repr*y)
   end
end

function Base.:*(y::GAPGroupElem, c::GroupCoset)
   @assert y in c.G "element not in the group"
   if c.side == :left
      return left_coset(c.H, y*c.repr)
   else
      return right_coset(c.H^(y^-1), y*c.repr)
   end
end

function Base.:*(c::GroupCoset, d::GroupCoset)
   if c.side != :right || d.side != :left
      throw(ArgumentError("Wrong input"))
   end
   return double_coset(c.H, c.repr*d.repr, d.H)
end

"""
    acting_domain(C::GroupCoset)
If `C` = `Hx` or `xH`, return `H`.
"""
acting_domain(C::GroupCoset) = C.H

"""
    representative(C::GroupCoset)
If `C` = `Hx` or `xH`, return `x`.
"""
representative(C::GroupCoset) = C.repr

function elements(C::GroupCoset)
   L = GAP.Globals.AsList(C.H.X)
   l = Vector{elem_type(C.G)}(undef, length(L))
   if C.side == :left
      for i = 1:length(l)
         l[i] = group_element(C.G, C.repr.X*L[i])
      end
   else
      for i = 1:length(l)
         l[i] = group_element(C.G, L[i]*C.repr.X)
      end
   end
   return l
end

"""
    isbicoset(C::GroupCoset)
Return whether `C` is simultaneously a right coset and a left coset for the same subgroup `H`. This is the case if, and only if, the coset representative normalizes `H`.
"""
isbicoset(C::GroupCoset) = C.H^(C.repr)==C.H

"""
    right_cosets(G::Group, H::Group)
Return the array of the right cosets of `H` in `G`.
"""
function right_cosets(G::GAPGroup, H::GAPGroup)
   L = GAP.Globals.RightTransversal(G.X,H.X)
   l = Vector{GroupCoset{typeof(G), elem_type(G)}}(undef, length(L))
   for i = 1:length(l)
     l[i] = _group_coset(G, H, group_element(G, L[i]), :right, H.X)
   end
   return l
end

"""
    left_cosets(G::Group, H::Group)
Return the array of the left cosets of `H` in `G`.
"""
function left_cosets(G::GAPGroup, H::GAPGroup)
   L = GAP.Globals.RightTransversal(G.X,H.X)
   l = Vector{GroupCoset{typeof(G), elem_type(G)}}(undef, length(L))
   for i = 1:length(l)
     l[i] = _group_coset(G, H, group_element(G, L[i])^-1, :left, H.X)
   end
   return l
end

"""
    right_transversal(G::T, H::T) where T<: Group
Return a vector containing a complete set of representatives for right cosets for `H`.
"""
function right_transversal(G::T, H::T) where T<: GAPGroup
   L = GAP.Globals.RightTransversal(G.X,H.X)
   l = Vector{elem_type(G)}(undef, length(L))
   for i in 1:length(l)
      l[i] = group_element(G,L[i])
   end
   return l
end

"""
    left_transversal(G::T, H::T) where T<: Group
Return a vector containing a complete set of representatives for left cosets for `H`.
"""
function left_transversal(G::T, H::T) where T<: GAPGroup
   return [x^-1 for x in right_transversal(G,H)]
end

function Base.iterate(C::GroupCoset)
  L=GAP.Globals.Iterator(C.X)
  if GAP.Globals.IsDoneIterator(L)
    return nothing
  end
  i = GAP.Globals.NextIterator(L)
  if C.side == :left
     h = group_element(C.G, C.repr.X*i)
  else
     h = group_element(C.G, i*C.repr.X)
  end
  return h, L
end

function Base.iterate(C::GroupCoset, state)
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  if C.side == :left
     h = group_element(C.G, C.repr.X*i)
  else
     h = group_element(C.G, i*C.repr.X)
  end
  return h, state
end

order(C::GroupCoset) = order(C.H)
Base.length(C::GroupCoset) = order(C.H)

function Base.rand(C::GroupCoset)
   x=GAP.Globals.Random(C.H.X)
   if C.side == :left
      return C.repr*group_element(C.G, x)
   else
      return group_element(C.G,x)*C.repr
   end
end




#######################################################################################

# Double cosets

#######################################################################################

"""
    GroupDoubleCoset{T<: Group, S <: GAPGroupElem}
Group double coset. It is displayed as `H * x * K`, where `H` and `K` are subgroups of a group `G` and `x` is an element of `G`. Two double cosets are equal if, and only if, they contain the same elements.
"""
# T=type of the group, S=type of the element
struct GroupDoubleCoset{T <: GAPGroup, S <: GAPGroupElem}
   G::T
   H::T
   K::T
   repr::S
   X::GapObj
end

Base.hash(x::GroupDoubleCoset) = 0 # FIXME

function ==(x::GroupDoubleCoset, y::GroupDoubleCoset)
   return x.X==y.X
end

function Base.show(io::IO, x::GroupDoubleCoset)
  print(io, GAP.gap_to_julia(GAP.Globals.StringViewObj(x.H.X)),
            " * ",
            GAP.gap_to_julia(GAP.Globals.StringViewObj(x.repr.X)),
            " * ",
            GAP.gap_to_julia(GAP.Globals.StringViewObj(x.K.X)))
end


"""
    double_coset(H::Group, x::GAPGroupElem, K::Group)
    *(H::Group, x::GAPGroupElem, K::Group)
Return the double coset `HxK`.
"""
function double_coset(G::T, g::GAPGroupElem{T}, H::T) where T<: GAPGroup
   # TODO: enforce that G, H have same type
   # TODO: enforce that G, H have common overgroup
   if !GAP.Globals.IsSubset(parent(g).X,G.X)
      throw(ArgumentError("G is not a subgroup of parent(g)"))
   end
   if !GAP.Globals.IsSubset(parent(g).X,H.X)
      throw(ArgumentError("H is not a subgroup of parent(g)"))
   end
   return GroupDoubleCoset(parent(g),G,H,g,GAP.Globals.DoubleCoset(G.X,g.X,H.X))
end

Base.:*(H::GAPGroup, g::GAPGroupElem, K::GAPGroup) = double_coset(H,g,K)

"""
    double_cosets(G::T, H::T, K::T; NC=false) where T<: GAPGroup
Return the array of all the double cosets `HxK` for `x` in `G`. If `NC` = `true`, do not check whether `H` and `K` are subgroups of `G`.
"""
function double_cosets(G::T, H::T, K::T; NC=false) where T<: GAPGroup
   if NC
      dcs = GAP.Globals.DoubleCosetsNC(G.X,H.X,K.X)
   else
      @assert issubgroup(G,H)[1] "H is not a subgroup of G"
      @assert issubgroup(G,K)[1] "K is not a subgroup of G"
      dcs = GAP.Globals.DoubleCosets(G.X,H.X,K.X)
   end
   res = Vector{GroupDoubleCoset{T,GAPGroupElem{T}}}(undef, length(dcs))
   for i = 1:length(res)
     dc = dcs[i]
     g = group_element(G, GAP.Globals.Representative(dc))
     res[i] = GroupDoubleCoset(G,H,K,g,dc)
   end
   return res
   #return [GroupDoubleCoset(G,H,K,group_element(G.X,GAP.Globals.Representative(dc)),dc) for dc in dcs]
end

"""
    elements(C::GroupDoubleCoset)
Return the array of all elements of the double coset `C`.
"""
function elements(C::GroupDoubleCoset)
  L = GAP.Globals.AsList(C.X)
  l = Vector{elem_type(C.G)}(undef, length(L))
  for i = 1:length(l)
    l[i] = group_element(C.G, L[i])
  end
  return l
end

"""
    order(C::Union{GroupCoset,GroupDoubleCoset})
Return the cardinality of the (double) coset `C`.
"""
order(C::GroupDoubleCoset) = GAP.Globals.Size(C.X)
Base.length(C::GroupDoubleCoset) = GAP.Globals.Size(C.X)

"""
    rand(C::Union{GroupCoset,GroupDoubleCoset})
Return a random element of the (double) coset `C`.
"""
Base.rand(C::GroupDoubleCoset) = group_element(C.G, GAP.Globals.Random(C.X))

"""
    representative(C::GroupDoubleCoset)
If `C` = `HxK`, return `x`.
"""
representative(C::GroupDoubleCoset) = C.repr

"""
    left_acting_group(C::GroupDoubleCoset)
If `C` = `HxK`, return `H`
"""
left_acting_group(C::GroupDoubleCoset) = C.H

"""
    left_acting_group(C::GroupDoubleCoset)
If `C` = `HxK`, return `K`
"""
right_acting_group(C::GroupDoubleCoset) = C.K

function Base.iterate(G::GroupDoubleCoset)
  L=GAP.Globals.Iterator(G.X)
  if GAP.Globals.IsDoneIterator(L)
    return nothing
  end
  i = GAP.Globals.NextIterator(L)
  return group_element(G.G, i), L
end

function Base.iterate(G::GroupDoubleCoset, state)
  if GAP.Globals.IsDoneIterator(state)
    return nothing
  end
  i = GAP.Globals.NextIterator(state)
  return group_element(G.G, i), state
end

"""
    intersect(V::AbstractVector{Union{T, GroupCoset, GroupDoubleCoset}}) where T <: GAPGroup
Return an array containing all elements belonging to all groups and cosets in the vector `V`.
"""
function intersect(V::AbstractVector{Union{T, GroupCoset, GroupDoubleCoset}}) where T <: GAPGroup
   if typeof(V[1]) <: GAPGroup
      G = V[1]
   else
      G = V[1].G
   end
   l = GAP.julia_to_gap([v.X for v in V])
   ints = GAP.Globals.Intersection(l)
   L = Vector{typeof(G)}(undef, length(ints))
   for i in 1:length(ints)
      L[i] = group_element(G,ints[i])
   end

   return L
end
