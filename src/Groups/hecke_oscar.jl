import Hecke:
    GrpAbFinGen,
    GrpAbFinGenElem,
    GrpAbFinGenMap,
    Map

export
    GenGroupHomomorphism,
    hecke_isomorphic_group,
    isomorphism,
    oscar_isomorphic_group


const AllGrp = Union{GAPGroup,GrpAbFinGen}

struct GenGroupHomomorphism{S<:AllGrp, T<:AllGrp} <: Map{S,T, GAPMap, GenGroupHomomorphism{S,T}}
   domain::S
   codomain::T
   f                                 # homomorphisms between the two groups
end


function fhom(g; genshecke=[0,0,0], y=1)
   G=parent(g)
   w = GAP.Globals.Factorization(G.X, g.X)
   arr = GAP.gap_to_julia(GAP.Globals.LetterRepAssocWord(w))
   for i in arr
      if i<0
         y -= genshecke[-i]
      else
         y += genshecke[i]
      end
   end

   return y
end
   

function hecke_isomorphic_group(G::PcGroup; PresGen=false)     # PresGen = preserve generators, i.e. f(G[i])=H[i] for every i
   isabelian(G) || throw(ArgumentError("G is not abelian"))
   Lgap=GAP.Globals.IndependentGeneratorsOfAbelianGroup(G.X)
   ng = length(Lgap)                                    # number of generators of the Hecke group
   Gnew = GAP.Globals.Group(Lgap)
   AbInv=[GAP.Globals.Order(Lgap[i]) for i in 1:ng]                         # abelian invariants
   G_hecke = abelian_group(AbInv)
   Words = [GAP.Globals.Factorization(Gnew, gens(G)[i].X) for i in 1:ngens(G)]
   id_Gh = id(G_hecke)
   gen_of_G = [id_Gh for i in 1:ngens(G)]                  # generating set for G_hecke isomorphic to the generating set for G
   for i in 1:ngens(G)
      arr=GAP.gap_to_julia(GAP.Globals.LetterRepAssocWord(Words[i]))
      global z = id_Gh
      for j in 1:length(arr)
         if arr[j] < 0
            global z -= G_hecke[-arr[j]]
         else
            global z += G_hecke[arr[j]]
         end
      end
      gen_of_G[i] = z
   end
   if PresGen
      G_hecke = sub(G_hecke,gen_of_G)[1]
      gen_of_G = gens(G_hecke)
      id_Gh = id(G_hecke)
   end
   homom(t::GAPGroupElem) = fhom(t; genshecke=gen_of_G, y=id_Gh)
   
   return G_hecke, GenGroupHomomorphism{PcGroup,GrpAbFinGen}(G,G_hecke,homom)
end


function ghom(x, G, f; pres_gen=false)
   if pres_gen
      y = x
   else
      y = haspreimage(f,x)[2]
   end
   z = one(G)
   for i in 1:length(y.coeff)
      z*=G[i]^(y.coeff[i])
   end
   return z
end

function oscar_isomorphic_group(G::GrpAbFinGen; PresGen=false)
   G1,f = snf(G)
   L = Int64[order(y) for y in gens(G1)]
   G_oscar = abelian_group(PcGroup,L)
   if PresGen
      gen_G = [prod([G_oscar[j]^(haspreimage(f,G[i])[2].coeff[j]) for j in 1:length(L)]) for i in 1:ngens(G)]
      G_oscar = sub(G_oscar, gen_G)[1]
   end
   homom(t::GrpAbFinGenElem) = ghom(t, G_oscar, f; pres_gen=PresGen)

   return G_oscar, GenGroupHomomorphism{GrpAbFinGen,PcGroup}(G,G_oscar,homom)
end

function isomorphism(G1::GrpAbFinGen, G2::PcGroup)
   Go,f = oscar_isomorphic_group(G1)
   vero,g = isisomorphic(Go,G2)
   vero || throw(ArgumentError("The groups are not isomorphic"))

   h(t::GrpAbFinGenElem) = g(f(t))
   return GenGroupHomomorphism{GrpAbFinGen,PcGroup}(G1,G2,h)
end

function isomorphism(G1::PcGroup, G2::GrpAbFinGen)
   Gh,f = hecke_isomorphic_group(G1)
   Sh,fh = snf(Gh)
   S2,f2 = snf(G2)
   filter(x -> x != 1, Sh.snf) == filter(x -> x != 1, S2.snf) || throw(ArgumentError("The groups are not isomorphic"))
   
   h(t::GAPGroupElem{PcGroup}) = f2(inv(fh)(f(t)))
   return GenGroupHomomorphism{PcGroup,GrpAbFinGen}(G1,G2,h)
end
   

##########################################################################################

# Functions for the isomorphisms between the two groups

##########################################################################################

const AllHom = Union{GAPGroupHomomorphism,GrpAbFinGenMap,GenGroupHomomorphism}
const AllExHom = Union{GAPGroupHomomorphism,GrpAbFinGenMap}

Base.show(io::IO, f::GenGroupHomomorphism) = print(io, "Isomorphism between \n", f.domain, " and \n", f.codomain)

function (f::GenGroupHomomorphism)(x::Union{GrpAbFinGenElem,GAPGroupElem})
   parent(x)==domain(f) || throw(ArgumentError("Element not in the domain of the map"))
   return f.f(x)
end

domain(f::GenGroupHomomorphism) = f.domain
codomain(f::GenGroupHomomorphism) = f.codomain

function compose(f::GenGroupHomomorphism{T,U}, g::GenGroupHomomorphism{S,T}) where {S,T,U}
   domain(f)==codomain(g) || throw(ArgumentError("Maps not composable"))
   h(t) = f(g(t))

   return GenGroupHomomorphism{S,U}(domain(g),codomain(f),h)
end

function compose(f::GenGroupHomomorphism{S,U}, g::T) where T<:AllExHom where {S,U}
   domain(f)==codomain(g) || throw(ArgumentError("Maps not composable"))
   h(t) = f(g(t))

   return GenGroupHomomorphism{typeof(domain(g)),U}(domain(g),codomain(f),h)
end

function compose(f::AllExHom, g::GenGroupHomomorphism{S,T}) where {S,T}
   domain(f)==codomain(g) || throw(ArgumentError("Maps not composable"))
   h(t) = f(g(t))

   return GenGroupHomomorphism{S,typeof(codomain(f))}(domain(g),codomain(f),h)
end

Base.:*(f::AllHom, g::AllHom) = compose(g,f)
(f::AllHom)(g::AllHom) = compose(f,g)


# change of type

function oscar_group_homomorphism(f::GenGroupHomomorphism{S,T}) where {S,T}
   if S<:GAPGroup && T<:GAPGroup
      G=domain(f)
      H=codomain(f)
      h=hom(G,H,gens(G),[f(y) for y in gens(G)])
   else
      throw(ArgumentError("The groups are not both GAPGroups"))
   end

   return h
end

function hecke_group_homomorphism(f::GenGroupHomomorphism{S,T}) where {S,T}
   if S<:GrpAbFinGen && T<:GrpAbFinGen
      G=domain(f)
      H=codomain(f)
      h=hom(G,H,[f(y) for y in gens(G)])
   else
      throw(ArgumentError("The groups are not both Hecke groups"))
   end

   return h
end

function Base.inv(f::GenGroupHomomorphism{S,T}) where {S,T}
   G=domain(f)
   H=codomain(f)

   if G isa GAPGroup && H isa GAPGroup
      g = oscar_group_homomorphism(f)
      return GenGroupHomomorphism{T,S}(H,G,inv(g))
   elseif G isa GrpAbFinGen && H isa GrpAbFinGen
      g = hecke_group_homomorphism(f)
      return GenGroupHomomorphism{T,S}(H,G,inv(g))
   elseif G isa GAPGroup && H isa GrpAbFinGen
      h = isomorphism(H,G)
      h1 = inv(oscar_group_homomorphism(f*h))
      return GenGroupHomomorphism{T,S}(H,G,h*h1)
   else
      h = isomorphism(H,G)
      h1 = inv(hecke_group_homomorphism(f*h))
      return GenGroupHomomorphism{T,S}(H,G,h*h1)
   end
end

# FIXME: change it?
function image(f::GenGroupHomomorphism{S,T}) where {S,T}
   G = domain(f)
   H = codomain(f)

   return sub(H,[f(x) for x in gens(G)])
end

function haspreimage(f::GenGroupHomomorphism{S,T}, x) where {S,T}
   G=domain(f)
   H=codomain(f)
   
   x in H || throw(ArgumentError("Element not in the codomain of the function"))
   if S==T
      if S==GrpAbFinGen
         return haspreimage(hecke_group_homomorphism(f),x)
      elseif T<:GAPGroup
         return haspreimage(oscar_group_homomorphism(f),x)
      end
   end      
   if T==GrpAbFinGen
      x!=id(H) || return true, one(G)
      H1,e=sub(H,[f(z) for z in gens(G)])
      vero,a=haspreimage(e,x)
      if vero
         return true, prod([G[i]^(a.coeff[i]) for i in 1:ngens(G)])
      else
         return false, one(G)
      end
   end
   if T<:GAPGroup
      imH=sub(H,[f(g) for g in gens(G)])[1]
      if x in imH
         x!=one(H) || return true, id(G)
         w=GAP.Globals.Factorization(H.X,x.X)
         arr=GAP.gap_to_julia(GAP.Globals.LetterRepAssocWord(w))
         z = id(G)
         for i in arr
            if i<0
               z -= G[-i]
            else
               z += G[i]
            end
         end
         return true, z
      else
         return false, one(G)
      end
   end
end


function kernel(f::GenGroupHomomorphism{S,T}) where {S,T}
   G=domain(f)
   H=codomain(f)

   if S==T
      if S==GrpAbFinGen
         return kernel(hecke_group_homomorphism(f))
      elseif T<:GAPGroup
         return kernel(oscar_group_homomorphism(f))
      end
   end
   if T==GrpAbFinGen
      H1,f1=hecke_isomorphic_group(G)
      f2=hecke_group_homomorphism(inv(f1)*f)
      K,e=kernel(f2)
      Snf,s=snf(K)         # TODO this passage serves to build a kernel with a minimal number of generators, but it can be skipped
      gen_K = [haspreimage(f1,e(s(g)))[2] for g in gens(Snf)]
      return sub(G,gen_K)
   end
   if T<:GAPGroup
      H1,f1=oscar_isomorphic_group(G)
      f2=oscar_group_homomorphism(inv(f1)*f)
      K,e=kernel(f2)
      gen_K = [haspreimage(f1,e(g))[2] for g in gens(K)]
      return sub(G,gen_K)
   end
end

