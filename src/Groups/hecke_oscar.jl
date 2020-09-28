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
   f                                 # homomorphisms between the two groups
   domain::S
   codomain::T
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
   

function hecke_isomorphic_group(G::PcGroup)
   isabelian(G) || throw(ArgumentError("G is not abelian"))
   Lgap=GAP.Globals.IndependentGeneratorsOfAbelianGroup(G.X)
   ng = length(Lgap)                                    # number of generators of the Hecke group
   Gnew = GAP.Globals.Group(Lgap)
   AbInv=[GAP.Globals.Order(Lgap[i]) for i in 1:ng]                         # abelian invariants
   G_hecke = abelian_group(AbInv)
   Words = [GAP.Globals.Factorization(Gnew, gens(G)[i].X) for i in 1:ngens(G)]
   id_Gh = G_hecke[1]*AbInv[1]                # TODO: how do I have access to the identity of G_hecke ?
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
   homom(t::GAPGroupElem) = fhom(t; genshecke=gen_of_G, y=id_Gh) 
   
   return G_hecke, GenGroupHomomorphism{PcGroup,GrpAbFinGen}(homom,G,G_hecke)
end


function ghom(x, G, f)
   y = haspreimage(f,x)[2]
   z = one(G)
   for i in 1:length(y.coeff)
      z*=G[i]^(y.coeff[i])
   end
   return z
end

function oscar_isomorphic_group(G::GrpAbFinGen)
   G1,f = snf(G)
   L = Int64[order(y) for y in gens(G1)]
   G_oscar = abelian_group(PcGroup,L)
   homom(t::GrpAbFinGenElem) = ghom(t, G_oscar, f)

   return G_oscar, GenGroupHomomorphism{GrpAbFinGen,PcGroup}(homom,G,G_oscar)
end

function isomorphism(G1::GrpAbFinGen, G2::PcGroup)
   Go,f = oscar_isomorphic_group(G1)
   vero,g = isisomorphic(Go,G2)
   vero || throw(ArgumentError("The groups are not isomorphic"))

   h(t::GrpAbFinGenElem) = g(f(t))
   return GenGroupHomomorphism{GrpAbFinGen,PcGroup}(h, G1,G2)
end

function isomorphism(G1::PcGroup, G2::GrpAbFinGen)
   Gh,f = hecke_isomorphic_group(G1)
   Sh,fh = snf(Gh)
   S2,f2 = snf(G2)
   filter(x -> x != 1, Sh.snf) == filter(x -> x != 1, S2.snf) || throw(ArgumentError("The groups are not isomorphic"))
   
   h(t::GAPGroupElem{PcGroup}) = f2(inv(fh)(f(t)))
   return GenGroupHomomorphism{PcGroup,GrpAbFinGen}(h,G1,G2)
end
   

##########################################################################################

# Functions for the isomorphisms between the two groups

##########################################################################################

const AllHom = Union{GAPGroupHomomorphism,GrpAbFinGenMap,GenGroupHomomorphism}

Base.show(io::IO, f::GenGroupHomomorphism) = print(io, "Isomorphism between \n", f.domain, " and \n", f.codomain)

(f::GenGroupHomomorphism)(x::Union{GrpAbFinGenElem,GAPGroupElem}) = f.f(x)

domain(f::GenGroupHomomorphism) = f.domain
codomain(f::GenGroupHomomorphism) = f.codomain

function compose(f::GenGroupHomomorphism{S,T}, g::GenGroupHomomorphism{T,U}) where {S,T,U}
   domain(g)==codomain(f) || throw(ArgumentError("Maps not composable"))
   h(t) = g(f(t))

   return GenGroupHomomorphism{S,U}(h,domain(f),codomain(g))
end

function compose(f::GenGroupHomomorphism{S,U}, g::T) where T<:AllHom where {S,U}
   domain(g)==codomain(f) || throw(ArgumentError("Maps not composable"))
   h(t) = g(f(t))

   return GenGroupHomomorphism{S,typeof(codomain(g))}(h,domain(f),codomain(g))
end

function compose(f::AllHom, g::GenGroupHomomorphism{S,T}) where {S,T}
   domain(g)==codomain(f) || throw(ArgumentError("Maps not composable"))
   h(t) = g(f(t))

   return GenGroupHomomorphism{typeof(domain(f)),T}(h,domain(f),codomain(g))
end

Base.:*(f::AllHom, g::AllHom) = compose(f,g)
(f::AllHom)(g::AllHom) = compose(g,f)
