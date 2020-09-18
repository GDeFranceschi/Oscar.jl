function hecke_isomorphic_group(G::PcGroup)
   isabelian(G) || throw(ArgumentError("G is not abelian"))
   Lgap=GAP.Globals.IndependentGeneratorsOfAbelianGroup(G)
   L=[group_elem(G,Lgap[i]) for i in 1:length(Lgap)]
   AbInv=[order(y) for y in L]                         # abelian invariants
   

# given an abelian group G with gens(G)=[x1,x2,x3,...] and g in G, return [m1,m2,m3,...] where g=prod(xi^mi)
function factorization_given_gen(G::GAPGroup)
   hom=GAP.Globals.EpimorphismFromFreeGroup(G.X)
   Lgap=[GAP.Globals.PreImagesRepresentative(hom,y.X) for y in gens(G)]
   # List(GeneratorsOfGroup(G),x->PreImagesRepresentative(hom,x));

