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
   
   return G_hecke, homom
end
