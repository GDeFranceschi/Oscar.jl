
@testset "Oscar vs Hecke" begin
   Gh=abelian_group(4,4,2)
   Go=abelian_group(PcGroup,[4,2,4])

   G,g=hecke_isomorphic_group(Go)
   @test G isa GrpAbFinGen
   @test g isa GenGroupHomomorphism
   @test domain(g)==Go
   @test codomain(g)==G
   for i=1:3
      @test g(Go[i]) in G
      @test order(Go[i])==order(g(Go[i]))
      for j=1:3
         @test g(Go[i]*Go[j])==g(Go[i])+g(Go[j])
      end
   end
   H,h=oscar_isomorphic_group(Gh)
   @test H isa PcGroup
   @test h isa GenGroupHomomorphism
   @test domain(h)==Gh
   @test codomain(h)==H
   for i=1:3
      @test h(Gh[i]) in H
      @test order(Gh[i])==order(h(Gh[i]))
      for j=1:3
         @test h(Gh[i]+Gh[j])==h(Gh[i])*h(Gh[j])
      end
   end

   f=isomorphism(Gh,Go)
   @test f isa GenGroupHomomorphism
   @test domain(f)==Gh
   @test codomain(f)==Go
   for i=1:3
      @test f(Gh[i]) in Go
      @test order(Gh[i])==order(f(Gh[i]))
      for j=1:3
         @test f(Gh[i]+Gh[j])==f(Gh[i])*f(Gh[j])
      end
   end

end
