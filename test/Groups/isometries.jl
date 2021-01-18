@testset "Definition forms" begin
   T,t = PolynomialRing(FiniteField(3),"t")
   F,z = FiniteField(t^2+1,"z")

   B = matrix(F,4,4,[0 1 0 0; 2 0 0 0; 0 0 0 z+2; 0 0 1-z 0])
   @test Oscar._is_skew_symmetric(B)
   f = alternating_form(B)
   @test f isa SesquilinearForm
   @test gram_matrix(f)==B
   @test f==alternating_form(gram_matrix(f))
   @test base_ring(B)==F
   @test !Oscar._is_symmetric(B)
   @test !Oscar._is_hermitian(B)
   @test isalternating_form(f)
   @test !isquadratic_form(f)
   @test !issymmetric_form(f)
   @test !ishermitian_form(f)
   @test_throws AssertionError f = symmetric_form(B)
   @test_throws AssertionError f = hermitian_form(B)
   @test hermitian_form(B; check=false) isa SesquilinearForm

   B = matrix(F,4,4,[0 1 0 0; 1 0 0 0; 0 0 0 z+2; 0 0 -1-z 0])   
   @test Oscar._is_hermitian(B)
   f = hermitian_form(B)
   @test f isa SesquilinearForm
   @test gram_matrix(f)==B
   @test ishermitian_form(f)
   @test_throws AssertionError f = symmetric_form(B)
   @test_throws AssertionError f = alternating_form(B)
   @test_throws ArgumentError corresponding_quadratic_form(f)

   B = matrix(F,4,4,[0 1 0 0; 1 0 0 0; 0 0 0 z+2; 0 0 z+2 0])
   @test Oscar._is_symmetric(B)
   f = symmetric_form(B)
   @test f isa SesquilinearForm
   @test gram_matrix(f)==B
   @test issymmetric_form(f)
   @test !ishermitian_form(f)
   @test_throws AssertionError f = alternating_form(B)
   Qf = corresponding_quadratic_form(f)
   R = PolynomialRing(F,4)[1]
   p = R[1]*R[2]+(z+2)*R[3]*R[4]
   Q = quadratic_form(p)
   @test Q==Qf
   @test corresponding_quadratic_form(corresponding_bilinear_form(Q))==Q
   @test corresponding_bilinear_form(Q)==f
   B1 = matrix(F,4,4,[0 0 0 0; 1 0 0 0; 0 0 0 0; 0 0 z+2 0])
   Q1 = quadratic_form(B1)
   @test Q1==Q
   @test gram_matrix(Q1)!=B1
   @test defining_polynomial(Q1) isa AbstractAlgebra.Generic.MPoly
   pf = defining_polynomial(Q1)
   @test defining_polynomial(Q1)==parent(pf)[1]*parent(pf)[2]+(z+2)*parent(pf)[3]*parent(pf)[4]
# I can't put pf==p, because it returns FALSE. The line 
 #      PolynomialRing(F,4)[1]==PolynomialRing(F,4)[1]
# returns FALSE.
   @test_throws ArgumentError corresponding_quadratic_form(Q)
   @test_throws ArgumentError corresponding_bilinear_form(f)

   T,t = PolynomialRing(FiniteField(2),"t")
   F,z = FiniteField(t^2+t+1,"z")
   R = PolynomialRing(F,4)[1]
   p = R[1]*R[2]+z*R[3]*R[4]
   Q = quadratic_form(p)
   @test isquadratic_form(Q)
   @test gram_matrix(Q)==matrix(F,4,4,[0 1 0 0; 0 0 0 0; 0 0 0 z; 0 0 0 0])
   f = corresponding_bilinear_form(Q)
   @test isalternating_form(f)
   @test gram_matrix(f)==matrix(F,4,4,[0 1 0 0; 1 0 0 0; 0 0 0 z; 0 0 z 0])
   @test_throws ArgumentError corresponding_quadratic_form(f)


end

@testset "Transform form" begin
   @testset "Unitary case" begin
      for q in [2,3]
      for n in [5,6]
         G = GL(n,q^2)
         for i in 1:5
            x = rand(G).elm
            x = x+conjugate_transpose(x)
            while rank(x) != n
               x = rand(G).elm
               x = x+conjugate_transpose(x)
            end
            y = rand(G).elm
            y = y+conjugate_transpose(y)
            while rank(y) != n
               y = rand(G).elm
               y = y+conjugate_transpose(y)
            end
            x = hermitian_form(x)
            y = hermitian_form(y)
            vero,C = iscongruent(x,y)
            @test vero
            @test C*x.matrix*conjugate_transpose(C)==y.matrix
         end
      end
      end
   end

   @testset "Symmetric case" begin
      for q in [3,9]
      for n in [5,6,8]
         G = GL(n,q^2)
         for i in 1:5
            x = rand(G).elm
            x = x+transpose(x)
            while rank(x) != n
               x = rand(G).elm
               x = x+transpose(x)
            end
            y = rand(G).elm
            y = y+transpose(y)
            while rank(y) != n
               y = rand(G).elm
               y = y+transpose(y)
            end
            x = symmetric_form(x)
            y = symmetric_form(y)
            vero,C = iscongruent(x,y)
            if vero @test C*x.matrix*transpose(C)==y.matrix end
         end
      end
      end
   end

   @testset "Alternating case" begin
      for q in [2,3,4,9]
      for n in [6,8]
         G = GL(n,q^2)
         for i in 1:5
            x = rand(G).elm
            x = x-transpose(x)
            while rank(x) != n
               x = rand(G).elm
               x = x-transpose(x)
            end
            y = rand(G).elm
            y = y-transpose(y)
            while rank(y) != n
               y = rand(G).elm
               y = y-transpose(y)
            end
            x = alternating_form(x)
            y = alternating_form(y)
            vero,C = iscongruent(x,y)
            @test vero
            @test C*x.matrix*transpose(C)==y.matrix
         end
      end
      end
   end
end
