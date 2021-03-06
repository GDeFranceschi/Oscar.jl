@testset "Polynomial Orderings" begin
  R, (x, y, z) = PolynomialRing(QQ, 3)
  f = x*y + 5*z^3

  @test collect(monomials(f, :lex)) == [ x*y, z^3 ]
  @test collect(monomials(f, :revlex)) == [ z^3, x*y ]
  @test collect(monomials(f, :deglex)) == [ z^3, x*y ]
  @test collect(monomials(f, :degrevlex)) == [ z^3, x*y ]
  @test collect(monomials(f, :neglex)) == [ z^3, x*y ]
  @test collect(monomials(f, :negrevlex)) == [ x*y, z^3 ]
  @test collect(monomials(f, :negdeglex)) == [ x*y, z^3 ]
  @test collect(monomials(f, :negdegrevlex)) == [ x*y, z^3 ]

  w = [ 1, 2, 1 ]
  @test collect(monomials(f, :weightlex, w)) == [ x*y, z^3 ]
  @test collect(monomials(f, :weightrevlex, w)) == [ x*y, z^3 ]

  M = [ 1 1 1; 1 0 0; 0 1 0 ]
  @test collect(monomials(f, M)) == collect(monomials(f, :deglex))

  @test collect(terms(f, :deglex)) == [ 5z^3, x*y ]
  @test collect(exponent_vectors(f, :deglex)) == [ [ 0, 0, 3 ], [ 1, 1, 0 ] ]
  @test collect(coeffs(f, :deglex)) == [ QQ(5), QQ(1) ]

  Fp = FiniteField(7)
  R, (x, y, z) = PolynomialRing(Fp, 3, ordering = :deglex)
  f = x*y + 5*z^3
  @test collect(monomials(f, :lex)) == [ x*y, z^3 ]
  @test Oscar.leading_monomial(f, :lex) == x*y
  @test Oscar.leading_coeff(f, :lex) == Fp(1)
  @test Oscar.leading_term(f, :lex) == x*y

end
