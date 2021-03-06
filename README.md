# Oscar.jl

Welcome to the OSCAR project, a visionary new computer algebra system
which combines the capabilities of four cornerstone systems: Gap,
Polymake, Antic and Singular.

## Installation

* Download Julia version 1.3 or higher
* In Julia type

```
julia> using Pkg
julia> Pkg.add("Oscar")
julia> using Oscar
```

or, to be at cutting edge:

```
julia> using Pkg
julia> Pkg.add(Oscar#master)
julia> using Oscar
```

Installation will take a couple of minutes.

## Current requirements

* we support Julia 1.3 and we try to get ready a.s.a.p. for the upcoming Julia 1.4
* target platforms: x86_64-linux-gnu, x86_64-apple-darwin14
* CxxWrap 0.9.0
* travis: distribution ubuntu bionic

Some explanations:

* Julia 1.0 is the latest LTS version, but it seems unlikely that GAP (with its garbage collection) can ever be supported
* Julia 1.1 and 1.2 do not receive any back ports any more
* Windows support only through WSL (i.e., this is covered by x86_64-linux-gnu)
* The travis configuration is important for building CxxWrap; otherwise the compiler is too old.


## Examples of usage

```
julia> using Oscar
...
 -----    -----    -----      -      -----   
|     |  |     |  |     |    | |    |     |  
|     |  |        |         |   |   |     |  
|     |   -----   |        |     |  |-----   
|     |        |  |        |-----|  |   |    
|     |  |     |  |     |  |     |  |    |   
 -----    -----    -----   -     -  -     -  

...combining (and extending) GAP, Hecke, Nemo, Polymake and Singular
Version 0.2.0 ... 
 ... which comes with absolutely no warranty whatsoever
Type: '?Oscar' for more information
(c) 2019-2020 by The Oscar Development Team


julia> k, a = quadratic_field(-5)
(Number field over Rational Field with defining polynomial x^2+5, sqrt(-5))

julia> zk = maximal_order(k)
Maximal order of Number field over Rational Field with defining polynomial x^2+5
with basis nf_elem[1, sqrt(-5)]

julia> factorisations(zk(6))
2-element Array{Fac{NfAbsOrdElem{AnticNumberField,nf_elem}},1}:
 -1 * (2) * (-3)
 -1 * (sqrt(-5)+1) * (sqrt(-5)-1)

julia> Qx, x = PolynomialRing(QQ, :x=>1:2)
(Multivariate Polynomial Ring in x1, x2 over Rational Field, fmpq_mpoly[x1, x2])

julia> R = grade(Qx, [1,2])
Multivariate Polynomial Ring in x1, x2 over Rational Field graded by 
        x1 -> [1]
        x2 -> [2]

julia> f = R(x[1]^2+x[2])
x1^2 + x2
julia> degree(f)
graded by [2]

julia> F = FreeModule(R, 1)
Free module of rank 1 over R, graded as R^1([0])

julia> s = sub(F, [f*F[1]])
Subquotient by Array of length 1
1 -> (x1^2 + x2)*e[1]

a> mH(H[1])
Map with following data
Domain:
=======
s
Codomain:
=========
Subquotient of Array of length 1
1 -> (1)*e[1]
 by Array of length 1
1 -> (x1^2 + x2)*e[1]
defined on the Singular side

julia> H, mH = hom(s, quo(F, s))
(hom of (s, Subquotient of Array of length 1
1 -> (1)*e[1]
 by Array of length 1
1 -> (x1^2 + x2)*e[1]
defined on the Singular side

), Map from
H to Set of all homomorphisms from Subquotient by Array of length 1
1 -> (x1^2 + x2)*e[1]
defined on the Singular side

 to Subquotient of Array of length 1
1 -> (1)*e[1]
 by Array of length 1
1 -> (x1^2 + x2)*e[1]
defined on the Singular side

 defined by a julia-function with inverse
)

julia> D = decoration(H)
GrpAb: Z

julia> homogenous_component(H, D[0])
(H_[0] of dim 2, Map from
H_[0] of dim 2 to H defined by a julia-function with inverse
)
```

Of course, the cornerstones are also available directly:

```
julia> C = Polymake.polytope.cube(3);

julia> C.F_VECTOR
pm::Vector<pm::Integer>
8 12 6

julia> RP2 = Polymake.topaz.real_projective_plane();

julia> RP2.HOMOLOGY
PropertyValue wrapping pm::Array<polymake::topaz::HomologyGroup<pm::Integer>>
({} 0)
({(2 1)} 0)
({} 0)

```
