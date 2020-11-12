
# return an element in the centralizer of x in GL(n,F) with determinant d

# first: brute force way
function _elem_given_det(x,d)
   C,e = centralizer(GL(x.parent.deg, x.parent.ring),x)
   U,fa = unit_group(x.parent.ring)
   GA,ea = sub(U, [preimage(fa,det(g)) for g in gens(C)])
   l = preimage(ea,preimage(fa,d))
   return prod([C[i]^Int(l[i]) for i in 1:ngens(C)])
end

# given an irreducible polynomial f in F[t], return the matrix that generates the field { g(C(f)) | g in F[t] }
# where C(f) = companion_matrix(f)
function _gen_matrix
end

