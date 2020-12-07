export preserved_quadratic_forms, preserved_sesquilinear_forms

GAP.Globals.LoadPackage(GAP.julia_to_gap("Forms"))

"""
    preserved_quadratic_forms(G::MatrixGroup)
Return a generating set for the vector space of quadratic forms preserved by `G`.

!!! warning "Note:"
    The process involves random procedures, so the function may return different outputs every time.

"""
function preserved_quadratic_forms(G::MatrixGroup)
   L = GAP.Globals.PreservedQuadraticrForms(G.X)
   return [G.mat_iso(GAP.Globals.GramMatrix(L[i])) for i in 1:length(L)]
end

"""
    preserved_sesquilinear_forms(G::MatrixGroup)
Return a generating set for the vector space of sesquilinear forms preserved by `G`.

!!! warning "Note:"
    The process involves random procedures, so the function may return different outputs every time.

"""
function preserved_sesquilinear_form(G::MatrixGroup)
   L = GAP.Globals.PreservedSesquilinearForms(G.X)
   return [G.mat_iso(GAP.Globals.GramMatrix(L[i])) for i in 1:length(L)]
end
