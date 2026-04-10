// ----------------------------------------------------------------------
// Hyperplanes and canonical representatives
// ----------------------------------------------------------------------
F := GF(3);
V5 := VectorSpace(F, 5);
ZZ := Integers();
// List of variables and degree-2 monomials, same order as in Sage.
P<x0, x1, x2, x3, x4> := PolynomialRing(F, 5);
vars  := [ P.1, P.2, P.3, P.4, P.5 ];
mons2 := [x0^2, x0*x1, x0*x2, x0*x3, x0*x4, x1^2, x1*x2, x1*x3, x1*x4, x2^2, x2*x3, x2*x4, x3^2, x3*x4, x4^2];

// Canonical representative for a nonzero row vector in V5*.
// Multiply by 2 if the first nonzero entry is 2, so that it becomes 1.
function Canonical(v)
    for a in Eltseq(v) do
        if a ne 0 then
            if a eq F!2 then
                v := 2*v;
            end if;
            return Eltseq(v);
        end if;
    end for;
    return [];
end function;

// Enumerate all hyperplanes in P^4(F_3) as canonical 5-tuples.
function AllHyperplanes()
    seen := {};
    H := [];
    for v in V5 do
        if v ne V5!0 then
            t := Canonical(v);
            if not t in seen then
                Include(~seen, t);
                Append(~H, t);
            end if;
        end if;
    end for;
    return H;
end function;

// Global list of hyperplanes, analogous to Sage's H.
Hyperplanes := AllHyperplanes();

changevars := AssociativeArray();  // keys will be the 5-tuples from Hyperplanes

// Build the 15x15 change-of-basis matrix M for a given hyperplane ell_seq.
function BuildChangeMatrix(ell_seq)
    ell := ell_seq;

    // Construct a 5x5 matrix A whose last row is ell, and the first four
    // rows are the standard basis, then tweak until full rank.
    rows := [
        [1, 0, 0, 0, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        ell
    ];
    A := Matrix(F, rows);

    while Rank(A) lt 5 do
        k := Random([1..4]);
        rows[k] := [ Random(F), Random(F), Random(F), Random(F), Random(F) ];
        A := Matrix(F, rows);
    end while;

    Ainv := A^-1;

    // New variables: x_i' = sum_j Ainv[i,j] * x_j
    new_vars := [ &+[ Ainv[i, j] * vars[j] : j in [1..5] ] : i in [1..5] ];

    // Build the 15x15 matrix M: columns are coords of transformed monomials
    // expressed in the basis mons2.
    M := ZeroMatrix(F, 15, 15);
    for j in [1..15] do
        mon := mons2[j];
        mon_new := Evaluate(mon, new_vars);  // substitute x -> Ainv * x

        // Extract coefficients of mon_new in the basis mons2.
        col := [ MonomialCoefficient(mon_new, mons2[k]) : k in [1..15] ];
        for i in [1..15] do
            M[i, j] := col[i];
        end for;
    end for;

    return M;
end function;

// Precompute changevars for ALL hyperplanes once, like in Sage.
for ell in Hyperplanes do
    changevars[ell] := BuildChangeMatrix(ell);
end for;