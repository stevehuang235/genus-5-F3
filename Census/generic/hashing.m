// hashing.m
// Magma translation of hashing.sage, with full precomputation of changevars
// and coefficient-level change of basis (15x15 matrices), mirroring Sage.
//
// Usage:
//   Attach("patterns.m");
//   Attach("hashing.m");
//   hashing1 := CountPatterns(quad_coeffs1);
//   hashing2 := CountPatterns(quad_coeffs2);
//   ...

//load "pattern.m";
// ----------------------------------------------------------------------
// Basic setup
// ----------------------------------------------------------------------
F := GF(3);
V5 := VectorSpace(F, 5);
ZZ := Integers();

P<x0, x1, x2, x3, x4> := PolynomialRing(F, 5);
P4<z0, z1, z2, z3> := PolynomialRing(F, 4);
P3<y1, y2, y3> := PolynomialRing(F, 3);
m := ideal< P3 | y1, y2, y3>;

// List of variables and degree-2 monomials, same order as in Sage.
vars  := [ P.1, P.2, P.3, P.4, P.5 ];
mons2 := [x0^2, x0*x1, x0*x2, x0*x3, x0*x4, x1^2, x1*x2, x1*x3, x1*x4, x2^2, x2*x3, x2*x4, x3^2, x3*x4, x4^2];

// For quotient-dimension (vector space dimension of P3 / I).
function QuotientDim(P, I)
    Q := quo< P | I >;
    card := #Q;             // cardinality of the finite ring
    /*
    q := #F;                // = 3 for GF(3)

    // Compute n such that card = q^n
    dim := 0;
    s := 1;
    while s lt card do
        s *:= q;
        dim +:= 1;
    end while;

    assert s eq card;

    return dim;
    */
    return Valuation(card, 3);
end function;

// ----------------------------------------------------------------------
// Quadrics and intersection patterns
// ----------------------------------------------------------------------

// Convert a list of 15 coefficients into a quadric in P.
function QuadricFromCoeffs(coeffs)
    poly := &+[ F!coeffs[i] * mons2[i] : i in [1..#mons2] ];
    return poly;
end function;

// Intersection with a hyperplane L:
//   Q1_coeffs, Q2_coeffs, Q3_coeffs : sequences of length 15
//   ell_seq : the 5-tuple hyperplane, exactly from Hyperplanes
// Returns an associative array mapping <deg, mult> -> count,
// where deg = degree of the radical component, mult = multiplicity.
function PolyIntersect(Q1_coeffs, Q2_coeffs, Q3_coeffs, ell_seq)
    // Look up the precomputed 15x15 change-of-basis matrix for ell_seq.    
    M := changevars[ell_seq];

    // Transform the coefficient vectors: new = M * old.
    v1 := M * Matrix(F, 15, 1, Q1_coeffs);
    v2 := M * Matrix(F, 15, 1, Q2_coeffs);
    v3 := M * Matrix(F, 15, 1, Q3_coeffs);

    // Rebuild the transformed quadrics as polynomials in P.
    Q1 := &+[ v1[i,1] * mons2[i] : i in [1..15] ];
    Q2 := &+[ v2[i,1] * mons2[i] : i in [1..15] ];
    Q3 := &+[ v3[i,1] * mons2[i] : i in [1..15] ];

    points := AssociativeArray();
    Q1 := Evaluate(Q1, [P4.1, P4.2, P4.3, P4.4, 0]);
    Q2 := Evaluate(Q2, [P4.1, P4.2, P4.3, P4.4, 0]);
    Q3 := Evaluate(Q3, [P4.1, P4.2, P4.3, P4.4, 0]);

    I := ideal< P4 | Q1, Q2, Q3>;
    Is, Is_rad := PrimaryDecomposition(I);

    for i in [1..#Is] do
        deg, _  := HilbertPolynomial(Is_rad[i]);   // degree (geom. multiplicity)
        dim, _  := HilbertPolynomial(Is[i]);      // length
        mult := dim div deg;

        key := <ZZ!deg, ZZ!mult>;
        if IsDefined(points, key) then
            points[key] +:= 1;
        else
            points[key] := 1;
        end if;
    end for;

    /*
    // Dehomogenise in the chart x0 = 1, and impose x4 = 0 (restrict to L).
    Q1_aff := Evaluate(Q1, [ F!1, P3.1, P3.2, P3.3, F!0 ]);
    Q2_aff := Evaluate(Q2, [ F!1, P3.1, P3.2, P3.3, F!0 ]);
    Q3_aff := Evaluate(Q3, [ F!1, P3.1, P3.2, P3.3, F!0 ]);

    I_aff := ideal< P3 | Q1_aff, Q2_aff, Q3_aff >;
    points := AssociativeArray();

    // This mirrors your Sage code: it only handles the case length = 8.
    if QuotientDim(P3, I_aff) eq 8 then 
        Is, Is_rad := PrimaryDecomposition(I_aff);
        for i in [1..#Is] do
            deg  := QuotientDim(P3, Is_rad[i]);   // degree (geom. multiplicity)
            dim  := QuotientDim(P3, Is[i]);      // length
            mult := dim div deg;

            key := <deg, mult>;
            if IsDefined(points, key) then
                points[key] +:= 1;
            else
                points[key] := 1;
            end if;
        end for;
    else
        Q1 := Evaluate(Q1, [P4.1, P4.2, P4.3, P4.4, 0]);
        Q2 := Evaluate(Q2, [P4.1, P4.2, P4.3, P4.4, 0]);
        Q3 := Evaluate(Q3, [P4.1, P4.2, P4.3, P4.4, 0]);

        I := ideal< P4 | Q1, Q2, Q3>;
        Is, Is_rad := PrimaryDecomposition(I);

        for i in [1..#Is] do
            deg, _  := HilbertPolynomial(Is_rad[i]);   // degree (geom. multiplicity)
            //assert deg eq HilbPoly(Is_rad[i]);
            dim, _  := HilbertPolynomial(Is[i]);      // length
            //assert dim eq HilbPoly(Is[i]);
            //deg := HilbPoly(Is_rad[i]);   // degree (geom. multiplicity)
            //dim := HilbPoly(Is[i]);      // length
            mult := dim div deg;

            key := <ZZ!deg, ZZ!mult>;
            if IsDefined(points, key) then
                points[key] +:= 1;
            else
                points[key] := 1;
            end if;
        end for;
    end if;
    */

    return points;
    //return [* <k, points[k]> : k in Keys(points) *];
end function;

S<u1,u2,u3,t> := PolynomialRing(F, 4);
S0<z1, z2, z3> := PolynomialRing(F, 3);
phi := hom< P -> S | [1,u1,u2,u3,0] >;
phi0 := hom< P -> P3 | [0, y1, y2, y3, 0] >;
function IntersectPoly(Q1_coeffs, Q2_coeffs, Q3_coeffs, ell_seq)
    // Look up the precomputed 15x15 change-of-basis matrix for ell_seq.
    M := changevars[ell_seq];

    // Transform the coefficient vectors: new = M * old.
    v1 := M * Matrix(F, 15, 1, Q1_coeffs);
    v2 := M * Matrix(F, 15, 1, Q2_coeffs);
    v3 := M * Matrix(F, 15, 1, Q3_coeffs);

    // Rebuild the transformed quadrics as polynomials in P.
    Q1 := &+[ v1[i,1] * mons2[i] : i in [1..15] ];
    Q2 := &+[ v2[i,1] * mons2[i] : i in [1..15] ];
    Q3 := &+[ v3[i,1] * mons2[i] : i in [1..15] ];

    Q1s := phi(Q1);
    Q2s := phi(Q2);
    Q3s := phi(Q3);

    maxTries := 5;

    flag := false;
    flag1 := false;
    f := 0;
    deg_f := -1;
    I_dehom_dim := QuotientDim(S, ideal< S | Q1s, Q2s, Q3s, t>);
    //print(I_dehom_dim);
    for tries in [1..maxTries] do
    //while not flag1 do 
        // random ell = sum c_i u_i, not all c_i = 0
        repeat
            coeffs := [ Random(F) : i in [1..3] ];
        until not &and[ c eq 0 : c in coeffs ];

        ell := coeffs[1]*u1 + coeffs[2]*u2 + coeffs[3]*u3;

        I := ideal< S | Q1s, Q2s, Q3s, t - ell >;
        J := EliminationIdeal(I, 3);      // eliminate u0..u4, keep t
        B := Basis(J);

        // Look for a non-constant polynomial in t
        polys := [ g : g in B | TotalDegree(g) gt 0 ];
        //print(polys);
        if #polys eq 0 then
            continue;
        end if;

        // Take the highest-degree one (should be the interesting one)
        g := polys[#polys];
        //print(g);
        d := Degree(g);
        //f := g / LeadingCoefficient(g);
        if d eq I_dehom_dim then
            f := g / LeadingCoefficient(g);
            deg_f := d;
            flag1 := true;
            break;
        end if;
    end for;
    //end while;

    if not flag1 then 
        return false, PolyIntersect(Q1_coeffs, Q2_coeffs, Q3_coeffs, ell_seq);
    end if;
    
    f_coeffs := <F!Coefficient(f, 4, i) : i in [0..8]>;
    //print(f);
    //print(f_coeffs);
    //points := AssociativeArray();
    points := PolyPattern[f_coeffs];

    if deg_f eq 8 then 
        flag := true;
    elif deg_f eq 7 then 
        key := <1, 1>;
        if IsDefined(points, key) then 
            points[key] +:= 1;
        else; 
            points[key] := 1;
        end if ;
        flag := true;
    else  
        Q10 := phi0(Q1);
        Q20 := phi0(Q2);
        Q30 := phi0(Q3);
        
        I := ideal< P3 | Q10, Q20, Q30>;
        I := Saturation(I, m);
        //print(I);
        Is, Is_rad := PrimaryDecomposition(I);
        //print(Is);
        if #Is eq 1 then 
            deg, _ := HilbertPolynomial(Is_rad[1]);
            mult := (8-deg_f) / deg;
            key := <ZZ!deg, ZZ!mult>; 
            if IsDefined(points, key) then
                points[key] +:= 1;
            else
                points[key] := 1;
            end if;
            flag := true;
        else 
            for i in [1..#Is] do
                deg, _  := HilbertPolynomial(Is_rad[i]);   // degree (geom. multiplicity)
                dim, _  := HilbertPolynomial(Is[i]);      // length
                deg_f +:= dim;
                if deg ne 0 then 
                    mult := dim div deg;
                    key := <ZZ!deg, ZZ!mult>;
                    if IsDefined(points, key) then
                        points[key] +:= 1;
                    else
                        points[key] := 1;
                    end if;
                end if;
            end for;
            //print(deg_f);
            if deg_f eq 8 then 
                flag := true;
            else 
                points := PolyIntersect(Q1_coeffs, Q2_coeffs, Q3_coeffs, ell_seq);
            end if;
        end if;
    end if;
        /*
        for i in [1..#Is] do
            deg, _  := HilbertPolynomial(Is_rad[i]);   // degree (geom. multiplicity)
            dim, _  := HilbertPolynomial(Is[i]);      // length

            print(<deg,dim>);
            if deg ne 0 then 
                mult := dim div deg;
                key := <ZZ!deg, ZZ!mult>;
                if IsDefined(points, key) then
                    points[key] +:= 1;
                else
                    points[key] := 1;
                end if;
            end if;
        end for;
        */
   

    //return f, flag, points;
    return flag, points;
    //return flag, [* <k, points[k]> : k in Keys(points) *];
end function;

// Given Qcoeffs := [Q1, Q2, Q3] and a hyperplane ell, return the pattern index.
function ComputePattern(Qcoeffs, ell)
    //multis := PolyIntersect(Qcoeffs[1], Qcoeffs[2], Qcoeffs[3], ell);
    _, multis := IntersectPoly(Qcoeffs[1], Qcoeffs[2], Qcoeffs[3], ell);
    return IndexOfPattern(multis);
end function;

// ----------------------------------------------------------------------
// Pattern counting over all hyperplanes
// ----------------------------------------------------------------------

// Count patterns over all hyperplanes in P^4(F_3).
// Returns a sequence of length #Patterns with counts.
function CountPatterns(Qcoeffs)
    //require assigned Patterns :
    //    "Patterns is not defined. Attach \"patterns.m\" first.";

    hashing := [ 0 : i in [1..#Patterns] ];
    for ell in Hyperplanes do
        idx := ComputePattern(Qcoeffs, ell);
        hashing[idx] +:= 1;
    end for;
    return hashing;
end function;

// ----------------------------------------------------------------------
// Example coefficient triples (same structure/order as in hashing.sage)
// ----------------------------------------------------------------------

// Order of coefficients matches mons2 = [x0^2, x0*x1, ..., x4^2].
quad_coeffs := [
    [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0 ],
    [ 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0 ],
    [ 2, 0, 2, 2, 2, 2, 0, 0, 1, 1, 0, 0, 0, 2, 1 ]
];

quad_coeffs1 := [
    [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0 ],
    [ 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0 ],
    [ 2, 0, 2, 2, 2, 2, 0, 0, 1, 1, 0, 0, 0, 2, 1 ]
];

quad_coeffs2 := [
    [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0 ],
    [ 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0 ],
    [ 1, 1, 0, 0, 2, 0, 1, 0, 2, 1, 0, 0, 1, 2, 1 ]
];

quad_coeffs3 := [
    [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0 ],
    [ 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0 ],
    [ 1, 2, 1, 1, 2, 2, 2, 0, 2, 2, 0, 0, 2, 2, 1 ]
];

quad_coeffs4 := [
    [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0 ],
    [ 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 2 ],
    [ 0, 0, 0, 0, 1, 2, 2, 0, 2, 2, 0, 2, 2, 1, 0 ]
]; // should have the same patterns as quad_coeffs2

quad_coeffs_prob := [
    [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0 ],
    [ 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0 ],
    [ 0, 2, 0, 1, 1, 1, 1, 1, 2, 2, 2, 0, 0, 0, 2 ]
];
