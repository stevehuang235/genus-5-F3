// ------------------------------------------------------------
// Precompute factorization patterns of all degree-8 polynomials
// over F_3, keyed by their coefficient tuples.
// Requires: patterns.m providing IndexOfPattern(D).
// ------------------------------------------------------------

// Load your pattern machinery (adjust path/filename as needed):
//load "pattern.m";   // defines IndexOfPattern, Patterns, PatternIndex, etc.

// Base field and polynomial ring
F := GF(3);
R<x> := PolynomialRing(F);

// Associative array: key = <c0,..,c8>, value = pattern index
PolyPattern := AssociativeArray();

// Total number of coefficient tuples (c0..c8) over F_3
total := 3^9;

// Enumerate all 9-tuples (c0..c8) via base-3 expansion of an integer.
for n in [0..total-1] do
    coeffs := [];
    m := n;

    // Extract coefficients in base 3: c0 + c1*x + ... + c8*x^8
    for i in [0..8] do
        c := F!(m mod 3);
        Append(~coeffs, c);
        m div:= 3;
    end for;

    // coeffs[i+1] corresponds to x^i; coeffs[9] is the coefficient of x^8
    //if coeffs[9] eq 0 then
        // Leading coefficient 0 => degree < 8, skip
    //    continue;
    //end if;

    // Build the polynomial f(x) = c0 + c1*x + ... + c8*x^8
    f := &+[ coeffs[i+1] * x^i : i in [0..8] ];

    if f eq 0 then 
        continue;
    end if;

    // Factor f over F_3
    fac := Factorization(f);

    // Build pattern D: Assoc from <deg, multiplicity> -> count
    D := AssociativeArray();
    for pf in fac do
        g := pf[1];
        e := pf[2];
        d := Degree(g);
        key := <d, e>;
        if IsDefined(D, key) then
            D[key] +:= 1;
        else
            D[key] := 1;
        end if;
    end for;

    // Map this pattern to its canonical index from patterns.m
    //idx := IndexOfPattern(D);

    // Use the coefficient tuple <c0,...,c8> as the key
    coeff_tuple := < coeffs[i] : i in [1..9] >;
    PolyPattern[coeff_tuple] := D;
end for;

// At this point, PolyPattern is filled.
// Example lookup:
//    f := x^8 + 2*x^3 + 1;
//    coeffs := < 1,0,0,2,0,0,0,0,1 >;  // 1 + 2x^3 + x^8
//    idx := PolyPattern[coeffs];
