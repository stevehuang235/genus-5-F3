
OutputFileName := "./data_trigonal.txt";
Input := "./data_filtered/" cat Input;
LinesOfInputFile := Split(Read(Input), "\n");

R<X,Y,Z> := PolynomialRing(GF(3), 3);
monos5 := [X^5, X^4*Y, X^4*Z, X^3*Y^2, X^3*Y*Z, X^3*Z^2, X^2*Y^3, X^2*Y^2*Z, X^2*Y*Z^2, X^2*Z^3, X*Y^4,
X*Y^3*Z, X*Y^2*Z^2, X*Y*Z^3, X*Z^4, Y^5, Y^4*Z, Y^3*Z^2, Y^2*Z^3, Y*Z^4, Z^5];
P2<X,Y,Z> := ProjectiveSpace(GF(3), 2);

// Quick construction of function field with IsIsomorphic functionality
function FFConstruction(fsupp)
    f := R!0;
    
    for i in [1..#fsupp] do
        f := f + fsupp[i] * monos5[i];
    end for;
    
    C_ := Scheme(P2,[f]);
    C := Curve(C_);
    F0 := FunctionField(C);

    return f, AlgorithmicFunctionField(F0);
end function;

function AutomorphismGroupCorrected(F)
    try
        return AutomorphismGroup(F);
    catch e
        /* When this error occurs, Automoprhisms(F) is still assumed to return a list
           of automorphisms which generates the full groups, but may be incomplete and/or
           include repetitions. */
        L := Automorphisms(F);
        G := FreeGroup(#L);
        rels := [];
        for i in [1..#L] do
            for j in [1..#L] do
                g := Composition(L[i], L[j]);
                for k in [1..#L] do
                    if Equality(g, L[k]) then
                        Append(~rels, G.i*G.j*G.k^(-1));
                    end if;
                end for;
            end for;
        end for;
        return quo<G|rels>;
    end try;
end function;

for MyLine in LinesOfInputFile do
    fsupp := eval(MyLine);
    poly, F := FFConstruction(fsupp);
    cpc := ([NumberOfPlacesOfDegreeOneECF(F,n) : n in [1..5]]);
    G := IdentifyGroup(AutomorphismGroupCorrected(F));
    fprintf OutputFileName, "[" cat "%o" cat "," cat "'" cat "%o" cat "'" cat "," cat "%o" cat "]" cat "\n", cpc, G, poly;
end for;

quit;