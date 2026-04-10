OutputFileName := "./data_generic/" cat Input;
Input := "./data_hashes/" cat Input;
LinesOfInputFile := Split(Read(Input), "\n");

F := GF(3);
R<x0,x1,x2,x3,x4> := PolynomialRing(F, 5);
monos2 := [x0^2,x0*x1,x0*x2,x0*x3,x0*x4,x1^2,x1*x2,x1*x3,x1*x4,x2^2,x2*x3,x2*x4,x3^2,x3*x4,x4^2];
P4<x0,x1,x2,x3,x4> := ProjectiveSpace(F, 4);

// Quick construction of function field with IsIsomorphic functionality
function FFConstruction(fsupp)
    quadrics := [];

    for tmp in fsupp do 
        f := R!0;
        for i in [1..#tmp] do
            f +:= tmp[i] * monos2[i];
        end for;
        Append(~quadrics, f);
    end for;
    
    C_ := Scheme(P4, quadrics);
    C := Curve(C_);
    F0 := FunctionField(C);

    return AlgorithmicFunctionField(F0);
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
    F := FFConstruction(fsupp);
    B := Parent(DefiningPolynomial(F));
    A := BaseRing(B);
    AssignNames(~A, ["x"]);
    AssignNames(~B, ["y"]);
    poly := DefiningPolynomial(F);
    cpc := ([NumberOfPlacesOfDegreeOneECF(F,n) : n in [1..5]]);
    G := IdentifyGroup(AutomorphismGroupCorrected(F));
    fprintf OutputFileName, "[" cat "%o" cat "," cat "'" cat "%o" cat "'" cat "," cat "%o" cat "]" cat "\n", cpc, G, poly;
end for;

fprintf "./database_gen_log.txt", Input cat "\n";

quit;