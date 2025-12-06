OutputFileName := "./data_hyperelliptic.txt";
load "genus5_Howe";

function FFConstruction(C)
    return AlgorithmicFunctionField(FunctionField(C));
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

// tmp := [];

stacky_count := 0;

for curve in genus5 do
    poly := curve[1]; aut_ord := curve[2];
    F := FFConstruction(HyperellipticCurve(poly));
    /*
    // check that there is no pairwise isomorphisms
    if #tmp eq 0 then 
        Append(~tmp, F);
    else 
        for F0 in tmp do 
            assert #Isomorphisms(F0, F) eq 0;
        end for;
        Append(~tmp, F);
    end if;
    */
    cpc := ([NumberOfPlacesOfDegreeOneECF(F,n) : n in [1..5]]);
    G := IdentifyGroup(AutomorphismGroupCorrected(F));
    assert G[1] eq aut_ord;  // verify that the order of the automorphism group agrees with Howe's data 
    stacky_count +:= 1/aut_ord;
    fprintf OutputFileName, "[" cat "%o" cat "," cat "'" cat "%o" cat "'" cat "," cat "%o" cat "]" cat "\n", cpc, G, poly;
end for;

printf "stacky count of the hyperelliptic locus is %o", stacky_count; 

quit;

