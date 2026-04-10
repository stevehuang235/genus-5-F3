
// ls ./data_unfiltered_genus_updated/batch_02/ | parallel -j50 "magma -b Input:={} isomclass_check.m&"

//OutputFileName := "./data_filtered_v2/batch_01/isomclass_" cat Input;
//InputFileName := "./data_unfiltered_genus_updated_2/batch_01/" cat Input;

load "pattern.m";
load "hyperplanes.m";
load "poly_pattern.m";
load "hashing.m";
InputFileName := "./data_unfiltered_genus_updated/" cat Input;
OutputFileName := "./data_unfiltered_genus_hashes/hash_" cat Input;
HashesFileName := "num_hashes_test.txt";
LinesOfInputFile := Split(Read(InputFileName), "\n");

// Count number of lines in text file
function LineCount(F)
    FP := Open(F, "r");
    count := 0;
    while true do
        line := Gets(FP);
        if IsEof(line) then
            break;
        end if;
        count +:= 1;
    end while;
    return count;
end function;

L := LineCount(InputFileName);

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

// Each line is a support ordered by point count, so we need to get starting and ending indices
function CountIndices(TxtFile, InitialPointCounts,StartingIndex)
    L := #TxtFile;
    for k in [StartingIndex..L] do
        tmp := eval(TxtFile[k]);
        if tmp[1] eq InitialPointCounts then
            continue;
        else
            return k-1;
        end if;
    end for;
    return L;
end function;

function CountAutoms(F)
    try
        return #AutomorphismGroup(F);
    catch e
        /* An error can occur when Automorphisms(F) returns a list with repetitions. */
        L := Automorphisms(F);
        L1 := [* *];
        for i in L do
            match := false;
            for j in L1 do
                if Equality(i, j) then
                    match := true;
                    break;
                end if;
            end for;
            if not match then
                Append(~L1, i);
            end if;
        end for;
        return #L1;
    end try;
end function;


OneLinePrint := procedure(list);
  fprintf OutputFileName, "[";
    for i in [1..#list] do
      if i eq #list then
        str := "%o]";
      else
        str := "%o, ";
      end if;
      fprintf OutputFileName, str, list[i]; 
  end for;
end procedure;


/*
for i in [1..L] do 
    lst := eval(LinesOfInputFile[i]);
    ct := lst[1];
    supp := lst[2];
    print(supp);
    try
        CountPatterns(supp);
    catch e 
        print("something went wrong");
    end try;
end for;
*/

// Main loop: check for pairwise isomorphism by varying over elements of the same point counts

procedure IsomCheck(supp, i, j)
    hyper_patterns := AssociativeArray();
    for ind in [i..j] do
        lst := eval(LinesOfInputFile[ind]);
        supp := lst[2];
        //try 
            pattern := CountPatterns(supp);
            if not IsDefined(hyper_patterns, pattern) then
                hyper_patterns[pattern] := [*supp*];
            else 
                Append(~hyper_patterns[pattern], supp); 
            end if;
        //catch e 
            //print("something went wrong");
            //print(e);
            //print(supp);
        //end try;
    end for;

    counter := 0;
    for key in Keys(hyper_patterns) do
        counter +:= 1;
        for supp in hyper_patterns[key] do 
            OneLinePrint(supp);
            fprintf OutputFileName, "," cat "%o" cat "\n", counter;
        end for;
        
        /*
        supps := hyper_patterns[key];
        F0 := FFConstruction(supps[1]);
        tmp := [F0];
        supptmp := [supps[1]];
        for i in [2..#supps] do 
            supp2 := supps[i];
            F02 := FFConstruction(supp2);
            if forall(u){m : m in tmp | #Isomorphisms(F02,m) eq 0} eq true then
                Append(~tmp,F02);
                Append(~supptmp,supp2);
            end if;
        end for;

        for i in [1..#supptmp] do 
            OneLinePrint(supptmp[i]);
            aut_ord := CountAutoms(tmp[i]);
            fprintf OutputFileName, "," cat "%o" cat "\n", aut_ord;
        end for;
        */
    end for;
    fprintf HashesFileName, "%o" cat "\n", counter;

end procedure;

i := 1;
while i le L do
    lst := eval(LinesOfInputFile[i]);
    ct := lst[1];
    supp := lst[2];
    
    //F0 := FFConstruction(supp);

    //tmp := [F0];
    //supptmp := [supp];
    j := CountIndices(LinesOfInputFile,ct,i);
    IsomCheck(supp, i, j);
    /*
    tmp := AssociativeArray();
    tmp[CountPatterns(supp)] := supp;
    
    for ind in [i..j] do
        lst2 := eval(LinesOfInputFile[ind]);
        supp2 := lst2[2];
        try 
            pattern := CountPatterns(supp2);
            if not IsDefined(tmp, pattern) then
                tmp[pattern] := supp2;
            end if;
        catch e 
            print("something went wrong");
        end try;
    end for;
    print(#Keys(tmp));
    for key in Keys(tmp) do
        eqn := tmp[key];
        //OneLinePrint(eqn);
        aut_ord := CountAutoms(FFConstruction(eqn));
        //fprintf OutputFileName, "," cat "%o" cat "\n", aut_ord;
    end for;
    */
    i := j + 1;
end while;

fprintf "./hashes_log.txt", Input cat "\n";



quit;