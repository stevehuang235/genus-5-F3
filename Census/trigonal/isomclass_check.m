OutputFileName := "./data_filtered/isomclass_" cat Input;
InputFileName := "./data_unfiltered_genus_updated/" cat Input;
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

// Main loop: check for pairwise isomorphism by varying over elements of the same point counts
i := 1;
while i le L do
    lst := eval(LinesOfInputFile[i]);
    ct := lst[1];
    supp := lst[2];
    F0 := FFConstruction(supp);

    tmp := [F0];
    supptmp := [supp];
    j := CountIndices(LinesOfInputFile,ct,i);
    for ind in [i..j] do
        lst2 := eval(LinesOfInputFile[ind]);
        supp2 := lst2[2];
        F02 := FFConstruction(supp2);
        if forall(u){m : m in tmp | #Isomorphisms(F02,m) eq 0} eq true then
            Append(~tmp,F02);
            Append(~supptmp,supp2);
        end if;
    end for;
    for eqn in supptmp do
        OneLinePrint(eqn);
        aut_ord := CountAutoms(FFConstruction(eqn));
        fprintf OutputFileName, "," cat "%o" cat "\n", aut_ord;
    end for;
    i := j + 1;
end while;

fprintf "./isomclasscheck_log.txt", Input cat "\n";

quit;