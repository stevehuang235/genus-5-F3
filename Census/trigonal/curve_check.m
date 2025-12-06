OutputFileName := "./data_unfiltered_genus/" cat Input;
InputFileName := "./data_unfiltered/" cat Input;

LinesOfInputFile := Split(Read(InputFileName), "\n");

R<X,Y,Z> := PolynomialRing(GF(3), 3);
monos5 := [X^5, X^4*Y, X^4*Z, X^3*Y^2, X^3*Y*Z, X^3*Z^2, X^2*Y^3, X^2*Y^2*Z, X^2*Y*Z^2, X^2*Z^3, X*Y^4,
X*Y^3*Z, X*Y^2*Z^2, X*Y*Z^3, X*Z^4, Y^5, Y^4*Z, Y^3*Z^2, Y^2*Z^3, Y*Z^4, Z^5];
P2<X,Y,Z> := ProjectiveSpace(GF(3), 2);


GenusCheck := function(_fsupp)
    fsupp := eval(_fsupp);
    f := R!0;
    
    for i in [1..#fsupp] do
        f := f+ fsupp[i] * monos5[i];
    end for;
    
    C_ := Scheme(P2,[f]);

    if Dimension(C_) eq 1 and IsIrreducible(C_) eq true then
        C := Curve(C_);
        F0 := FunctionField(C);
        F := AlgorithmicFunctionField(F0);

        if Genus(F) eq 5 then
            ct := "[";
            for n in [1..5] do
                ct :=  ct cat IntegerToString(NumberOfPlacesOfDegreeOneECF(F,n)) cat ",";
            end for;
            Prune(~ct);
            ct := ct cat "]"; 
            return true, ct;
        else
            return false, "";
        end if; 
    else
        return false, "";
    end if;
end function;
    

for MyLine in LinesOfInputFile do
    boo,ct := GenusCheck(MyLine);
    if boo eq true then
        to_print := "[*" cat ct cat "," cat MyLine cat "*]" cat "\n";
        fprintf OutputFileName, to_print;
    end if;
end for;

fprintf "./curve_check_log.txt", Input cat "\n";

quit;