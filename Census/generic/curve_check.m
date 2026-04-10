OutputFileName := "./data_unfiltered_genus/" cat Input;
InputFileName := "./data_unfiltered/" cat Input;
print(Input);

LinesOfInputFile := Split(Read(InputFileName), "\n");

F := GF(3);
R<x0,x1,x2,x3,x4> := PolynomialRing(F, 5);
monos2 := [x0^2,x0*x1,x0*x2,x0*x3,x0*x4,x1^2,x1*x2,x1*x3,x1*x4,x2^2,x2*x3,x2*x4,x3^2,x3*x4,x4^2];
P4<x0,x1,x2,x3,x4> := ProjectiveSpace(F, 4);

GenusCheck := function(_fsupp)
    fsupp := eval(_fsupp);
    quadrics := [];

    for tmp in fsupp do 
        f := R!0;
        for i in [1..#tmp] do
            f +:= tmp[i] * monos2[i];
        end for;
        Append(~quadrics, f);
    end for;
    
    C_ := Scheme(P4, quadrics);
    
    if Dimension(C_) eq 1 and IsReduced(C_) eq true and IsIrreducible(C_) eq true then
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

//line_count := 0;
for MyLine in LinesOfInputFile do
    //line_count +:= 1;
    //print(line_count);
    boo,ct := GenusCheck(MyLine);
    if boo eq true then
        to_print := "[*" cat ct cat "," cat MyLine cat "*]" cat "\n";
        fprintf OutputFileName, to_print;
    end if;
end for;

fprintf "./curve_check_log.txt", Input cat "\n";

quit;