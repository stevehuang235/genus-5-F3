import ast
from collections import defaultdict

F = GF(3)
V = Matrix(F,0,15).right_kernel()
custom_function = """
F := GF(3);
R<x0,x1,x2,x3,x4> := PolynomialRing(F, 5);
monos2 := [x0^2,x0*x1,x0*x2,x0*x3,x0*x4,x1^2,x1*x2,x1*x3,x1*x4,x2^2,x2*x3,x2*x4,x3^2,x3*x4,x4^2];


P4<x0,x1,x2,x3,x4> := ProjectiveSpace(F, 4);

function Process(eqns) 
    quadrics := [];

    for tmp in eqns do 
        f := R!0;
        for i in [1..#tmp] do
            f +:= tmp[i] * monos2[i];
        end for;
        Append(~quadrics, f);
    end for;
    X := Scheme(P4, quadrics);
    
    return IsIrreducible(X) and (Dimension(SingularSubscheme(X)) lt 1);
end function;
"""
magma.eval(custom_function)

with open('./2-tuples-f3.txt', 'r') as file:
    content = file.read().strip()
    tuples_list = ast.literal_eval(content)

    with open(f"./generic_unfiltered.txt", 'w') as file_out:
        for tup in tuples_list: 
            eqns = [list(tmp) for tmp in tup]
            magma.eval("eqns := {};".format(eqns))
            bool = magma("Process(eqns)").sage()
            if bool:
                subspace = V.subspace(eqns)
                V0 = V.quotient(subspace)
                for v in V0:
                    tmp = eqns + [list(V0.lift(v))]
                    file_out.write(str(tmp))
                    file_out.write("\n")

        file_out.close()
file.close()