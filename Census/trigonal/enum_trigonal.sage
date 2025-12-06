load("./find_orbits.sage")

F3 = GF(3)
P.<X,Y,Z> = PolynomialRing(F3, 3, order='lex')
vars = [X,Y,Z]
monos5 = P.monomials_of_degree(5)  # 21 monomials of degree 5
monos5.reverse()
G = GL(3,F3)

# enforce singular point at [0:0:1]
a14 = [0]*14 + [1] + [0]*6  # a14 = 0
a19 = [0]*19 + [1] + [0]  # a19 = 0
a20 = [0]*20 + [1]  # a20 = 0 

# coefficients of the tangent cone a9x^2 + 2a13xy + a18y^2 = 0
a13 = [0]*13 + [1] + [0]*7  
a9 = [0]*9  + [1] + [0]*11  
a18 = [0]*18 + [1] + [0]*2  

M = matrix(F3, [a14, a19, a20, a13, a9, a18])
V = M.right_kernel()
V_basis = V.basis()
V_lst = list(V)

# Case 1: split node, tangent cone in affine coordinate at (0,0) is x^2 - y^2
b1 = M.solve_right(vector(F3, [0, 0, 0, 0, 1, 2])) # a9 = 1, a18=2
G1 = []
for g in G: 
    #if act(g, b1) == list(b1) and all(vector(F3, act(g, v)) in V for v in V.basis()):
    if (vector(F3, act(g,b1)) - b1) in V and all(vector(F3, act(g, v)) in V for v in V.basis()):
        G1.append(g)
G1 = G.subgroup(G1)

orbit_reps_1 = compute_orbits(V_lst, V, G1, b1)  
print(f"The number of orbit reps in Case 1 is {len(orbit_reps_1)}")
with open("./split_node_unfiltered.txt", 'w') as f1:
    for orbit in orbit_reps_1:
        tmp = list(orbit)
        f1.write(str(tmp))
        f1.write("\n")
f1.close()

# Case 2: non-split node, tangent cone in affine coordinate at (0,0) is x^2 + y^2
b2 = M.solve_right(vector(F3, [0, 0, 0, 0, 1, 1]))  # now a18 = 1
G2 = []
for g in G: 
    if (vector(F3, act(g,b2)) - b2) in V and all(vector(F3, act(g, v)) in V for v in V.basis()):
        G2.append(g)
G2 = G.subgroup(G2)
orbit_reps_2 = compute_orbits(V_lst, V, G2, b2)  
print(f"The number of orbit reps in Case 2 is {len(orbit_reps_2)}")
with open("./nonsplit_node_unfiltered.txt", 'w') as f2:
    for orbit in orbit_reps_2:
        tmp = list(orbit)
        f2.write(str(tmp))
        f2.write("\n")
f2.close()

# Case 3: potential cusp, tangent cone in affine coordinate at (0,0) is y^2
b3 = M.solve_right(vector(F3, [0, 0, 0, 0, 0, 1]))  # now a9=a13=0 a18 = 1
G3 = []
for g in G: 
    if (vector(F3, act(g,b3)) - b3) in V and all(vector(F3, act(g, v)) in V for v in V.basis()):
        G3.append(g)
G3 = G.subgroup(G3)
orbit_reps_3 = compute_orbits(V_lst, V, G3, b3)  
print(f"The number of orbit reps in Case 3 is {len(orbit_reps_3)}")
with open("./potential_cusp_unfiltered.txt", 'w') as f3:
    for orbit in orbit_reps_3:
        tmp = list(orbit)
        f3.write(str(tmp))
        f3.write("\n")
f3.close()