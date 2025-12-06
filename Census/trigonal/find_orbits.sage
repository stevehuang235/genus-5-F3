# Functions that computes the orbit of a group acting on a collection of certain vector spaces 
# This code is specialized to the rational g26 case 

def minimal_generators(G):
    gapG = G.gap()
    gens = gapG.MinimalGeneratingSet()
    return [G(g) for g in gens]

def act(g, V0):
    """
    Returns the output from the action of an element g in G0 on a given basis element 
    """
    g = g.matrix()
    V0_new = []
    
    f = P(sum(V0[i]*monos5[i] for i in range(21)))
    f_new = f.subs({X: sum(c*m for c,m in zip(vars,g.row(0))), Y: sum(c*m for c,m in zip(vars,g.row(1))),
    Z: sum(c*m for c,m in zip(vars,g.row(2)))})
        
    for mon in monos5:
        V0_new.append(f_new.monomial_coefficient(mon))

    return V0_new

# Function to compute the orbit of an element
def compute_orbit_ele(V0, gens):
    return [tuple(act(g, V0)) for g in gens]

def compute_orbits(L, V, G, b): 
    """
    INPUTS: 
    * "L" -- list, list of vector spaces on which G0 acts 
    * "V" -- vector space, ambient space where everything lives in 
    """
    gens = minimal_generators(G) if G.cardinality() != 1 else [G.one()]

    graph = Graph()
    for V0 in L:
        V0_mod = V0 + b
        graph.add_vertex(tuple(V0_mod))
        orbit = compute_orbit_ele(V0_mod, gens)
        for V0_new in orbit:
            if V0_new != tuple(V0_mod):
                graph.add_vertex(V0_new)
                graph.add_edge(tuple(V0_mod), V0_new)
    
    return [component[-1] for component in graph.connected_components(sort=True)]