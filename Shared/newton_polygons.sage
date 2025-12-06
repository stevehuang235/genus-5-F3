def newton_polygons(g):
    """
    Returns the set of all possible Newton polygons of g-dimensional abelian varieties. 

    INPUT: 
    
    * "g" -- a positive integer 
    """

    polygons = set()
    possible_slopes = [(0, 1)]
    for s in range(2, g+1):
        for r in range(1, int(floor(s/2)) + 1):
            num = QQ(r/s)
            if num == QQ(1/2):
                possible_slopes.append((QQ(1/2), 1))
            else:
                possible_slopes.append((num, num.denominator()))

    # (each r/s is paired with (s-r)/s)
    unique_slopes = sorted(set((frac, 1-frac, denom) for frac, denom in possible_slopes))

    for weights in WeightedIntegerVectors(g, [denom for _, _, denom in unique_slopes]):
        slopes = []

        for (slope1, slope2, s), weight in zip(unique_slopes, weights):
            slopes.extend([slope1, slope2] *  s * weight)
            if len(slopes) == 2*g and sum(slopes) == g: 
                polygons.add(tuple(sorted(slopes)))

    return sorted(polygons)

    
