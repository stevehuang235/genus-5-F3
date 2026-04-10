// patterns.m
// Magma translation of patterns.sage
// Generates all patterns of (degree, multiplicity) for total weight N=8.

// Return all pairs <d, m> with d * m = w.
function FactorPairs(w)
    pairs := [ <d, w div d> : d in [1..w] | w mod d eq 0 ];
    return pairs;
end function;

// Normalize a pattern D (Assoc from <d,m> to count) to a sorted sequence
// of <<d,m>, count>, sorted lex in (d,m).
function NormalizePattern(D)
    keys := SetToSequence(Keys(D));

    // Comparator must return -1,0,1 (not Bool)
    Sort(~keys, func< a, b |
        a[1] lt b[1] select -1 else
        a[1] gt b[1] select 1 else
        a[2] lt b[2] select -1 else
        a[2] gt b[2] select 1 else
        0
    >);

    seq := [ < keys[i], Integers()!D[keys[i]] >
             : i in [1..#keys] | Integers()!D[keys[i]] ne 0 ];
    return seq;
end function;

// Convert a pattern to a string id, used as a hash key.
function PatternId(D)
    return Sprintf("%o", NormalizePattern(D));
end function;

// Key for sorting patterns: multiset of (w = d*m, d, m).
function PatternKey(D)
    items := [* *];
    for k in Keys(D) do
        d := k[1];
        m := k[2];
        c := Integers()!D[k];
        for i in [1..c] do
            Append(~items, <d*m, d, m>);
        end for;
    end for;
    // Default Sort is fine here (no custom comparator).
    Sort(~items);
    return items;
end function;

// Integer partitions of N with nonincreasing parts, each <= maxPart.
function PartitionsOf(N, maxPart)
    if N eq 0 then
        return [ [] ];
    end if;
    res := [];
    limit := (N lt maxPart) select N else maxPart;
    for p in [limit..1 by -1] do
        subs := PartitionsOf(N - p, p);
        for s in subs do
            Append(~res, [p] cat s);
        end for;
    end for;
    return res;
end function;

// All r-tuples of nonnegative integers summing to c.
function IntegerVectorsOfSum(c, r)
    if r eq 1 then
        return [ [c] ];
    end if;
    res := [];
    for k in [0..c] do
        subs := IntegerVectorsOfSum(c - k, r - 1);
        for v in subs do
            Append(~res, [k] cat v);
        end for;
    end for;
    return res;
end function;

// ----------------------------------------------------------------------
// DFS helper, top-level: we pass patterns and seen BY REFERENCE.
// ----------------------------------------------------------------------
//   w_list, mult_w, factors, N : data describing current partition
//   idx                         : which weight index we're at
//   accum                       : list of <w, alloc> so far
//   ~patterns                   : global list of patterns (Assoc arrays)
//   ~seen                       : AssociativeArray of ids we’ve already added
procedure DFS(w_list, mult_w, factors, N, idx, accum, ~patterns, ~seen)
    if idx gt #w_list then
        // Build the pattern D from accum.
        D := AssociativeArray();
        total := 0;
        for t in accum do
            w := t[1];
            alloc := t[2];
            pairs := factors[w];
            for k in [1..#alloc] do
                cnt := alloc[k];
                if cnt gt 0 then
                    pair := pairs[k];      // <d, m>
                    key := pair;
                    if IsDefined(D, key) then
                        D[key] +:= cnt;
                    else
                        D[key] := cnt;
                    end if;
                    total +:= cnt * key[1] * key[2];
                end if;
            end for;
        end for;

        // Sanity check: total degree should be N.
        if total ne N then
            return;
        end if;

        id := PatternId(D);
        if IsDefined(seen, id) then
            return;
        end if;

        seen[id] := true;
        Append(~patterns, D);
        return;
    end if;

    w := w_list[idx];
    c := mult_w[w];
    pairs := factors[w];
    vecs := IntegerVectorsOfSum(c, #pairs);

    for alloc in vecs do
        DFS(w_list, mult_w, factors, N, idx + 1, accum cat [ <w, alloc> ],
            ~patterns, ~seen);
    end for;
end procedure;

// ----------------------------------------------------------------------
// Main pattern generator. Returns:
//   patterns  : sequence of patterns D (Assoc from <d,m> to count)
//   indexMap  : AssociativeArray from PatternId(D) to index.
// ----------------------------------------------------------------------
function GeneratePatterns(N)
    weights := [1..N];

    // Precompute factor pairs for each weight.
    factors := AssociativeArray(Integers());
    for w in weights do
        factors[w] := FactorPairs(w);
    end for;

    patterns := [];
    seen := AssociativeArray(Strings());

    parts := PartitionsOf(N, N);
    for p in parts do
        // multiplicity of each weight w in the partition p
        mult_w := AssociativeArray(Integers());
        for w in p do
            if IsDefined(mult_w, w) then
                mult_w[w] +:= 1;
            else
                mult_w[w] := 1;
            end if;
        end for;

        w_list := SetToSequence(Keys(mult_w));
        // Default Sort is fine here (integers).
        Sort(~w_list);

        accum := [];
        DFS(w_list, mult_w, factors, N, 1, accum, ~patterns, ~seen);
    end for;

    // Sort patterns in a deterministic order using their string ids.
    ids := [ PatternId(D) : D in patterns ];
    perm := [1..#patterns];

    Sort(~perm, func< i, j |
        ids[i] lt ids[j] select -1 else
        ids[i] gt ids[j] select 1 else
        0
    >);

    patterns := [ patterns[i] : i in perm ];

    // Build index map (PatternId -> index).
    indexMap := AssociativeArray(Strings());
    for i in [1..#patterns] do
        indexMap[ PatternId(patterns[i]) ] := i;
    end for;

    return patterns, indexMap;
end function;

// Global data for N = 8, mirroring patterns.sage
Patterns, PatternIndex := GeneratePatterns(8);

// Public helper: given a pattern D, return its index.
function IndexOfPattern(D)
    return PatternIndex[ PatternId(D) ];
end function;
