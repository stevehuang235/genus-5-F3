This collection of code computes all isomorphism classes of generic curves of genus 5 over $$\mathbb{F}_{3}$$ and stores all the data, including the ones used in the intermediate steps. Our computation uses the fact that a generic curve of genus 5 arises as a complete intersection of three quadrics in $$\mathbf{P}^4$$. *We do not store the intermediate data files due to space limit.*

```2-tuples-f3.txt``` - representatives of $$\operatorname{PGL}_5(\mathbf{F}_3)$$-orbits of two-dimensional subspaces of the space of quadrics in $$\mathbf{P}^4$$, equivalently, $$PGL_5(\mathbf{F}_3)$$-orbits of $$Gr(2,15)$$.

```enum_generic.sage``` - SageMath code that enumerates candidate curves. For each two-dimensional space of quadrics that has a singular subscheme of dimension less than 2, we extend it to a three-dimensional space. 

```curve_check.m``` - Magma code that checks whether each of the candidate curves generated in the previous step is of genus 5 and computes their point count over $$\mathbf{F}_{3^i}$$ for i=1,2,3,4,5. 

For each isogeny class, we compute an isomorphism invariant described in Section 2.3 of our paper using `isom_invariant.m`. The script there uses the precomputations handled by `hashing.m`, `hyperplanes.m`, `pattern.m`, and `poly_pattern.m` to speed up the computation.

The final data is stored in the folder ```./data_generic```. In each line of this table, we store the point counts of the curve over $$\mathbf{F}_{3^i}$$ for i=1,2,3,4,5, the label of its $$\mathbf{F}_3$$-automorphism group, and the singular plane model. These files are computed using ```database_gen_generic.m```. 

We obtain the following (nonstacky, stacky) count for the number of isomorphism classes of trigonal curves of genus 5 over $$\mathbf{F}_3$$: (?, ?), computed using `moduli_count.py` by iterating through all the files in the folder `./data_generic`.
