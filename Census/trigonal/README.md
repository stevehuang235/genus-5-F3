This collection of code computes all isomorphism classes of trigonal curves of genus 5 over $$\mathbb{F}_{3}$$ and stores all the data, including the ones used in the intermediate steps. Our computation crucially relies on a well-known result that a curve of genus 5 is trigonal if and only if it can be represented as a plane quintic with exactly one singularity of delta invariant one. 

```enum_trigonal.sage``` - SageMath code that enumerates candidate curves. Our computation is divided into three separate cases depending on the type of singularity: a split node, a non-split node, and a potential cusp. In each case, we use the functions in ```find_orbits.sage``` to compute the orbit representatives of quintics under the action of the subgroup of $$\operatorname{PGL}_3(\mathbf{F}_3)$$ that preserves the type of singularity in question. 

```curve_check.m``` - Magma code that checks whether each of the candidate curves generated in the previous step is of genus 5 and computes their point count over $$\mathbf{F}_{3^i}$$ for i=1,2,3,4,5. The filtered data are then stored in ./data_unfiltered_genus/. 

After running ```curve_check.m```, we already have the complete list of trigonal curves of genus 5 over $$\mathbf{F}_3$$ with one representative for each isomorphism class. As an extra consistency check, and to get the weighted stacky point-count for this stratum, we have the following script: 
```isom_class_check.m``` - Magma code that computes the distinct isomorphism classes of curves in ./data_unfiltered_genus_updated/ and stores each isomorphism class in ./data_filtered. For each isomorphism class, it also computes the order of its automorphism group over $$\mathbf{F}_3$$. 

The final data is stored in the file ```data_trigonal.txt```. In each line of this table, we store the point counts of the curve over $$\mathbf{F}_{3^i}$$ for i=1,2,3,4,5, the label of its $$\mathbf{F}_3$$-automorphism group, and the defining equation in the form of a plane quintic. 

We obtain the following (nonstacky, stacky) count for the number of isomorphism classes of trigonal curves of genus 5 over $$\mathbf{F}_3$$: (230327, 229636). In particular, the weighted stacky point-count for this stratum agrees with the theoretical result of Wennink (https://arxiv.org/abs/1701.00375). 
