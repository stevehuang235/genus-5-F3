## A Census of Curves of Genus 5 over $\mathbb{F}_3$

This repository contains the code (written in SageMath and Magma) and data used to build a census of (smooth, projective, and geometrically irreducible) curves of genus 5 over the finite field $\mathbb{F}_3$.

We conduct our census according to the classical classification of genus 5 curves, whereby a curve is of genus 5 if and only if it is one of the following: 
- a hyperelliptic curve (c.f. `/Census/hyperelliptic`),
- a trigonal curve (c.f. `/Census/trigonal`), or
- an intersection of three quadrics in $$\mathbb{P}^4$$ (c.f. `/Census/generic`).

The data in the hyperelliptic stratum was computed by Everett Howe using the algorithm described in https://arxiv.org/abs/2401.15255. The table of all hyperelliptic curves of genus 5 over $\mathbb{F}_3$ can be found in `/Census/hyperelliptic/data_hyperelliptic.txt`.

Our tabulation of trigonal curves of genus 5 over $\mathbb{F}_3$ is also complete, and a complete list of distinct isomorphism classes can be found in `/Census/trigonal/data_trigonal.txt`. 

Our tabulation of generic curves is in progess. 