# Add missing Charmm36 terms

- Curently most Amber terms are scales.
- CMAP has been added
- Additional terms to scale for charmm36m


# Allow intra molecular solute scaling

The gREST (generalized REST) framework introduces a new approach where the solute region is defined not only as a whole solute molecule but also as a part of a solute molecule. Crucially, l_i is defined as the maximum number of particles involved in the i-th solute-solvent interaction:
* 4 for dihedrals or improper torsions
\* 8 for CMAP.

This is the key parameter governing how boundary-crossing terms are scaled.

The general scaling rule for a bonded term crossing the boundary is:
scale(term) = λ^(k/l)
Where:
* λ is the scaling parameter (= T₀/T_solute)
* k = number of atoms in the term that are inside the solute region
* l = total number of atoms in the term (2, 3, or 4)

* Done for dihedrals


# Try to solve slow weight convergence

Exponential Moving Average (EMA) 

Instead of a simple cumulative average, use a recency-weighted average:

```python
g_i(t+1) = (1 - α) * g_i(t) + α * U_i(t)
```
-
* Recent (more native-like) energy samples dominate
* Old unfolded-state energies are progressively forgotten
* α ~ 0.01–0.05 is a good starting range
* No structural change to the algorithm, just replace the average

Gradually reduce α as statistics accumulate (1/t schedule)


# Implement REST3 algorithm ?

seems difficult as SST2 works by scaling epsilon solute parameters. 