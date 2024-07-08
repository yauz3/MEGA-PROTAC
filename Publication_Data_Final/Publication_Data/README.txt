There are three data directories: Input_data, Prerefinement_Publication_data and Refinement_Publication_data.

Input_data contains 22 input structures for the BOTCP (together with native structures).

Prerefinement_Publication_data contains 22 complexes before refinement. 
Every complex contains 100 best dockq structures from 100 best clusters (ranked by Vina Energy).
Clarification: we find the RRTs with the best vina energy in each cluster and ranked clusters according to it.
Afterwards, out of the first-ranked 100 clusters we save the best dockq structure for each of the clusters. 
There is also a csv file "before_refinement_dataframe" which contains various scores for each of the structures.

Refinement_Publication_data contains 22 complexes after refinement.
Every complex contains different number of structures according to the number of cluster that specific complex has.
From every cluster we save the best dockq structure.
There is also a csv file "refinement_dataframe" which stores various features regarding that specific structure.

**Reference**
[1] Rao, Arjun, et al. "Bayesian optimization for ternary complex prediction (BOTCP)." Artificial Intelligence in the Life Sciences 3 (2023): 100072.
