# Using McSplit to solve induced subgraph isomorphism

This is a slightly modified version of the sequential C++ program used in our IJCAI 2017 paper.
In this version, we backtrack as soon as the calculated upper bound is less than the order
of the pattern graph.  In addition, there is a small optimisation that stops splitting early if 
possible.
