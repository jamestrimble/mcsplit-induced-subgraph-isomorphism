# Using McSplit to solve induced subgraph isomorphism

This is a slightly modified version of the sequential C++ program used in our IJCAI 2017 paper.
In this version, we backtrack as soon as the calculated upper bound is less than the order
of the pattern graph.  In addition, there is a small optimisation
that stops splitting early if possible
(a250c4550680eda61ef2c78922a27e64e8715f90 and 23a9ec068db6f6ce80d0fd86d675808858580461);
this sometimes reduces the amount of work per search node.
