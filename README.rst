======
Murphy
======
Generic tools to identify overlapping genomic regions
-----------------------------------------------------
Murphy provides a set of tools based on an implementation of the interval
tree algorithm based on red-black trees from Introduction to Algorithms by
Cormen, Leiserson, Rivest and Stein (2001) 2nd Edition, The MIT Press

The intervalTree class makes a nice self-balancing tree.  The alrogithm
may be more efficient placing items in the tree if items are randomly
selected for insertion, instead of in sort order.

The CLRS algorithm has been modified to report all overlapping nodes instead of
only the first node.  This involved adding a min value instead of only a max
value in order to speed up the search by checking if the subtree min,max
overlaps with the search interval.  The search routine is also recursive.

:Authors:
    Suzy M. Stiegelmeyer
:Version: 0.1
:Dedication: to Murphy with love
