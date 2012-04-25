Chan-Vese Segmentation using Branch-and-Mincut
==========================

Branch-and-Mincut Chan-Vese Segmentation COMP790-098 Course Project

This program derives from the paper [1] by Lempitsky et. al. that using Branch-and-Mincut 
framework on segmenting images using Chan-Vese model, by discretizing the foreground and
background mean colors $c_f$ and $c_b$ into 255 levels and do tree-search on each quad-interval.
The maxflow algorithm used in calculting graph cuts in this program comes from Boykov and 
Kolmogorov's paper[2], and Kohli and Torr's paper[3].

The aim of this project is to accelarate this algorithm from several aspects. 

1) By calculating the mean color c_m of the entire image we see that 
              $min{c_f, c_b} <= c_m <= max{c_f,c_b}$.
   This can be used on pruning at least one half of the branches at the start of searching.

2) The maxflow algorithm currently uses preflow pushing algorithm using LIFO queues. However
   if we pick active nodes according to their labeling (i.e. DIST to terminal), and maintain
   a priority queue to achieve that, the time complexity can be reduced from O(V^3) to 
   O(V^{2.5}) in this problem.
   
   
   
=====Reference=====

[1] Lempitsky, V.S., Blake, A., Rother, C. Image segmentation by branch-and-mincut. In: 
    ECCV (4), pp. 15–29 (2008)
    
[2] Boykov, Y. and Kolmogorov, V. An Experimental Comparison of Min-Cut/Max-Flow 
    Algorithms for Energy Minimization in Vision. PAMI 26(9) (2004)
    
[3] Kohli, P., Torr, P.: Effciently Solving Dynamic Markov Random Fields Using Graph Cuts.
    ICCV 2005 (2005)
