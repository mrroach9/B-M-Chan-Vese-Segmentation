#Fast Branch-and-Mincut based Chan-Vese Segmentation

COMP790-098 Course Project

### Introduction
This program derives from the paper [1] by Lempitsky et. al. that using Branch-and-Mincut 
framework on segmenting images using Chan-Vese model, by discretizing the foreground and
background mean colors c_f and c_b into 255 levels and do tree-search on each quad-interval.
The maxflow algorithm used in calculting graph cuts in this program comes from Boykov and 
Kolmogorov's paper[2], and Kohli and Torr's paper[3].

### Our Work

The aim of this project is to accelarate this algorithm from several aspects. 

1) Mean Pruning: by calculating the mean color c_m of the entire image we see that 

              min{c_f, c_b} <= c_m <= max{c_f,c_b}.
              
   This can be used on pruning at least one half of the branches at the start of searching.
Another point is that since background and foreground are actually symmetric under the 
energy function of C-V model, we can add another constrint that c_b < c_f, therefore we have

              c_b < c_m < c_f
              
   This constraint furthur prunes the param space into 1/4 of the original space.

2) Thumbsnail： The algorithm running time increases dramatically with the increasing of pixels since the
   Preflow Pushing Algorithm runs in super-quadratic time. However for a thumbsnail of the 
   original image, the Chan-Vese segmentation result will not change much (c_f and c_b). Therefore
   we can run the Branch-and-Mincut on the thumbsnail with the same energy function (in practice we 
   reduce the smooth term lambda to half so that the smoothness is invariant under scaling),
   then segment the original image in a pretty small range of (c_b, c_f).

###Future Works

1) The network flow algorithm can be further improved by using priority queue instead of normal
FIFO queue, picking the highset labeled node every time when augmenting, yielding an O(V^2* sqrt(E))
run time instead of the current O(V^2*E) time.

2) For extremely large image, the thumbsnail method has a dilemma between large thumbsnail and small
one: large thumbsnail will give better estimation but cost too long time before doing the real segmentation,
too small thumbsnail will give bad estimation such that the B&M has either to search in a large 
space, or giving local optimum. Therefore we can apply an image-pyramid method on these pictures:
build different levels of images, start searching from the smallest, and apply the estimation on the
second smallest, then apply the second estimation on the third one, and so on. It is believed that
in practice, for even huge image like 4000*4000, 3~4 levels will be sufficient (5%, 20%, 40%, 100%).
  
###Reference

[1] Lempitsky, V.S., Blake, A., Rother, C. Image segmentation by branch-and-mincut. In: 
    ECCV (4), pp. 15–29 (2008)
    
[2] Boykov, Y. and Kolmogorov, V. An Experimental Comparison of Min-Cut/Max-Flow 
    Algorithms for Energy Minimization in Vision. PAMI 26(9) (2004)
    
[3] Kohli, P., Torr, P.: Effciently Solving Dynamic Markov Random Fields Using Graph Cuts.
    ICCV 2005 (2005)
