/*
This software contains the C++ implementation of the "branch-and-mincut" framework for image segmentation
with various high-level priors as described in the paper:

V. Lempitsky, A. Blake, C. Rother. Image Segmentation by Branch-and-Mincut. 
In proceedings of European Conference on Computer Vision (ECCV), October 2008.

The software contains the core algorithm and an example of its application (globally-optimal 
segmentations under Chan-Vese functional).

Implemented by Victor Lempitsky, 2008
*/

#ifndef BRANCH_AND_MINCUT_H
#define BRANCH_AND_MINCUT_H

#include "maxflow\graph.h"

typedef int gtype; //working type, can be int, double or integer
const gtype INFTY = 1 << 29; //a large value
const gtype EPSILON = 1; //a small value
typedef Graph<gtype,gtype,gtype> GraphT;

extern int imWidth, imHeight; //global variables corresponding to image dimensions


//main class, implements a branch, i.e. a node in the tree
class Branch
{
public:
	gtype bound;

	virtual bool IsLeaf() = 0; //needs to be defined. Should return true if the node
	
	virtual void BranchFurther(Branch **br1_, Branch **br2_) = 0; //needs to be defined. Should return two branches of teh same type corresponding to the subtrees

	virtual void Clone(Branch **br) = 0; //needs to be defined. Should create a copy of the branch
	
	virtual gtype GetConstant() {return 0; } //can be redefined. Should return a constant term C_w for the branch

	virtual bool SkipEvaluation() { return false; } //can be redefined. If returns true, the bound is not evaluated and is assumed -infinity
	
	virtual void GetUnaries(gtype *bgUnaries, gtype *fgUnaries) = 0; //needs to be defined. Should fill in the arrays of aggregated unary potentials
																	//for the background and for the foreground
	virtual bool prune() { return false; }
};

//these two functions should be called before and after all procedures (or in the case if the image size changes)
void PrepareGraph(int imwidth, int imheight);
void ReleaseGraph();

//main function

Branch * //output - globally optimal branch(node) of the tree. Should be deleted afterwards.
	BranchAndMincut(int imwidth, int imheight, //input image sizes (should be the same as in the call to PrepareGraphs
					  Branch *root, //root branch
					  int *segmentation,  //output: globally optimal segmentation. For each pixel either 1(foreground) or 0(background).
					  bool bestFirst, Branch *initialGuess, //branch-and-bound variations bestFirst/depthFirst, initialGuess needed only if the bestFirst=false
					  gtype *pairwise, //pairwise terms. For each pixel (including boundary) - 4 edge-strength values: top-right, right, bottom-right, bottom. Edges going outside the grid are simply ignored.
					  gtype *commonUnaries,//foreground unaries independent on the branch. For each pixel - a value.
					  int *nCalls //output: number of calls to the lower bound evaluation (including leaf branch-nodes)
					  ); 

#endif