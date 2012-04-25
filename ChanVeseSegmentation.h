/*
This software contains the C++ implementation of the "branch-and-mincut" framework for image segmentation
with various high-level priors as described in the paper:

V. Lempitsky, A. Blake, C. Rother. Image Segmentation by Branch-and-Mincut. 
In proceedings of European Conference on Computer Vision (ECCV), October 2008.

The software contains the core algorithm and an example of its application (globally-optimal 
segmentations under Chan-Vese functional).

Implemented by Victor Lempitsky, 2008
*/

#ifndef CHAN_VESE_SEGMENTATION_H
#define CHAN_VESE_SEGMENTATION_H

#include "BranchAndMincut.h"

//EXAMPLE(section 5.1 of the paper): segmentation under discretized Chan-Vese functional.
//Original continuous functional was suggested in:
//T. Chan, L. Vese: Active contours without edges. Trans. Image Process., 10(2), 2001.


class ChanVeseBranch: public Branch
{
public:
	static gtype mu; //bias
	static gtype lambda; //smoothness
	//the branch is defined by minimal and maximal bounds on the parameters c_b and c_f (corresponding to the average intensities of the foreground and the background)
	int minb;
	int maxb;
	int minf;
	int maxf;
	int mean;

	virtual bool IsLeaf()
	{
		if(minf >= maxf && minb >= maxb) 
			return true; //we are at the finest level of discretization
		return false;
	}
	
	virtual void BranchFurther(Branch **br1, Branch **br2); //see cpp file

	virtual void Clone(Branch **br_)
	{
		*br_ = new ChanVeseBranch;
		ChanVeseBranch *br = (ChanVeseBranch *)*br_;
		br->bound = bound;
		br->mean = mean;
		br->minb = minb;
		br->maxb = maxb;
		br->minf = minf;
		br->maxf = maxf;
	}

	virtual gtype GetConstant()
	{
		if(!mu && minb > maxf) return INFTY;  //with mu=0 the energy becomes symmetric with respect c_f <-> c_b. This line add a constraint c_b <= c_f.
		return 0;
	}
	
	virtual void GetUnaries(gtype *bgUnaries, gtype *fgUnaries); //see cpp file

	virtual bool prune() {
		if (maxb < mean && maxf < mean) {
			return true;
		} else if (minb > mean && minf > mean) {
			return true;
		}
		return false;
	}
};





#endif