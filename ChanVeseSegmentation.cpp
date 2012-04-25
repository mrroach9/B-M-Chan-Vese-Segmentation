/*
This software contains the C++ implementation of the "branch-and-mincut" framework for image segmentation
with various high-level priors as described in the paper:

V. Lempitsky, A. Blake, C. Rother. Image Segmentation by Branch-and-Mincut. 
In proceedings of European Conference on Computer Vision (ECCV), October 2008.

The software contains the core algorithm and an example of its application (globally-optimal 
segmentations under Chan-Vese functional).

Implemented by Victor Lempitsky, 2008
*/


#include "ChanVeseSegmentation.h"
#include "image.h"
#include <time.h>

int *image;

gtype ChanVeseBranch::mu; //bias
gtype ChanVeseBranch::lambda; //smoothness

//splitting the branch
void ChanVeseBranch::BranchFurther(Branch **br1_, Branch **br2_)
{
	*br1_ = new ChanVeseBranch;
	*br2_ = new ChanVeseBranch;

	ChanVeseBranch *br1 = (ChanVeseBranch *)*br1_;
	ChanVeseBranch *br2 = (ChanVeseBranch *)*br2_;

	//splitting into halves either the range of c_f or the range of c_b
	br1->minf = minf;
	br1->minb = minb;
	br2->maxf = maxf;
	br2->maxb = maxb;
	br1->mean = mean;
	br2->mean = mean;

/*	if(maxf > mean && minf < mean) {
		br1->maxf = mean;
		br2->minf = mean + 1;
		br1->maxb = maxb;
		br2->minb = minb;
	} else if (maxb > mean && minb < mean) {
		br1->maxb = mean;
		br2->minb = mean + 1;
		br1->maxf = maxf;
		br2->minf = minf;
	}else */if(maxf-minf > maxb-minb) {
		br1->maxf = (maxf+minf)/2;
		br2->minf = br1->maxf+1;
		br1->maxb = maxb;
		br2->minb = minb;
	} else {
		br1->maxb = (maxb+minb)/2;
		br2->minb = br1->maxb+1;
		br1->maxf = maxf;
		br2->minf = minf;
	}
}



inline gtype dist2segment(gtype val, gtype minSegm, gtype maxSegm)
{
	if(val <= minSegm) return minSegm-val;
	if(val <= maxSegm) return 0;
	return val-maxSegm;
}

//computing aggregated unary potentials for each pixel
void ChanVeseBranch::GetUnaries(gtype *bgUnaries, gtype *fgUnaries)
{
	for(int i = 0; i < imWidth*imHeight; i++)
	{
		bgUnaries[i] = dist2segment(image[i], minb, maxb);
		bgUnaries[i] *= bgUnaries[i];
		fgUnaries[i] = dist2segment(image[i], minf, maxf);
		fgUnaries[i] *= fgUnaries[i];
	}
}



int main()
{
	double totalTime = -clock();
	const char *path = "google.png";

	int w,h;

//// IMAGE LOADING (your preferred image I/O tool can be used instead)
	printf("Loading image...");
	image = LoadImage8bpp<gtype>(path, w, h); //loading the image into the array (row-wise). w and h contain image dimensions after the call
	if(!image)
	{
		puts("Invalid path to the test image!");
		return -1;
	}
	printf("done.\n");
//////////////////////////////////

	printf("Preparing graph...");
	int *segm = new int[w*h];

	ChanVeseBranch::lambda = 10000; //smoothness in the Chan-Vese functional
	ChanVeseBranch::mu = 0; //bias in the Chan-Vese functional 

	gtype *unaries = new gtype[w*h]; //array for branch independent unary terms
	gtype *pairwise = new gtype[w*h*4]; //array for pairwise terms

	for(int i = 0; i < w*h; i++)
	{
		unaries[i] = gtype(ChanVeseBranch::mu);
		//creating contrast-independent (Euclidean-regularization) edge links
		pairwise[4*i] = pairwise[4*i+2] = gtype(ChanVeseBranch::lambda); //horizonta and vertical edges
		pairwise[4*i+1] = pairwise[4*i+3] = gtype(ChanVeseBranch::lambda/sqrt(2.0)); //diagonal edges
	}

	//initialiizing the root branch
	ChanVeseBranch root;
	root.minb = 0;
	root.maxb = 255;
	root.minf = 0;
	root.maxf = 255;
	root.mean = 86;

	//preparing the graph for mincuts
	PrepareGraph(w, h);
	printf("done.\n");

	printf("Start Branch-and-Mincut searching!\n");
	int nCalls;
	double time = -clock();
	ChanVeseBranch *resultLeaf = (ChanVeseBranch *)BranchAndMincut(w, h, &root, segm, true, NULL, pairwise, unaries, &nCalls); //main function call
	time += clock();
	time /= CLOCKS_PER_SEC;

	printf("PROPORTION = %lf, TIME = %lf\n", (ChanVeseBranch::mu? 2.0 : 1.0)*double((root.maxb-root.minb)*(root.maxf-root.minf))/nCalls, time);
	totalTime += clock();
	totalTime /= CLOCKS_PER_SEC;
	printf("Total Time = %lf.\n", totalTime);

	printf("c_b = %lf, c_f = %lf\n", double(resultLeaf->minb), double(resultLeaf->minf));


//// VISUALIZATION (your preferred image I/O tool can be used instead)
	double *imageColor = LoadImage24bpp<double>(path, w, h);
	DrawSegmentation24bpp<double,int>(imageColor, segm, w, h);
	ShowImage24bpp<double>(imageColor, w, h, 0, "result");
////////////////////////////////

	delete resultLeaf;
	delete[] imageColor;
	delete[] image;
	delete[] pairwise;
	delete[] unaries;
	delete[] segm;
	return 0;
}