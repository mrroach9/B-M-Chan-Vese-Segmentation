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
	br1->image = image;
	br2->image = image;
	if(maxf-minf > maxb-minb) {
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

double calcMean(int* image, int w, int h) {
	double total = 0;
	for (int i = 0;i < w; ++i){
		for (int j = 0; j < h; ++j){
			total += image[i*h+j];
		}
	}
	return total / (w*h);
}

ChanVeseBranch* runBranchAndMinCut(int* image, int w, int h, gtype lambda, gtype mu, int** segm) {
	double mean = calcMean(image,w,h);
	ChanVeseBranch root;
	root.minb = 0;
	root.maxb = (int)mean;
	root.minf = (int)mean + 1;
	root.maxf = 255;
	root.image = image;
	*segm = new int[w*h];

	ChanVeseBranch::lambda = lambda; //smoothness in the Chan-Vese functional
	ChanVeseBranch::mu = mu; //bias in the Chan-Vese functional 

	gtype *unaries = new gtype[w*h]; //array for branch independent unary terms
	gtype *pairwise = new gtype[w*h*4]; //array for pairwise terms

	for(int i = 0; i < w*h; i++)
	{
		unaries[i] = gtype(ChanVeseBranch::mu);
		//creating contrast-independent (Euclidean-regularization) edge links
		pairwise[4*i] = pairwise[4*i+2] = gtype(ChanVeseBranch::lambda); //horizonta and vertical edges
		pairwise[4*i+1] = pairwise[4*i+3] = gtype(ChanVeseBranch::lambda/sqrt(2.0)); //diagonal edges
	}
	PrepareGraph(w, h);
	int nCalls;
	ChanVeseBranch *resultLeaf = (ChanVeseBranch *)BranchAndMincut(
		w, h, &root, *segm, true, NULL, pairwise, unaries, &nCalls); //main function call
	delete[] pairwise;
	delete[] unaries;

	return resultLeaf;
}

ChanVeseBranch* thumbsnailEstimate(const char* path, gtype lambda, gtype mu, int** segm) {
	int* image;
	int w, h;
	image = LoadImage8bpp<gtype>(path, w, h); 
	if(!image)
	{
		puts("Invalid path to the test image!");
		return NULL;
	}
	return runBranchAndMinCut(image, w, h, lambda/2, mu, segm);
}

void visualize(const char* path, int* segm){
	int w,h;
	double *imageColor = LoadImage24bpp<double>(path, w, h);
	DrawSegmentation24bpp<double,int>(imageColor, segm, w, h);
	ShowImage24bpp<double>(imageColor, w, h, 0, "result");
	delete[] imageColor;
}

int main()
{
	double totalTime = -clock();
	const char *rPath = "test.png";
	const char *path = "lake_google_maps.png";
	int** segm = new (int*);
	*segm = NULL;

	ChanVeseBranch* resultLeaf = thumbsnailEstimate(rPath, 10000, 0, segm);
	totalTime += clock();
	totalTime /= CLOCKS_PER_SEC;
	printf("Total Time = %lf.\n", totalTime);

	printf("c_b = %lf, c_f = %lf\n", double(resultLeaf->minb), double(resultLeaf->minf));
	
	visualize(rPath, *segm);


	delete resultLeaf;
	delete[] *segm;
	delete	 segm;
	return 0;
}