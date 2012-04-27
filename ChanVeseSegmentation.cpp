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
#include <math.h>
#include <algorithm>
#include <time.h>
#include <stdlib.h>

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

ChanVeseBranch* runBranchAndMincut(int* image, int w, int h, gtype lambda, gtype mu, 
								   int** segm, ChanVeseBranch root) {
	int* segment = new int[w*h];

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
		w, h, &root, segment, true, NULL, pairwise, unaries, &nCalls); //main function call
	delete[] pairwise;
	delete[] unaries;

	if (segm == NULL) {
		delete[] segment;
	} else {
		*segm = segment;
	}
	return resultLeaf;
}

ChanVeseBranch* thumbsnailEstimate(const char* path, gtype lambda, gtype mu, int** segm = NULL) {
	int* image;
	int w, h;
	image = LoadImage8bpp<gtype>(path, w, h); 
	if(!image)
	{
		puts("Invalid path to the test image!");
		return NULL;
	}
	double mean = calcMean(image,w,h);
	ChanVeseBranch root;
	root.minb = 0;
	root.maxb = (int)mean;
	root.minf = (int)mean + 1;
	root.maxf = 255;
	root.image = image;
	return runBranchAndMincut(image, w, h, lambda/2, mu, segm, root);
}

bool calcSSD(int* image, int w, int h, int b, int f, int bound){
	int total = 0;
	for (int i = 0; i < w*h; ++i){
		if (abs(image[i] - f) > abs(image[i] - b)) {
			total += (image[i] - b) * (image[i] - b);
		} else {
			total += (image[i] - f) * (image[i] - f);
		}
		if (total >= bound) return false;
	}
	return true;
}

ChanVeseBranch* calcFeasibleRegion(int* image, int w, int h, int bound) {
	double mean = calcMean(image, w, h);
	ChanVeseBranch root;
	root.image = image;
	root.minb = -1;
	root.maxb = -1;
	root.minf = -1;
	root.maxf = -1;
	for (int b = 0; b <= (int)mean; ++b){
		if (root.minb != -1) break;
		for (int f = (int)mean + 1; f <= 255; ++f) {	
			if (calcSSD(image, w, h, b, f, bound)) {
				root.minb = b;
				break;
			}
		}
	}
	for (int b = (int)mean; b >= 0; --b){
		if (root.maxb != -1) break;
		for (int f = (int)mean + 1; f <= 255; ++f) {	
			if (calcSSD(image, w, h, b, f, bound)) {
				root.maxb = b;
				break;
			}
		}
	}
	for (int f = (int)mean + 1; f <= 255; ++f){
		if (root.minf != -1) break;
		for (int b = root.minb; b <= root.maxb; ++b) {		
			if (calcSSD(image, w, h, b, f, bound)) {
				root.minf = f;
				break;
			}
		}
	}
	for (int f = 255; f >= (int)mean + 1; --f){
		if (root.maxf != -1) break;
		for (int b = root.minb; b <= root.maxb; ++b) {		
			if (calcSSD(image, w, h, b, f, bound)) {
				root.maxf = f;
				break;
			}
		}
	}
	return &root;
}

ChanVeseBranch* origImageSeg(const char* path, int lambda, int mu, int** segm,
							 int est_cf, int est_cb){
	int* image;
	int w, h;
	image = LoadImage8bpp<gtype>(path, w, h); 
	if(!image)
	{
		puts("Invalid path to the test image!");
		return NULL;
	}
	ChanVeseBranch root;
	root.maxb = std::min(est_cb + 10, 255);
	root.minb = std::max(0, est_cb - 10);
	root.maxf = std::min(255, est_cf + 10);
	root.minf = std::max(0, est_cf - 10);
	root.image = image;
/*
	printf("Estimating lower bound...");

	ChanVeseBranch* resultLeaf = runBranchAndMincut(
								image, w, h, lambda, mu, segm, root);
	int bound = resultLeaf->bound;

	printf("done.\nLower bound = %d.\n\n", bound);
	printf("Calculating feasible region...");

	root = *calcFeasibleRegion(image, w, h, bound);

	printf("done.\n");
	printf("Feasible region: c_b in [%d, %d], c_f in [%d, %d].\n\n",
		root.minb, root.maxb, root.minf, root.maxf);
*/
	ChanVeseBranch* resultLeaf = runBranchAndMincut(image, w, h, lambda, mu, segm, root);
	return resultLeaf;
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
	const char *thumbPath = "lake_google_4.png";
	const char *origPath  = "lake_google_1.png";
	int lambda = 10000;
	int mu = 0;

	int** segm = new (int*);
	*segm = NULL;

	double totalTime = -clock();

	printf("Running thumbsnail estimator...");

	ChanVeseBranch* resultLeaf = thumbsnailEstimate(thumbPath, lambda, mu);

	printf("done.\n");

	delete[] *segm;
	int est_cf = resultLeaf->maxf;
	int est_cb = resultLeaf->maxb;
	delete resultLeaf;

	printf("Estimated c_b = %d, c_f = %d.\n\n", est_cb, est_cf);
	
	printf("Segmenting original image...");

	resultLeaf = origImageSeg(origPath, lambda, mu, segm, est_cf, est_cb);

	totalTime += clock();
	totalTime /= CLOCKS_PER_SEC;

	printf("done.\n");
	printf("Total Time = %lf.\n", totalTime);
	printf("Energy = %d, c_b = %d, c_f = %d\n", 
		resultLeaf->bound, resultLeaf->minb, resultLeaf->minf);
	
	visualize(origPath, *segm);

	delete resultLeaf;
	delete[] *segm;
	delete	 segm;
	return 0;
}