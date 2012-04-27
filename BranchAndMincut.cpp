/*
This software contains the C++ implementation of the "branch-and-mincut" framework for image segmentation
with various high-level priors as described in the paper:

V. Lempitsky, A. Blake, C. Rother. Image Segmentation by Branch-and-Mincut. 
In proceedings of European Conference on Computer Vision (ECCV), October 2008.

The software contains the core algorithm and an example of its application (globally-optimal 
segmentations under Chan-Vese functional).

Implemented by Victor Lempitsky, 2008
*/


#include "BranchAndMincut.h"
#include <stdio.h>
#include <time.h>
#include <float.h>

//using stl for the queue in the min
#include <queue>
#include <deque>
#include <vector>
#include <functional>


//////////////////////////////////


int imWidth = 0;
int imHeight = 0;
static int *bestSegm;
static Branch *bestBranch = NULL;
static gtype upperBound; //best leaf energy found so far

//////////////////////////////////////////////

static int statFlowCalls; //counting calls to lower bound/energy evaluations

//////////////////////////////////////////////


struct 
{
	GraphT *graph;
	gtype *bgUnaries;
	gtype *fgUnaries;
	bool maxflowWasCalled;

	void Reset(gtype *pairwise, gtype *commonUnaries)
	{
		maxflowWasCalled = false;
		graph->reset();
		graph->add_node(imWidth*imHeight);

		int x,y,i;

		if(commonUnaries)
			for(i = 0; i < imWidth*imHeight; i++)
			{
				if(commonUnaries[i] > 0)
					graph->add_tweights(i, commonUnaries[i], 0);
				else
					graph->add_tweights(i, 0, -commonUnaries[i]);
			}

		for(y = 0, i = 0; y < imHeight; y++)
			for(x = 0; x < imWidth; x++, i++)
			{
				if(y && x < imWidth-1)	graph->add_edge(i, i-imWidth+1, pairwise[i*4], pairwise[i*4]);
				if(x < imWidth-1)	graph->add_edge(i, i+1, pairwise[i*4+1], pairwise[i*4+1]);
				if(y < imHeight-1 && x < imWidth-1)	graph->add_edge(i, i+imWidth+1, pairwise[i*4+2], pairwise[i*4+2]);
				if(y < imHeight-1)	graph->add_edge(i, i+imWidth, pairwise[i*4+3], pairwise[i*4+3]);
			}
		memset(fgUnaries, 0, sizeof(gtype)*imWidth*imHeight);
		memset(bgUnaries, 0, sizeof(gtype)*imWidth*imHeight);
	}
} reusable;

void PrepareGraph(int imwidth, int imheight)
{
	imWidth = imwidth;
	imHeight = imheight;

	reusable.graph = new GraphT(imwidth*imheight, imwidth*imheight*4);
	reusable.bgUnaries = new gtype[imwidth*imheight];
	reusable.fgUnaries = new gtype[imwidth*imheight];
}

void ReleaseGraph()
{
	delete reusable.graph;
	delete reusable.bgUnaries;
	delete reusable.fgUnaries;
}
//////////////////////////////////////////

//STL stuff

struct BranchWrapper
{
	Branch *br;
	BranchWrapper(Branch *b): br(b) {}
};

using namespace std;
bool operator<(const BranchWrapper& a, const BranchWrapper& b)
{
	return a.br->bound < b.br->bound;
}
bool operator>(const BranchWrapper& a, const BranchWrapper& b)
{
	return a.br->bound > b.br->bound;
}
typedef std::priority_queue<BranchWrapper, vector<BranchWrapper>, greater<vector<BranchWrapper>::value_type>> FRONT_QUEUE;

FRONT_QUEUE frontQueue;

///////////////////////////////////////////////


///////////////////////////////////////////////////////

gtype *currentBgUnaries = NULL;
gtype *currentFgUnaries = NULL;

////////////////////////////////////////////

bool BestFirstSearch();
void DepthFirstSearch(Branch *br);
gtype EvaluateBound(Branch *br);

Branch *BranchAndMincut(int imwidth, int imheight, 
					  Branch *root, int *segmentation, 
					  bool bestFirst, Branch *initialGuess, 
					  gtype *pairwise, 
					  gtype *commonUnaries,
					  int *nCalls)
{
	assert(imwidth == imWidth && imheight == imHeight);
	clock_t start = clock();

	bestSegm = segmentation;
	bestBranch = NULL;

	statFlowCalls = 0;

	if(currentBgUnaries)
		delete currentBgUnaries;
	if(currentFgUnaries)
		delete currentFgUnaries;

	currentBgUnaries = new gtype[imWidth*imHeight];
	currentFgUnaries = new gtype[imWidth*imHeight];

	reusable.Reset(pairwise, commonUnaries);

	upperBound = INFTY;

	if(!bestFirst && initialGuess)
		upperBound = EvaluateBound(initialGuess);

	Branch *root_;
	root->Clone(&root_);

	EvaluateBound(root_);

	if(bestFirst)
	{
		frontQueue.push(BranchWrapper(root_));
		while(BestFirstSearch());
		while(!frontQueue.empty())
		{
			Branch *br = frontQueue.top().br;
			delete br;
			frontQueue.pop();
		}
	}
	else
		DepthFirstSearch(root_);


	delete currentBgUnaries;
	delete currentFgUnaries;
	currentBgUnaries = NULL;
	currentFgUnaries = NULL;

	bestBranch->bound = upperBound;

	if(nCalls)
		*nCalls = statFlowCalls;

//	printf("Time spent in Branch-And-Mincut is %lf sec\n", double(clock()-start)/CLOCKS_PER_SEC);
	return bestBranch;
}

////////////////////////////////////////////


gtype EvaluateBound(Branch *br)
{
	statFlowCalls++;
	int i, x, y, imsize = imWidth*imHeight;

	if(br->SkipEvaluation())
	{
		br->bound = -INFTY;
		return -INFTY;
	}

//working with the constant term
	gtype boundVal = 0;
	gtype constant = br->GetConstant();

	gtype flow_limit = upperBound-constant;
	if(flow_limit < 0)
	{
		br->bound = upperBound+EPSILON;
		return upperBound+EPSILON;
	}

//updating unary terms in the graph
	br->GetUnaries(currentBgUnaries, currentFgUnaries);
	
	for(y = 0, i = 0; y < imHeight; y++)
		for(x = 0; x < imWidth; x++,i++)
		{
			gtype unaryUpdateBg = currentBgUnaries[i]-reusable.bgUnaries[i];
			gtype unaryUpdateFg = currentFgUnaries[i]-reusable.fgUnaries[i];
			reusable.bgUnaries[i] = currentBgUnaries[i];
			reusable.fgUnaries[i] = currentFgUnaries[i];
			if(unaryUpdateBg || unaryUpdateFg)
			{
				reusable.graph->add_tweights(i, unaryUpdateFg, unaryUpdateBg);
				
				if(reusable.maxflowWasCalled)
					reusable.graph->mark_node(i);
			}
		}

//evaluating lower bound by pushing flow
	boundVal = reusable.graph->maxflow(reusable.maxflowWasCalled, NULL)+constant;
	reusable.maxflowWasCalled = true;
	br->bound = boundVal;
	
	if(br->IsLeaf() && boundVal < upperBound)
	{
		//the new candidate for a global minimum
		upperBound = boundVal;
		if(bestBranch)
			delete bestBranch;
		br->Clone(&bestBranch);
		for(i = 0; i < imsize; i++)
			bestSegm[i] = (int)reusable.graph->what_segment(i);
//		printf("Bound value = %lf\n", double(boundVal));
	}

	return boundVal;
}

//////////////////////////////////////////////

bool BestFirstSearch()
{
	Branch *br = frontQueue.top().br;
//	printf("%d\t%d\n", br->bound, frontQueue.size());
	frontQueue.pop();

	if(br->IsLeaf())
	{
//		printf("Minimum found!\n");
		return false;
	}
	
	Branch *br1, *br2;
	br->BranchFurther(&br1, &br2);
	delete br;

	EvaluateBound(br1);
	frontQueue.push(BranchWrapper(br1));

	EvaluateBound(br2);
	frontQueue.push(BranchWrapper(br2));

	return true;
}



void DepthFirstSearch(Branch *br)
{
	if(br->IsLeaf())
		return;

	Branch *br1, *br2;
	br->BranchFurther(&br1, &br2);

	delete br;
	
	EvaluateBound(br1);
	EvaluateBound(br2);

	if(br1->bound < br2->bound)
	{
		if(br1->bound < upperBound)
		{
			DepthFirstSearch(br1);
			if(br2->bound < upperBound)
				DepthFirstSearch(br2);
		}
	}
	else
	{
		if(br2->bound < upperBound)
		{
			DepthFirstSearch(br2);
			if(br1->bound < upperBound)
				DepthFirstSearch(br1);
		}
	}	
}


