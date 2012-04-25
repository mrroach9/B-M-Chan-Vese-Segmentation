/*
This software contains the C++ implementation of the "branch-and-mincut" framework for image segmentation
with various high-level priors as described in the paper:

V. Lempitsky, A. Blake, C. Rother. Image Segmentation by Branch-and-Mincut. 
In proceedings of European Conference on Computer Vision (ECCV), October 2008.

The software contains the core algorithm and an example of its application (globally-optimal 
segmentations under Chan-Vese functional).

Implemented by Victor Lempitsky, 2008
*/

#ifndef IMAGE_H
#define IMAGE_H

#include <cv.h>
#include <highgui.h>
#include <stdio.h>

//a number of simple (and inefficient) wrappers around Intel OPEN_CV library

template<class T> T* LoadImage8bpp(const char *filename, int& width, int& height)
{
	IplImage *im = cvLoadImage(filename, 0);
	if(!im) return NULL;
	
	width = im->width;
	height = im->height;
	T *image = new T[width*height];

	for(int y = 0, i = 0; y < height; y++)
		for(int x = 0; x < width; x++, i++)
			image[i] = (T)(((uchar*)(im->imageData + im->widthStep*y))[x]);

	cvReleaseImage(&im);
	return image;
}

template<class T> T* LoadImage24bpp(const char *filename, int& width, int& height)
{
	IplImage *im = cvLoadImage(filename, 1);
	if(!im) return NULL;
	
	width = im->width;
	height = im->height;
	T *image = new T[width*height*3];

	for(int y = 0, i = 0; y < height; y++)
		for(int x = 0; x < width; x++,i++)
			for(int c = 0; c < 3; c++)
				image[3*i+c] = (T)(((uchar*)(im->imageData + im->widthStep*y))[3*x+c]);

	cvReleaseImage(&im);
	return image;
}


template<class T> void ShowImage24bpp(T *image, int width, int height, int pause, const char *caption, const char *outFile = NULL)
{

	IplImage *out = cvCreateImage(cvSize(width, height), IPL_DEPTH_8U, 3);

	cvNamedWindow(caption, CV_WINDOW_AUTOSIZE);
	

	for(int y = 0, i =0; y < height; y++)
		for(int x = 0; x < width; x++, i++)
		{
			CvScalar s;
			s.val[0] = image[3*i];
			s.val[1] = image[3*i+1];
			s.val[2] = image[3*i+2];
			//s.val[0] = s.val[1] = s.val[2] = double(image[i]-minVal)/double(maxVal-minVal+0.0001)*255;
			cvSet2D(out, y, x, s);
		}


	if(outFile)
		cvSaveImage(outFile, out);

	cvShowImage(caption, out);
	cvWaitKey(pause);

	cvReleaseImage(&out);

}

template<class T, class U> void DrawSegmentation24bpp(T *im, U *mask, int w, int h)
{
	for(int y = 0, i = 0; y < h; y++)
		for(int x = 0; x < w; x++, i++)
		{
			if(y >= 0  && x >= 0 && x < w && y < h &&
				(y > 0 && mask[i] != mask[i-w] || y < h - 1 && mask[i] != mask[i+w] ||
				 x > 0 && mask[i] != mask[i-1] || x < w - 1 && mask[i] != mask[i+1] ))
			{
				im[3*i+0] = 0;
				im[3*i+1] = 0;
				im[3*i+2] = 255;
			}
		}
}


#endif