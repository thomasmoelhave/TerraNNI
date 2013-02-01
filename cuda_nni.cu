// -*- mode: c++; tab-width: 4; indent-tabs-mode: t; eval: (progn (c-set-style "stroustrup") (c-set-offset 'innamespace 0)); -*-
// vi:set ts=4 sts=4 sw=4 noet :

/*
 * TerraNNI: Natural neighbor interpolation in 2D and 3D for large datasets
 * Copyright (C) 2010, 2011: Pankaj K. Agarwal, Alex Beutel & Thomas Mølhave
 *
 * This file is part of TerraNNI.
 *
 * TerraNNI is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TerraNNI is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with TerraNNI.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "nni.h"
#include "cuda_nni.h"

using namespace std;
//using namespace nni;
using namespace cuda_nni;

#include <iostream>
#include <cassert>

#include <cmath>
#include <map>
#include <algorithm>
#include <string.h>

#include </usr/include/GL/glew.h>

//#include <helper_cuda.h>
//#include <helper_cuda_gl.h>

//#include "cutil_inline.h"
//#include <cutil_gl_inline.h>

//#include <cudaGL.h>

#include <cuda_gl_interop.h>
#include <cuda_runtime.h>
//#include <cutil_gl_error.h>

#include <cstdlib>

#include <GL/glut.h>
#include <GL/glext.h>
//#include <GL/gl.h>
//#include <GL/glu.h>
#include <GL/glx.h>


#define GL_ERROR2() CheckGLError2(__FILE__, __LINE__)
bool CheckGLError2(char* acSourceFile, int iLine)
{
	GLenum eErr;
	bool bError = false;
	eErr = glGetError();
	while (eErr != GL_NO_ERROR) {
		fprintf(stderr, "OpenGL: %s, errno %d, source file %s, source line %d\n", gluErrorString(eErr), eErr, acSourceFile, iLine);
		bError = true;
		eErr = glGetError();
	}
	return bError;
}

size_t checkMemory() {

	size_t theFree, theTotal;
	cudaError_t res = cudaMemGetInfo( &theFree, &theTotal );
	if(res != cudaSuccess)
		cerr << "cudaMemGetInfo Failed\n";
	cerr << "Memory: " << theFree << " :: " << theTotal << "\n";

	return theFree;

}

texture<float4, 2, cudaReadModeElementType> siteTex;
texture<float4, 2, cudaReadModeElementType> queryTex;

/* -------------------------------
Offscreen rendering
------------------------------- */

GLuint vboNum, vboDenom;
struct cudaGraphicsResource *vbo_res_num;
struct cudaGraphicsResource *vbo_res_denom;
bool vboInitialized = false;

struct cudaGraphicsResource* siteBuffer_CUDA;
struct cudaGraphicsResource* queryBuffer_CUDA;

struct cudaGraphicsResource* siteBuffer_CUDA_depth;

cudaGraphicsResource ** allSiteBuffers_Color;
cudaGraphicsResource ** allSiteBuffers_Depth;

// Initialize vertex buffer objects (VBOs) to store the numerator and
// denominator during NNI computation
void initVBO(int w_width, int w_height) {
	fprintf(stderr, "VBO Initialized: %d x %d\n", w_width, w_height);
	vboInitialized = true;

	unsigned int size = w_width * w_height * sizeof(int);	
	glGenBuffers(1, &vboNum);
	glBindBuffer(GL_ARRAY_BUFFER, vboNum);
	glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	cudaGraphicsGLRegisterBuffer(&vbo_res_num, vboNum, cudaGraphicsMapFlagsNone);

	glGenBuffers(1, &vboDenom);
	glBindBuffer(GL_ARRAY_BUFFER, vboDenom);
	glBufferData(GL_ARRAY_BUFFER, size, 0, GL_DYNAMIC_DRAW);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
	cudaGraphicsGLRegisterBuffer(&vbo_res_denom, vboDenom, cudaGraphicsMapFlagsNone);
}

int * rb_voronois_indices;
int rb_query_index;

void cuda_nni::cudaRegisterAllBuffers(int tC, int * rb_v, int rb_q) {

	rb_voronois_indices=rb_v; 
	rb_query_index=rb_q; 

	/*

	//unsigned int theFree, theTotal,theFree2;
	size_t theFree, theTotal,theFree2;
	cudaError_t res = cudaMemGetInfo( &theFree, &theTotal );

	//cudaGraphicsResource * allSiteBuffers_Color[tC];
	allSiteBuffers_Color = new cudaGraphicsResource*[tC];

	for(int i = 0; i < tC; i++) {
		if(true) continue;
		//if(i == 1) continue;

		fprintf(stderr, "Register Voronoi render buffer %d\n", rb_voronois[i]);

		res = cudaMemGetInfo( &theFree, &theTotal );
		if(res != CUDA_SUCCESS)
			fprintf(stderr, "cudaMemGetInfo FAILED\n");
		fprintf(stderr, "Memory: %u :: %u\n", theFree, theTotal);


		//cutilSafeCall(cudaGraphicsGLRegisterImage(&allSiteBuffers_Color[i], rb_voronois[i], GL_RENDERBUFFER, cudaGraphicsMapFlagsReadOnly));
		cudaError_t reg_res = cudaGraphicsGLRegisterImage(&allSiteBuffers_Color[i], rb_voronois[i], GL_RENDERBUFFER, cudaGraphicsMapFlagsReadOnly);
		if(reg_res != CUDA_SUCCESS) {
			cout << "Error registering image\n\n";
		}


		//cudaFreeArray(allSiteBuffers_Color[i]);
		//cutilSafeCall(cudaGraphicsGLUnregisterImage(allSiteBuffers_Color[i]));
		cudaError_t res2 = cudaGraphicsUnregisterResource(allSiteBuffers_Color[i]);
		if(res2 != cudaSuccess) {
			fprintf(stderr, "*******error unregistering resource*******\n\n");
		}

		res = cudaMemGetInfo( &theFree2, &theTotal );
		if(res != CUDA_SUCCESS)
			fprintf(stderr, "cudaMemGetInfo FAILED\n");
		fprintf(stderr, "Memory: %u :: %u\n", theFree2, theTotal);
		fprintf(stderr, "Memory used: %d\n\n", (theFree-theFree2));
	}

	fprintf(stderr, "Register Query render buffer %d\n", rb_query);
	res = cudaMemGetInfo( &theFree, &theTotal );
	if(res != CUDA_SUCCESS)
		fprintf(stderr, "cudaMemGetInfo FAILED\n");
	fprintf(stderr, "Memory: %u :: %u\n", theFree, theTotal);

	cutilSafeCall(cudaGraphicsGLRegisterImage(&queryBuffer_CUDA, rb_query, GL_RENDERBUFFER, cudaGraphicsMapFlagsReadOnly));

	res = cudaMemGetInfo( &theFree, &theTotal );
	if(res != CUDA_SUCCESS)
		fprintf(stderr, "cudaMemGetInfo FAILED\n");
	fprintf(stderr, "Memory: %u :: %u\n", theFree, theTotal);
     */



	cudaSwitchVoronoiRB(0);
	//cudaSwitchVoronoiRB(rb_voronois[0]);
}


// Switch to renderbuffer at index rbIndex
void cuda_nni::cudaSwitchVoronoiRB(int rbIndex) {
	//siteBuffer_CUDA = allSiteBuffers_Color[rbIndex];


	//Unregister previous renderbuffers if they exist
	cudaError_t res2 = cudaGraphicsUnregisterResource(siteBuffer_CUDA);
	if(res2 != cudaSuccess) {
		fprintf(stderr, "*******error unregistering resource*******\n\n");
	}

	res2 = cudaGraphicsUnregisterResource(queryBuffer_CUDA);
	if(res2 != cudaSuccess) {
		fprintf(stderr, "*******error unregistering query resource*******\n\n");
	}


	cout << "Switch to renderbuffer " << rb_voronois_indices[rbIndex] << "\n";

	size_t theFree, theFree2;

	theFree = checkMemory();

	cudaError_t reg_res = cudaGraphicsGLRegisterImage(&siteBuffer_CUDA, rb_voronois_indices[rbIndex], GL_RENDERBUFFER, cudaGraphicsMapFlagsReadOnly);
	if(reg_res != cudaSuccess) {
		cout << "Error registering image\n\n";
	}

	theFree2 = checkMemory();
	cerr << "Memory used: " << (theFree-theFree2) << "\n";



	fprintf(stderr, "Register Query render buffer %d\n", rb_query_index);
	theFree = checkMemory();
	cudaGraphicsGLRegisterImage(&queryBuffer_CUDA, rb_query_index, GL_RENDERBUFFER, cudaGraphicsMapFlagsReadOnly);
	theFree2 = checkMemory();
	cerr << "Memory used: " << (theFree-theFree2) << "\n";

}


void init() {
	glClearColor(0.0, 0.0, 0.0, 0.0);
	glEnable(GL_DEPTH_TEST);
}

void reshape(int w, int h) {
	glViewport(0,0, (GLsizei) w, (GLsizei) h);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(-5.0, 20.0, -5.0, 20.0, -50.0, 50.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

void fakeDisplayFunc() { }


// Start GLUT
void startWithGLUT() {
	int argc =0;	
	char * argv = "";

	glutInit ( &argc, &argv );
	glutInitDisplayMode ( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize ( 800, 800 );
	glutCreateWindow ( "TerraNNI" );
	init();
	glutDisplayFunc ( fakeDisplayFunc );
	glutReshapeFunc ( reshape );
	
	glewInit();

}

// Attempt to not create a window (without GLUT) for true offscreen rendering
// Has never quite worked
void startOffscreen2() {	
	Display *display = XOpenDisplay(0);
	XVisualInfo *vinfo;
	int attrList[20];
	int indx=0;
	GLXContext util_glctx;

//   Colormap cmap;
//   XSetWindowAttributes swa;
//   Window win;
//   XEvent event;

	if(!display) exit (1);
	
	attrList[indx] = GLX_USE_GL; indx++;
	attrList[indx] = GLX_RGBA; indx++;
	attrList[indx] = GLX_DEPTH_SIZE; indx++;
	attrList[indx] = 8; indx++;
	attrList[indx] = GLX_RGBA; indx++;
	attrList[indx] = GLX_RED_SIZE; indx++;
	attrList[indx] = 8; indx++;
	attrList[indx] = GLX_GREEN_SIZE; indx++;
	attrList[indx] = 8; indx++;
	attrList[indx] = GLX_BLUE_SIZE; indx++;
	attrList[indx] = 8; indx++;
	attrList[indx] = None;

	vinfo = glXChooseVisual(display, DefaultScreen(display), attrList);
	if (vinfo == NULL) {
		printf ("ERROR: Can't open window\n");
		exit (1);
	}

	Pixmap pixmap = XCreatePixmap(display, DefaultRootWindow(display), 1,1,vinfo->depth);
	GLXPixmap glxpix = glXCreateGLXPixmap(display, vinfo, pixmap);

	util_glctx = glXCreateContext(display, vinfo, NULL, False);
	if (util_glctx == NULL) {
		printf("glXCreateContext failed \n");
		return;
	}

	if (!glXMakeCurrent(display, glxpix, util_glctx)) {
		printf("glXMakeCurrent failed \n");
		exit(1);
	}
	
	glewInit();

	if(glXGetCurrentContext() == NULL) {
		fprintf(stderr, "Uh oh, no context is currently set...\n");
	}
	
}


/* -------------------------------
End offscreen rendering
------------------------------- */

/* -------------------------------
CUDA functions
------------------------------- */


// Set all values to zero in the array
// Used for zeroing the VBOs
__global__ void zeroValues(int1 * arrayPoint, int N) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
	if(i < N) arrayPoint[i].x = 0;
}




__device__ void cudaCheckBits(uint t, // value in query buffer (representing the bitwise-OR of rendered queries 
		int1 * denom,  
		int1 * numer, 
		uint nearestPointVal, // value in color buffer, which is the value (usually elevation of the nearest neighbor to the pixel
		int x_pixel, 
		int y_pixel, 
		int ncols, 
		int nrows, 
		int blockWidth, 
		int scaling, 
		int pixel_offset, 
		int curRound, 
		int spacingFactor) 
{

	// Only continue if both a query cone AND a input point covered the pixel
    if(t != 0 && nearestPointVal > 0) {
		int rBlockWidth = blockWidth * spacingFactor;
		int nrounds = spacingFactor*spacingFactor;

		// This code assumes we only have 32 bits in the color buffer
		// Loop through each bit in the depth buffer
		for(int ti = 0; ti < 32; ti++) {
			// If the bit is set to one then a query covered it
			if( (t >> ti) & 1 == 1) {
				
				// Find the query that was rendered
				int offset_i = (int)(scaling * (0.5 + ((ti*nrounds + curRound) % rBlockWidth))) + pixel_offset;
				int offset_j = (int)(scaling * (0.5 + ((ti*nrounds + curRound) / rBlockWidth))) + pixel_offset;

				int block_i = round((x_pixel - offset_i) / (float)(rBlockWidth * scaling));
				int block_j = round((y_pixel - offset_j) / (float)(rBlockWidth * scaling));
				
				int p_i = block_i * rBlockWidth + ((ti*nrounds + curRound) % rBlockWidth);
				int p_j = block_j * rBlockWidth + ((ti*nrounds + curRound) / rBlockWidth);
				
				int whole_index = p_j * ncols + p_i;
				
				// Increment the numerator and denominator for that query
				if(whole_index < ncols * nrows && whole_index >= 0) {
					atomicAdd(&(numer[whole_index].x),nearestPointVal );
					atomicAdd(&(denom[whole_index].x),1);
				}		
			}
		}
    }

}



// General function to perform Buffer Analysis
__global__ void calcNNI(int1 * numer, int1 * denom, int width, int height, int scaling, int blockWidth, int pixel_offset, int curRound, int spacingFactor) {
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;

	// (x,y) is a pixel location in the buffers
    // confirm that the pixel is in bounds	
	if(x < ((width*scaling)+2*pixel_offset) && y < ((height*scaling)+2*pixel_offset)) {

		// Read RGBA values from the two buffers
		float4 voronoi = tex2D(siteTex, x, y);
		float4 query = tex2D(queryTex, x, y);

		// Combine the RGBA values
		int vx = *((int *)(&(voronoi.x)));
		int vy = *((int *)(&(voronoi.y)));
		int vz = *((int *)(&(voronoi.z)));
		int vw = *((int *)(&(voronoi.w)));
		uint v = vx | (vy << 8) | (vz << 16) | (vw << 24); // do i need to cast to uint before the shfit -amb79

		int tx = *((int *)(&(query.x)));
		int ty = *((int *)(&(query.y)));
		int tz = *((int *)(&(query.z)));
		int tw = *((int *)(&(query.w)));
		uint t = tx | (ty << 8) | (tz << 16) | (tw << 24);  // do i need to cast to uint before the shfit -amb79

		// Use read values
		cudaCheckBits(t,denom, numer,v,x,y, width,height, blockWidth, scaling, pixel_offset, curRound, spacingFactor);
    }

}

// Given the numerator and denominator, compute the final values
__global__ void divide(int1 * numer, int1 * denom, int width, int height, int zMult) {
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;

	if(y < height && x < width) {	
		// If the query did not reach any colored pixels then set the result to NODATA
		if(denom[y*width+x].x == 0)	
			denom[y*width+x].x = (-9999 * zMult);
		else 
			denom[y*width + x].x = numer[y*width+x].x/denom[y*width+x].x;
	}
	
}


/* -------------------------------
End CUDA functions
------------------------------- */


int1 *dptrNum;
int1 *dptrDenom;


// Set up a run in CUDA
void cuda_nni::cudaSetup(int nrows, int ncols) {

	size_t num_bytes; 
	size_t num_bytes2; 

	cerr << "Map VBO resuorces\n";
	cudaGraphicsMapResources(1, &vbo_res_num, 0);
	cudaGraphicsMapResources(1, &vbo_res_denom, 0);
	cudaGraphicsResourceGetMappedPointer((void **)&dptrNum, &num_bytes, vbo_res_num);
	cudaGraphicsResourceGetMappedPointer((void **)&dptrDenom, &num_bytes2, vbo_res_denom);

	cudaGraphicsMapResources(1, &siteBuffer_CUDA, 0); 
	cudaArray * siteArray;
	cudaGraphicsSubResourceGetMappedArray(&siteArray, siteBuffer_CUDA, 0,0);
	cudaBindTextureToArray(siteTex, siteArray);

	cudaGraphicsMapResources(1, &queryBuffer_CUDA, 0); 
	cudaArray * queryArray;
	cudaGraphicsSubResourceGetMappedArray(&queryArray, queryBuffer_CUDA, 0,0);
	cudaBindTextureToArray(queryTex, queryArray);

	// Zero VBO values
	int N = nrows * ncols;
	int threadsPerBlock = 1000;
	int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
	fprintf(stderr, "Zero values\n");
	zeroValues<<<blocksPerGrid, threadsPerBlock>>>(dptrNum, N);
	zeroValues<<<blocksPerGrid, threadsPerBlock>>>(dptrDenom, N);
	cudaThreadSynchronize();

	cudaGraphicsUnmapResources(1, &queryBuffer_CUDA, 0); 
	cudaGraphicsUnmapResources(1, &siteBuffer_CUDA, 0); 
}

// Perform buffer analysis
__host__ void cudaBufferAnalysis(int nrows, int ncols, int scaling, int blockWidth, int pixel_offset, int curRound, int spacingFactor) {

	int p_height = scaling * nrows + 2 * pixel_offset;
	int p_width = scaling * ncols + 2 * pixel_offset;


	// Map renderbuffers for use
	cudaGraphicsMapResources(1, &siteBuffer_CUDA, 0); 
	cudaArray * siteArray;
	cudaGraphicsSubResourceGetMappedArray(&siteArray, siteBuffer_CUDA, 0,0);
	cudaBindTextureToArray(siteTex, siteArray);

	cudaGraphicsMapResources(1, &queryBuffer_CUDA, 0); 
	cudaArray * queryArray;
	cudaGraphicsSubResourceGetMappedArray(&queryArray, queryBuffer_CUDA, 0,0);
	cudaBindTextureToArray(queryTex, queryArray);

	//TODO: Do this value for nthreads need to be taken from CUDA?
	int nthreads = 32;
	dim3 blocks((p_width + nthreads - 1) / nthreads,(p_height + nthreads - 1) / nthreads);//(250,250);
	dim3 threads(nthreads,nthreads);//(20,20);
	calcNNI<<<blocks, threads>>>(dptrNum, dptrDenom, ncols, nrows, scaling, blockWidth, pixel_offset, curRound, spacingFactor);
	cudaThreadSynchronize();

	/*
	cutilSafeCall(cudaUnbindTexture(siteTex));
	cutilSafeCall(cudaUnbindTexture(queryTex));
	*/

	//cutilSafeCall(cudaGraphicsUnmapResources(1, &siteArray, 0)); 
	//cutilSafeCall(cudaGraphicsUnmapResources(1, &queryArray, 0)); 


	// Unmap renderbuffers
	// If this is not done then future OpenGL calls will produce unexpected results
	cudaGraphicsUnmapResources(1, &queryBuffer_CUDA, 0); 
	cudaGraphicsUnmapResources(1, &siteBuffer_CUDA, 0); 

}


// Complete buffer analysis
__host__ void cudaCompleteAnalysis(int nrows, int ncols, int ** res, int zMult) {

	//TODO: Do this value for nthreads need to be taken from CUDA?
	int nthreads = 32;

	// Divide the values in numerator and denominator
	dim3 blocks2((ncols + nthreads - 1) / nthreads,(nrows + nthreads - 1) / nthreads);
	dim3 threads2(nthreads,nthreads);
	divide<<<blocks2, threads2>>>(dptrNum, dptrDenom, ncols, nrows, zMult);
	cudaThreadSynchronize();

	// Unmap numerator VBO
	cudaGraphicsUnmapResources(1, &vbo_res_num, 0);
	cudaGraphicsUnregisterResource(vbo_res_num);
	glBindBuffer(GL_ARRAY_BUFFER_ARB, vboNum );

	// Unmap denominator VBO
	cudaGraphicsUnmapResources(1, &vbo_res_denom, 0);
	cudaGraphicsUnregisterResource(vbo_res_denom);
	glBindBuffer(GL_ARRAY_BUFFER_ARB, vboDenom );
	// Read values from denominator
	int * data = (int*) glMapBuffer(GL_ARRAY_BUFFER, GL_READ_ONLY);
	
	*res = data;

}

void cuda_nni::completeAnalysis(int nrows, int ncols, int ** res, int zMult) {
	cudaCompleteAnalysis(nrows,ncols,res,zMult);
}

void cuda_nni::cleanupCUDA() {
	cerr << "Clean up CUDA\n";

	glBindBuffer(GL_ARRAY_BUFFER, vboDenom );
	glUnmapBuffer(GL_ARRAY_BUFFER);

	glDeleteBuffers(1, &vboDenom);
	glDeleteBuffers(1, &vboNum);
}

void cuda_nni::bufferAnalysis(const nni::WorkPackage& wp, int scaling, int blockWidth, 
		int pixel_offset, int curRound, int spacingFactor) { 

	cudaBufferAnalysis(wp.nrows, wp.ncols, scaling, blockWidth, pixel_offset,curRound,spacingFactor);
}

// Initialize graphics card environment and CUDA
void cuda_nni::init(const nni::Settings& s) {	
	
	// Setup OpenGL with GLUT
	startWithGLUT();
	//startOffscreen2();



	//int deviceCount = 0;



	// Set CUDA
	//cerr << "Max gflops: "<< cutGetMaxGflopsDeviceId() << "\n";
	//cudaSetDevice(cutGetMaxGflopsDeviceId());
	//cudaGLSetGLDevice(cutGetMaxGflopsDeviceId());

	cudaSetDevice(0);
	cudaGLSetGLDevice(0);

	cudaDeviceProp deviceProps;
    //cudaGetDeviceProperties(&deviceProps, cutGetMaxGflopsDeviceId());
    cudaGetDeviceProperties(&deviceProps, 0);
	cerr << "CUDA device [" << deviceProps.name << "] has " << deviceProps.multiProcessorCount << "Multi-Processors\n";

	checkMemory();

}

void cuda_nni::setup_run(int nrows, int ncols) {
	cout << "CUDA setup run \n";
	initVBO(ncols, nrows);
}

