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

#include <iostream>
#include <cassert>
#include "tpie_config.h"
#include <tpie/stream.h>
#include <cmath>
#include <map>
#include <algorithm>
#include <string.h>

#include <tpie/sort.h>
#include <tpie/progress_indicator_arrow.h>
#include <tpie/mm_manager.h>
#include <tpie/cpu_timer.h>

#include <cstdlib>

#include </usr/include/GL/glew.h>
#include <GL/glut.h>
#include <GL/glext.h>
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glx.h>

using namespace std;
using namespace nni;

#define PI 3.14159265
#define CONE_ANGLE 45.0

// Amount of segmentation around a strip in a hyperboloid or cone
#define CONE_SEGMENTS 50  

// Number of strips in rendering a hyperboloid
#define HYPERBOLOID_STRIPS 5 

// Value to offset z of points so that all heights are positive
int zOffset = 0; 

// Multiplier for z-values offering higher precision in CUDA calculations
int zMultiplier = 1.0;

float site_radius = 5.0;
float query_radius;

// Total number of rows and columns in all computation
int total_nrows, total_ncols; 

// 0 means render planes, anything else means render cones
int drawing_algorithm = 0; 

// Number of pixels between two adjacent queries
int scaling = 5;  

// Distance between two adjacent queries
double query_spacing = 2.0; 

const int bitsPerColor = 8;
const int bitsPower = 255; 
const int bitsMod = 256; 
const int totalBits = 4 * bitsPerColor;
int blockWidth = 5; 

int * i_map;
int * j_map;

#define GL_ERROR() CheckGLError(__FILE__, __LINE__)
bool CheckGLError(std::string acSourceFile, int iLine)
{
	GLenum eErr;
	bool bError = false;
	eErr = glGetError();
	while (eErr != GL_NO_ERROR) {
		fprintf(stderr, "OpenGL: %s, errno %d, source file %s, source line %d\n", gluErrorString(eErr), eErr, acSourceFile.c_str(), iLine);
		bError = true;
		eErr = glGetError();
	}
	return bError;
}


/* -------------------------------
   Rendering functionality
   ------------------------------- */

#define QUERY_INDEX 99
#define SITE_INDEX 100


// Draw a hyperboloid with maximum radius r 
// where time is the final term in the distance function that has been set constant
void displayListDrawHyperboloid(double r, double time) {
	double t2 = pow(time,2);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	double degInc = 2.0 * M_PI / CONE_SEGMENTS;
	double height = sqrt(pow(r,2) + t2) - abs(time); 

	double sheight = height / HYPERBOLOID_STRIPS;

	double curHeight = abs(time);
	for(int j = 0; j < HYPERBOLOID_STRIPS; j++) {
		glBegin(GL_TRIANGLE_STRIP);
		double curDeg = 0;
		for(int i = 0; i <= CONE_SEGMENTS; i++) {	
			double th = curHeight;
			double tr = sqrt(pow(th,2)-t2); 
			float x1 =tr * cos(curDeg);
			float y1 =tr * sin(curDeg);
			glVertex3f(x1,y1, -1.0f * th);

			th = curHeight+sheight;
			tr = sqrt(pow(th,2)-t2); 
			x1 = tr * cos(curDeg);
			y1 = tr * sin(curDeg);
			glVertex3f(x1,y1, -1.0f * th);
			curDeg += degInc;
		}
		glEnd();
		curHeight += sheight;
	}
}

// Create a display list for a render hyperboloid given a specific radius and time offset
GLuint getDrawHyperboloid(double radius, float timeOffset) {
	GLuint displayList = glGenLists(1);
	glNewList(displayList, GL_COMPILE);
	displayListDrawHyperboloid(radius,timeOffset);
	glEndList();
	return displayList;
}

//Store display list indices for hyperboloids
// Key is the time offset
// Query and site cones only differ in radius
map<float, GLuint> * queryCones = new map<float,GLuint>();
map<float, GLuint> * siteCones = new map<float,GLuint>();


float plane_min_z =std::numeric_limits<float>::infinity();
float plane_max_z =-std::numeric_limits<float>::infinity();


// Render plane
void drawPlane(int coneIndex,  // determines radius used (query vs. site radius)
		const Point3D& site, // point being drawn 
		const WorkPackage& wp,
		double x_offset = 0, 
		double y_offset = 0, 
		float timeOffset = 0) 
{
	double r = site_radius;
	if(coneIndex == QUERY_INDEX) {
		r = query_radius;
	}

	r = r / sqrt(2.0); // make it inscribed in circular radius of influence

	// Half the width of the buffer in real space
	float halfwidth = ((wp.nrows-1) * query_spacing)/2.0; 

	// Location in the buffer in real space
	float px = site.p[0] - x_offset;
	float py = site.p[1] - y_offset;

	// Location int he buffer if the origin was at the center of the buffer
	// Also scaled such that the buffer ranges approximately from -1 to 1 in both x and y
	float zpx = (site.p[0] - x_offset + halfwidth)/(halfwidth);
	float zpy = (site.p[1] - y_offset + halfwidth)/(halfwidth);

	glTranslatef(0,0,0);

	float z_offset =  zpx*zpx+ zpy*zpy + timeOffset*timeOffset/(halfwidth*halfwidth);
	float a = -2.0 * zpx / (halfwidth);
	float b = -2.0 * zpy / (halfwidth);


	glBegin(GL_TRIANGLE_STRIP);

	// Draw the plane as a square composed of two triangles
	for(int i = -1; i <= 1; i += 2){
		for(int j = -1; j <= 1; j += 2){
			float x = px + i*r;
			float y = py + j*r;

			float zx = x + halfwidth;
			float zy = y + halfwidth;

			//Height of the plane at x,y given coordinate system centered in the buffer
			float z = a*zx + b*zy + z_offset;  


			// Keep track of range of planes drawn to make sure glOrtho was set correctly
			plane_min_z = std::min(plane_min_z,z);
			plane_max_z = std::max(plane_max_z,z);

			glVertex3f(x,y,-1.0*z);
		}
	}

	glEnd();

	glLoadIdentity();
}

// Render a hyperboloid
void drawCone(int coneIndex, //determines radius of influence used 
		const Point3D& site, 
		double x_offset = 0, 
		double y_offset = 0, 
		float timeOffset = 0) 
{

	map<float, GLuint>* m = siteCones;
	double r = site_radius;
	if(coneIndex == QUERY_INDEX) {
		m = queryCones;
		r = query_radius;
	}

	glTranslatef((site.p[0] - x_offset), (site.p[1] - y_offset),0);

	GLuint displayList;
	if(m->count(timeOffset) > 0)
		displayList = (*m)[timeOffset];
	else { 
		displayList = getDrawHyperboloid(r,timeOffset);
		m->insert(pair<float,GLuint>(timeOffset,displayList));
	}
	glCallList(displayList);
	glLoadIdentity();
}


/* -------------------------------
   End rendering functionality
   ------------------------------- */

/* -------------------------------
   Offscreen rendering
   ------------------------------- */


int window_width = 1;
int window_height = 1;

// Render buffer objects currently being used
GLuint fb, rb_voronoi = -1, rb_query = -1, rb_depth = -1;

// Arrays of names of renderbuffers used for the color and depth buffers when
// drawing the Voronoi diagrams for input points
GLuint * rbVColors;
GLuint * rbVDepths;

// Array linking time of Voronoi diagram rendered in each buffer in rbVColors and rbVDepths
// Initialized to -1 when not set
// This only works because intensity and thus time must be a positive integer
int * rbTimes;

// Index of current render buffer
int curRB = 0;

// Number of renderbuffers needed based on how many times we render for
int timeCount = 0;


int getMinIndex(int * a) {
	int min = 0;
	for(int i = 1; i < timeCount; i++){
		if(a[i] < a[min])
			min = i;
	}
	return min;
}

void deleteAllRenderbuffers(){
	glDeleteRenderbuffers(1, &rb_query);
	glDeleteRenderbuffers(timeCount, rbVColors);
	glDeleteRenderbuffers(timeCount, rbVDepths);
}


// If FBO generated is not complete, print out information and exit program
void fboError() {
	cerr << "FBO not complete on renderbuffer switch\n"
		<< "Color renderbuffer: " << rb_voronoi << "\n"
		<< "Depth renderbuffer: " << rb_depth << "\n"
		;
	exit(1);
}

// Switch to renderbuffer for Voronoi diagram at time
// Return true if the Voronoi diagram has been generated for the given time,
// and false otherwise
bool switchRenderbuffers(int time) {
	bool found = false;

	for(int i = 0; i < timeCount; i++){

		// If the Voronoi diagram has already been generated and is still in a render buffer, use that one
		if(time == rbTimes[i]) {
			found =true;
			curRB = i;
			cerr << "Switch to saved renderbuffer " << i << "\n";
			break;
		}
	}

	// If the Voronoi diagram hasn't been generated, use the oldest renderbuffer
	// This assume times are processed from least to greatest
	if(!found){ 
		curRB = getMinIndex(rbTimes);
		rbTimes[curRB] = time;
	}

	// Set current renderbuffer
	rb_voronoi = rbVColors[curRB];
	rb_depth = rbVDepths[curRB];

	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rb_voronoi);GL_ERROR();
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb_depth);GL_ERROR();

	GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(status != GL_FRAMEBUFFER_COMPLETE) {
		fboError();
	}

	glFinish();

	cuda_nni::cudaSwitchVoronoiRB(curRB);

	return found;
}


void nni::clearRenderbuffers(){
	cout << "Clear render buffers, tC: " << timeCount << "\n\n";
	for(int i = 0; i < timeCount; i++){
		rbTimes[i] = -1; 
	}

	//Go backwards so rb_voronoi and rb_depth are both set to 0
	for(int i = timeCount-1; i >= 0; i--) {

		rb_voronoi = rbVColors[i];
		rb_depth = rbVDepths[i];
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rb_voronoi);GL_ERROR();
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb_depth);GL_ERROR();

		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	
	}

}

void nni::setRenderbuffers(int w_width, int w_height, int tC) {

	// Don't remake renderbuffers unnecessarily
	if(w_width == window_width && w_height == window_height && tC == timeCount) 
		return;

	cerr << "Actually set renderbuffers\n";
	timeCount = tC;

	if(rb_voronoi != -1 && rb_depth != -1 && rb_query != -1) {
		cerr << "Delete old renderbuffers\n";
		deleteAllRenderbuffers();
	}

	// Create renderbuffer for color buffer for query step
	glGenRenderbuffers(1, &rb_query); GL_ERROR();
	glBindRenderbuffer(GL_RENDERBUFFER, rb_query);  GL_ERROR();
	glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA, w_width, w_height);  GL_ERROR();
	//glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rb_query);

	// Initialize arrays indexing the buffers
	rbTimes = new int[timeCount];
	rbVColors = new GLuint[timeCount];
	rbVDepths = new GLuint[timeCount];

	// Reset times on renderbuffers
	for(int i = 0; i < timeCount; i++){
		rbTimes[i] = -1; 
	}

	// Generate renderbuffers for depth and color at each time
	GLuint * bufs = new GLuint[timeCount*2];
	glGenRenderbuffers(timeCount*2, bufs);	GL_ERROR();

	for(int i = 0; i < timeCount; i++) {

		rbVColors[i] = bufs[2*i];
		rbVDepths[i] = bufs[2*i+1];

		glBindRenderbuffer(GL_RENDERBUFFER, rbVColors[i]);GL_ERROR();
		glRenderbufferStorage(GL_RENDERBUFFER, GL_RGBA, w_width, w_height);GL_ERROR();

		glBindRenderbuffer(GL_RENDERBUFFER, rbVDepths[i]);GL_ERROR();
		glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT32, w_width, w_height);GL_ERROR();

	}

	curRB = 0;
	rb_voronoi = rbVColors[0];
	rb_depth = rbVDepths[0];	

	// amb79 - Need to figure out why but apparently registering buffers up front makes it better for memory

	cuda_nni::cudaRegisterAllBuffers(tC,(int*)rbVColors, rb_query);

	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rb_voronoi);GL_ERROR();
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, rb_depth);GL_ERROR();


	GLenum status = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(status != GL_FRAMEBUFFER_COMPLETE) {
		fboError();
	}

	// Set initial OpenGL setup
	glEnable(GL_DEPTH_TEST);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();

	// glOrtho is updated later
	glOrtho(0.0, 2000.0, 0.0, 2000.0, -50.0, 50.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	//Set viewport to be the full size
	glViewport(0,0,w_width, w_height);

	window_width = w_width;
	window_height = w_height;

}

/* -------------------------------
   End offscreen rendering
   ------------------------------- */


// Break unsigned int index into four bytes, each to fill RGBA of color buffer
void setUIntColor(uint index) {

	int i1 = index%bitsMod;
	int i2 = (index >> bitsPerColor)%bitsMod;
	int i3 = (index >> (bitsPerColor*2))%bitsMod;
	int i4 = (index >> (bitsPerColor*3))%bitsMod;

	glColor4ub(i1,i2,i3,i4);

}


// Draw query cone for query i,j with color set to have a 1 only in position pos
void drawBitCone(int i, int j, int pos, const WorkPackage& wp, float timeOffset = 0) {
	uint index = (1 << pos);

	setUIntColor(index);

	Point3D pt(i*query_spacing,j*query_spacing,0);

	if(drawing_algorithm == 0)
		drawPlane(QUERY_INDEX, pt,wp,0,0, timeOffset);
	else 
		drawCone(QUERY_INDEX, pt,0,0, timeOffset);
}


// Make query blocks and render points appropriately
void makeBlocks(const WorkPackage& wp, 
		int curRound, // If spacingFactor > 1 then multiple rounds are required to perform all queries in each query block 
		int spacingFactor, // Multiplier for width of query blocks 
		float timeOffset = 0) 
{

	int rBlockWidth = spacingFactor * blockWidth;

	// Number of query blocks across and down
	int width = ceil(1.0 * wp.ncols / rBlockWidth);
	int height = ceil(1.0 * wp.nrows / rBlockWidth);

	int nrounds = spacingFactor * spacingFactor; // Total number of rounds required

	for(int i = 0; i < wp.ncols; i++) {
		for(int j = 0; j < wp.nrows; j++) {

			int block = width * floor(1.0 * j / rBlockWidth) + floor(1.0 * i / rBlockWidth);

			// Index within the given block
			int i_b = i % rBlockWidth;
			int j_b = j % rBlockWidth;
			int index = j_b * rBlockWidth + i_b;

			int whole_index = j * wp.ncols + i;
			if(index % nrounds == curRound){

				drawBitCone(i,j,(int)floor(1.0 * index / nrounds),wp, timeOffset);

				assert((whole_index < (wp.nrows * wp.ncols) && whole_index >= 0) && "WRITE ERROR FOR I_MAP");
				i_map[whole_index] = i;
				j_map[whole_index] = j;

			}
		}
	}
}


// Draw input point pnt from perspective of timeOffset
void addPoint(const WorkPackage& wp, const Point3D& pnt, float timeOffset){
	uint index = floor((pnt.p[2] + zOffset + 1) * zMultiplier);

	assert((index != 0) && "Error - Point height is 0");

	setUIntColor(index);

	if(drawing_algorithm == 0)
		drawPlane(SITE_INDEX, pnt, wp, wp.origin_x, wp.origin_y, timeOffset);
	else 
		drawCone(SITE_INDEX, pnt, wp.origin_x, wp.origin_y, timeOffset);
}


tpie::cpu_timer timeSiteDraw, timeQueryVoronoi, timeQueryAnalyze, timeQueryWrite;

void nni::printTimes() {
	fprintf(stdout, "Time Results:\n");
	fprintf(stdout, "\tDraw site voronoi : \t\t%f\n", timeSiteDraw.wall_time());
	fprintf(stdout, "\tDraw query voronoi: \t\t%f\n", timeQueryVoronoi.wall_time());
	fprintf(stdout, "\tAnalayze query res: \t\t%f\n", timeQueryAnalyze.wall_time());
	fprintf(stdout, "\tWrite query points: \t\t%f\n", timeQueryWrite.wall_time());
}


void nni::cleanup() {
	deleteAllRenderbuffers();
	glDeleteFramebuffers(1, &fb);
}

// Initialize NNI system
int nni::init(const Settings& s) {	

	// Initialize parameters
	query_spacing = s.cell_size;
	scaling = s.scale;
	site_radius = s.site_radius;
	query_radius = s.query_radius;
	drawing_algorithm = s.algorithm;
	zOffset = -1 * s.pointMin;
	zMultiplier = s.multiplier;
	total_nrows = s.nrows;
	total_ncols = s.ncols;

	// Initialize GLUT and CUDA
	cuda_nni::init(s);

	// Generate framebuffer object
	glGenFramebuffers(1, &fb);
	glBindFramebuffer(GL_FRAMEBUFFER, fb);


	//Should figure out maximum size texture
	GLint max_size;
	glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max_size);

	cout << "Max size from OpenGL: " << GL_MAX_RENDERBUFFER_SIZE << "\n";
	cout << "Max size from CUDA: " << max_size << "\n";

	return max_size;
}

// Set up NNI system for individual run on a work package
void nni::setupRun(const WorkPackage& wp, int pixel_offset, int timeCount) {
	cerr << "Setup run for new work package\n"
		<< "\tOrigin: " << wp.origin_x << ", " << wp.origin_y << "\n"
		<< "\tRows, cols: " << wp.nrows << ", " << wp.ncols << "\n"
		<< "\tPixel offset: " << pixel_offset << "\n"
		<< "\tScaling: " << scaling << "\n"
		;

	int x_offset = wp.origin_x;
	int y_offset = wp.origin_y;

	double pixel_size = query_spacing / scaling;

	int x_pixels = wp.ncols * scaling + 2 * pixel_offset;
	int y_pixels = wp.nrows * scaling + 2 * pixel_offset;

	setRenderbuffers(x_pixels, y_pixels, timeCount);

	// Initialize vertex buffer objects for buffer analysis
	cuda_nni::setup_run(wp.nrows, wp.ncols);

	i_map = new int[(wp.ncols * wp.nrows)];
	j_map = new int[(wp.ncols * wp.nrows)];


	// Set up OpenGL system appropriately
	glDepthMask(true);
	glDepthFunc(GL_LEQUAL);
	glDisable(GL_COLOR_LOGIC_OP);
	glLogicOp(GL_COPY);

	glAlphaFunc(GL_ALWAYS, 1);
	glDisable(GL_ALPHA_TEST);


	glViewport(0,0, x_pixels, y_pixels);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();



	// Set up z-range for glOrtho
	// TODO: the z area should be calculated more dynamically for the case of planes
	float ortho_zmin = -200.0;
	float ortho_zmax = 200.0;
	if(drawing_algorithm != 0) {
		float height = site_radius / tan(CONE_ANGLE * PI / 180.0);
		ortho_zmin = -1.0;
		ortho_zmax = height+1;
	}

	// Set up glOrtho with border areas
	// Each border area includes half a query spacing length along with an
	// addition pixel_offset to cover points that fall near by
	GL_ERROR();
	glOrtho(	0 - (query_spacing / 2.0) - (pixel_offset * pixel_size), 
			query_spacing * (wp.ncols-1) + (query_spacing/2.0) + (pixel_offset * pixel_size), 
			0 - (query_spacing/2.0) - (pixel_offset * pixel_size), 
			query_spacing * (wp.nrows-1) + (query_spacing/2.0) + (pixel_offset * pixel_size), 
			ortho_zmin, ortho_zmax); GL_ERROR();


	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glPolygonMode(GL_FRONT, GL_FILL);

	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rb_voronoi);GL_ERROR();
}



// Output arbitrary color buffer as an ASCII grid file
// glT must match bufname and type paramters (see Color and Depth buffer examples below)
template <typename glT>
void outputBuffer(int x_pixels, int y_pixels, double cellsize, float nodata,
		GLenum bufname, GLenum type, std::string filename = "buffer.asc") {

	glT * buf;
	buf = new glT[x_pixels*y_pixels];
	glReadPixels(0,0, x_pixels,y_pixels, bufname, type, buf);

	ofstream outfile(filename.c_str());

	outfile 
		<< "NROWS " << x_pixels << "\n"
		<< "NCOLS " << y_pixels << "\n"
		<< "XLLCORNER " << 0 << "\n"
		<< "YLLCORNER " << 0 << "\n"
		<< "CELLSIZE " << cellsize << "\n"
		<< "NODATA_VALUE " << nodata << "\n";

	for(int i = 0; i < x_pixels; i++) {
		for(int j =0; j < y_pixels; j++) {
			//if(buf[i*x_pixels+j] != nodata)
			//cout << i << ", " << j << ", " << buf[i*x_pixels+j] << "\n";
			outfile << buf[i*x_pixels+j] << " ";
		}
		outfile << "\n";
	}

}


void outputColorBuffer(int x_pixels, int y_pixels, double cellsize,
		std::string filename = "color_buffer.asc") {

	outputBuffer<GLuint>(x_pixels,y_pixels,cellsize,0,GL_RGBA,GL_UNSIGNED_INT_8_8_8_8_REV,filename);

}

void outputDepthBuffer(int x_pixels, int y_pixels, double cellsize,
		std::string filename = "depth_buffer.asc") {

	outputBuffer<GLfloat>(x_pixels,y_pixels,cellsize,1,GL_DEPTH_COMPONENT,GL_FLOAT,filename);

}



// Render Voronoi diagram for time t from set of input points
// When done seek back to the appropriate position in the TPIE stream
void renderInputPoints(float t, //Time we are rendering for 
		const WorkPackage& wp, // Current workpackage
		tpie::ami::stream<nni::Point4D>& vec, //Input points (already at correct starting position) 
		float timeRadius, 
		bool useDelta, // if true use Delta algorithm, otherwise use Delta^2 algorithm 
		TPIE_OS_OFFSET& finalSeekTo, 
		bool& firstSeek, 
		bool isFinalSeek) 
{

	// Prepare OpenGL for rendering
	glDepthMask(true);
	glDepthFunc(GL_LEQUAL);
	glDisable(GL_COLOR_LOGIC_OP);
	glLogicOp(GL_COPY);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	


	int cnt = 0;
	TPIE_OS_OFFSET seekto = vec.tell();
	bool foundSeek = false;
	nni::Point4D * p_ptr;

	cerr 
		<< "Draw Voronoi at time " << t << "\n"
		<< "Start TPIE stream at " <<  vec.tell() << " out of " << vec.stream_len() << "\n"
		;

	// While points are available, read points
	while (vec.read_item(&p_ptr) == tpie::ami::NO_ERROR) {
		const nni::Point4D& p = *p_ptr;

		// Check if we need to set seek position based on the new point
		if(!foundSeek && p.t >= t + 1 - timeRadius){ 
			seekto = vec.tell()-1;   // Can I subtract one? - amb79 
			foundSeek = true;
			if(firstSeek){
				firstSeek = false;
				finalSeekTo = vec.tell()-1;
			}
		}

		//cout << p.t << "\n";

		// If the point is within the time radius of influence, render it
		// If not and the point is after the current time, we are done
		// rendering for this time, since the points are sorted, so
		// break loop
		if(abs(t-p.t) <= timeRadius + 0.000001) { /// UGLY HACK - amb79
			cnt++;
			addPoint(wp,p,abs(t - p.t));
		}
		else if(p.t > t + timeRadius)
			break;

	}

	cerr << cnt << " points rendered\n";

	// Seek in TPIE stream to appropriate point
	if(!useDelta && isFinalSeek) {
		vec.seek(finalSeekTo);
		cerr << "Final Seek to " << finalSeekTo << "\n";
	} else {
		vec.seek(seekto);
		cerr << "Seek to " << seekto << "\n";
	}

}


// Perform NNI queries - render query blocks and perform buffer analysis
void performNNIQueries(const nni::WorkPackage& wp, int scaling, int blockWidth, 
		int pixel_offset, int spacingFactor, int timeOffset){

	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_RENDERBUFFER, rb_query);GL_ERROR();

	cerr << "Spacing factor " << spacingFactor << "\n";
	for(int i = 0; i < (int)pow(spacingFactor,2); i++){

		//Set up buffer
		glDepthMask(false);
		glDepthFunc(GL_LEQUAL);
		glEnable(GL_COLOR_LOGIC_OP);
		glLogicOp(GL_OR);
		glClear(GL_COLOR_BUFFER_BIT);	

		// Render queries
		timeQueryVoronoi.start();
		makeBlocks(wp,i,spacingFactor,timeOffset);	
		glFinish();
		timeQueryVoronoi.stop();

		//Analyze buffers
		timeQueryAnalyze.start();
		cuda_nni::bufferAnalysis(wp, scaling, blockWidth,pixel_offset,i,spacingFactor);
		timeQueryAnalyze.stop();
	}
}


// Read results from CUDA VBO and add results to TPIE stream
void completeBufferAnalysis(const nni::WorkPackage& wp, tpie::ami::stream<nni::OutPoint>& out, 
		int * i_map, int * j_map) { 
	int * res;

	// Get points from CUDA
	timeQueryAnalyze.start();
	cuda_nni::completeAnalysis(wp.nrows, wp.ncols, &res,zMultiplier);
	timeQueryAnalyze.stop();

	timeQueryWrite.start();
	int cnt1 = 0;
	int cnt2 = 0;
	for(int i = 0; i < wp.ncols; i++) {
		for(int j = 0; j < wp.nrows; j++) {			
			int whole_index = j * wp.ncols + i;

			// Undo scaling 
			double z = 1.0 * res[whole_index] / zMultiplier; 
			if(z != -9999) {
				cnt1++;
				z = z - zOffset - 1;
			} else { 
				cnt2++; 
			}

			nni::OutPoint threeDPt(i_map[whole_index] + wp.grid_origin_x, j_map[whole_index] + wp.grid_origin_y, z);
			//nni::OutPoint threeDPt(i + wp.grid_origin_x, j + wp.grid_origin_y, z);

			assert( (i_map[whole_index] >= 0 && j_map[whole_index] >=0 && i_map[whole_index] < wp.ncols && j_map[whole_index] < wp.nrows) && "Error - output point out of bounds");

			//cout << "Error: " << i << ", " << j << ", " << i_map[whole_index] << ", " << j_map[whole_index] <<"\n";

			// Only enter point in TPIE stream if we desired a result
			if(threeDPt.i < total_ncols && threeDPt.j < total_nrows) 
				out.write_item(threeDPt);
		}
	}
	timeQueryWrite.stop();

	cerr << "NODATA for " << cnt2 << " of " << (cnt1+cnt2) << "\n";

	cuda_nni::cleanupCUDA();
}


// General function to perform NNI
void nni::performNNI(const WorkPackage& wp, tpie::ami::stream<OutPoint>& out,
		tpie::ami::stream<nni::Point4D>& vec, float timeRadius, float curTime,
		bool useDelta, bool debug) {


	// Derive initial values based on passed parameters

	// Set spacing factor such that queries don't overlap
	// If query radius was set automatically then spacingFactor will be 1
	int spacingFactor = (int)ceil((query_radius * 2 + 0.001) / (query_spacing * blockWidth)); 

	double pixel_size = query_spacing / scaling;

	int pixel_offset = ceil((query_radius-floor(scaling / 2.0)*pixel_size) / pixel_size) + 1;

	int x_pixels = wp.ncols * scaling + 2 * pixel_offset;
	int y_pixels = wp.nrows * scaling + 2 * pixel_offset;


	cerr 
		<< "Perform NNI for " << curTime << "\n"
		<< "Cols, rows: " << wp.ncols << ", " << wp.nrows << "\n"
		<< "Pixel offset: " << pixel_offset << "\n"
		<< "Query spacing: " << query_spacing << "\n"
		<< "Block width: " << blockWidth << "\n"
		<< "Spacing factor: " << spacingFactor << "\n"
		<< "Time radius: " << timeRadius << "\n";
	;


	// Number of renderbuffers required
	// TODO: Verify this is correct
	int numRB = useDelta ? (timeRadius*2+1) : 1;
	setupRun(wp,pixel_offset, numRB);
	cuda_nni::cudaSetup(wp.nrows, wp.ncols);

	TPIE_OS_OFFSET finalSeekTo;
	bool firstSeek = true;

	// TODO: time resolution rather than just 1 (j++)
	for(int j = -1 * (int)timeRadius; j <= timeRadius; j++) {

		float t = curTime + j;

		// If useDelta is set to false then we are using the Delta^2 algorithm, switchRenderbuffers is not run
		// If useDelta then checks switchRenderbuffers to see if we have already generated the Voronoi diagram
		timeSiteDraw.start();
		if(!useDelta || !switchRenderbuffers(t)){
			bool isFinalSeek = (j + 1 > timeRadius);
			renderInputPoints(t,wp,vec,timeRadius,useDelta, finalSeekTo, firstSeek, isFinalSeek);
		} else { 
			cerr << "Saved renderbuffer used for time " << t << "\n";
		}
		glFinish(); 
		timeSiteDraw.stop();

		// Useful for debugging
		//outputColorBuffer(x_pixels, y_pixels, query_spacing);
		//outputDepthBuffer(x_pixels, y_pixels, query_spacing);

		// Perform queries
		performNNIQueries(wp,scaling,blockWidth,pixel_offset,spacingFactor,j);

	}

	// Complete the NNI
	completeBufferAnalysis(wp, out, i_map, j_map);


	delete[] i_map;
	delete[] j_map;

	cerr << "Plane min/max: " << plane_min_z << ", " << plane_max_z << "\n";
	printTimes();
}

