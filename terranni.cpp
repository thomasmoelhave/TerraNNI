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


#include "tpie_config.h"
#include <iostream>
#include <tpie/tpie.h>
#include <tpie/sort.h>
#include <tpie/stream.h>
#include <tpie/cpu_timer.h>
#include <tpie/sort.h>
#include <cmath>
#include <map>
#include <algorithm>
#include <algorithm>
#include <string.h>
#include <boost/filesystem.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/environment_iterator.hpp>
#include <boost/program_options/eof_iterator.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/program_options/option.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/positional_options.hpp>
#include <boost/program_options/value_semantic.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/version.hpp>

#include "nni.h"

#include "las_reader.h"

using namespace boost::filesystem;
using namespace std;

struct tile {
	path lasfile;
};

nni::Settings s;


/**************************
 * Get File Data
 *************************/

// Store information about all tiles
struct gridheader {
	float minx,miny,maxx,maxy,minz,maxz;
};

void update_super(const liblas::LASHeader &newtile, gridheader& super) {
	super.minx = std::min(super.minx,(float)newtile.GetMinX());
	super.miny = std::min(super.miny,(float)newtile.GetMinY());
	super.minz = std::min(super.minz,(float)newtile.GetMinZ());
	super.maxx = std::max(super.maxx,(float)newtile.GetMaxX());
	super.maxy = std::max(super.maxy,(float)newtile.GetMaxY());
	super.maxz = std::max(super.maxz,(float)newtile.GetMaxZ());
}


// Recursively find tiles in directory
void find_tiles(
		const path & dir_path,         // in this directory,
		const std::string & ext, // search for this extension,
		vector<tile>& tilevec,
		gridheader& super_grid
		)   // add results here
{
	if ( !exists( dir_path ) ) return ;

	directory_iterator end_itr; // default construction yields past-the-end
	for ( directory_iterator itr( dir_path );
			itr != end_itr;
			++itr )
	{
		if ( is_directory(itr->status()) )
		{
			find_tiles( itr->path(), ext, tilevec,super_grid);
		} else if ( extension(itr->leaf()) == ext ) {
			tile t;
			t.lasfile=itr->path();
			tilevec.push_back(t);

			nni_io::las_reader<nni::Point4D,double> pr;
			pr.open(t.lasfile.string());
			update_super(pr.get_header(), super_grid);
		}
	}
}



// Simple struct to sort points in TPIE stream by time (ascending)
struct Timecmp{
	Timecmp()
	{}
	int compare(const nni::Point4D& p1, const nni::Point4D& p2){
		return p1.t - p2.t;
	}
};

// Simple struct to sort points by rows and columns,
// used for outputting ASCII files
struct XYcmp{
	XYcmp()
	{}
	int compare(const nni::OutPoint& p1, const nni::OutPoint& p2){
		if(p1.j != p2.j)
			return p2.j - p1.j;
		return p1.i - p2.i;
	}
};

// Calculate the number of streams that can be handled at once given available memory
// Used to determine if the binning procedure must be recursive
int streamsPerPass() {
	tpie::ami::stream<nni::Point3D> temp;
	tpie::ami::err ae;	
	TPIE_OS_SIZE_T mmBytesPerStream;	
	TPIE_OS_SIZE_T mmBytesAvail = tpie::MM_manager.consecutive_memory_available();
	if ((ae = temp.main_memory_usage(&mmBytesPerStream,tpie::mem::STREAM_USAGE_MAXIMUM)) != tpie::ami::NO_ERROR) {
		fprintf(stderr, "TPIE MM_manager error\n");
		return -1;
	}
	return (int)min(floor(mmBytesAvail/mmBytesPerStream), 1.0*temp.available_streams())-1;
}

// Generic class to read points from tiles
class tile_reader {
	public:
		tile_reader(vector<tile>& tv): tilevec(tv), cur_tile(-1) {}
		void setTiles(vector<tile>& tv) {tilevec = tv;}
		bool next_point(nni::Point4D& p) {
			if(cur_tile == -1){
				cur_tile++;
				pr.open(tilevec[cur_tile].lasfile.string());
			}
			//Fix this up!
			if(pr.next_point(p)) {
				return true;
			} else if ( cur_tile + 1 < tilevec.size()) {
				cur_tile++;
				pr.open(tilevec[cur_tile].lasfile.string());
				pr.next_point(p);
				return true;
			} else {
				return false;
			}
		}
	private:
		vector<tile>& tilevec;
		int cur_tile;
		nni_io::las_reader<nni::Point4D,double> pr;
};

// Generic class to read points from TPIE stream
class tpie_stream_reader {
	public:
		tpie_stream_reader(tpie::ami::stream<nni::Point4D>& nwp) : wp(nwp) {}
		bool next_point(nni::Point4D& p){
			nni::Point4D * p_ptr;
			if (wp.read_item(&p_ptr) == tpie::ami::NO_ERROR) {
				p = *p_ptr;
				return true;
			} else {
				return false;
			}
		}

	private:
		tpie::ami::stream<nni::Point4D>& wp;
};


// Bin data in x,y space such that each work package tile has its own stream
// containing all necessary points 
template <typename Reader>
void bin(const int width, const int height, const float border, const
		gridheader& super, const std::string wppath,const int maxStreams,
		Reader& r, int level, int offset_x, int offset_y, bool twoD){

	int cols = ceil(1.0*(super.maxx-super.minx)/width);
	int rows = ceil(1.0*(super.maxy-super.miny)/height);

	int maxWidth = (int)sqrt(maxStreams);

	int p_width = (int)(width * pow(maxWidth,level));
	int p_height = (int)(height * pow(maxWidth,level));
	int c_width = (int)(width * pow(maxWidth,level-1));
	int c_height = (int)(height * pow(maxWidth,level-1));

	int num_cols = maxWidth;
	int num_rows = maxWidth;

	//NEED TO IMPLEMENT CUSTOM WIDTH AND HEIGHTS BASED ON OVERALL BOUNDS

	if(level == 1) {
		fprintf(stderr, "cur width and height: %d, %d\n", c_width, c_height);
		fprintf(stderr, "offset: %d, %d\n", offset_x, offset_y );
	}	


	// Create all TPIE streams for this level of binning
	std::stringstream ss;
	tpie::ami::stream<nni::Point4D>** wps = new tpie::ami::stream<nni::Point4D>*[maxStreams];
	for(int i = 0; i < maxWidth; i++) {
		for(int j = 0; j < maxWidth; j++) {
			ss.str("");
			if(level == 1)
				ss << wppath << "/" << "r" << i+offset_y << "c" << j+offset_x << ".tpie";
			else
				ss << wppath << "/l" << level << "r" << i+offset_y << "c" << j+offset_x << ".tpie";
			wps[i*maxWidth + j] = new tpie::ami::stream<nni::Point4D>(ss.str());
			wps[i*maxWidth + j]->truncate(0);
		}
	}

	int cnt =0;
	nni::Point4D p;
	while(r.next_point(p)) {
		cnt++;

		if(twoD) {
			p.t = 0;
		}

		float tx = (p.p[0]-super.minx - offset_x * width);
		float ty = (p.p[1]-super.miny - offset_y * height);
		int col = (int)(tx/c_width);
		int row = (int)(ty/c_height);

		//Add point p to the bin it falls directly within
		int index = row*maxWidth + col;
		if(row >= 0 && col >= 0 && row < maxWidth && col < maxWidth){
			if(!(level==1 && (row + offset_y >= rows || col + offset_x >= cols)) )
				wps[index]->write_item(p);
		}


		// Check if point p also falls within neighboring work packages border region
		// Border region makes sure that points are included if their region of
		// influence can effect the interpolation of a work package, even if
		// the point itself is outside the work package

		int row_offset = 0;
		int col_offset = 0;
		// Check if the point falls within adjacent work packages on the left or right
		if(tx - (floor(tx/c_width)*c_width) <= border) { //left
			col_offset = -1;
		} else if (tx - (floor(tx/c_width)*c_width) >= (c_width-border)) {//right
			col_offset = 1;
		}
		// If so and it is a valid work package to add points to, add the point
		if(col_offset != 0 && row >= 0 && col+col_offset >= 0 && row < maxWidth && col+col_offset < maxWidth)
			if(!(level==1 && (row + offset_y >= rows || col + col_offset + offset_x >= cols)) )
				wps[row*maxWidth + (col+col_offset)]->write_item(p);

		// Check if the point falls within adjacent work packages on the top or bottom
		if(ty - (floor(ty/c_height)*c_height) <= border) { //bottom
			row_offset = 1;
		} else if (ty - (floor(ty/c_height)*c_height) >= (c_height-border)) {//top
			row_offset = -1;
		}

		// Again, if so and it is a valid work package, add it
		if(row_offset != 0 && row+row_offset >= 0 && col >= 0 && row+row_offset < maxWidth && col < maxWidth)
			if(!(level==1 && (row + row_offset + offset_y >= rows || col + offset_x >= cols)) )
				wps[(row+row_offset)*maxWidth + col]->write_item(p);

		// Check if the point is in the corner region such that it should be
		// added to a diagonally adjacent work package, and if so add it
		if(row_offset != 0 && col_offset != 0) {
			if((row+row_offset) >= 0 && (col+col_offset) >= 0 && (row+row_offset) < maxWidth && (col+col_offset) < maxWidth)
				if(!(level==1 && (row + row_offset + offset_y >= rows || col + col_offset + offset_x >= cols)) )
					wps[(row+row_offset)*maxWidth + (col+col_offset)]->write_item(p);
		}
	}

	tpie::cpu_timer timeDataSort;

	// Sort points in time (necessary for interpolation and seeking order in NNI)
	// After sorted delete pointer to the stream
	timeDataSort.start();
	Timecmp cmp;
	for(int i = 0; i < maxStreams; i++) {
		if(level == 1) {
			tpie::ami::sort(wps[i], &cmp);
			wps[i]->seek(0);
		}
		delete wps[i];
	}
	delete[] wps;
	timeDataSort.stop();
	fprintf(stdout, "\tTime to sort data (by time): \t\t%f\n", timeDataSort.wall_time());


	// If necessary recurse in each area, one at a time
	if(level != 1){ 
		for(int i = 0; i < maxWidth; i++) { //row
			for(int j = 0; j < maxWidth; j++) { //col
				if(offset_x + (j * pow(maxWidth,level-1)) > cols || offset_y + (i * pow(maxWidth, level-1)) > rows)
					continue;
				ss.str("");
				ss << wppath << "/l" << level << "r" << i+offset_y << "c" << j+offset_x << ".tpie";

				tpie::ami::stream<nni::Point4D> wp(ss.str());
				wp.seek(0);
				tpie_stream_reader temp_reader(wp);
				bin<tpie_stream_reader>(width,height,border,super, wppath, maxStreams, 
						temp_reader, level-1, offset_x + (j * pow(maxWidth,level-1)), offset_y + (i*pow(maxWidth,level-1)),twoD);

				wp.truncate(0);
			}
		}
	}
}


// Set up binning procedure and then bin input points
void binSetup(const gridheader& super, vector<tile>& tilevec,const int width, const int height, const float border, const std::string wppath, bool twoD){
	fprintf(stderr, "Set up bins...\n");

	int maxStreams = (int)pow(floor(sqrt(streamsPerPass())),2);

	fprintf(stderr, "Maximum of %d streams.\n", maxStreams);	
	int cols = ceil(1.0*(super.maxx-super.minx)/width);
	int rows = ceil(1.0*(super.maxy-super.miny)/height);
	fprintf(stderr, "Binning rows,cols: %d,%d\n", rows, cols);

	int levels = 1;
	while((rows*cols) >= pow(maxStreams,levels)) levels++;
	fprintf(stderr, "Optimal: %d level(s) required.\n", (levels));	

	levels = 1;
	while((rows) >= pow(sqrt(maxStreams),levels)) levels++;
	while((cols) >= pow(sqrt(maxStreams),levels)) levels++;
	fprintf(stderr, "%d level(s) required.\n", (levels));	

	tile_reader tr(tilevec);
	bin<tile_reader>(width,height,border,super, wppath, maxStreams, tr, levels,0,0, twoD);

}



nni::WorkPackage get_workpackage(const int width, const int height, const float origin_x, const float origin_y, const gridheader& super) {

	nni::WorkPackage wp;
	wp.origin_x = origin_x;
	wp.origin_y = origin_y;

	wp.ncols= (int)(width/s.cell_size);
	wp.nrows= (int)(height/s.cell_size);

	wp.grid_origin_x=(wp.origin_x-super.minx)/s.cell_size;
	wp.grid_origin_y=(wp.origin_y-super.miny)/s.cell_size;

	cerr << "New work package: \n"
		<< "\tOrigin: " << origin_x << ", " << origin_y << "\n"
		<< "\tCols,Rows: " << wp.ncols << ", " << wp.nrows << "\n"
		<< "\tGrid origin: " << wp.grid_origin_x << ", " << wp.grid_origin_y << "\n"
		;

	return wp;
}

// Run required interpolation on each work package
void binRun(gridheader& super, int width, int height, std::string wppath, std::string outfile, float timeStart, float timeLength, float timeStep, int timeRadius, bool useDelta) {
	int cols = ceil(1.0*(super.maxx-super.minx)/width);
	int rows = ceil(1.0*(super.maxy-super.miny)/height);

	std::stringstream ss;
	bool firstBox = true;
	for(int row = 0; row < rows; row++) {
		for(int col = 0; col < cols; col++) {
			ss.str("");
			ss << wppath << "/" << "r" << row << "c" << col << ".tpie";
			tpie::ami::stream<nni::Point4D> wps(ss.str());

			cerr << "NNI on Tile: " << col << ", " << row << "\n";
			nni::WorkPackage wp = get_workpackage(width, height, col*width+super.minx, row*height+super.miny,super);

			for(float t = timeStart; t <= timeStart + timeLength; t += timeStep) {
				ss.str("");
				ss << outfile << "/" << setw(5) << left  << showpoint << t << ".tpie";

				tpie::ami::stream<nni::OutPoint> out_stream(ss.str(),tpie::ami::APPEND_STREAM);
				out_stream.persist();
				if(firstBox){
					out_stream.truncate(0);
				}

				nni::performNNI(wp,out_stream,wps,timeRadius,t,useDelta);
			}
			firstBox = false;
			nni::clearRenderbuffers();
		}
	}

}


// Output all interpolated data as ASCII grids (one per time slice interpolated)
void outputASC(float timeStart, float timeLength, float timeStep, int
		timeRadius, int nrrows, int nrcols, float xlcorner, float ylcorner,
		float cellsize, std::string outfile) { 

	cout << "Output ASCII\n";
	cout << "Cols, Rows: " << nrcols << ", " << nrrows << "\n";
	tpie::cpu_timer timeDataOut;
	std::stringstream ss;
	XYcmp cmp;

	timeDataOut.start();
	for(float t = timeStart; t <= timeStart + timeLength; t += timeStep) {
		cout << "Output time " << t << "\n";
		ss.str("");
		ss << outfile << "/" << setw(5) << left  << showpoint << t << ".tpie";

		tpie::ami::stream<nni::OutPoint> out_stream(ss.str(),tpie::ami::APPEND_STREAM);

		out_stream.seek(0);
		cout << "Stream length: " << out_stream.stream_len() << "\n";
		tpie::ami::sort(&out_stream, &cmp);
		out_stream.seek(0);

		ss << ".asc";

		ofstream outfile(ss.str().c_str());

		cout << "Output to " << ss.str() << "\n";

		outfile 
			<< "NROWS " << nrrows << "\n"
			<< "NCOLS " << nrcols << "\n"
			<< "XLLCORNER " << xlcorner << "\n"
			<< "YLLCORNER " << ylcorner << "\n"
			<< "CELLSIZE " << cellsize << "\n"
			<< "NODATA_VALUE -9999\n";

		nni::OutPoint * p_ptr;
		for(int j = nrrows-1; j >= 0; j--) {
			for(int i = 0; i < nrcols; i++) {
				if( out_stream.read_item(&p_ptr) == tpie::ami::NO_ERROR){
					const nni::OutPoint& p = *p_ptr;
					if(p.i != i || p.j != j) 
						fprintf(stderr, "Error in output order: %d,%d -- %d,%d\n",i,j,p.i,p.j);
					outfile << p.h << " "; 
					//outfile << p.i << " " << p.j << " " << p.h << "\n"; 
				} else {
					cerr << "Stream appears to be incorrect\n";
					break;
				}
			}
			outfile << "\n";
		}

	}
	timeDataOut.stop();

	fprintf(stderr, "Total Print Output time: %f\n", timeDataOut.wall_time());
}


int main (int argc, char** argv) {
	
	std::cout
		<< "TerraNNI Copyright (C) 2010, 2011 Pankaj K. Agarwal, Alex Beutel & Thomas Mølhave\n"
		<< "This program comes with ABSOLUTELY NO WARRANTY.\n"
		<< "This is free software, and you are welcome to redistribute it\n"
		<< "under certain conditions.\n\n";

	tpie::tpie_init();

	tpie::MM_manager.set_memory_limit(2000*1024*1024);
	tpie::MM_manager.warn_memory_limit();

	namespace po = boost::program_options;


	float timeStep = 1;
	float timeLength = 1;
	float timeStart = 0;
	int timeRadius = 0;

	int algo = 0;

	bool twoD = false;

	std::string tmppath=".";
	std::string outfile = "./";
	std::string wppath=".";
	std::string tilepath;

	s.site_radius = 5.0;

	int wp_width;
	int wp_height;

	po::options_description desc("Allowed options");
	desc.add_options()
		("tilepath", po::value<std::string>(&tilepath), "set path to tiles (Required)")
		("output", po::value<std::string>(&outfile)->default_value("./"), "required: set output folder")
		("scaling", po::value<int>(&(s.scale))->default_value(5), "set scaling factor")
		("cell-size", po::value<float>(&(s.cell_size))->default_value(2.0), "set cell size")
		("site-radius", po::value<float>(&(s.site_radius))->default_value(5.0), "set cone radius for input points")
		("query-radius", po::value<float>(), "set cone radius for queries, default is largest possible that optimizes speed")
		("wp-width", po::value<int>(&wp_width), "set tile width")
		("wp-height", po::value<int>(&wp_height), "set tile height")
		("tmp-path", po::value<std::string>(&tmppath)->default_value("."), "set path for temporary files")
		("wp-path", po::value<std::string>(&wppath)->default_value("."), "set path for work package files (TPIE streams)")
		("wp-computed", "use precomputed TPIE streams for work package construction")
		("origin-x", po::value<float>(), "set origin for grid on x-axis")
		("origin-y", po::value<float>(), "set origin for grid on y-axis")
		("grid-cols", po::value<int>(), "number of grid columns")
		("grid-rows", po::value<int>(), "number of grid rows")
#ifndef TERRANNI_2DONLY
		("time-radius", po::value<int>(&timeRadius)->default_value(0), "time radius")
		("time-start", po::value<float>(&timeStart)->default_value(0), "time to begin interpolation")
		("time-length", po::value<float>(&timeLength)->default_value(1), "length of time to interpolate over")
		("time-step", po::value<float>(&timeStep)->default_value(1), "time increments to interpolate over")
		("2d", "Ignore the time for all points and render as 2D for one time slice.  time-radius, time-start, time-length, and time-step will be ignored.")
		("delta-squared", "use delta squared algorithm rather than delta")
#endif
		("draw-cones", "Draw Cones rather than Planes (using the lifting transform.")
		("help", "produce help message")
		;

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);    

	tpie::cpu_timer timeFindTile;
	tpie::cpu_timer timeToTile;
	tpie::cpu_timer timeToNNI;


	timeLength -= 1; // doing only 1 year means you step 0 times

	bool valid = true;
	if(!(vm.count("tilepath") || vm.count("help"))) {
		cout << "No tile path specified.\n";
		valid = false;
	}

	if (vm.count("help") || !valid) {
		cout << desc << "\n";
		return EXIT_SUCCESS;
	}

#ifdef TERRANNI_2DONLU
	twoD = true;
#else
	twoD = vm.count("2d");
#endif 

	if(twoD){
		cerr << "Running 2D version.\n";
		timeRadius = 0;
		timeStart = 0;
		timeLength = 0;
		timeStep = 1;
	}


	std::vector<tile> tilevec;

	tpie::tempname::set_default_path(tmppath);
	cout << "Storing temporary files in: " << tmppath << "\n";

	timeFindTile.start();

	gridheader super_grid;
	super_grid.minx = super_grid.miny = super_grid.minz = std::numeric_limits<int>::max();
	super_grid.maxx = super_grid.maxy = super_grid.maxz = std::numeric_limits<int>::min();

	find_tiles(tilepath,".las",tilevec,super_grid);

	if (tilevec.empty()) {
		cout << tilepath << " contains no .las files, aborting.";
		return EXIT_SUCCESS;
	}

	if(vm.count("origin-x")) super_grid.minx = vm["origin-x"].as<float>();
	if(vm.count("origin-y")) super_grid.miny = vm["origin-y"].as<float>();

	int nrcols = (super_grid.maxx-super_grid.minx)/s.cell_size;
	int nrrows = (super_grid.maxy-super_grid.miny)/s.cell_size;

	if(vm.count("grid-cols")){
		super_grid.maxx = vm["grid-cols"].as<int>() * s.cell_size + super_grid.minx;
		nrcols = vm["grid-cols"].as<int>();
	}
	if(vm.count("grid-rows")) {
		super_grid.maxy = vm["grid-rows"].as<int>() * s.cell_size + super_grid.miny;
		nrrows = vm["grid-rows"].as<int>();
	}

	cerr << "Z [min,max] = " << super_grid.minz << ", " << super_grid.maxz << "\n";
	s.pointMin = (int)floor(super_grid.minz);

	s.ncols = nrcols;
	s.nrows = nrrows;

	cout << "BBOX: x:["<<super_grid.minx <<", " << super_grid.maxx << "]. y:[" << super_grid.miny << ", " << super_grid.maxy << "]\n";

	timeFindTile.stop();


	// TODO: SHOULD VERIFY THAT THIS IS A FOLDER
	ofstream infofile((outfile+"/tpie.info").c_str());
	if(!infofile) {
		cerr << "Output path is not a folder.\n";
		exit(1);
	}

	infofile 
		<< "origin-x=" << super_grid.minx << "\n"
		<< "origin-y=" << super_grid.miny+s.cell_size*nrrows << "\n"
		<< "nrrows=" << s.nrows << "\n"
		<< "nrcols=" << s.ncols << "\n"
		<< "grid-length=" << s.cell_size << "\n";


	int blockWidth = 5;
	float border = blockWidth * s.cell_size + s.site_radius;

	//Set query radius to be half the width of one query block
	//This is as large as possible without performing multiple rounds of querying 
	double qr = blockWidth * s.cell_size / 2.0 - 0.01; // amb79: not - 0.01?
	if(vm.count("query-radius")){
		qr = vm["query-radius"].as<float>();
	}
	s.query_radius = qr;

	s.algorithm = 0;
	if(vm.count("draw-cones")) s.algorithm = 1;

	bool delta = !vm.count("delta-squared");


	int pw = (int)ceil(2.0*s.site_radius * s.scale / s.cell_size);  // pixel width of region of influence
	double max_sum = pw*pw * (super_grid.maxz - s.pointMin);

	s.multiplier = 1;
	if(max_sum < std::numeric_limits<int>::max()) {
		s.multiplier = (int)floor(1.0 * std::numeric_limits<int>::max() / max_sum);
		s.multiplier = (int)pow(10,floor(log10(s.multiplier)));  //To match with idea of decimal precision
	} else { 
		cerr << "Sizes so large that there may be integer overflow in CUDA resulting in faulty results.\n\n";
	}

	cerr 
		<< "Point minimum: " << s.pointMin << "\n"
		<< "Multiplier: " << s.multiplier << "\n"
		;


	// FIX THIS (amb79) - should automatically find best size
	cout << "Initializing NNI Library...\n";
	int max_size = nni::init(s);
	max_size = 8000;
	max_size = 5000;

	float max_meters = floor((max_size - (s.scale/s.cell_size) * 2 * border - 100) / s.scale) * s.cell_size; 
	fprintf(stderr, "Max size: %d pixels --> %f meters\n", max_size, max_meters);

	if(!vm.count("wp-width") || wp_width > max_meters) 
		wp_width = max_meters;

	if(!vm.count("wp-height") || wp_height > max_meters) 
		wp_height = max_meters;

	fprintf(stderr, "Work Package width and height: %d,%d\n", wp_width, wp_height);

	if(!vm.count("wp-computed")){
		timeToTile.start();
		binSetup(super_grid, tilevec, wp_width, wp_height, border, wppath,twoD);
		timeToTile.stop();
	}

	cout << "Outputting to: " << outfile << "\n";

	timeToNNI.start();
	binRun(super_grid,wp_width, wp_height, wppath, outfile, timeStart, timeLength, timeStep, timeRadius, delta);
	timeToNNI.stop();

	nni::cleanup();

	outputASC(timeStart,timeLength,timeStep,timeRadius,nrrows,nrcols,super_grid.minx,super_grid.miny,s.cell_size,outfile);

	cerr << "Total Find Tile time: " << timeFindTile.wall_time() << "\n"
		<< "Total Tile time: " << timeToTile.wall_time() << "\n"
		<< "Total NNI time: " << timeToNNI.wall_time() << "\n"
		;

	return EXIT_SUCCESS;
}


