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

#ifndef LIBNNI_H
#define LIBNNI_H

#include "tpie_config.h"
#include <tpie/stream.h>
#include <vector>
#include <tpie/cpu_timer.h>

namespace nni {

	struct Point3D {
		double p[3]; //x,y,z
		Point3D(double x,double y,double z, double t = 0) {
			p[0]=x;
			p[1]=y;
			p[2]=z;
		}
		Point3D() {}
	};

	struct Point4D : public Point3D {
		int t;
		Point4D(double x, double y, double z, int time =0) {
			p[0]=x;
			p[1]=y;
			p[2]=z;
			t = time;
		}
		Point4D() {}
	};

	struct OutPoint {
		int i,j;
		float h;
		OutPoint(int x, int y, float z) {
			i=x;
			j=y;
			h=z;
		}
		OutPoint() {}
	};

	struct WorkPackage {
		float origin_x,origin_y; //grid offset
		float grid_origin_x, grid_origin_y;
		int nrows,ncols;
	};

	struct Settings {
		int scale;
		float cell_size; //cell dimension
		int nrows,ncols;
		float site_radius;
		double query_radius;
		int algorithm;
		int pointMin;
		int multiplier;
	};

	void cleanup();

	int init(const Settings& s);

	void setRenderbuffers(int w_width, int w_height, int tC = 1);

	void printTimes();

	void setupRun(const WorkPackage& wp, int pixel_offset = -1, int timeCount = 1);

	void performNNI(const WorkPackage& wp, tpie::ami::stream<OutPoint>& out,
			tpie::ami::stream<nni::Point4D>& vec, float timeRadius, float curTime,
			bool useDelta, bool debug = false);

	void clearRenderbuffers();

}

#endif
