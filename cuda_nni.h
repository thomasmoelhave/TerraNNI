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

#ifndef LIBCUDA_NNI_H
#define LIBCUDA_NNI_H

#include "nni.h"

namespace cuda_nni {

	void init(const nni::Settings& s);
	void setup_run(int nrows, int ncols);
	void cleanupCUDA();

	void cudaSwitchVoronoiRB(int rbIndex);
	void cudaSetup(int nrows, int ncols);
	void bufferAnalysis(const nni::WorkPackage& wp, int scaling, int
			blockWidth, int pixel_offset, int curRound, int rounds);
	void completeAnalysis(int nrows, int ncols, int ** res, int zMult);
	void cudaRegisterAllBuffers(int tC, int * rb_voronois, int rb_query);

}

#endif
