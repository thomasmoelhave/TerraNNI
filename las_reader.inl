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

template<typename P, typename Z>
void las_reader<P, Z>::do_open(const std::string &filename){
	this->stream = new std::ifstream(filename.c_str(), std::ios::binary | std::ios::in);
	if(!this->stream){
		throw std::runtime_error("Could not open file " + filename);
	}
	try{
		this->reader = new liblas::LASReader(*stream);
	} catch(const std::runtime_error &err){
		throw std::runtime_error("Could not open file " + filename + " an error happened in libLAS: " + err.what());
	}

	const liblas::LASHeader &head = reader->GetHeader();
}

template<typename P, typename Z>
void las_reader<P, Z>::do_close(){
	delete this->reader;
	delete this->stream;
}

template<typename P, typename Z>
bool las_reader<P, Z>::ignore(const liblas::LASPoint& p) const {
	const liblas::uint8_t c = p.GetClassification();
	//return (c != 2);
	//never ignore, regardless of classification
	return false;
}

template<typename P, typename Z>
bool las_reader<P, Z>::get_next_point(P &point){
	try {
		while(reader->ReadNextPoint()){
			const liblas::LASPoint &p = reader->GetPoint();
			if (!ignore(p)) {
				point = P(static_cast<double>(p.GetX()), static_cast<double>(p.GetY()), static_cast<double>(p.GetZ()), (int)p.GetIntensity());
				return true;
			}
		}
		return false;
	} catch (...) {
		//liblas threw an exception, skipping the rest of the file.
		std::cerr << "failure occured while reading points from LAS file.\n";
		return false;
	}
}

template<typename P, typename Z>
void las_reader<P, Z>::do_reset()
{
	this->stream->seekg (0, std::ios::beg);
	delete reader;
	this->reader = new liblas::LASReader(*stream);
	const liblas::LASHeader &head = reader->GetHeader();
}
