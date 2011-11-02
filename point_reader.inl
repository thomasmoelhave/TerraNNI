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
void point_reader<P, Z>::open(const std::string &filename){
	this->close();

	this->do_open(filename);
	is_open=true;
}

template<typename P, typename Z>
void point_reader<P,Z>::close(){
	if(is_open)
		this->do_close();
	is_open=false;
}

template<typename P, typename Z>
bool point_reader<P, Z>::next_point(P &point){
	if(!is_open)
		return false;

	if(this->get_next_point(point)){
		return true;
	}
	return false;
}
