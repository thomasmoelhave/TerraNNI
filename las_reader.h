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

#ifndef NNI_IO_LAS_READER
#define NNI_IO_LAS_READER

#include <tpie_config.h>
#include "point_reader.h"

#include <fstream>
#include <liblas/reader.hpp>
#include <liblas/header.hpp>
#include <boost/cstdint.hpp>

namespace nni_io {
	template<typename P, typename Z>
	class las_reader : public point_reader<P, Z>{
		public:
			virtual ~las_reader(){this->close();}

			const liblas::Header& get_header(){
				assert(reader!=NULL);
				const liblas::Header &head = reader->GetHeader();
				return head;
			}

		protected:
			virtual void do_open(const std::string &filename);
			virtual void do_close();
			virtual bool get_next_point(P &point);
			virtual void do_reset();

		private:
			//ignore point (use classification to decide)
			bool ignore(const liblas::Point& p) const;
			std::ifstream *stream;
			liblas::Reader *reader;
	};
#include "las_reader.inl"
}
#endif
