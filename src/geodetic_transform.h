//
// geodetic_transform.h
// 
// Copyright (C) 2005-2015 Simon E. Ching
// 
// This file is part of the denise suite.
//
// libdenise is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// libdenise is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with libdenise.  If not, see <http://www.gnu.org/licenses/>.

#ifndef DENISE_SRC_GEODETIC_TRANSFORM_H
#define DENISE_SRC_GEODETIC_TRANSFORM_H

#include <map>
#include <denise/geodesy.h>
#include "package.h"

using namespace std;
using namespace denise;

namespace andrea
{

   class Geodetic_Transform_Package : public Package
   {

      protected:

         map<Dstring, Dstring>
         geodetic_transform_str_map;

         Geodetic_Transform_Package (Andrea& andrea);

         void
         geodetic_transform_assign (const Dstring& identifier,
                                    const Dstring& str);

         void
         geodetic_transform_print (const Dstring& identifier,
                                   const Tokens& arguments) const;

         void
         geodetic_transform_parse (const Tokens& tokens);

      public:

         const map<Dstring, Dstring>&
         get_geodetic_transform_str_map () const;

         const Geodetic_Transform*
         get_geodetic_transform_ptr (const Dstring& identifier,
                                     const Point_2D& point) const;

   };

};

#endif /* DENISE_SRC_GEODETIC_TRANSFORM_H */
