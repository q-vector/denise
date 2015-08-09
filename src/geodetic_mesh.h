//
// geodetic_mesh.h
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

#ifndef DENISE_SRC_GEODETIC_MESH_H
#define DENISE_SRC_GEODETIC_MESH_H

#include <map>
#include <denise/geodesy.h>
#include "package.h"

using namespace std;
using namespace denise;

namespace andrea
{

   class Geodetic_Mesh_Package : public Package
   {

      protected:

         map<string, Geodetic_Mesh>
         geodetic_mesh_map;

         Geodetic_Mesh_Package (Andrea& andrea);

         void
         geodetic_mesh_assign (const string& identifier,
                               const Size_2D& size_2d,
                               const Domain_2D& domain_2d);

         void
         geodetic_mesh_add (const string& identifier,
                            const Tokens& arguments);

         void
         geodetic_mesh_print (const string& identifier,
                              const Tokens& arguments) const;

         void
         geodetic_mesh_parse (const Tokens& tokens);

      public:

         const map<string, Geodetic_Mesh>&
         get_geodetic_mesh_map () const;

         const Geodetic_Mesh&
         get_geodetic_mesh (const string& identifier) const;

         void
         surface_geodetic_mesh (const string& surface_identifier,
                                const string& geodetic_transform_identifier,
                                const string& geodetic_mesh_identifier);

   };

};

#endif /* DENISE_SRC_GEODETIC_MESH_H */
