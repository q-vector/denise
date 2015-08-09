//
// andrea.h
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

#ifndef DENISE_SRC_ANDREA_H
#define DENISE_SRC_ANDREA_H

#include <map>
#include <denise/dstring.h>
#include "package.h"
#include "geodesy.h"
#include "geodetic_mesh.h"
#include "geodetic_transform.h"
#include "gshhs.h"
#include "journey.h"
#include "sounding.h"
#include "surface.h"

using namespace std;
using namespace denise;

namespace andrea
{

   class Entity : public string
   {

      public:

         Entity (const string& str);

         Real
         value () const;

   };

   class Andrea : public Geodesy_Package,
                  public Geodetic_Mesh_Package,
                  public Geodetic_Transform_Package,
                  public Gshhs_Package,
                  public Journey_Package,
                  public Sounding_Package,
                  public Surface_Package
   {

      private:

         string
         prompt;

         void
         wind_shear (const Tokens& arguments) const;

         void
         print (const Entity& entity) const;

      public:

         Andrea (const string& prompt);

         virtual void
         parse (const Tokens& tokens);

         void
         loop ();

   };

};

#endif /* DENISE_SRC_ANDREA_H */

