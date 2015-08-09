//
// gshhs.h
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

#ifndef DENISE_SRC_GSHHS_H
#define DENISE_SRC_GSHHS_H

#include <map>
#include <denise/gis.h>
#include "package.h"

using namespace std;
using namespace denise;

namespace andrea
{

   class Gshhs_Package : public Package
   {

      protected:

         map<string, Gshhs*>
         gshhs_ptr_map;

         Gshhs_Package (Andrea& andrea);

         ~Gshhs_Package ();

         void
         gshhs_load (const string& identifier,
                     const string& file_path);

         void
         gshhs_print (const string& identifier,
                      const Tokens& arguments) const;

         void
         gshhs_parse (const Tokens& tokens);

      public:

         const map<string, Gshhs*>&
         get_gshhs_ptr_map () const;

         const Gshhs*
         get_gshhs_ptr (const string& identifier) const;

         void
         surface_gshhs (const string& surface_identifier,
                        const string& geodetic_transform_identifier,
                        const string& gshhs_identifier,
                        const Tokens& arguments);

   };

};

#endif /* DENISE_SRC_GSHHS_H */
