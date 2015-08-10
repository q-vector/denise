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

         map<Dstring, Gshhs*>
         gshhs_ptr_map;

         Gshhs_Package (Andrea& andrea);

         ~Gshhs_Package ();

         void
         gshhs_load (const Dstring& identifier,
                     const Dstring& file_path);

         void
         gshhs_print (const Dstring& identifier,
                      const Tokens& arguments) const;

         void
         gshhs_parse (const Tokens& tokens);

      public:

         const map<Dstring, Gshhs*>&
         get_gshhs_ptr_map () const;

         const Gshhs*
         get_gshhs_ptr (const Dstring& identifier) const;

         void
         surface_gshhs (const Dstring& surface_identifier,
                        const Dstring& geodetic_transform_identifier,
                        const Dstring& gshhs_identifier,
                        const Tokens& arguments);

   };

};

#endif /* DENISE_SRC_GSHHS_H */
