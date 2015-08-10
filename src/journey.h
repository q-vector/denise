//
// journey.h
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

#ifndef DENISE_SRC_JOURNEY_H
#define DENISE_SRC_JOURNEY_H

#include <map>
#include <denise/geodesy.h>
#include "package.h"

using namespace std;
using namespace denise;

namespace andrea
{

   class Journey_Package : protected Package
   {

      protected:

         map<Dstring, Journey>
         journey_map;

         Journey_Package (Andrea& andrea);

         void
         journey_assign (const Dstring& identifier,
                         const Tokens& arguments);

         void
         journey_print (const Dstring& identifier) const;

         void
         journey_parse (const Tokens& tokens);

      public: 

         const map<Dstring, Journey>&
         get_journey_map () const;

         const Journey&
         get_journey (const Dstring& identifier) const;

         void
         surface_journey (const Dstring& surface_identifier,
                          const Dstring& geodetic_transform_identifier,
                          const Dstring& journey_identifier);

   };

};

#endif /* DENISE_SRC_JOURNEY_H */
