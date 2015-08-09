//
// sounding.h
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

#ifndef DENISE_SRC_SOUNDING_H
#define DENISE_SRC_SOUNDING_H

#include <map>
#include <denise/thermo.h>
#include "package.h"

using namespace std;
using namespace denise;

namespace andrea
{

   class Sounding_Package : public Package
   {

      protected:

         map<string, Sounding>
         sounding_map;

         Tephigram
         tephigram;

         Sounding_Package (Andrea& andrea);

         void
         sounding_load (const string& identifier,
                        const string& file_path);

         void
         sounding_print (const string& identifier,
                         const Tokens& arguments) const;

         void
         sounding_parse (const Tokens& tokens);

      public:

         const map<string, Sounding>&
         get_sounding_map () const;

         const Sounding&
         get_sounding (const string& identifier) const;

         void
         surface_sounding (const Tokens& tokens) const;

         void
         surface_sounding_tephigram (const Tokens& tokens) const;

         void
         surface_sounding_chart (const Tokens& tokens) const;

         void
         surface_sounding_chart (const RefPtr<Context>& cr,
                                 const Transform_2D& transform,
                                 const bool is_p,
                                 const Mesh_2D& mesh_2d,
                                 const string& fmt_x,
                                 const string& fmt_y,
                                 const string& unit_x,
                                 const string& unit_y,
                                 const Sounding& sounding,
                                 const Real_Profile& real_profile,
                                 const Symbol& symbol,
                                 const Color& color) const;


   };

};

#endif /* DENISE_SRC_SOUNDING_H */
