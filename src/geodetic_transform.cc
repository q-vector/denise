//
// geodetic_transform.cc
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

#include "geodetic_transform.h"

using namespace std;
using namespace denise;
using namespace andrea;

Geodetic_Transform_Package::Geodetic_Transform_Package (Andrea& andrea)
   : Package (andrea)
{
}

void
Geodetic_Transform_Package::geodetic_transform_assign (const Dstring& identifier,
                                                       const Dstring& str)
{
   geodetic_transform_str_map[identifier] = str;
}

void
Geodetic_Transform_Package::geodetic_transform_print (const Dstring& identifier,
                                                      const Tokens& arguments) const
{
   auto iterator = geodetic_transform_str_map.find (identifier);
   const bool is_present = (iterator != geodetic_transform_str_map.end ());
   if (is_present)
   {
      cout << "geodetic_transform " << identifier << " is present" << endl;
   }
}

void
Geodetic_Transform_Package::geodetic_transform_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "assign")
   {
      const Dstring& identifier = tokens[1];
      const Dstring& str = tokens[2];
      geodetic_transform_assign (identifier, str);
   }
   else
   if (tokens[0] == "print")
   {
      const Dstring& identifier = tokens[1];
      geodetic_transform_print (identifier, tokens.subtokens (2));
   }

}

const map<Dstring, Dstring>&
Geodetic_Transform_Package::get_geodetic_transform_str_map () const
{
   return geodetic_transform_str_map;
}

const Geodetic_Transform*
Geodetic_Transform_Package::get_geodetic_transform_ptr (const Dstring& identifier,
                                                        const Point_2D& point) const
{
   Exception e ("geodetic_transform not found: " + identifier);
   try
   {
      typedef Geodetic_Transform Gt;
      const Dstring& str = geodetic_transform_str_map.at (identifier);
      return Gt::get_transform_ptr (str, point);
   }
   catch (const std::out_of_range& oor) { throw e; }
}

