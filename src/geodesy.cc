//
// geodesy.h
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

#include "geodesy.h"

using namespace std;
using namespace denise;
using namespace andrea;

void
Geodesy_Package::geodesy_assign (const Dstring& identifier,
                                 const Dstring& str)
{
   geodesy_map[identifier] = Geodesy (str);
}

void
Geodesy_Package::geodesy_print (const Dstring& identifier) const
{

   auto iterator = geodesy_map.find (identifier);
   const bool is_present = (iterator != geodesy_map.end ());
   if (is_present)
   {
      wcout << L"geodesy " << identifier << L" is present" << endl;
   }

}

void
Geodesy_Package::geodesy_distance (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 2) { throw Exception (L"geodesy_distance"); }

   const Lat_Long origin (tokens[0]);
   const Lat_Long destination (tokens[1]);
   const Real distance = Geodesy::get_distance (origin, destination);

   wcout << distance << endl;

}

void
Geodesy_Package::geodesy_azimuth (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 2) { throw Exception (L"geodesy_azimuth"); }

   const Lat_Long origin (tokens[0]);
   const Lat_Long destination (tokens[1]);
   const Real azimuth_f = Geodesy::get_azimuth_forward (origin, destination);
   const Real azimuth_b = Geodesy::get_azimuth_backward (origin, destination);

   wcout << azimuth_f << L" " << azimuth_b << endl;

}

void
Geodesy_Package::geodesy_destination (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 3) { throw Exception (L"geodesy_destination"); }

   const Lat_Long origin (tokens[0]);
   const Real distance = stof (tokens[1]);
   const Real azimuth_forward = stof (tokens[2]);
   const Lat_Long& destination =
      Geodesy::get_destination (origin, distance, azimuth_forward);

   wcout << destination.get_string (lat_long_dp) << endl;

}

Geodesy_Package::Geodesy_Package (Andrea& andrea)
   : Package (andrea),
     lat_long_dp (4)
{
}

void
Geodesy_Package::geodesy_parse (const Tokens& tokens)
{


   const Integer n = tokens.size ();

   if (tokens[0] == L"assign")
   {
      const Dstring& identifier = tokens[1];
      geodesy_assign (identifier, tokens[2]);
   }
   else
   if (tokens[0] == L"print")
   {
      const Dstring& identifier = tokens[1];
      geodesy_print (identifier);
   }
   else
   if (tokens[0] == L"distance")
   {
      geodesy_distance (tokens.subtokens (1));
   }
   else
   if (tokens[0] == L"azimuth")
   {
      geodesy_azimuth (tokens.subtokens (1));
   }
   else
   if (tokens[0] == L"destination")
   {
      geodesy_destination (tokens.subtokens (1));
   }
   else
   {
      throw Exception (L"geodesy_parse");
   }

}

const map<Dstring, Geodesy>&
Geodesy_Package::get_geodesy_map () const
{
   return geodesy_map;
}

const Geodesy&
Geodesy_Package::get_geodesy (const Dstring& identifier) const
{
   Exception e (L"geodetsy not found: " + identifier);
   try { return geodesy_map.at (identifier); }
   catch (const std::out_of_range& oor) { throw e; }
}

