//
// journey.cc
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

#include "andrea.h"
#include "journey.h"

using namespace std;
using namespace denise;
using namespace andrea;

Journey_Package::Journey_Package (Andrea& andrea)
   : Package (andrea)
{
}

void
Journey_Package::journey_assign (const Dstring& identifier,
                                 const Tokens& arguments)
{

   Journey journey;

   for (auto iterator = arguments.begin ();
        iterator != arguments.end (); iterator++)
   {
      const Lat_Long lat_long (*(iterator));
      journey.push_back (lat_long);
   }

   journey_map[identifier] = journey;

}

void
Journey_Package::journey_print (const Dstring& identifier) const
{
   auto iterator = journey_map.find (identifier);
   const bool is_present = (iterator != journey_map.end ());
   if (is_present)
   {
      wcout << L"journey " << identifier << L" is present" << endl;
   }
}

void
Journey_Package::journey_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == L"assign")
   {
      const Dstring& identifier = tokens[1];
      journey_assign (identifier, tokens.subtokens (2));
   }
   else
   if (tokens[0] == L"print")
   {
      const Dstring& identifier = tokens[1];
      journey_print (identifier);
   }

}

const map<Dstring, Journey>&
Journey_Package::get_journey_map () const
{
   return journey_map;
}

const Journey&
Journey_Package::get_journey (const Dstring& identifier) const
{
   Exception e (L"journey not found: " + identifier);
   try { return journey_map.at (identifier); }
   catch (const std::out_of_range& oor) { throw e; }
}

void
Journey_Package::surface_journey (const Dstring& surface_identifier,
                                  const Dstring& geodetic_transform_identifier,
                                  const Dstring& journey_identifier)
{

   const RefPtr<Surface>& surface = andrea.get_surface (surface_identifier);
   const RefPtr<Context> cr = andrea.get_cr (surface_identifier);
   const Size_2D& size_2d = andrea.get_size_2d (surface_identifier);
   const Point_2D centre (Real (size_2d.i) / 2, Real (size_2d.j) / 2);

   const Geodetic_Transform* geodetic_transform_ptr =
      andrea.get_geodetic_transform_ptr (geodetic_transform_identifier, centre);
   const Geodetic_Transform& geodetic_transform = *geodetic_transform_ptr;

   const Journey& journey = andrea.get_journey (journey_identifier);

   cr->save ();
   Color::black ().cairo (cr);
   journey.cairo (cr, geodetic_transform);
   cr->restore ();

   delete geodetic_transform_ptr;

}

