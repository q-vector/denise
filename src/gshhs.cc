//
// gshhs.cc
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
#include "gshhs.h"

using namespace std;
using namespace denise;
using namespace andrea;

Gshhs_Package::Gshhs_Package (Andrea& andrea)
   : Package (andrea)
{
}

Gshhs_Package::~Gshhs_Package ()
{
   for (auto iterator = gshhs_ptr_map.begin ();
        iterator != gshhs_ptr_map.end (); iterator++)
   {
      Gshhs* gshhs_ptr = iterator->second;
      delete gshhs_ptr;
   }
}

void
Gshhs_Package::gshhs_load (const Dstring& identifier,
                           const Dstring& file_path)
{

   auto iterator = gshhs_ptr_map.find (identifier);
   const bool is_present = (iterator != gshhs_ptr_map.end ());
   if (is_present) { delete iterator->second; }

   Gshhs* gshhs_ptr = new Gshhs (file_path);
   gshhs_ptr_map[identifier] = gshhs_ptr;

}

void
Gshhs_Package::gshhs_print (const Dstring& identifier,
                            const Tokens& arguments) const
{
   auto iterator = gshhs_ptr_map.find (identifier);
   const bool is_present = (iterator != gshhs_ptr_map.end ());
   if (is_present)
   {
      cout << "gshhs " << identifier << " is present" << endl;
   }
}

void
Gshhs_Package::gshhs_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "load")
   {
      const Dstring& identifier = tokens[1];
      const Dstring& file_path = tokens[2];
      gshhs_load (identifier, file_path);
   }
   else
   if (tokens[0] == "print")
   {
      const Dstring& identifier = tokens[1];
      gshhs_print (identifier, tokens.subtokens (2));
   }

}

const map<Dstring, Gshhs*>&
Gshhs_Package::get_gshhs_ptr_map () const
{
   return gshhs_ptr_map;
}

const Gshhs*
Gshhs_Package::get_gshhs_ptr (const Dstring& identifier) const
{
   Exception e ("gshhs not found: " + identifier);
   try { return gshhs_ptr_map.at (identifier); }
   catch (const std::out_of_range& oor) { throw e; }
}

void
Gshhs_Package::surface_gshhs (const Dstring& surface_identifier,
                              const Dstring& geodetic_transform_identifier,
                              const Dstring& gshhs_identifier,
                              const Tokens& arguments)
{

   const RefPtr<Surface>& surface = andrea.get_surface (surface_identifier);
   const RefPtr<Context> cr = andrea.get_cr (surface_identifier);
   const Size_2D& size_2d = andrea.get_size_2d (surface_identifier);
   const Point_2D centre (Real (size_2d.i) / 2, Real (size_2d.j) / 2);

   const Geodetic_Transform* geodetic_transform_ptr =
      andrea.get_geodetic_transform_ptr (geodetic_transform_identifier, centre);
   const Geodetic_Transform& geodetic_transform = *geodetic_transform_ptr;

   const Gshhs& gshhs = *(andrea.get_gshhs_ptr (gshhs_identifier));

   bool is_fill = false;
   for (auto iterator = arguments.begin ();
        iterator != arguments.end (); iterator++)
   {
      const Tokens tokens (iterator->get_lower_case (), "=");
      const Dstring& option = tokens[0];
      const Dstring& value = tokens[1];
      if (option == "fill")
      {
         is_fill = (value == "yes" || value == "y" || value == "true");
      }
   }

   gshhs.cairo (cr, geodetic_transform);
   if (is_fill) { cr->fill (); } else { cr->stroke (); }

   delete geodetic_transform_ptr;

}

