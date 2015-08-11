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

#include "andrea.h"
#include "geodetic_mesh.h"

using namespace std;
using namespace denise;
using namespace andrea;

Geodetic_Mesh_Package::Geodetic_Mesh_Package (Andrea& andrea)
   : Package (andrea)
{
}

void
Geodetic_Mesh_Package::geodetic_mesh_assign (const Dstring& identifier,
                                             const Size_2D& size_2d,
                                             const Domain_2D& domain_2d)
{
   const Geodetic_Mesh geodetic_mesh (size_2d, domain_2d);
   geodetic_mesh_map[identifier] = geodetic_mesh;
}

void
Geodetic_Mesh_Package::geodetic_mesh_add (const Dstring& identifier,
                                          const Tokens& arguments)
{

   const Integer n = arguments.size ();
   Geodetic_Mesh& geodetic_mesh = geodetic_mesh_map.at (identifier);

   switch (n)
   {
      case 2:
      {
         const Real interval = stof (arguments[0]); 
         const Color color (arguments[1]);
         geodetic_mesh.add (Simple_Mesh_2D (interval, color));
         break;
      }

      case 3:
      {
         const Real interval_lat = stof (arguments[0]); 
         const Real interval_long = stof (arguments[1]); 
         const Color color (arguments[2]);
         const Simple_Mesh_2D sm (interval_lat, interval_long, color);
         geodetic_mesh.add (sm);
         break;
      }

   } 

}

void
Geodetic_Mesh_Package::geodetic_mesh_print (const Dstring& identifier,
                                            const Tokens& arguments) const
{
   auto iterator = geodetic_mesh_map.find (identifier);
   const bool is_present = (iterator != geodetic_mesh_map.end ());
   if (is_present)
   {
      cout << "geodetic_mesh " << identifier << " is present" << endl;
   }
}

void
Geodetic_Mesh_Package::geodetic_mesh_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "assign")
   {
      const Dstring& identifier = tokens[1];
      const Size_2D size_2d (tokens[2]);
      const Domain_2D domain_2d (tokens[3]);
      geodetic_mesh_assign (identifier, size_2d, domain_2d);
   }
   else
   if (tokens[0] == "add")
   {
      const Dstring& identifier = tokens[1];
      geodetic_mesh_add (identifier, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "print")
   {
      const Dstring& identifier = tokens[1];
      geodetic_mesh_print (identifier, tokens.subtokens (2));
   }

}

const map<Dstring, Geodetic_Mesh>&
Geodetic_Mesh_Package::get_geodetic_mesh_map () const
{
   return geodetic_mesh_map;
}

const Geodetic_Mesh&
Geodetic_Mesh_Package::get_geodetic_mesh (const Dstring& identifier) const
{
   Exception e ("geodetic_mesh not found: " + identifier);
   try { return geodetic_mesh_map.at (identifier); }
   catch (const std::out_of_range& oor) { throw e; }
}

void
Geodetic_Mesh_Package::surface_geodetic_mesh (const Dstring& surface_identifier,
                                              const Dstring& geodetic_transform_identifier,
                                              const Dstring& geodetic_mesh_identifier)
{

   const RefPtr<Surface>& surface = andrea.get_surface (surface_identifier);
   const RefPtr<Context> cr = andrea.get_cr (surface_identifier);
   const Size_2D& size_2d = andrea.get_size_2d (surface_identifier);
   const Point_2D centre (Real (size_2d.i) / 2, Real (size_2d.j) / 2);

   const Geodetic_Transform* geodetic_transform_ptr =
      andrea.get_geodetic_transform_ptr (geodetic_transform_identifier, centre);
   const Geodetic_Transform& geodetic_transform = *geodetic_transform_ptr;

   const Geodetic_Mesh& geodetic_mesh =
      andrea.get_geodetic_mesh (geodetic_mesh_identifier);

   cr->save ();
   Color::black ().cairo (cr);
   geodetic_mesh.cairo (cr, geodetic_transform);
   cr->stroke ();
   cr->restore ();

   delete geodetic_transform_ptr;

}

