//
// andrea.cc
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

#include <cairomm/xlib_surface.h>
#include <X11/Xlib.h>
#include "package.h"
#include "andrea.h"
#include "geodetic_transform.h"
#include "gshhs.h"
#include "surface.h"

using namespace std;
using namespace denise;
using namespace andrea;

Surface_Package::Surface_Package (Andrea& andrea)
   : Package (andrea)
{
}

void
Surface_Package::surface_init (const string& identifier,
                               const string& file_path,
                               const string& geometry)
{

   const Tokens geometry_tokens (geometry, "x");
   const Integer ww = stof (geometry_tokens[0]);
   const Integer hh = stof (geometry_tokens[1]);
   const Size_2D size_2d (ww, hh);

   const Tokens file_path_tokens (file_path, ".");
   const string& extension = file_path_tokens.back ();
   const string& e = get_lower_case (extension);
   const RefPtr<Surface> surface = denise::get_surface (size_2d, e, file_path);
   const RefPtr<Context> cr = denise::get_cr (surface);

   surface_map[identifier] = surface;
   size_2d_map[identifier] = size_2d;
   extension_map[identifier] = e;
   cr_map[identifier] = cr;

}

void
Surface_Package::surface_show (const string& identifier) const
{

   const RefPtr<Surface>& surface = surface_map.at (identifier);
   const Size_2D& size_2d = size_2d_map.at (identifier);
   const Integer ww = size_2d.i;
   const Integer hh = size_2d.j;

   XWindowAttributes xwindow_attributes_return;
   Display* d = XOpenDisplay (NULL);
   if (d == NULL)
   {
      cerr << "Cannot open display" << endl;
      return;
   }

   int s = DefaultScreen (d);

   Window w = XCreateSimpleWindow (d, RootWindow (d, s), 10, 10,
      ww, hh, 1, BlackPixel (d, s), WhitePixel (d, s));
   
   XSelectInput (d, w, ExposureMask | KeyPressMask | ButtonPressMask);
   XMapWindow (d, w);

   XGetWindowAttributes (d, w, &xwindow_attributes_return);
   Visual* v = xwindow_attributes_return.visual;
   RefPtr<XlibSurface> xlib_surface = XlibSurface::create (d, w, v, ww, hh);
   RefPtr<Context> cr = Context::create (xlib_surface);

   Atom WM_DELETE_WINDOW = XInternAtom (d, "WM_DELETE_WINDOW", False); 
   XSetWMProtocols (d, w, &WM_DELETE_WINDOW, 1);

   for (bool done = false; !done; )
   {

      XEvent e;
      XNextEvent (d, &e);

      switch (e.type)
      {

         case Expose:
         {
            if (surface != 0)
            {
               Checkered (Color::gray (0.55), Color::gray (0.45)).cairo (cr);
               cr->paint();
               cr->set_source (surface, 0, 0);
               cr->paint ();
            }
            break;
         }

         case KeyPress:
         case ClientMessage:
         {
            done = true;
            break;
         }

      }

    }

    
    XCloseDisplay (d);

}

void
Surface_Package::surface_finish (const string& identifier)
{
   surface_map.erase (identifier);
   size_2d_map.erase (identifier);
   extension_map.erase (identifier);
   cr_map.erase (identifier);
}

void
Surface_Package::surface_write (const string& identifier,
                                const string& file_path) const
{
   const RefPtr<Surface>& surface = andrea.get_surface (identifier);
   const string& e = andrea.get_extension (identifier);
   surface->write_to_png (file_path);
}

void
Surface_Package::surface_save (const string& identifier) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->save ();
}

void
Surface_Package::surface_restore (const string& identifier) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->restore ();
}

void
Surface_Package::surface_stroke (const string& identifier) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->stroke ();
}

void
Surface_Package::surface_stroke_preserve (const string& identifier) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->stroke_preserve ();
}

void
Surface_Package::surface_fill (const string& identifier) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->fill ();
}

void
Surface_Package::surface_fill_preserve (const string& identifier) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->fill_preserve ();
}

void
Surface_Package::surface_clip (const string& identifier) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->clip ();
}

void
Surface_Package::surface_clip_preserve (const string& identifier) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->clip_preserve ();
}

void
Surface_Package::surface_color (const string& identifier,
                                const Color& color) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   color.cairo (cr);
}

void
Surface_Package::surface_line_width (const string& identifier,
                                     const Real line_width) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->set_line_width (line_width);
}

void
Surface_Package::surface_font_size (const string& identifier,
                                    const Real font_size) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->set_font_size (font_size);
}

void
Surface_Package::surface_font_face (const string& identifier,
                                    const string& font_face) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->select_font_face (font_face, FONT_SLANT_NORMAL, FONT_WEIGHT_NORMAL);

}

void
Surface_Package::surface_paint (const string& identifier) const
{
   const RefPtr<Context> cr = get_cr (identifier);
   cr->paint ();
}

void
Surface_Package::surface_edge (const string& identifier,
                               const Tokens& arguments) const
{

   const RefPtr<Context> cr = get_cr (identifier);
   const Point_2D point_a (arguments[0]);
   const Point_2D point_b (arguments[1]);
   const Edge edge (point_a, point_b);

   edge.cairo (cr);

}

void
Surface_Package::surface_circle (const string& identifier,
                                 const Tokens& arguments) const
{

   const RefPtr<Context> cr = get_cr (identifier);
   const Point_2D centre (arguments[0]);
   const Real radius = stof (arguments[1]);
   const Circle circle (centre, radius);

   circle.cairo (cr);

}

void
Surface_Package::surface_ellipse (const string& identifier,
                                  const Tokens& arguments) const
{

   const RefPtr<Context> cr = get_cr (identifier);
   const Point_2D centre (arguments[0]);
   const Real major_radius = stof (arguments[1]);
   const Real minor_radius = stof (arguments[2]);
   const Real slant = stof (arguments[3]);
   const Ellipse ellipse (centre, major_radius, minor_radius, slant);

   ellipse.cairo (cr);

}

void
Surface_Package::surface_label (const string& identifier,
                                const Tokens& arguments) const
{

   const RefPtr<Context> cr = get_cr (identifier);
   const string& text = arguments[0];
   const Point_2D point (arguments[1]);
   const char justify_h = arguments[2][0];
   const char justify_v = arguments[3][0];
   const Real padding = stof (arguments[4]);
   Label label (text, point, justify_h, justify_v, padding);

   label.cairo (cr);

}

void
Surface_Package::surface_title (const string& identifier,
                                const Tokens& tokens) const
{

   const RefPtr<Surface>& surface = andrea.get_surface (identifier);
   const RefPtr<Context> cr = andrea.get_cr (identifier);
   const Size_2D& size_2d = andrea.get_size_2d (identifier);

   Title title (size_2d);
   title.set (tokens);
   title.cairo (cr);

}

void
Surface_Package::surface_range_circle (const string& surface_identifier,
                                       const string& geodetic_transform_identifier,
                                       const Lat_Long& lat_long,
                                       const Real distance) const
{

   const RefPtr<Surface>& surface = andrea.get_surface (surface_identifier);
   const RefPtr<Context> cr = andrea.get_cr (surface_identifier);
   const Size_2D& size_2d = andrea.get_size_2d (surface_identifier);
   const Point_2D centre (Real (size_2d.i) / 2, Real (size_2d.j) / 2);

   const Geodetic_Transform* geodetic_transform_ptr =
      andrea.get_geodetic_transform_ptr (geodetic_transform_identifier, centre);
   const Geodetic_Transform& geodetic_transform = *geodetic_transform_ptr;

cout << distance << " " << lat_long << endl;
   Range_Circle range_circle (lat_long, distance);
   range_circle.cairo (cr, geodetic_transform);
   cr->fill ();

   delete geodetic_transform_ptr;

}

void
Surface_Package::surface_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "init")
   {
      const string& identifier = tokens[1];
      const string& file_path = tokens[2];
      const string& geometry = tokens[3];
      surface_init (identifier, file_path, geometry);
   }
   else
   if (tokens[0] == "show")
   {
      const string& identifier = tokens[1];
      surface_show (identifier);
   }
   else
   if (tokens[0] == "finish")
   {
      const string& identifier = tokens[1];
      surface_finish (identifier);
   }
   else
   if (tokens[0] == "write")
   {
      const string& identifier = tokens[1];
      const string& file_path = tokens[2];
      surface_write (identifier, file_path);
   }
   else
   if (tokens[0] == "save")
   {
      const string& identifier = tokens[1];
      surface_save (identifier);
   }
   else
   if (tokens[0] == "restore")
   {
      const string& identifier = tokens[1];
      surface_restore (identifier);
   }
   else
   if (tokens[0] == "stroke")
   {
      const string& identifier = tokens[1];
      surface_stroke (identifier);
   }
   else
   if (tokens[0] == "stroke_preserve")
   {
      const string& identifier = tokens[1];
      surface_stroke_preserve (identifier);
   }
   else
   if (tokens[0] == "fill")
   {
      const string& identifier = tokens[1];
      surface_fill (identifier);
   }
   else
   if (tokens[0] == "fill_preserve")
   {
      const string& identifier = tokens[1];
      surface_fill_preserve (identifier);
   }
   else
   if (tokens[0] == "clip")
   {
      const string& identifier = tokens[1];
      surface_clip (identifier);
   }
   else
   if (tokens[0] == "clip_preserve")
   {
      const string& identifier = tokens[1];
      surface_clip_preserve (identifier);
   }
   else
   if (tokens[0] == "color")
   {
      const string& identifier = tokens[1];
      const Color& color (tokens[2]);
      surface_color (identifier, color);
   }
   else
   if (tokens[0] == "line_width")
   {
      const string& identifier = tokens[1];
      const Real line_width = stof (tokens[2]);
      surface_line_width (identifier, line_width);
   }
   else
   if (tokens[0] == "font_size")
   {
      const string& identifier = tokens[1];
      const Real font_size = stof (tokens[2]);
      surface_font_size (identifier, font_size);
   }
   else
   if (tokens[0] == "font_face")
   {
      const string& identifier = tokens[1];
      const string& font_face = tokens[2];
      surface_font_face (identifier, font_face);
   }
   else
   if (tokens[0] == "paint")
   {
      const string& identifier = tokens[1];
      surface_paint (identifier);
   }
   else
   if (tokens[0] == "title")
   {
      const string& surface_identifier = tokens[1];
      surface_title (surface_identifier, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "paint")
   {
      const string& surface_identifier = tokens[1];
      surface_paint (surface_identifier);
   }
   else
   if (tokens[0] == "edge")
   {
      const string& surface_identifier = tokens[1];
      surface_edge (surface_identifier, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "circle")
   {
      const string& surface_identifier = tokens[1];
      surface_circle (surface_identifier, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "ellipse")
   {
      const string& surface_identifier = tokens[1];
      surface_ellipse (surface_identifier, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "label")
   {
      const string& surface_identifier = tokens[1];
      surface_label (surface_identifier, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "range_circle")
   {
      const string& surface_identifier = tokens[1];
      const string& geodetic_transform_identifier = tokens[2];
      const Lat_Long& lat_long (tokens[3]);
      const Real distance = stof (tokens[4]);
      surface_range_circle (surface_identifier,
         geodetic_transform_identifier, lat_long, distance);
   }
   else
   if (tokens[0] == "sounding")
   {
      andrea.surface_sounding (tokens.subtokens (1));
   }
   else
   if (tokens[0] == "journey")
   {
      const string& surface_identifier = tokens[1];
      const string& geodetic_transform_identifier = tokens[2];
      const string& journey_identifier = tokens[3];
      andrea.surface_journey (surface_identifier,
         geodetic_transform_identifier, journey_identifier);
   }
   else
   if (tokens[0] == "geodetic_mesh")
   {
      const string& surface_identifier = tokens[1];
      const string& geodetic_transform_identifier = tokens[2];
      const string& geodetic_mesh_identifier = tokens[3];
      andrea.surface_geodetic_mesh (surface_identifier,
         geodetic_transform_identifier, geodetic_mesh_identifier);
   }
   else
   if (tokens[0] == "gshhs")
   {
      const string& surface_identifier = tokens[1];
      const string& geodetic_transform_identifier = tokens[2];
      const string& gshhs_identifier = tokens[3];
      andrea.surface_gshhs (surface_identifier,
         geodetic_transform_identifier, gshhs_identifier,
         tokens.subtokens (4));
   }

}

const map<string, RefPtr<Surface> >&
Surface_Package::get_surface_map () const
{
   return surface_map;
}

const map<string, Size_2D>&
Surface_Package::get_size_2d_map () const
{
   return size_2d_map;
}

const map<string, string>&
Surface_Package::get_extension_map () const
{
   return extension_map;
}

const map<string, RefPtr<Context> >&
Surface_Package::get_cr_map () const
{
   return cr_map;
}

const RefPtr<Surface>&
Surface_Package::get_surface (const string& identifier) const
{
   Exception e ("surface not found: " + identifier);
   try { return surface_map.at (identifier); }
   catch (const std::out_of_range& oor) { throw e; }
}

const Size_2D&
Surface_Package::get_size_2d (const string& identifier) const
{
   Exception e ("surface not found: " + identifier);
   try { return size_2d_map.at (identifier); }
   catch (const std::out_of_range& oor) { throw e; }
}

const string&
Surface_Package::get_extension (const string& identifier) const
{
   Exception e ("surface not found: " + identifier);
   try { return extension_map.at (identifier); }
   catch (const std::out_of_range& oor) { throw e; }
}

const RefPtr<Context>&
Surface_Package::get_cr (const string& identifier) const
{
   Exception e ("surface not found: " + identifier);
   try { return cr_map.at (identifier); }
   catch (const std::out_of_range& oor) { throw e; }
}

