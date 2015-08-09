//
// surface.h
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

#ifndef DENISE_SRC_SURFACE_H
#define DENISE_SRC_SURFACE_H

#include <map>
#include <denise/dstring.h>
#include "package.h"

using namespace std;
using namespace denise;

namespace andrea
{

   class Surface_Package : public Package
   {

      protected:

         map<string, RefPtr<Surface> >
         surface_map;

         map<string, Size_2D>
         size_2d_map;

         map<string, string>
         extension_map;

         map<string, RefPtr<Context> >
         cr_map;

         Surface_Package (Andrea& andrea);

         void
         surface_init (const string& identifier,
                       const string& format,
                       const string& geometry);

         void
         surface_show (const string& identifier) const;


         void
         surface_finish (const string& identifier);

         void
         surface_write (const string& identifier,
                        const string& file_path) const;

         void
         surface_save (const string& identifier) const;

         void
         surface_restore (const string& identifier) const;

         void
         surface_stroke (const string& identifier) const;

         void
         surface_stroke_preserve (const string& identifier) const;

         void
         surface_fill (const string& identifier) const;

         void
         surface_fill_preserve (const string& identifier) const;

         void
         surface_clip (const string& identifier) const;

         void
         surface_clip_preserve (const string& identifier) const;

         void
         surface_color (const string& identifier,
                        const Color& color) const;

         void
         surface_line_width (const string& identifier,
                             const Real line_width) const;

         void
         surface_font_size (const string& identifier,
                            const Real font_size) const;

         void
         surface_font_face (const string& identifier,
                            const string& font_face) const;

         void
         surface_paint (const string& identifier) const;

         void
         surface_edge (const string& identifier,
                       const Tokens& arguments) const;

         void
         surface_circle (const string& identifier,
                         const Tokens& arguments) const;

         void
         surface_ellipse (const string& identifier,
                          const Tokens& arguments) const;

         void
         surface_label (const string& identifier,
                        const Tokens& arguments) const;

         void
         surface_title (const string& identifier,
                        const Tokens& tokens) const;

         void
         surface_range_circle (const string& surface_identifier,
                               const string& geodetic_transform_identifier,
                               const Lat_Long& lat_long,
                               const Real distance) const;

         void
         surface_parse (const Tokens& tokens);

      public:

         const map<string, RefPtr<Surface> >&
         get_surface_map () const;

         const map<string, Size_2D>&
         get_size_2d_map () const;

         const map<string, string>&
         get_extension_map () const;

         const map<string, RefPtr<Context> >&
         get_cr_map () const;

         const RefPtr<Surface>&
         get_surface (const string& identifier) const;

         const Size_2D&
         get_size_2d (const string& identifier) const;

         const string&
         get_extension (const string& identifier) const;

         const RefPtr<Context>&
         get_cr (const string& identifier) const;

   };

};

#endif /* DENISE_SRC_SURFACE_H */

