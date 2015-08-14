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

         map<Dstring, RefPtr<Surface> >
         surface_map;

         map<Dstring, Size_2D>
         size_2d_map;

         map<Dstring, Dstring>
         extension_map;

         map<Dstring, RefPtr<Context> >
         cr_map;

         Surface_Package (Andrea& andrea);

         void
         surface_init (const Dstring& identifier,
                       const Dstring& format,
                       const Dstring& geometry);

         void
         surface_show (const Dstring& identifier) const;


         void
         surface_finish (const Dstring& identifier);

         void
         surface_write (const Dstring& identifier,
                        const Dstring& file_path) const;

         void
         surface_save (const Dstring& identifier) const;

         void
         surface_restore (const Dstring& identifier) const;

         void
         surface_stroke (const Dstring& identifier) const;

         void
         surface_stroke_preserve (const Dstring& identifier) const;

         void
         surface_fill (const Dstring& identifier) const;

         void
         surface_fill_preserve (const Dstring& identifier) const;

         void
         surface_clip (const Dstring& identifier) const;

         void
         surface_clip_preserve (const Dstring& identifier) const;

         void
         surface_color (const Dstring& identifier,
                        const Color& color) const;

         void
         surface_line_width (const Dstring& identifier,
                             const Real line_width) const;

         void
         surface_font_size (const Dstring& identifier,
                            const Real font_size) const;

         void
         surface_font_face (const Dstring& identifier,
                            const Dstring& font_face) const;

         void
         surface_paint (const Dstring& identifier) const;

         void
         surface_edge (const Dstring& identifier,
                       const Tokens& arguments) const;

         void
         surface_circle (const Dstring& identifier,
                         const Tokens& arguments) const;

         void
         surface_ellipse (const Dstring& identifier,
                          const Tokens& arguments) const;

         void
         surface_label (const Dstring& identifier,
                        const Tokens& arguments) const;

         void
         surface_title (const Dstring& identifier,
                        const Tokens& tokens) const;

         void
         surface_range_circle (const Dstring& surface_identifier,
                               const Dstring& geodetic_transform_identifier,
                               const Lat_Long& lat_long,
                               const Real distance) const;

         void
         surface_domain (const Dstring& surface_identifier,
                         const Dstring& geodetic_transform_identifier,
                         const Domain_2D& domain) const;

         void
         surface_gtopo30 (const Dstring& surface_identifier,
                          const Dstring& geodetic_transform_identifier,
                          const Tokens& tokens) const;

         void
         surface_blue_marble (const Dstring& surface_identifier,
                              const Dstring& geodetic_transform_identifier) const;

         void
         surface_parse (const Tokens& tokens);

      public:

         const map<Dstring, RefPtr<Surface> >&
         get_surface_map () const;

         const map<Dstring, Size_2D>&
         get_size_2d_map () const;

         const map<Dstring, Dstring>&
         get_extension_map () const;

         const map<Dstring, RefPtr<Context> >&
         get_cr_map () const;

         const RefPtr<Surface>&
         get_surface (const Dstring& identifier) const;

         const Size_2D&
         get_size_2d (const Dstring& identifier) const;

         const Dstring&
         get_extension (const Dstring& identifier) const;

         const RefPtr<Context>&
         get_cr (const Dstring& identifier) const;

   };

};

#endif /* DENISE_SRC_SURFACE_H */

