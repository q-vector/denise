//
// graphics.h
// 
// Copyright (C) 2005-2013 Simon E. Ching
// 
// This file is part of libdenise.
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

#ifndef DENISE_GRAPHICS_H
#define DENISE_GRAPHICS_H

#include <cairomm/context.h>
#include <denise/basics.h>
#include <denise/cairoable.h>
#include <denise/geometry.h>

using namespace std;
using namespace Cairo;
using namespace denise;

namespace denise
{

   class Lat_Long;

   enum Color_Mode
   {
      ALPHA,
      HUE,
      SATURATION,
      BRIGHTNESS
   };

   enum Color_Gradient_Type
   {
      LINEAR,
      RADIAL
   };

   enum Dash_Mode
   {
      SOLID,
      DOTTED,
      DOT_DASHED,
      SHORT_DASHED,
      LONG_DASHED,
      DOT_DOT_DASHED,
      DOT_DOT_DOT_DASHED
   };

   enum Number_Type
   {
      NUMBER_REAL,
      NUMBER_TIME,
      NUMBER_LATITUDE,
      NUMBER_LONGITUDE
   };

   Real
   get_hue (const Integer i,
            const Integer n);

   RefPtr<ImageSurface>
   get_surface (const Size_2D& size_2d);

   RefPtr<Context>
   get_cr (const RefPtr<Surface> surface);

   class Dashes : public Tuple
   {

      public:

         Dashes (const string& str,
                 const string& delimiter = string (":"));

         Dashes (const Tuple& tuple);

         void
         cairo (const RefPtr<Context>& cr,
                const Real offset = 0.0) const;

   };

   class Paintable
   {

      public:

         virtual void
         cairo (const RefPtr<Context>& cr) const = 0;

   };

   /// The Color class is used to encapsulate colors in the default
   /// sRGB color space.  Each color component assumes a Real value
   /// in the range [0.0, 1.0].
   class Color : public Paintable
   {

      public:

         /// The red component
         Real
         r;

         /// The green component
         Real
         g;

         /// The blue component
         Real
         b;

         /// The alpha component
         Real
         a;

         /// Creates an sRGB color withthe specified values.
         ///
         /// If rgb is set true, the constructor expects HSV values
         /// instead of RGB values.
         Color (const Real r = 0,
                const Real g = 0,
                const Real b = 0,
                const Real a = 1.0);

         Color (const RefPtr<Context>& cr);

         static Color
         hsb (const Real h,
              const Real s,
              const Real b,
              const Real a = 1.0);

         static Color
         get_nac ();

         static Color
         transparent ();

         static Color
         black (const Real a = 1.0);

         static Color
         white (const Real a = 1.0);

         static Color
         gray (const Real b = 0.5,
               const Real a = 1.0);

         static Color
         red (const Real a = 1.0);

         static Color
         green (const Real a = 1.0);

         static Color
         blue (const Real a = 1.0);

         static Color
         cyan (const Real a = 1.0);

         static Color
         yellow (const Real a = 1.0);

         static Color
         magenta (const Real a = 1.0);

         Color&
         operator = (const Color& color);

         bool
         operator == (const Color& color) const;

         bool
         operator != (const Color& color) const;

         void
         set (const Real r,
              const Real g,
              const Real b,
              const Real a = 1.0);

         void
         set_hsb (const Real h,
                  const Real s,
                  const Real b,
                  const Real a = 1.0);

         /// Scales the brightness of this Color
         ///
         /// @param percentage  the percentage difference
         void
         scale_brightness (const Real percentage);

         /// Scales the saturation of this Color
         ///
         /// @param percentage  the percentage difference
         void
         scale_saturation (const Real percentage);

         /// Scales the opacity of this Color
         ///
         /// @param percentage  the percentage difference
         void
         scale_opacity (const Real percentage);

         /// Returns the hue value of this Color
         Real
         get_hue () const;

         /// Returns the saturation value of this Color
         Real
         get_saturation () const;

         /// Returns the brightness value of this Color
         Real
         get_brightness () const;

         /// Returns the brightness value of this Color
         Real
         get_luminosity () const;

         /// Returns a Color that is a brighter version of this Color.
         ///
         /// @param percentage  the percentage amount in brightening
         Color
         get_brighter_color (const Real percentage) const;

         /// Returns a Color that is a brighter version of this Color.
         ///
         /// @param percentage  the percentage amount in darkening
         Color
         get_darker_color (const Real percentage) const;

         /// Returns true if any of the component assums the value of nan.
         bool
         is_nac () const;

         void
         cairo (const RefPtr<Context>& cr) const;

   };

   const Color
   transparent (0, 0, 0, 0);

   class Color_Gradient : public Paintable
   {

      public:

         Color_Gradient_Type
         type;

         Color
         color_a;

         Color
         color_b;

         Point_2D
         point_a;

         Point_2D
         point_b;

         Real
         radius_a;

         Real
         radius_b;

         Color_Gradient (const Color& color_a,
                         const Color& color_b,
                         const Point_2D& point_a,
                         const Point_2D& point_b);

         Color_Gradient (const Color& color_a,
                         const Color& color_b,
                         const Point_2D& point_a,
                         const Point_2D& point_b,
                         const Real radius_a,
                         const Real radius_b);

   };

   class Motif : public Paintable
   {

      protected:

         RefPtr<SurfacePattern>
         pattern;

         virtual void
         init (const Paintable& paintable_0,
               const Paintable& paintable_1,
               const Real width) = 0;

      public:

         void
         cairo (const RefPtr<Context>& cr) const;

   };

   class Checkered : public Motif
   {

      private:

         void
         init (const Paintable& paintable_0,
               const Paintable& paintable_1,
               const Real width);

      public:

         Checkered (const Paintable& paintable_0,
                    const Paintable& paintable_1,
                    const Real width = 10);

   };

   class Stripped : public Motif
   {

      private:

         void
         init (const Paintable& paintable_0,
               const Paintable& paintable_1,
               const Real width);

      public:

         Stripped (const Paintable& paintable_0,
                   const Paintable& paintable_1,
                   const Real width = 10);

   };

   /// The Color_Chooser class maps from a Real value to a Color.
   ///
   /// This is an abstract base class that takes a Real input value
   /// and returns a Color.  This can be thought of a (Real -> Color)
   /// function.
   class Color_Chooser
   {

      protected:

         Color
         invalid_color;

      public:

         Color_Chooser (const Color& invalid_color = transparent);

         /// Returns a Color
         ///
         /// @param value  input value
         virtual Color
         get_color (const Real value) const = 0;

         void
         cairo (const RefPtr<Context>& cr,
                const Real value) const;

   };

   /// Constantly return a Color regardless of the input value.
   class Constant_Color_Chooser : public Color_Chooser
   {

      private:

         Color
         color;

      public:

         /// Creates a Constant_Color_Chooser which returns the given Color.
         ///
         /// @param color  Color that Constant_Color_Chooser returns.
         Constant_Color_Chooser (const Color& color,
                                 const Color& invalid_color = transparent);

         Color
         get_color (const Real value) const;

   };

   class Rgb_Color_Chooser : public Color_Chooser
   {

      public:

         class Component
         {

            private:

               const Rgb_Color_Chooser&
               rgb_color_chooser;

               Real
               start_fraction;

               Real
               fraction_range;

               Real
               gamma;

            public:

               Component (const Rgb_Color_Chooser& rcc,
                          const string& config_str);

               Real
               get_fraction (const Real value) const;

         };

      protected:

         const Real
         start_value;

         const Real
         end_value;

         const Real
         value_range;

         const Component
         component_r;

         const Component
         component_g;

         const Component
         component_b;

         const Component
         component_a;

         virtual Color
         color (const Real r,
                const Real g,
                const Real b,
                const Real a) const;

      public:

         Rgb_Color_Chooser (const Real start_value,
                            const Real end_value,
                            const string& str_r,
                            const string& str_g,
                            const string& str_b,
                            const string& str_a = "1",
                            const Color& invalid_color = transparent);

         Color
         get_color (const Real value) const;

         const Real
         get_start_value () const;

         const Real
         get_end_value () const;

         const Real
         get_value_range () const;

         bool
         within_range (const Real value) const;

   };

   class Hsb_Color_Chooser : public Rgb_Color_Chooser
   {


      protected:

         virtual Color
         color (const Real hue,
                const Real saturation,
                const Real brightness,
                const Real alpha) const;

      public:

         Hsb_Color_Chooser (const Real start_value,
                            const Real end_value,
                            const string& hue_str = "0",
                            const string& saturation_str = "1",
                            const string& brightness_str = "1",
                            const string& alpha_str = "1",
                            const Color& invalid_color = transparent);

   };

   class Hue_Color_Chooser : public Hsb_Color_Chooser
   {

      public:

         Hue_Color_Chooser (const Real start_value,
                            const Real end_value,
                            const string& hue_str,
                            const Real alpha = 1,
                            const Color& invalid_color = transparent);

         Hue_Color_Chooser (const Real start_value,
                            const Real end_value,
                            const Real saturation,
                            const Real brightness,
                            const Real alpha = 1,
                            const Real gamma = 1,
                            const Color& invalid_color = transparent);

         Hue_Color_Chooser (const Real start_value,
                            const Real end_value,
                            const Real start_hue,
                            const Real end_hue,
                            const Real saturation,
                            const Real brightness,
                            const Real alpha = 1,
                            const Real gamma = 1,
                            const Color& invalid_color = transparent);

   };

   class Gray_Color_Chooser : public Hsb_Color_Chooser
   {

      public:

         Gray_Color_Chooser (const Real start_value,
                             const Real end_value,
                             const Real gamma = 1,
                             const Color& invalid_color = transparent);

         Gray_Color_Chooser (const Real start_value,
                             const Real end_value,
                             const string& brightness_str,
                             const Color& invalid_color = transparent);

   };

   /// Mono_Color_Chooser is a subclass of Color_Chooser that returns
   /// a Color with either varying hue, saturation, or brightness.
   ///
   /// Creating Mono_Color_Chooser requires a specification of:
   /// \li end_value, start_value which specifies the value range
   /// \li color_mode which specifies where hue, saturation or brightness
   ///     would be varying
   /// \li default hue, saturation and brightness values
   class Mono_Color_Chooser : public Color_Chooser
   {

      private:

         Color_Mode
         color_mode;

         Real
         alpha;

         Real
         hue;

         Real
         saturation;

         Real
         brightness;

         Real
         start_value;

         Real
         end_value;

         Real
         value_range;

         Real
         start_extent;

         Real
         end_extent;

         Real
         extent_range;

      public:

         Mono_Color_Chooser (const Real end_value,
                             const Real start_value = 0,
                             const Color_Mode color_mode = SATURATION,
                             const Real alpha = 1,
                             const Real hue = 0,
                             const Real saturation = 1,
                             const Real brightness = 1,
                             const Color& invalid_color = transparent);

         Mono_Color_Chooser (const Real end_value,
                             const Real start_value,
                             const Real end_extent,
                             const Real start_extent,
                             const Color_Mode color_mode = SATURATION,
                             const Real alpha = 1,
                             const Real hue = 0,
                             const Real saturation = 1,
                             const Real brightness = 1,
                             const Color& invalid_color = transparent);

         Color
         get_color (const Real value) const;

   };

   /// Bi_Color_Chooser is a subclass of Color_Chooser that returns
   /// a Color with varying saturation. Colors with different signs
   /// are given different hues.
   ///
   /// Creating Bi_Color_Chooser requires a specification of:
   /// \li max_value which specifies the maximum absolute value
   /// \li positive_hue, negative_hue which specifies the hue for
   ///     positive or negative input values respectively
   /// \li max_saturation which specifies the maximum saturation
   /// \li default brightness
   class Bi_Color_Chooser : public Color_Chooser
   {

      private:

         Real
         positive_hue;

         Real
         negative_hue;

         Real
         alpha;

         Real
         brightness;

         Real
         max_saturation;

         Real
         max_value;

      public:

         Bi_Color_Chooser (const Real max_value,
                           const Real positive_hue = 0,
                           const Real negative_hue = 0.667,
                           const Real alpha = 1,
                           const Real brightness = 1,
                           const Real max_saturation = 1,
                           const Color& invalid_color = transparent);

         Color
         get_color (const Real value) const;

   };

   class Multi_Color_Chooser : public Color_Chooser,
                               public map<Real, Color>
   {

      private:

         bool
         staircase;

      public:

         Multi_Color_Chooser (const Real value_a,
                              const Real value_b,
                              const Color& color_a,
                              const Color& color_b,
                              const bool staircase = false);

         Multi_Color_Chooser (const Real value_a,
                              const Real value_b,
                              const Real value_c,
                              const Color& color_a,
                              const Color& color_b,
                              const Color& color_c,
                              const bool staircase = false);

         Multi_Color_Chooser (const Real value_a,
                              const Real value_b,
                              const Real value_c,
                              const Real value_d,
                              const Color& color_a,
                              const Color& color_b,
                              const Color& color_c,
                              const Color& color_d,
                              const bool staircase = false);

         Multi_Color_Chooser (const Real value_a,
                              const Real value_b,
                              const Real value_c,
                              const Real value_d,
                              const Real value_e,
                              const Color& color_a,
                              const Color& color_b,
                              const Color& color_c,
                              const Color& color_d,
                              const Color& color_e,
                              const bool staircase = false);

         void
         set_staircase (const bool staircase);

         Color
         get_color (const Real value) const;

   };

   class Multi_Rgb_Color_Chooser : public Color_Chooser,
                                   public list<Rgb_Color_Chooser>
   {

      public:

         Multi_Rgb_Color_Chooser (const Color& invalid_color = transparent);

         void
         add (const Rgb_Color_Chooser& ccc);

         Color
         get_color (const Real value) const;

   };

   class Level_And_Color
   {

      public:

         Real
         level;

         Color
         color;

         Level_And_Color (const Real level,
                          const Color& color);

   };

   class Custom_Color_Chooser : public Color_Chooser
   {

      private:

         const Color
         base_color;

         vector<Level_And_Color>
         level_and_color_vector;

      public:

         Custom_Color_Chooser ();

         Custom_Color_Chooser (const Color& base_color);

         void
         append (const Real level,
                 const Color& color);

         Color
         get_color (const Real value) const;

   };

   class Label : public Rect
   {

      private:

         string
         text;

         Point_2D
         point_2d;

         char
         justify_h;

         char
         justify_v;

         Point_2D
         offset;

         Real
         text_angle;

         Rect
         cairo_text (const RefPtr<Context>& cr,
                     const Point_2D& point,
                     const bool align,
                     const Real margin_x,
                     const Real margin_y,
                     const bool outline) const;

         static string
         get_string (const Real number,
                     const Number_Type number_type,
                     const string& format_str,
                     const Real multiplier = 1,
                     const Real offset = 0);
                     
         void
         set_offset (const char justify_h,
                     const char justify_v,
                     const Real padding);

      public:

         Label (const string& text,
                const Point_2D& point_2d,
                const char justify_h,
                const char justify_v);

         Label (const string& text,
                const Point_2D& point_2d,
                const char justify_h,
                const char justify_v,
                const Real padding);

         Label (const Real number,
                const Number_Type number_type,
                const string& format_str,
                const Point_2D& point_2d,
                const char justify_h,
                const char justify_v);

         Label (const Real number,
                const Number_Type number_type,
                const string& format_str,
                const Point_2D& point_2d,
                const char justify_h,
                const char justify_v,
                const Real padding);

         Label (const Real number,
                const Number_Type number_type,
                const string& format_str,
                const Real multiplier,
                const Real offset,
                const Point_2D& point_2d,
                const char justify_h,
                const char justify_v);

         Label (const Real number,
                const Number_Type number_type,
                const string& format_str,
                const Real multiplier,
                const Real offset,
                const Point_2D& point_2d,
                const char justify_h,
                const char justify_v,
                const Real padding);

         Label (const Label& label);

         void
         set_offset (const Point_2D& offset);

         void
         set_text_angle (const Real text_angle);

         void
         cairo (const RefPtr<Context>& cr,
                const bool align = false,
                const Real margin_x = 3,
                const Real margin_y = 3,
                const bool outline = false);

         void
         cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d,
                const bool align = false,
                const Real margin_x = 3,
                const Real margin_y = 3,
                const bool outline = false);

         void
         cairo (const RefPtr<Context>& cr,
                const Color& fg_color,
                const Color& bg_color,
                const Point_2D& offset,
                const bool align = false,
                const Real margin_x = 3,
                const Real margin_y = 3,
                const bool outline = false);

         void
         cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d,
                const Color& fg_color,
                const Color& bg_color,
                const Point_2D& offset,
                const bool align = false,
                const Real margin_x = 3,
                const Real margin_y = 3,
                const bool outline = false);

   };

   class Simple_Mesh_2D
   {

      private:

         Tuple
         get_tuple (const Domain_1D domain_1d,
                    const Real interval,
                    const Real multiplier,
                    const Real offset) const;

         Tuple
         get_tuple_x (const Domain_1D& domain_x) const;

         Tuple
         get_tuple_y (const Domain_1D& domain_y) const;

      public:

         Color
         color;

         Tuple
         tuple_x;

         Tuple
         tuple_y;

         Real
         interval_x;

         Real
         interval_y;

         Real
         multiplier_x;

         Real
         offset_x;

         Real
         multiplier_y;

         Real
         offset_y;

         Real
         line_width;

         Simple_Mesh_2D (const Color& color,
                         const Tuple& tuple_x,
                         const Tuple& tuple_y);

         Simple_Mesh_2D (const Color& color,
                         const Tuple& tuple_x,
                         const Real interval_y,
                         const Real multiplier_y = 1,
                         const Real offset_y = 0);

         Simple_Mesh_2D (const Color& color,
                         const Real interval_x,
                         const Tuple& tuple_y,
                         const Real multiplier_x = 1,
                         const Real offset_x = 0);

         Simple_Mesh_2D (const Color& color,
                         const Real interval_x,
                         const Real interval_y,
                         const Real multiplier_x = 1,
                         const Real offset_x = 0,
                         const Real multiplier_y = 1,
                         const Real offset_y = 0);

         void
         set_line_width (const Real line_width);

         void
         render (const RefPtr<Context>& cr,
                 const Transform_2D& transform,
                 const Size_2D& size_2d,
                 const Domain_2D& domain_2d) const;

         void
         render_label_x (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Domain_1D& domain_1d,
                         const Real position_y,
                         const string& format_str,
                         const Number_Type number_type = NUMBER_REAL,
                         const char justify_h = 'c',
                         const char justify_v = 'c',
                         const Real padding = 0) const;

         void
         render_label_y (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Domain_1D& domain_1d,
                         const Real position_x,
                         const string& format_str,
                         const Number_Type number_type = NUMBER_REAL,
                         const char justify_h = 'c',
                         const char justify_v = 'c',
                         const Real padding = 0) const;

         void
         render_label_lat_long (const RefPtr<Context>& cr,
                                const Transform_2D& transform,
                                const Domain_2D& domain_2d,
                                const Lat_Long& lat_long,
                                const string& format_str,
                                const char justify_h = 'c',
                                const char justify_v = 'c',
                                const Real padding = 0) const;

   };

   class Mesh_2D : protected vector<Simple_Mesh_2D>
   {

      protected:

         Size_2D
         size_2d;

         Domain_2D
         domain_2d;

      public:

         Mesh_2D ();

         Mesh_2D (const Size_2D& size_2d,
                  const Domain_2D& domain_2d);

         Mesh_2D (const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Simple_Mesh_2D& simple_mesh_2d);

         Mesh_2D (const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Simple_Mesh_2D& simple_mesh_2d_a,
                  const Simple_Mesh_2D& simple_mesh_2d_b);

         Mesh_2D (const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Simple_Mesh_2D& simple_mesh_2d_a,
                  const Simple_Mesh_2D& simple_mesh_2d_b,
                  const Simple_Mesh_2D& simple_mesh_2d_c);

         Mesh_2D (const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Simple_Mesh_2D& simple_mesh_2d_a,
                  const Simple_Mesh_2D& simple_mesh_2d_b,
                  const Simple_Mesh_2D& simple_mesh_2d_c,
                  const Simple_Mesh_2D& simple_mesh_2d_d);

         Mesh_2D (const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const vector<Simple_Mesh_2D>& simple_mesh_2d_vector);

         void
         set_size_2d (const Size_2D& size_2d);

         void
         set_domain_2d (const Domain_2D& domain_2d);

         void
         set_domain_2d (const Domain_1D& domain_x,
                        const Domain_1D& domain_y);

         void
         add (const Simple_Mesh_2D& simple_mesh_2d);

         const Domain_2D&
         get_domain_2d () const;

         void
         render (const RefPtr<Context>& cr,
                 const Transform_2D& transform) const;

         void
         render_label_x (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Integer index,
                         const Real position_y,
                         const string& format_str,
                         const Number_Type number_type = NUMBER_REAL,
                         const char justify_h = 'c',
                         const char justify_v = 'c',
                         const Real padding = 0) const;

         void
         render_label_y (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Integer index,
                         const Real position_x,
                         const string& format_str,
                         const Number_Type number_type = NUMBER_REAL,
                         const char justify_h = 'c',
                         const char justify_v = 'c',
                         const Real padding = 0) const;

         void
         render_label_lat_long (const RefPtr<Context>& cr,
                                const Transform_2D& transform,
                                const Integer index,
                                const Lat_Long& lat_long,
                                const string& format_str,
                                const char justify_h = 'c',
                                const char justify_v = 'c',
                                const Real padding = 0) const;

   };

   class Raster : public Box_2D
   {

      protected:

         uint32_t*
         buffer;

      public:

         Raster (const Size_2D& size_2d);

         Raster (const Box_2D& Box_2D);

         ~Raster ();

         void
         blend_pixel (const Integer i,
                      const Integer j,
                      const Color& color);

         void
         set_pixel (const Integer i,
                    const Integer j,
                    const Color& color);

         Color
         get_pixel (const Integer i,
                    const Integer j) const;

         RefPtr<ImageSurface>
         get_image_surface () const;

         void
         blit (const RefPtr<Context>& cr,
               const Real alpha = 1) const;

   };

   class Title : public Tokens,
                 private Size_2D
   {

      protected:

         Size_2D
         size_2d;

         Color
         bg_color;

         Color
         fg_color;

         Color
         shadow_color;

         virtual void
         cairo (const RefPtr<Context>& cr,
                const string& string_l,
                const string& string_c,
                const string& string_r) const;

         virtual void
         cairo (const RefPtr<Context>& cr,
                const string& string_ul,
                const string& string_ll,
                const string& string_c,
                const string& string_ur,
                const string& string_lr) const;

      public:

         Title (const Size_2D& size_2d,
                const Color& bg_color = Color (0, 0, 0, 0.5),
                const Color& fg_color = Color (1, 1, 1),
                const Color& shadow_color = Color (0, 0, 0, 0.25));

         void
         set_size_2d (const Size_2D& size_2d);

         Real
         get_height () const;

         void
         set (const string& string_l,
              const string& string_c,
              const string& string_r);

         void
         set (const string& string_ul,
              const string& string_ll,
              const string& string_c,
              const string& string_ur,
              const string& string_lr);

         void
         cairo (const RefPtr<Context>& cr);

   };

   ostream&
   operator << (ostream &out_file,
                const Color& color);

}

#endif /* DENISE_GRAPHICS_H */

