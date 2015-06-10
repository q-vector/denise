//
// graphics.cc
// 
// Copyright (C) 2010 Simon E. Ching
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

#include <denise/geodesy.h>
#include <denise/graphics.h>

namespace denise
{

   Real
   get_hue (const Integer i,
            const Integer n = 3)
   {

      if (i < n) { return Real (i) / Real (n); }

      Integer nn;
      for (nn = n; i >= nn; nn += nn);

      Integer ii = i - nn/2;
      return Real (ii + ii + 1) / Real (nn);

   }

}

Dashes::Dashes (const string& str,
                const string& delimiter)
   : Tuple (str, delimiter)
{
}

Dashes::Dashes (const Tuple& tuple)
   : Tuple (tuple)
{
}

void
Dashes::cairo (const RefPtr<Context>& cr,
               const Real offset) const
{
   std::vector<double> d;
   for (Tuple::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const double x = *(iterator);
      d.push_back (x);
   }
   cr->set_dash (d, offset);
}

Color::Color (const Real r,
              const Real g,
              const Real b,
              const Real a)
         : r (r),
           g (g),
           b (b),
           a (a)
{
}

Color::Color (const RefPtr<Context>& cr)
{

   double r = GSL_NAN, g = GSL_NAN, b = GSL_NAN, a = GSL_NAN;
   RefPtr<Cairo::Pattern> pattern = cr->get_source ();

   if (pattern->get_type () == PATTERN_TYPE_SOLID)
   {
      RefPtr<SolidPattern>::cast_dynamic (pattern)->get_rgba (r, g, b, a);
   } 

   this->r = r;
   this->g = g;
   this->b = b;
   this->a = a;

}

Color
Color::hsb (const Real h,
            const Real s,
            const Real b,
            const Real a)
{
   Color color;
   color.set_hsb (h, s, b, a);
   return color;
}

Color
Color::get_nac ()
{
   return Color (GSL_NAN, GSL_NAN, GSL_NAN, GSL_NAN);
}

Color
Color::black (const Real a)
{
   return Color (0, 0, 0, a);
}

Color
Color::white (const Real a)
{
   return Color (1, 1, 1, a);
}

Color
Color::gray (const Real saturation,
             const Real a)
{
   return Color (saturation, saturation, saturation, a);
}

Color&
Color::operator = (const Color& color)
{
   r = color.r;
   g = color.g;
   b = color.b;
   a = color.a;
   return *this;
}

bool
Color::operator == (const Color& color) const
{
   return (is_nac () && color.is_nac ()) || ((r == color.r) &&
      (g == color.g) && (b == color.b) && (a == color.a)) ;
}

bool
Color::operator != (const Color& color) const
{
   return (!is_nac () || !color.is_nac ()) && ((r != color.r) ||
      (g != color.g) || (b != color.b) || (a != color.a));
}

void
Color::set (const Real r,
            const Real g,
            const Real b,
            const Real a)
{
   this->a = a;
   this->r = r;
   this->g = g;
   this->b = b;
}

void
Color::set_hsb (const Real h,
                const Real s,
                const Real b,
                const Real a)
{

   this->a = a;

   if (s == 0)
   {
      this->r = b;
      this->g = b;
      this->b = b;
   }
   else
   {

      Real hh = (h == 1 ? 0 : h*6);
      Integer i = Integer (floor (hh));

      Real f = hh - i;
      Real x = b * (1 - s);
      Real y = b * (1 - s * f);
      Real z = b * (1 - s * (1 - f));

      switch (i)
      {
         case 0: this->r = b; this->g = z; this->b = x; break;
         case 1: this->r = y; this->g = b; this->b = x; break;
         case 2: this->r = x; this->g = b; this->b = z; break;
         case 3: this->r = x; this->g = y; this->b = b; break;
         case 4: this->r = z; this->g = x; this->b = b; break;
         case 5: this->r = b; this->g = x; this->b = y; break;
      }

   }

}

void
Color::scale_brightness (const Real percentage)
{
   Real scaled_brightness = get_brightness () * percentage;
   if (scaled_brightness > 1) { scaled_brightness = 1; }
   if (scaled_brightness < 0) { scaled_brightness = 0; }

   set_hsb (get_hue (), get_saturation (), scaled_brightness, a);
}

void
Color::scale_saturation (const Real percentage)
{
   Real scaled_saturation = get_saturation () * percentage;
   if (scaled_saturation > 1) { scaled_saturation = 1; }
   if (scaled_saturation < 0) { scaled_saturation = 0; }

   set_hsb (get_hue (), scaled_saturation, get_brightness (), a);
}

void
Color::scale_opacity (const Real percentage)
{
   Real scaled_opacity = a * percentage;
   if (scaled_opacity > 1) { scaled_opacity = 1; }
   if (scaled_opacity < 0) { scaled_opacity = 0; }

   set_hsb (get_hue (), get_saturation (), get_brightness (), scaled_opacity);
}

Real
Color::get_hue () const
{

   Real A = ((r - g) + (r - b)) / 2;
   Real B = (r - g) * (r - g);
   Real C = (r - b) * (g - b);
   Real hue = acos (A / sqrt (B + C));
   Real brightness = get_brightness ();

   if (b / brightness > g / brightness) { hue = M_2_TIMES_PI - hue; }
   return hue / (M_2_TIMES_PI);

}

Real
Color::get_saturation () const
{

   Real min;
   min = (r < g ? r : g);
   min = (min < b ? min : b);

   return (1 - 3 / (r + g + b) * min);

}

Real
Color::get_brightness () const
{
   return max (r, max (g, b));
}

Real
Color::get_luminosity () const
{
   return r * 0.299 + g * 0.587 + b * 0.144;
}

Color
Color::get_brighter_color (const Real percentage) const
{
   Real brightness = get_brightness () * (1 + percentage / 100);
   if (brightness > 1) { brightness = 1; }
   return Color::hsb (get_hue (), get_saturation (), brightness, a);
}

Color
Color::get_darker_color (const Real percentage) const
{
   Real brightness = get_brightness () * (1 - percentage / 100);
   if (brightness < 0) { brightness = 0; }
   return Color::hsb (get_hue (), get_saturation (), brightness, a);
}

bool
Color::is_nac () const
{
   return gsl_isnan (r) || gsl_isnan (g) || gsl_isnan (b) || gsl_isnan (a);
}

void
Color::cairo (const RefPtr<Context>& cr) const
{
   cr->set_source_rgba (r, g, b, a);
}

Color_Gradient::Color_Gradient (const Color& color_a,
                                const Color& color_b,
                                const Point_2D& point_a,
                                const Point_2D& point_b)
   : type (LINEAR),
     color_a (color_a),
     color_b (color_b),
     point_a (point_a),
     point_b (point_b)
{
}

Color_Gradient::Color_Gradient (const Color& color_a,
                                const Color& color_b,
                                const Point_2D& point_a,
                                const Point_2D& point_b,
                                const Real radius_a,
                                const Real radius_b)
   : type (RADIAL),
     color_a (color_a),
     color_b (color_b),
     point_a (point_a),
     point_b (point_b),
     radius_a (radius_a),
     radius_b (radius_b)
{
}

void
Checkered::init (const Paintable& paintable_0,
                 const Paintable& paintable_1,
                 const Real width)
{

   const Integer w = Integer (round (width));
   const Integer w2 = Integer (round (width * 2));

   RefPtr<ImageSurface> surface = ImageSurface::create (FORMAT_ARGB32, w2, w2);
   RefPtr<Context> cr = Context::create (surface);
   cr->select_font_face ("Verdana", FONT_SLANT_NORMAL, FONT_WEIGHT_NORMAL);

   Color (1, 1, 1, 0.0).cairo (cr);
   cr->paint ();

   paintable_0.cairo (cr);
   cr->move_to (0, 0);
   cr->line_to (w, 0);
   cr->line_to (w, w2);
   cr->line_to (w2, w2);
   cr->line_to (w2, w);
   cr->line_to (0, w);
   cr->fill ();

   paintable_1.cairo (cr);
   cr->move_to (w, 0);
   cr->line_to (w2, 0);
   cr->line_to (w2, w);
   cr->line_to (0, w);
   cr->line_to (0, w2);
   cr->line_to (w, w2);
   cr->fill ();

   pattern = SurfacePattern::create (surface);
   pattern->set_extend (EXTEND_REPEAT);

}

Checkered::Checkered (const Paintable& paintable_0,
                      const Paintable& paintable_1,
                      const Real width)
{
   init (paintable_0, paintable_1, width);
}

void
Stripped::init (const Paintable& paintable_0,
                const Paintable& paintable_1,
                const Real width)
{

   const Integer w = Integer (round (width));
   const Integer w2 = Integer (round (width * 2));

   RefPtr<ImageSurface> surface = ImageSurface::create (FORMAT_ARGB32, w2, w2);
   RefPtr<Context> cr = Context::create (surface);
   cr->select_font_face ("Verdana", FONT_SLANT_NORMAL, FONT_WEIGHT_NORMAL);

   Color (1, 1, 1, 0.0).cairo (cr);
   cr->paint ();

   paintable_0.cairo (cr);
   cr->move_to (0, 0);
   cr->line_to (w, 0);
   cr->line_to (w2, w);
   cr->line_to (w2, w2);
   cr->fill ();
   cr->move_to (0, w);
   cr->line_to (w, w2);
   cr->line_to (0, w2);
   cr->fill ();

   paintable_1.cairo (cr);
   cr->move_to (0, 0);
   cr->line_to (w2, w2);
   cr->line_to (w, w2);
   cr->line_to (0, w);
   cr->fill ();
   cr->move_to (w, 0);
   cr->line_to (w2, 0);
   cr->line_to (w2, w);
   cr->fill ();

   pattern = SurfacePattern::create (surface);
   pattern->set_extend (EXTEND_REPEAT);

}

Stripped::Stripped (const Paintable& paintable_0,
                    const Paintable& paintable_1,
                    const Real width)
{
   init (paintable_0, paintable_1, width);
}

void
Motif::cairo (const RefPtr<Context>& cr) const
{
   cr->set_source (pattern);
}

Color_Chooser::Color_Chooser (const Color& invalid_color)
             : invalid_color (invalid_color)
{
}

void
Color_Chooser::cairo (const RefPtr<Context>& cr,
                      const Real value) const
{
   const Color& color = get_color (value);
   color.cairo (cr);
}

Constant_Color_Chooser::Constant_Color_Chooser (const Color& color,
                                                const Color& invalid_color)
                               : Color_Chooser (invalid_color),
                                         color (color)
{
}

Color
Constant_Color_Chooser::get_color (Real value) const
{
   return color;
}

Rgb_Color_Chooser::Component::Component (const Rgb_Color_Chooser& rcc,
                                         const string& config_str)
   : rgb_color_chooser (rcc)
{

   const Tokens tokens (config_str, ":");
   const Integer n = tokens.size ();

   start_fraction = tokens.real (0);
   fraction_range = ((n > 1) ? (tokens.real (1) - tokens.real (0)) : 0);
   this->gamma = ((n > 2) ? tokens.real (2) : 1);

}

Real
Rgb_Color_Chooser::Component::get_fraction (const Real value) const
{

   Real f;
   const Real& f_range = fraction_range;
   const Real& start_f = start_fraction;

   if (fabs (f_range) < 0.0001)
   {
      return start_f + f_range / 2;
   }
   else
   {
      const Real& start_v = rgb_color_chooser.start_value;
      const Real& v_range = rgb_color_chooser.value_range;
      Real delta_f = (((value - start_v) / v_range) * f_range);
      if (f_range > 0)
      {
         if (delta_f < 0) { delta_f = 0; }
         if (delta_f > f_range) { delta_f = f_range; }
      }
      else
      {
         if (delta_f > 0) { delta_f = 0; }
         if (delta_f < f_range) { delta_f = f_range; }
      }
      f = start_f + delta_f;
   }

   f = pow (f, gamma);
   if (f > 0.9999) { return 0.9999; }
   if (f < 0.0001) { return 0.0001; }
   return f;

}

Color
Rgb_Color_Chooser::color (const Real r,
                          const Real g,
                          const Real b,
                          const Real a) const
{
   return Color (r, g, b, a);
}

Rgb_Color_Chooser::Rgb_Color_Chooser (const Real start_value,
                                      const Real end_value,
                                      const string& str_r,
                                      const string& str_g,
                                      const string& str_b,
                                      const string& str_a,
                                      const Color& invalid_color)
   : Color_Chooser (invalid_color),
     start_value (start_value),
     end_value (end_value),
     value_range (end_value - start_value),
     component_r (*this, str_r),
     component_g (*this, str_g),
     component_b (*this, str_b),
     component_a (*this, str_a)
{
}

Color
Rgb_Color_Chooser::get_color (const Real value) const
{

   if (!gsl_finite (value)) { return invalid_color; }

   const Real r = component_r.get_fraction (value);
   const Real g = component_g.get_fraction (value);
   const Real b = component_b.get_fraction (value);
   const Real a = component_a.get_fraction (value);

   return color (r, g, b, a);

}

const Real
Rgb_Color_Chooser::get_start_value () const
{
   return start_value;
}

const Real
Rgb_Color_Chooser::get_end_value () const
{
   return end_value;
}

const Real
Rgb_Color_Chooser::get_value_range () const
{
   return value_range;
}

bool
Rgb_Color_Chooser::within_range (const Real value) const
{
   return (((value - start_value) * (value - end_value)) <= 0);
}

Color
Hsb_Color_Chooser::color (const Real hue,
                          const Real saturation,
                          const Real brightness,
                          const Real alpha) const
{
   return Color::hsb (hue, saturation, brightness, alpha);
}

Hsb_Color_Chooser::Hsb_Color_Chooser (const Real start_value,
                                      const Real end_value,
                                      const string& hue_str,
                                      const string& saturation_str,
                                      const string& brightness_str,
                                      const string& alpha_str,
                                      const Color& invalid_color)
   : Rgb_Color_Chooser (start_value, end_value, hue_str,
        saturation_str, brightness_str, alpha_str, invalid_color)
{
}

Hue_Color_Chooser::Hue_Color_Chooser (const Real start_value,
                                      const Real end_value,
                                      const string& hue_str,
                                      const Real alpha,
                                      const Color& invalid_color)
   : Hsb_Color_Chooser (start_value, end_value, hue_str, "1",
        "1", string_render ("%f", alpha), invalid_color)
{
}

Hue_Color_Chooser::Hue_Color_Chooser (const Real start_value,
                                      const Real end_value,
                                      const Real saturation,
                                      const Real brightness,
                                      const Real alpha,
                                      const Real gamma,
                                      const Color& invalid_color)
   : Hsb_Color_Chooser (start_value, end_value,
        string_render ("0:0.8333:%f", gamma),
        string_render ("%f", saturation),
        string_render ("%f", brightness),
        string_render ("%f", alpha),
        invalid_color)
{
}

Hue_Color_Chooser::Hue_Color_Chooser (const Real start_value,
                                      const Real end_value,
                                      const Real start_hue,
                                      const Real end_hue,
                                      const Real saturation,
                                      const Real brightness,
                                      const Real alpha,
                                      const Real gamma,
                                      const Color& invalid_color)
   : Hsb_Color_Chooser (start_value, end_value,
        string_render ("%f:%f:%f", start_hue, end_hue, gamma),
        string_render ("%f", saturation),
        string_render ("%f", brightness),
        string_render ("%f", alpha),
        invalid_color)
{
}

Gray_Color_Chooser::Gray_Color_Chooser (const Real start_value,
                                        const Real end_value,
                                        const string& brightness_str,
                                        const Color& invalid_color)
   : Hsb_Color_Chooser (start_value, end_value, "0", "0",
        brightness_str, "1", invalid_color)
{
}

Gray_Color_Chooser::Gray_Color_Chooser (const Real start_value,
                                        const Real end_value,
                                        const Real gamma,
                                        const Color& invalid_color)
   : Hsb_Color_Chooser (start_value, end_value, "0", "0",
        string_render ("0:1:%f", gamma), "1", invalid_color)
{
}

Mono_Color_Chooser::Mono_Color_Chooser (const Real end_value,
                                        const Real start_value,
                                        const Color_Mode color_mode,
                                        const Real alpha,
                                        const Real hue,
                                        const Real saturation,
                                        const Real brightness,
                                        const Color& invalid_color)
                       : Color_Chooser (invalid_color),
                             end_value (end_value),
                           start_value (start_value),
                            end_extent (1),
                          start_extent (0),
                            color_mode (color_mode),
                                 alpha (alpha),
                                   hue (hue),
                            saturation (saturation),
                            brightness (brightness)
{

   if (color_mode == HUE)
   {
      this->start_extent = 0.833;
      this->end_extent = 0;
   }

   this->value_range = end_value - start_value;
   this->extent_range = end_extent - start_extent;
}

Mono_Color_Chooser::Mono_Color_Chooser (const Real end_value,
                                        const Real start_value,
                                        const Real end_extent,
                                        const Real start_extent,
                                        const Color_Mode color_mode,
                                        const Real alpha,
                                        const Real hue,
                                        const Real saturation,
                                        const Real brightness,
                                        const Color& invalid_color)
                       : Color_Chooser (invalid_color),
                             end_value (end_value),
                           start_value (start_value),
                            end_extent (end_extent),
                          start_extent (start_extent),
                            color_mode (color_mode),
                                 alpha (alpha),
                                   hue (hue),
                            saturation (saturation),
                            brightness (brightness)
{
   this->value_range = end_value - start_value;
   this->extent_range = end_extent - start_extent;
}

Color
Mono_Color_Chooser::get_color (const Real value) const
{

   if (!gsl_finite (value)) { return invalid_color; }

   Color color;
   Real fraction = ((value - start_value) / value_range);
   fraction = (fraction * extent_range) + start_extent;

        if (fraction > 1) { fraction = 1; }
   else if (fraction < 0) { fraction = 0; }

   switch (color_mode)
   {

      case ALPHA:
         color.set_hsb (hue, saturation, brightness, fraction);
         break;

      case HUE:
         color.set_hsb (fraction, saturation, brightness, alpha);
         break;

      case SATURATION:
         color.set_hsb (hue, fraction, brightness, alpha);
         break;

      case BRIGHTNESS:
         color.set_hsb (hue, saturation, fraction, alpha);
         break;

   }

   return color;

}

Bi_Color_Chooser::Bi_Color_Chooser (const Real max_value,
                                    const Real positive_hue,
                                    const Real negative_hue,
                                    const Real alpha,
                                    const Real brightness,
                                    const Real max_saturation,
                                    const Color& invalid_color)
                   : Color_Chooser (invalid_color),
                         max_value (max_value),
                      positive_hue (positive_hue),
                      negative_hue (negative_hue),
                             alpha (alpha),
                        brightness (brightness),
                    max_saturation (max_saturation)
{
}

Color
Bi_Color_Chooser::get_color (const Real value) const
{

   if (!gsl_finite (value)) { return invalid_color; }

   Real v = value;
   v = (v >  max_value ?  max_value : v);
   v = (v < -max_value ? -max_value : v);

   const Real& hue = (v >= 0 ? positive_hue : negative_hue);
   Real saturation = (fabs (v / (2*max_value) * max_saturation));

   return Color::hsb (hue, saturation, brightness, alpha);

}

Multi_Color_Chooser::Multi_Color_Chooser (const Real value_a,
                                          const Real value_b,
                                          const Color& color_a,
                                          const Color& color_b,
                                          const bool staircase)
   : staircase (staircase)
{
   insert (make_pair (value_a, color_a));
   insert (make_pair (value_b, color_b));
}

Multi_Color_Chooser::Multi_Color_Chooser (const Real value_a,
                                          const Real value_b,
                                          const Real value_c,
                                          const Color& color_a,
                                          const Color& color_b,
                                          const Color& color_c,
                                          const bool staircase)
   : staircase (staircase)
{
   insert (make_pair (value_a, color_a));
   insert (make_pair (value_b, color_b));
   insert (make_pair (value_c, color_c));
}

Multi_Color_Chooser::Multi_Color_Chooser (const Real value_a,
                                          const Real value_b,
                                          const Real value_c,
                                          const Real value_d,
                                          const Color& color_a,
                                          const Color& color_b,
                                          const Color& color_c,
                                          const Color& color_d,
                                          const bool staircase)
   : staircase (staircase)
{
   insert (make_pair (value_a, color_a));
   insert (make_pair (value_b, color_b));
   insert (make_pair (value_c, color_c));
   insert (make_pair (value_d, color_d));
}

Multi_Color_Chooser::Multi_Color_Chooser (const Real value_a,
                                          const Real value_b,
                                          const Real value_c,
                                          const Real value_d,
                                          const Real value_e,
                                          const Color& color_a,
                                          const Color& color_b,
                                          const Color& color_c,
                                          const Color& color_d,
                                          const Color& color_e,
                                          const bool staircase)
   : staircase (staircase)
{
   insert (make_pair (value_a, color_a));
   insert (make_pair (value_b, color_b));
   insert (make_pair (value_c, color_c));
   insert (make_pair (value_d, color_d));
   insert (make_pair (value_e, color_e));
}

void
Multi_Color_Chooser::set_staircase (const bool staircase)
{
   this->staircase = staircase;
}

Color
Multi_Color_Chooser::get_color (const Real value) const
{

   Multi_Color_Chooser::const_iterator lower_iterator = lower_bound (value);
   Multi_Color_Chooser::const_iterator upper_iterator = upper_bound (value);

   if (lower_iterator == upper_iterator) // Not a node
   {
      if (lower_iterator == lower_bound (GSL_NEGINF))
      {
         // Before first node
         return begin ()->second;
      }
      else
      if (lower_iterator == end ())
      {
         // After last node
         return rbegin ()->second;
      }
      else
      {
         // Ordinary segment
         lower_iterator--;
      }
   }
   else
   {
      if (upper_iterator == end ())
      {
         lower_iterator--;
         upper_iterator--;
      }
   }

   if (staircase) { return lower_iterator->second; }

   const Real value_a = lower_iterator->first;
   const Real value_b = upper_iterator->first;
   const Color& color_a = lower_iterator->second;
   const Color& color_b = upper_iterator->second;

   const Real r_a = color_a.r;
   const Real g_a = color_a.g;
   const Real b_a = color_a.b;
   const Real a_a = color_a.a;

   const Real r_b = color_b.r;
   const Real g_b = color_b.g;
   const Real b_b = color_b.b;
   const Real a_b = color_b.a;

   const Real u = (value - value_a) / (value_b - value_a);

   const Real r = u * (r_b - r_a) + r_a;
   const Real g = u * (g_b - g_a) + g_a;
   const Real b = u * (b_b - b_a) + b_a;
   const Real a = u * (a_b - a_a) + a_a;

   return Color (r, g, b, a);

}

Multi_Rgb_Color_Chooser::Multi_Rgb_Color_Chooser (const Color& invalid_color)
   : Color_Chooser (invalid_color)
{
}

void
Multi_Rgb_Color_Chooser::add (const Rgb_Color_Chooser& rgb_color_chooser)
{
   push_back (rgb_color_chooser);
}

Color
Multi_Rgb_Color_Chooser::get_color (const Real value) const
{

   for (Multi_Rgb_Color_Chooser::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Rgb_Color_Chooser& rcc = *(iterator);
      if (rcc.within_range (value)) { return rcc.get_color (value); }
   }

   return invalid_color;

}

Level_And_Color::Level_And_Color (const Real level,
                                  const Color& color)
   : level (level),
     color (color)
{
}

Custom_Color_Chooser::Custom_Color_Chooser ()
   : base_color (0, 0, 0, 0)
{
}

Custom_Color_Chooser::Custom_Color_Chooser (const Color& base_color)
   : base_color (base_color)
{
}

void
Custom_Color_Chooser::append (const Real level,
                              const Color& color)
{
   const Level_And_Color level_and_color (level, color);
   level_and_color_vector.push_back (level_and_color);
}

Color
Custom_Color_Chooser::get_color (const Real value) const
{

   Integer index = -1;

   for (Integer i = 0; i < level_and_color_vector.size (); i++)
   {
      if (value < level_and_color_vector[i].level) { break; }
      index = i;
   }

   return ((index < 0) ? base_color : level_and_color_vector[index].color);

}

Rect
Label::cairo_text (const RefPtr<Context>& cr,
                   const Point_2D& point,
                   const bool align,
                   const Real margin_x,
                   const Real margin_y,
                   const bool outline) const
{

   Real dx, dy;
   FontExtents fe;
   TextExtents te;
   cr->get_font_extents (fe);
   cr->get_text_extents (text, te);

   const Real text_width = (align ? te.x_advance : te.width);
   const Real text_height = (align ? -te.y_bearing : te.height);
   //const Real text_height = (align ? (fe.descent + fe.ascent) : te.height);

   switch (justify_h)
   {
      default:
      case 'c': dx = -text_width / 2; break;
      case 'l': dx =  0;              break;
      case 'r': dx = -text_width;     break;
   }

   switch (justify_v)
   {
      default:
      case 'c': dy = text_height / 2; break;
      case 't': dy = text_height;     break;
      case 'b': dy = 0;               break;
   }


   cr->save ();

   if (text_angle != 0)
   {
      cr->translate (point.x, point.y);
      cr->rotate (text_angle);
      cr->translate (-point.x, -point.y);
   }

   cr->move_to (point.x + offset.x + dx, point.y + offset.y + dy); 

   if (outline) { cr->text_path (text); }
   else { cr->show_text (text); }

   cr->restore ();

   Rect rect (point + offset, text_width, text_height,
      justify_h, justify_v, text_angle);

   rect.grow (margin_x, margin_y);
   return rect;

}

string
Label::get_string (const Real number,
                   const Number_Type number_type,
                   const string& format_str,
                   const Real multiplier,
                   const Real offset)
{

   const Real x = number * multiplier + offset;

   switch (number_type)
   {

      default:
      case NUMBER_REAL:
      {
         return string_render (format_str.c_str (), x);
      }

      case NUMBER_TIME:
      {
         return Dtime (x).get_string (format_str);
      }

      case NUMBER_LATITUDE:
      {
         const bool small = fabs (x) < 0.0001;
         const bool positive = x > 0;
         const string ns = (small ? "" : (positive ? "N" : "S"));
         const string& fs = format_str + "\u00b0" + ns;
         return string_render (fs.c_str (), fabs (x));
      }

      case NUMBER_LONGITUDE:
      {
         Lat_Long ll (0, x);
         ll.standardize ();
         const Real longitude = ll.longitude;
         const bool small = fabs (longitude) < 0.0001;
         const bool dateline = fabs (longitude - 180) < 0.0001;
         const bool positive = longitude > 0;
         const string ew = ((small || dateline) ? "" : (positive ? "E" : "W"));
         const string& fs = format_str + "\u00b0" + ew;
         return string_render (fs.c_str (), fabs (longitude));
      }

   }

}

void
Label::set_offset (const char justify_h,
                   const char justify_v,
                   const Real padding)
{

   switch (justify_h)
   {
      default:
      case 'c': offset.x = 0;        break;
      case 'l': offset.x = padding;  break;
      case 'r': offset.x = -padding; break;
   }

   switch (justify_v)
   {
      default:
      case 'c': offset.y = 0;        break;
      case 'b': offset.y = -padding; break;
      case 't': offset.y = padding;  break;
   }

}

Label::Label (const string& text,
              const Point_2D& point_2d,
              const char justify_h,
              const char justify_v)
   : text (text),
     point_2d (point_2d),
     justify_h (justify_h),
     justify_v (justify_v),
     offset (Point_2D (0, 0)),
     text_angle (0)
{
}

Label::Label (const string& text,
              const Point_2D& point_2d,
              const char justify_h,
              const char justify_v,
              const Real padding)
   : text (text),
     point_2d (point_2d),
     justify_h (justify_h),
     justify_v (justify_v),
     text_angle (0)
{
   set_offset (justify_h, justify_v, padding);
}

Label::Label (const Real number,
              const Number_Type number_type,
              const string& format_str,
              const Point_2D& point_2d,
              const char justify_h,
              const char justify_v)
   : text (get_string (number, number_type, format_str)),
     point_2d (point_2d),
     justify_h (justify_h),
     justify_v (justify_v),
     offset (Point_2D (0, 0)),
     text_angle (0)
{
}

Label::Label (const Real number,
              const Number_Type number_type,
              const string& format_str,
              const Point_2D& point_2d,
              const char justify_h,
              const char justify_v,
              const Real padding)
   : text (get_string (number, number_type, format_str)),
     point_2d (point_2d),
     justify_h (justify_h),
     justify_v (justify_v),
     text_angle (0)
{
   set_offset (justify_h, justify_v, padding);
}

Label::Label (const Real number,
              const Number_Type number_type,
              const string& format_str,
              const Real multiplier,
              const Real offset,
              const Point_2D& point_2d,
              const char justify_h,
              const char justify_v)
   : text (get_string (number, number_type, format_str, multiplier, offset)),
     point_2d (point_2d),
     justify_h (justify_h),
     justify_v (justify_v),
     offset (Point_2D (0, 0)),
     text_angle (0)
{
}

Label::Label (const Real number,
              const Number_Type number_type,
              const string& format_str,
              const Real multiplier,
              const Real offset,
              const Point_2D& point_2d,
              const char justify_h,
              const char justify_v,
              const Real padding)
   : text (get_string (number, number_type, format_str, multiplier, offset)),
     point_2d (point_2d),
     justify_h (justify_h),
     justify_v (justify_v),
     text_angle (0)
{
   set_offset (justify_h, justify_v, padding);
}


Label::Label (const Label& label)
   : Rect (label),
     text (label.text),
     point_2d (label.point_2d),
     justify_h (label.justify_h),
     justify_v (label.justify_v),
     offset (label.offset),
     text_angle (label.text_angle)
{
}

void
Label::set_offset (const Point_2D& offset)
{
   this->offset = offset;
}

void
Label::set_text_angle (const Real text_angle)
{
   this->text_angle = text_angle;
}

void
Label::cairo (const RefPtr<Context>& cr,
              const bool align,
              const Real margin_x,
              const Real margin_y,
              const bool outline)
{
   const Point_2D& p = point_2d;
   set (cairo_text (cr, p, align, margin_x, margin_y, outline));
}

void
Label::cairo (const RefPtr<Context>& cr,
              const Transform_2D& transform_2d,
              const bool align,
              const Real margin_x,
              const Real margin_y,
              const bool outline)
{
   const Point_2D& p = transform_2d.transform (point_2d);
   set (cairo_text (cr, p, align, margin_x, margin_y, outline));
}

void
Label::cairo (const RefPtr<Context>& cr,
              const Color& fg_color,
              const Color& bg_color,
              const Point_2D& offset,
              const bool align,
              const Real margin_x,
              const Real margin_y,
              const bool outline)
{

   const Point_2D& p = point_2d;

   cr->save ();

   bg_color.cairo (cr);
   cairo_text (cr, p + offset, align, margin_x, margin_y, outline);
   fg_color.cairo (cr);
   set (cairo_text (cr, p, align, margin_x, margin_y, outline));

   cr->restore ();

}

void
Label::cairo (const RefPtr<Context>& cr,
              const Transform_2D& transform_2d,
              const Color& fg_color,
              const Color& bg_color,
              const Point_2D& offset,
              const bool align,
              const Real margin_x,
              const Real margin_y,
              const bool outline)
{

   const Point_2D& p = transform_2d.transform (point_2d);

   cr->save ();

   bg_color.cairo (cr);
   cairo_text (cr, p + offset, align, margin_x, margin_y, outline);
   fg_color.cairo (cr);
   set (cairo_text (cr, p, align, margin_x, margin_y, outline));

   cr->restore ();

}

Tuple
Simple_Mesh_2D::get_tuple (const Domain_1D domain_1d,
                           const Real interval,
                           const Real multiplier,
                           const Real offset) const
{

   Tuple tuple;
   if (interval <= 0) { return tuple; }

   const bool r = domain_1d.is_reverse ();

   Real s = domain_1d.start * multiplier + offset;
   Real e = domain_1d.end * multiplier + offset;
   const Real& h = interval;

   Integer start_i = r ? Integer (floor (s / h)) : Integer (ceil (s / h));
   Integer end_i = r ? Integer (ceil (e / h)) : Integer (floor (e / h));
   if (start_i > end_i) { std::swap (start_i, end_i); }

   for (Integer i = start_i; i <= end_i; i++)
   {
      const Real x = (i * h - offset) / multiplier;
      tuple.push_back (x);
   }

/*
   Domain_1D d = domain_1d;
   d.swap_if_reverse ();
   const Real start = d.start;
   const Real end = d.end;

   const Real offset = 0;

   Real first = start;
   first += fmod (interval - fmod (start + offset, interval), interval);

   Real x = first;
   Integer i = 0;

   while (x >= start && x <= end)
   {
      tuple.push_back (x);
      i++;
      x = first + i * interval;
   }
*/

}

Tuple
Simple_Mesh_2D::get_tuple_x (const Domain_1D& domain_x) const
{
   if (gsl_isnan (interval_x) ||
       gsl_isnan (multiplier_x) ||
       gsl_isnan (offset_x))
   {
      return tuple_x;
   }
   return get_tuple (domain_x, interval_x, multiplier_x, offset_x);
}

Tuple
Simple_Mesh_2D::get_tuple_y (const Domain_1D& domain_y) const
{
   if (gsl_isnan (interval_y) ||
       gsl_isnan (multiplier_y) ||
       gsl_isnan (offset_y))
   {
      return tuple_y;
   }
   return get_tuple (domain_y, interval_y, multiplier_y, offset_y);
}

Simple_Mesh_2D::Simple_Mesh_2D (const Color& color,
                                const Tuple& tuple_x,
                                const Tuple& tuple_y)
   : color (color),
     tuple_x (tuple_x),
     tuple_y (tuple_y),
     interval_x (GSL_NAN),
     multiplier_x (GSL_NAN),
     offset_x (GSL_NAN),
     interval_y (GSL_NAN),
     multiplier_y (GSL_NAN),
     offset_y (GSL_NAN),
     line_width (1)
{
}

Simple_Mesh_2D::Simple_Mesh_2D (const Color& color,
                                const Tuple& tuple_x,
                                const Real interval_y,
                                const Real multiplier_y,
                                const Real offset_y)
   : color (color),
     tuple_x (tuple_x),
     interval_x (GSL_NAN),
     multiplier_x (GSL_NAN),
     offset_x (GSL_NAN),
     interval_y (interval_y),
     multiplier_y (multiplier_y),
     offset_y (offset_y),
     line_width (1)
{
}

Simple_Mesh_2D::Simple_Mesh_2D (const Color& color,
                                const Real interval_x,
                                const Tuple& tuple_y,
                                const Real multiplier_x,
                                const Real offset_x)
   : color (color),
     tuple_y (tuple_y),
     interval_x (interval_x),
     multiplier_x (multiplier_x),
     offset_x (offset_x),
     interval_y (GSL_NAN),
     multiplier_y (GSL_NAN),
     offset_y (GSL_NAN),
     line_width (1)
{
}

Simple_Mesh_2D::Simple_Mesh_2D (const Color& color,
                                const Real interval_x,
                                const Real interval_y,
                                const Real multiplier_x,
                                const Real offset_x,
                                const Real multiplier_y,
                                const Real offset_y)
   : color (color),
     interval_x (interval_x),
     interval_y (interval_y),
     multiplier_x (multiplier_x),
     offset_x (offset_x),
     multiplier_y (multiplier_y),
     offset_y (offset_y),
     line_width (1)
{
}

void
Simple_Mesh_2D::set_line_width (const Real line_width)
{
   this->line_width = line_width;
}

void
Simple_Mesh_2D::render (const RefPtr<Context>& cr,
                        const Transform_2D& transform,
                        const Size_2D& size_2d,
                        const Domain_2D& domain_2d) const
{

   Point_2D point;
   Domain_1D domain_x = domain_2d.domain_x;
   Domain_1D domain_y = domain_2d.domain_y;
   domain_x.swap_if_reverse ();
   domain_y.swap_if_reverse ();

   const Real& start_x = domain_x.start;
   const Real& end_x   = domain_x.end;
   const Real& start_y = domain_y.start;
   const Real& end_y   = domain_y.end;
   const Real delta_x = domain_x.get_span () / (size_2d.i - 1);
   const Real delta_y = domain_y.get_span () / (size_2d.j - 1);

   const Tuple& tuple_x = get_tuple_x (domain_x);
   const Tuple& tuple_y = get_tuple_y (domain_y);

   color.cairo (cr);
   cr->set_line_width (line_width);

   for (Tuple::const_iterator iterator = tuple_x.begin ();
        iterator != tuple_x.end (); iterator++)
   {

      point.x = *(iterator);
      Simple_Polyline simple_polyline;

      for (Integer j = 0; j < size_2d.j; j++)
      {
         point.y = start_y + j * delta_y;
         simple_polyline.add (point);
      }

      transform.cairo (cr, simple_polyline);
      cr->stroke ();

   }

   for (Tuple::const_iterator iterator = tuple_y.begin ();
        iterator != tuple_y.end (); iterator++)
   {

      point.y = *(iterator);
      Simple_Polyline simple_polyline;

      for (Integer i = 0; i < size_2d.i; i++)
      {
         point.x = start_x + i * delta_x;
         simple_polyline.add (point);
      }

      transform.cairo (cr, simple_polyline);
      cr->stroke ();

   }

}

void
Simple_Mesh_2D::render_label_x (const RefPtr<Context>& cr,
                                const Transform_2D& transform,
                                const Domain_1D& domain_1d,
                                const Real position,
                                const string& format_str,
                                const Number_Type number_type,
                                const char justify_h,
                                const char justify_v,
                                const Real padding) const
{

   Domain_1D d = domain_1d;
   d.swap_if_reverse ();

   const Real offset = 0;

   Point_2D point (0, position);
   const Tuple& tuple = get_tuple_x (d);

   const bool non_regular = (gsl_isnan (multiplier_x) || gsl_isnan (offset_x));
   const Real m_x = (non_regular ? 1 : multiplier_x);
   const Real o_x = (non_regular ? 1 : offset_x);

   for (Tuple::const_iterator iterator = tuple.begin ();
        iterator != tuple.end (); iterator++)
   {

      point.x = *(iterator);

      Label label (point.x, number_type, format_str, m_x,
         o_x, point, justify_h, justify_v, padding);
      label.cairo (cr, transform);

   }

}

void
Simple_Mesh_2D::render_label_y (const RefPtr<Context>& cr,
                                const Transform_2D& transform,
                                const Domain_1D& domain_1d,
                                const Real position,
                                const string& format_str,
                                const Number_Type number_type,
                                const char justify_h,
                                const char justify_v,
                                const Real padding) const
{

   Domain_1D d = domain_1d;
   d.swap_if_reverse ();

   const Real offset = 0;

   Point_2D point (position, 0);
   const Tuple& tuple = get_tuple_y (d);

   const bool non_regular = (gsl_isnan (multiplier_y) || gsl_isnan (offset_y));
   const Real m_y = (non_regular ? 1 : multiplier_y);
   const Real o_y = (non_regular ? 1 : offset_y);

   for (Tuple::const_iterator iterator = tuple.begin ();
        iterator != tuple.end (); iterator++)
   {

      point.y = *(iterator);

      Label label (point.y, number_type, format_str, m_y,
         o_y, point, justify_h, justify_v, padding);
      label.cairo (cr, transform);

   }

}

void
Simple_Mesh_2D::render_label_lat_long (const RefPtr<Context>& cr,
                                       const Transform_2D& transform,
                                       const Domain_2D& domain_2d,
                                       const Lat_Long& lat_long,
                                       const string& format_str,
                                       const char justify_h,
                                       const char justify_v,
                                       const Real padding) const
{

   render_label_x (cr, transform, domain_2d.domain_x,
      lat_long.longitude, format_str, NUMBER_LATITUDE, justify_h,
      justify_v, padding);

   render_label_y (cr, transform, domain_2d.domain_y,
      lat_long.latitude, format_str, NUMBER_LONGITUDE, justify_h,
      justify_v, padding);

}

Mesh_2D::Mesh_2D ()
{
}

Mesh_2D::Mesh_2D (const Domain_2D& domain_2d)
   : domain_2d (domain_2d)
{
}

Mesh_2D::Mesh_2D (const Domain_2D& domain_2d,
                  const Simple_Mesh_2D& simple_mesh_2d)
   : domain_2d (domain_2d)
{
   add (simple_mesh_2d);
}

Mesh_2D::Mesh_2D (const Domain_2D& domain_2d,
                  const Simple_Mesh_2D& simple_mesh_2d_a,
                  const Simple_Mesh_2D& simple_mesh_2d_b)
   : domain_2d (domain_2d)
{
   add (simple_mesh_2d_a);
   add (simple_mesh_2d_b);
}

Mesh_2D::Mesh_2D (const Domain_2D& domain_2d,
                  const Simple_Mesh_2D& simple_mesh_2d_a,
                  const Simple_Mesh_2D& simple_mesh_2d_b,
                  const Simple_Mesh_2D& simple_mesh_2d_c)
   : domain_2d (domain_2d)
{
   add (simple_mesh_2d_a);
   add (simple_mesh_2d_b);
   add (simple_mesh_2d_c);
}

Mesh_2D::Mesh_2D (const Domain_2D& domain_2d,
                  const Simple_Mesh_2D& simple_mesh_2d_a,
                  const Simple_Mesh_2D& simple_mesh_2d_b,
                  const Simple_Mesh_2D& simple_mesh_2d_c,
                  const Simple_Mesh_2D& simple_mesh_2d_d)
   : domain_2d (domain_2d)
{
   add (simple_mesh_2d_a);
   add (simple_mesh_2d_b);
   add (simple_mesh_2d_c);
   add (simple_mesh_2d_d);
}

Mesh_2D::Mesh_2D (const Domain_2D& domain_2d,
                  const vector<Simple_Mesh_2D>& simple_mesh_2d_vector)
   : domain_2d (domain_2d)
{

   typedef vector<Simple_Mesh_2D>::const_iterator Iterator;

   for (Iterator iterator = simple_mesh_2d_vector.begin ();
        iterator != simple_mesh_2d_vector.end (); iterator++)
   {
      const Simple_Mesh_2D& simple_mesh_2d = *(iterator);
      add (simple_mesh_2d);
   }

}

void
Mesh_2D::set_domain_2d (const Domain_2D& domain_2d)
{
   this->domain_2d = domain_2d;
}

void
Mesh_2D::set_domain_2d (const Domain_1D& domain_x,
                        const Domain_1D& domain_y)
{
   this->domain_2d.domain_x = domain_x;
   this->domain_2d.domain_y = domain_y;
}

void
Mesh_2D::add (const Simple_Mesh_2D& simple_mesh_2d)
{
   push_back (simple_mesh_2d);
}

const Domain_2D&
Mesh_2D::get_domain_2d () const
{
   return domain_2d;
}

void
Mesh_2D::render (const RefPtr<Context>& cr,
                 const Transform_2D& transform,
                 const Size_2D& size_2d) const
{

   cr->save ();

   for (Mesh_2D::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Simple_Mesh_2D& simple_mesh_2d = *(iterator);
      simple_mesh_2d.render (cr, transform, size_2d, domain_2d);
   }

   cr->restore ();

}

void
Mesh_2D::render_label_x (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Integer index,
                         const Real position_y,
                         const string& format_str,
                         const Number_Type number_type,
                         const char justify_h,
                         const char justify_v,
                         const Real padding) const
{
   if (index >= size ()) { return; }
   const Simple_Mesh_2D& simple_mesh_2d = at (index);
   simple_mesh_2d.render_label_x (cr, transform, domain_2d.domain_x,
      position_y, format_str, number_type, justify_h, justify_v, padding);
}

void
Mesh_2D::render_label_y (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Integer index,
                         const Real position_x,
                         const string& format_str,
                         const Number_Type number_type,
                         const char justify_h,
                         const char justify_v,
                         const Real padding) const
{
   if (index >= size ()) { return; }
   const Simple_Mesh_2D& simple_mesh_2d = at (index);
   simple_mesh_2d.render_label_y (cr, transform, domain_2d.domain_y,
      position_x, format_str, number_type, justify_h, justify_v, padding);
}

void
Mesh_2D::render_label_lat_long (const RefPtr<Context>& cr,
                                const Transform_2D& transform,
                                const Integer index,
                                const Lat_Long& lat_long,
                                const string& format_str,
                                const char justify_h,
                                const char justify_v,
                                const Real padding) const
{
   render_label_x (cr, transform, index, lat_long.longitude,
      format_str, NUMBER_LATITUDE, justify_h, justify_v, padding);
   render_label_y (cr, transform, index, lat_long.latitude,
      format_str, NUMBER_LONGITUDE, justify_h, justify_v, padding);
}

Raster::Raster (const Size_2D& size_2d)
   : Box_2D (size_2d)
{
   buffer = new uint32_t[size_2d.i * size_2d.j];
}

Raster::Raster (const Box_2D& box_2d)
   : Box_2D (box_2d)
{
   buffer = new uint32_t[box_2d.size_2d.i * box_2d.size_2d.j];
}

Raster::~Raster ()
{
   delete[] buffer;
}

void
Raster::blend_pixel (const Integer i,
                     const Integer j,
                     const Color& color)
{

   const Integer offset = j * size_2d.i + i;
   uint32_t& pixel = *(buffer + offset);

   const Real r_orig = Real ((pixel & 0x00ff0000) >> 16) / 255;
   const Real g_orig = Real ((pixel & 0x0000ff00) >> 8) / 255;
   const Real b_orig = Real ((pixel & 0x000000ff) / 255);

   const Real alpha = color.a;
   const Real alpha_com = 1 - alpha;

   uint32_t a = 255;
   uint32_t r = uint32_t (rint ((r_orig * alpha_com + color.r * alpha) * 255));
   uint32_t g = uint32_t (rint ((g_orig * alpha_com + color.g * alpha) * 255));
   uint32_t b = uint32_t (rint ((b_orig * alpha_com + color.b * alpha) * 255));

   pixel = (a << 24) + (r << 16) + (g << 8) + b;

}

void
Raster::set_pixel (const Integer i,
                   const Integer j,
                   const Color& color)
{

   const Integer offset = j * size_2d.i + i;
   uint32_t& pixel = *(buffer + offset);

   uint32_t a = uint32_t (rint (color.a * 255));
   uint32_t r = uint32_t (rint (color.a * color.r * 255));
   uint32_t g = uint32_t (rint (color.a * color.g * 255));
   uint32_t b = uint32_t (rint (color.a * color.b * 255));

   pixel = (a << 24) + (r << 16) + (g << 8) + b;

}

Color
Raster::get_pixel (const Integer i,
                   const Integer j) const
{

   const Integer offset = j * size_2d.i + i;
   uint32_t& pixel = *(buffer + offset);

   const Real a = Real ((pixel & 0xff000000) >> 24) / 255;
   const Real r = Real ((pixel & 0x00ff0000) >> 16) / 255;
   const Real g = Real ((pixel & 0x0000ff00) >> 8) / 255;
   const Real b = Real ((pixel & 0x000000ff) / 255);

   return Color (r / a, g / a, b / a, a);

}

RefPtr<ImageSurface>
Raster::get_image_surface () const
{
   return ImageSurface::create ((unsigned char*)buffer,
      FORMAT_ARGB32, size_2d.i, size_2d.j, size_2d.i * 4);
}

void
Raster::blit (const RefPtr<Context>& cr,
              const Real alpha) const
{
   RefPtr<ImageSurface> surface = get_image_surface ();
   cr->set_source (surface, index_2d.i, index_2d.j);
   cr->paint_with_alpha (alpha);
}

void
Title::cairo (const RefPtr<Context>& cr,
              const string& string_l,
              const string& string_c,
              const string& string_r) const
{

   const Real width = Real (i);
   const Real title_height = get_height ();
   const Real font_size = title_height / 2;
   const Real y = title_height * 0.875;
   const Point_2D origin (0, 0);
   const Point_2D shadow_offset (2, -2);
   
   cr->save ();
   cr->set_font_size (font_size);
   
   Label label_l (string_l, Point_2D (10, y), 'l', 'b');
   Label label_c (string_c, Point_2D (width / 2, y), 'c', 'b');
   Label label_r (string_r, Point_2D (width - 10, y), 'r', 'b');
   
   shadow_color.cairo (cr);
   label_l.set_offset (shadow_offset);
   label_l.cairo (cr);
   label_c.set_offset (shadow_offset);
   label_c.cairo (cr);
   label_r.set_offset (shadow_offset);
   label_r.cairo (cr); 
   
   fg_color.cairo (cr);
   label_l.set_offset (origin);
   label_l.cairo (cr); 
   label_c.set_offset (origin);
   label_c.cairo (cr); 
   label_r.set_offset (origin);
   label_r.cairo (cr);

   cr->restore ();

}

void
Title::cairo (const RefPtr<Context>& cr,
              const string& string_ul,
              const string& string_ll,
              const string& string_c,
              const string& string_ur,
              const string& string_lr) const
{

   const Real width = Real (i);
   const Real title_height = get_height ();
   const Real large_font_size = title_height * 0.500;
   const Real small_font_size = title_height * 0.400;
   const Real upper_y = title_height * 0.125;
   const Real lower_y = title_height * 0.875;
   const Point_2D origin (0, 0);
   const Point_2D shadow_offset (2, -2);

   cr->save ();

   Label label_ul (string_ul, Point_2D (10, upper_y), 'l', 't');
   Label label_ll (string_ll, Point_2D (10, lower_y), 'l', 'b');
   Label label_c (string_c, Point_2D (width / 2, lower_y), 'c', 'b');
   Label label_ur (string_ur, Point_2D (width - 10, upper_y), 'r', 't');
   Label label_lr (string_lr, Point_2D (width - 10, lower_y), 'r', 'b');

   shadow_color.cairo (cr);
   cr->set_font_size (large_font_size);
   label_c.set_offset (shadow_offset);
   label_c.cairo (cr);
   cr->set_font_size (small_font_size);
   label_ul.set_offset (shadow_offset);
   label_ul.cairo (cr);
   label_ll.set_offset (shadow_offset);
   label_ll.cairo (cr);
   label_ur.set_offset (shadow_offset);
   label_ur.cairo (cr);
   label_lr.set_offset (shadow_offset);
   label_lr.cairo (cr);

   fg_color.cairo (cr);
   cr->set_font_size (large_font_size);
   label_c.set_offset (origin);
   label_c.cairo (cr);
   cr->set_font_size (small_font_size);
   label_ul.set_offset (origin);
   label_ul.cairo (cr);
   label_ll.set_offset (origin);
   label_ll.cairo (cr);
   label_ur.set_offset (origin);
   label_ur.cairo (cr);
   label_lr.set_offset (origin);
   label_lr.cairo (cr);

   cr->restore ();

}

Title::Title (const Size_2D& size_2d,
              const Color& bg_color,
              const Color& fg_color,
              const Color& shadow_color)
   : Size_2D (size_2d),
     bg_color (bg_color),
     fg_color (fg_color),
     shadow_color (shadow_color)
{
}

void
Title::set_size_2d (const Size_2D& size_2d)
{
   this->i = size_2d.i;
   this->j = size_2d.j;
}

Real
Title::get_height () const
{
   const Real h = j * 0.042;
   return std::max (std::min (h, 45.0), 35.0);
}

void
Title::set (const string& string_l,
            const string& string_c,
            const string& string_r)
{
   clear ();
   push_back (string_l);
   push_back (string_c);
   push_back (string_r);
}

void
Title::set (const string& string_ul,
            const string& string_ll,
            const string& string_c,
            const string& string_ur,
            const string& string_lr)
{
   clear ();
   push_back (string_ul);
   push_back (string_ll);
   push_back (string_c);
   push_back (string_ur);
   push_back (string_lr);
}

void
Title::cairo (const RefPtr<Context>& cr)
{

   const Real width = Real (i);
   const Real height = Real (j);
   const Real title_height = get_height ();

   cr->save ();
   bg_color.cairo (cr);
   Rect (Point_2D (0, 0), width, title_height).cairo (cr);
   cr->fill ();
   cr->set_line_width (6);
   Rect (Point_2D (0, 0), width, height).cairo (cr);
   cr->stroke ();

   const Tokens& tokens = *(this);

   if (size () == 3)
   {
      cairo (cr, tokens[0], tokens[1], tokens[2]);
   }
   else
   if (size () == 5)
   {
      cairo (cr, tokens[0], tokens[1], tokens[2], tokens[3], tokens[4]);
   }

   cr->restore ();

}

namespace denise
{

   ostream&
   operator << (ostream &out_file,
                const Color& color)
   {
      out_file << "Color (" << color.r << ", " << color.g
               << ", " << color.b << ", " << color.a << ")";
      return out_file;
   }

}

