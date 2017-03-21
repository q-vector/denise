//
// met.cc
// 
// Copyright (C) 2005-2015 Simon E. Ching
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

#include "met.h"
#include "thermo.h"

using namespace denise;

namespace denise
{

   Tc_Category
   get_tc_category (const Real max_wind)
   {

      Tc_Category tc_category = LOW;

      if (max_wind >= 32.7) { tc_category = TC_T; }

      else
      if (max_wind >= 24.5) { tc_category = TC_STS; }

      else
      if (max_wind >= 17.2) { tc_category = TC_TS; }

      else
      if (max_wind >= 10.8) { tc_category = TC_TD; }

      return tc_category;

   }

   Tc_Category
   get_tc_category (const Dstring& tc_category_string)
   {

      Tc_Category tc_category = UNKNOWN_TC_CATEGORY;

      if ((tc_category_string == "AL") ||
          (tc_category_string == "ALL_TC_CATEGORIES"))
      {
         tc_category = ALL_TC_CATEGORIES;
      }

      else
      if ((tc_category_string == "LOW"))
      {
         tc_category = LOW;
      }

      else
      if (((tc_category_string == "TD")) ||
          ((tc_category_string == "TROPICAL_DEPRESSION")))
      {
         tc_category = TC_TD;
      }

      else
      if (((tc_category_string == "TS")) ||
          ((tc_category_string == "TROPICAL_STORM")))
      {
         tc_category = TC_TS;
      }

      else
      if (((tc_category_string == "STS")) ||
          ((tc_category_string == "SEVERE_TROPICAL_STORM")))
      {
         tc_category = TC_STS;
      }

      else
      if (((tc_category_string == "T")) ||
          ((tc_category_string == "TY")) ||
          ((tc_category_string == "TYPHOON")))
      {
         tc_category = TC_T;
      }

      return tc_category;

   }

   Real
   get_maximal_intensity (const Tc_Category tc_category)
   {

      Real intensity = GSL_NAN;

      switch (tc_category)
      {
         case TC_TD:  intensity = 17.0; break;
         case TC_TS:  intensity = 24.3; break;
         case TC_STS: intensity = 32.5; break;
         case TC_T:   intensity = 99.9; break;
      }

      return intensity;

   }

   Real
   get_minimal_intensity (const Tc_Category tc_category)
   {

      Real intensity = 0;

      switch (tc_category)
      {
         case TC_TD:  intensity =  0.0; break;
         case TC_TS:  intensity = 17.4; break;
         case TC_STS: intensity = 24.7; break;
         case TC_T:   intensity = 32.9; break;
      }

      return intensity;

   }

   Dstring
   get_tc_category_string (const Tc_Category tc_category)
   {

      Dstring s ("XXX");

      switch (tc_category)
      {
         case LOW:    s.assign ("LOW"); break;
         case TC_TD:  s.assign ("TD");  break;
         case TC_TS:  s.assign ("TS");  break;
         case TC_STS: s.assign ("STS"); break;
         case TC_T:   s.assign ("T");   break;
         case ALL_TC_CATEGORIES: s.append ("ALL_TC_CATEGORIES"); break;
      }

      return s;

   }

}

void
Tc_Symbol::init_north (const Point_2D& point,
                       const Point_2D& left,
                       const Point_2D& right,
                       const Point_2D& up,
                       const Point_2D& down,
                       const Real r,
                       const Real r_o,
                       const Real r_i,
                       const Real s_o,
                       const Real s_i,
                       const Real w,
                       const Real w_2,
                       const Real theta,
                       const Real phi,
                       const Real lambda)
{

   add_arc (true, right, r_o, M_PI, M_PI - lambda);
   add_arc (false, up, w/2, M_PI - lambda, -lambda);
   add_arc (false, right, r_i, M_PI - lambda, M_PI - phi);
   add_arc (false, point, s_o, M_PI - theta, 0);
   add_arc (false, left, r_o, 0, -lambda);
   add_arc (false, down, w/2, -lambda, -lambda - M_PI);
   add_arc (false, left, r_i, -lambda, -phi);
   add_arc (false, point, s_o, -theta, -M_PI);

   add_circle (true, point, s_i);

}

void
Tc_Symbol::init_south (const Point_2D& point,
                       const Point_2D& left,
                       const Point_2D& right,
                       const Point_2D& up,
                       const Point_2D& down,
                       const Real r,
                       const Real r_o,
                       const Real r_i,
                       const Real s_o,
                       const Real s_i,
                       const Real w,
                       const Real w_2,
                       const Real theta,
                       const Real phi,
                       const Real lambda)
{

   add_arc (true, left, r_o, 0, lambda);
   add_arc (false, up, w/2, lambda, M_PI + lambda);
   add_arc (false, left, r_i, lambda, phi);
   add_arc (false, point, s_o, theta, M_PI);
   add_arc (false, right, r_o, M_PI, M_PI + lambda);
   add_arc (false, down, w/2, lambda - M_PI, lambda);
   add_arc (false, right, r_i, M_PI + lambda, M_PI + phi);
   add_arc (false, point, s_o, M_PI + theta, M_2_TIMES_PI);

   add_circle (true, point, s_i);

}

Tc_Symbol::Tc_Symbol (const bool southern_hemisphere,
                      const Real size,
                      const Real width)
            : Symbol (size)
{

   Real s_o = size;
   Real w = (gsl_finite (width) ? width : size / 3);
   Real w_2 = w / 2;
   Real s_i = s_o - w;
   Real r_o = s_o * 2.75;
   Real r_i = r_o - w;
   Real c = r_o - s_o;
   Real r = r_o - w_2;

   Real x = -(s_o*s_o - r_i*r_i + c*c) / (2*c);
   Real theta = acos (x / s_o);
   Real phi = acos ((r_o - s_o + x) / r_i);
   Real lambda = acos ((r_o - s_o) / r);

   Point_2D left (-c, 0);
   Point_2D right (c, 0);
   Point_2D up (0, r * sin (lambda));
   Point_2D down (0, -r * sin (lambda));

   Point_2D point (0, 0);

   if (southern_hemisphere)
   {
      init_south (point, left, right, up, down, r, r_o,
         r_i, s_o, s_i, w, w_2, theta, phi, lambda);
   }
   else
   {
      init_north (point, left, right, up, down, r, r_o,
         r_i, s_o, s_i, w, w_2, theta, phi, lambda);
   }

}

Octa::Octa (const Octa::Number octa_number,
            const Real size,
            const Real width)
  : Symbol (size)
{

   Point_2D point (0, 0);

   Real s_o = size;
   Real w = (gsl_finite (width) ? width : size / 3);
   Real s_i = s_o - w;

   switch (octa_number)
   {

      case OCTA_ZERO:
      {
         add_circle (true, point, s_o);
         add_circle (true, point, s_i);
         break;
      }

      case OCTA_ONE:
      {
         Real theta = atan ((w/2) / s_i);
         add_circle (true, point, s_o);
         add_arc (true, point, s_i, -M_PI/2 + theta, M_PI/2 - theta);
         add_arc (true, point, s_i, M_PI/2 + theta, M_PI/2*3 - theta);
         break;
      }

      case OCTA_TWO:
      {
         add_circle (true, point, s_o);
         add_sector (true, point, s_i, M_PI/2, M_2_TIMES_PI);
         break;
      }

      case OCTA_THREE:
      {
         Real theta = atan ((w/2) / s_i);
         add_circle (true, point, s_o);
         add_arc (true, point, s_i, M_PI/2 + theta, M_PI/2*3 - theta);
         add_arc (true, point, s_i, M_PI/2*3 + theta, M_2_TIMES_PI - theta);
         add (Point_2D (point.x + w/2, point.y - w/2));
         break;
      }

      case OCTA_FOUR:
      {
         add_circle (true, point, s_o);
         add_arc (true, point, s_i, M_PI/2, M_PI/2*3);
         break;
      }

      case OCTA_FIVE:
      {
         Real theta = atan ((w/2) / s_i);
         add_circle (true, point, s_o);
         add_arc (true, point, s_i, M_PI/2 + theta, M_PI - theta);
         add (Point_2D (point.x - w/2, point.y + w/2));
         add_arc (true, point, s_i, M_PI + theta, M_PI/2*3 - theta);
         add (Point_2D (point.x - w/2, point.y - w/2));
         break;
      }

      case OCTA_SIX:
      {
         add_circle (true, point, s_o);
         add_sector (true, point, s_i, M_PI/2, M_PI);
         break;
      }

      case OCTA_SEVEN:
      {
         add_circle (true, point, s_o);
         add_arc (true, point, s_i, M_PI*0.425, M_PI*0.575);
         add_arc (false, point, s_i, M_PI*1.425, M_PI*1.575);
         break;
      }

      case OCTA_EIGHT:
      {
         add_circle (true, point, s_o);
         break;
      }

      case OCTA_OBSCURED:
      {
         Real theta = atan ((w/2) / s_i);
         add_circle (true, point, s_o);
         add_arc (true, point, s_i, M_PI/4 + theta, M_PI/4*3 - theta);
         add (Point_2D (point.x, point.y + w/M_SQRT2));
         add_arc (true, point, s_i, M_PI/4*3 + theta, M_PI/4*5 - theta);
         add (Point_2D (point.x - w/M_SQRT2, point.y));
         add_arc (true, point, s_i, M_PI/4*5 + theta, M_PI/4*7 - theta);
         add (Point_2D (point.x, point.y - w/M_SQRT2));
         add_arc (true, point, s_i, M_PI/4*7 + theta, M_PI/4*9 - theta);
         add (Point_2D (point.x + w/M_SQRT2, point.y));
         break;
      }

      case OCTA_GARBLED:
      {
         Real sin_theta = sin (M_PI/15);
         Real theta_minus = asin ((s_i * sin_theta - w/2) / s_i);
         Real theta_plus = asin ((s_i * sin_theta + w/2) / s_i);
         add_circle (true, point, s_o);
         add_arc (true, point, s_i, theta_plus, M_PI - theta_plus);
         add_arc (true, point, s_i, -M_PI + theta_plus, -theta_plus);
         add_arc (true, point, s_i, -theta_minus, theta_minus);
         add_arc (false, point, s_i, M_PI - theta_minus, M_PI + theta_minus);
         break;
      }

      default:
         throw Exception ("Bad Octa Type");
         break;

   }

}

Wind::Wind (const Real u,
            const Real v)
   : Motion (u, v)
{
}

Wind::Wind (const Wind& wind)
   : Motion (wind)
{
}

Wind
Wind::direction_speed (const Real direction,
                       const Real speed)
{
   Real u = -speed * sin (direction * DEGREE_TO_RADIAN);
   Real v = -speed * cos (direction * DEGREE_TO_RADIAN);
   return Wind (u, v);
}

void
Wind::set_from_direction_speed (const Real direction,
                                const Real speed)
{
   u = -speed * sin (direction * DEGREE_TO_RADIAN);
   v = -speed * cos (direction * DEGREE_TO_RADIAN);
}

void
Wind::normalize (const Real magnitude)
{
   const Real multiplier = magnitude / get_speed ();
   this->u *= multiplier;
   this->v *= multiplier;
}

void
Wind::veer (const Real degree)
{
   Wind temp (u, v);
   Real theta = -degree * DEGREE_TO_RADIAN;
   u = temp.u * cos (theta) - temp.v * sin (theta);
   v = temp.u * sin (theta) + temp.v * cos (theta);
}

void
Wind::back (const Real degree)
{
   veer (-degree);
}

Wind
Wind::get_normalized_wind (const Real magnitude) const
{
   Wind normalized (u, v);
   normalized.normalize (magnitude);
   return normalized;
}

Wind
Wind::get_veered_wind (const Real degree) const
{
   Wind veered (u, v);
   veered.veer (degree);
   return veered;
}

Wind
Wind::get_backed_wind (const Real degree) const
{
   return get_veered_wind (-degree);
}

Real
Wind::get_direction () const
{
   if (u == 0 && v == 0) { return 0; }
   return fmod (atan2 (-u, -v) * RADIAN_TO_DEGREE + 360, 360);
}

Dstring
Wind::get_direction_string (const Integer n,
                            const Dstring& format) const
{

   const Real d = get_direction ();
   const Real delta_d = 360.0 / n;
   const Integer dd = Integer (floor ((d + delta_d / 2) / delta_d));

   const bool no_format = (format == "");

   if (n == 4 && no_format)
   {
      switch (dd)
      {
         case 0: return "N";
         case 1: return "E";
         case 2: return "S";
         case 3: return "W";
      }
   }

   if (n == 8 && no_format)
   {
      switch (dd)
      {
         case 0: return "N";
         case 1: return "NE";
         case 2: return "E";
         case 3: return "SE";
         case 4: return "S";
         case 5: return "SW";
         case 6: return "W";
         case 7: return "NW";
      }
   }

   if (n == 16 && no_format)
   {
      switch (dd)
      {
         case 0: return "N";
         case 1: return "NNE";
         case 2: return "NE";
         case 3: return "ENE";
         case 4: return "E";
         case 5: return "ESE";
         case 6: return "SE";
         case 7: return "SSE";
         case 8: return "S";
         case 9: return "SSW";
         case 10: return "SW";
         case 11: return "WSW";
         case 12: return "W";
         case 13: return "WNW";
         case 14: return "NW";
         case 15: return "NNW";
      }
   }

   if (n == 32 && no_format)
   {
      switch (dd)
      {
         case 0: return "N";
         case 1: return "NbE";
         case 2: return "NNE";
         case 3: return "NEbE";
         case 4: return "NE";
         case 5: return "NEbN";
         case 6: return "ENE";
         case 7: return "EbN";
         case 8: return "E";
         case 9: return "EbS";
         case 10: return "ESE";
         case 11: return "SEbE";
         case 12: return "SE";
         case 13: return "SEbS";
         case 14: return "SSE";
         case 15: return "SbE";
         case 16: return "S";
         case 17: return "SbW";
         case 18: return "SSW";
         case 19: return "SWbS";
         case 20: return "SW";
         case 21: return "SWbW";
         case 22: return "WSW";
         case 23: return "WbS";
         case 24: return "W";
         case 25: return "WbN";
         case 26: return "WNW";
         case 27: return "NWbW";
         case 28: return "NW";
         case 29: return "NWbN";
         case 30: return "NNW";
         case 31: return "NbW";
      }
   }

   const Real cardinal_direction = (dd == 0 ? n : dd) * delta_d;
   const Dstring& fmt = (no_format ? "%03.0f\u00b0" : format);
   return Dstring::render (fmt, cardinal_direction);

}

bool
Wind::is_naw () const
{
   return gsl_isnan (u) || gsl_isnan (v);
}

Dstring
Wind::get_string (const Real speed_multiplier,
                  const Dstring& format) const
{
   const Real direction = get_direction ();
   const Real speed = get_speed () * speed_multiplier;
   return Dstring::render (format, direction, speed);
}

Dstring
Wind::get_string (const Dstring& format,
                  const Real speed_multiplier) const
{
   const Real direction = get_direction ();
   const Real speed = get_speed () * speed_multiplier;
   return Dstring::render (format, direction, speed);
}

Wind::Category
Wind::get_category (const Real speed)
{
        if (speed >= 32.7) { return HURRICANE; }
   else if (speed >= 24.5) { return STORM; }
   else if (speed >= 17.2) { return GALE; }
   else if (speed >= 10.8) { return STRONG; }
   else if (speed >= 8.0)  { return FRESH; }
   else if (speed >= 3.4)  { return MODERATE; }
   else if (speed >= 0.2)  { return LIGHT; }
   else                    { return CALM; }
}

Dstring
Wind::get_category_string (const Wind::Category category)
{
   switch (category)
   {
      case CALM:      return "Calm";      break;
      case LIGHT:     return "Light";     break;
      case MODERATE:  return "Moderate";  break;
      case FRESH:     return "Fresh";     break;
      case STRONG:    return "Strong";    break;
      case GALE:      return "Gale";      break;
      case STORM:     return "Storm";     break;
      case HURRICANE: return "Hurricane"; break;
   }
}

void
Wind_Barb::analysis (Integer& pennants,
                     Integer& barbs,
                     Real& residue,
                     Real speed)
{

   pennants = 0;
   barbs = 0;

   while (speed >= 25) { pennants++; speed -= 25; }
   while (speed >= 5)  { barbs++; speed -= 5; }

   residue = speed;
//   if (residue < 2 && pennants == 0 && barbs == 0) { residue = 2; }

//   if (barbs == 0 && residue > 0) { barbs = 1; }

}

void
Wind_Barb::add_pennant (const Real theta,
                        const Point_2D& point,
                        const Real length)
{

   Real phi = (north ? theta + M_PI_2 : theta - M_PI_2);
   Real lambda = length * M_SQRT1_3;

   add (point);
   add (Point_2D (point.x + length*cos (phi), point.y + length * sin (phi)));
   add (Point_2D (point.x - lambda*cos (theta), point.y - lambda * sin (theta)));

}

void
Wind_Barb::add_barb (const Real theta,
                     const Point_2D& point,
                     const Real length,
                     const Real width,
                     const bool invisible_barb)
{

   Real chi = M_PI_2 - M_PI/6;
   Real psi = (north ? theta + chi : theta - chi);

   const Real& l = length;

   Point_2D point_a (point.x + l*cos (psi), point.y + l*sin (psi));
   Point_2D point_c (point.x - width*cos (theta), point.y - width*sin (theta));
   Point_2D point_b (point_c.x + l*cos (psi), point_c.y + l*sin (psi));

   add (point);

   if (!invisible_barb)
   {
      add (point_a);
      add (point_b);
   }

   add (point_c);

}

void
Wind_Barb::barb (const Wind& wind,
                 const Real size,
                 const Real width)
{

   const Real speed = std::min (wind.get_speed (), Real (400.0));
//   Real theta = (450 - wind.get_direction ()) * DEGREE_TO_RADIAN;
   const Real theta = (270 + wind.get_direction ()) * DEGREE_TO_RADIAN;

   Real offset = size;
   Real length = size * 0.3;
   const Real interval = width;
   const Real w_2 = width / 2;

   const Real theta_r = theta + M_PI_2;
   const Real theta_l = theta - M_PI_2;
   Point_2D anchor_l (w_2 * cos (theta_r), w_2 * sin (theta_r));
   Point_2D anchor_r (w_2 * cos (theta_l), w_2 * sin (theta_l));

   Real residue;
   Integer pennants, barbs;
   analysis (pennants, barbs, residue, speed);

   bool no_pennants = (pennants == 0);
   bool reduced_offset = (pennants && (speed < 5));

   if (reduced_offset)
   {
      offset -= interval;
   }

   if (north)
   {
      std::swap (anchor_l, anchor_r);
   }

   while (pennants > 0)
   {

      Point_2D point (anchor_r.x + offset * cos (theta),
                      anchor_r.y + offset * sin (theta));
      add_pennant (theta, point, length);

      offset -= (interval + length * M_SQRT1_3);
      pennants--;

   }

   while (barbs > 0)
   {

      Real l = (reduced_offset && l < length * 0.4) ? length * 0.4 : length;
      l = (l * 2 / M_SQRT3);
      if (l <= w_2) { l = w_2; }

      Point_2D point (anchor_r.x + offset * cos (theta),
                      anchor_r.y + offset * sin (theta));
      add_barb (theta, point, l, width, (speed < 5));

      offset -= (interval + width);
      barbs--;

   }

   if (residue > 0)
   {

      length *= (residue / 5);
      Real l = (reduced_offset && l < length * 0.4) ? length * 0.4 : length;
      l = (l * 2 / M_SQRT3);
      if (l <= w_2) { l = w_2; }

      Point_2D point (anchor_r.x + offset * cos (theta),
                      anchor_r.y + offset * sin (theta));
      add_barb (theta, point, l, width, false);

      offset -= (interval + width);

   }

   add (anchor_r);
   add (anchor_l);

   Real s = size;
   if (no_pennants && !reduced_offset) { s -= (width / M_SQRT3); }
   add (Point_2D (anchor_l.x + s * cos (theta), anchor_l.y + s * sin (theta)));

   if (reduced_offset)
   {
      add (Point_2D (anchor_r.x + size * cos (theta),
                     anchor_r.y + size * sin (theta)));
   }

}

Wind_Barb::Wind_Barb (const Wind& wind,
                      const Real size,
                      const bool north,
                      const Real width,
                      const Real calm_threshold)
            : Symbol (size),
               north (north)
{

   Real speed = wind.get_speed ();
   Real w = (gsl_isnan (width) ? size / 20 : width);

   if (speed < calm_threshold)
   {
      Point_2D point (0, 0);
      add_circle (true, point, size / 8 + w / 2, true);
      add_circle (true, point, size / 8 - w / 2, false);
   }
   else
   {
      barb (wind, size, w);
   }

}

Wind_Transform::Wind_Transform ()
   : Polar_Transform_2D ()
{
}

Wind_Transform::Wind_Transform (const Point_2D& origin,
                                const Real scale)
   : Polar_Transform_2D (origin, scale)
{
}

void
Wind_Transform::transform (Real& x,
                           Real& y,
                           const Real direction,
                           const Real speed) const
{
   const Real theta = (direction - 90) * M_PI/180;
   Polar_Transform_2D::transform (x, y, speed, theta);
}

void
Wind_Transform::reverse (Real& direction,
                         Real& speed,
                         const Real x,
                         const Real y) const
{
   Real theta;
   Polar_Transform_2D::reverse (speed, theta, x, y);
   direction = (theta * 180/M_PI + 90);
}

Wind_Rose::Threshold::Threshold (const Real threshold,
                                 const Dstring& label_str)
   : value (threshold),
     label_str (label_str)
{
}

Wind_Rose::Record::Record (const Real threshold_value,
                           const Real percentage,
                           const Dstring& label_str)
   : threshold (threshold_value, label_str),
     percentage (percentage)
{
}

Wind_Rose::Record::Record (const Wind_Rose::Threshold& threshold,
                           const Real percentage)
   : threshold (threshold),
     percentage (percentage)
{
}

void
Wind_Rose::Arm::render (const RefPtr<Context> cr,
                        const Color& pen_color,
                        const Color_Chooser& color_chooser,
                        const Point_2D& point,
                        const Real direction,
                        const Real percentage_size,
                        const Real delta_width,
                        const bool show_label) const
{

   Polygon polygon;

   Real width = 1;
   Real radius = 0;
   Real theta = direction * DEGREE_TO_RADIAN;

   Real tl = theta + M_PI_2;
   Real tr = theta - M_PI_2;

   Real s = sin (theta), c = cos (theta);
   Real sin_tl = sin (tl), cos_tl = cos (tl);
   Real sin_tr = sin (tr), cos_tr = cos (tr);

   for (Integer i = 0; i < size (); i++)
   {

      const Wind_Rose::Record& record = at (i);
      const Real& threshold = record.threshold.value;
      const Real& percentage = record.percentage;

      if (percentage > 0)
      {

         Real w2 = width / 2;
         Real radius_h = radius * percentage_size;
         Real radius_t = (radius + percentage) * percentage_size;

         Point_2D ph (point.x + radius_h * s, point.y + radius_h * c);
         Point_2D pt (point.x + radius_t * s, point.y + radius_t * c);
         Point_2D pc ((ph.x + pt.x) / 2 + (w2 + 5) * sin_tl,
                      (ph.y + pt.y) / 2 + (w2 + 5) * cos_tl);

         if (width <= 1)
         {
            pen_color.cairo (cr);
            cr->move_to (ph.x, ph.y);
            cr->line_to (pt.x, pt.y);
            cr->stroke ();
         }
         else
         {

            polygon.clear ();
            polygon.add (Point_2D (ph.x + w2*sin_tl, ph.y + w2*cos_tl));
            polygon.add (Point_2D (ph.x + w2*sin_tr, ph.y + w2*cos_tr));
            polygon.add (Point_2D (pt.x + w2*sin_tr, pt.y + w2*cos_tr));
            polygon.add (Point_2D (pt.x + w2*sin_tl, pt.y + w2*cos_tl));


            const Color& color = color_chooser.get_color (threshold);
            color.cairo (cr);
            polygon.cairo (cr);
            cr->fill ();

            polygon.cairo (cr);
            pen_color.cairo (cr);
            cr->stroke ();

         }

         if (show_label)
         {
            Label label (record.threshold.label_str, pc, 'r', 'c');
            label.set_text_angle (M_PI - theta);
            label.cairo (cr);
         }

         radius += percentage;

      }

      width += delta_width;

   }

}

void
Wind_Rose::init (const Dstring& unit_string)
{
   this->unit_string = unit_string;
   this->multiplier = get_multiplier (unit_string);
}

Wind_Rose::Thresholds
Wind_Rose::get_thresholds (const Tuple& threshold_tuple) const
{

   Wind_Rose::Thresholds thresholds;

   for (const auto& value : threshold_tuple)
   {
      const Dstring format = "%g" + unit_string;
      const Dstring label_str = Dstring::render (format, value);
      const Wind_Rose::Threshold threshold (value, label_str);
      thresholds.push_back (threshold);
   }

   return thresholds;

}

Real
Wind_Rose::get_multiplier (const Dstring& unit_string)
{
   if (unit_string == "kt") { return 0.51444444; }
   if (unit_string == "m/s") { return 1; }
   if (unit_string == "km/h") { return 0.277777778; }
   return GSL_NAN;
}

Wind_Rose::Wind_Rose (const Integer number_of_directions,
                      const Wind_Rose::Thresholds& thresholds,
                      const Dstring& unit_string)
{
   reset (number_of_directions, thresholds, unit_string);
}

Wind_Rose::Wind_Rose (const Integer number_of_directions,
                      const Tuple& threshold_tuple,
                      const Dstring& unit_string)
{
   reset (number_of_directions, threshold_tuple, unit_string);
}

Wind_Rose::Wind_Rose (const Integer number_of_directions,
                      const Dstring& threshold_vector_string,
                      const Dstring& unit_string)
{

   if (threshold_vector_string == "beaufort")
   {

      const Real multiplier = get_multiplier (unit_string);

      thresholds.push_back (Wind_Rose::Threshold (0.51 / multiplier, "LT"));
      thresholds.push_back (Wind_Rose::Threshold (3.60 / multiplier, "MD"));
      thresholds.push_back (Wind_Rose::Threshold (8.75 / multiplier, "FS"));
      thresholds.push_back (Wind_Rose::Threshold (11.32 / multiplier, "SG"));
      thresholds.push_back (Wind_Rose::Threshold (17.49 / multiplier, "GA"));
      thresholds.push_back (Wind_Rose::Threshold (24.69 / multiplier, "SM"));
      thresholds.push_back (Wind_Rose::Threshold (32.92 / multiplier, "HC"));

      reset (number_of_directions, thresholds, unit_string);

   }
   else
   if (threshold_vector_string == "full_beaufort")
   {

      Wind_Rose::Thresholds thresholds;
      const Real multiplier = get_multiplier (unit_string);

      thresholds.push_back (Wind_Rose::Threshold (0.51 / multiplier, "F1"));
      thresholds.push_back (Wind_Rose::Threshold (1.54 / multiplier, "F2"));
      thresholds.push_back (Wind_Rose::Threshold (3.60 / multiplier, "F3"));
      thresholds.push_back (Wind_Rose::Threshold (5.66 / multiplier, "F4"));
      thresholds.push_back (Wind_Rose::Threshold (8.75 / multiplier, "F5"));
      thresholds.push_back (Wind_Rose::Threshold (11.32 / multiplier, "F6"));
      thresholds.push_back (Wind_Rose::Threshold (13.89 / multiplier, "F7"));
      thresholds.push_back (Wind_Rose::Threshold (17.49 / multiplier, "F8"));
      thresholds.push_back (Wind_Rose::Threshold (21.09 / multiplier, "F9"));
      thresholds.push_back (Wind_Rose::Threshold (24.69 / multiplier, "F10"));
      thresholds.push_back (Wind_Rose::Threshold (28.81 / multiplier, "F11"));
      thresholds.push_back (Wind_Rose::Threshold (32.92 / multiplier, "F12"));

      reset (number_of_directions, thresholds, unit_string);

   }
   else
   {
      const Tuple threshold_tuple (threshold_vector_string);
      reset (number_of_directions, threshold_tuple, unit_string);
   }

}

Wind_Rose::~Wind_Rose ()
{
   delete count_data;
}

void
Wind_Rose::reset (const Integer number_of_directions,
                  const Wind_Rose::Thresholds& thresholds,
                  const Dstring& unit_string)
{

   Integer n = number_of_directions;
   Integer m = thresholds.size ();

   this->number_of_directions = n;
   this->thresholds = thresholds;
   this->delta_direction = Real (360) / number_of_directions;

   //count_data_ptr = new Data_2D<Integer> (n, m);
   //count_data_ptr->init (Integer (0));

   count_data = new Integer[n*m];
   for (Integer i = 0; i < n*m; i++)
   {
      count_data[i] = 0;
   }

   calm_count = 0;
   total_count = 0;

   init (unit_string);

}

void
Wind_Rose::reset (const Integer number_of_directions,
                  const Tuple& threshold_tuple,
                  const Dstring& unit_string)
{
   const Wind_Rose::Thresholds& thresholds = get_thresholds (threshold_tuple);
   reset (number_of_directions, thresholds, unit_string);
}

const Integer&
Wind_Rose::get_number_of_directions () const
{
   return number_of_directions;
}

const Wind_Rose::Thresholds&
Wind_Rose::get_thresholds () const
{
   return thresholds;
}

const Wind_Rose::Threshold&
Wind_Rose::get_max_threshold () const
{
   return thresholds.back ();
}

Tuple
Wind_Rose::get_threshold_tuple () const
{

   Tuple tuple;

   for (const Threshold& threshold : thresholds)
   {
      tuple.push_back (threshold.value);
   }

   return tuple;

}

const Integer&
Wind_Rose::get_total_count () const
{
   return total_count;
}

Integer&
Wind_Rose::get_count (const Integer direction_index,
                      const Integer speed_index)
{
   if (direction_index == -1 && speed_index == -1) { return calm_count; }
   else
   {
      const Integer d = direction_index;
      const Integer s = speed_index;
      return count_data[d * thresholds.size () + s];
   }
}

const Integer&
Wind_Rose::get_count (const Integer direction_index,
                      const Integer speed_index) const
{
   if (direction_index == -1 && speed_index == -1) { return calm_count; }
   else
   {
      const Integer d = direction_index;
      const Integer s = speed_index;
      return count_data[d * thresholds.size () + s];
   }
}

Integer
Wind_Rose::get_count_stronger (const Integer direction_index,
                               const Integer speed_index) const
{
   if (direction_index == -1 && speed_index == -1)
   {
      return total_count - calm_count;
   }
   else
   {
      Integer count = 0;
      for (Integer si = speed_index; si < thresholds.size (); si++)
      {
         count += get_count (direction_index, si);
      }
      return count;
   }
}

Integer
Wind_Rose::get_count_d (const Integer direction_index) const
{
   return get_count_stronger (direction_index, 0);
}

Real
Wind_Rose::get_percentage (const Integer direction_index,
                           const Integer speed_index) const
{
   return Real (get_count (direction_index, speed_index)) / total_count * 100;
}

Real
Wind_Rose::get_percentage_stronger (const Integer direction_index,
                                    const Integer speed_index) const
{
   Integer count= get_count_stronger (direction_index, speed_index);
   return Real (count) / total_count * 100;
}

Real
Wind_Rose::get_percentage_d (const Integer direction_index) const
{
   return Real (get_count_d (direction_index)) / total_count * 100;
}

Real
Wind_Rose::get_max_percentage () const
{

   Real max_percentage = GSL_NEGINF;

   for (Integer i = 0; i < number_of_directions; i++)
   {
      Real percentage = get_percentage_d (i);
      if (percentage > max_percentage) { max_percentage = percentage; }
   }

   return max_percentage;

}

Integer
Wind_Rose::get_calm_count () const
{
   return calm_count;
}

Real
Wind_Rose::get_calm_percentage () const
{
   return Real (calm_count) / total_count * 100;
}

Index_2D
Wind_Rose::get_index (const Wind& wind) const
{

   const Real speed = wind.get_speed ();
   const Real calm_threshold = thresholds.front ().value * multiplier;
   if (speed < calm_threshold) { return Index_2D (-1, -1); }

   const Real direction = wind.get_direction ();

   const Integer n = number_of_directions;
   const Real dd = 360.0 / n;
   const Integer d = Integer (round (direction / dd)) % n;

   Integer s = 0;
   typedef vector<Wind_Rose::Threshold>::const_iterator Iterator;
   for (Iterator iterator = thresholds.begin ();
        iterator != thresholds.end (); iterator++)
   {
      if (iterator == thresholds.begin ()) { continue; }
      const Real threshold = iterator->value * multiplier;
      if (speed < threshold) { break; }
      s++;
   }

   return Index_2D (d, s);

/*
   Real direction = wind.get_direction ();
   Real aligned_direction = fmod (direction + delta_direction/2, 360);
   Integer d = Integer (floor (aligned_direction / delta_direction));
   d %= number_of_directions;

   for (Integer s = thresholds.size ()-1; s >= 0; s--)
   {

      const Real& threshold = thresholds[s].value * multiplier;
      if (speed >= threshold)
      {
         get_count (d, s)++;
         break;
      }

   }
*/

}

void
Wind_Rose::clear ()
{
   total_count = 0;
   calm_count = 0;
   for (Integer d = 0; d < number_of_directions; d++)
   {
      for (Integer s = 0; s < thresholds.size (); s++)
      {
         get_count (d, s) = 0;
      }
   }
}

void
Wind_Rose::add_wind (const Wind& wind)
{
   total_count++;
   const Index_2D& index_2d = get_index (wind);
   get_count (index_2d.i, index_2d.j)++;

   //Real speed = wind.get_speed ();
   //const Real calm_threshold = thresholds.front ().value * multiplier;
   //if (speed < calm_threshold) { calm_count++; } 
   //else
   //{
   //   const Index_2D& index_2d = get_index (wind);
   //   get_count (index_2d.i, index_2d.j)++;
   //}

}

void
Wind_Rose::add_wind (const vector<Wind>& wind_vector)
{

   for (vector<Wind>::const_iterator iterator = wind_vector.begin ();
        iterator != wind_vector.end (); iterator++)
   {
      const Wind& wind = *(iterator);
      add_wind (wind);
   }

}

Wind_Rose::Arm
Wind_Rose::get_arm (const Integer direction_index) const
{

   Wind_Rose::Arm arm;

   for (Integer s = 0; s < thresholds.size (); s++)
   {
      const Wind_Rose::Threshold& threshold = thresholds[s];
      Real percentage = get_percentage (direction_index, s);
      arm.push_back (Wind_Rose::Record (threshold, percentage));
   }

   return arm;

}

Wind_Rose::Arm
Wind_Rose::get_unit_arm () const
{

   Wind_Rose::Arm arm;

   for (Integer s = 0; s < thresholds.size (); s++)
   {
      const Wind_Rose::Threshold& threshold = thresholds[s];
      Wind_Rose::Record record (threshold, 1.0);
      arm.push_back (record);
   }

   return arm;

}

void
Wind_Rose::render (const RefPtr<Context>& cr,
                   const Point_2D& point,
                   const Color& pen_color,
                   const Color& calm_color,
                   const Color& ring_color,
                   const Color_Chooser& color_chooser,
                   const Real calm_ring_size,
                   const Real percentage_size,
                   const Real delta_width,
                   const Tuple& ring_tuple) const
{

   Real radius_h, radius_t;

   const Integer& number_of_directions = get_number_of_directions ();
   Real delta_direction = Real (360) / Real (number_of_directions);
   Real calm_percentage = get_calm_percentage ();

   const Wind_Rose::Thresholds& thresholds = get_thresholds ();

   cr->save ();
   cr->move_to (0, 0);
   cr->stroke ();
   Tuple clone_ring_tuple = ring_tuple;

   if (clone_ring_tuple.empty ())
   {

      Real max_percentage = get_max_percentage ();
      Real delta_p = get_delta (0, max_percentage, 5);

      for (Real p = delta_p; p <= max_percentage; p += delta_p)
      {
         clone_ring_tuple.push_back (p);
      }

   }

   { // plot contentric rings

      ring_color.cairo (cr);

      for (Tuple::const_iterator iterator = clone_ring_tuple.begin ();
           iterator != clone_ring_tuple.end (); iterator++)
      {

         const Real& p = *(iterator);
         if (gsl_isnan (p) || p <= 0) { continue; }

         const Dstring& label_str = Dstring::render ("%g%%", p);
         Real r = p * percentage_size + calm_ring_size;
         Point_2D label_point (point.x + 0.643*r, point.y - 0.766*r);

         cr->arc (point.x, point.y, r, 0, 2*M_PI);

         Label label (label_str, label_point, 'r', 'b');
         label.cairo (cr);

         cr->stroke ();

      } 

   }

   calm_color.cairo (cr);
   cr->arc (point.x, point.y, calm_ring_size, 0, 2*M_PI);
   cr->fill ();

   pen_color.cairo (cr);
   cr->arc (point.x, point.y, calm_ring_size, 0, 2*M_PI);
   cr->stroke ();

   const Dstring& calm_percentage_str =
      Dstring::render ("%.1f%%", calm_percentage);

   Label label (calm_percentage_str, point, 'c', 'c');
   label.cairo (cr);

   for (Integer d = 0; d < number_of_directions; d++)
   {

      Real width = 1;
      Real radius = 0;
      Real direction = (180 - d * delta_direction);
      Real theta = direction * DEGREE_TO_RADIAN;

      Point_2D p (point.x + calm_ring_size * sin (theta),
                  point.y + calm_ring_size * cos (theta));

      Wind_Rose::Arm arm = get_arm (d);

      arm.render (cr, pen_color, color_chooser, p,
         direction, percentage_size, delta_width, false);

   }

   cr->restore ();

}

void
Wind_Rose::render (const RefPtr<Context>& cr,
                   const Transform_2D& transform_2d,
                   const Point_2D& point,
                   const Color& pen_color,
                   const Color& calm_color,
                   const Color& ring_color,
                   const Color_Chooser& color_chooser,
                   const Real calm_ring_size,
                   const Real percentage_size,
                   const Real delta_width,
                   const Tuple& ring_tuple) const
{
   render (cr, transform_2d.transform (point), pen_color,
      calm_color, ring_color, color_chooser, calm_ring_size,
      percentage_size, delta_width, ring_tuple);
}

Wind_Disc::Transform::Transform (const Wind_Disc& wind_disc,
                                 const Point_2D& origin,
                                 const Real max_speed,
                                 const Real max_radius,
                                 const Real calm_radius,
                                 const Real label_height)
   : wind_disc (wind_disc),
     origin (origin),
     max_speed (max_speed),
     max_radius (max_radius),
     calm_radius (calm_radius),
     label_height (label_height)
{
   const Wind_Rose::Thresholds& thresholds = wind_disc.get_thresholds ();
   const Real calm_threshold = thresholds.front ().value;
   const Real speed_span = max_speed - calm_threshold;
   const Real r_span = (max_radius - calm_radius - label_height);
   this->scale = (r_span / speed_span);
}

bool
Wind_Disc::Transform::out_of_domain (const Real x,
                                     const Real y) const
{
   const Real dx = (x - origin.x);
   const Real dy = -(y - origin.y);
   const Real radius = sqrt (dx*dx + dy*dy);
   return (radius > max_radius);
}

void
Wind_Disc::Transform::transform (Real& x,
                                 Real& y,
                                 const Real direction,
                                 const Real speed) const
{
   const Real r = get_radius (speed);
   const Real theta = direction * DEGREE_TO_RADIAN;
   x = origin.x + r * sin (theta);
   y = origin.y - r * cos (theta);
}
 
void
Wind_Disc::Transform::reverse (Real& direction,
                               Real& speed,
                               const Real x, 
                               const Real y) const
{
   const Real dx = (x - origin.x);
   const Real dy = -(y - origin.y);
   const Real radius = sqrt (dx*dx + dy*dy);
   const Real d = atan2 (dx, dy) * RADIAN_TO_DEGREE;
   direction = (d < 0 ? d + 360 : d);
   speed = get_speed (radius);
}

const Point_2D&
Wind_Disc::Transform::get_origin () const
{
   return origin;
}

const Real
Wind_Disc::Transform::get_max_speed () const
{
   return max_speed;
}

const Real
Wind_Disc::Transform::get_max_radius () const
{
   return max_radius;
}

const Real
Wind_Disc::Transform::get_calm_radius () const
{
   return calm_radius;
}

const Real
Wind_Disc::Transform::get_label_height () const
{
   return label_height;
}

Real
Wind_Disc::Transform::get_radius (const Real speed) const
{
   const Wind_Rose::Thresholds& thresholds = wind_disc.get_thresholds ();
   const Real calm_threshold = thresholds.front ().value;
   return calm_radius + scale * (speed - calm_threshold);
}

Real
Wind_Disc::Transform::get_speed (const Real radius) const
{
   if (radius < calm_radius) { return 0; }
   else
   {
      const Wind_Rose::Thresholds& thresholds = wind_disc.get_thresholds ();
      const Real calm_threshold = thresholds.front ().value;
      return (radius - calm_radius) / scale + calm_threshold;
   }
}

void
Wind_Disc::render_background (const RefPtr<Context> cr) const
{
   cr->set_fill_rule (FILL_RULE_EVEN_ODD);
   render_directions (cr);
   render_mesh (cr);
}

void
Wind_Disc::render_directions (const RefPtr<Context> cr) const
{

   if (is_even (number_of_directions))
   {
      render_directions_even (cr);
   }
   else
   {
      render_directions_odd (cr);
   }

   render_direction_labels (cr);

}

void
Wind_Disc::render_direction_labels (const RefPtr<Context> cr) const
{

   const Color& color = major_color;

   cr->save ();
   color.cairo (cr);

   const Real max_speed = transform_ptr->get_max_speed ();
   const Real label_height = transform_ptr->get_label_height ();
   const Real delta_theta = 2*M_PI / number_of_directions;

   for (Integer i = 0; i < number_of_directions; i++)
   {

      const Real theta = i * delta_theta;
      const Real d = theta * RADIAN_TO_DEGREE;

      const Point_2D p (d, max_speed);
      const Wind& wind = Wind::direction_speed (d, max_speed);
      const Dstring& str = wind.get_direction_string (number_of_directions);

      Label label (str, p, 'c', 'c');
      label.set_text_angle (theta);
      label.set_offset (Point_2D (0, -label_height / 2));
      label.cairo (cr, *transform_ptr);

   }

   cr->restore ();

}

void
Wind_Disc::render_directions_even (const RefPtr<Context> cr) const
{

   const Real delta_theta = 2*M_PI / number_of_directions;
   const Real dt = delta_theta / 2;

   const Point_2D& origin = transform_ptr->get_origin ();
   const Real max_speed = transform_ptr->get_max_speed ();
   const Real label_height = transform_ptr->get_label_height ();

   const Real outer_r = transform_ptr->get_radius (max_speed) + label_height;
   const Real middle_r = transform_ptr->get_radius (max_speed);
   const Real inner_r = transform_ptr->get_calm_radius ();

   const Color& color = shade_color;

   cr->save ();
   color.cairo (cr);

   for (Integer i = 0; i < number_of_directions / 2; i++)
   {

      const Real theta = i * 2*delta_theta - M_PI/2;

      const Real x = origin.x + middle_r * cos (theta - dt);
      const Real y = origin.y + middle_r * sin (theta - dt);

      cr->move_to (x, y);
      cr->arc (origin.x, origin.y, middle_r, theta - dt, theta + 3*dt);
      cr->arc_negative (origin.x, origin.y, outer_r, theta + 3*dt, theta + dt);
      cr->arc_negative (origin.x, origin.y, inner_r, theta + dt, theta - dt);
      cr->fill ();

   }

   cr->restore ();

}

void
Wind_Disc::render_directions_odd (const RefPtr<Context> cr) const
{

   const Real delta_theta = 2*M_PI / number_of_directions;
   const Real dt = delta_theta / 2;

   const Point_2D& origin = transform_ptr->get_origin ();
   const Real max_speed = transform_ptr->get_max_speed ();
   const Real label_height = transform_ptr->get_label_height ();

   const Real outer_r = transform_ptr->get_radius (max_speed) + label_height;
   const Real inner_r = transform_ptr->get_calm_radius ();

   const Color& color = shade_color;

   cr->save ();
   color.cairo (cr);

   for (Integer i = 0; i < number_of_directions; i++)
   {

      const Real theta = (i + 0.5) * delta_theta - M_PI/2;

      const Real x_0 = origin.x + inner_r * cos (theta);
      const Real y_0 = origin.y + inner_r * sin (theta);
      const Real x_1 = origin.x + outer_r * cos (theta);
      const Real y_1 = origin.y + outer_r * sin (theta);

      const Edge edge (Point_2D (x_0, y_0), Point_2D (x_1, y_1));
      edge.cairo (cr);
      cr->stroke ();

   }

   cr->restore ();

}

Polygon*
Wind_Disc::render_ring_label (const RefPtr<Context> cr) const
{

   const Color& color = major_color;

   Point_2D point;

   Polygon* clip_polygon_ptr = new Polygon ();
   Polygon& clip_polygon = *clip_polygon_ptr;

   const Point_2D& origin = transform_ptr->get_origin ();
   const Real max_speed = transform_ptr->get_max_speed ();
   const Real label_height = transform_ptr->get_label_height ();
   const Real outer_r = transform_ptr->get_radius (max_speed) + label_height;

   for (Integer i = 0; i < 360; i++)
   {
      const Real theta = Real (i) * DEGREE_TO_RADIAN;
      point.x = origin.x + outer_r * cos (theta);
      point.y = origin.y + outer_r * sin (theta);
      clip_polygon.add (point);
   }

   const Real delta_d = 360.0 / number_of_directions;

   cr->save ();
   cr->set_fill_rule (FILL_RULE_EVEN_ODD);
   cr->set_font_size (label_height * 0.5 * 0.9);
   color.cairo (cr);
   //cr->select_font_face ("AR PL Mingti2L Big5", FONT_SLANT_NORMAL, FONT_WEIGHT_BOLD);

   const Tuple& threshold_tuple = get_threshold_tuple ();

   const Integer number_of_labels = 8;
   const Real delta_d_label = 360.0 / number_of_labels;
   const Dstring& format = "%.0f" + unit_string;

   for (Integer j = 0; j < number_of_labels; j++)
   {

      const Real direction = j * delta_d_label + delta_d / 2;
      const Real theta = direction * DEGREE_TO_RADIAN;

      for (Integer i = 0; i < speed_label_tuple.size (); i++)
      { 

         const Real speed = speed_label_tuple[i];
         const Point_2D p (direction, speed);

         const Dstring& str = Dstring::render (format, speed);
         Label label (str, p, 'c', 'c');
         label.set_text_angle (theta);
         label.cairo (cr, *transform_ptr);

         const Rect& rect = (const Rect&)label;
         clip_polygon.add (rect);

      }

      if (speed_label_tuple.size () == 0)
      {

         typedef Wind_Rose::Threshold Wrt;

         // no need to draw calm circle
         for (Integer i = 1; i < thresholds.size (); i++)
         { 

            const Wrt& wind_rose_threshold = thresholds[i];
            const Real speed = wind_rose_threshold.value;
            const Dstring& str = wind_rose_threshold.label_str;

            const Point_2D p (direction, speed);

            Label label (str, p, 'c', 'c');
            label.set_text_angle (theta);
            label.cairo (cr, *transform_ptr);

            const Rect& rect = (const Rect&)label;
            clip_polygon.add (rect);

         }

      }

   }

   cr->restore ();

   return clip_polygon_ptr;

}

void
Wind_Disc::render_mesh (const RefPtr<Context> cr) const
{

   const Point_2D& origin = transform_ptr->get_origin ();
   const Real label_height = transform_ptr->get_label_height ();
   const Real max_speed = transform_ptr->get_max_speed ();
   const Real middle_r = transform_ptr->get_radius (max_speed);
   const Real outer_r = middle_r + label_height;

   Polygon* clip_polygon_ptr = render_ring_label (cr);
   const Polygon& clip_polygon = *clip_polygon_ptr;

   cr->save ();

   clip_polygon.cairo (cr);
   cr->clip ();

   const Real calm_threshold = thresholds.front ().value;
   const Integer minor_i_start = Integer (ceil (calm_threshold));
   const Integer minor_i_end = Integer (floor (max_speed));
   const Integer minor_n = minor_i_end - minor_i_start + 1;;
   const Tuple minor_tuple (minor_n, minor_i_start, minor_i_end);

   render_mesh (cr, minor_color, minor_tuple, 1);
   render_mesh (cr, middle_color, get_threshold_tuple (), 1.5);

   major_color.cairo (cr);

   for (Integer i = 0; i < speed_label_tuple.size (); i++)
   { 
      const Real speed = speed_label_tuple[i];
      const Real radius = transform_ptr->get_radius (speed);
      const Ring ring (radius);
      ring.cairo (cr, origin);
      cr->stroke ();
   }

   cr->set_line_width (5);
   Ring (outer_r).cairo (cr, origin);
   Ring (middle_r).cairo (cr, origin);
   cr->stroke ();

   cr->restore ();

   delete clip_polygon_ptr;

}

void
Wind_Disc::render_mesh (const RefPtr<Context> cr,
                        const Color& color,
                        const Tuple& speed_tuple,
                        const Real line_width) const
{

   const Point_2D& origin = transform_ptr->get_origin ();

   cr->save ();
   color.cairo (cr);
   cr->set_line_width (line_width);

   for (Integer i = 0; i < speed_tuple.size (); i++)
   {

      const Real& speed = speed_tuple[i];
      const Real radius = transform_ptr->get_radius (speed);

      const Ring ring (radius);
      ring.cairo (cr, origin);
      cr->stroke ();
   }

   cr->restore ();


}

void
Wind_Disc::render_scatter_plot (const RefPtr<Context> cr,
                                const Real hue,
                                const Real dir_scatter) const
{

   // All speed units here in designated unit_string unit

   Real alpha = 50.0 / get_total_count ();
   if (alpha < 0.04) { alpha = 0.04; }
   if (alpha > 0.30) { alpha = 0.30; }

   const Ring ring (scatter_ring_size);
   const Color& color = Color::hsb (hue, 0.4, 0.8, alpha);
   const Color& color_a2 = Color::hsb (hue, 0.4, 0.8, alpha * 2);

   const Transform_2D& t = *transform_ptr;
   const Real max_speed = transform_ptr->get_max_speed ();
   const Real calm_threshold = thresholds.front ().value;

   bool calm = false;
   bool out_of_bounds = false;

   for (const Wind& wind : wind_vector)
   {

      const Real speed = wind.get_speed () / multiplier;

      if (wind.is_naw ()) { continue; }
      if (speed < calm_threshold) { calm = true; continue; }
      if (speed > max_speed) { out_of_bounds = true; continue; }

      const Real r = random (dir_scatter, -dir_scatter);
      const Real direction = wind.get_direction () + r;

      ring.cairo (cr, t.transform (Point_2D (direction, speed)));
      color.cairo (cr);
      cr->fill_preserve ();
      color_a2.cairo (cr);
      cr->stroke ();

   }

   if (calm)
   {
      const Real r = transform_ptr->get_radius (calm_threshold);
      render_scatter_ring (cr, r, scatter_ring_size);
   }

   if (out_of_bounds)
   {
      const Real r = transform_ptr->get_radius (max_speed);
      render_scatter_ring (cr, r, scatter_ring_size);
   }

}

void
Wind_Disc::render_scatter_ring (const RefPtr<Context> cr,
                                const Real radius,
                                const Real line_width) const
{

   cr->save ();
   cr->set_line_width (line_width);
   cr->set_line_cap (LINE_CAP_ROUND);

   const Point_2D& origin = transform_ptr->get_origin ();

   Dashes (Tuple (2, line_width * 2, line_width * 1.5)).cairo (cr);
   Ring (radius).cairo (cr, origin);
   cr->stroke ();

   cr->restore ();

}

void
Wind_Disc::render_percentages (const RefPtr<Context> cr) const
{

   const Transform_2D& t = *transform_ptr;
   const Real max_speed = transform_ptr->get_max_speed ();
   const Real d_direction = 360.0 / number_of_directions;

   for (Integer i = 0; i < number_of_directions; i++)
   {

      const Real direction = i * d_direction;

      for (Integer j = 0; j < thresholds.size (); j++)
      {

         Real lower_speed = thresholds[j].value;
         if  (lower_speed >= max_speed) { break; }

         Real upper_speed = thresholds[j+1].value;
         const Real truncated = (max_speed < upper_speed);
         const Real strongest_speed = (j == thresholds.size () - 1);
         if (truncated || strongest_speed)
         {
            upper_speed = max_speed;
         }

         const Real speed = (lower_speed + upper_speed) / 2;

         const Point_2D& p = t.transform (Point_2D (direction, speed));
         const Real percentage = truncated ?
            get_percentage_stronger (i, j) : get_percentage (i, j);

         render_percentage (cr, p, percentage);

      }

   }

   const Point_2D& origin = transform_ptr->get_origin ();
   const Real calm_percentage = get_calm_percentage ();
   render_percentage (cr, origin, calm_percentage);

}

void
Wind_Disc::render_percentage (const RefPtr<Context> cr,
                              const Point_2D& point,
                              const Real percentage) const
{

   const bool b = gsl_isnan (percentage);
   const Dstring& str = (b ? "-" : Dstring::render ("%.0f%%", percentage));
   Label label (str, point, 'c', 'c');

   const Color& color_j = major_color;
   const Color& color_n = minor_color;

   color_n.cairo (cr);
   label.set_offset (Point_2D (2, -2));
   label.cairo (cr);

   color_j.cairo (cr);
   label.set_offset (Point_2D (0, 0));
   label.cairo (cr);

}

void
Wind_Disc::render_percentage_d (const RefPtr<Context> cr,
                                const Real hue) const
{

   Real max_p = GSL_NEGINF;
   Real* percentages = new Real[number_of_directions];

   for (Integer i = 0; i < number_of_directions; i++)
   {
      const Real p = get_percentage_d (i);
      percentages[i] = p;
      if (p > max_p) { max_p = p; }
   }

   max_p *= 1.1;

   const Point_2D& origin = transform_ptr->get_origin ();
   const Real max_speed = transform_ptr->get_max_speed ();
   const Real max_r = transform_ptr->get_radius (max_speed);
   const Real min_r = transform_ptr->get_calm_radius ();
   const Real m = (max_r - min_r) / max_p;

   const Transform_2D& t = *transform_ptr;
//   spline.set_interpolation (gsl_interp_cspline_periodic);

   const Color& color = Color::hsb (hue, 0.8, 0.5, 0.5);

   cr->save ();
   cr->set_line_width (6);
   cr->set_fill_rule (FILL_RULE_EVEN_ODD);
   color.cairo (cr);

   const Real delta_p = 360.0 / number_of_directions;
   for (Integer i = 0; i < number_of_directions; i++)
   {

      const Real d = i * delta_p - 90;
      const Real dl = d - delta_p/2;
      const Real dr = d + delta_p/2;
      const Real theta_l = dl * DEGREE_TO_RADIAN;
      const Real theta_r = dr * DEGREE_TO_RADIAN;

      const Real r = percentages[i] * m + min_r;
      cr->arc (origin.x, origin.y, r, theta_l, theta_r);

   }

   delete[] percentages;

   cr->close_path ();
   cr->stroke ();
   cr->restore ();

}

Wind_Disc::Wind_Disc (const Integer number_of_directions,
                      const Tuple& threshold_tuple,
                      const Point_2D& origin,
                      const Real max_radius,
                      const Tuple& speed_label_tuple,
                      const Real max_speed,
                      const Color& major_color,
                      const Color& middle_color,
                      const Color& minor_color,
                      const Color& shade_color,
                      const Real font_size,
                      const Real scatter_ring_size,
                      const Real calm_radius,
                      const Real label_height,
                      const Dstring& unit_string)
   : Wind_Rose (number_of_directions, threshold_tuple, unit_string),
     speed_label_tuple (speed_label_tuple),
     major_color (major_color),
     middle_color (middle_color),
     minor_color (minor_color),
     shade_color (shade_color),
     font_size (font_size),
     scatter_ring_size (scatter_ring_size),
     transform_ptr (NULL)
{

   const Real ms_given = !gsl_isnan (max_speed);
   const Real ms = (ms_given ? max_speed : get_max_threshold ().value * 1.15);

   transform_ptr = new Wind_Disc::Transform (*this, origin, ms,
      max_radius, calm_radius, label_height);

}

Wind_Disc::~Wind_Disc ()
{
   delete transform_ptr;
}

void
Wind_Disc::set (const Integer number_of_directions,
                const Tuple& threshold_tuple,
                const Tuple& speed_label_tuple,
                const Real max_speed)
{
   //Wind_Rose::set (number_of_directions, thresdhold_tuple); 
   //this->speed_label_tuple = speed_lable_tuple;
   //this->max_speed = max_speed;
}

void
Wind_Disc::set_position (const Point_2D& origin,
                         const Real max_radius)
{

   const Real max_speed = transform_ptr->get_max_speed ();
   const Real calm_radius = transform_ptr->get_calm_radius ();
   const Real label_height = transform_ptr->get_label_height ();

   if (transform_ptr != NULL) { delete transform_ptr; }

   transform_ptr = new Wind_Disc::Transform (*this, origin, max_speed,
      max_radius, calm_radius, label_height);

}

void
Wind_Disc::clear ()
{
   Wind_Rose::clear ();
   wind_vector.clear ();
}

void
Wind_Disc::add_wind (const Wind& wind)
{
   wind_vector.push_back (wind);
   Wind_Rose::add_wind (wind);
}

Wind
Wind_Disc::get_wind (const Point_2D& point) const
{

   if (transform_ptr->out_of_domain (point.x, point.y))
   {
      throw Exception ("Wind_Disc::get_wind out of domain");
   }

   Real direction, speed;
   transform_ptr->reverse (direction, speed, point.x, point.y);
   return Wind::direction_speed (direction, speed * multiplier);

}

void
Wind_Disc::render (const RefPtr<Context> cr,
                   const Real hue,
                   const bool outline,
                   const Real dir_scatter) const
{

   const Point_2D& origin = transform_ptr->get_origin ();
   const Real max_radius = transform_ptr->get_max_radius ();

   cr->save ();
   cr->set_font_size (font_size);
   Color (1, 1, 1, 1).cairo (cr);
   Ring (max_radius).cairo (cr, origin);
   cr->fill ();
   
   render_background (cr);
   render_scatter_plot (cr, hue, dir_scatter);

   if (outline) { render_percentage_d (cr, hue); }
   render_percentages (cr);

   cr->restore ();

}

void
Wind_Disc::render_bg (const RefPtr<Context> cr) const
{

   const Point_2D& origin = transform_ptr->get_origin ();
   const Real max_radius = transform_ptr->get_max_radius ();

   cr->save ();
   cr->set_font_size (font_size);
   Color (1, 1, 1, 1).cairo (cr);
   Ring (max_radius).cairo (cr, origin);
   cr->fill ();
   
   render_background (cr);

   cr->restore ();

}

void
Wind_Disc::render_index (const RefPtr<Context> cr,
                         const Index_2D& index) const
{

   const Integer s = index.j;

   const Point_2D& origin = transform_ptr->get_origin ();
   const Color& color_0 = Color::hsb (0.0, 0.6, 0.6, 0.3);
   const Color& color_1 = Color::hsb (0.1, 0.5, 0.6, 0.3);

   if (index.i == -1 && index.j == -1)
   {
      cr->save ();
      //Color (0.5, 0.5, 1, 0.3).cairo (cr);
      Stripped (color_0, color_1, 20).cairo (cr);
      const Real calm_radius = transform_ptr->get_calm_radius ();
      Circle (origin, calm_radius).cairo (cr); 
      cr->fill ();
      cr->restore ();
      return;
   }

   const bool top_speed = (s >= thresholds.size () - 1);
   const Real max_speed = transform_ptr->get_max_speed ();
   const Real inner_speed = thresholds[s].value;
   const Real outer_speed = (top_speed ? max_speed : thresholds[s + 1].value);
   const Real inner_r = transform_ptr->get_radius (inner_speed);
   const Real outer_r = transform_ptr->get_radius (outer_speed);

   const Real d_direction = 360.0 / number_of_directions;
   const Real direction = index.i * d_direction;
   const Real left_direction = direction - d_direction / 2;
   const Real right_direction = direction + d_direction / 2;


   Point_2D lo, ri;
   transform_ptr->transform (lo.x, lo.y, left_direction, outer_speed);
   transform_ptr->transform (ri.x, ri.y, right_direction, inner_speed);

   const Real theta_l = atan2 (lo.y - origin.y, lo.x - origin.x);
   const Real theta_r = atan2 (ri.y - origin.y, ri.x - origin.x);

   cr->save ();
   Stripped (color_0, color_1, 20).cairo (cr);
   cr->move_to (lo.x, lo.y);
   cr->arc (origin.x, origin.y, outer_r, theta_l, theta_r);
   cr->arc_negative (origin.x, origin.y, inner_r, theta_r, theta_l);
   //cr->line_to (p.x, p.y);
   cr->fill ();
   cr->restore ();

}

const Point_2D&
Wind_Disc::get_origin () const
{
   return transform_ptr->get_origin ();
}

const Real
Wind_Disc::get_max_radius () const
{
   return transform_ptr->get_max_radius ();
}

const Wind_Disc::Transform&
Wind_Disc::get_transform () const
{
   return *transform_ptr;
}

Axisymmetric_Vortex::Axisymmetric_Vortex (const Real inflow)
   : center (Point_2D (0, 0)),
     inflow (inflow)
{
}

Axisymmetric_Vortex::Axisymmetric_Vortex (const Point_2D& center,
                                          const Real inflow)
   : center (center),
     inflow (inflow)
{
}

const Point_2D&
Axisymmetric_Vortex::get_center () const
{
   return center;
}

const Real&
Axisymmetric_Vortex::get_inflow () const
{
   return inflow;
}

Wind
Axisymmetric_Vortex::get_wind (const Real x,
                               const Real y) const
{

   const Real dx = x - center.x;
   const Real dy = y - center.y;

   const Real r = sqrt (dx*dx + dy*dy);
   const Real theta = atan2 (dy, dx) + inflow + M_PI/2;

   const Real V = get_V (r);
   return Wind (V * cos (theta), V * sin (theta));

}

Wind
Axisymmetric_Vortex::get_wind (const Point_2D& point) const
{
   return get_wind (point.x, point.y);
}

Real
Axisymmetric_Vortex::get_radial_wind (const Point_2D& radar,
                                      const Real x,
                                      const Real y) const
{

   const Real dx = x - radar.x;
   const Real dy = y - radar.y;
   const Real theta = atan2 (dy, dx);

   const Wind& wind = get_wind (x, y);
   return wind.u * cos (theta) + wind.v * sin (theta);

}

Real
Axisymmetric_Vortex::get_radial_wind (const Point_2D& radar,
                         const Point_2D& point) const
{
   return get_radial_wind (radar, point.x, point.y);
}

Circle
Axisymmetric_Vortex::get_isodop (const Point_2D& radar_center) const
{

   const Real m = 1.0 / tan (get_inflow ());
   const Real D = m * (radar_center.y - center.y) - radar_center.x - center.x;
   const Real E = m * (center.x - radar_center.x) - radar_center.y - center.y;
   const Real F = radar_center.x * center.x + radar_center.y * center.y +
             m * (radar_center.x * center.y - radar_center.y * center.x);

   return Circle (1.0, D, E, F);

//   const Real t = tan (get_infow () + M_PI_2);
//   const Real x_v_tilde = centre.x - radar_center.x;
//   const Real y_v_tilde = centre.y - radar_center.y;
//   const Real A = (x_v_tilde - t * y_v_tilde) / 2;
//   const Real B = (t * x_v_tilde + y_v_tilde) / 2;
//   const Real radius = sqrt (A*A + B*B);
//   const Point_2D centre (radar_center.x + A, radar_centre.y + B);
//   return Circle (centre, radius);
   

}

Exponential_Vortex::Exponential_Vortex (const Real max_wind,
                                        const Real decay,
                                        const Real inflow)
   : Axisymmetric_Vortex (inflow),
     max_wind (max_wind),
     decay (decay)
{
}

Exponential_Vortex::Exponential_Vortex (const Point_2D& center,
                                        const Real max_wind,
                                        const Real decay,
                                        const Real inflow)
     : Axisymmetric_Vortex (center, inflow),
       max_wind (max_wind),
       decay (decay)
{
}

Real
Exponential_Vortex::get_V (const Real r) const
{
   return max_wind * exp (-fabs (r) / decay);
}

const Real&
Exponential_Vortex::get_max_wind () const
{
   return max_wind;
}

const Real&
Exponential_Vortex::get_decay () const
{
   return decay;
}

Rankine_Vortex::Rankine_Vortex (const Real max_wind,
                                const Real max_radius,
                                const Real inflow)
   : Axisymmetric_Vortex (inflow),
     k_in (max_wind / max_radius),
     k_out (max_wind * max_radius),
     max_radius (max_radius)
{
}

Rankine_Vortex::Rankine_Vortex (const Point_2D& center,
                                const Real max_wind,
                                const Real max_radius,
                                const Real inflow)
   : Axisymmetric_Vortex (center, inflow),
     k_in (max_wind / max_radius),
     k_out (max_wind * max_radius),
     max_radius (max_radius)
{
}

Real
Rankine_Vortex::get_V (const Real r) const
{
   Real fabs_r = fabs (r);
   return (fabs_r < max_radius ? k_in * fabs_r : k_out / fabs_r);
}

Real
Rankine_Vortex::get_max_wind () const
{
   return k_in * max_radius;
}

const Real&
Rankine_Vortex::get_max_radius () const
{
   return max_radius;
}

Real
Fire::get_gfdi (const Real t,
                const Real rh,
                const Real kph,
                const Real curing)
{
   const Real a = 0.009254;
   const Real b = -0.004096 * pow (double (100 - curing), double (1.536));
   const Real c = 0.01201 * t;
   const Real d = 0.2789 * sqrt (kph);
   const Real e = -0.09577 * sqrt (rh);
   return pow (10, a + b + c + d + e);
}

Real
Fire::get_gfdi_si (const Real temperature,
                   const Real dew_point,
                   const Real speed,
                   const Real curing)
{
   const Real t = temperature - K;
   const Real rh = Moisture::get_rh (t, dew_point - denise::K) * 100;
   const Real kph = speed * 3.6;
   return get_gfdi (t, rh, kph, curing);
}

Real
Fire::get_ffdi (const Real t,
                const Real rh,
                const Real kph,
                const Real df)
{
   const Real a = -0.45;
   const Real b = 0.987 * log (df);
   const Real c = -0.0345 * rh;
   const Real d = 0.0338 * t;
   const Real e = 0.0234 * kph;
   return 2 * exp (a + b + c + d + e);
}

Real
Fire::get_ffdi_si (const Real temperature,
                   const Real dew_point,
                   const Real speed,
                   const Real df)
{
   const Real t = temperature - K;
   const Real rh = Moisture::get_rh (t, dew_point - denise::K) * 100;
   const Real kph = speed * 3.6;
   return get_ffdi (t, rh, kph, df);
}

Level::Level ()
   : type (NAL),
     value (GSL_NAN),
     value_ (GSL_NAN)
{
}

Level::Level (const Level& level)
   : type (level.type),
     value (level.value),
     value_ (level.value_)
{
}

Level::Level (const Dstring& str)
   : value (GSL_NAN),
     value_ (GSL_NAN)
{
   if (str == "Screen")
   {
      type = SCREEN;
   }
   else
   if (str == "10m")
   {
      type = TEN_METRE;
   }
   else
   if (str == "50m")
   {
      type = FIFTY_METRE;
   }
   else
   if (str == "MS")
   {
      type = MEAN_SEA;
   }
   else
   if (str == "Nil")
   {
      type = NIL;
   }
   else
   if (str == "Surface")
   {
      type = SURFACE;
   }
   else
   if (str.find ("magl") != Dstring::npos)
   {
      type = MAGL;
      value = stof (str);
   }
   else
   if (str.find ("m") != Dstring::npos)
   {
      type = HEIGHT;
      value = stof (str);
   }
   else
   if (str.find ("hPa") != Dstring::npos)
   {
      type = PRESSURE;
      value = stof (str) * 1e2;
   }
   else
   if (str.find ("K") != Dstring::npos)
   {
      type = THETA;
      value = stof (str);
   }
   else
   if (str.find ("M") != Dstring::npos)
   {
      type = MODEL;
      value = stof (str.substr (1));
   }
   else
   {
      type = SIGMA;
      value = stof (str);
   }
}

Level::Level (const Level::Type type,
              const Real value)
   : type (type),
     value (value),
     value_ (GSL_NAN)
{
}

Level::Level (const Level::Type type,
              const Real value_0,
              const Real value_1)
   : type (type),
     value (value_0),
     value_ (value_1)
{
}

bool
Level::is_layer () const
{
   return !gsl_isnan (value_);
}

void
Level::order ()
{
   if (gsl_isnan (value) || gsl_isnan (value_)) { return; }
   if (value > value_) { std::swap (value, value_); }
}

Level
Level::theta_level (const Real theta)
{
   return Level (THETA, theta, GSL_NAN);
}

Level
Level::sigma_level (const Real sigma)
{
   return Level (SIGMA, sigma, GSL_NAN);
}

Level
Level::pressure_level (const Real pressure)
{
   return Level (PRESSURE, pressure, GSL_NAN);
}

Level
Level::z_level (const Real z)
{
   return Level (HEIGHT, z, GSL_NAN);
}

Level
Level::model_level (const Real m)
{
   return Level (MODEL, m, GSL_NAN);
}

Level
Level::screen_level ()
{
   return Level (SCREEN);
}
Level
Level::fifty_metre_level ()
{
   return Level (FIFTY_METRE);
}

Level
Level::ten_metre_level ()
{
   return Level (TEN_METRE);
}

Level
Level::mean_sea_level ()
{
   return Level (MEAN_SEA);
}

Level
Level::nil_level ()
{
   return Level (NIL);
}

Level
Level::surface_level ()
{
   return Level (SURFACE);
}

void
Level::set_height (const Real z)
{
   this->type = HEIGHT;
   this->value = z;
}

void
Level::set_pressure (const Real pressure)
{
   this->type = PRESSURE;
   this->value = pressure;
}

void
Level::set_theta (const Real theta)
{
   this->type = THETA;
   this->value = theta;
}

void
Level::set_sigma (const Real sigma)
{
   this->type = SIGMA;
   this->value = sigma;
}

void
Level::set_model (const Real m)
{
   this->type = MODEL;
   this->value = m;
}

void
Level::set_screen ()
{
   this->type = SCREEN;
}
void
Level::set_fifty_metre ()
{
   this->type = FIFTY_METRE;
}

void
Level::set_ten_metre ()
{
   this->type = TEN_METRE;
}

void
Level::set_mean_sea ()
{
   this->type = MEAN_SEA;
}

void
Level::set_nil ()
{
   this->type = NIL;
}

void
Level::set_surface ()
{
   this->type = SURFACE;
}

void
Level::set_value (const Real value)
{
   this->value = value;
}

Real
Level::get_value () const
{
   return value;
}

Dstring
Level::get_string () const
{

   if (gsl_isnan (value_))
   {

      switch (type)
      {
         case THETA:
            return Dstring::render ("%.0fK", round (value));
         case SIGMA:
            return Dstring::render ("%.4f", value);
         case PRESSURE:
            return Dstring::render ("%.0fhPa", round (value * 1e-2));
         case HEIGHT:
            return Dstring::render ("%.0fm", round (value));
         case MAGL:
            return Dstring::render ("%.0fmagl", round (value));
         case MODEL:
            return Dstring::render ("M%.0f", round (value));
         case SCREEN:
            return "Screen";
         case FIFTY_METRE:
            return "50m";
         case TEN_METRE:
            return "10m";
         case MEAN_SEA:
            return "MS";
         case NIL:
            return "";
         case SURFACE:
            return "Surface";
         case NAL:
            return "";
      }

   }
   else
   {

      switch (type)
      {

         case THETA:
         {
            const Real theta_0 = round (value);
            const Real theta_1 = round (value_);
            return Dstring::render ("%.0fK to %.0fK", theta_0, theta_1);
         }

         case SIGMA:
         {
            return Dstring::render ("%.4f to %.4f", value, value_);
         }

         case PRESSURE:
         {
            const Real p_0 = round (value * 1e-2);
            const Real p_1 = round (value_ * 1e-2);
            return Dstring::render ("%.0fhPa to %.0fhPa", p_0, p_1);
         }

         case HEIGHT:
         {
            return Dstring::render ("%.0fm to %.0fm", value, value_);
         }

         case MODEL:
         {
            return Dstring::render ("M%.0f to M%.0f", value, value_);
         }

      }

   }

}

Level
Level::get_level_0 () const
{
   return Level (type, value);
}

Level
Level::get_level_1 () const
{
   return Level (type, value_);
}

void
Level::set_nal ()
{
   type = NAL;
}

bool
Level::is_nal () const
{
   return (type == NAL);
}

ostream&
Level::operator<< (ostream& o) const
{
   o << get_string ();
   return o;
}

Layer::Layer (const Level::Type type,
              const Real value_0,
              const Real value_1)
   : Level (type, value_0, value_1)
{
}

void
Layer::set (const Real value_0,
            const Real value_1)
{
   this->value = value;
   this->value_ = value_1;
}

Z_Layer::Z_Layer (const Z_Layer& z_layer)
   : Level (HEIGHT, z_layer.get_start_z (), z_layer.get_end_z ())
{
   order ();
}

Z_Layer::Z_Layer (const Real start_z,
                  const Real end_z)
   : Level (HEIGHT, start_z, end_z)
{
   order ();
}

void
Z_Layer::set (const Real start_z,
              const Real end_z)
{
   this->value = start_z;
   this->value_ = end_z;
   order ();
}

void
Z_Layer::set_start_z (const Real start_z)
{
   this->value = start_z;
}

void
Z_Layer::set_end_z (const Real end_z)
{
   this->value_ = end_z;
}

Real&
Z_Layer::get_start_z ()
{
   return value;
}

const Real&
Z_Layer::get_start_z () const
{
   return value;
}

Real&
Z_Layer::get_end_z ()
{
   return value_;
}

const Real&
Z_Layer::get_end_z () const
{
   return value_;
}

Real
Z_Layer::get_span_z () const
{
   return fabs (value_ - value);
}

P_Layer::P_Layer (const P_Layer& p_layer)
   : Level (PRESSURE, p_layer.get_start_p (), p_layer.get_end_p ())
{
   order ();
}

P_Layer::P_Layer (const Real start_p,
                  const Real end_p)
   : Level (PRESSURE, start_p, end_p)
{
   order ();
}

void
P_Layer::set (const Real start_p,
              const Real end_p)
{
   this->value = start_p;
   this->value_ = end_p;
   order ();
}

void
P_Layer::set_start_p (const Real start_p)
{
   this->value = start_p;
}

void
P_Layer::set_end_p (const Real end_p)
{
   this->value_ = end_p;
}

Real&
P_Layer::get_start_p ()
{
   return value;
}

const Real&
P_Layer::get_start_p () const
{
   return value;
}

Real&
P_Layer::get_end_p ()
{
   return value_;
}

const Real&
P_Layer::get_end_p () const
{
   return value_;
}

Real
P_Layer::get_span_p () const
{
   return fabs (value_ - value);
}

Real
P_Layer::get_middle_p () const
{
   return (value + value_) * 0.5;
}

Tuple
P_Layer::get_tuple_p_specify_p (const Real specified_p,
                                const Real approx_delta_p) const
{

   const Real& start_p = get_start_p ();
   const Real& end_p = get_end_p ();

   const bool out_of_bounds = (specified_p < start_p) || (specified_p > end_p);
   if (out_of_bounds) { return get_tuple_p (approx_delta_p); }

   const Real upper_span_p = fabs (specified_p - start_p);
   const Real lower_span_p = fabs (end_p - specified_p);
   const Integer upper_n = Integer (round (upper_span_p) / approx_delta_p) + 2;
   const Integer lower_n = Integer (round (lower_span_p) / approx_delta_p) + 2;
   const Real upper_delta_p = upper_span_p / (upper_n - 1);
   const Real lower_delta_p = lower_span_p / (lower_n - 1);

   Tuple tuple_p;

   for (Integer i = 0; i < upper_n - 1; i++)
   {
      const Real p = start_p + i * upper_delta_p;
      tuple_p.push_back (p);
   }

   for (Integer i = 0; i < lower_n - 1; i++)
   {
      const Real p = specified_p + i * lower_delta_p;
      tuple_p.push_back (p);
   }

   tuple_p.push_back (end_p);
   return tuple_p;

}

Tuple
P_Layer::get_tuple_p (const Real approx_delta_p) const
{

   const Real& start_p = get_start_p ();
   const Real& end_p = get_end_p ();
   const Real span_p = end_p - start_p;
   const Integer n = Integer (round (fabs (span_p) / approx_delta_p)) + 2;
   const Real delta_p = span_p / (n - 1);

   Tuple tuple_p;

   for (Integer i = 0; i < n-1; i++)
   {
      const Real p = start_p + i * delta_p;
      tuple_p.push_back (p);
   }

   tuple_p.push_back (end_p);

   return tuple_p;

}

