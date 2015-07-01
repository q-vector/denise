//
// geodesy.cc
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

#include "geodesy.h"

using namespace std;
using namespace denise;

Lat_Long::Lat_Long (const Real latitude,
                    const Real longitude)
{
   this->latitude = latitude;
   this->longitude = longitude;
}

Lat_Long::Lat_Long (const string& lat_long_string)
{

   char* end;
   const Reg_Exp south ("[Ss]$");
   const Reg_Exp west ("[Ww]$");

   const Tokens tokens (lat_long_string, ",");
   const string& latitude_string = tokens[0];
   const string& longitude_string = tokens[1];

   this->latitude = strtod (latitude_string.c_str (), &end);
   this->longitude = strtod (longitude_string.c_str (), &end);

   if (south.match (latitude_string)) { latitude *= -1; }
   if (west.match (longitude_string)) { longitude *= -1; }

}

Lat_Long::Lat_Long (const string& latitude_string,
                    const string& longitude_string)
{

   char* end;
   const Reg_Exp south ("[Ss]$");
   const Reg_Exp west ("[Ww]$");

   this->latitude = strtod (latitude_string.c_str (), &end);
   this->longitude = strtod (longitude_string.c_str (), &end);

   if (south.match (latitude_string)) { latitude *= -1; }
   if (west.match (longitude_string)) { longitude *= -1; }

}

Lat_Long::Lat_Long (const Point_2D& point)
{
   this->latitude = point.x;
   this->longitude = point.y;
}

void
Lat_Long::standardize (Real& latitude,
                       Real& longitude,
                       const Lat_Long_Genre lat_long_genre)
{

   switch (lat_long_genre)
   {

      case LAT_LONG_STANDARD:
      {
         standardize (latitude, longitude, 0);
         break;
      }

      case LAT_LONG_PACIFIC:
      {
         standardize (latitude, longitude, 180);
         break;
      }

   }

}

void
Lat_Long::standardize (Real& latitude,
                       Real& longitude,
                       const Real standard_longitude)
{

   latitude = modulo (latitude, Domain_1D (-180, 180));

   if (fabs (latitude) > 90)
   {
      longitude += 180;
      latitude = (latitude > 0 ? 180 : -180) - latitude;
   }

   const Real start_longitude = standard_longitude - 180;
   const Real end_longitude = standard_longitude + 180;
   const Domain_1D domain_longitude (start_longitude, end_longitude);
   longitude = modulo (longitude, domain_longitude);

}

void
Lat_Long::standardize (const Lat_Long_Genre lat_long_genre)
{
   Lat_Long::standardize (latitude, longitude, lat_long_genre);
}

void
Lat_Long::standardize (const Real standard_longitude)
{
   Lat_Long::standardize (latitude, longitude, standard_longitude);
}

string
Lat_Long::get_string (const Integer decimal_places,
                      const bool plain) const
{
   const string& number_format = string_render ("%%.%df\u00b0", decimal_places);
   return get_string (plain, number_format);
}

string
Lat_Long::get_string (const bool plain,
                      const string& number_format) const
{

   if (is_nall ()) { return ""; }

   Lat_Long ll = *(this);
   ll.standardize (LAT_LONG_STANDARD);
   const Real latitude = ll.latitude;
   const Real longitude = ll.longitude;

   const string& latitude_string = ((latitude >= 0) ?
      string_render ((number_format + "N").c_str (), latitude) :
      string_render ((number_format + "S").c_str (), -latitude));

   const string& longitude_string = ((longitude >= 0) ?
      string_render ((number_format + "E").c_str (), longitude) :
      string_render ((number_format + "W").c_str (), -longitude));

   if (plain) { return latitude_string + longitude_string; }
   else { return "(" + latitude_string + ", " + longitude_string + ")"; }

}

bool
Lat_Long::operator == (const Lat_Long& lat_long) const
{
   return (latitude == lat_long.latitude) && (longitude == lat_long.longitude);
}

bool
Lat_Long::operator != (const Lat_Long& lat_long) const
{
   return (latitude != lat_long.latitude) || (longitude != lat_long.longitude);
}

Lat_Long
Lat_Long::operator = (const Lat_Long& lat_long)
{
   latitude = lat_long.latitude;
   longitude = lat_long.longitude;

   return *this;
}

Lat_Long
Lat_Long::operator + (const Lat_Long& lat_long) const
{
   return Lat_Long (latitude + lat_long.latitude,
                    longitude + lat_long.longitude);
}

void
Lat_Long::operator += (const Lat_Long& lat_long)
{
   latitude += lat_long.latitude;
   longitude += lat_long.longitude;
}

Lat_Long
Lat_Long::operator - ()
{
   return Lat_Long (-latitude, -longitude);
}

Lat_Long
Lat_Long::operator - (const Lat_Long& lat_long) const
{
   return Lat_Long (latitude - lat_long.latitude,
                    longitude - lat_long.longitude);
}

void
Lat_Long::operator -= (const Lat_Long& lat_long)
{
   latitude -= lat_long.latitude;
   longitude -= lat_long.longitude;
}

Lat_Long
Lat_Long::operator * (const Real scalar) const
{
   return Lat_Long (latitude * scalar, longitude * scalar);
}

void
Lat_Long::operator *= (const Real scalar)
{
   latitude *= scalar;
   longitude *= scalar;
}

Lat_Long
Lat_Long::operator / (const Real scalar) const
{
   return Lat_Long (latitude / scalar, longitude / scalar);
}

void
Lat_Long::operator /= (const Real scalar)
{
   latitude /= scalar;
   longitude /= scalar;
}

Lat_Long::operator Point_2D () const
{
   return Point_2D (latitude, longitude);
}

bool
Lat_Long::is_nall () const
{
   return gsl_isnan (latitude) || gsl_isnan (longitude);
}

/*
Location::Location (const string& identifier,
                    const Lat_Long& lat_long)
   : identifier (identifier),
     Lat_Long (lat_long)
{
}

Location::Location (const Lat_Long& lat_long)
   : identifier (""),
     Lat_Long (lat_long)
{
}

Location::Location (const Location& location)
   : identifier (location.identifier),
     Lat_Long (location)
{
}

void
Location::set_identifier (const string& identifier)
{
   this->identifier = identifier;
}

const string&
Location::get_identifier () const
{
   return identifier;
}
 
Location_Set::const_iterator 
Location_Set::nearest (const Lat_Long& lat_long) const
{

   Real min_distance = GSL_POSINF;
   Location_Set::const_iterator nearest = end ();
 
   for (Location_Set::const_iterator iterator = begin (); 
        iterator != end (); iterator++)
   {
 
      const Location& location = *(iterator);
      const Real d = Geodesy::get_distance (lat_long, location);
      
      if (d < min_distance)
      {
         min_distance = d;
         nearest = iterator;
      }

   }

   return nearest;

}
*/

Real
Geodesy::asin (const Real sine)
{
   if (sine > 1) { return M_PI/2; }
   else if (sine < -1) { return -M_PI/2; }
   else { return std::asin (sine); }
}

Real
Geodesy::acos (const Real cosine)
{
   if (cosine > 1) { return 0; }
   else if (cosine < -1) { return M_PI; }
   else { return std::acos (cosine); }
}

void
Geodesy::spherical_inverse (Journey& journey) const
{

   if (journey.get_origin () == journey.get_destination ())
   {
      journey.set_distance_azimuth (0, GSL_NAN, GSL_NAN);
   }

   else
   {

      const Lat_Long& origin = journey.get_origin ();
      const Lat_Long& destination = journey.get_destination ();

      const Real phi_o  = origin.latitude * DEGREE_TO_RADIAN;
      const Real phi_d  = destination.latitude * DEGREE_TO_RADIAN;
      const Real lambda_o = origin.longitude * DEGREE_TO_RADIAN;
      const Real lambda_d = destination.longitude * DEGREE_TO_RADIAN;
      const Real spo = sin (phi_o), cpo = cos (phi_o);
      const Real spd = sin (phi_d), cpd = cos (phi_d);

      const Real d_lambda = lambda_d - lambda_o;
      const Real sdl = sin (d_lambda);
      const Real cdl = cos (d_lambda);

      const Real d = Geodesy::acos (spo*spd + cpo*cpd * cdl);
      const Real sd = sin (d), cd = cos (d);

      Real bf = Geodesy::acos ((spd - spo * cd) / (sd * cpo));
      Real bb = Geodesy::acos ((spo - spd * cd) / (sd * cpd));

      if (sdl < 0) { bf = M_2_TIMES_PI - bf; }
      else         { bb = M_2_TIMES_PI - bb; }


      const Real distance = d * a;
      const Real azimuth_forward = bf * RADIAN_TO_DEGREE;
      const Real azimuth_backward = bb * RADIAN_TO_DEGREE;

      journey.set_distance_azimuth (distance,
         azimuth_forward, azimuth_backward);

   }

}

void
Geodesy::spherical_direct (Journey& journey) const
{

   const Lat_Long& origin = journey.get_origin ();

   const Real phi_o  = origin.latitude * DEGREE_TO_RADIAN;
   const Real lambda_o = origin.longitude * DEGREE_TO_RADIAN;

   const Real d = journey.get_distance () / a;
   const Real bf = journey.get_azimuth_forward () * DEGREE_TO_RADIAN;

   const Real spo = sin (phi_o), cpo = cos (phi_o);
   const Real sd = sin (d), cd = cos (d);

   Real lambda_d;
   const Real phi_d = asin (spo * cd + cpo * sd * cos (bf));

   const Real spd = sin (phi_d), cpd = cos (phi_d);

   if (fabs (cos (phi_d)) < epsilon_v) { lambda_d = 0; }
   else
   {
      const Real d_lambda = atan2 (sin (bf) * sd * cpo, cd - spo * spd);
      lambda_d = modulo (lambda_o + d_lambda + M_PI, M_2_TIMES_PI) - M_PI;
   }

   Real bb = Geodesy::acos ((spo - spd * cd) / (sd * cpd));
   if (sin (lambda_o - lambda_d) < 0) { bb = M_2_TIMES_PI - bb; }

   Lat_Long destination;
   destination.latitude = phi_d * RADIAN_TO_DEGREE;
   destination.longitude = lambda_d * RADIAN_TO_DEGREE;

   journey.set_destination_azimuth (destination, bb * RADIAN_TO_DEGREE);

}

void
Geodesy::vincenty_inverse (Journey& journey) const
{

   const Lat_Long& origin = journey.get_origin ();
   const Lat_Long& destination = journey.get_destination ();

   Real phi_o    = origin.latitude * DEGREE_TO_RADIAN;
   Real lambda_o = origin.longitude * DEGREE_TO_RADIAN;
   Real phi_d    = destination.latitude * DEGREE_TO_RADIAN;
   Real lambda_d = destination.longitude * DEGREE_TO_RADIAN;

   if ((phi_o + phi_d == 0) && (fabs (lambda_o - lambda_d) == M_PI))
   {
      // Antipodal points. Adjust origin slightly.
      phi_o = phi_o + epsilon_v;
   }

   if (phi_o == phi_d && (lambda_o == lambda_d ||
       fabs (fabs (lambda_o - lambda_d) - M_2_TIMES_PI) < epsilon_v))
   {
      // Identical points.
      journey.set_distance_azimuth (0, GSL_NAN, GSL_NAN);
      return;
   }

   const Real r = 1 - f;
   Real tu1 = r * tan (phi_o);
   Real tu2 = r * tan (phi_d);
   const Real cu1 = 1 / sqrt (1 + tu1*tu1);
   const Real su1 = cu1 * tu1;
   const Real cu2 = 1 / sqrt (1 + tu2*tu2);
   const Real s1 = cu1 * cu2;
   const Real b1 = s1 * tu2;
   const Real f1 = b1 * tu1;

   Real x = lambda_d - lambda_o, y;
   Real xx, sx, sy, cx, cy, sa, c2a, cz, e, c, d;

   do
   {

      sx = sin (x);
      cx = cos (x);
      tu1 = cu2*sx;
      tu2 = b1 - su1*cu2*cx;
      sy = sqrt (tu1*tu1 + tu2*tu2);
      cy = s1*cx + f1;
      y = atan2 (sy, cy);
      sa = s1 * sx / sy;
      c2a = 1 - sa * sa;
      cz = f1 + f1;
      if (c2a > 0) { cz = cy - cz / c2a; }
      e = 2 * cz*cz - 1;
      c = ((-3 * c2a + 4) * f + 4) * c2a * f / 16;
      xx = x;
      x = ((e * cy * c + cz) * sy * c + y) * sa;
      x = (1. - c) * x * f + lambda_d - lambda_o;

   }
   while (fabs (xx - x) > epsilon_v);

   const Real theta_f = atan2 (tu1, tu2);
   const Real theta_b = atan2 (cu1*sx, b1*cx - su1*cu2) + M_PI;

   x = sqrt ((1 / (r*r) - 1) * c2a + 1);
   x += 1;
   x = (x - 2) / x;

   c = 1 - x;
   c = (x * x/4 + 1) / c;

   d = (0.375*x*x - 1) * x;
   x = e * cy;

   const Real nautical_miles = ((((4*sy*sy-3) *
      (1-e-e)*cz*d/6 - x) * d/4 + cz) * sy*d + y) * c*a*r;

   const Real distance = nautical_miles;
   const Real azimuth_forward = fmod (theta_f * RADIAN_TO_DEGREE, 360);
   const Real azimuth_backward = fmod (theta_b * RADIAN_TO_DEGREE, 360);
   journey.set_distance_azimuth (distance, azimuth_forward, azimuth_backward);

}
          
void
Geodesy::vincenty_direct (Journey& journey) const
{

   const Lat_Long& origin = journey.get_origin ();

   const Real phi_o    = origin.latitude * M_PI/180;
   const Real lambda_o = origin.longitude * M_PI/180;

   const Real s = journey.get_distance ();
   const Real faz = journey.get_azimuth_forward () * DEGREE_TO_RADIAN;

   const Real r = 1 - f;
   Real tu = r * tan (phi_o);
   const Real sf = sin (faz);
   const Real cf = cos (faz);
   Real b = (cf == 0) ? 0 : 2 * atan2 (tu, cf);
   const Real cu = 1 / sqrt (1 + tu*tu);
   const Real su = tu*cu;
   const Real sa = cu*sf;
   const Real c2a = 1 - sa*sa;

   Real c, d, x, y, sy, cy, cz, e;

   x = 1 + sqrt (1 + c2a * (1/(r*r) - 1));
   x = (x - 2) / x;

   c = 1 - x;
   c = (x*x / 4 + 1) / c;

   d = (0.375*x*x - 1) * x;

   tu = s / (r*a*c);
   y = tu;
   c = y + 1;

   while (fabs (y - c) > epsilon_v)
   {
      sy = sin (y);
      cy = cos (y);
      cz = cos (b + y);
      e = 2 * cz * cz - 1;
      c = y;
      x = e * cy;
      y = e + e - 1;
      y = (((4*sy*sy - 3) *y*cz*d/6 + x) * d/4 - cz) * sy*d + tu;
   }

   b = cu*cy*cf - su*sy;
   c = r * sqrt (sa*sa + b*b);
   d = su*cy + cu*sy*cf;
   const Real phi_d = fmod (atan2 (d, c) - M_PI_2, M_PI) - M_PI_2;

   c = cu*cy - su*sy*cf;
   x = atan2 (sy*sf, c);
   c = ((-3*c2a + 4) * f + 4) *c2a*f/16;
   d = ((e*cy*c + cz) * sy*c + y) * sa;
   Real lambda_d = lambda_o + x - (1 - c) * d*f;
   lambda_d = modulo (lambda_d - M_PI, M_2_TIMES_PI) - M_PI;

   const Real baz = atan2 (sa, b) + M_PI;

   Lat_Long destination;
   destination.latitude = phi_d * RADIAN_TO_DEGREE;
   destination.longitude = lambda_d * RADIAN_TO_DEGREE;

   journey.set_destination_azimuth (destination,
      fmod (baz * RADIAN_TO_DEGREE, 360));

}

void
Geodesy::inverse (Journey& journey) const
{
   if (f == 0) { spherical_inverse (journey); }
   else        { vincenty_inverse (journey); }
}

void
Geodesy::direct (Journey& journey) const
{
   if (f == 0) { spherical_direct (journey); }
   else        { vincenty_direct (journey); }
}

Geodesy::Geodesy (const Geodesy_Model geodesy_model,
                  const Real epsilon_v)
     : epsilon_v (epsilon_v)
{

   switch (geodesy_model)
   {

      case AIRY_1930:             a = 6377563.396; f = 1/299.3249646;   break;
      case MODIFIED_AIRY:         a = 6377340.189; f = 1/299.3249646;   break;
      case AUSTRALIAN_NATIONAL:   a = 6378160.0;   f = 1/298.25;        break;
      case BESSEL_1841:           a = 6377397.155; f = 1/299.1528128;   break;
      case BESSEL_NAMIBIA_1841:   a = 6377483.865; f = 1/299.1528128;   break;
      case CLARKE_1866:           a = 6378206.4;   f = 1/294.9786982;   break;
      case CLARKE_1880:           a = 6378249.145; f = 1/293.465;       break;
      case EVEREST_INDIA_1830:    a = 6377276.345; f = 1/300.8017;      break;
      case EVEREST_INDIA_1956:    a = 6377301.243; f = 1/300.8017;      break;
      case EVEREST_SABAH:         a = 6377298.556; f = 1/300.8017;      break;
      case EVEREST_MALAYSIA_1948: a = 6377304.063; f = 1/300.8017;      break;
      case EVEREST_MALAYSIA_1969: a = 6377295.664; f = 1/300.8017;      break;
      case MODIFIED_FISCHER_1960: a = 6378155.0;   f = 1/298.3;         break;
      case GRS_1967:              a = 6378160.0;   f = 1/298.247167427; break;
      case GRS_1980:              a = 6378137.0;   f = 1/298.257222101; break;
      case HELMERT_1906:          a = 6378200.0;   f = 1/298.3;         break;
      case HOUGH_1960:            a = 6378270.0;   f = 1/297.0;         break;
      case INDONESIAN_1974:       a = 6378160.0;   f = 1/298.247;       break;
      case INTERNATIONAL_1924:    a = 6378388.0;   f = 1/297.0;         break;
      case KRASSOVSKY_1940:       a = 6378245.0;   f = 1/298.3;         break;
      case SOUTH_AMERICAN_1969:   a = 6378160.0;   f = 1/298.25;        break;
      case WGS60:                 a = 6378165.0;   f = 1/298.3;         break;
      case WGS66:                 a = 6378145.0;   f = 1/298.25;        break;
      case WGS72:                 a = 6378135.0;   f = 1/298.26;        break;
      case WGS84:                 a = 6378137.0;   f = 1/298.257223563; break;

      default:
      case SPHERE:                a = EARTH_RADIUS; f = 0; break;

   }

}
          
Geodesy::Geodesy (const Real semi_major_axis,
                  const Real flattening,
                  const Real epsilon_v)
             : a (semi_major_axis),
               f (flattening),
       epsilon_v (epsilon_v)
{
}

Real
Geodesy::get_semi_major_axis () const
{
   return a;
}

Real
Geodesy::get_semi_minor_axis () const
{
   return a * (1 - f);
}

Real
Geodesy::get_flattening () const
{
   return f;
}
          
Real
Geodesy::get_eccentricity () const
{
   return f * (2 - f);
}

Real
Geodesy::get_epsilon_v () const
{
   return epsilon_v;
}

void
Geodesy::complete (Journey& journey) const
{

   if (!journey.get_origin ().is_nall () &&
       !gsl_isnan (journey.get_distance ()) &&
       !gsl_isnan (journey.get_azimuth_forward ()))
   {
      direct (journey);
      return;
   }

   if (!journey.get_origin ().is_nall () &&
       !journey.get_destination ().is_nall ())
   {
      inverse (journey);
      return;
   }

   if (!journey.get_destination ().is_nall () &&
       !gsl_isnan (journey.get_distance ()) &&
       !gsl_isnan (journey.get_azimuth_backward ()))
   {
      journey.swap ();
      direct (journey);
      journey.swap ();
      return;
   }

   if (!journey.get_origin ().is_nall () &&
       (journey.get_distance () < 1))
   {
      journey.set_destination (journey.get_origin ());
      return;
   }

   throw Exception ("Insufficient info to complete journey.");

}

Real
Geodesy::get_angle (const Lat_Long& lat_long_a,
                    const Lat_Long& lat_long_b,
                    const Geodesy& geodesy)
{
   Journey journey (lat_long_a, lat_long_b);
   geodesy.complete (journey);
   return journey.get_distance () / EARTH_RADIUS;
}

Real
Geodesy::get_distance (const Lat_Long& lat_long_a,
                       const Lat_Long& lat_long_b,
                       const Geodesy& geodesy)
{
   Journey journey (lat_long_a, lat_long_b);
   geodesy.complete (journey);
   return journey.get_distance ();
}

Lat_Long
Geodesy::get_destination (const Lat_Long& origin,
                          const Real distance,
                          const Real azimuth_forward,
                          const Geodesy& geodesy)
{
   Journey journey (origin, distance, azimuth_forward);
   geodesy.complete (journey);
   return journey.get_destination ();
}

Real
Geodesy::get_azimuth_forward (const Lat_Long& origin,
                              const Lat_Long& destination,
                              const Geodesy& geodesy)
{
   Journey journey (origin, destination);
   geodesy.complete (journey);
   return journey.get_azimuth_forward ();
}

Real
Geodesy::get_azimuth_backward (const Lat_Long& origin,
                               const Lat_Long& destination,
                               const Geodesy& geodesy)
{
   Journey journey (origin, destination);
   geodesy.complete (journey);
   return journey.get_azimuth_backward ();
}

Journey::Journey (const Lat_Long& origin,
                  const Lat_Long& destination)
   : Edge (origin, destination),
     distance (GSL_NAN),
     azimuth_forward (GSL_NAN),
     azimuth_backward (GSL_NAN)
{
}

Journey::Journey (const Lat_Long& origin,
                  const Lat_Long& destination,
                  const Geodesy& geodesy)
   : Edge (origin, destination),
     distance (GSL_NAN),
     azimuth_forward (GSL_NAN),
     azimuth_backward (GSL_NAN)
{
   geodesy.complete (*this);
}

Journey::Journey (const Lat_Long& origin,
                  const Real distance,
                  const Real azimuth_forward)
   : Edge (origin, Point_2D (GSL_NAN, GSL_NAN)),
     distance (distance),
     azimuth_forward (azimuth_forward),
     azimuth_backward (GSL_NAN)
{
}

Journey::Journey (const Lat_Long& origin,
                  const Real distance,
                  const Real azimuth_forward,
                  const Geodesy& geodesy)
   : Edge (origin, Point_2D (GSL_NAN, GSL_NAN)),
     distance (distance),
     azimuth_forward (azimuth_forward),
     azimuth_backward (GSL_NAN)
{
   geodesy.complete (*this);
}

Journey::Journey (const Journey& journey)
   : Edge (journey),
     distance (journey.distance),
     azimuth_forward (journey.azimuth_forward),
     azimuth_backward (journey.azimuth_backward)
{
}

void
Journey::standardize (const Lat_Long_Genre lat_long_genre)
{

   Lat_Long origin (point_a);
   Lat_Long destination (point_b);

   origin.standardize (lat_long_genre);
   destination.standardize (lat_long_genre);

   point_a = origin;
   point_b = destination;

}

void
Journey::standardize (const Real standard_longitude)
{

   Lat_Long origin (point_a);
   Lat_Long destination (point_b);

   origin.standardize (standard_longitude);
   destination.standardize (standard_longitude);

   point_a = origin;
   point_b = destination;

}

string
Journey::get_cardinal_direction (const Real azimuth)
{

   const Integer a = Integer (round (azimuth / 22.5)) % 16;

   switch (a)
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

   return "";

}

void
Journey::swap ()
{
   Edge::swap ();
   std::swap (azimuth_forward, azimuth_backward);
}

void
Journey::translate (const Real distance,
                    const Real azimuth)
{

   const Geodesy geodesy;
   const Journey journey_a (point_a, distance, azimuth, geodesy);
   const Journey journey_b (point_b, distance, azimuth, geodesy);

   point_a = journey_a.get_destination ();
   point_b = journey_b.get_destination ();

   geodesy.complete (*this);

}

void
Journey::set_origin (const Lat_Long& origin)
{
   this->point_a = origin;
   this->distance = GSL_NAN;
   this->azimuth_forward = GSL_NAN;
   this->azimuth_backward = GSL_NAN;
}

void
Journey::set_destination (const Lat_Long& destination)
{
   this->point_b = destination;
   this->distance = GSL_NAN;
   this->azimuth_forward = GSL_NAN;
   this->azimuth_backward = GSL_NAN;
}

void
Journey::set_destination_azimuth (const Lat_Long& destination,
                                  const Real azimuth_backward)
{
   this->point_b = destination;
   this->azimuth_backward = azimuth_backward;
}

void
Journey::set_distance_azimuth (const Real distance,
                               const Real azimuth_forward,
                               const Real azimuth_backward)
{
   this->distance = distance;
   this->azimuth_forward = azimuth_forward;
   this->azimuth_backward = azimuth_backward;
}

Lat_Long
Journey::get_origin () const
{
   return Lat_Long (point_a);
}

Lat_Long
Journey::get_destination () const
{
   return Lat_Long (point_b);
}

const Real&
Journey::get_distance () const
{
   return distance;
}
   
const Real&
Journey::get_azimuth_forward () const
{
   return azimuth_forward;
}
   
const Real&
Journey::get_azimuth_backward () const
{
   return azimuth_backward;
}

Lat_Long
Journey::get_middle_lat_long (const Geodesy& geodesy)
{
   return get_lat_long (0.5, geodesy);
}

Lat_Long
Journey::get_lat_long (const Real fraction,
                       const Geodesy& geodesy)
{
   geodesy.complete (*this);
   const Lat_Long& origin = get_origin ();
   const Journey j (origin, distance/2, azimuth_forward, geodesy);
   return j.get_destination ();
}

void
Journey::fill_lat_long_list (list<Lat_Long>& lat_long_list,
                             const Real approx_d,
                             const Integer max_n_per_leg,
                             const bool first_leg,
                             const Geodesy& geodesy) const
{

   Integer n = Integer (round (distance / (approx_d))) + 1;
   if (n < 2) { n = 2; } else if (n > max_n_per_leg) { n = max_n_per_leg; }

   const Real d = distance / (n - 1);

   for (Integer i = (first_leg ? 0 : 1); i < n; i++)
   {

      Journey j (get_origin (), i * d, azimuth_forward);
      geodesy.complete (j);

      const Lat_Long& ll = j.get_destination ();
      lat_long_list.push_back (ll);

   }

}

list<Lat_Long>*
Journey::get_lat_long_list_ptr (const Real approx_d,
                                const Integer max_n_per_leg,
                                const Geodesy& geodesy) const
{
   list<Lat_Long>* lat_long_list_ptr = new list<Lat_Long>;
   fill_lat_long_list (*lat_long_list_ptr, approx_d, max_n_per_leg, true, geodesy);
   return lat_long_list_ptr;
}

Multi_Journey::Multi_Journey (const bool closed)
   : Simple_Polyline (closed)
{
}

Multi_Journey::Multi_Journey (const Journey& journey,
                              const Geodesy& geodesy,
                              const Real d_distance,
                              const Integer max_n)
   : Simple_Polyline (false)
{

   Journey j (journey);
   geodesy.complete (j);

   const Lat_Long& origin = j.get_origin ();
   const Real distance = j.get_distance ();
   const Real azimuth_forward = j.get_azimuth_forward ();

   Integer n = Integer (round (distance / (d_distance))) + 1;
   if (n < 2) { n = 2; } else if (n > max_n) { n = max_n; }
   const Real d = distance / (n - 1);

   for (Integer i = 0; i < n; i++)
   {
      Journey jj (origin, i * d, azimuth_forward);
      geodesy.complete (jj);
      push_back (jj.get_destination ());
   }

}

Multi_Journey::Multi_Journey (const Multi_Journey& multi_journey,
                              const Geodesy& geodesy,
                              const Real d_distance,
                              const Integer max_n)
   : Simple_Polyline (multi_journey.closed)
{

   for (Multi_Journey::const_iterator iterator = multi_journey.begin ();
        iterator != multi_journey.end (); iterator++)
   {

      // Skip the node if this is the last node of an unclosed multi_journey
      if (!multi_journey.closed && multi_journey.is_last (iterator))
      {
         continue;
      }

      Journey journey = multi_journey.get_journey (iterator);
      geodesy.complete (journey);

      const Lat_Long& origin = journey.get_origin ();
      const Real distance = journey.get_distance ();
      const Real azimuth_forward = journey.get_azimuth_forward ();

      Integer n = Integer (round (distance / (d_distance))) + 1;
      if (n < 2) { n = 2; } else if (n > max_n) { n = max_n; }
      const Real d = distance / (n - 1);

      for (Integer i = 0; i < n; i++)
      {

         if (i == 0)
         {
            if (iterator == multi_journey.begin ()) { push_back (*(iterator)); }
            continue;
         }

         Journey j (origin, i * d, azimuth_forward);
         geodesy.complete (j);
         push_back (j.get_destination ());

      }

   }

}

Multi_Journey::Multi_Journey (const Simple_Polyline& simple_polyline)
   : Simple_Polyline (simple_polyline)
{
}

void
Multi_Journey::standardize (const Lat_Long_Genre lat_long_genre)
{
   for (Multi_Journey::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Point_2D& point = *(iterator);
      Lat_Long lat_long (point);
      lat_long.standardize (lat_long_genre);
      point = lat_long;
   }
}

void
Multi_Journey::standardize (const Real standard_longitude)
{
   for (Multi_Journey::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Point_2D& point = *(iterator);
      Lat_Long lat_long (point);
      lat_long.standardize (standard_longitude);
      point = lat_long;
   }
}

void
Multi_Journey::translate (const Real distance,
                          const Real azimuth)
{

   const Geodesy geodesy;

   for (Multi_Journey::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Point_2D& point = *(iterator);
      const Journey journey (point, distance, azimuth, geodesy);
      point = journey.get_destination ();
   }

}

Journey
Multi_Journey::get_journey (const Real x,
                            const Geodesy& geodesy) const
{

   if (x < 0) { throw Exception ("Before Lat Point"); }
   Real distance = 0;

   for (Multi_Journey::const_iterator i = begin (); i != end (); i++)
   {
      if (!closed && is_last (i)) { throw Exception ("Beyond Lat Point"); }
      Journey journey = get_journey (i);
      geodesy.complete (journey);
      distance += journey.get_distance ();
      if (x < distance) { return journey; }
   }

}

Journey
Multi_Journey::get_journey (Multi_Journey::iterator iterator) const
{
   if (is_last (iterator))
   {
      if (closed) { return Journey (*(iterator), *(begin ())); }
      else { throw Exception ("Last Node in Multi_Journey"); }
   }
   else
   {
      Multi_Journey::iterator next = iterator;
      return Journey (*(iterator), *(++next));
   }
}

Journey
Multi_Journey::get_journey (Multi_Journey::const_iterator iterator) const
{
   if (is_last (iterator))
   {
      if (closed) { return Journey (*(iterator), *(begin ())); }
      else { throw Exception ("Last Node in Multi_Journey"); }
   }
   else
   {
      Multi_Journey::const_iterator next = iterator;
      return Journey (*(iterator), *(++next));
   }
}

Real
Multi_Journey::get_distance (const Geodesy& geodesy) const
{

   if (size () < 2) { return GSL_NAN; }

   Real distance = 0;

   for (Multi_Journey::const_iterator i = begin (); i != end (); i++)
   {
      if (!closed && is_last (i)) { continue; }
      Journey journey = get_journey (i);
      geodesy.complete (journey);
      distance += journey.get_distance ();
   }

   return distance;

}

Tuple
Multi_Journey::get_tuple_x (const Geodesy& geodesy) const
{

   Tuple tuple_x;
   tuple_x.push_back (0);

   Real distance = 0;

   for (Multi_Journey::const_iterator i = begin (); i != end (); i++)
   {
      if (!closed && is_last (i)) { continue; }
      Journey journey = get_journey (i);
      geodesy.complete (journey);
      distance += journey.get_distance ();
      tuple_x.push_back (distance);
   }

   return tuple_x;

}

Lat_Long
Multi_Journey::get_lat_long (const Real x,
                             const Geodesy& geodesy) const
{

   Real distance = 0;
   if (x < 0) { return Lat_Long (GSL_NAN, GSL_NAN); }

   for (Multi_Journey::const_iterator i = begin (); i != end (); i++)
   {

      if (!closed && is_last (i)) { continue; }
      Journey journey = get_journey (i);
      geodesy.complete (journey);
      const Real d_distance = journey.get_distance ();

      if (x >= distance && x <= distance + d_distance)
      {
         const Real azimuth = journey.get_azimuth_forward ();
         const Lat_Long& origin = journey.get_origin ();
         return geodesy.get_destination (origin, x - distance, azimuth);
      }

      distance += d_distance;

   }

   return Lat_Long (GSL_NAN, GSL_NAN);

}

Real
Multi_Journey::get_azimuth_forward (const Real x,
                                    const Geodesy& geodesy) const
{
   const Lat_Long lat_long = get_lat_long (x, geodesy);
   const Lat_Long& destination = get_journey (x, geodesy).get_destination ();
   return Geodesy::get_azimuth_forward (lat_long, destination, geodesy);
}

Real
Multi_Journey::get_azimuth_forward (Multi_Journey::const_iterator iterator,
                                    const Geodesy& geodesy) const
{

   if (size () < 2) { throw Exception ("Multi_Journey less than 2 nodes"); }

   if (is_last (iterator))
   {
      if (closed)
      {
         const Journey journey (*(iterator), front (), geodesy);
         return journey.get_azimuth_forward ();
      }
      else
      {
         Multi_Journey::const_iterator prev = iterator;
         prev--;
         const Journey journey (*(prev), *(iterator), geodesy);
         return modulo ((journey.get_azimuth_backward () + 180), 360);
      }
   }
   else
   {
      Multi_Journey::const_iterator next = iterator;
      next++;
      const Journey journey (*(iterator), *(next), geodesy);
      return journey.get_azimuth_forward ();
   }

}

Real
Multi_Journey::get_azimuth_forward (Multi_Journey::iterator iterator,
                                    const Geodesy& geodesy) const
{

   if (size () < 2) { throw Exception ("Multi_Journey not long enough"); }

   if (is_last (iterator))
   {
      if (closed)
      {
         const Journey journey (*(iterator), front (), geodesy);
         return journey.get_azimuth_forward ();
      }
      else
      {
         Multi_Journey::iterator prev = iterator;
         prev--;
         const Journey journey (*(prev), *(iterator), geodesy);
         return modulo ((journey.get_azimuth_backward () + 180), 360);
      }
   }
   else
   {
      Multi_Journey::iterator next = iterator;
      next++;
      const Journey journey (*(iterator), *(next), geodesy);
      return journey.get_azimuth_forward ();
   }

}

void
Multi_Journey::cairo (const RefPtr<Context> cr,
                      const Transform_2D& transform) const
{

   const Geodesy geodesy;
   const Real distance = get_distance (geodesy);
   if (distance < 1) { return; }

   Real d_distance = 100e3;
   if (distance < 1000e3) { d_distance = 50e3; }
   if (distance < 400e3) { d_distance = 20e3; }
   if (distance < 200e3) { d_distance = 10e3; }
   Multi_Journey multi_journey (*this, geodesy, d_distance);
   multi_journey.standardize (LAT_LONG_PACIFIC);

   const Real node_size = 16;
   const Color& bg_color = Color::hsb (0.167, 0.2, 0.5, 0.7);
   const Color& fg_color = Color::hsb (0.167, 0.2, 1.0, 1.0);

   // Simple_Polyline is needed because it is "curved"
   Simple_Polyline simple_polyline;
   for (Multi_Journey::const_iterator iterator = multi_journey.begin ();
        iterator != multi_journey.end (); iterator++)
   {
      const Lat_Long& ll = *(iterator);
      const Point_2D& p = transform.transform (ll);
      simple_polyline.push_back (p);
   }

   cr->save ();
   cr->set_line_cap (LINE_CAP_ROUND);
   cr->set_line_join (LINE_JOIN_ROUND);
   cr->set_line_width (node_size);
   bg_color.cairo (cr);
   simple_polyline.cairo (cr);
   cr->stroke ();

   cr->set_line_width (2);
   cr->set_font_size (7);
   fg_color.cairo (cr);

   // draw the nominal numbers
   for (auto iterator = simple_polyline.begin ();
        iterator != simple_polyline.end (); iterator++)
   {
      const Point_2D& p = *(iterator);
      const Integer d = std::distance (simple_polyline.begin (), iterator);
      const string& str = string_render ("%d", d);
      Label (str, p, 'c', 'c').cairo (cr);
   }

   // Two things:
   // - draw dotted circles at each stop
   // - draw the distances of each leg
   for (auto i = begin (); i != end (); i++)
   {

      const Lat_Long& lat_long = *(i);
      const Point_2D& point = transform.transform (lat_long);

      cr->save ();
      Dashes ("1:2").cairo (cr);
      cr->set_line_width (1);
      Ring (node_size / 2).cairo (cr, point);
      fg_color.cairo (cr);
      cr->stroke ();
      cr->restore ();

      if (!closed && is_last (i)) { continue; }

      cr->save ();
      cr->set_font_size (12);

      Journey journey = this->get_journey (i);
      geodesy.complete (journey);

      const Real distance = journey.get_distance ();
      const Lat_Long& origin = journey.get_origin ();
      const Lat_Long& destination = journey.get_destination ();
      const Point_2D& origin_point = transform.transform (origin);
      const Point_2D& destination_point = transform.transform (destination);
      const Real dx = destination_point.x - origin_point.x;
      const Real dy = destination_point.y - origin_point.y;
      const Real theta = atan (dy / dx);

      Lat_Long ll = journey.get_middle_lat_long (geodesy);
      ll.standardize (LAT_LONG_PACIFIC);

      const Point_2D& p = transform.transform (ll);
      const string& str = string_render ("%.0fkm", distance / 1e3);

      Label label (str, p, 'c', 'b', 12);
      label.set_text_angle (theta);
      //label.cairo (cr, fg_color, bg_color, Point_2D (-2, 2));
      label.cairo (cr, bg_color, fg_color, Point_2D (-2, 2));

      cr->restore ();

   }

   cr->restore ();

}

Multi_Journey::iterator
Multi_Journey::get_iterator (const Transform_2D& transform,
                             const Point_2D& point_2d,
                             const Geodesy& geodesy,
                             const Real threshold,
                             const Real standard_longitude,
                             const Real d_distance,
                             const Integer max_n)
{

   Point_2D p_a, p_b;

   for (Multi_Journey::iterator i = begin (); i != end (); i++)
   {

      if (!closed && is_last (i)) { continue; }

      Journey journey = get_journey (i);
      Multi_Journey mj (journey, geodesy, d_distance, max_n);

      for (Multi_Journey::const_iterator j = mj.begin (); j != mj.end (); j++)
      {

         if (!mj.closed && mj.is_last (j)) { continue; }

         // This is original,... but should be wrong!
         //Journey jj = get_journey (j);
         Journey jj = mj.get_journey (j);
         Lat_Long ll_a = jj.get_origin ();
         Lat_Long ll_b = jj.get_destination ();

         ll_a.standardize (standard_longitude);
         ll_b.standardize (standard_longitude);

         transform.transform (p_a.x, p_a.y, ll_a.latitude, ll_a.longitude);
         transform.transform (p_b.x, p_b.y, ll_b.latitude, ll_b.longitude);

         if (Edge::distance (p_a, p_b, point_2d) < threshold)
         {
            return i;
         }

      }

   }

   return end ();

}

Multi_Journey::const_iterator
Multi_Journey::get_iterator (const Transform_2D& transform,
                             const Point_2D& point_2d,
                             const Geodesy& geodesy,
                             const Real threshold,
                             const Real standard_longitude,
                             const Real d_distance,
                             const Integer max_n) const
{

   Point_2D p_a, p_b;

   for (Multi_Journey::const_iterator i = begin (); i != end (); i++)
   {

      if (!closed && is_last (i)) { continue; }

      Journey journey = get_journey (i);
      Multi_Journey mj (journey, geodesy, d_distance, max_n);

      for (Multi_Journey::const_iterator j = mj.begin (); j != mj.end (); j++)
      {

         if (!mj.closed && mj.is_last (j)) { continue; }

         // This is original,... but should be wrong!
         //Journey jj = get_journey (j);
         Journey jj = mj.get_journey (j);
         Lat_Long ll_a = jj.get_origin ();
         Lat_Long ll_b = jj.get_destination ();

         ll_a.standardize (standard_longitude);
         ll_b.standardize (standard_longitude);

         transform.transform (p_a.x, p_a.y, ll_a.latitude, ll_a.longitude);
         transform.transform (p_b.x, p_b.y, ll_b.latitude, ll_b.longitude);

         if (Edge::distance (p_a, p_b, point_2d) < threshold)
         {
            return i;
         }

      }

   }

   return end ();

}

Multi_Journey::iterator
Multi_Journey::implant (const Transform_2D& transform,
                        const Point_2D& point_2d,
                        const Geodesy& geodesy, 
                        const Real threshold,
                        const Real standard_longitude,
                        const Real d_distance, 
                        const Integer max_n)
{

   // if point is not near existing Multi_Journey
   Multi_Journey::iterator iterator = get_iterator (transform,
      point_2d, geodesy, threshold, standard_longitude, d_distance, max_n);
   if (iterator == end ()) { return iterator; }

   const Lat_Long lat_long = transform.reverse (point_2d);
   return insert (++iterator, lat_long);

}

Journey_List::Journey_List (const Lat_Long& lat_long)
   : geodesy (SPHERE)
{
   Journey journey (lat_long, lat_long);
   geodesy.complete (journey);
   push_back (journey);
}

Journey_List::Journey_List (const Lat_Long& lat_long_f,
                            const Lat_Long& lat_long_t)
   : geodesy (SPHERE)
{
   Journey journey (lat_long_f, lat_long_t);
   geodesy.complete (journey);
   push_back (journey);
}

void
Journey_List::add_via_lat_long (Journey_List::iterator iterator,
                                const Lat_Long& lat_long)
{

   Journey& journey = *(iterator);

   Journey j (lat_long, journey.get_destination ());
   journey.set_destination (lat_long);

   geodesy.complete (journey);
   geodesy.complete (j);

   insert (iterator, j);

}

void
Journey_List::set_origin (Journey_List::iterator iterator,
                          const Lat_Long& lat_long)
{
   Journey& journey = *(iterator);
   journey.set_origin (lat_long);
   geodesy.complete (journey);
}

void
Journey_List::set_destination (Journey_List::iterator iterator,
                               const Lat_Long& lat_long)
{
   Journey& journey = *(iterator);
   journey.set_destination (lat_long);
   geodesy.complete (journey);
}

list<Lat_Long>*
Journey_List::get_lat_long_list_ptr (const Real approx_d,
                                     const Integer max_n_per_leg) const
{

   list<Lat_Long>* lat_long_list_ptr = new list<Lat_Long>;

   for (Journey_List::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Journey& journey = *(iterator);
      journey.fill_lat_long_list (*lat_long_list_ptr, approx_d,
         max_n_per_leg, (iterator == begin ()), geodesy);
   }

   return lat_long_list_ptr;

}

pair<string, Lat_Long>
Geodetic_Attractor::nearest (const Lat_Long& lat_long) const
{
   return make_pair ("", lat_long);
}

Degree_Geodetic_Attractor::Degree_Geodetic_Attractor (const Integer n)
   : n (n)
{
}

void
Degree_Geodetic_Attractor::set_n (const Integer n)
{
   this->n = n;
}

const Integer
Degree_Geodetic_Attractor::get_n ()
{
   return n;
}

pair<string, Lat_Long>
Degree_Geodetic_Attractor::nearest (const Lat_Long& lat_long) const
{
   const Real latitude = round (lat_long.latitude * n) / n;
   const Real longitude = round (lat_long.longitude * n) / n;
   return make_pair ("", Lat_Long (latitude, longitude));
}

Range_Circle::Range_Circle (const Lat_Long& lat_long,
                            const Real range,
                            const Real standard_longitude,
                            const Integer n,
                            const Geodesy_Model geodesy_model,
                            const Real epsilon_v)
{

   Lat_Long last_ll;
   const Real da = 360.0 / n;

   const Real start_longitude = standard_longitude - 180;
   const Real end_longitude = standard_longitude + 180;

   for (Integer i = 0; i < n; i++)
   {

      const Real azimuth = 90 + i * da;
      Lat_Long ll = Geodesy::get_destination (lat_long, range, azimuth);
      ll.standardize (standard_longitude);

      if (ll.longitude - last_ll.longitude > 350)
      {
         add (Point_2D (ll.latitude, start_longitude));
         add (Point_2D (90, start_longitude));
         add (Point_2D (90, end_longitude));
         add (Point_2D (last_ll.latitude, end_longitude));
      }
      else
      if (last_ll.longitude - ll.longitude > 350)
      {
         add (Point_2D (last_ll.latitude, end_longitude));
         add (Point_2D (-90, end_longitude));
         add (Point_2D (-90, start_longitude));
         add (Point_2D (ll.latitude, start_longitude));
      }

      add (Point_2D (ll.latitude, ll.longitude));
      last_ll.latitude = ll.latitude;
      last_ll.longitude = ll.longitude;
   }

}

Range_Sector::Range_Sector (const Lat_Long& lat_long,
                            const Real range,
                            const Real start_azimuth,
                            const Real end_azimuth,
                            const Real standard_longitude,
                            const Integer n,
                            const Geodesy_Model geodesy_model,
                            const Real epsilon_v)
{

   Lat_Long last_ll;
   const Real azimuth_span = end_azimuth - start_azimuth;
   const Integer nn = Integer (round (azimuth_span / 360.0 * n));
   const Real da = azimuth_span / nn;

   const Real start_longitude = standard_longitude - 180;
   const Real end_longitude = standard_longitude + 180;

   add (Point_2D (lat_long.latitude, lat_long.longitude));

   for (Integer i = 0; i < nn; i++)
   {

      const Real azimuth = start_azimuth + i * da;
      Lat_Long ll = Geodesy::get_destination (lat_long, range, azimuth);
      ll.standardize (standard_longitude);

      if (ll.longitude - last_ll.longitude > 350)
      {
         add (Point_2D (ll.latitude, start_longitude));
         add (Point_2D (90, start_longitude));
         add (Point_2D (90, end_longitude));
         add (Point_2D (last_ll.latitude, end_longitude));
      }
      else
      if (last_ll.longitude - ll.longitude > 350)
      {
         add (Point_2D (last_ll.latitude, end_longitude));
         add (Point_2D (-90, end_longitude));
         add (Point_2D (-90, start_longitude));
         add (Point_2D (ll.latitude, start_longitude));
      }

      add (Point_2D (ll.latitude, ll.longitude));
      last_ll.latitude = ll.latitude;
      last_ll.longitude = ll.longitude;

   }

}

Geodetic_Transform::Data::Data (const Genre genre,
                                const Real scale,
                                const Lat_Long& lat_long)
   : genre (genre),
     scale (scale),
     lat_long (lat_long)
{
}

Geodetic_Transform::Data::Data (const string& str)
{

   const Tokens tokens (str, ":");

   if (tokens[0] == "MERCATOR")
   {
      genre = MERCATOR;
   }
   else
   if (tokens[0] == "LAMBERT_CONIC_NORTH")
   {
      genre = LAMBERT_CONIC_NORTH;
   }
   else
   if (tokens[0] == "LAMBERT_CONIC_SOUTH")
   {
      genre = LAMBERT_CONIC_SOUTH;
   }
   else
   if (tokens[0] == "POLAR_STEREOGRAPHIC_NORTH")
   {
      genre = POLAR_STEREOGRAPHIC_NORTH;
   }
   else
   if (tokens[0] == "POLAR_STEREOGRAPHIC_SOUTH")
   {
      genre = POLAR_STEREOGRAPHIC_SOUTH;
   }
   else
   if (tokens[0] == "GEOS")
   {
      genre = GEOS;
   }
   else
   if (tokens[0] == "MOLLWEIDE")
   {
      genre = MOLLWEIDE;
   }
   else
   {
      genre = MERCATOR;
   }

   scale = atof (tokens[1].c_str ());
   lat_long.latitude = atof (tokens[2].c_str ());
   lat_long.longitude = atof (tokens[3].c_str ());

}

Geodetic_Transform::Data::Data (const Data& data)
   : genre (data.genre),
     scale (data.scale),
     lat_long (data.lat_long)
{
}

void
Geodetic_Transform::Data::standardize (Lat_Long& lat_long) const
{
   const Real standard_longitude = this->lat_long.longitude;
   lat_long.standardize (standard_longitude);
}

string
Geodetic_Transform::Data::get_string () const
{

   string str;

   switch (genre)
   {
      case MERCATOR:                  str = "MERCATOR"; break;
      case LAMBERT_CONIC_NORTH:       str = "LAMBERT_CONIC_NORTH"; break;
      case LAMBERT_CONIC_SOUTH:       str = "LAMBERT_CONIC_SOUTH"; break;
      case POLAR_STEREOGRAPHIC_NORTH: str = "POLAR_STEREOGRAPHIC_NORTH"; break;
      case POLAR_STEREOGRAPHIC_SOUTH: str = "POLAR_STEREOGRAPHIC_SOUTH"; break;
      case GEOS:                      str = "GOES"; break;
      case MOLLWEIDE:                 str = "MOLLWEIDE"; break;
   }

   const Real latitude = lat_long.latitude;
   const Real longitude = lat_long.longitude;
   str += string_render (":%.0f:%.2f:%.2f", scale, latitude, longitude);
   return str;

}

Geodetic_Transform::Geodetic_Transform (const Geodetic_Transform::Genre genre,
                                        const Real scale,
                                        const Lat_Long& lat_long)
   : data (genre, scale, lat_long)
{
}

Geodetic_Transform::Geodetic_Transform (const string& str)
   : data (str)
{
}

Geodetic_Transform::Geodetic_Transform (const Geodetic_Transform& gt)
   : data (gt.data)
{
}

bool
Geodetic_Transform::is_geodetic (const string& str)
{
   const Tokens tokens (str, ":");
   return ((tokens[0] == "MERCATOR") ||
           (tokens[0] == "LAMBERT_CONIC_NORTH") ||
           (tokens[0] == "LAMBERT_CONIC_SOUTH") ||
           (tokens[0] == "POLAR_STEREOGRAPHIC_NORTH") ||
           (tokens[0] == "POLAR_STEREOGRAPHIC_SOUTH") ||
           (tokens[0] == "GEOS") ||
           (tokens[0] == "MOLLWEIDE"));
}

Geodetic_Transform*
Geodetic_Transform::get_transform_ptr (const string& str,
                                       const Point_2D& point)
{
   const Geodetic_Transform::Data data (str);
   return get_transform_ptr (data.genre, data.scale, data.lat_long, point);
}

Geodetic_Transform*
Geodetic_Transform::get_transform_ptr (const Geodetic_Transform::Genre genre,
                                       const Real scale,
                                       const Lat_Long& lat_long,
                                       const Point_2D& point)
{

   typedef Polar_Stereographic_Transform Pst;

   switch (genre)
   {

      default:
      case MERCATOR:
         return new Mercator_Transform (scale, lat_long, point);

      case LAMBERT_CONIC_NORTH:
         return new Lambert_Conic_Transform (scale, lat_long, point, false);

      case LAMBERT_CONIC_SOUTH:
         return new Lambert_Conic_Transform (scale, lat_long, point, true);

      case POLAR_STEREOGRAPHIC_NORTH:
         return new Pst (scale, lat_long.longitude, point, false);

      case POLAR_STEREOGRAPHIC_SOUTH:
         return new Pst (scale, lat_long.longitude, point, true);

      case MOLLWEIDE:
         return new Mollweide_Transform (scale, point, lat_long);

   }

}

Geodetic_Transform*
Geodetic_Transform::get_transform_ptr (const Geodetic_Transform::Data& data,
                                       const Point_2D& point)
{
   return get_transform_ptr (data.genre, data.scale, data.lat_long, point);
}

bool
Geodetic_Transform::is_out_of_domain (const Lat_Long& lat_long) const
{
   const Real latitude = lat_long.latitude;
   const Real longitude = lat_long.longitude;
   return out_of_domain (latitude, longitude);
}

Real
Geodetic_Transform::get_scale () const
{
   return data.scale;
}

void
Geodetic_Transform::reverse (Lat_Long& lat_long,
                             const Point_2D& point) const
{
   Real& latitude = lat_long.latitude;
   Real& longitude = lat_long.longitude;
   reverse (latitude, longitude, point.x, point.y);
   standardize (lat_long);
}

const Lat_Long&
Geodetic_Transform::get_lat_long () const
{
   return data.lat_long;
}

Lat_Long&
Geodetic_Transform::get_lat_long ()
{
   return data.lat_long;
}

Lat_Long
Geodetic_Transform::get_lat_long (const Point_2D& point) const
{
   Lat_Long lat_long;
   reverse (lat_long.latitude, lat_long.longitude, point.x, point.y);
   standardize (lat_long);
   return lat_long;
}

void
Geodetic_Transform::standardize (Lat_Long& lat_long) const
{
   data.standardize (lat_long);
}

Equidistant_Cylindrical_Transform::Equidistant_Cylindrical_Transform (const Domain_1D& domain_latitude,
                                                                      const Domain_1D& domain_longitude,
                                                                      const Real width,
                                                                      const Real height,
                                                                      const Point_2D& anchor)
   : Geodetic_Transform (EQUIDISTANT_CYLINDRICAL, GSL_NAN, Lat_Long (GSL_NAN, GSL_NAN))
{

   const Real scale_x = width / domain_longitude.get_span ();
   const Real scale_y = height / domain_latitude.get_span ();

   affine_transform.translate (anchor.x, anchor.y);
   affine_transform.scale (scale_x, scale_y);
   affine_transform.translate (-domain_longitude.start, -domain_latitude.start);
   affine_transform.rotate (-M_PI_2);
   //affine_transform.scale (1, 1); 
   ////affine_transform.rotate (M_PI_2);
   ////affine_transform.scale (1, -1); 

}

Geodetic_Transform*
Equidistant_Cylindrical_Transform::clone () const
{
   return new Equidistant_Cylindrical_Transform (*this);
}

bool
Equidistant_Cylindrical_Transform::out_of_domain (const Real latitude,
                                                  const Real longitude) const
{
   return false;
}

void
Equidistant_Cylindrical_Transform::transform (Real& x,
                                              Real& y,
                                              const Real latitude,
                                              const Real longitude) const
{
   affine_transform.transform (x, y, latitude, longitude);
}

void
Equidistant_Cylindrical_Transform::reverse (Real& latitude,
                                            Real& longitude,
                                            const Real x,
                                            const Real y) const
{
   affine_transform.reverse (latitude, longitude, x, y);
   Lat_Long::standardize (latitude, longitude, data.lat_long.longitude);
}

Real
Lambert_Conic_Transform::get_cone_constant (const Real latitude_1,
                                            const Real latitude_2) const
{

   const Real colatitude_1 = (90 - latitude_1) * DEGREE_TO_RADIAN;
   const Real colatitude_2 = (90 - latitude_2) * DEGREE_TO_RADIAN;

   return ((log (sin (colatitude_1)) - log (sin (colatitude_2))) /
      (log (tan (colatitude_1 / 2)) - log (tan (colatitude_2 / 2))));

}

Real
Lambert_Conic_Transform::get_r (const Real latitude) const
{
   const Real colatitude = 90 - latitude;
   const Real tan_half_colatitude = tan (colatitude / 2 * DEGREE_TO_RADIAN);
   const Real p = pow (tan_half_colatitude / tan_half_colatitude_1, cone_constant);
   return (EARTH_RADIUS / data.scale) / cone_constant * sin_colatitude_1 * p;
}

Real
Lambert_Conic_Transform::get_latitude (const Real r) const
{
   const Real scale = data.scale;
   const Real p = r * cone_constant * (scale / EARTH_RADIUS) / sin_colatitude_1;
   const Real q = pow (p, 1 / cone_constant);
   return 90 - (2 * atan (q * tan_half_colatitude_1) * RADIAN_TO_DEGREE);
}

Lambert_Conic_Transform::Lambert_Conic_Transform (const Real scale,
                                                  const Lat_Long& ref_lat_long,
                                                  const Point_2D& ref_point,
                                                  const bool southern_hemisphere)
   : Geodetic_Transform (southern_hemisphere ? LAMBERT_CONIC_SOUTH : LAMBERT_CONIC_NORTH, scale, ref_lat_long),
     ref_lat_long (ref_lat_long),
     ref_point (ref_point)
{

   const Real true_latitude_1 = (southern_hemisphere ? -30 : 30);
   const Real true_latitude_2 = (southern_hemisphere ? -60 : 60);

   cone_constant = get_cone_constant (true_latitude_1, true_latitude_2);

   const Real true_colatitude_1 = 90 - true_latitude_1;
   sin_colatitude_1 = sin (true_colatitude_1 * DEGREE_TO_RADIAN);
   tan_half_colatitude_1 = tan (true_colatitude_1 * M_PI/360);

   // Calculation of r_ref_lat must be done after initialization
   //   of cone_constant of sin_colatitude_1 and tan_half_colatitude
   //   because get_r () needs these variables
   r_ref_lat = get_r (ref_lat_long.latitude);

}

Lambert_Conic_Transform::Lambert_Conic_Transform (const Real scale,
                                                  const Lat_Long& ref_lat_long,
                                                  const Point_2D& ref_point,
                                                  const Real true_latitude_1,
                                                  const Real true_latitude_2)
   : Geodetic_Transform (true_latitude_1 < 0 ? LAMBERT_CONIC_SOUTH : LAMBERT_CONIC_NORTH, scale, ref_lat_long),
     ref_lat_long (ref_lat_long),
     ref_point (ref_point)
{

   cone_constant = get_cone_constant (true_latitude_1, true_latitude_2);

   const Real true_colatitude_1 = 90 - true_latitude_1;
   sin_colatitude_1 = sin (true_colatitude_1 * DEGREE_TO_RADIAN);
   tan_half_colatitude_1 = tan (true_colatitude_1 * M_PI/360);

   // Calculation of r_ref_lat must be done after initialization
   //   of cone_constant of sin_colatitude_1 and tan_half_colatitude
   //   because get_r () needs these variables
   r_ref_lat = get_r (ref_lat_long.latitude);

}

Lambert_Conic_Transform::Lambert_Conic_Transform (const Lambert_Conic_Transform& transform)
   : Geodetic_Transform (transform.data.genre, transform.data.scale, transform.data.lat_long),
     ref_point (transform.ref_point),
     ref_lat_long (transform.ref_lat_long),
     cone_constant (transform.cone_constant),
     sin_colatitude_1 (transform.sin_colatitude_1),
     tan_half_colatitude_1 (transform.tan_half_colatitude_1),
     r_ref_lat (transform.r_ref_lat)
{
}

Geodetic_Transform*
Lambert_Conic_Transform::clone () const
{
   return new Lambert_Conic_Transform (*this);
}

bool
Lambert_Conic_Transform::out_of_domain (const Real latitude,
                                        const Real longitude) const
{
   return false;
}

void
Lambert_Conic_Transform::transform (Real& x,
                                    Real& y,
                                    const Real latitude,
                                    const Real longitude) const
{

   const Real& n = cone_constant;
   const Real r = get_r (latitude);
   const Real nl = (longitude - ref_lat_long.longitude) * n * DEGREE_TO_RADIAN;

   const Real dx = r * sin (nl);
   const Real dy = -r * cos (nl);

   x = ref_point.x + dx;
   y = ref_point.y - (r_ref_lat + dy);
   //y = ref_point.y + r_ref_lat + dy;

}

void
Lambert_Conic_Transform::reverse (Real& latitude,
                                  Real& longitude,
                                  const Real x,
                                  const Real y) const
{

   const Real dx = x - ref_point.x;
   const Real dy = (ref_point.y - y) - r_ref_lat;
   //Real dy = y - (ref_point.y + r_ref_lat);

   Real r = sqrt (dx * dx + dy * dy);

   //// cone_constant negative : southern hemisphere
   //if (cone_constant < 0)
   //{
   //   r *= -1; nl += 180;
   //   const Real domain_start = -180 - ref_lat_long.longitude;
   //   nl = modulo (nl, Domain_1D (domain_start, domain_start+360));
   //}

   if (cone_constant < 0)
   {
      const Real nl = (atan2 (-dx, dy) * RADIAN_TO_DEGREE);
      longitude = nl / cone_constant + ref_lat_long.longitude;
      latitude = get_latitude (-r);
   }
   else
   {
      const Real nl = (atan2 (dx, -dy) * RADIAN_TO_DEGREE);
      longitude = nl / cone_constant + ref_lat_long.longitude;
      latitude = get_latitude (r);
   }

   Lat_Long::standardize (latitude, longitude, data.lat_long.longitude);

}

void
Lambert_Conic_Transform::transform_uv (Real& u,
                                       Real& v,
                                       const Real latitude,
                                       const Real longitude) const
{
   const Real theta = cone_constant *
      (longitude - ref_lat_long.longitude) * DEGREE_TO_RADIAN;
   const Real c = cos (theta);
   const Real s = sin (theta);
   const Real uu = u * c - v * s;
   const Real vv = u * s + v * c;
   u = uu;
   v = vv;
}

Real
Lambert_Conic_Transform::get_theta (const Real u,
                                    const Real v,
                                    const Real latitude,
                                    const Real longitude) const
{
   const Real atan2_u_v = (u == 0 && v == 0) ? 0 : -atan2 (v, u);
   return atan2_u_v - cone_constant *
          (longitude - ref_lat_long.longitude) * DEGREE_TO_RADIAN;
}

Real
Mercator_Transform::project_x (const Real longitude) const
{
   return a_e * (longitude - ref_lat_long.longitude) *
          DEGREE_TO_RADIAN * cos_true_latitude;
}

Real
Mercator_Transform::project_y (const Real latitude) const
{
   return a_e * log (tan ((45 + latitude/2) * DEGREE_TO_RADIAN)) *
          cos_true_latitude;
}

Real
Mercator_Transform::reverse_latitude (const Real y) const
{
   return ((atan (exp (y / (a_e * cos_true_latitude))) *
           RADIAN_TO_DEGREE) - 45) * 2;
}

Real
Mercator_Transform::reverse_longitude (const Real x) const
{
   return x / (a_e * DEGREE_TO_RADIAN * cos_true_latitude) +
          ref_lat_long.longitude;
}

Mercator_Transform::Mercator_Transform (const Real scale,
                                        const Lat_Long& ref_lat_long,
                                        const Point_2D& ref_point,
                                        const Real true_latitude)
   : Geodetic_Transform (MERCATOR, scale, ref_lat_long),
     ref_lat_long (ref_lat_long)
{

   a_e = EARTH_RADIUS;
   cos_true_latitude = cos (true_latitude * DEGREE_TO_RADIAN);

   ////offset_x = -ref_point.x;
   ////offset_y = project_y (ref_lat_long.latitude) / scale - ref_point.y;
   offset_x = ref_point.x;
   offset_y = ref_point.y + project_y (ref_lat_long.latitude) / scale;

}

Geodetic_Transform*
Mercator_Transform::clone () const
{
   return new Mercator_Transform (*this);
}

bool
Mercator_Transform::out_of_domain (const Real latitude,
                                   const Real longitude) const
{
   return false;
}

void
Mercator_Transform::transform (Real& x,
                               Real& y,
                               const Real latitude,
                               const Real longitude) const
{
   //yy = project_y (latitude) / data.scale - offset_y;
   //yy = offset_y - project_y (latitude) / data.scale;
   x = offset_x + project_x (longitude) / data.scale;
   y = offset_y - project_y (latitude) / data.scale;
}

void
Mercator_Transform::reverse (Real& latitude,
                             Real& longitude,
                             const Real x,
                             const Real y) const
{
   //latitude = reverse_latitude ((y + offset_y) * data.scale);
   latitude = reverse_latitude ((offset_y - y) * data.scale);
   longitude = reverse_longitude ((x - offset_x) * data.scale);
   Lat_Long::standardize (latitude, longitude, data.lat_long.longitude);
}

Real
Mercator_Transform::get_theta (const Real u,
                               const Real v,
                               const Real latitude,
                               const Real longitude) const
{
   return -atan2 (v, u);
}

Real
Polar_Stereographic_Transform::get_r (const Real latitude) const
{
   const Real phi = latitude * DEGREE_TO_RADIAN;
   const Real q = cos (phi) / (1 + (southern_hemisphere ? -1 : +1) * sin (phi));
   return (EARTH_RADIUS / data.scale) * q * 2;
}
         
Real
Polar_Stereographic_Transform::get_latitude (const Real r) const
{
   Real q = r * (data.scale / EARTH_RADIUS) / 2;
   if (southern_hemisphere) { q = 1 / q; }
   return ((M_PI/2) - 2 * atan (q)) * RADIAN_TO_DEGREE;
}
         
Polar_Stereographic_Transform::Polar_Stereographic_Transform (const Real scale,
                                                              const Real ref_longitude,
                                                              const Point_2D& ref_point,
                                                              const bool southern_hemisphere)
   : Geodetic_Transform (southern_hemisphere ? POLAR_STEREOGRAPHIC_SOUTH : POLAR_STEREOGRAPHIC_NORTH, scale, Lat_Long (0, ref_longitude)),
     ref_point (ref_point),
     ref_longitude (ref_longitude),
     southern_hemisphere (southern_hemisphere)
{
}
                            
Geodetic_Transform*
Polar_Stereographic_Transform::clone () const
{
   return new Polar_Stereographic_Transform (*this);
}

bool
Polar_Stereographic_Transform::out_of_domain (const Real latitude,
                                              const Real longitude) const
{
   return false;
}

void
Polar_Stereographic_Transform::transform (Real& x,
                                          Real& y,
                                          const Real latitude,
                                          const Real longitude) const
{

   Real r = get_r (latitude);
   Real nl = (longitude - ref_longitude) * DEGREE_TO_RADIAN;

   if (southern_hemisphere) { r *= -1; nl *= -1; }

   const Real dx = r * sin (nl);
   const Real dy = r * cos (nl);

   x = ref_point.x + dx;
   y = ref_point.y + dy;

}

void
Polar_Stereographic_Transform::reverse (Real& latitude,
                                        Real& longitude,
                                        const Real x,
                                        const Real y) const
{

   const Real dx = x - ref_point.x;
   const Real dy = y - ref_point.y;

   Real r = sqrt (dx * dx + dy * dy);
   Real nl = (atan2 (dx, dy) * RADIAN_TO_DEGREE);

   if (southern_hemisphere) { nl += 180; }

   const Real domain_start = (southern_hemisphere ?
      ref_longitude - 180 : -180 - ref_longitude);

   nl = modulo (nl, Domain_1D (domain_start, domain_start + 360));

   longitude = ref_longitude + (southern_hemisphere ? -1 : 1) * nl;
   latitude = get_latitude (r);

   Lat_Long::standardize (latitude, longitude, data.lat_long.longitude);

}

void
Polar_Stereographic_Transform::transform_uv (Real& u,
                                             Real& v,
                                             const Real latitude,
                                             const Real longitude) const
{
   Real theta = (longitude - ref_longitude) * DEGREE_TO_RADIAN;
   if (southern_hemisphere) { theta *= -1; }
   const Real c = cos (theta);
   const Real s = sin (theta);
   const Real uu = u * c - v * s;
   const Real vv = u * s + v * c;
   u = uu;
   v = vv;
}

Real
Polar_Stereographic_Transform::get_theta (const Real u,
                                          const Real v,
                                          const Real latitude,
                                          const Real longitude) const
{
   const Real atan2_u_v = (u == 0 && v == 0) ? 0 : -atan2 (v, u);
   const Real rel_lambda = (longitude - ref_longitude) * DEGREE_TO_RADIAN;
   return atan2_u_v - rel_lambda * (southern_hemisphere ? -1 : 1);
}

Perspective_Transform::Perspective_Transform (const Real scale,
                                              const Lat_Long& ref_lat_long,
                                              const Point_2D& ref_point,
                                              const Real height)
   : Geodetic_Transform (PERSPECTIVE, scale, ref_lat_long),
     geodesy (SPHERE),
     ref_lat_long (ref_lat_long),
     ref_point (ref_point),
     height (height),
     range_circle (ref_lat_long, EARTH_RADIUS * acos
        (EARTH_RADIUS / (EARTH_RADIUS + height)), 180, 360),
     ttt (Domain_1D (-91, 91), Domain_1D (-5, 365), 1100, 935, Point_2D (0, 0))

{
   const Real a = EARTH_RADIUS;
   const Real max_range = a * acos (a / (a + height));
   domain_polygon_ptr = new Range_Circle (ref_lat_long, max_range, 180, 360);
}

Perspective_Transform::~Perspective_Transform ()
{
   delete domain_polygon_ptr;
}

Geodetic_Transform*
Perspective_Transform::clone () const
{
   throw Exception ("Perspective_Transform::clone () not implemented");
   return new Perspective_Transform (*this);
}

const Range_Circle&
Perspective_Transform::get_domain () const
{
   return range_circle;
}

bool
Perspective_Transform::out_of_domain (const Real latitude,
                                      const Real longitude) const
{
   return !range_circle.contains (Point_2D (latitude, longitude));
}

void
Perspective_Transform::transform (Real& x,
                                  Real& y,
                                  const Real latitude, 
                                  const Real longitude) const
{

//   ttt.transform (x, y, latitude, longitude);

   Journey journey (ref_lat_long, Lat_Long (latitude, longitude));
   geodesy.complete (journey);

   const Real azimuth_forward = journey.get_azimuth_forward ();
   if (!gsl_finite (azimuth_forward))
   {
      x = ref_point.x;
      y = ref_point.y;
      return;
   }

   const Real alpha = journey.get_azimuth_forward () * DEGREE_TO_RADIAN;
   const Real chi = journey.get_distance () / EARTH_RADIUS;
   const Real tan_eta = sin (chi) / (1 - cos (chi) + height / EARTH_RADIUS);
   const Real r = height * tan_eta;

   x = ref_point.x + r * sin (alpha) / data.scale;
   y = ref_point.y - r * cos (alpha) / data.scale;

}

void
Perspective_Transform::reverse (Real& latitude,
                                Real& longitude,
                                const Real x,
                                const Real y) const
{

   if (x == ref_point.x && y == ref_point.y)
   {
      latitude = ref_lat_long.latitude;
      longitude = ref_lat_long.longitude;
      return;
   }

   const Real xx = ((x - ref_point.x) * data.scale);
   const Real yy = ((ref_point.y - y) * data.scale);

   const Real r = sqrt (xx * xx + yy * yy);
   const Real alpha = atan2 (xx, yy);

   const Real eta = atan (r / height);
   const Real theta = asin (sin (eta) / EARTH_RADIUS * (EARTH_RADIUS + height));
   const Real chi = theta - eta;

   Journey journey (ref_lat_long, chi * EARTH_RADIUS, alpha * RADIAN_TO_DEGREE);
   geodesy.complete (journey);

   const Lat_Long& destination = journey.get_destination ();
   latitude = destination.latitude;
   longitude = destination.longitude;

   Lat_Long::standardize (latitude, longitude, data.lat_long.longitude);

}

void    
Perspective_Transform::cairo (const RefPtr<Context> cr,
                              const Simple_Polyline& simple_polyline) const
{

   Real x, y;
   Lat_Long lat_long;
   bool first_point = true;

   const Real a = EARTH_RADIUS;
   const Real r_m = a * acos (a / (a + height));

   for (Simple_Polyline::const_iterator iterator = simple_polyline.begin ();
        iterator != simple_polyline.end (); iterator++)
   {

      const Point_2D& point = *(iterator);
      lat_long.latitude = point.x;
      lat_long.longitude = point.y;

      if (Geodesy::get_distance (ref_lat_long, lat_long) <= r_m)
      {

         transform (x, y, lat_long.latitude, lat_long.longitude);

         if (first_point)
         {
            cr->move_to (x, y);
            first_point = false;
         }
         else
         {
            cr->line_to (x, y);
         }

      }
      else
      {
         first_point = true;
      }

   }

}

void
Perspective_Transform::cairo (const RefPtr<Context> cr,
                              const Polygon& polygon) const
{

   const Integer n = 360;
   const Real a = EARTH_RADIUS;

   Polygon* p_ptr = Polygon::boolean_op (INTERSECTION, range_circle, polygon);
   Transform_2D::cairo (cr, *p_ptr);
   delete p_ptr;

}

Geos_Transform::Geos_Transform (const Real nadir_longitude,
                                const Real coff,
                                const Real loff,
                                const Real cfac,
                                const Real lfac)
   : Geodetic_Transform (GEOS, GSL_NAN, Lat_Long (0, nadir_longitude)),
     nadir_longitude (nadir_longitude),
     coff (coff),
     loff (loff),
     cfac (cfac),
     lfac (lfac)
{
}

Geodetic_Transform*
Geos_Transform::clone () const
{
   throw Exception ("Geos_Transform::clone () not implemented");
   return new Geos_Transform (*this);
}

bool
Geos_Transform::out_of_domain (const Real latitude,
                               const Real longitude) const
{
   return false;
}

void
Geos_Transform::transform (Real& x,
                           Real& y,
                           const Real latitude,
                           const Real longitude) const
{

   const Real lambda = (longitude - nadir_longitude) * DEGREE_TO_RADIAN;
   const Real phi = atan (0.993243 * tan (latitude * DEGREE_TO_RADIAN));

   const Real cos_phi = cos (phi);
   const Real sin_phi = sin (phi);
   const Real cos_lambda = cos (lambda);
   const Real sin_lambda = sin (lambda);

   const Real rl = 6356.5838 / sqrt (1 - 0.00675701 * cos_phi*cos_phi);
   const Real r1 = 42164 - rl * cos_phi * cos_lambda;
   const Real r2 = -rl * cos_phi * sin_lambda;
   const Real r3 = rl * sin_phi;
   const Real rn = sqrt (r1*r1 + r2*r2 + r3*r3);

   x = coff + (atan (-r2 / r1) * RADIAN_TO_DEGREE / 65536 * cfac);
   y = loff + (asin (-r3 / rn) * RADIAN_TO_DEGREE / 65536 * lfac);

}

void
Geos_Transform::reverse (Real& latitude,
                         Real& longitude,
                         const Real x,
                         const Real y) const
{

   const Real xx = (x - coff) / cfac * 65536 * DEGREE_TO_RADIAN;
   const Real yy = (y - loff) / lfac * 65536 * DEGREE_TO_RADIAN;

   const Real cos_x = cos (xx);
   const Real sin_x = sin (xx);
   const Real cos_y = cos (yy);
   const Real sin_y = sin (yy);

   const Real a = 42164 * cos_x * cos_y;
   const Real b = cos_y*cos_y + 1.006803 * sin_y*sin_y;

   const Real sd = sqrt (a*a - b * 1737121856);
   const Real sn = (a - sd) / b;

   const Real s1 = 42164 - sn * cos_x * cos_y;
   const Real s2 = sn * sin_x * cos_y;
   const Real s3 = -sn * sin_y;
   const Real sxy = sqrt (s1*s1 + s2*s2);

   latitude = atan (1.006803 * (s3 / sxy));
   longitude = atan (s2 / s1) * RADIAN_TO_DEGREE + nadir_longitude;

}

Mollweide_Transform::Mollweide_Transform (const Real scale,
                                          const Point_2D& anchor,
                                          const Lat_Long& lat_long)
   : Geodetic_Transform (MOLLWEIDE, scale, lat_long),
     scale (scale),
     anchor (anchor),
     k (M_SQRT2 * EARTH_RADIUS / scale),
     epsilon (1e-4)
{
}

Geodetic_Transform*
Mollweide_Transform::clone () const
{
   throw Exception ("Mollweide_Transform::clone () not implemented");
   return new Mollweide_Transform (*this);
}

bool
Mollweide_Transform::out_of_domain (const Real latitude,
                                    const Real longitude) const
{
   return false;
}

Real
Mollweide_Transform::get_theta (const Real phi) const
{

   if (fabs (fabs (phi) - M_PI_2) < epsilon) { return phi; }

   Real theta = phi;
   Real step = GSL_POSINF;
   const Real pi_sin_phi = M_PI * sin (phi);

   while (fabs (step) > epsilon)
   {
      const Real two_theta = 2 * theta;
      const Real f = two_theta + sin (two_theta) - pi_sin_phi;
      const Real f_prime = 2 * (1 + cos (two_theta));
      step = -f / f_prime;
      theta += step;
   }

   return theta;

}

Real
Mollweide_Transform::get_phi (const Real theta) const
{
   const Real two_theta = 2 * theta;
   return asin ((two_theta + sin (two_theta)) / M_PI);
}

void
Mollweide_Transform::transform (Real& x,
                                Real& y,
                                const Real latitude,
                                const Real longitude) const
{
   const Real theta = get_theta (latitude * DEGREE_TO_RADIAN);
   const Real lambda = longitude * DEGREE_TO_RADIAN;
   x = anchor.x + k * M_2_PI * cos (theta) * lambda;
   y = anchor.y - k * sin (theta);
}

void
Mollweide_Transform::reverse (Real& latitude,
                              Real& longitude,
                              const Real x,
                              const Real y) const
{
   const Real theta = asin ((anchor.y - y) / k);
   const Real lambda = (x - anchor.x) / (k * M_2_PI * cos (theta));
   const Real phi = get_phi (theta);
   latitude = phi * RADIAN_TO_DEGREE;
   longitude = lambda * RADIAN_TO_DEGREE;
}

Real
Geodetic_Vector_Data_2D::get_dmagnitude_dx (const Integer vector_element_u,
                                            const Integer vector_element_v,
                                            const Integer node_x,
                                            const Integer node_y,
                                            const Real magnitude) const
{
   const Real latitude = get_coordinate (0, node_y);
   const Real eval = Vector_Data_2D::get_dmagnitude_dy (
      vector_element_u, vector_element_v, node_x, node_y, magnitude);
   return eval / (LATITUDE_LENGTH * cos (latitude * DEGREE_TO_RADIAN));
}

Real
Geodetic_Vector_Data_2D::get_dmagnitude_dy (const Integer vector_element_u,
                                            const Integer vector_element_v,
                                            const Integer node_x,
                                            const Integer node_y,
                                            const Real magnitude) const
{
   const Real eval = Vector_Data_2D::get_dmagnitude_dx (
      vector_element_u, vector_element_v, node_x, node_y, magnitude);
   return eval / LATITUDE_LENGTH;
}

Real
Geodetic_Vector_Data_2D::get_dmagnitude_dx (const Integer vector_element_u,
                                            const Integer vector_element_v,
                                            const Real x,
                                            const Real y,
                                            const Real magnitude) const
{
   const Real latitude = y;
   const Real eval = Vector_Data_2D::get_dmagnitude_dy (
      vector_element_u, vector_element_v, x, y, magnitude);
   return eval / (LATITUDE_LENGTH * cos (latitude * DEGREE_TO_RADIAN));
}

Real
Geodetic_Vector_Data_2D::get_dmagnitude_dy (const Integer vector_element_u,
                                            const Integer vector_element_v,
                                            const Real x,
                                            const Real y,
                                            const Real magnitude) const
{
   const Real eval = Vector_Data_2D::get_dmagnitude_dx (
      vector_element_u, vector_element_v, x, y, magnitude);
   return eval / LATITUDE_LENGTH;
}

Geodetic_Vector_Data_2D::Geodetic_Vector_Data_2D (const Integer vector_size,
                                                  const Size_2D& size_2d,
                                                  const Domain_2D& domain_2d,
                                                  const bool periodic_longitude)
   : Vector_Data_2D (vector_size, size_2d, domain_2d, false, periodic_longitude)
{
}

Geodetic_Vector_Data_2D::Geodetic_Vector_Data_2D (const Integer vector_size,
                                                  const Tuple tuple_latitude,
                                                  const Tuple tuple_longitude,
                                                  const bool periodic_longitude)
   : Vector_Data_2D (vector_size, tuple_latitude, tuple_longitude,
                     false, periodic_longitude)
{
}

Geodetic_Vector_Data_2D::~Geodetic_Vector_Data_2D ()
{
}

Real
Geodetic_Vector_Data_2D::get_f (const Real latitude)
{
   return 2 * EARTH_ROTATION * sin (latitude * DEGREE_TO_RADIAN);
}

Real
Geodetic_Vector_Data_2D::get_f_hat (const Real latitude)
{
   return 2 * EARTH_ROTATION * cos (latitude * DEGREE_TO_RADIAN);
}

Real
Geodetic_Vector_Data_2D::get_latitude (const Integer i) const
{
   return get_coordinate (0, i);
}

Real
Geodetic_Vector_Data_2D::get_longitude (const Integer j) const
{
   return get_coordinate (1, j);
}

Real
Geodetic_Vector_Data_2D::evaluate (const Integer vector_element,
                                   const Integer i,
                                   const Integer j,
                                   const Evaluate_Op evaluate_op) const
{

   const Integer& ve = vector_element;

   switch (evaluate_op)
   {

      default:
      {
         return Vector_Data_2D::evaluate (ve, i, j, VALUE);
      }

      case DX:
      {
         const Real latitude = get_coordinate (0, i);
         const Real eval = Vector_Data_2D::evaluate (ve, i, j, DY);
         return eval / (LATITUDE_LENGTH * cos (latitude * DEGREE_TO_RADIAN));
      }

      case DY:
      {
         const Real eval = Vector_Data_2D::evaluate (ve, i, j, DX);
         return eval / LATITUDE_LENGTH;
      }

      case DX2:
      {
         const Real latitude = get_coordinate (0, i);
         const Real c = cos (latitude * DEGREE_TO_RADIAN);
         const Real longitude_length = LATITUDE_LENGTH * c;
         const Real eval = Vector_Data_2D::evaluate (ve, i, j, DY2);
         return eval / (longitude_length * longitude_length);
      }

      case DY2:
      {
         const Real eval = Vector_Data_2D::evaluate (ve, i, j, DX2);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH);
      }

      case DXY:
      {
         const Real latitude = get_coordinate (0, i);
         const Real c = cos (latitude * DEGREE_TO_RADIAN);
         const Real eval = Vector_Data_2D::evaluate (ve, i, j, DXY);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH * c);
      }

      case LAPLACIAN:
      {
         const Real latitude = get_coordinate (0, i);
         const Real c = cos (latitude * DEGREE_TO_RADIAN);
         const Real dx2 = Vector_Data_2D::evaluate (ve, i, j, DY2);
         const Real dy2 = Vector_Data_2D::evaluate (ve, i, j, DX2);
         return (dx2 / (c * c) + dy2) / (LATITUDE_LENGTH * LATITUDE_LENGTH);
      }

   }

}

Real
Geodetic_Vector_Data_2D::evaluate (const Integer vector_element,
                                   const Real latitude,
                                   const Real longitude,
                                   const Evaluate_Op evaluate_op) const
{

   const Integer& ve = vector_element;

   switch (evaluate_op)
   {

      default:
      {
         return Vector_Data_2D::evaluate (ve, latitude, longitude, VALUE);
      }

      case DX:
      {
         Real eval = Vector_Data_2D::evaluate (ve, latitude, longitude, DY);
         return eval / (LATITUDE_LENGTH * cos (latitude * DEGREE_TO_RADIAN));
      }

      case DY:
      {
         Real eval = Vector_Data_2D::evaluate (ve, latitude, longitude, DX);
         return eval / LATITUDE_LENGTH;
      }

      case DX2:
      {
         Real c = cos (latitude * DEGREE_TO_RADIAN);
         Real longitude_length = LATITUDE_LENGTH * c;
         Real eval = Vector_Data_2D::evaluate (ve, latitude, longitude, DY2);
         return eval / (longitude_length * longitude_length);
      }

      case DY2:
      {
         Real eval = Vector_Data_2D::evaluate (ve, latitude, longitude, DX2);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH);
      }

      case DXY:
      {
         Real c = cos (latitude * DEGREE_TO_RADIAN);
         Real eval = Vector_Data_2D::evaluate (ve, latitude, longitude, DXY);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH * c);
      }

      case LAPLACIAN:
      {
         Real c = cos (latitude * DEGREE_TO_RADIAN);
         Real dx2 = Vector_Data_2D::evaluate (ve, latitude, longitude, DY2);
         Real dy2 = Vector_Data_2D::evaluate (ve, latitude, longitude, DX2);
         return (dx2 / (c * c) + dy2) / (LATITUDE_LENGTH * LATITUDE_LENGTH);
      }

   }

}

Real
Geodetic_Vector_Data_2D::evaluate (const Integer vector_element,
                                   const Lat_Long& lat_long,
                                   const Evaluate_Op evaluate_op) const
{
   return evaluate (vector_element, lat_long.latitude,
       lat_long.longitude, evaluate_op);
}

Real
Geodetic_Vector_Data_2D::get_shear_vorticity (const Integer vector_element_u,
                                              const Integer vector_element_v,
                                              const Integer i,
                                              const Integer j) const
{

   return GSL_NAN;
/*
   const Tuple& tuple_x = coordinate_tuples[0];
   const Tuple& tuple_y = coordinate_tuples[1];

   const Size_2D& size_2d = get_size_2d ();
   const Real latitude = tuple_x[i];

   const bool i_start = (i == 0);
   const bool i_end = (i == size_2d.i - 1);
   const Integer im1 = (i_start ? 0 : i - 1);
   const Integer ip1 = (i_end ? size_2d.i - 1 : i + 1);
   const Real x_im1 = tuple_x[im1];
   const Real x_ip1 = tuple_x[ip1];
   const Real h_x = x_ip1 - x_im1;

   const bool j_start = (j == 0);
   const bool j_end = (j == size_2d.j - 1);
   const Integer jm1 = (j_start ? 0 : j - 1);
   const Integer jp1 = (j_end ? size_2d.j - 1 : j + 1);
   const Real y_jm1 = tuple_y[jm1];
   const Real y_jp1 = tuple_y[jp1];
   const Real h_y = y_jp1 - y_jm1;

   const Real u = get_datum (vector_element_u, i, j);
   const Real v = get_datum (vector_element_v, i, j);
   const Real speed = sqrt (u*u + v*v);

   const Real u_im1 = get_datum (vector_element_u, im1, j);
   const Real v_im1 = get_datum (vector_element_v, im1, j);
   const Real speed_im1 = sqrt (u_im1*u_im1 + v_im1*v_im1);

   const Real u_ip1 = get_datum (vector_element_u, ip1, j);
   const Real v_ip1 = get_datum (vector_element_v, ip1, j);
   const Real speed_ip1 = sqrt (u_ip1*u_ip1 + v_ip1*v_ip1);

   const Real u_jm1 = get_datum (vector_element_u, i, jm1);
   const Real v_jm1 = get_datum (vector_element_v, i, jm1);
   const Real speed_jm1 = sqrt (u_jm1*u_jm1 + v_jm1*v_jm1);

   const Real u_jp1 = get_datum (vector_element_u, i, jp1);
   const Real v_jp1 = get_datum (vector_element_v, i, jp1);
   const Real speed_jp1 = sqrt (u_jp1*u_jp1 + v_jp1*v_jp1);

   const Real dspeed_dx = (speed_ip1 - speed_im1) / h_x;
   const Real dspeed_dy = (speed_jp1 - speed_jm1) / h_y;
   return (v * dspeed_dx - u * dspeed_dy) / speed;
*/

}

void
Geodetic_Vector_Data_2D::subtract_zonal_mean (const Integer vector_element)
{
   subtract_y_mean (vector_element);
}

void
Geodetic_Vector_Data_2D::subtract_meridional_mean (const Integer vector_element)
{
   subtract_x_mean (vector_element);
}

Geodetic_Vector_Data_3D::Geodetic_Vector_Data_3D (const Integer vector_size,
                                                  const Size_3D& size_3d,
                                                  const Domain_3D& domain_3d,
                                                  const bool periodic_longitude)
   : Vector_Data_3D (vector_size, size_3d, domain_3d,
                     false, false, periodic_longitude)
{
}

Geodetic_Vector_Data_3D::Geodetic_Vector_Data_3D (const Integer vector_size,
                                                  const Tuple tuple_z,
                                                  const Tuple tuple_latitude,
                                                  const Tuple tuple_longitude,
                                                  const bool periodic_longitude)
   : Vector_Data_3D (vector_size, tuple_z, tuple_latitude, tuple_longitude,
                     false, false, periodic_longitude)
{
}

Geodetic_Vector_Data_3D::Geodetic_Vector_Data_3D (const Integer vector_size,
                                                  const Tuple tuple_z,
                                                  const Size_2D& size_2d,
                                                  const Domain_2D& domain_2d,
                                                  const bool periodic_longitude)
   : Vector_Data_3D (vector_size, tuple_z, size_2d, domain_2d,
                     false, false, periodic_longitude)
{
}

Geodetic_Vector_Data_3D::~Geodetic_Vector_Data_3D ()
{
}

Real
Geodetic_Vector_Data_3D::get_f (const Real latitude)
{
   return 2 * EARTH_ROTATION * sin (latitude * DEGREE_TO_RADIAN);
}

Real
Geodetic_Vector_Data_3D::get_f_hat (const Real latitude)
{
   return 2 * EARTH_ROTATION * cos (latitude * DEGREE_TO_RADIAN);
}

Real
Geodetic_Vector_Data_3D::get_latitude (const Integer i) const
{
   return get_coordinate (1, i);
}

Real
Geodetic_Vector_Data_3D::get_longitude (const Integer j) const
{
   return get_coordinate (2, j);
}

Real
Geodetic_Vector_Data_3D::evaluate (const Integer vector_element,
                                   const Real z,
                                   const Integer i,
                                   const Integer j,
                                   const Evaluate_Op evaluate_op) const
{

   const Integer& ve = vector_element;

   switch (evaluate_op)
   {

      default:
      case VALUE:
      {
         return Vector_Data_3D::evaluate (ve, z, i, j, evaluate_op);
      }

      case DX:
      {
         const Real latitude = get_coordinate (1, i);
         const Real eval = Vector_Data_3D::evaluate (ve, z, i, j, DY);
         return eval / (LATITUDE_LENGTH * cos (latitude * DEGREE_TO_RADIAN));
      }

      case DY:
      {
         const Real eval = Vector_Data_3D::evaluate (ve, z, i, j, DX);
         return eval / LATITUDE_LENGTH;
      }

      case DX2:
      {
         const Real latitude = get_coordinate (1, i);
         const Real c = cos (latitude * DEGREE_TO_RADIAN);
         const Real longitude_length = LATITUDE_LENGTH * c;
         const Real eval = Vector_Data_3D::evaluate (ve, z, i, j, DY2);
         return eval / (longitude_length * longitude_length);
      }

      case DY2:
      {
         const Real eval = Vector_Data_3D::evaluate (ve, z, i, j, DX2);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH);
      }

      case DXY:
      {
         const Real latitude = get_coordinate (1, i);
         const Real c = cos (latitude * DEGREE_TO_RADIAN);
         const Real eval = Vector_Data_3D::evaluate (ve, z, i, j, DXY);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH * c);
      }

      case LAPLACIAN:
      {
         const Real latitude = get_coordinate (1, i);
         const Real c = cos (latitude * DEGREE_TO_RADIAN);
         const Real dz2 = Vector_Data_3D::evaluate (ve, z, i, j, DY2);
         const Real dx2 = Vector_Data_3D::evaluate (ve, z, i, j, DY2);
         const Real dy2 = Vector_Data_3D::evaluate (ve, z, i, j, DX2);
         return (dx2 / (c * c) + dy2) / (LATITUDE_LENGTH * LATITUDE_LENGTH);
      }

   }

}

Real
Geodetic_Vector_Data_3D::evaluate (const Integer vector_element,
                                   const Integer k,
                                   const Integer i,
                                   const Integer j,
                                   const Evaluate_Op evaluate_op) const
{

   const Integer& ve = vector_element;

   switch (evaluate_op)
   {

      default:
      case VALUE:
      {
         return Vector_Data_3D::evaluate (ve, k, i, j, evaluate_op);
      }

      case DX:
      {
         const Real latitude = get_coordinate (1, i);
         const Real eval = Vector_Data_3D::evaluate (ve, k, i, j, DY);
         return eval / (LATITUDE_LENGTH * cos (latitude * DEGREE_TO_RADIAN));
      }

      case DY:
      {
         const Real eval = Vector_Data_3D::evaluate (ve, k, i, j, DX);
         return eval / LATITUDE_LENGTH;
      }

      case DX2:
      {
         const Real latitude = get_coordinate (1, i);
         const Real c = cos (latitude * DEGREE_TO_RADIAN);
         const Real longitude_length = LATITUDE_LENGTH * c;
         const Real eval = Vector_Data_3D::evaluate (ve, k, i, j, DY2);
         return eval / (longitude_length * longitude_length);
      }

      case DY2:
      {
         const Real eval = Vector_Data_3D::evaluate (ve, k, i, j, DX2);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH);
      }

      case DXY:
      {
         const Real latitude = get_coordinate (1, i);
         const Real c = cos (latitude * DEGREE_TO_RADIAN);
         const Real eval = Vector_Data_3D::evaluate (ve, k, i, j, DXY);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH * c);
      }

      case LAPLACIAN:
      {
         const Real latitude = get_coordinate (1, i);
         const Real c = cos (latitude * DEGREE_TO_RADIAN);
         const Real dz2 = Vector_Data_3D::evaluate (ve, k, i, j, DY2);
         const Real dx2 = Vector_Data_3D::evaluate (ve, k, i, j, DY2);
         const Real dy2 = Vector_Data_3D::evaluate (ve, k, i, j, DX2);
         return (dx2 / (c * c) + dy2) / (LATITUDE_LENGTH * LATITUDE_LENGTH);
      }

   }

}

Real
Geodetic_Vector_Data_3D::evaluate (const Integer vector_element,
                                   const Real z,
                                   const Real latitude,
                                   const Real longitude,
                                   const Evaluate_Op evaluate_op) const
{

   const Integer& ve = vector_element;

   switch (evaluate_op)
   {

      default:
      case VALUE:
      {
         return Vector_Data_3D::evaluate (ve, z, latitude, longitude, evaluate_op);
      }

      case DX:
      {
         Real eval = Vector_Data_3D::evaluate (ve, z, latitude, longitude, DY);
         return eval / (LATITUDE_LENGTH * cos (latitude * DEGREE_TO_RADIAN));
      }

      case DY:
      {
         Real eval = Vector_Data_3D::evaluate (ve, z, latitude, longitude, DX);
         return eval / LATITUDE_LENGTH;
      }

      case DX2:
      {
         Real c = cos (latitude * DEGREE_TO_RADIAN);
         Real longitude_length = LATITUDE_LENGTH * c;
         Real eval = Vector_Data_3D::evaluate (ve, z, latitude, longitude, DY2);
         return eval / (longitude_length * longitude_length);
      }

      case DY2:
      {
         Real eval = Vector_Data_3D::evaluate (ve, z, latitude, longitude, DX2);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH);
      }

      case DXY:
      {
         Real c = cos (latitude * DEGREE_TO_RADIAN);
         Real eval = Vector_Data_3D::evaluate (ve, z, latitude, longitude, DXY);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH * c);
      }

      case LAPLACIAN:
      {
         Real c = cos (latitude * DEGREE_TO_RADIAN);
         Real dz2 = Vector_Data_3D::evaluate (ve, z, latitude, longitude, DY2);
         Real dx2 = Vector_Data_3D::evaluate (ve, z, latitude, longitude, DY2);
         Real dy2 = Vector_Data_3D::evaluate (ve, z, latitude, longitude, DX2);
         return (dx2 / (c * c) + dy2) / (LATITUDE_LENGTH * LATITUDE_LENGTH);
      }

   }

}

Real
Geodetic_Vector_Data_3D::evaluate (const Integer vector_element,
                                   const Integer k,
                                   const Real latitude,
                                   const Real longitude,
                                   const Evaluate_Op evaluate_op) const
{

   const Integer& ve = vector_element;

   switch (evaluate_op)
   {

      default:
      case VALUE:
      {
         return Vector_Data_3D::evaluate (ve, k, latitude, longitude, evaluate_op);
      }

      case DX:
      {
         Real eval = Vector_Data_3D::evaluate (ve, k, latitude, longitude, DY);
         return eval / (LATITUDE_LENGTH * cos (latitude * DEGREE_TO_RADIAN));
      }

      case DY:
      {
         Real eval = Vector_Data_3D::evaluate (ve, k, latitude, longitude, DX);
         return eval / LATITUDE_LENGTH;
      }

      case DXY:
      {
         Real c = cos (latitude * DEGREE_TO_RADIAN);
         Real eval = Vector_Data_3D::evaluate (ve, k, latitude, longitude, DXY);
         return eval / (LATITUDE_LENGTH * LATITUDE_LENGTH * c);
      }

   }

}

void
Geodetic_Vector_Data_3D::subtract_zonal_mean (const Integer vector_element)
{
   subtract_y_mean (vector_element);
}

void
Geodetic_Vector_Data_3D::subtract_meridional_mean (const Integer vector_element)
{
   subtract_x_mean (vector_element);
}

void
Geodetic_Vector_Data_3D::subtract_horizontal_mean (const Integer vector_element)
{
   subtract_xy_mean (vector_element);
}

Geodetic_Scalar_Data_2D::Geodetic_Scalar_Data_2D (const Size_2D& size_2d,
                                                  const Domain_2D& domain_2d,
                                                  const bool periodic_longitude)
   : Vector_Data_2D (1, size_2d, domain_2d, false, periodic_longitude),
     Scalar_Data_2D (size_2d, domain_2d, false, periodic_longitude),
     Geodetic_Vector_Data_2D (1, size_2d, domain_2d, periodic_longitude)
{
}

Geodetic_Scalar_Data_2D::Geodetic_Scalar_Data_2D (const Tuple tuple_latitude,
                                                  const Tuple tuple_longitude,
                                                  const bool periodic_longitude)
   : Vector_Data_2D (1, tuple_latitude, tuple_longitude,
                     false, periodic_longitude),
     Scalar_Data_2D (tuple_latitude, tuple_longitude,
                     false, periodic_longitude),
     Geodetic_Vector_Data_2D (1, tuple_latitude, tuple_longitude,
                              periodic_longitude)
{
}

Geodetic_Scalar_Data_2D::~Geodetic_Scalar_Data_2D ()
{
}

void
Geodetic_Scalar_Data_2D::set_datum (const Integer node_latitude,
                                    const Integer node_longitude,
                                    const Real datum)
{
   Geodetic_Vector_Data_2D::set_datum (0, node_latitude, node_longitude, datum);
}

const Real&
Geodetic_Scalar_Data_2D::get_datum (const Integer node_latitude,
                                    const Integer node_longitude) const
{
   return Geodetic_Vector_Data_2D::get_datum (0, node_latitude, node_longitude);
}

Real
Geodetic_Scalar_Data_2D::evaluate (const Real latitude,
                                   const Real longitude,
                                   const Evaluate_Op evaluate_op) const
{
   return Geodetic_Vector_Data_2D::evaluate (
      0, latitude, longitude, evaluate_op);
}

Geodetic_Scalar_Data_3D::Geodetic_Scalar_Data_3D (const Size_3D& size_3d,
                                                  const Domain_3D& domain_3d,
                                                  const bool periodic_longitude)
   : Vector_Data_3D (1, size_3d, domain_3d, false, false, periodic_longitude),
     Scalar_Data_3D (size_3d, domain_3d, false, false, periodic_longitude),
     Geodetic_Vector_Data_3D (1, size_3d, domain_3d, periodic_longitude)
{
}

Geodetic_Scalar_Data_3D::Geodetic_Scalar_Data_3D (const Tuple tuple_z,
                                                  const Tuple tuple_latitude,
                                                  const Tuple tuple_longitude,
                                                  const bool periodic_longitude)
   : Vector_Data_3D (1, tuple_z, tuple_latitude, tuple_longitude,
                     false, false, periodic_longitude),
     Scalar_Data_3D (tuple_z, tuple_latitude, tuple_longitude,
                     false, false, periodic_longitude),
     Geodetic_Vector_Data_3D (1, tuple_z, tuple_latitude,
                              tuple_longitude, periodic_longitude)
{
}

Geodetic_Scalar_Data_3D::Geodetic_Scalar_Data_3D (const Tuple tuple_z,
                                                  const Size_2D& size_2d,
                                                  const Domain_2D& domain_2d,
                                                  const bool periodic_longitude)
   : Vector_Data_3D (1, tuple_z, size_2d, domain_2d,
                     false, false, periodic_longitude),
     Scalar_Data_3D (tuple_z, size_2d, domain_2d,
                     false, false, periodic_longitude),
     Geodetic_Vector_Data_3D (1, tuple_z, size_2d, domain_2d,
                              periodic_longitude)
{
}

Geodetic_Scalar_Data_3D::~Geodetic_Scalar_Data_3D ()
{
}


void
Geodetic_Scalar_Data_3D::set_datum (const Integer node_z,
                                    const Integer node_latitude,
                                    const Integer node_longitude,
                                    const Real datum)
{
   Geodetic_Vector_Data_3D::set_datum (
      0, node_z, node_latitude, node_longitude, datum);
}

const Real&
Geodetic_Scalar_Data_3D::get_datum (const Integer node_z,
                                    const Integer node_latitude,
                                    const Integer node_longitude) const
{
   return Geodetic_Vector_Data_3D::get_datum (
      0, node_z, node_latitude, node_longitude);
}

Real
Geodetic_Scalar_Data_3D::evaluate (const Real z,
                                   const Real latitude,
                                   const Real longitude,
                                   const Evaluate_Op evaluate_op) const
{
   return Geodetic_Vector_Data_3D::evaluate (
      0, z, latitude, longitude, evaluate_op);
}

Real
Geodetic_Scalar_Data_3D::evaluate (const Integer k,
                                   const Real latitude,
                                   const Real longitude,
                                   const Evaluate_Op evaluate_op) const
{
   return Geodetic_Vector_Data_3D::evaluate (
      0, k, latitude, longitude, evaluate_op);
}

Track::Track ()
   : lat_long_spline_ptr (NULL)
{
}

Track::~Track ()
{
   if (lat_long_spline_ptr != NULL) { delete lat_long_spline_ptr; }
}

void
Track::add (const Real& tau,
            const Lat_Long& lat_long)
{
   insert (make_pair (tau, lat_long));
}

void
Track::okay ()
{

   const Integer n = size ();
   const Real last_tau = rbegin ()->first;

   if (n < 2) { return; }
   if (lat_long_spline_ptr != NULL) { delete lat_long_spline_ptr; }

   Tuple tau_tuple;
   for (Track::iterator i = begin (); i != end (); i++)

   lat_long_spline_ptr = new Vector_Data_1D (2, n, Domain_1D (0, last_tau));

   for (Track::iterator i = begin (); i != end (); i++)
   {

      const Real tau = i->first;
      const Lat_Long& lat_long = i->second;

      const Integer node = distance (begin (), i);
      const Real& latitude = lat_long.latitude;
      const Real& longitude = lat_long.longitude;

      lat_long_spline_ptr->modify_coordinate_tuple (node, tau);
      lat_long_spline_ptr->set_datum (0, node, latitude);
      lat_long_spline_ptr->set_datum (1, node, longitude);

   }

   const gsl_interp_type* interp_type = (n > 3 ?
      gsl_interp_cspline : gsl_interp_linear);
   lat_long_spline_ptr->set_interpolation (interp_type);

}

Lat_Long
Track::get_lat_long (const Real tau,
                     const bool forbid_extrapolate) const
{

   if (forbid_extrapolate)
   {

      const Real min_tau = begin ()->first;
      const Real max_tau = rbegin ()->first;

      if (tau < min_tau || tau > max_tau)
      {
         return Lat_Long (GSL_NAN, GSL_NAN);
      }

   }

   const Real latitude = lat_long_spline_ptr->evaluate (0, tau);
   const Real longitude = lat_long_spline_ptr->evaluate (1, tau);

   return Lat_Long (latitude, longitude);

}

Geodetic_Mesh::Geodetic_Mesh (const Size_2D& size_2d,
                              const Domain_2D& domain_2d)
   : Mesh_2D (size_2d, domain_2d)
{
}

Geodetic_Mesh::Geodetic_Mesh (const Simple_Mesh_2D& simple_mesh_2d,
                              const Size_2D& size_2d,
                              const Domain_2D& domain_2d)
   : Mesh_2D (size_2d, domain_2d, simple_mesh_2d)
{
}

Geodetic_Mesh::Geodetic_Mesh (const Simple_Mesh_2D& simple_mesh_2d_a,
                              const Simple_Mesh_2D& simple_mesh_2d_b,
                              const Size_2D& size_2d,
                              const Domain_2D& domain_2d)
   : Mesh_2D (size_2d, domain_2d, simple_mesh_2d_a, simple_mesh_2d_b)
{
}

Geodetic_Mesh::Geodetic_Mesh (const Simple_Mesh_2D& simple_mesh_2d_a,
                              const Simple_Mesh_2D& simple_mesh_2d_b,
                              const Simple_Mesh_2D& simple_mesh_2d_c,
                              const Size_2D& size_2d,
                              const Domain_2D& domain_2d)
   : Mesh_2D (size_2d, domain_2d, simple_mesh_2d_a,
              simple_mesh_2d_b, simple_mesh_2d_c)
{
}

Geodetic_Mesh::Geodetic_Mesh (const Simple_Mesh_2D& simple_mesh_2d_a,
                              const Simple_Mesh_2D& simple_mesh_2d_b,
                              const Simple_Mesh_2D& simple_mesh_2d_c,
                              const Simple_Mesh_2D& simple_mesh_2d_d,
                              const Size_2D& size_2d,
                              const Domain_2D& domain_2d)
   : Mesh_2D (size_2d, domain_2d, simple_mesh_2d_a,
              simple_mesh_2d_b, simple_mesh_2d_c, simple_mesh_2d_d)
{
}

void
Geodetic_Mesh::cairo (const RefPtr<Context> cr,
                      const Geodetic_Transform& transform) const
{

   const Real start_latitude = domain_2d.domain_x.start;
   const Real end_latitude = domain_2d.domain_x.end;
   const Real start_longitude = domain_2d.domain_y.start;
   const Real end_longitude = domain_2d.domain_y.end;
   const Real latitude_span = domain_2d.domain_x.get_span ();
   const Real longitude_span = domain_2d.domain_y.get_span ();
   const Real anchor_lat_a = start_latitude + latitude_span * 0.2;
   const Real anchor_long_a = start_longitude + longitude_span * 0.2;
   const Real anchor_lat_b = start_latitude + latitude_span * 0.8;
   const Real anchor_long_b = start_longitude + longitude_span * 0.8;
   const Lat_Long anchor_lat_long_a (anchor_lat_a, anchor_long_a);
   const Lat_Long anchor_lat_long_b (anchor_lat_b, anchor_long_b);

   Mesh_2D::render (cr, transform);
   render_label_lat_long (cr, transform, 1, anchor_lat_long_a, "%.0f");
   render_label_lat_long (cr, transform, 1, anchor_lat_long_b, "%.0f");

}

