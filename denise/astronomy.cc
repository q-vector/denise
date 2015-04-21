//
// astronomy.cc
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

#include "astronomy.h"

using namespace std;
using namespace denise;

Real
Cosmos::get_ascending_node (const Dtime& dtime) const
{
   return 0;
}

Real
Cosmos::get_inclination (const Dtime& dtime) const
{
   return 0;
}

Real
Cosmos::get_eccentric_anormaly (const Real first_guess,
                                const Real mean_anormaly,
                                const Real eccentricity,
                                const Real epsilon) const
{

   const Real M = mean_anormaly;
   const Real e = eccentricity;
   const Real E_0 = first_guess;
   const Real E_1 = E_0 - (E_0 - e * sin (E_0) - M) / (1 - e * cos (E_0));

   if (fabs (E_0 - E_1) < epsilon) { return E_1; }
   else { return get_eccentric_anormaly (E_1, M, e, epsilon); }

}

Real
Cosmos::get_eccentric_anormaly (const Dtime& dtime,
                                const Real epsilon) const
{
   const Real M = get_mean_anormaly (dtime);
   const Real e = get_eccentricity (dtime);
   if (e < 0.05) { return M + e * sin (M) * (1 + e * cos (M)); }
   else { return get_eccentric_anormaly (M, M, e, epsilon); }
}

pair<Real, Real>
Cosmos::get_r_v (const Dtime& dtime) const
{

   const Real E = get_eccentric_anormaly (dtime);
   const Real e = get_eccentricity (dtime);

   const Real x = a * (cos (E) - e);
   const Real y = a * (sqrt (1 - e*e) * sin (E));
   const Real r = sqrt (x*x + y*y);
   const Real v = atan2 (y, x);

   return make_pair (r, v);

}

Point_3D
Cosmos::get_ecliptic_heliocentric (const Dtime& dtime) const
{

   const Real N = get_ascending_node (dtime);
   const Real i = get_inclination (dtime);
   const Real w = get_perihelion (dtime);

   const pair<Real, Real> r_v = get_r_v (dtime);
   const Real& r = r_v.first;
   const Real& v = r_v.second;
   const Real vw = v + w;

   const Real x = r * (cos (N) * cos (vw) - sin (N) * sin (vw) * cos (i));
   const Real y = r * (sin (N) * cos (vw) + cos (N) * sin (vw) * cos (i));
   const Real z = r * sin (vw) * sin (i);

   return Point_3D (z, x, y);

}

Point_3D
Cosmos::get_ecliptic_geocentric (const Dtime& dtime) const
{

   const Sun sun;
   const Point_3D& sun_pos = sun.get_equatorial_geocentric (dtime);
   const Point_3D& pos = get_ecliptic_geocentric (dtime);

   const Real x = pos.x + sun_pos.x;
   const Real y = pos.y + sun_pos.y;

   return Point_3D (pos.z, pos.x + sun_pos.x, pos.x + sun_pos.y);

}

Point_3D
Cosmos::get_equatorial_geocentric (const Dtime& dtime) const
{

   const Point_3D& point_3d = get_ecliptic_geocentric (dtime);
   const Real ecliptic = Sun::get_ecliptic (dtime);

   const Real& z = point_3d.z;
   const Real& y = point_3d.y;

   const Real yy = y * cos (ecliptic) - z * sin (ecliptic);
   const Real zz = y * sin (ecliptic) + z * cos (ecliptic);

   return (zz, point_3d.x, yy);

}

Cosmos::Cosmos (const Real a)
   : a (a)
{
}

Cosmos::Cosmos (const Cosmos& cosmos)
   : a (cosmos.a)
{
}

Real
Cosmos::get_day_number (const Dtime& dtime)
{

   const Integer yy = dtime.get_year ();
   const Integer mm = dtime.get_month ();
   const Integer dd = dtime.get_day ();

   const Integer d = 367*yy - 7*(yy+(mm+9)/12)/4 + 275*mm/9 + dd - 730530;

   const Real ut = modulo (dtime.t, 24);
   return Real (d) + ut / 24;

}

Lat_Long
Cosmos::get_lat_long (const Dtime& dtime) const
{

/*
cout << "ha" << endl;
   const Real day = Cosmos::get_day_number (dtime);

   const Real ecliptic = (23.4393 - 3.563e-7 * day) * DEGREE_TO_RADIAN;
   const Real w = (282.9404 + 4.70935e-5 * day) * DEGREE_TO_RADIAN;
   const Real e = (0.016709 - 1.151e-9 * day) * DEGREE_TO_RADIAN;
   const Real M = (356.0470 + 0.9856002585 * day) * DEGREE_TO_RADIAN;
   const Real E = M + e * sin (M) * (1 + e * cos (M));

   const Sun sun;
   const Real Ms = sun.get_mean_anormaly (dtime);
   const Real ws = sun.get_perihelion (dtime);
   const Real Ls = Ms + ws;

//   const Real E = get_eccentric_anormaly (dtime);
//   const Real w = get_perihelion (dtime);
//   const Real e = get_eccentricity (dtime);
//   const Real ecliptic = Sun::get_ecliptic (dtime);

   const Real xv = cos (E) - e;
   const Real yv = sqrt (1 - e*e) * sin (E);

   const Real v = atan2 (yv, xv);
   const Real r = sqrt (xv*xv + yv*yv);

   const Real vw = v + w;
   const Real xs = r * cos (vw);
   const Real ys = r * sin (vw);

   const Real xe = xs;
   const Real ye = ys * cos (ecliptic);
   const Real ze = ys * sin (ecliptic);

   const Real gmst0 = ((Ls * RADIAN_TO_DEGREE) + 180) / 15;
   const Real lst = (gmst0 + modulo (dtime.t, 24) * 15);
   const Real d_lambda = lst * 15 * DEGREE_TO_RADIAN;

   const Real local_hour = modulo (dtime.t, 24);
   const Real lambda = atan2 (ye , xe) - d_lambda;
   const Real phi = atan2 (ze, sqrt (xe*xe + ye*ye));

   Lat_Long lat_long (phi * RADIAN_TO_DEGREE, lambda * RADIAN_TO_DEGREE);
cout << "from" << lat_long << endl;
   lat_long.standardize ();
cout << "to" << lat_long << endl;

   return lat_long;
*/

   const Point_3D& point_3d = get_equatorial_geocentric (dtime);
   const Real& z = point_3d.z;
   const Real& x = point_3d.x;
   const Real& y = point_3d.y;

   const Sun sun;
   const Real Ms = sun.get_mean_anormaly (dtime);
   const Real ws = sun.get_perihelion (dtime);
   const Real Ls = Ms + ws;
   const Real gmst0 = (Ls + M_PI);
   const Real lst = gmst0 + modulo (dtime.t, 24) * 15 * DEGREE_TO_RADIAN;
   const Real d_lambda = lst;

   const Real lambda = atan2 (y , x) - d_lambda;
   const Real phi = atan2 (z, sqrt (x*x + y*y));

   Lat_Long ll (phi * RADIAN_TO_DEGREE, lambda * RADIAN_TO_DEGREE);
   ll.standardize ();

   return ll;

}

Real
Cosmos::get_distance (const Dtime& dtime) const
{
   const Point_3D& point_3d = get_equatorial_geocentric (dtime);
   const Real& z = point_3d.z;
   const Real& x = point_3d.x;
   const Real& y = point_3d.y;
   return sqrt (x*x + y*y + z*z);
}

Real
Cosmos::get_zenith_angle (const Dtime& dtime,
                          const Lat_Long& lat_long) const
{

/*
//   const Real day = Cosmos::get_day_number (dtime);
//
//   const Point_3D& point_3d = get_equatorial_geocentric (dtime);
//   const Real& z = point_3d.z;
//   const Real& x = point_3d.x;
//   const Real& y = point_3d.y;
//
//   const Real lambda = atan2 (point_3d.y , point_3d.x);
//   const Real phi = atan2 (point_3d.z, sqrt (x*x + y*y));
//   const Lat_Long ll (phi * RADIAN_TO_DEGREE, lambda * RADIAN_TO_DEGREE);


   const Real day = Cosmos::get_day_number (dtime);

   const Real ecliptic = (23.4393 - 3.563e-7 * day) * DEGREE_TO_RADIAN;
   const Real w = (282.9404 + 4.70935e-5 * day) * DEGREE_TO_RADIAN;
   const Real e = (0.016709 - 1.151e-9 * day) * DEGREE_TO_RADIAN;
   const Real M = (356.0470 + 0.9856002585 * day) * DEGREE_TO_RADIAN;
   const Real E = M + e * sin (M) * (1 + e * cos (M));

   const Real xv = cos (E) - e;
   const Real yv = sqrt (1 - e*e) * sin (E);

   const Real v = atan2 (yv, xv);
   const Real r = sqrt (xv*xv + yv*yv);

   const Real vw = v + w;
   const Real xs = r * cos (vw);
   const Real ys = r * sin (vw);

   const Real xe = xs;
   const Real ye = ys * cos (ecliptic);
   const Real ze = ys * sin (ecliptic);

   const Real local_hour = modulo (dtime.t, 24);
   const Real lambda = atan2 (ye , xe) - local_hour * 15 * DEGREE_TO_RADIAN;
   const Real phi = atan2 (ze, sqrt (xe*xe + ye*ye));
*/

   const Point_3D& point_3d = get_equatorial_geocentric (dtime);
   const Real& z = point_3d.z;
   const Real& x = point_3d.x;
   const Real& y = point_3d.y;

   const Sun sun;
   const Real Ms = sun.get_mean_anormaly (dtime);
   const Real ws = sun.get_perihelion (dtime);
   const Real Ls = Ms + ws;
   const Real gmst0 = (Ls + M_PI);
   const Real lst = gmst0 + modulo (dtime.t, 24) * 15 * DEGREE_TO_RADIAN;
   const Real d_lambda = lst;

   const Real lambda = atan2 (y , x) - d_lambda;
   const Real phi = atan2 (z, sqrt (x*x + y*y));

   Lat_Long ll (phi * RADIAN_TO_DEGREE, lambda * RADIAN_TO_DEGREE);
   ll.standardize ();

   if (ll == lat_long) { return 90; }

   const Real d = sqrt (x*x + y*y + z*z) * 149597870700.0;

   const Real a = EARTH_RADIUS;
   const Real theta_star = Geodesy::get_distance (lat_long, ll);
   const Real theta = theta_star / LATITUDE_LENGTH * DEGREE_TO_RADIAN;

   const Real d_star = sqrt (a*a + d*d - 2*a*d* cos (theta));
   Real zenith_angle = asin (d / d_star * sin (theta)) * RADIAN_TO_DEGREE;

   if (d_star*d_star > a*a + d*d) { zenith_angle = 180 - zenith_angle; }
   return zenith_angle;

}

Dtime
Cosmos::get_rise_or_set (const Lat_Long& lat_long,
                         const Dtime& dtime_0,
                         const Dtime& dtime_1,
                         const Real zenith_angle,
                         const Real epsilon) const
{

   const Dtime& dt0 = dtime_0;
   const Dtime& dt1 = dtime_1;
   const Real zenith = zenith_angle;
   const Dtime dtm ((dtime_0.t + dtime_1.t) / 2);
   const Real time_span = (dtime_1.t - dtime_0.t);

const string format ("%Y.%m.%d %H:%M:%S");
const string& str_0 = dtime_0.get_string (format, true);
const string& str_1 = dtime_1.get_string (format, true);
const string& str_m = dtm.get_string (format, true);

cout << endl;
   if (time_span < epsilon) { return dtm; }
   else
   {
      
      const Real zenith_0 = get_zenith_angle (dtime_0, lat_long);
      const Real zenith_1 = get_zenith_angle (dtime_1, lat_long);
      const Real zenith_m = get_zenith_angle (dtm, lat_long);

      const Real b_0 = (zenith_0 - zenith_angle);
      const Real b_1 = (zenith_1 - zenith_angle);
      const Real b_m = (zenith_m - zenith_angle);

cout << str_0 << " " << zenith_0 << " " << b_0 << endl;
cout << str_m << " " << zenith_m << " " << b_m << endl;
cout << str_1 << " " << zenith_1 << " " << b_1 << endl;

      if ((b_0 * b_m) < 0)
      {
         return get_rise_or_set (lat_long, dt0, dtm, zenith, epsilon);
      }
      else
      {
         return get_rise_or_set (lat_long, dtm, dt1, zenith, epsilon);
      }

   }

}

Dtime
Cosmos::get_next_rise (const Dtime& dtime,
                       const Lat_Long& lat_long,
                       const Real zenith_angle,
                       const Real epsilon) const
{

   const Real longitude = lat_long.longitude;
   const Real h = 12 - (longitude / 15);
   const Dtime dtime_noon_approx (modulo (dtime.t, 24) + h);

   const Real earth_tilt = 23.439291111;
   const Real t_0 = Dtime ("197003210056").t;
   const Real& t = dtime.t;
   const Real over_head_latitude = earth_tilt * sin ((t - t_0) / 8765.812776);
   

   for (Dtime dt = Dtime (dtime.t); dt.t <= (dtime.t + 30); dt.t += 0.5)
   {

      const Dtime dtime_0 (dt.t);
      const Real zenith_0 = get_zenith_angle (dtime_0, lat_long);
      if (zenith_0 == zenith_angle) { return dtime_0; }

      const Dtime dtime_1 (dt.t + 0.5);
      const Real zenith_1 = get_zenith_angle (dtime_1, lat_long);
      if (zenith_1 == zenith_angle) { return dtime_1; }

      if (zenith_0 > zenith_angle && zenith_1 < zenith_angle)
      {
         return get_rise_or_set (lat_long, dtime_0,
            dtime_1, zenith_angle, epsilon);
      }

   }

   return Dtime (GSL_NAN);

}

Dtime
Cosmos::get_next_set (const Dtime& dtime,
                      const Lat_Long& lat_long,
                      const Real zenith_angle,
                      const Real epsilon) const
{

   for (Dtime dt = Dtime (dtime.t); dt.t <= (dtime.t + 30); dt.t += 0.5)
   {

      const Dtime dtime_0 (dt.t);
      const Real zenith_0 = get_zenith_angle (dt, lat_long);
      if (zenith_0 == zenith_angle) { return dtime_0; }

      const Dtime dtime_1 (dt.t + 0.5);
      const Real zenith_1 = get_zenith_angle (dtime_1, lat_long);
      if (zenith_1 == zenith_angle) { return dtime_1; }

      if (zenith_0 < zenith_angle && zenith_1 > zenith_angle)
      {
         return get_rise_or_set (lat_long, dtime_0,
            dtime_1, zenith_angle, epsilon);
      }

   }

   return Dtime (GSL_NAN);

}

Real
Sun::get_eccentricity (const Dtime& dtime) const
{
   const Real d = Cosmos::get_day_number (dtime);
   return 0.016709 - 1.515e-9 * d;
}

Sun::Sun ()
   : Cosmos (1)
{
}

Real
Sun::get_ecliptic (const Dtime& dtime)
{
   const Real d = Cosmos::get_day_number (dtime);
   return (23.4393 - 3.563e-7 * d) * DEGREE_TO_RADIAN;
}

Real
Sun::get_mean_anormaly (const Dtime& dtime) const
{
   const Real d = Cosmos::get_day_number (dtime);
   return (356.0470 + 0.9856002585 * d) * DEGREE_TO_RADIAN;
}

Real
Sun::get_perihelion (const Dtime& dtime) const
{
   const Real d = Cosmos::get_day_number (dtime);
   return (282.9404 + 4.70935e-5 * d) * DEGREE_TO_RADIAN;
}

Point_3D
Sun::get_equatorial_geocentric (const Dtime& dtime) const
{

   const Real w = get_perihelion (dtime);
   const Real ecliptic = get_ecliptic (dtime);

   const pair<Real, Real> r_v = get_r_v (dtime);
   const Real& r = r_v.first;
   const Real& v = r_v.second;
   const Real sun_longitude = v + w;

   const Real xs = r * cos (sun_longitude);
   const Real ys = r * sin (sun_longitude);

   const Real y = ys * cos (ecliptic);
   const Real z = ys * sin (ecliptic);

   return Point_3D (z, xs, y);

}

Point_3D
Perturbed_Cosmos::get_ecliptic_heliocentric (const Dtime& dtime) const
{

   const Point_3D& point_3d = Cosmos::get_ecliptic_heliocentric (dtime);
   const Real& x = point_3d.x;
   const Real& y = point_3d.y;
   const Real& z = point_3d.z;

   Real r = sqrt (x*x + y*y + z*z);
   Real ecliptic_latitude = atan2 (y, x);
   Real ecliptic_longitude = atan2 (z, sqrt (x*x + y*y));

   apply_corrections (dtime, r, ecliptic_latitude, ecliptic_longitude);

   const Real xx = r * cos (ecliptic_latitude) * cos (ecliptic_longitude);
   const Real yy = r * cos (ecliptic_latitude) * sin (ecliptic_longitude);
   const Real zz = r * sin (ecliptic_latitude);

   return Point_3D (zz, xx, yy);

}

Real
Moon::get_ascending_node (const Dtime& dtime) const
{
   const Real d = Cosmos::get_day_number (dtime);
   return (125.1228 - 0.0529538083 * d) * DEGREE_TO_RADIAN;
}

Real
Moon::get_inclination (const Dtime& dtime) const
{
   return 5.1454;
}

Real
Moon::get_eccentricity (const Dtime& dtime) const
{
   return 0.054900;
}

Real
Moon::get_mean_anormaly (const Dtime& dtime) const
{
   const Real d = Cosmos::get_day_number (dtime);
   return (115.3654 + 13.0649929509 * d) * DEGREE_TO_RADIAN;
}

Real
Moon::get_perihelion (const Dtime& dtime) const
{
   const Real d = Cosmos::get_day_number (dtime);
   return (318.0634 + 0.1643573223 * d) * DEGREE_TO_RADIAN;
}

void
Moon::apply_corrections (const Dtime& dtime,
                         Real& r,
                         Real& ecliptic_latitude,
                         Real& ecliptic_longitude) const
{

   const Sun sun;
   const Real Ms = sun.get_mean_anormaly (dtime);
   const Real ws = sun.get_perihelion (dtime);
   const Real Mm = get_mean_anormaly (dtime);
   const Real wm = get_perihelion (dtime);
   const Real Nm = get_ascending_node (dtime);

   const Real Ls = Ms + ws;
   const Real Lm = Mm + wm + Nm;
   const Real D = Lm - Ls;
   const Real F = Lm - Nm;

   const Real d_latitude = (-0.173 * sin (F - 2*D) +
      -0.055 * sin (Mm - F - 2*D) + -0.046 * sin (Mm + F - 2*D) +
      0.033 * sin (F + 2*D) + 0.017 * sin (2*Mm + F));
   const Real d_longitude = (-1.274 * sin (Mm - 2*D) + 0.658 * sin (2*D) +
      -0.186 * sin (Ms) + -0.059 * sin (2*Mm - 2*D) +
      -0.057 * sin (Mm - 2*D + Ms) + 0.053 * sin (Mm + 2*D) +
      0.046 * sin (2*D - Ms) + 0.041 * sin (Mm - Ms) + -0.035 * sin (D) +
      -0.031 * sin (Mm + Ms) + -0.015 * sin (2*F - 2*D) +
      0.011 * sin (Mm - 4*D));

   ecliptic_latitude += d_latitude * DEGREE_TO_RADIAN;
   ecliptic_longitude += d_longitude * DEGREE_TO_RADIAN;
   r += (-0.58 * cos (Mm - 2*D) + -0.46 * cos (2*D));

}

Point_3D
Moon::get_ecliptic_geocentric (const Dtime& dtime) const
{
   return Perturbed_Cosmos::get_ecliptic_heliocentric (dtime);
}

Point_3D
Moon::get_equatorial_geocentric (const Dtime& dtime) const
{
   return Perturbed_Cosmos::get_ecliptic_geocentric (dtime);
}

Zenith_Field::Zenith_Field (const Cosmos& cosmos,
                            const Dtime& dtime)
   : cosmos (cosmos),
     dtime (dtime)
{
   if (dtime.is_nat ()) { this->dtime = Dtime (); }
}

const Dtime&
Zenith_Field::get_time () const
{
   return dtime;
}

Real
Zenith_Field::evaluate (const Integer vector_element,
                        const Real latitude,
                        const Real longitude,
                        const Evaluate_Op evaluate_op) const
{
   const Lat_Long lat_long (latitude, longitude);
   return cosmos.get_zenith_angle (dtime, lat_long);
}

/*
Real
Solar::get_earth_tilt (const Dtime& time) const
{

   switch (earth_tilt_model)
   {

      default:
      EARTH_TILT_CONSTANT:
      {
         return 23.439291111;
      }

      EARTH_TILT_IAU:
      {
         const Dtime epoch ("2000010100");
         const Real c = ((time.t - epoch.t) / 876600);
         return (84381.448 + (-46.84024 + (-5.9e-4  + 1.813e-3*c)*c)*c) / 3600;
      }

      EARTH_TILT_WITTMANN:
      {
         const Dtime epoch ("2000010100");
         const Real c = ((time.t - epoch.t) / 876600);
         return 23.496932 + -0.86 * sin (0.01532 * (c + 4.4));
      }

   }

}

Real
Solar::get_rise_set_t (const Dtime& start_time,
                       const Dtime& end_time,
                       const Real altitude,
                       const Real epsilon_t)
{

   const Real start_altitude = get_altitude (lat_long, start_time);
   const Real end_altitude = get_altitude (lat_long, end_time);
   const bool rising = (end_altitude > start_altitude);

   if ((start_altitude - altitude) - (end_altitude - altitude) > 0)
   {
      throw Exception ("start/end_latitude same sign");
   }

   if (start_altitude == altitude)
   {
      return start_time.t;
   }

   if (end_altitude == altitude)
   {
      return end_time.t;
   }

   const Dtime middle_time ((start_time.t + end_time.t) / 2);
   const Real middle_altitude = get_altitude (lat_long, middle_time);

   if ((fabs (end_time.t - start_time.t) < epsilon_t) ||
       (middle_altitude == altitude))
   {
      return middle_time.t;
   }

   const Real start_d_altitude = start_altitude - altitude;
   const Real middle_d_altitude = middle_altitude - altitude;

   if (middle_d_altitude - start_d_altitude <= 0)
   {
      return get_rise_set_t (start_time, middle_time, altitude, epsilon_t);
   }
   else
   {
      return get_rise_set_t (middle_time, end_time, altitude, epsilon_t);
   }

}

Solar::Solar (const Earth_Tilt_Model earth_tilt_model)
   : earth_tilt_model (earth_tilt_model)
{
}

void
Solar::set_lat_long (const Dtime& dtime)
{

   const Real earth_tilt = get_earth_tilt (dtime);
   const Real t_0 = Dtime ("197003210056").t;
   const Real& t = dtime.t;
   
   lat_long.latitude = earth_tilt * sin ((t - t_0) / 8765.812776);
   lat_long.longitude = modulo (180 - (t * 15), 0, 360);

}

const Lat_Long&
Solar::get_lat_long (const Dtime& dtime)
{
   set_lat_long (dtime);
   return lat_long;
}

const Lat_Long&
Solar::get_lat_long () const
{
   return lat_long;
}

void
Solar::acquire_alt_az (Real& altitude,
                       Real& azimuth,
                       const Lat_Long& lat_long,
                       const Dtime& dtime)
{

   set_lat_long (dtime);
   const Lat_Long& sun_lat_long = this->lat_long;

   Geodesy geodesy (SPHERE);
   Journey journey (lat_long, sun_lat_long);
   geodesy.complete (journey);

   altitude = (90 - journey.distance / EARTH_RADIUS) * DEGREE_TO_RADIAN;
   azimuth = journey.azimuth_forward;

}

Real
Solar::get_altitude (const Lat_Long& lat_long,
                     const Dtime& dtime)
{

   set_lat_long (dtime);
   const Lat_Long& sun_lat_long = this->lat_long;

   Geodesy geodesy (SPHERE);
   Journey journey (lat_long, sun_lat_long);
   geodesy.complete (journey);

   return (90 - journey.distance / EARTH_RADIUS * RADIAN_TO_DEGREE);

}

pair<Dtime, Dtime>
Solar::get_rise_set_time_pair (const Lat_Long& lat_long,
                               const Dtime& dtime,
                               const Real altitude,
                               const Real epsilon_t)
{

   set_lat_long (dtime);

   const Real& l = lat_long.longitude;
   const Real& sl = this->lat_long.longitude;

   const Real local_hour = modulo ((l - sl) / 15, 0, 24);

   const Dtime start_time (dtime.t - local_hour);
   const Dtime midday (start_time.t + 12);
   const Dtime end_time (start_time.t + 24);

   Dtime rise_time (get_rise_set_t (start_time, midday, altitude, epsilon_t));
   Dtime set_time (get_rise_set_t (midday, end_time, altitude, epsilon_t));

   return make_pair (rise_time, set_time);

}

Sun_Altitude_Field::Sun_Altitude_Field (const Dtime& dtime,
                                        const Earth_Tilt_Model earth_tilt_model)
   : dtime (dtime)
{
}

void
Sun_Altitude_Field::set_time (const Dtime& dtime)
{
   this->dtime.t = dtime.t;
}

const Dtime&
Sun_Altitude_Field::get_time () const
{
   return dtime;
}

Real
Sun_Altitude_Field::evaluate (const Real latitude,
                              const Real longitude,
                              const Evaluate_Op evaluate_op) const
{
   Solar solar (EARTH_TILT_IAU);
   return solar.get_altitude (Lat_Long (latitude, longitude), dtime);
}
*/

