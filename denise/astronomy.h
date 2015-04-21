//
// astronomy.h
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

#ifndef DENISE_ASTRONOMY_H
#define DENISE_ASTRONOMY_H

#include <cmath>
#include <denise/basics.h>
#include <denise/dtime.h>
#include <denise/exception.h>
#include <denise/geodesy.h>
#include <denise/geometry.h>
#include <denise/util.h>

using namespace std;

namespace denise
{

   class Cosmos
   {

      private:

         const Real
         a;

      protected:

         virtual Real
         get_ascending_node (const Dtime& dtime) const;

         virtual Real
         get_inclination (const Dtime& dtime) const;

         virtual Real
         get_eccentricity (const Dtime& dtime) const = 0;

         Real
         get_eccentric_anormaly (const Real first_guess,
                                 const Real mean_anormaly,
                                 const Real eccentricity,
                                 const Real epsilon = 1e-8) const;

         Real
         get_eccentric_anormaly (const Dtime& dtime,
                                 const Real epsilon = 1e-8) const;

         pair<Real, Real>
         get_r_v (const Dtime& dtime) const; 

         virtual Point_3D
         get_ecliptic_heliocentric (const Dtime& dtime) const;

         virtual Point_3D
         get_ecliptic_geocentric (const Dtime& dtime) const;

         virtual Point_3D
         get_equatorial_geocentric (const Dtime& dtime) const;

      public:

         Cosmos (const Real a);

         Cosmos (const Cosmos& cosmos);

         static Real
         get_day_number (const Dtime& dtime);

         virtual Real
         get_mean_anormaly (const Dtime& dtime) const = 0;

         virtual Real
         get_perihelion (const Dtime& dtime) const = 0;

         Lat_Long
         get_lat_long (const Dtime& dtime) const;

         Real
         get_distance (const Dtime& dtime) const;

         Real
         get_zenith_angle (const Dtime& dtime,
                           const Lat_Long& lat_long) const;

         Dtime
         get_rise_or_set (const Lat_Long& lat_long,
                          const Dtime& dtime_0,
                          const Dtime& dtime_1,
                          const Real zenith_angle,
                          const Real epsilon) const;

         Dtime
         get_next_rise (const Dtime& dtime,
                        const Lat_Long& lat_long,
                        const Real zenith_angle = 90.6,
                        const Real epsilon = 0.0001) const;

         Dtime
         get_next_set (const Dtime& dtime,
                       const Lat_Long& lat_long,
                       const Real zenith_angle = 90.6,
                       const Real epsilon = 0.0001) const;

   };

   class Sun : public Cosmos
   {

      private:

         Real
         get_eccentricity (const Dtime& dtime) const;

      public:

         Sun ();

         static Real
         get_ecliptic (const Dtime& dtim);

         Real
         get_mean_anormaly (const Dtime& dtime) const;

         Real
         get_perihelion (const Dtime& dtime) const;

         Point_3D
         get_equatorial_geocentric (const Dtime& dtime) const;

   };

   class Perturbed_Cosmos : public Cosmos
   {

      protected:

         virtual void
         apply_corrections (const Dtime& dtime,
                            Real& r,
                            Real& ecliptic_latitude,
                            Real& ecliptic_longitude) const = 0;

      public:

         Point_3D
         get_ecliptic_heliocentric (const Dtime& dtime) const;

   };

   class Moon : public Perturbed_Cosmos
   {

      private:

         Real
         get_ascending_node (const Dtime& dtime) const;

         Real
         get_inclination (const Dtime& dtime) const;

         Real
         get_eccentricity (const Dtime& dtime) const;

         Real
         get_mean_anormaly (const Dtime& dtime) const;

         Real
         get_perihelion (const Dtime& dtime) const;

         void
         apply_corrections (const Dtime& dtime,
                            Real& r,
                            Real& ecliptic_latitude,
                            Real& ecliptic_longitude) const;

         Point_3D
         get_ecliptic_geocentric (const Dtime& time) const;

         Point_3D
         get_equatorial_geocentric (const Dtime& time) const;

      public:

         

   };

   class Zenith_Field : public Vector_Field_2D
   {

      private:

         const Cosmos&
         cosmos;

         Dtime
         dtime;

      public:

         Zenith_Field (const Cosmos& cosmos,
                       const Dtime& dtime = Dtime (GSL_NAN));

         const Dtime&
         get_time () const;

         Real
         evaluate (const Integer vector_element,
                   const Real latitude,
                   const Real longitude,
                   const Evaluate_Op evaluate_op = VALUE) const;

   };

/*
   enum Earth_Tilt_Model
   {
      EARTH_TILT_CONSTANT,
      EARTH_TILT_IAU,
      EARTH_TILT_WITTMANN
   };

   class Sun_Elevation_Field;

   class Solar
   {

      private:

         const Earth_Tilt_Model
         earth_tilt_model;

         Lat_Long
         lat_long;

         Real
         get_earth_tilt (const Dtime& dtime) const;

         Real
         get_rise_set_t (const Dtime& start_time,
                         const Dtime& end_time,
                         const Real altitude,
                         const Real epsilon_t);

      public:

         Solar (const Earth_Tilt_Model earth_tilt_model = EARTH_TILT_IAU);

         void
         set_lat_long (const Dtime& dtime);

         const Lat_Long&
         get_lat_long (const Dtime& dtime);

         const Lat_Long&
         get_lat_long () const;

         void
         acquire_alt_az (Real& altitude,
                         Real& azimuth,
                         const Lat_Long& lat_long,
                         const Dtime& dtime);

         Real
         get_altitude (const Lat_Long& lat_long,
                       const Dtime& dtime);

         pair<Dtime, Dtime>
         get_rise_set_time_pair (const Lat_Long& lat_long,
                                 const Dtime& dtime,
                                 const Real altitude,
                                 const Real epsilon_t = 0.001);

   };

   class Sun_Altitude_Field : public Scalar_Field_2D
   {

      private:

         Dtime
         dtime;

      public:

         Sun_Altitude_Field (const Dtime& dtime,
                             const Earth_Tilt_Model earth_tilt_model = EARTH_TILT_IAU);

         void
         set_time (const Dtime& dtime);

         const Dtime&
         get_time () const;

         virtual Real
         evaluate (const Real latitude,
                   const Real longitude,
                   const Evaluate_Op evaluate_op = VALUE) const;

   };
*/

}

#endif /* DENISE_ASTRONOMY_H */ 
