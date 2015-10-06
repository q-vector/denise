//
// geodesy.h
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

#ifndef DENISE_GEODESY_H
#define DENISE_GEODESY_H

#include <cmath>
#include <cairomm/context.h>
#include <denise/analysis.h>
#include <denise/basics.h>
#include <denise/dtime.h>
#include <denise/exception.h>
#include <denise/geometry.h>
#include <denise/transform.h>
#include <denise/util.h>

using namespace std;
using namespace Cairo;

namespace denise
{

   class Lat_Long;
   class Geodesy;
   class Journey;
   class Geodetic_Mesh;
   class Geodetic_Transform;

   /// The FAI sphere radius
   const Real
   EARTH_RADIUS = 6371009;

   /// Length of one degree of latitude
   const Real
   LATITUDE_LENGTH = EARTH_RADIUS * DEGREE_TO_RADIAN;

   const Real
   EARTH_ROTATION = 7.2921159e-5;

   enum Lat_Long_Genre
   {
      LAT_LONG_STANDARD,
      LAT_LONG_PACIFIC
   };

   class Motion
   {

      public:

         Real
         u;

         Real
         v;

         Motion (const Real u = 0,
                 const Real v = 0);

         Motion (const Motion& motion);

         static Motion
         direction_speed (const Real direction,
                          const Real speed);

         virtual void
         set_from_direction_speed (const Real direction,
                                   const Real speed);
         virtual Real
         get_direction () const;

         Real
         get_speed () const;

         Motion
         operator+ (const Motion& motion) const;

         void
         operator+= (const Motion& motion);

         Motion
         operator- (const Motion& motion) const;

         void
         operator-= (const Motion& motion);

         Motion
         operator* (const Real a) const;

         void
         operator*= (const Real a);

         Motion
         operator/ (const Real a) const;

         void
         operator/= (const Real a);

   };

   /// Represents the Geodesy.
   ///
   /// This is an ellipsoidal geodesy model.  The calculations
   /// utilize Vincenty's (1975) iterative algorithm.  If flattening
   /// is 0 the geodesy is spherical and the calculations are based
   /// on direct methods of spherical trignometry.
   class Geodesy
   {

       private:

          Real
          a;

          Real
          f;

          Real
          epsilon_v;

          static Real
          asin (const Real cosine);

          static Real
          acos (const Real cosine);

          void
          spherical_inverse (Real& distance,
                             Real& azimuth_forward,
                             Real& azimuth_backward,
                             const Lat_Long& origin,
                             const Lat_Long& destination) const;

          void
          spherical_direct (Lat_Long& destination,
                            Real& azimuth_backward,
                            const Lat_Long& origin,
                            const Real distance,
                            const Real azimuth_forward) const;

          void
          vincenty_inverse (Real& distance,
                            Real& azimuth_forward,
                            Real& azimuth_backward,
                            const Lat_Long& origin,
                            const Lat_Long& destination) const;

          void
          vincenty_direct (Lat_Long& destination,
                           Real& azimuth_backward,
                           const Lat_Long& origin,
                           const Real distance,
                           const Real azimuth_forward) const;

       public:

          enum Model
          {
             AIRY_1930,
             MODIFIED_AIRY,
             AUSTRALIAN_NATIONAL,
             BESSEL_1841,
             BESSEL_NAMIBIA_1841,
             CLARKE_1866,
             CLARKE_1880,
             EVEREST_INDIA_1830,
             EVEREST_INDIA_1956,
             EVEREST_SABAH,
             EVEREST_MALAYSIA_1948,
             EVEREST_MALAYSIA_1969,
             MODIFIED_FISCHER_1960,
             GRS_1967,
             GRS_1980,
             HELMERT_1906,
             HOUGH_1960,
             INDONESIAN_1974,
             INTERNATIONAL_1924,
             KRASSOVSKY_1940,
             SOUTH_AMERICAN_1969,
             WGS60,
             WGS66,
             WGS72,
             WGS84,
             SPHERE
          };

          Geodesy (const Dstring& str);          

          Geodesy (const Model model,
                   const Real epsilon_v = 5e-12);

          Geodesy (const Real semi_major_axis = EARTH_RADIUS,
                   const Real flattening = 0,
                   const Real epsilon_v = 5e-12);

          void
          set (const Model model,
               const Real epsilon_v = 5e-12);

          Real
          get_semi_major_axis () const;

          Real
          get_semi_minor_axis () const;

          Real
          get_flattening () const;

          Real
          get_eccentricity () const;

          Real
          get_epsilon_v () const;

          void
          inverse (Real& distance,
                   Real& azimuth_forward,
                   Real& azimuth_backward,
                   const Lat_Long& origin,
                   const Lat_Long& destination) const;

          void
          direct (Lat_Long& destination,
                  Real& azimuth_backward,
                  const Lat_Long& origin,
                  const Real distance,
                  const Real azimuth_forward) const;

          static Real
          get_angle (const Lat_Long& lat_long_a,
                     const Lat_Long& lat_long_b,
                     const Geodesy& geodesy = Geodesy ());

          static Real
          get_distance (const Lat_Long& lat_long_a,
                        const Lat_Long& lat_long_b,
                        const Geodesy& geodesy = Geodesy ());

          static Lat_Long
          get_destination (const Lat_Long& origin,
                           const Real distance,
                           const Real azimuth_forward,
                           const Geodesy& geodesy = Geodesy ());

          static Real
          get_azimuth_forward (const Lat_Long& lat_long_a,
                               const Lat_Long& lat_long_b,
                               const Geodesy& geodesy = Geodesy ());

          static Real
          get_azimuth_backward (const Lat_Long& lat_long_a,
                                const Lat_Long& lat_long_b,
                                const Geodesy& geodesy = Geodesy ());

   };

   /// Represents a point on the Geodesy.
   ///
   /// Lat_Long consists of 2 Real's in (latitude, longitude)
   class Lat_Long
   {

      public:

         Real
         latitude;

         Real
         longitude;

         class List : public list<Lat_Long>
         {

            public:

               //void
               //add (const Journey::Simple& simple_journey,
               //     const Geodesy& geodesy = Geodesy (),
               //     const Real d_distance = 10e3,
               //     const Integer max_n = 1000);

               void
               add (const Journey& journey,
                    const Geodesy& geodesy = Geodesy (),
                    const Real d_distance = 10e3,
                    const Integer max_n = 1000);

         };

         Lat_Long (const Real latitude = 0.0,
                   const Real longitude = 0.0);

         Lat_Long (const Dstring& lat_long_string);
         
         Lat_Long (const Dstring& latitude_string,
                   const Dstring& longitude_string);
         
         Lat_Long (const Point_2D& point);

         static void
         standardize (Real& latitude,
                      Real& longitude,
                      const Lat_Long_Genre lat_long_genre = LAT_LONG_STANDARD);

         static void
         standardize (Real& latitude,
                      Real& longitude,
                      const Real standard_longitude);

         void
         standardize (const Lat_Long_Genre lat_long_genre = LAT_LONG_STANDARD);

         void
         standardize (const Real standard_longitude);

         Dstring
         get_string (const Integer decimal_places,
                     const bool nsew = false,
                     const bool with_symbol = false,
                     const bool with_parenthesis = false) const;

         Dstring
         get_string (const bool with_parenthesis = false,
                     const Dstring& number_format = "-%f") const;

         bool
         operator == (const Lat_Long& lat_long) const;
         
         bool
         operator != (const Lat_Long& lat_long) const;
         
         Lat_Long 
         operator = (const Lat_Long& lat_long);
         
         Lat_Long 
         operator + (const Lat_Long& lat_long) const;
         
         void
         operator += (const Lat_Long& lat_long);
         
         Lat_Long 
         operator - ();
         
         Lat_Long 
         operator - (const Lat_Long& lat_long) const;
         
         void
         operator -= (const Lat_Long& lat_long);
         
         Lat_Long 
         operator * (const Real scalar) const;

         void
         operator *= (const Real scalar);

         Lat_Long
         operator / (const Real scalar) const;

         void
         operator /= (const Real scalar);

         operator Point_2D () const;

         bool
         is_nall () const;

   };

/*
   class Location : public Lat_Long
   {

      protected:

         string
         identifier;

      public:

         Location (const string& identifier,
                   const Lat_Long& lat_long);

         Location (const Lat_Long& lat_long);

         Location (const Location& location);

         void
         set_identifier (const string& identifier);

         const string&
         get_identifier () const;

   };

   class Location_Set : public set<Location>
   {

      public:

         Location_Set::const_iterator
         nearest (const Lat_Long& lat_long) const;

   };
*/


   class Journey : public Simple_Polyline
   {

      public:

         /// Represents a Simple Journey on the Geodesy.
         ///
         /// Journey consists of the Lat_Long's of the origin and destination;
         /// the distance in between as well as the course back and fro between
         /// the 2 locations.
         class Simple : public Edge
         {

            protected:

               Real
               distance;

               Real
               azimuth_forward;

               Real
               azimuth_backward;

            public:

               Simple ();

               Simple (const Lat_Long& origin,
                       const Lat_Long& destination,
                       const bool complete = true);

               Simple (const Lat_Long& origin,
                       const Real distance,
                       const Real azimuth_forward,
                       const Geodesy& geodesy);

               Simple (const Lat_Long& origin,
                       const Real distance,
                       const Real azimuth_forward,
                       const bool complete = true);

               Simple (const Lat_Long& origin,
                       const Lat_Long& destination,
                       const Geodesy& geodesy);

               Simple (const Simple& simple_journey);

               void
               complete (const Geodesy& geodesy = Geodesy ());

               void
               standardize (const Lat_Long_Genre lat_long_genre = LAT_LONG_STANDARD);

               void
               standardize (const Real standard_longitude);

               static Dstring
               get_cardinal_direction (const Real azimuth);

               /// Swaps the origin and destination; as well as the
               /// azimuth_forward and azimuth_backward
               void
               swap ();

               virtual void
               translate (const Real distance,
                          const Real azimuth);

               void
               set_origin (const Lat_Long& origin);

               void
               set_destination (const Lat_Long& destination);

               void
               set_destination_azimuth (const Lat_Long& destination,
                                        const Real azimuth_backward);

               void
               set_distance_azimuth (const Real distance,
                                     const Real azimuth_forward,
                                     const Real azimuth_backward);

               Lat_Long
               get_origin () const;

               Lat_Long
               get_destination () const;

               const Real&
               get_distance () const;

               const Real&
               get_azimuth_forward () const;

               const Real&
               get_azimuth_backward () const;

               Lat_Long
               get_middle_lat_long (const Geodesy& geodesy = Geodesy ());

               Lat_Long
               get_lat_long (const Real fraction,
                             const Geodesy& geodesy = Geodesy ());

               void
               fill_lat_long_list (list<Lat_Long>& lat_long_list,
                                   const Real approx_d,
                                   const Integer max_n_per_leg,
                                   const bool first_leg,
                                   const Geodesy& geodesy = Geodesy ()) const;

               list<Lat_Long>*
               get_lat_long_list_ptr (const Real approx_d,
                                      const Integer max_n_per_leg,
                                      const Geodesy& geodesy = Geodesy ()) const;

         };

         Journey (const bool closed = false);

         Journey (const Simple_Polyline& simple_polyline);

         Journey (const Simple& simple_journey,
                  const Geodesy& geodesy = Geodesy ());

         Journey (const Journey& journey);

         Journey (const Dstring& str);

         static Journey
         fine_journey (const Simple& simple_journey,
                       const Geodesy& geodesy = Geodesy (),
                       const Real d_distance = 100e3,
                       const Integer max_n = 50);

         static Journey
         fine_journey (const Journey& journey,
                       const Geodesy& geodesy = Geodesy (),
                       const Real d_distance = 100e3,
                       const Integer max_n = 50);

         void
         standardize (const Lat_Long_Genre lat_long_genre = LAT_LONG_STANDARD);

         void
         standardize (const Real standard_longitude);

         virtual void
         translate (const Real distance,
                    const Real azimuth);

         Journey::Simple
         get_simple_journey (const Real x,
                             const Geodesy& geodesy = Geodesy ()) const;

         Journey::Simple
         get_simple_journey (Journey::iterator iterator) const;

         Journey::Simple
         get_simple_journey (Journey::const_iterator iterator) const;

         Real
         get_distance (const Geodesy& geodesy = Geodesy ()) const;

         Tuple
         get_tuple_x (const Geodesy& geodesy = Geodesy ()) const;

         Lat_Long
         get_lat_long (const Real x,
                       const Geodesy& geodesy = Geodesy ()) const;

         Lat_Long::List
         get_lat_long_list (const Geodesy& geodesy = Geodesy (),
                            const Real d_distance = 100e3,
                            const Integer max_n = 50) const;

         Real
         get_azimuth_forward (const Real x,
                              const Geodesy& geodesy = Geodesy ()) const;

         Real
         get_azimuth_forward (Journey::const_iterator iterator,
                              const Geodesy& geodesy = Geodesy ()) const;

         Real
         get_azimuth_forward (Journey::iterator iterator,
                              const Geodesy& geodesy = Geodesy ()) const;

         void
         cairo (const RefPtr<Context> cr,
                const Transform_2D& transform) const;

         Journey::iterator
         get_iterator (const Transform_2D& transform,
                       const Point_2D& point_2d,
                       const Geodesy& geodesy,
                       const Real threshold,
                       const Real standard_longitude = 0,
                       const Real d_distance = 100e3,
                       const Integer max_n = 50);

         Journey::const_iterator
         get_iterator (const Transform_2D& transform,
                       const Point_2D& point_2d,
                       const Geodesy& geodesy,
                       const Real threshold,
                       const Real standard_longitude = 0,
                       const Real d_distance = 100e3,
                       const Integer max_n = 50) const;

         Journey::iterator
         implant (const Transform_2D& transform,
                  const Point_2D& point_2d,
                  const Geodesy& geodesy,
                  const Real threshold,
                  const Real standard_longitude = 0,
                  const Real d_distance = 100e3,
                  const Integer max_n = 50);

   };

/*
   class Journey_List : public list<Journey>
   {

      public:

         Geodesy
         geodesy;

      public:

         Journey_List (const Lat_Long& lat_long);

         Journey_List (const Lat_Long& lat_long_f,
                       const Lat_Long& lat_long_t);

         void
         add_via_lat_long (Journey_List::iterator iterator,
                           const Lat_Long& lat_long);

         void
         set_origin (Journey_List::iterator iterator,
                     const Lat_Long& lat_long);

         void
         set_destination (Journey_List::iterator iterator,
                          const Lat_Long& lat_long);

         list<Lat_Long>*
         get_lat_long_list_ptr (const Real approx_d,
                                const Integer max_n_per_leg = 30) const;

   };
*/

   class Geodetic_Attractor
   {

      public:

         virtual pair<Dstring, Lat_Long>
         nearest (const Lat_Long& lat_long) const = 0;

   };

   class Degree_Geodetic_Attractor : public Geodetic_Attractor
   {

      private:

         Integer
         n;

      public:

         Degree_Geodetic_Attractor (const Integer n = 1);

         void
         set_n (const Integer n);

         const Integer
         get_n ();

         virtual pair<Dstring, Lat_Long>
         nearest (const Lat_Long& lat_long) const;

   };

   class Range_Circle : public Polygon
   {

      public:

         Range_Circle (const Lat_Long& lat_long,
                       const Real range,
                       const Real standard_longitude = 0,
                       const Integer n = 360,
                       const Geodesy& geodesy = Geodesy ());

   };   

   class Range_Sector : public Polygon
   {

      public:

         Range_Sector (const Lat_Long& lat_long,
                       const Real range,
                       const Real start_azimuth,
                       const Real end_azimuth,
                       const Real standard_longitude = 0,
                       const Integer n = 360,
                       const Geodesy& geodesy = Geodesy ());

   };   

   class Geodetic_Transform : public Transform_2D
   {

      public:

         enum Genre
         { 
            MERCATOR,
            LAMBERT_CONIC_NORTH,
            LAMBERT_CONIC_SOUTH,
            EQUIDISTANT_CYLINDRICAL,
            POLAR_STEREOGRAPHIC_NORTH,
            POLAR_STEREOGRAPHIC_SOUTH,
            PERSPECTIVE,
            GEOS,
            MOLLWEIDE,
            UNDEFINED
         };

         class Data
         {

            public:

               Genre
               genre;

               Real
               scale;

               Lat_Long
               lat_long;

               Data (const Genre genre,
                     const Real scale,
                     const Lat_Long& lat_long);

               Data (const Dstring& str);

               Data (const Data& data);

               void
               standardize (Lat_Long& lat_long) const;

               Dstring
               get_string () const;

         };

         Data
         data;

         Geodetic_Transform (const Geodetic_Transform::Genre genre,
                             const Real scale,
                             const Lat_Long& lat_long);

         Geodetic_Transform (const Dstring& str);

         Geodetic_Transform (const Geodetic_Transform& geodetic_transform);

         static bool
         is_geodetic (const Dstring& str);

         virtual Geodetic_Transform*
         clone () const = 0;

         static Geodetic_Transform*
         get_transform_ptr (const Dstring& str,
                            const Point_2D& point);

         static Geodetic_Transform*
         get_transform_ptr (const Geodetic_Transform::Genre genre,
                            const Real scale,
                            const Lat_Long& lat_long,
                            const Point_2D& point);

         static Geodetic_Transform*
         get_transform_ptr (const Geodetic_Transform::Data& data,
                            const Point_2D& point);

         virtual bool
         is_out_of_domain (const Lat_Long& lat_long) const;

         Real
         get_scale () const;

         virtual void
         reverse (Real& latitude,
                  Real& longitude,
                  const Real x,
                  const Real y) const = 0;

         void
         reverse (Lat_Long& lat_long,
                  const Point_2D& point) const;

         const Lat_Long&
         get_lat_long () const;

         Lat_Long&
         get_lat_long ();

         Lat_Long
         get_lat_long (const Point_2D& point) const;

         void
         standardize (Lat_Long& lat_long) const;

         Domain_2D
         get_domain_2d (const Size_2D& size_2d) const;

         Geodetic_Mesh
         get_mesh (const Size_2D& size_2d) const;

   };

   class Equidistant_Cylindrical_Transform : public Geodetic_Transform
   {

      private:

         Affine_Transform_2D
         affine_transform;

      public:

         Equidistant_Cylindrical_Transform (const Domain_1D& domain_latitude,
                                            const Domain_1D& domain_longitude,
                                            const Real width,
                                            const Real height,
                                            const Point_2D& anchor = Point_2D (0, 0));

         virtual Geodetic_Transform*
         clone () const;

         virtual bool
         out_of_domain (const Real latitude,
                        const Real longitude) const;

         void
         transform (Real& x,
                    Real& y,
                    const Real latitude,
                    const Real longitude) const;

         void
         reverse (Real& latitude,
                  Real& longitude,
                  const Real x,
                  const Real y) const;
   
   }; 

   class Lambert_Conic_Transform : public Geodetic_Transform
   {

      private:

         Point_2D
         ref_point;
      
         Lat_Long
         ref_lat_long;
      
         Real
         cone_constant;

         Real
         sin_colatitude_1;

         Real
         tan_half_colatitude_1;

         Real
         r_ref_lat;

         Real
         get_cone_constant (const Real latitude_1,
                            const Real latitude_2) const;

         Real
         get_r (const Real latitude) const;

         Real
         get_latitude (const Real r) const;

      public:

         Lambert_Conic_Transform (const Real scale,
                                  const Lat_Long& ref_lat_long,
                                  const Point_2D& ref_point,
                                  const bool southern_hemisphere = false);

         Lambert_Conic_Transform (const Real scale,
                                  const Lat_Long& ref_lat_long,
                                  const Point_2D& ref_point,
                                  const Real true_latitude_1 = 30,
                                  const Real true_latitude_2 = 60);

         Lambert_Conic_Transform (const Lambert_Conic_Transform& transform);

         virtual Geodetic_Transform*
         clone () const;

         virtual bool
         out_of_domain (const Real latitude,
                        const Real longitude) const;

         void
         transform (Real& x,
                    Real& y,
                    const Real latitude,
                    const Real longitude) const;

         void
         reverse (Real& latitude,
                  Real& longitude,
                  const Real x,
                  const Real y) const;

         void
         transform_uv (Real& u,
                       Real& v,
                       const Real latitude,
                       const Real longitude) const;

         Real
         get_theta (const Real u,
                    const Real v,
                    const Real latitude,
                    const Real longitude) const;

   };

   class Mercator_Transform : public Geodetic_Transform
   {

      private:

         Real
         cos_true_latitude;

         Lat_Long
         ref_lat_long;

         Real
         offset_x;

         Real
         offset_y;

         Real
         a_e;

         Real
         project_x (const Real longitude) const;

         Real
         project_y (const Real latitude) const;

         Real
         reverse_latitude (const Real y) const;

         Real
         reverse_longitude (const Real x) const;

      public:

         Mercator_Transform (const Real scale,
                             const Lat_Long& ref_lat_long,
                             const Point_2D& ref_point,
                             const Real true_latitude = 0);

         virtual Geodetic_Transform*
         clone () const;

         virtual bool
         out_of_domain (const Real latitude,
                        const Real longitude) const;

         void
         transform (Real& x,
                    Real& y,
                    const Real latitude,
                    const Real longitude) const;

         void
         reverse (Real& latitude,
                  Real& longitude,
                  const Real x,
                  const Real y) const;

         Real
         get_theta (const Real u,
                    const Real v,
                    const Real latitude = GSL_NAN,
                    const Real longitude = GSL_NAN) const;

   };

   class Polar_Stereographic_Transform : public Geodetic_Transform
   {
   
      private:

         bool
         southern_hemisphere;
      
         Real
         ref_longitude;
         
         Real
         m_0;

         Point_2D
         ref_point;
         
         Real
         get_r (const Real latitude) const;
         
         Real
         get_latitude (const Real r) const;
         
      public:
      
         Polar_Stereographic_Transform (const Real scale,
                                        const Real ref_longitude,
                                        const Point_2D& ref_point,
                                        const bool southern_hemisphere = false);
                            
         virtual Geodetic_Transform*
         clone () const;

         virtual bool
         out_of_domain (const Real latitude,
                        const Real longitude) const;

         void
         transform (Real& x,
                    Real& y,
                    const Real latitude,
                    const Real longitude) const;

         void
         reverse (Real& latitude,
                  Real& longitude,
                  const Real x,
                  const Real y) const;

         void
         transform_uv (Real& u,
                       Real& v,
                       const Real latitude,
                       const Real longitude) const;

         Real
         get_theta (const Real u,
                    const Real v,
                    const Real latitude,
                    const Real longitude) const;

   };

   class Perspective_Transform : public Geodetic_Transform
   {

      private:
 
         Equidistant_Cylindrical_Transform
         ttt;

         const Geodesy
         geodesy;

         const Range_Circle
         range_circle;

         const Lat_Long
         ref_lat_long;

         const Point_2D
         ref_point;

         const Real
         height;

      public:

         Perspective_Transform (const Real scale,
                                const Lat_Long& ref_lat_long,
                                const Point_2D& ref_point,
                                const Real height);

         ~Perspective_Transform ();

         virtual Geodetic_Transform*
         clone () const;

         const Range_Circle&
         get_domain () const;

         virtual bool
         out_of_domain (const Real latitude,
                        const Real longitude) const;

         void
         transform (Real& transformed_x,
                    Real& transformed_y,
                    const Real x, 
                    const Real y) const;
         
         void
         reverse (Real& reversed_x,
                  Real& reversed_y,
                  const Real x,
                  const Real y) const;

         void
         cairo (const RefPtr<Context> cr,
                const Polygon& polygon) const;

         void
         cairo (const RefPtr<Context> cr,
                const Simple_Polyline& simple_polyline) const;

   };

   class Geos_Transform : public Geodetic_Transform
   {

      private:

         const Real
         nadir_longitude;

         const Real
         coff;

         const Real
         loff;

         const Real
         cfac;

         const Real
         lfac;

      public:

         Geos_Transform (const Real nadir_longitude,
                         const Real coff,
                         const Real loff,
                         const Real cfac,
                         const Real lfac);

         virtual Geodetic_Transform*
         clone () const;

         virtual bool
         out_of_domain (const Real latitude,
                        const Real longitude) const;

         void
         transform (Real& transformed_x,
                    Real& transformed_y,
                    const Real x, 
                    const Real y) const;
         
         void
         reverse (Real& reversed_x,
                  Real& reversed_y,
                  const Real x,
                  const Real y) const;

   };

   class Mollweide_Transform : public Geodetic_Transform
   {

      private:

         const Real
         epsilon;

         const Real
         scale;

         const Real
         k;

         const Point_2D
         anchor;

         const Lat_Long
         lat_long;

         Real
         get_theta (const Real phi) const;

         Real
         get_phi (const Real theta) const;

      public:

         Mollweide_Transform (const Real scale,
                              const Point_2D& anchor,
                              const Lat_Long& lat_long = Lat_Long (0, 0));

         virtual Geodetic_Transform*
         clone () const;

         virtual bool
         out_of_domain (const Real latitude,
                        const Real longitude) const;

         void
         transform (Real& transformed_x,
                    Real& transformed_y,
                    const Real x, 
                    const Real y) const;
         
         void
         reverse (Real& reversed_x,
                  Real& reversed_y,
                  const Real x,
                  const Real y) const;

   };

   class Geodetic_Vector_Data_2D : public virtual Vector_Data_2D
   {

      protected:

         virtual Real
         get_dmagnitude_dx (const Integer vector_element_u,
                            const Integer vector_element_v,
                            const Integer node_x,
                            const Integer node_y,
                            const Real magnitude) const;

         virtual Real
         get_dmagnitude_dy (const Integer vector_element_u,
                            const Integer vector_element_v,
                            const Integer node_x,
                            const Integer node_y,
                            const Real magnitude) const;

         virtual Real
         get_dmagnitude_dx (const Integer vector_element_u,
                            const Integer vector_element_v,
                            const Real x,
                            const Real y,
                            const Real magnitude) const;

         virtual Real
         get_dmagnitude_dy (const Integer vector_element_u,
                            const Integer vector_element_v,
                            const Real x,
                            const Real y,
                            const Real magnitude) const;

      public:

         Geodetic_Vector_Data_2D (const Integer vector_size,
                                  const Size_2D& size_2d,
                                  const Domain_2D& domain_2d,
                                  const bool periodic_longitude = false);

         Geodetic_Vector_Data_2D (const Integer vector_size,
                                  const Tuple tuple_latitude,
                                  const Tuple tuple_longitude,
                                  const bool periodic_longitude = false);

         virtual
         ~Geodetic_Vector_Data_2D ();

         static Real
         get_f (const Real latitude);

         static Real
         get_f_hat (const Real latitude);

         Real
         get_latitude (const Integer i) const;

         Real
         get_longitude (const Integer j) const;

         Real
         evaluate (const Integer vector_element,
                   const Integer i,
                   const Integer j,
                   const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate (const Integer vector_element,
                   const Real latitude,
                   const Real longitude,
                   const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate (const Integer vector_element,
                   const Lat_Long& lat_long,
                   const Evaluate_Op evaluate_op = VALUE) const;

         Real
         get_shear_vorticity (const Integer vector_element_u,
                              const Integer vector_element_v,
                              const Integer i,
                              const Integer j) const;

         void
         subtract_zonal_mean (const Integer vector_element);

         void
         subtract_meridional_mean (const Integer vector_element);

   };

   class Geodetic_Vector_Data_3D : public virtual Vector_Data_3D
   {

      public:

         Geodetic_Vector_Data_3D (const Integer vector_size,
                                  const Size_3D& size_3d,
                                  const Domain_3D& domain_3d,
                                  const bool periodic_longitude = false);

         Geodetic_Vector_Data_3D (const Integer vector_size,
                                  const Tuple tuple_z,
                                  const Tuple tuple_latitude,
                                  const Tuple tuple_longitude,
                                  const bool periodic_longitude = false);

         Geodetic_Vector_Data_3D (const Integer vector_size,
                                  const Tuple tuple_z,
                                  const Size_2D& size_2d,
                                  const Domain_2D& domain_2d,
                                  const bool periodic_longitude = false);

         virtual
         ~Geodetic_Vector_Data_3D ();

         static Real
         get_f (const Real latitude);

         static Real
         get_f_hat (const Real latitude);

         Real
         get_latitude (const Integer i) const;

         Real
         get_longitude (const Integer j) const;

         Real
         evaluate (const Integer vector_element,
                   const Real z,
                   const Integer i,
                   const Integer j,
                   const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate (const Integer vector_element,
                   const Integer k,
                   const Integer i,
                   const Integer j,
                   const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate (const Integer vector_element,
                   const Real z,
                   const Real latitude,
                   const Real longitude,
                   const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate (const Integer vector_element,
                   const Integer k,
                   const Real latitude,
                   const Real longitude,
                   const Evaluate_Op evaluate_op = VALUE) const;

         void
         subtract_zonal_mean (const Integer vector_element);

         void
         subtract_meridional_mean (const Integer vector_element);

         void
         subtract_horizontal_mean (const Integer vector_element);

   };

   class Geodetic_Scalar_Data_2D : public Scalar_Data_2D,
                                   public Geodetic_Vector_Data_2D
   {

      public:

         Geodetic_Scalar_Data_2D (const Size_2D& size_2d,
                                  const Domain_2D& domain_2d,
                                  const bool periodic_longitude = false);

         Geodetic_Scalar_Data_2D (const Tuple tuple_latitude,
                                  const Tuple tuple_longitude,
                                  const bool periodic_longitude = false);

         ~Geodetic_Scalar_Data_2D ();

         void
         set_datum (const Integer node_latitude,
                    const Integer node_longitude,
                    const Real datum);

         const Real&
         get_datum (const Integer node_latitude,
                    const Integer node_longitude) const;

         Real
         evaluate (const Real latitude,
                   const Real longitude,
                   const Evaluate_Op evaluate_op = VALUE) const;

   };

   class Geodetic_Scalar_Data_3D : public Scalar_Data_3D,
                                   public Geodetic_Vector_Data_3D
   {

      public:

         Geodetic_Scalar_Data_3D (const Size_3D& size_3d,
                                  const Domain_3D& domain_3d,
                                  const bool periodic_longitude = false);

         Geodetic_Scalar_Data_3D (const Tuple tuple_z,
                                  const Tuple tuple_latitude,
                                  const Tuple tuple_longitude,
                                  const bool periodic_longitude = false);

         Geodetic_Scalar_Data_3D (const Tuple tuple_z,
                                  const Size_2D& size_2d,
                                  const Domain_2D& domain_2d,
                                  const bool periodic_longitude = false);

         ~Geodetic_Scalar_Data_3D ();

         void
         set_datum (const Integer node_z,
                    const Integer node_latitude,
                    const Integer node_longitude,
                    const Real datum);

         const Real&
         get_datum (const Integer node_z,
                    const Integer node_latitude, 
                    const Integer node_longitude) const;

         Real
         evaluate (const Real z,
                   const Real latitude,
                   const Real longitude,
                   const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate (const Integer k,
                   const Real latitude,
                   const Real longitude,
                   const Evaluate_Op evaluate_op = VALUE) const;

   };

   class Track_Data : public map<Real, Real>
   {

      private:

         Scalar_Data_1D*
         spline_ptr;

      public:

         Track_Data ();

         Track_Data (const Track_Data& track_data);

         ~Track_Data ();

         void
         okay (const bool cubic = false);

         set<Real>
         get_tau_set () const;

         Real
         get_datum (const Real tau,
                    const bool forbid_extrapolate) const;

         Real
         get_dt (const Real tau,
                 const bool forbid_extrapolate) const;

         Domain_1D
         get_domain_1d (const Real dt = 0.1) const;

         friend ostream&
         operator<< (ostream &o,
                     const Track_Data& track_data);

   };

   ostream&
   operator<< (ostream &out,
               const Track_Data& track_data);

   class Track : public map<Dstring, Track_Data>
   {

      protected:

         Dtime
         epoch;

         set<Dstring>
         element_set;

         Tokens
         cubic_tokens;

      public:

         Track (const Dtime& epoch = Dtime (0.0));

         ~Track ();

         void
         set_epoch (const Dtime& epoch);

         const Dtime&
         get_epoch () const;

         Dtime
         get_dtime (const Real tau) const;

         Dtime
         get_start_time () const;

         Dtime
         get_end_time () const;

         Real
         get_start_tau () const;

         Real
         get_end_tau () const;

         const set<Dstring>&
         get_element_set () const;

         set<Real>
         get_tau_set (const Dstring& element = "latitude") const;

         set<Dtime>
         get_dtime_set (const Dstring& element = "latitude") const;

         void
         add (const Real tau,
              const Lat_Long& lat_long);

         void
         add (const Dstring& element,
              const Real tau,
              const Real datum);

         void
         add (const Dtime& dtime,
              const Lat_Long& lat_long);

         void
         add (const Dstring& element,
              const Dtime& dtime,
              const Real datum);

         void
         okay ();

         void
         okay (const Dstring& element);

         Lat_Long
         get_lat_long (const Dtime& dtime,
                       const bool forbid_extrapolate = true) const;

         Lat_Long
         get_lat_long (const Real tau,
                       const bool forbid_extrapolate = true) const;

         Motion
         get_motion (const Dtime& dtime,
                     const bool forbid_extrapolate = true) const;

         Motion
         get_motion (const Real tau,
                     const bool forbid_extrapolate = true) const;

         Real
         get_datum (const Dstring& element,
                    const Dtime& dtime,
                    const bool forbid_extrapolate = true) const;

         Real
         get_datum (const Dstring& element,
                    const Real tau,
                    const bool forbid_extrapolate = true) const;

         Domain_2D
         get_domain_2d (const Real dt = 0.1) const;

         bool
         trespass (const Domain_2D& domain_2d,
                   const Real dt = 0.1) const; 

         void
         write (ostream& o,
                const Dstring& id) const;

         friend ostream&
         operator<< (ostream& o,
                     const Track& track);

         class Map : public map<Dstring, Track>
         {

            public:

               Map ();

               Map (igzstream& i);

               void
               ingest (igzstream& i);

               void
               write (ostream& o) const;

         };

   };

   ostream&
   operator<< (ostream &out,
               const Track& track);

   class Geodetic_Cairoable
   {

      public:

         virtual void
         cairo (const RefPtr<Context> cr,
                const Geodetic_Transform& transform) const = 0;

   };

   class Geodetic_Mesh : public Mesh_2D,
                         public Geodetic_Cairoable
   {

      public:

         Geodetic_Mesh (const Size_2D& s = Size_2D (100, 100),
                        const Domain_2D& d = Domain_2D (-89, 89, -180, 179.9));

         Geodetic_Mesh (const Real interval_x,
                        const Real interval_y,
                        const Color& color,
                        const Size_2D& s = Size_2D (100, 100),
                        const Domain_2D& d = Domain_2D (-89, 89, -180, 179.9));

         Geodetic_Mesh (const Real interval,
                        const Color& color,
                        const Size_2D& s = Size_2D (100, 100),
                        const Domain_2D& d = Domain_2D (-89, 89, -180, 179.9));

         Geodetic_Mesh (const Real interval_x_0,
                        const Real interval_y_0,
                        const Color& color_0,
                        const Real interval_x_1,
                        const Real interval_y_1,
                        const Color& color_1,
                        const Size_2D& s = Size_2D (100, 100),
                        const Domain_2D& d = Domain_2D (-89, 89, -180, 179.9));

         Geodetic_Mesh (const Real interval_0,
                        const Color& color_0,
                        const Real interval_1,
                        const Color& color_1,
                        const Size_2D& s = Size_2D (100, 100),
                        const Domain_2D& d = Domain_2D (-89, 89, -180, 179.9));

         Geodetic_Mesh (const Real interval_x_0,
                        const Real interval_y_0,
                        const Color& color_0,
                        const Real interval_x_1,
                        const Real interval_y_1,
                        const Color& color_1,
                        const Real interval_x_2,
                        const Real interval_y_2,
                        const Color& color_2,
                        const Size_2D& s = Size_2D (100, 100),
                        const Domain_2D& d = Domain_2D (-89, 89, -180, 179.9));

         Geodetic_Mesh (const Real interval_0,
                        const Color& color_0,
                        const Real interval_1,
                        const Color& color_1,
                        const Real interval_2,
                        const Color& color_2,
                        const Size_2D& s = Size_2D (100, 100),
                        const Domain_2D& d = Domain_2D (-89, 89, -180, 179.9));

         Geodetic_Mesh (const Real interval_x_0,
                        const Real interval_y_0,
                        const Color& color_0,
                        const Real interval_x_1,
                        const Real interval_y_1,
                        const Color& color_1,
                        const Real interval_x_2,
                        const Real interval_y_2,
                        const Color& color_2,
                        const Real interval_x_3,
                        const Real interval_y_3,
                        const Color& color_3,
                        const Size_2D& s = Size_2D (100, 100),
                        const Domain_2D& d = Domain_2D (-89, 89, -180, 179.9));

         Geodetic_Mesh (const Real interval_0,
                        const Color& color_0,
                        const Real interval_1,
                        const Color& color_1,
                        const Real interval_2,
                        const Color& color_2,
                        const Real interval_3,
                        const Color& color_3,
                        const Size_2D& s = Size_2D (100, 100),
                        const Domain_2D& d = Domain_2D (-89, 89, -180, 179.9));

         void
         cairo (const RefPtr<Context> cr,
                const Geodetic_Transform& transform) const;

   };

   ostream&
   operator << (ostream &out,
                const Lat_Long& lat_long);

}

#endif /* DENISE_GEODESY_H */ 
