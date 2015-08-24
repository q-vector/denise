//
// tc.h
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

#ifndef DENISE_TC_H
#define DENISE_TC_H

#include <denise/analysis.h>
#include <denise/cluster.h>
#include <denise/graphics.h>
#include <denise/geodesy.h>
#include <denise/dtime.h>
#include <denise/stat.h>

using namespace std;
using namespace Cairo;
using namespace denise;

namespace denise
{

   class Tc_Track : public Track
   {

      protected:

         const Dstring
         name;

      public:

         Tc_Track (const Dstring& name = "",
                   const Dtime& dtime = Dtime (0.0));

         const Dstring&
         get_name () const;

   };

   class Best_Tracks : public map<Dstring, Tc_Track>
   {

      public:

         class Id_Set : public set<Dstring>
         {

            private:

               const Best_Tracks&
               best_tracks;

            public:

               Id_Set (const Id_Set& id_set);

               Id_Set (const Best_Tracks& best_tracks,
                       const bool fill = false);

               Id_Set
               get_id_set (const Integer year) const;

               Id_Set
               get_id_set (const Dstring& name) const;

               Id_Set
               get_id_set (const Domain_2D& domain_2d,
                           const Real dt = 0.1) const;

               Id_Set
               get_id_set (const Integer day_of_year,
                           const Integer delta_days,
                           const Domain_2D& domain_2d,
                           const Real dt = 0.1) const;

         };

      private:

         Id_Set
         id_set;

      public:

         Best_Tracks ();

         void
         ingest_jma (const Dstring& file_path);

         Id_Set
         get_id_set (const Integer year) const;

         Id_Set
         get_id_set (const Dstring& name) const;

         Id_Set
         get_id_set (const Domain_2D& domain_2d,
                     const Real dt = 0.1) const;

         Id_Set
         get_id_set (const Integer day_of_year,
                     const Integer delta_days,
                     const Domain_2D& domain_2d,
                     const Real dt = 0.1) const;

   };

   class Forecast
   {

      private:

         static Real
         get_pressure (const Dstring& pressure_string);

         static Real
         get_max_wind (const Dstring& max_wind_string);

      public:

         Lat_Long
         lat_long;

         Real
         pressure;

         Real
         max_wind;

         Forecast () { }

         Forecast (const Lat_Long& lat_long,
                   const Real pressure,
                   const Real max_wind);

         Forecast (const Lat_Long& lat_long,
                   const Dstring& pressure_string,
                   const Dstring& max_wind_string);

   };

   class Advisory : public map<Integer, Forecast>
   {

      protected:

         Track
         track;

         Dtime
         get_time (const Dstring& time_string,
                   const Dtime& time_stamp) const;

      public:

         Dstring
         forecast_centre;

         Dstring
         icon_string;

         Dstring
         tc_name;

         Dstring
         tc_id;

         Dtime
         initial_time;

         Dstring
         initial_time_string;

         Advisory (const Dstring& forecast_centre);

         Advisory (const Dstring& forecast_centre,
                   const Dstring& icon_string);

         void
         make_track ();

         Dstring
         get_key () const;

         bool
         operator == (const Advisory& advisory) const;

         bool
         operator > (const Advisory& advisory) const;

         bool
         operator < (const Advisory& advisory) const;

         Lat_Long
         get_lat_long (const Real tau) const;

         void
         cairo (const RefPtr<Context> cr,
                const Geodetic_Transform& transform,
                const Real intensity = 1.0) const;

         ostream&
         operator << (ostream& out_file) const;

   };

   class Wtss20_Vhhh : public Advisory
   {

      public:

         Wtss20_Vhhh (const Tokens& content,
                      const Dtime& time_stamp);

   };

   class Wtpq20_Babj : public Advisory
   {

      public:

         Wtpq20_Babj (const Tokens& content,
                      const Dtime& time_stamp);

   };

   class Wtpq20_Rjtd : public Advisory
   {

      public:

         Wtpq20_Rjtd (const Tokens& content,
                      const Dtime& time_stamp);

   };

   class Wtko20_Rksl : public Advisory
   {

      public:

         Wtko20_Rksl (const Tokens& content,
                      const Dtime& time_stamp);

   };

   class Wtpn31_Pgtw : public Advisory
   {

      public:

         Wtpn31_Pgtw (const Tokens& content,
                      const Dtime& time_stamp);

   };

   class Wtph20_Rpmm : public Advisory
   {

      private:

         static Real
         get_real (const Dstring& digit_string);

         static Lat_Long
         get_lat_long (const Dstring& lat_long_string);

         static Real
         get_central_pressure (const Tokens& tokens);

         static Real
         get_meters_per_second (const Tokens& tokens);

      public:

         Wtph20_Rpmm (const Tokens& content,
                      const Dtime& time_stamp);

   };

   class Wtnt80_Egrr : public Advisory
   {

      public:

         Wtnt80_Egrr (const Tokens& content,
                      const Dtime& time_stamp);

   };

   class Wtth20_Vtbb : public Advisory
   {

      public:

         Wtth20_Vtbb (const Tokens& content,
                      const Dtime& time_stamp);

   };

   class Knhc_4 : public Advisory
   {

      public:

         Knhc_4 (const Tokens& content,
                 const Dtime& time_stamp);

   };

   class Advisory_Store : public map<Dstring, Advisory>
   {

      public:

         class Cluster_Info
         {

            private:

               const Advisory_Store&
               advisory_store;

               const Tokens
               key_tokens;

               Cluster*
               cluster_ptr;

               Cluster::Multimap
               cluster_multimap;

            public:

               Cluster_Info (const Advisory_Store& advisory_store,
                             const Tokens& key_tokens,
                             const Dataset& dataset,
                             const Real cluster_distance = 2);

               ~Cluster_Info ();

               set<Integer>
               get_cluster_index_set () const;

               Dstring
               get_tc_name (const Integer cluster_index) const;

               Lat_Long
               get_lat_long (const Integer cluster_index) const;

         };

      private:

         void
         insert (const Advisory& advisory);

         void
         parse (const Tokens& content,
                const Dtime& time_stamp);

      public:

         Advisory_Store ();

         void
         ingest_dir (const Dstring& dir_path,
                     const Reg_Exp& file_reg_exp);

         void
         ingest_file (const Dstring& file_path);

         set<Dtime>
         get_initial_time_set (const Area* area_ptr = NULL,
                               const bool synoptic_hours_only = true,
                               const bool sub_synoptic_okay = false) const;

         Tokens
         get_key_tokens (const Area* area_ptr = NULL) const;

         Tokens
         get_key_tokens (const Dstring& initial_time_string,
                         const Area* area_ptr = NULL) const;

         Tokens
         get_key_tokens (const set<Dstring>& initial_time_string_set,
                         const Area* area_ptr = NULL) const;

         Tokens
         get_key_tokens (const Dtime& initial_time,
                         const Area* area_ptr = NULL) const;

         Tokens
         get_key_tokens (const set<Dtime>& initial_time_set,
                         const Area* area_ptr = NULL) const;

         Cluster_Info
         get_cluster_info (const Tokens& key_tokens,
                           const Real criteria_distance = 2) const;

         Cluster::Multimap
         get_cluster_multimap (const Tokens& key_tokens,
                               const Real criteria_distance) const;

         vector<Point_2D>
         get_position_vector (const Tokens& key_tokens,
                              const Real tau) const;

         Bivariate_Gaussian_Distribution
         get_bivariate_gaussian_distribution (const Tokens& key_tokens,
                                              const Real tau) const;

         Ellipse
         get_ellipse (const Tokens& key_tokens,
                      const Real tau,
                      const Real probability) const;

   };

};

#endif /* DENISE_TC_H */
