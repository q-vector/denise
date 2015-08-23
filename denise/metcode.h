//
// metcode.h
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

#ifndef DENISE_METCODE_H
#define DENISE_METCODE_H

#include <denise/dstring.h>
#include <denise/exception.h>
#include <denise/geodesy.h>
#include <denise/met.h>

using namespace std;

#define METAR_TIME_TOLERANCE 0.01

namespace denise
{

   class Location
   {

      public:

         Dstring
         id;

         Dstring
         name;

         Lat_Long
         lat_long;

         Real
         height;

         Location ();

         Location (const Dstring& str,
                   const Dstring& separator,
                   const Integer id_index,
                   const Integer name_index,
                   const Integer latitude_index,
                   const Integer longitude_index);

         Location (const Dstring& str,
                   const Dstring& separator,
                   const Integer id_index,
                   const Integer name_index,
                   const Integer latitude_index,
                   const Integer longitude_index,
                   const Integer height_index);

         Location (const Dstring& id,
                   const Dstring& name,
                   const Lat_Long& lat_long);

         Location (const Dstring& id,
                   const Dstring& name,
                   const Lat_Long& lat_long,
                   const Real height);

         bool
         operator == (const Location& location) const;

         bool
         operator < (const Location& location) const;

   };

   class Location_Map : public map<Dstring, Location>
   {

      public:

         const Location&
         get_nearest (const Lat_Long& lat_long) const;

   };

   class Location_Multimap : public multimap<Dstring, Location>
   {

      public:

         const Location&
         get_nearest (const Lat_Long& lat_long) const;

   };

   class Metars;
   class Station_Metars;

   class Metar_Wind : public Wind
   {

      private:

         bool
         variable_direction;

         Real
         gust;

      public:

         Metar_Wind ();

         Metar_Wind (const Dstring& metar_wind_string);

         Real
         get_gust () const;

   };

   class Metar
   {

      private:

         class Key
         {

            public:

               Dtime
               time;

               Dstring
               station_name;

               Key ();

               Key (const Dtime& time,
                    const Dstring& station_name);

               bool
               operator == (const Key& key) const;

               bool
               operator > (const Key& key) const;

               bool
               operator < (const Key& key) const;

         };

         Key
         key;

         bool
         speci;

         Metar_Wind
         metar_wind;

         Real
         temperature;

         Real
         dew_point;

         Real
         qnh;

         Real
         visibility;

         Dstring
         cloud;

         Dstring
         weather;

         Real
         rain_10_minute;

         Real
         rain_hour;

         Real
         rain_since_9am;

         Real
         visibility_auto;

         Dstring
         cloud_auto;

         void
         read_t_td (const Dstring& token);

         void
         read_qnh (const Dstring& token);

         void
         read_rainfall (const Dstring& token);

      public:

         Metar (const Dstring& metar_string);

         Dstring
         get_string () const;

         const Dtime&
         get_time () const;

         const Dstring&
         get_station_name () const;

         const Real&
         get_temperature () const;

         const Real&
         get_dew_point () const;

         const Real
         get_rh () const;

         const Metar_Wind&
         get_metar_wind () const;

         const Real&
         get_qnh () const;

         const Real&
         get_rain_10_minute () const;

         const Real&
         get_rain_hour () const;

         const Real&
         get_rain_since_9am () const;

         const Real
         get_visibility () const;

         const Real
         get_gfdi () const;

         const Real
         get_ffdi () const;

         bool
         operator == (const Metar& metar) const;

         bool
         operator > (const Metar& metar) const;

         bool
         operator < (const Metar& metar) const;

   };

   class Metars : public set<Metar>,
                  public Attractor
   {

      private:

         Dstring
         dir_path;

         Dstring
         file_format;

         map<Dstring, Lat_Long>
         station_name_lat_long_map;

         set<Dtime>
         time_set;

         map<Dstring, Station_Metars*>
         station_metars_ptr_map;

         void
         clear_station_metars_ptr_map ();

      public:

         Metars ();

         Metars (const Dstring& dir_path,
                 const Dstring& file_format);

         Metars (const Dstring& dir_path,
                 const Dstring& file_format,
                 const Location_Map& location_map);

         Metars (const Dstring& dir_path,
                 const Dstring& file_format,
                 const Location_Multimap& location_multimap);

         ~Metars ();

         void
         setup (const Dstring& dir_path,
                const Dstring& file_format);

         Dstring
         get_status () const;

         void
         reload ();

         void
         reload (const Location_Map& location_map);

         void
         reload (const Location_Multimap& location_multimap);

         void
         read (const Dstring& dir_path,
               const Dstring& file_format);

         void
         read (const Dstring& dir_path,
               const Dstring& file_format,
               const Location_Map& location_map);

         void
         read (const Dstring& dir_path,
               const Dstring& file_format,
               const Location_Multimap& location_multimap);

         void
         read (const Dstring& file_path);

         void
         read (const Dstring& file_path,
               const Location_Map& location_map);

         void
         read (const Dstring& file_path,
               const Location_Multimap& location_multimap);

         void
         construct_station_metars_ptr_map ();

         const map<Dstring, Lat_Long>&
         get_station_name_lat_long_map () const;

         void
         attract (Real& latitude,
                  Real& longitude) const;

         map<Dstring, Lat_Long>::const_iterator
         get_nearest_station_name_lat_long (const Lat_Long& lat_long) const;

         const Metar&
         get_nearest_metar (const Dtime& time,
                            const Lat_Long& lat_long) const;

         const Lat_Long&
         get_lat_long (const Dstring& station_name) const;

         const set<Dtime>&
         get_time_set () const;

         const Station_Metars&
         get_station_metars (const Dstring& station_name) const;

         Station_Metars&
         get_station_metars (const Dstring& station_name);

         const Station_Metars&
         get_station_metars (const Lat_Long& lat_long) const;

         Station_Metars&
         get_station_metars (const Lat_Long& lat_long);

         const Metar&
         get_metar (const Dtime& time,
                    const Dstring& station_name) const;

         vector<const Metar*>
         get_metar_ptr_vector (const Dtime& time) const;

         pair<Dstring, Lat_Long>
         nearest (const Lat_Long& lat_long) const;

   };

   class Station_Metars : public map<Dtime, const Metar*>
   {

      private:

         const Dstring
         station_name;

      public:

         Station_Metars (const Metars& metars,
                         const Dstring& station_name);

         Station_Metars (const Metars& metars,
                         const Dstring& station_name,
                         const Dtime& start_time,
                         const Dtime& end_time);

         const Metar&
         get_metar (const Dtime& time) const;

   };

};

#endif /* DENISE_METCODE_H */

