//
// metcode.cc
// 
// Copyright (C) 2012 Simon E. Ching
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

#include "metcode.h"
#include <denise/thermo.h>

Location::Location ()
{
}

Location::Location (const string& str,
                    const string& separator,
                    const Integer id_index,
                    const Integer name_index,
                    const Integer latitude_index,
                    const Integer longitude_index)
{

   const Tokens tokens (str, separator.c_str ());

   id = tokens[id_index];
   name = tokens[name_index];
   lat_long.latitude = atof (tokens[latitude_index].c_str ());
   lat_long.longitude = atof (tokens[longitude_index].c_str ());
   height = GSL_NAN;

   trim (id);
   trim (name);

}

Location::Location (const string& str,
                    const string& separator,
                    const Integer id_index,
                    const Integer name_index,
                    const Integer latitude_index,
                    const Integer longitude_index,
                    const Integer height_index)
{

   const Tokens tokens (str, separator.c_str ());

   id = tokens[id_index];
   name = tokens[name_index];
   lat_long.latitude = atof (tokens[latitude_index].c_str ());
   lat_long.longitude = atof (tokens[longitude_index].c_str ());
   height = atof (tokens[height_index].c_str ());

   trim (id);
   trim (name);

}

Location::Location (const string& id,
                    const string& name,
                    const Lat_Long& lat_long)
   : id (id),
     name (name),
     lat_long (lat_long)
{
}

Location::Location (const string& id,
                    const string& name,
                    const Lat_Long& lat_long,
                    const Real height)
   : id (id),
     name (name),
     lat_long (lat_long),
     height (height)
{
}

bool
Location::operator == (const Location& location) const
{
   return id == location.id &&
          fabs (lat_long.longitude - location.lat_long.longitude < 0.000001) &&
          fabs (lat_long.latitude - location.lat_long.latitude < 0.000001) &&
          fabs (height - location.height < 0.000001);
}

bool
Location::operator < (const Location& location) const
{
   if (id != location.id)
   {
      return id < location.id;
   }
   else
   if (lat_long.longitude != location.lat_long.longitude)
   {
      return lat_long.longitude < location.lat_long.longitude;
   }
   else
   if (lat_long.latitude != location.lat_long.latitude)
   {
      return lat_long.latitude < location.lat_long.latitude;
   }
   else
   if (height != location.height)
   {
      return height < location.height;
   }
}

const Location&
Location_Map::get_nearest (const Lat_Long& lat_long) const
{

   Real min_distance = GSL_POSINF;
   typedef Location_Multimap::const_iterator Iterator;

   Iterator nearest_iterator = end ();

   for (Iterator iterator = begin (); iterator != end (); iterator++)
   {

      const Location& location = iterator->second;
      const Lat_Long& ll = location.lat_long;
      if (ll.is_nall ()) { continue; }

      const Real distance = Geodesy::get_distance (lat_long, ll);

      if (distance < min_distance)
      {
         min_distance = distance;
         nearest_iterator = iterator;
      }

   }

   if (nearest_iterator == end ()) { throw Exception (); }
   return nearest_iterator->second;

}

const Location&
Location_Multimap::get_nearest (const Lat_Long& lat_long) const
{

   Real min_distance = GSL_POSINF;
   typedef Location_Multimap::const_iterator Iterator;

   Iterator nearest_iterator = end ();

   for (Iterator iterator = begin (); iterator != end (); iterator++)
   {

      const Location& location = iterator->second;
      const Lat_Long& ll = location.lat_long;
      if (ll.is_nall ()) { continue; }

      const Real distance = Geodesy::get_distance (lat_long, ll);

      if (distance < min_distance)
      {
         min_distance = distance;
         nearest_iterator = iterator;
      }

   }

   if (nearest_iterator == end ()) { throw Exception (); }
   return nearest_iterator->second;

}


Metar_Wind::Metar_Wind ()
   : Wind ()
{
}

Metar_Wind::Metar_Wind (const string& metar_wind_string)
{

   if (metar_wind_string.substr (0, 3) == "///")
   {
      this->u = GSL_NAN;
      this->v = GSL_NAN;
      this->gust = GSL_NAN;
      this->variable_direction = false;
      return;
   }

   const string& str = metar_wind_string;
   const string& dir_str = str.substr (0, 3);
   const bool no_gust = (metar_wind_string.size () < 10);

   const Integer unit_pos = (no_gust ? 5 : 8);
   const string unit_str = str.substr (unit_pos, 2);

   Real direction = atof (dir_str.c_str ());
   Real speed = atof (str.substr (3, 2).c_str ());
   Real gust = (no_gust ? GSL_NAN : atof (str.substr (6, 2).c_str ()));

   if (unit_str == "KT")
   {
      const Real kt_to_ms = 1.852 / 3.6;
      speed *= kt_to_ms;
      gust *= kt_to_ms;
   }

   set_from_direction_speed (direction, speed);
   this->gust = gust;
   this->variable_direction = (str.substr (0, 3) == "000");

}

Real
Metar_Wind::get_gust () const
{
   return gust;
}

using namespace std;
using namespace denise;

Metar::Key::Key ()
{
}

Metar::Key::Key (const Dtime& time,
                 const string& station_name)
   : time (time),
     station_name (station_name)
{
}

bool
Metar::Key::operator == (const Key& key) const
{
   return (fabs (time.t - key.time.t) <= METAR_TIME_TOLERANCE) &&
      (station_name == key.station_name);
}

bool
Metar::Key::operator > (const Key& key) const
{
   const bool bt = (fabs (time.t - key.time.t) > METAR_TIME_TOLERANCE);
   if (bt) { return time > key.time; }
   else { return station_name > key.station_name; }
}

bool
Metar::Key::operator < (const Key& key) const
{
   const bool bt = (fabs (time.t - key.time.t) > METAR_TIME_TOLERANCE);
   if (bt) { return time.t < key.time.t; }
   else { return station_name < key.station_name; }
}

void
Metar::read_t_td (const string& token)
{

   if (token == "/////////")
   {
      this->temperature = GSL_NAN;
      this->dew_point = GSL_NAN;
      return;
   }

   string::size_type slash = token.find_first_of ('/');
   const string& t_str = token.substr (0, slash);
   const string& td_str = token.substr (slash + 1);

   if (t_str.substr (0, 2) == "//")
   {
      this->temperature = GSL_NAN;
   }
   else
   if (t_str.substr (0, 2) == "MS")
   {
      this->temperature = -1 * atof (t_str.substr (2).c_str ());
   }
   else
   {
      this->temperature = atof (t_str.c_str ());
   }

   if (this->temperature < -40 || this->temperature > 60)
   {
      this->temperature = GSL_NAN;
   }

   if (td_str.substr (0, 2) == "//")
   {
      this->dew_point = GSL_NAN;
   }
   else
   if (td_str.substr (0, 2) == "MS")
   {
      this->dew_point = -1 * atof (td_str.substr (2).c_str ());
   }
   else
   {
      this->dew_point = atof (td_str.c_str ());
   }

   if (this->dew_point < -40 || this->dew_point > 40)
   {
      this->dew_point = GSL_NAN;
   }

}

void
Metar::read_qnh (const string& token)
{

   if (token.substr (0, 4) == "////")
   {
      this->qnh = GSL_NAN;
      return;
   }

   this->qnh = atof (token.c_str ());

   if (this->qnh < 880 || this->qnh > 1060)
   {
      this->qnh = GSL_NAN;
   }

}

void
Metar::read_rainfall (const string& token)
{

   if (token.size () == 12)
   {
      this->rain_10_minute = atof (token.substr (2, 4).c_str ());
      this->rain_hour = GSL_NAN;
      this->rain_since_9am = atof (token.substr (7, 5).c_str ());
   }
   else
   if (token.size () == 18)
   {
      this->rain_10_minute = atof (token.substr (2, 4).c_str ());
      this->rain_hour = atof (token.substr (7, 5).c_str ());
      this->rain_since_9am = atof (token.substr (13, 5).c_str ());
   }
   else
   {
      this->rain_10_minute = GSL_NAN;
      this->rain_hour = GSL_NAN;
      this->rain_since_9am = GSL_NAN;
   }

}

Metar::Metar (const string& metar_string)
   : visibility (GSL_NAN),
     visibility_auto (GSL_NAN)
{

   const string wx_str_a ("MI|DR|BL|SH|TS|FZ|DZ|RA|SN|SG|IC|PL|GR");
   const string wx_str_b ("|GS|BR|FG|FU|VA|DI|SA|HZ|PO|SQ|FC|SS|DS");
   const string wx_str = wx_str_a + wx_str_b;

   const Reg_Exp wx_regexp (wx_str);
   const Reg_Exp colon_regexp (":");
   const Reg_Exp cloud_regexp ("[0-9][A-Z].[0-9]..");

   vector<string> tokens = tokenize (metar_string);
   const Integer n = tokens.size ();
   const string date_format ("%Y%m%d%H%M");

   this->speci = (tokens[0] == "SPECIAWS");
   this->key = Metar::Key (Dtime (tokens[2], date_format), tokens[1]);

   if (n < 4) { return; }
   this->metar_wind = Metar_Wind (tokens[3]);

   // CAVOK or NOVIS
   if (n < 5) { return; }
   bool cavok = (tokens[4] == "CAVOK");
   bool no_vis = (tokens[4] == "////");

   if (!cavok && !no_vis)
   {
      const string& vis_str = tokens[4];
      visibility = atof (vis_str.c_str ());
   }

   if (n < 6) { return; }
   Integer i = 5;

   if (!cavok && !no_vis)
   {

      const string& vis_str = tokens[i];
      i++;

      while (wx_regexp.match (tokens[i]))
      {
         weather += tokens[i] + " ";
         i++;
      }

      while (cloud_regexp.match (tokens[i]))
      {
         cloud += tokens[i] + " ";
         i++;
      }

   }

   // Temperature and Dew Point
   read_t_td (tokens[i]);
   i++;
   if (i >= n) { return; }

   // QNH
   read_qnh (tokens[i]);

   // RMK
   while (i < n)
   {

      i++;
      if (i >= n) { return; }

      // Rainfall
      if (tokens[i].substr (0, 2) == "RF")
      {
         read_rainfall (tokens[i]);
         continue;
      }

/*
      // cloud_auto group
      bool first_cloud_auto = (tokens[i].substr (0, 3) == "CLD");
      if (first_cloud_auto)
      {
         while ((i < n) && (first_cloud_auto || !colon_regexp.match (tokens[i])))
         {
            cloud_auto += tokens[i] + " ";
            i++;
            first_cloud_auto = false;
         }
         continue;
      }

      // visibility_auto group
      const string& vis_auto_str = tokens[i];
      if (vis_auto_str.substr (0, 3) == "VIS")
      {
         visibility_auto = atof (vis_auto_str.substr (4, 4).c_str ());
         continue;
      }
*/

   }

}

string
Metar::get_string () const
{

   const string& station_name = get_station_name ();
   const Dtime& dtime = get_time ();
   const Real& t = get_temperature ();
   const Real& td = get_dew_point ();
   const Real& qnh = get_qnh ();
   const Real& rain_10_minutes = get_rain_10_minute ();
   const Real& rain_hour = get_rain_hour ();
   const Real& rain_since_9am = get_rain_since_9am ();
   const Real vis = get_visibility ();

   const Metar_Wind& metar_wind = get_metar_wind ();
   const Real wind_dir = metar_wind.get_direction ();
   const Real wind_speed = metar_wind.get_speed ();
   const Real wind_gust = metar_wind.get_gust ();

   const string& time_str = dtime.get_string ("%d%H%M");
   const string& ttd_str = string_render ("%04.1f/%04.1f", t, td);
   const string& qnh_str = string_render ("Q%06.1f", qnh);
   const string& vis_str = (gsl_isnan (vis) ? "////" : string_render ("%04.0f", vis));
   const string& wind_str = string_render ("%03.0f%02.0f/%02.0fKT",
      wind_dir, wind_speed / 0.514444, wind_gust / 0.514444);
   const string& rf9am_str = string_render ("RF%.1f", rain_since_9am);

   return (station_name + " " + time_str + " " + wind_str + " "
      + vis_str + " " + ttd_str + " " + qnh_str + " " + rf9am_str);

}

const Dtime&
Metar::get_time () const
{
   return key.time;
}

const string&
Metar::get_station_name () const
{
   return key.station_name;
}

const Real&
Metar::get_temperature () const
{
   return temperature;
}

const Real&
Metar::get_dew_point () const
{
   return dew_point;
}

const Real
Metar::get_rh () const
{
   return Moisture::get_rh (temperature, dew_point);
}

const Metar_Wind&
Metar::get_metar_wind () const
{
   return metar_wind;
}

const Real&
Metar::get_qnh () const
{
   return qnh;
}

const Real&
Metar::get_rain_10_minute () const
{
   return rain_10_minute;
}

const Real&
Metar::get_rain_hour () const
{
   return rain_hour;
}

const Real&
Metar::get_rain_since_9am () const
{
   return rain_since_9am;
}

const Real
Metar::get_visibility () const
{
   return gsl_isnan (visibility) ? visibility_auto : visibility;
}

const Real
Metar::get_gfdi () const
{
   const Real rh = get_rh () * 100;
   const Real kph = metar_wind.get_speed () * 3.6;
   return Fire::get_gfdi (temperature, rh, kph);
}

const Real
Metar::get_ffdi () const
{
   const Real rh = get_rh () * 100;
   const Real kph = metar_wind.get_speed () * 3.6;
   return Fire::get_ffdi (temperature, rh, kph);
}

bool
Metar::operator == (const Metar& metar) const
{
   return (key == metar.key);
}

bool
Metar::operator > (const Metar& metar) const
{
   return (key > metar.key);
}

bool
Metar::operator < (const Metar& metar) const
{
   return (key < metar.key);
}

void
Metars::clear_station_metars_ptr_map ()
{

   typedef map<string, Station_Metars*>::const_iterator Iterator;
   station_metars_ptr_map.clear ();

   for (Iterator iterator = station_metars_ptr_map.begin ();
        iterator != station_metars_ptr_map.end (); iterator++)
   {
      Station_Metars* station_metars_ptr = iterator->second;
      delete station_metars_ptr;
   }

}

Metars::Metars ()
{
}

Metars::Metars (const string& dir_path,
                const string& file_format)
   : dir_path (dir_path),
     file_format (file_format)
{
   try
   {
      read (dir_path, file_format);
   }
   catch (const IO_Exception& ioe)
   {
      cerr << ioe << endl;
   }
}

Metars::Metars (const string& dir_path,
                const string& file_format,
                const Location_Map& location_map)
   : dir_path (dir_path),
     file_format (file_format)
{
   try
   {
      read (dir_path, file_format, location_map);
   }
   catch (const IO_Exception& ioe)
   {
      cerr << ioe << endl;
   }

}

Metars::Metars (const string& dir_path,
                const string& file_format,
                const Location_Multimap& location_multimap)
   : dir_path (dir_path),
     file_format (file_format)
{
   try
   {
      read (dir_path, file_format, location_multimap);
   }
   catch (const IO_Exception& ioe)
   {
      cerr << ioe << endl;
   }

}

Metars::~Metars ()
{
   station_metars_ptr_map.clear ();
}

void
Metars::setup (const string& dir_path,
               const string& file_format)
{
   this->dir_path = dir_path;
   this->file_format = file_format;
}

string
Metars::get_status () const
{
   if (time_set.size () == 0) { return "Unloaded"; }
   set<Dtime>::const_reverse_iterator i = time_set.rbegin ();
   return i->get_string ("Latest %Y.%m.%d %H:%M");
}

void
Metars::reload ()
{
   if ((dir_path == "") || file_format == "") { return; }
   read (dir_path, file_format);
   construct_station_metars_ptr_map ();
}

void
Metars::reload (const Location_Map& location_map)
{
   if ((dir_path == "") || file_format == "") { return; }
   read (dir_path, file_format, location_map);
   construct_station_metars_ptr_map ();
}

void
Metars::reload (const Location_Multimap& location_multimap)
{
   if ((dir_path == "") || file_format == "") { return; }
   read (dir_path, file_format, location_multimap);
   construct_station_metars_ptr_map ();
}

void
Metars::read (const string& dir_path,
              const string& file_format)
{

   typedef vector<string> Tokens;
   const Tokens& dir_listing = get_dir_listing (dir_path, file_format);

//   for (Tokens::const_iterator iterator = dir_listing.begin ();
//        iterator != dir_listing.end (); iterator++)
//   {
//
//      const string& file_name = *(iterator);

   for (Integer i = 0; i < dir_listing.size (); i++)
   {
      const string& file_name = dir_listing[i];
      const string& file_path = dir_path + "/" + file_name;

cout << "metar reading " << file_path << endl;
      read (file_path);

   }

}

void
Metars::read (const string& dir_path,
              const string& file_format,
              const Location_Map& location_map)
{

   typedef vector<string> Tokens;
   const Tokens& dir_listing = get_dir_listing (dir_path, file_format);

   for (Integer i = 0; i < dir_listing.size (); i++)
   {
      const string& file_name = dir_listing[i];
      const string& file_path = dir_path + "/" + file_name;

cout << "metar reading " << file_path << endl;
      read (file_path, location_map);

   }

}

void
Metars::read (const string& dir_path,
              const string& file_format,
              const Location_Multimap& location_multimap)
{

   typedef vector<string> Tokens;
   const Tokens& dir_listing = get_dir_listing (dir_path, file_format);

   for (Integer i = 0; i < dir_listing.size (); i++)
   {
      const string& file_name = dir_listing[i];
      const string& file_path = dir_path + "/" + file_name;

cout << "metar reading " << file_path << endl;
      read (file_path, location_multimap);

   }

}

void
Metars::read (const string& file_path)
{

   char input_line_str[1024];
   gzFile file = get_gzfile (file_path);

   while (gz_readline (input_line_str, 1024, file) != NULL)
   {

      const string metar_string (input_line_str);
      const Metar metar (metar_string);
      const string& station_name = metar.get_station_name ();
      const Dtime& dtime = metar.get_time ();

      insert (metar);

      const Lat_Long lat_long (GSL_NAN, GSL_NAN);
      station_name_lat_long_map.insert (make_pair (station_name, lat_long));
      time_set.insert (dtime);

   }

   gzclose (file);

}

void
Metars::read (const string& file_path,
              const Location_Map& location_map)
{

   char input_line_str[1024];
   gzFile file = get_gzfile (file_path);

   while (gz_readline (input_line_str, 1024, file) != NULL)
   {

      const string metar_string (input_line_str);
      const Metar metar (metar_string);
      const string& station_name = metar.get_station_name ();
      const Dtime& dtime = metar.get_time ();

      insert (metar);

      typedef Location_Map::const_iterator Iterator;
      Iterator i = location_map.find (station_name);
      const Lat_Long lat_long = ((i == location_map.end ()) ?
         Lat_Long (GSL_NAN, GSL_NAN) : (i->second).lat_long);
      station_name_lat_long_map.insert (make_pair (station_name, lat_long));
      time_set.insert (dtime);

   }

   gzclose (file);

}

void
Metars::read (const string& file_path,
              const Location_Multimap& location_multimap)
{

   char input_line_str[1024];
   gzFile file = get_gzfile (file_path);

   while (gz_readline (input_line_str, 1024, file) != NULL)
   {

      const string metar_string (input_line_str);
      const Metar metar (metar_string);
      const string& station_name = metar.get_station_name ();
      const Dtime& dtime = metar.get_time ();

      insert (metar);

      typedef Location_Multimap::const_iterator Iterator;
      Iterator i = location_multimap.find (station_name);
      const Lat_Long lat_long = ((i == location_multimap.end ()) ?
         Lat_Long (GSL_NAN, GSL_NAN) : (i->second).lat_long);
      station_name_lat_long_map.insert (make_pair (station_name, lat_long));
      time_set.insert (dtime);

   }

   gzclose (file);

}

void
Metars::construct_station_metars_ptr_map ()
{

   station_metars_ptr_map.clear ();
   typedef map<string, Lat_Long>::const_iterator Iterator;

   for (Iterator iterator = station_name_lat_long_map.begin ();
        iterator != station_name_lat_long_map.end (); iterator++)
   {
      const string& sn = (iterator->first);
      Station_Metars* station_metars_ptr = new Station_Metars (*this, sn);
      station_metars_ptr_map.insert (make_pair (sn, station_metars_ptr));
   }

}

const map<string, Lat_Long>&
Metars::get_station_name_lat_long_map () const
{
   return station_name_lat_long_map;
}

void
Metars::attract (Real& latitude,
                 Real& longitude) const
{

   map<string, Lat_Long>::const_iterator iterator =
      get_nearest_station_name_lat_long (Lat_Long (latitude, longitude));

   const Lat_Long& ll = iterator->second; 
   latitude = ll.latitude;
   longitude = ll.longitude;

}

map<string, Lat_Long>::const_iterator
Metars::get_nearest_station_name_lat_long (const Lat_Long& lat_long) const
{

   Real min_distance = GSL_POSINF;
   typedef map<string, Lat_Long>::const_iterator Iterator;

   Iterator nearest_iterator = station_name_lat_long_map.end ();

   for (Iterator i = station_name_lat_long_map.begin ();
        i != station_name_lat_long_map.end (); i++)
   {

      const Lat_Long& ll = (i->second);
      if (ll.is_nall ()) { continue; }
      const Real distance = Geodesy::get_distance (lat_long, ll);

      if (distance < min_distance)
      {
         min_distance = distance;
         nearest_iterator = i;
      }

   }

   if (nearest_iterator == station_name_lat_long_map.end ())
   {
      throw Exception ();
   }

   return nearest_iterator;

}

const Metar&
Metars::get_nearest_metar (const Dtime& time,
                           const Lat_Long& lat_long) const
{

   Real min_distance = GSL_POSINF;
   typedef map<string, Lat_Long>::const_iterator Iterator;

   Iterator nearest_iterator = station_name_lat_long_map.end ();

   for (Iterator i = station_name_lat_long_map.begin ();
        i != station_name_lat_long_map.end (); i++)
   {

      const Lat_Long& ll = (i->second);
      if (ll.is_nall ()) { continue; }
      const Real distance = Geodesy::get_distance (lat_long, ll);

      if (distance < min_distance)
      {
         min_distance = distance;
         nearest_iterator = i;
      }

   }

   if (nearest_iterator == station_name_lat_long_map.end ())
   {
      throw Exception ();
   }

   const string& station_name = nearest_iterator->first;
   return get_metar (time, station_name);

}

const Lat_Long&
Metars::get_lat_long (const string& station_name) const
{
   typedef map<string, Lat_Long>::const_iterator Iterator;
   Iterator i = station_name_lat_long_map.find (station_name);
   if (i == station_name_lat_long_map.end ()) { throw Exception (); }
   return (i->second);
}

const set<Dtime>&
Metars::get_time_set () const
{
   return time_set;
}

const Station_Metars&
Metars::get_station_metars (const string& station_name) const
{
   typedef map<string, Station_Metars*> Station_Metars_Ptr_Map;
   typedef Station_Metars_Ptr_Map::const_iterator Iterator;
   const Station_Metars_Ptr_Map& smpm = station_metars_ptr_map;
   Iterator iterator = smpm.find (station_name);
   if (iterator == smpm.end ()) { throw Exception ("No Such Station"); }
   return *(iterator->second);
}

Station_Metars&
Metars::get_station_metars (const string& station_name)
{
   typedef map<string, Station_Metars*> Station_Metars_Ptr_Map;
   typedef Station_Metars_Ptr_Map::iterator Iterator;
   Station_Metars_Ptr_Map& smpm = station_metars_ptr_map;
   Iterator iterator = smpm.find (station_name);
   if (iterator == smpm.end ()) { throw Exception ("No Such Station"); }
   return *(iterator->second);
}

const Station_Metars&
Metars::get_station_metars (const Lat_Long& lat_long) const
{
   typedef map<string, Lat_Long>::const_iterator Iterator;
   Iterator i = get_nearest_station_name_lat_long (lat_long);
   const string& station_name = i->first;
   return get_station_metars (station_name);
}

Station_Metars&
Metars::get_station_metars (const Lat_Long& lat_long)
{
   typedef map<string, Lat_Long>::const_iterator Iterator;
   Iterator i = get_nearest_station_name_lat_long (lat_long);
   const string& station_name = i->first;
   return get_station_metars (station_name);
}

const Metar&
Metars::get_metar (const Dtime& time,
                   const string& station_name) const
{
   const Station_Metars& sm = get_station_metars (station_name);
   return sm.get_metar (time);
}

vector<const Metar*>
Metars::get_metar_ptr_vector (const Dtime& time) const
{

   vector<const Metar*> metar_ptr_vector;
   typedef map<string, Lat_Long>::const_iterator Iterator;

   for (Iterator iterator = station_name_lat_long_map.begin ();
        iterator != station_name_lat_long_map.end (); iterator++)
   {
      const string& station_name = (iterator->first);
      try
      {
         const Metar& metar = get_metar (time, station_name);
         metar_ptr_vector.push_back (&metar);
      }
      catch (const Exception& e)
      {
      }
   }

   return metar_ptr_vector;

}

pair<string, Lat_Long>
Metars::nearest (const Lat_Long& lat_long) const
{
   try
   {
      return *(get_nearest_station_name_lat_long (lat_long));
   }
   catch (const Exception& e)
   {
      return make_pair (string (""), lat_long);
   }
}

Station_Metars::Station_Metars (const Metars& metars,
                                const string& station_name)
   : station_name (station_name)
{

   for (Metars::const_iterator iterator = metars.begin ();
        iterator != metars.end (); iterator++)
   {

      const Metar& metar = *(iterator);
      const string& sn = metar.get_station_name ();

      if (sn == station_name)
      {
         const Dtime& time = metar.get_time ();
         insert (make_pair (time, &metar));
      }
   }

}

const Metar&
Station_Metars::get_metar (const Dtime& time) const
{

   Station_Metars::const_iterator lb = lower_bound (time);
   Station_Metars::const_iterator ub = upper_bound (time);

   if (lb == ub) // Not at a node
   {

      if (ub == lower_bound (Dtime (GSL_NEGINF))) // Before first node
      {
         throw Exception ("Out of bounds:");
      }
      else
      {
         lb--;
         return *(lb->second);
      }
   }
   else // At a node
   {
      return *(lb->second);

   }

}


