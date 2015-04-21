//
// gts.cc
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
// but WITHOUT ANY WARRANTY), without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with libdenise.  If not, see <http://www.gnu.org/licenses/>.

#include <fstream>
#include "gts.h"

using namespace std;
using namespace denise;

Gts::Gts (const string& message,
          const Dtime& time_hint)
   : Tokens (message),
     time_hint (time_hint)
{
   if (size () < 3) { throw Exception ("Empty GTS"); }
//   if (((*this)[0]).substr (0, 3) == "NIL") { throw Exception ("Empty GTS"); }
   if (((*this)[1]).substr (0, 3) == "NIL") { throw Exception ("Empty GTS"); }
   if (((*this)[2]).substr (0, 3) == "NIL") { throw Exception ("Empty GTS"); }
}

string
Gts::get_yygg () const
{
   const string& yyggi = (*this)[1];
   const Integer yy = atoi (yyggi.substr (0, 2).c_str ());
   const Integer gg = atoi (yyggi.substr (2, 2).c_str ());
   return string_render ("%02d%02d", yy - (yy >= 50 ? 50 : 0), gg);
}

string
Gts::get_time_string () const
{

   const string& yygg = get_yygg ();
   const string yyyymmdd = time_hint.get_string ("%Y%m%d");
   const bool same_day = (yygg.substr (0, 2) == yyyymmdd.substr (6, 2));

   const Dtime yyyymm_time (time_hint.t - (same_day ? 0 : 24));
   const string& yyyymm = yyyymm_time.get_string ("%Y%m");

   return yyyymm + yygg;

}

const string&
Gts::get_wmo_id () const
{
//   return ((*this)[2]).substr (0, 5);
   return ((*this)[2]);
}

string
Gts::get_key () const
{
   return get_wmo_id () + ":" + get_time_string ();
}

Tokens*
Gts::parse_file (const string& file_path)
{

   bool first_sticker;
   string message, input_string, sticker;

   Reg_Exp block_delimiter ("[\001\003]");
   Reg_Exp block_header_1 ("^[0-9]..");
   Reg_Exp block_header_2 ("^[A-Z0-9]..... [A-Z0-9]... [A-Z0-9].....");
   Reg_Exp empty_line ("^[ \015]*$");
   Reg_Exp white_head ("^[^A-Z0-9/=]+", true);
   Reg_Exp white_tail ("[^A-Z0-9/=]+$", true);
   Reg_Exp end_of_message ("=", true);
   Reg_Exp control_char ("[\001-\037]+", true);
   Reg_Exp aaxx ("^AAXX *[0-9]....[^A-Z0-9/=]*$");
   Reg_Exp bbxx ("^BBXX[^A-Z0-9/=]*$");

   ifstream file (file_path.c_str ());

   Tokens* message_tokens_ptr = new Tokens ();

   int i;
   for (i = 0; std::getline (file, input_string); i++)
   {

      if (block_delimiter.match (input_string))
      {
         i = 0;
         sticker = "";
         first_sticker = true;
      }

      control_char.replace (input_string, "");
      white_head.replace (input_string, "");
      white_tail.replace (input_string, "");

      if (i == 1)
      {
         if (!block_header_1.match (input_string))
         {
            i = -1;
         }
         continue;
      }

      if (i == 2)
      {
         if (!block_header_2.match (input_string))
         {
            i = -1;
         }
         continue;
      }
      if (i == 3)
      {
         if (aaxx.match (input_string) || bbxx.match (input_string))
         {
            sticker = input_string;
         }
      }

      if (empty_line.match (input_string))
      {
         if (!empty_line.match (message))
         {
            if (sticker != "")
            {
               if (first_sticker) { first_sticker = false; }
               else { message = sticker + " " + message; }
            }
            end_of_message.replace (message, "");
            message_tokens_ptr->push_back (message);
         }
         message.clear ();
      }
      else
      if (end_of_message.match (input_string))
      {
         message += input_string;
         if (sticker != "")
         {
            if (first_sticker) { first_sticker = false; }
            else { message = sticker + " " + message; }
         }
         end_of_message.replace (message, "");
         message_tokens_ptr->push_back (message);
         message.clear ();
      }
      else
      {
         message += input_string + " ";
      }

   }

   file.close ();
   return message_tokens_ptr;

}

Gts_Ptr_Store::~Gts_Ptr_Store ()
{
   for (Gts_Ptr_Store::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Gts* gts_ptr = iterator->second;
//      delete gts_ptr;
   }
}

void
Gts_Ptr_Store::ingest (const string& message,
                       const Dtime& time_hint)
{

   if (Reg_Exp ("NIL").match (message)) { return; }

   Gts* gts_ptr = NULL;
   const string& header = message.substr (0, 5);

   try
   {
           if (header == "TTAA ") { gts_ptr = (new Ttaa (message, time_hint)); }
      else if (header == "TTBB ") { gts_ptr = (new Ttbb (message, time_hint)); }
      else if (header == "TTCC ") { gts_ptr = (new Ttcc (message, time_hint)); }
      else if (header == "TTDD ") { gts_ptr = (new Ttdd (message, time_hint)); }
   }
   catch (const Exception& e) { }

   if (gts_ptr == NULL) { return; }

   const string gts_key = gts_ptr->get_key ();
   insert (make_pair (gts_key, gts_ptr));

}

const Gts&
Gts_Ptr_Store::get_gts (const string& gts_key) const
{
   Gts_Ptr_Store::const_iterator i = find (gts_key);
   if (i == end ()) { throw No_Match_Exception ("GTS with this key N/A"); }
   return *(i->second);
}

void
Gts_Ptr_Grandstore::ingest (const string& genre,
                            const string& message,
                            const Dtime& time_hint)
{

   Gts_Ptr_Grandstore::iterator i = find (genre);
   if (i == end ()) { i = insert (i, make_pair (genre, Gts_Ptr_Store ())); }

   Gts_Ptr_Store& gps = (i->second);
   gps.ingest (message, time_hint);

}

void
Gts_Ptr_Grandstore::ingest (const string& message,
                            const Dtime& time_hint)
{

   const string& header = message.substr (0, 5); 

        if (header == "TTAA ") { ingest ("TTAA", message, time_hint); }
   else if (header == "TTBB ") { ingest ("TTBB", message, time_hint); }
   else if (header == "TTCC ") { ingest ("TTCC", message, time_hint); }
   else if (header == "TTDD ") { ingest ("TTDD", message, time_hint); }

}

void
Gts_Ptr_Grandstore::ingest_file (const string& file_path)
{

   Reg_Exp yymmddhh_re ("[0-9].......");
   Tokens* message_tokens_ptr = Gts::parse_file (file_path);
   const string file_name = Tokens (file_path, "/").back ();

   const bool match = (yymmddhh_re.match (file_name));

   Dtime time_hint = Dtime::now ();
   if (yymmddhh_re.match (file_name))
   {
      const string& yymmddhh = file_name.substr (0, 8);
      time_hint = Dtime (yymmddhh, string ("%y%m%d%H"));
   }

   for (Tokens::const_iterator iterator = message_tokens_ptr->begin ();
        iterator != message_tokens_ptr->end (); iterator++)
   {
      const string& message = *(iterator);
      ingest (message, time_hint);
   }

   delete message_tokens_ptr;

}

void
Gts_Ptr_Grandstore::ingest_dir (const string& dir_path,
                                const Reg_Exp& file_reg_exp)
{

   const Tokens& dir_listing = get_dir_listing (dir_path, file_reg_exp);

   for (Tokens::const_iterator iterator = dir_listing.begin ();
        iterator != dir_listing.end (); iterator++)
   {
      const string& file_name = *(iterator);
      const string file_path = dir_path + "/" + file_name;
      ingest_file (file_path);
   }

}

const Gts_Ptr_Store&
Gts_Ptr_Grandstore::get_gts_ptr_store (const string& genre) const
{
   Gts_Ptr_Grandstore::const_iterator i = find (genre);
   if (i == end ()) { throw Exception ("Invalid GTS genre"); }
   return i->second;
}

const Gts&
Gts_Ptr_Grandstore::get_gts (const string& genre,
                             const string& gts_key) const
{
   const Gts_Ptr_Store& gts_ptr_store = get_gts_ptr_store (genre);
   return gts_ptr_store.get_gts (gts_key);
}

pair<Real, Real>
Temp::parse_tttdd (const string& tttdd) const
{

   Real temperature = GSL_NAN;
   Real dew_point = GSL_NAN;

   const string& ttt = tttdd.substr (0, 3);
   if (ttt != "///")
   {
      const bool minus = (atoi (tttdd.substr (2, 1).c_str ()) % 2 == 1);
      const Real sign = (minus ? -1 : 1);
      temperature = sign * atof (ttt.c_str ()) * 0.1;
   }

   const string& dd = tttdd.substr (3, 2);
   if (dd  != "//")
   {
      const Real dep = atof (dd.c_str ());
      const bool dry = dep > 55;
      const Real depression = (dry ? dep - 50 : dep * 0.1);
      dew_point = temperature - depression;
   }

   return make_pair (temperature, dew_point);

}

Wind
Temp::parse_ddfff (const string& ddfff) const
{

   Real direction = GSL_NAN;
   Real speed = GSL_NAN;

   if (ddfff != "/////")
   {

      const Integer centre_digit = atoi (ddfff.substr (2, 1).c_str ());
      const bool five = (centre_digit >= 5);

      direction = atof (ddfff.substr (0, 2).c_str ()) * 10;
      speed = atof (ddfff.substr (2, 3).c_str ());

      if (five) { direction += 5; speed -= 500; }
      if (use_knots ()) { speed *= 1.852 / 3.6; }

   }

   return Wind::direction_speed (direction, speed);

}

bool
Temp::use_knots () const
{
   const string& yyggi = (*this)[1];
   const Integer yy = atoi (yyggi.substr (0, 2).c_str ());
   return (yy >= 50);
}

Temp::Temp (const string& message,
            const Dtime& time_hint)
   : Gts (message, time_hint)
{
   for (Integer i = 1; i < size (); i++)
   {
      const string& group = at (i);
      const Integer n = group.size ();
      if (n != 5) { throw Exception ("Invalid TEMP message"); }
   }
}

const char
Ttac::get_highest_wind_code () const
{
   const string& yyggi = (*this)[1];
   return yyggi[4];
}

Integer
Ttac::parse_special_level (Sounding& sounding,
                           Integer& index) const
{

   if (index >= size ()) { return index; }

   const string& nnppp = (*this)[index++];
   const Integer nn = atoi (nnppp.substr (0, 2).c_str ());


   switch (nn)
   {

      case 88:
      case 99:
      {

         const string& ppp = nnppp.substr (2, 3);
         if (ppp == "999") { return index; }

         const Real p = atof (ppp.c_str ()) * 1e2;
         const Real pressure = (p < 500e2 ? p + 1000e2 : p);
         if (!gsl_finite (pressure)) { return index; }

         const string& tttdd = (*this)[index++];
         const string& ddfff = (*this)[index++];

         const bool station_level = (nn == 99);

         if (station_level)
         {
            Height_Profile& height_profile = sounding.get_height_profile ();
            height_profile.add (pressure, 0);
         }

         const pair<Real, Real>& tttdd_pair = parse_tttdd (tttdd);
         const Real temperature = tttdd_pair.first;
         const Real dew_point = tttdd_pair.second;

         if (gsl_finite (temperature))
         {
            Thermo_Line& t_line = sounding.get_t_line ();
            t_line.add (pressure, temperature);
         }

         if (gsl_finite (dew_point))
         {
            Thermo_Line& t_d_line = sounding.get_t_d_line ();
            t_d_line.add (pressure, dew_point);
         }

         const Wind& wind = parse_ddfff (ddfff);
         if (!wind.is_naw ())
         {
            Wind_Profile& wind_profile = sounding.get_wind_profile ();
            wind_profile.add (pressure, wind);
         }

         break;
      }

   }

   return index;

}

Integer
Ttac::parse_standard_levels (Sounding& sounding,
                             Integer& index) const
{

   //const Integer n = get_number_of_standard_levels ();
   const char highest_wind_code = get_highest_wind_code ();

   while (index < size ())
   {

      const string& pphhh = (*this)[index++];
      const char& p = pphhh[0];
      if (p != '0' && (p == pphhh[1])) { return --index; }

      const string& tttdd = (*this)[index++];

      // only 2 groups
      const bool with_wind_group = ((highest_wind_code != '/')  &&
         ((p >= highest_wind_code) || (p == '0')));


      const pair<Real, Real>& pphhh_pair = parse_pphhh (pphhh);
      const Real pressure = pphhh_pair.first;
      const Real geopotential_height = pphhh_pair.second;

      if (!gsl_finite (pressure)) { continue; }
      Height_Profile& height_profile = sounding.get_height_profile ();
      height_profile.add (pressure, geopotential_height);

      const pair<Real, Real>& tttdd_pair = parse_tttdd (tttdd);
      const Real temperature = tttdd_pair.first;
      const Real dew_point = tttdd_pair.second;

      if (gsl_finite (temperature))
      {
         Thermo_Line& t_line = sounding.get_t_line ();
         t_line.add (pressure, temperature);
      }

      if (gsl_finite (dew_point))
      {
         Thermo_Line& t_d_line = sounding.get_t_d_line ();
         t_d_line.add (pressure, dew_point);
      }

      if (with_wind_group)
      {

         const string& ddfff = (*this)[index++];
         const Wind& wind = parse_ddfff (ddfff);

         if (!wind.is_naw ())
         {
            Wind_Profile& wind_profile = sounding.get_wind_profile ();
            wind_profile.add (pressure, wind);
         }

      }

   }

   return index;

}

pair<Real, Real>
Ttac::parse_pphhh (const string& pphhh) const
{

   Real pressure = GSL_NAN;
   Real geopotential_height = GSL_NAN;

   const Integer pp = atoi (pphhh.substr (0, 2).c_str ());
   const string& hhh = pphhh.substr (2, 3);
   if (hhh != "///") { geopotential_height = atof (hhh.c_str ()); }

   interpret_p_z (pp, pressure, geopotential_height);
   return make_pair (pressure, geopotential_height);

}

Ttac::Ttac (const string& message,
            const Dtime& time_hint)
   : Temp (message, time_hint)
{
}

bool
Ttbd::parse_significant_temperature_levels (Sounding& sounding,
                                            Integer& index) const
{

   while (index < size ())
   {

      const string& nnppp = (*this)[index++];
      // reached end of message for temperatures
      if (nnppp[0] != nnppp[1]) { return false; }
      if (nnppp == "21212") { return false; }

      const string& tttdd = (*this)[index++];

      const Real pressure = parse_nnppp (nnppp);
      if (!gsl_finite (pressure)) { continue; }

      const pair<Real, Real>& tttdd_pair = parse_tttdd (tttdd);
      const Real temperature = tttdd_pair.first;
      const Real dew_point = tttdd_pair.second;

      if (gsl_finite (temperature))
      {
         Thermo_Line& t_line = sounding.get_t_line ();
         t_line.add (pressure, temperature);
      }

      if (gsl_finite (dew_point))
      {
         Thermo_Line& t_d_line = sounding.get_t_d_line ();
         t_d_line.add (pressure, dew_point);
      }

   }

   return true;

}

bool
Ttbd::parse_significant_wind_levels (Sounding& sounding,
                                     Integer& index) const
{

   while (index < size ())
   {

      const string& nnppp = (*this)[index++];
      if (nnppp[0] != nnppp[1]) { return false; }
      if (nnppp == "31313") { return false; }

      const string& ddfff = (*this)[index++];

      // reached end of message for temperatures

      const Real pressure = parse_nnppp (nnppp);
      if (!gsl_finite (pressure)) { continue; }

      const Wind& wind = parse_ddfff (ddfff);
      if (!wind.is_naw ())
      {
         Wind_Profile& wind_profile = sounding.get_wind_profile ();
         wind_profile.add (pressure, wind);
      }

   }

   return true;

}

Ttbd::Ttbd (const string& message,
            const Dtime& time_hint)
   : Temp (message, time_hint)
{
}

void
Ttaa::interpret_p_z (const Integer pp,
                     Real& pressure,
                     Real& geopotential_height) const
{

   Real& p = pressure;
   Real& z = geopotential_height;

   switch (pp)
   {
      case 0:  p = 1000e2; if (z > 300) { z -= 500; } break;
      case 92: p = 925e2; break;
      case 85: p = 850e2; z += 1000; break;
      case 70: p = 700e2; z += 3000; break;
      case 50: p = 500e2; z *= 10; break;
      case 40: p = 400e2; z *= 10; break;
      case 30: p = 300e2; z *= 10; break;
      case 25: p = 250e2; z = z*10 + 10000; break;
      case 20: p = 200e2; z = z*10 + 10000; break;
      case 15: p = 150e2; z = z*10 + 10000; break;
      case 10: p = 100e2; z = z*10 + 10000; break;
   }

}

Integer
Ttaa::get_number_of_standard_levels () const
{

   const char highest_wind_code = get_highest_wind_code ();
   switch (highest_wind_code)
   {
      case '0': return 1;
      case '1': return 11;
      case '2': return 9;
      case '3': return 7;
      case '4': return 6;
      case '5': return 5;
      case '7': return 4;
      case '8': return 3;
      case '9': return 2;
      case '/': return 0;
      default: throw Exception ("TEMP Parse Error");
   }

}

Ttaa::Ttaa (const string& message,
            const Dtime& time_hint)
:   Ttac (message, time_hint)
{
}

void
Ttaa::parse_to (Sounding& sounding) const
{

   Integer index = 3;
   if (size () < 4) { return; }

   parse_special_level (sounding, index); // Station Level
   parse_standard_levels (sounding, index);
   //parse_special_level (sounding, index); // Tropopause Level

}

Real
Ttbb::parse_nnppp (const string& nnppp) const
{
   const string& ppp = nnppp.substr (2, 3);
   if (ppp == "///") { return GSL_NAN; }
   else
   {
      Real pressure = atof (ppp.c_str ()) * 1e2;
      if (pressure < 100e2) { pressure += 1000e2; }
      return pressure;
   }
}

Ttbb::Ttbb (const string& message,
            const Dtime& time_hint)
   : Ttbd (message, time_hint)
{
}

void
Ttbb::parse_to (Sounding& sounding) const
{

   Integer index = 3;
   if (size () < 4) { return; }

   // if not reaching end of message by completion of
   //    parse_significant_temperature_levels,
   //    proceed to parse_significant_wind_levels
   if (!parse_significant_temperature_levels (sounding, index))
   {
      parse_significant_wind_levels (sounding, index);
   }

}

void
Ttcc::interpret_p_z (const Integer pp,
                     Real& pressure,
                     Real& geopotential_height) const
{

   Real& p = pressure;
   Real& z = geopotential_height;

   switch (pp)
   {
      case 70: p = 70e2; z = z*10 + 10000; break;
      case 50: p = 50e2; z = z*10 + 20000; break;
      case 30: p = 30e2; z = z*10 + 20000; break;
      case 20: p = 20e2; z = z*10 + 20000; break;
      case 10: p = 10e2; z = z*10 + 20000; break;
   }

}

Integer
Ttcc::get_number_of_standard_levels () const
{
   const char highest_wind_code = get_highest_wind_code ();
   switch (highest_wind_code)
   {
      case '1': return 5;
      case '2': return 4;
      case '3': return 3;
      case '5': return 2;
      case '7': return 1;
      case '/': return 0;
      default: throw Exception ("TEMP Parse Error");
   }
}

Ttcc::Ttcc (const string& message,
            const Dtime& time_hint)
   : Ttac (message, time_hint)
{
}

void
Ttcc::parse_to (Sounding& sounding) const
{

   Integer index = 3;
   if (size () < 4) { return; }
   if (((*this)[index]) == "NIL") { return; }

   parse_standard_levels (sounding, index);
   //parse_special_level (sounding, index); // Tropopause Level

}

Real
Ttdd::parse_nnppp (const string& nnppp) const
{
   const string& ppp = nnppp.substr (2, 3);
   if (ppp == "///") { return GSL_NAN; }
   else { return atof (ppp.c_str ()) * 10; }
}

Ttdd::Ttdd (const string& message,
            const Dtime& time_hint)
   : Ttbd (message, time_hint)
{
}

void
Ttdd::parse_to (Sounding& sounding) const
{

   Integer index = 3;
   if (size () < 4) { return; }
   if (((*this)[index]) == "NIL") { return; }

   // if not reaching end of message by completionof
   //    parse_significant_temperature_levels, proceed to
   //    parse_significant_wind_levels
   if (!parse_significant_temperature_levels (sounding, index))
   {
      parse_significant_wind_levels (sounding, index);
   }
}

namespace denise
{

   ostream&
   operator << (ostream &out_file,
                const Gts& gts)
   {
      for (Gts::const_iterator iterator = gts.begin ();
           iterator != gts.end (); iterator++)
      {
         const string& group = *(iterator);
         out_file << group << " ";
      }
      return out_file;
   }

}


