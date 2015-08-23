//
// tc.cc
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

#include <denise/gzstream.h>
#include <denise/tc.h>

using namespace std;
using namespace Cairo;
using namespace denise;

Tc_Track::Tc_Track (const Dstring& name,
                    const Dtime& dtime)
   : Track (dtime),
     name (name)
{
}

const Dstring&
Tc_Track::get_name () const
{
   return name;
}

Best_Tracks::Best_Tracks ()
{
}

void
Best_Tracks::ingest_jma (const Dstring& file_path)
{

   igzstream file (file_path);
   Best_Tracks::iterator i = end ();

   for (Dstring is; getline (file, is); )
   {

      if (is.substr (0, 5) == "66666")
      {
         if (i != end ()) { i->second.okay (); }
         const bool post_49 = stoi (is.substr (6, 2)) > 49;
         const Dstring& id = (post_49 ? "19" : "20") + is.substr (6, 4);
         const Dstring& name = Dstring (is.substr (30, 20)).get_trimmed ();
         i = insert (make_pair (id, Tc_Track (name))).first;
         continue;
      }

      if (i != end ())
      {

         Tc_Track& tc_track = i->second;

         const bool post_49 = stoi (is.substr (0, 2)) > 49;
         const Dtime dtime ((post_49 ? "19" : "20") + is.substr (0, 2));
         const Real latitude = stof (is.substr (15, 3)) * 0.1;
         const Real longitude = stof (is.substr (19, 4)) * 0.1;
         const Real pressure = stof (is.substr (24, 4));
         const Real max_wind = stof (is.substr (33, 3));

         tc_track.add (dtime.t, Lat_Long (latitude, longitude));
         tc_track.add (dtime.t, "pressure", pressure);
         tc_track.add (dtime.t, "max_wind", max_wind);

      }

   }

   file.close ();

}

set<Dstring>
Best_Tracks::get_subset (const Integer year) const
{

   set<Dstring> subset;

   for (auto iterator = begin (); iterator != end (); iterator++)
   {
      const Dstring& id = iterator->first;
      const Tc_Track& tc_track = iterator->second;
      const Integer y = tc_track.get_start_time ().get_year ();
      if (y == year) { subset.insert (id); }
   }

   return subset;

}

set<Dstring>
Best_Tracks::get_subset (const Domain_2D& domain_2d,
                         const Real dt) const
{

   set<Dstring> subset;

   for (auto iterator = begin (); iterator != end (); iterator++)
   {
      const Dstring& id = iterator->first;
      const Tc_Track& tc_track = iterator->second;
      if (!tc_track.trespass (domain_2d, dt)) { continue;}
      subset.insert (id);
   }

   return subset;

}

set<Dstring>
Best_Tracks::get_subset (const Integer day_of_year,
                         const Integer delta_days,
                         const Domain_2D& domain_2d,
                         const Real dt) const
{

   set<Dstring> subset;
   const Integer sjd = (day_of_year - delta_days);
   const Integer ejd = (day_of_year + delta_days);

   for (auto iterator = begin (); iterator != end (); iterator++)
   {

      const Dstring& id = iterator->first;
      const Tc_Track& tc_track = iterator->second;

      const Integer s = tc_track.get_start_time ().get_day_of_year ();
      const Integer e = tc_track.get_end_time ().get_day_of_year ();

      if ((sjd < s && ejd < s) || (sjd > e && ejd < e)) { continue; }
      if (!tc_track.trespass (domain_2d, dt)) { continue;}

      subset.insert (id);

   }

   return subset;

}

Real
Forecast::get_pressure (const Dstring& pressure_string)
{

   //const Reg_Exp mmhg ("MMHG");
   //const Reg_Exp hpa ("(HPA)|(MB)");

   Real pressure = stof (pressure_string) * 1e2;
   return pressure;

}

Real
Forecast::get_max_wind (const Dstring& max_wind_string)
{

   const Reg_Exp ms ("(MS)|(M/S)");
   const Reg_Exp kph ("(KPH)|(KM/H)");
   //const Reg_Exp knot ("(KT)|(KNOT)|(KNOTS)");

   Real max_wind = stof (max_wind_string);

   if (ms.match (max_wind_string))
   {
      max_wind *= 1.943844492;
   }
   else
   if (kph.match (max_wind_string))
   {
      max_wind /= 1.852;
   }

   return max_wind;

}

Forecast::Forecast (const Lat_Long& lat_long,
                    const Real pressure,
                    const Real max_wind)
   : lat_long (lat_long),
     pressure (pressure),
     max_wind (max_wind)
{
}

Forecast::Forecast (const Lat_Long& lat_long,
                    const Dstring& pressure_string,
                    const Dstring& max_wind_string)
   : lat_long (lat_long),
     pressure (get_pressure (pressure_string)),
     max_wind (get_max_wind (max_wind_string))
{
}

Dtime
Advisory::get_time (const Dstring& time_string,
                    const Dtime& time_stamp) const
{

   const Integer d = stoi (time_string.substr (0, 2));
   const Integer h = stoi (time_string.substr (2, 2));

   Integer y_ = time_stamp.get_year ();
   Integer m_ = time_stamp.get_month ();
   Integer d_ = time_stamp.get_day ();

   const bool last_day_of_month =
      (m_ == 1 && d_ == 31) ||
      (m_ == 2 && d_ == 28) ||
      (m_ == 3 && d_ == 31) ||
      (m_ == 4 && d_ == 30) ||
      (m_ == 5 && d_ == 31) ||
      (m_ == 6 && d_ == 30) ||
      (m_ == 7 && d_ == 31) ||
      (m_ == 8 && d_ == 31) ||
      (m_ == 9 && d_ == 30) ||
      (m_ == 10 && d_ == 31) ||
      (m_ == 11 && d_ == 30) ||
      (m_ == 12 && d_ == 31);

   const bool next_month = (d < d_) && (last_day_of_month);
   const bool next_year = (!next_month) && (m_ == 12);

   if (next_year) { y_++; }
   if (next_month) { m_ += (m_ == 12 ? -11 : 1); }

   const Dstring& yyyymm = Dstring::render ("%04d%02d", y_, m_);
   return Dtime (yyyymm + time_string.substr (0, 4));


}

void
Advisory::make_track ()
{

   track.set_dtime (initial_time);

   for (Advisory::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Real tau = Real (iterator->first);
      const Lat_Long& lat_long = iterator->second.lat_long;
      track.add (tau, lat_long);
   }

   track.okay ();

}

Advisory::Advisory (const Dstring& forecast_centre)
   : forecast_centre (forecast_centre),
     icon_string (forecast_centre)
{
}

Advisory::Advisory (const Dstring& forecast_centre,
                    const Dstring& icon_string)
   : forecast_centre (forecast_centre),
     icon_string (icon_string)
{
}

Dstring
Advisory::get_key () const
{
   const Lat_Long& ll = begin ()->second.lat_long;
   const Real& latitude = ll.latitude;
   const Real& longitude = ll.longitude;
   const Dstring& ll_str = Dstring::render ("%.1f %.1f", latitude, longitude);
   const Dstring& initial_time_str = initial_time.get_string ();
   return forecast_centre + " " + initial_time_str + " " + ll_str;
}

bool
Advisory::operator == (const Advisory& advisory) const
{
   return (get_key () == advisory.get_key ());
}

bool
Advisory::operator > (const Advisory& advisory) const
{
   return (get_key () > advisory.get_key ());
}

bool
Advisory::operator < (const Advisory& advisory) const
{
   return (get_key () < advisory.get_key ());
}

Lat_Long
Advisory::get_lat_long (const Real tau) const
{
   if (size () == 1)
   {
      const Integer i_tau = Integer (round (tau));
      Advisory::const_iterator iterator = find (i_tau);
      if (iterator == end ()) { return Lat_Long (GSL_NAN, GSL_NAN); }
      const Forecast& forecast = iterator->second;
      return forecast.lat_long;
   }
   return track.get_lat_long (tau);
}

void
Advisory::cairo (const RefPtr<Context> cr,
                 const Geodetic_Transform& transform,
                 const Real intensity) const
{

   const Real start_tau = track.begin ()->first;
   const Real end_tau = track.rbegin ()->first;
   const Integer n = Integer (round (end_tau - start_tau)) + 1;

   Real font_size;
   Real alpha;
   for (Integer i = 0; i < n; i++)
   {

      if (i % 24 == 0)
      {
         font_size = 24;
         alpha = 1;
      }
      else
      if (i % 6 == 0)
      {
         font_size = 12;
         alpha = 0.5;
      }
      else
      {
         font_size = 6;
         alpha = 0.2;
      }

      const Real tau = Real (i);
      const Lat_Long& lat_long = get_lat_long (tau);
      const Point_2D& point_2d = transform.transform (lat_long);

      const Real hue = (tau - 24) / 144;
      const Color& color = Color::hsb (hue, 1, 1, alpha * intensity);

      cr->set_font_size (font_size);
      cr->set_source_rgba (color.r, color.g, color.b, color.a);
      Label (icon_string, point_2d, 'c', 'c').cairo (cr, true);

   }

}

ostream&
Advisory::operator << (ostream &out_file) const
{

   ostream& o = out_file;

   o << "forecast_centre " << forecast_centre << endl;
   o << "tc_name " << tc_name << endl;
   o << "tc_id " << tc_id << endl;
   o << "initial_time_string " << initial_time_string << endl;

   for (auto iterator = begin (); iterator != end (); iterator++)
   {

      const Real tau = iterator->first;
      const Forecast& forecast = iterator->second;

      o << tau << " " << forecast.lat_long << " ";
      o << forecast.pressure << " " << forecast.max_wind << endl;

   }

   return o;

}

Wtss20_Vhhh::Wtss20_Vhhh (const Tokens& content,
                          const Dtime& time_stamp)
   : Advisory ("VHHH", "港")
//   : Advisory ("VHHH", "HK")
{

   const Reg_Exp blank ("^ *$");

   const Reg_Exp analysis ("WITH CENTRAL PRESSURE");
   const Reg_Exp forecast ("^FORECAST POSITION AND INTENSITY AT");
   const Reg_Exp max_wind_analysis ("^MAXIMUM WINDS NEAR THE CENTRE");

   const Reg_Exp tcid_re ("\\([0-9]...\\)", true);
   const Reg_Exp utc_re ("[0-9]..... UTC", true);
   const Reg_Exp knots_re ("([0-9].|[0-9]..) KNOTS", true);
   const Reg_Exp latitude_re ("([0-9]|[0-9].)\\.[0-9] (N|S)", true);
   const Reg_Exp longitude_re ("([0-9]|[0-9].|[0-9]..)\\.[0-9] (E|W)", true);
   const Reg_Exp hpa_re ("(9[0-9].|1[0-9]..) HECTOPASCALS", true);
   const Reg_Exp dissipate_re ("DISSIPATED", true);
   const Reg_Exp extratropical_re ("EXTRATROPICA", true);

   Dstring paragraph;
   Tokens paragraphs;

   // Grouping into paragraphs
   for (Tokens::const_iterator iterator = content.begin ();
        iterator != content.end (); iterator++)
   {

      if (iterator == content.begin ()) { continue; }

      const Dstring& line = *(iterator);
      const bool is_blank = blank.match (line);

      if (is_blank)
      {
         if (paragraph.size () != 0) { paragraphs.push_back (paragraph); }
         paragraph.clear ();
      }
      else
      {
         paragraph += line + " ";
      }

   }

   for (Tokens::const_iterator iterator = paragraphs.begin ();
        iterator != paragraphs.end (); iterator++)
   {

      const Dstring& paragraph = *(iterator);

      // First paragrah, Initial Time, Analysis Position
      if (analysis.match (paragraph))
      {

         const Integer utc_i = utc_re.get_match (paragraph).first;
         const Integer latitude_i = latitude_re.get_match (paragraph).first;
         const Integer longitude_i = longitude_re.get_match (paragraph).first;
         const Integer hpa_i = hpa_re.get_match (paragraph).first;
         //const Integer tcid_i = tcid_re.get_match (paragraph).first;

         const Real latitude = stof (paragraph.substr (latitude_i));
         const Real longitude = stof (paragraph.substr (longitude_i));
         const Real hpa = stof (paragraph.substr (hpa_i)) * 1e2;

         this->initial_time_string = paragraph.substr (utc_i, 6);
         this->initial_time = get_time (initial_time_string, time_stamp);
         //this->tc_id = paragraph.substr (tcid_i + 1, 4);

         const Tokens tokens (paragraph);
         for (Integer i = 0; i < tokens.size (); i++)
         {
            if (tcid_re.match (tokens[i]))
            {
               this->tc_id = tokens[i].substr (1, 4);
               this->tc_name = tokens[i - 1];
               break;
            }
         }

         const Lat_Long ll (latitude, longitude);
         (*this)[0] = Forecast (ll, hpa, 0);

         continue;

      }

      // Max Wind Analysis
      if (max_wind_analysis.match (paragraph))
      {
         const Integer knots_i = knots_re.get_match (paragraph).first;
         const Real knots = stof (paragraph.substr (knots_i));

         Forecast& forecast = begin ()->second;
         forecast.max_wind = knots;
         continue;
      }

      // Forecast positions
      if (forecast.match (paragraph))
      {

         // Do nothing for dissipated TC
         if (dissipate_re.match (paragraph)) { continue; }

         // Do nothing for TC becoming extratropical
         if (extratropical_re.match (paragraph)) { continue; }

         const Integer utc_i = utc_re.get_match (paragraph).first;
         const Integer knots_i = knots_re.get_match (paragraph).first;
         const Integer latitude_i = latitude_re.get_match (paragraph).first;
         const Integer longitude_i = longitude_re.get_match (paragraph).first;

         const Dstring utc = paragraph.substr (utc_i, 6);
         const Real knots = stof (paragraph.substr (knots_i));
         const Real latitude = stof (paragraph.substr (latitude_i));
         const Real longitude = stof (paragraph.substr (longitude_i));

         const Dstring& its = this->initial_time_string;
         const Dstring& its_dd = its.substr (0, 2);
         const Dstring& its_hh = its.substr (2, 2);
         const Integer initial_dd = stoi (its.substr (0, 2));
         const Integer initial_hh = stoi (its.substr (2, 2));
         const Integer utc_dd = stoi (utc.substr (0, 2));
         const Integer utc_hh = stoi (utc.substr (2, 2));
         const Integer dd = utc_dd - initial_dd;
         const Integer hh = utc_hh - initial_hh;

         const Dtime& forecast_time = get_time (utc, time_stamp);
         const Integer tau = Integer (round (forecast_time.t - initial_time.t));

         const Lat_Long lat_long (latitude, longitude);
         (*this)[tau] = Forecast (lat_long, GSL_NAN, knots);
         continue;

      }

   }

}

Wtpq20_Babj::Wtpq20_Babj (const Tokens& content,
                          const Dtime& time_stamp)
   : Advisory ("BABJ", "中")
//   : Advisory ("BABJ", "BJ")
{

   const Reg_Exp reg_exp ("P\\+[0-9].HR");
   const Reg_Exp time_re ("INITIAL +TIME");
   const Reg_Exp utc_re ("[0-9]..... *UTC", true);

   for (Tokens::const_iterator iterator = content.begin ();
        iterator != content.end (); iterator++)
   {

      const Dstring& input_string = *(iterator);
      const Tokens tokens (input_string);
      const Integer n = tokens.size ();

      if (time_re.match (input_string))
      {
         this->tc_name = tokens[1];
         this->tc_id = tokens[2];
         const Integer utc_i = utc_re.get_match (input_string).first;
         this->initial_time_string = input_string.substr (utc_i, 6);
         this->initial_time = get_time (initial_time_string, time_stamp);
         continue;
      }

      if (n > 0 && tokens[0] == "00HR")
      {
         const Integer tau = 0;
         const Lat_Long ll (tokens[1], tokens[2]);
         (*this)[0] = Forecast (ll, tokens[3], tokens[4]);
         continue;
      }

      if (n > 0 && reg_exp.match (tokens[0]))
      {
         const Integer tau = stoi (tokens[0].substr (2));
         const Lat_Long ll (tokens[1], tokens[2]);
         (*this)[tau] = Forecast (ll, tokens[3], tokens[4]);
         continue;
      }

   }

}

Wtpq20_Rjtd::Wtpq20_Rjtd (const Tokens& content,
                          const Dtime& time_stamp)
   : Advisory ("RJTD", "日")
//   : Advisory ("RJTD", "JMA")
{

   const Reg_Exp name_re ("^NAME");
   const Reg_Exp pstn_re ("^PSTN");
   const Reg_Exp hrfc_re ("^[0-9].HF");
   const Reg_Exp utc_re ("[0-9].....UTC", true);
   const Reg_Exp tcid_re ("\\([0-9]...\\)");

   for (Tokens::const_iterator iterator = content.begin ();
        iterator != content.end (); iterator++)
   {

      const Dstring& input_string = *(iterator);
      const Tokens tokens (input_string);
      const Integer n = tokens.size ();

      if (name_re.match (input_string))
      {
         for (Integer i = 0; i < n; i++)
         {
            if (tcid_re.match (tokens[i]))
            {
               this->tc_id = tokens[i].substr (1, 4);
               this->tc_name = tokens[i - 1];
               break;
            }
         }
         continue;
      }

      if (pstn_re.match (input_string))
      {
         const Integer utc_i = utc_re.get_match (input_string).first;
         this->initial_time_string = input_string.substr (utc_i, 6);
         this->initial_time = get_time (initial_time_string, time_stamp);
         const Integer tau = 0;
         const Lat_Long ll (tokens[2], tokens[3]);
         (*this)[0] = Forecast (ll, GSL_NAN, GSL_NAN);
         continue;
      }

      if (hrfc_re.match (input_string))
      {

         const Integer utc_i = utc_re.get_match (input_string).first;
         const Dstring utc = input_string.substr (utc_i, 6);

         const Dstring& its = initial_time_string;
         const Integer initial_dd = stoi (its.substr (0, 2));
         const Integer initial_hh = stoi (its.substr (2, 2));
         const Integer utc_dd = stoi (utc.substr (0, 2));
         const Integer utc_hh = stoi (utc.substr (2, 2));
         const Integer dd = utc_dd - initial_dd;
         const Integer hh = utc_hh - initial_hh;

         const Dtime& forecast_time = get_time (utc, time_stamp);
         const Integer tau = Integer (round (forecast_time.t - initial_time.t));

         const Lat_Long ll (tokens[2], tokens[3]);
         (*this)[tau] = Forecast (ll, GSL_NAN, GSL_NAN);

         continue;

      }

      if (n > 0 && tokens[0] == "PRES")
      {
         const Real pres = stof (tokens[1]) * 1e2;
         rbegin ()->second.pressure = pres;
         continue;
      }

      if (n > 0 && tokens[0] == "MXWD")
      {
         const Real mxwd = stof (tokens[1]);
         rbegin ()->second.max_wind = mxwd;
         continue;
      }

   }

}

Wtko20_Rksl::Wtko20_Rksl (const Tokens& content,
                          const Dtime& time_stamp)
   : Advisory ("RKS", "韓")
//   : Advisory ("RKS", "KMA")
{

   const Reg_Exp name_re ("^NAME");
   const Reg_Exp position_re ("^POSITION");
   const Reg_Exp hrfc_re ("^[0-9].HF");
   const Reg_Exp utc_re ("[0-9].....UTC", true);
   const Reg_Exp tcid_re ("[0-9]...", true);

   for (Tokens::const_iterator iterator = content.begin ();
        iterator != content.end (); iterator++)
   {

      const Dstring& input_string = *(iterator);
      const Tokens tokens (input_string);
      const Integer n = tokens.size ();

      if (name_re.match (input_string))
      {
         for (Integer i = 0; i < n; i++)
         {
            if (tcid_re.match (tokens[i]))
            {
               this->tc_id = tokens[i];
               if (i <= n - 1) { this->tc_name = tokens[i + 1]; }
               break;
            }
         }
         continue;
      }

      if (position_re.match (input_string))
      {

         const Integer utc_i = utc_re.get_match (input_string).first;
         const Dstring utc = input_string.substr (utc_i, 6);

         if (size () == 0)
         {
            this->initial_time_string = utc;
            this->initial_time = get_time (initial_time_string, time_stamp);
         }

         const Dstring& its = initial_time_string;
         const Integer initial_dd = stoi (its.substr (0, 2));
         const Integer initial_hh = stoi (its.substr (2, 2));
         const Integer utc_dd = stoi (utc.substr (0, 2));
         const Integer utc_hh = stoi (utc.substr (2, 2));
         const Integer dd = utc_dd - initial_dd;
         const Integer hh = utc_hh - initial_hh;

         const Dtime& forecast_time = get_time (utc, time_stamp);
         const Integer tau = Integer (round (forecast_time.t - initial_time.t));

         const Lat_Long ll (tokens[2], tokens[3]);
         (*this)[tau] = Forecast (ll, GSL_NAN, GSL_NAN);

         continue;

      }

      if (n > 0 && tokens[0] == "PRES")
      {
         const Real pres = stof (tokens[1]) * 1e2;
         rbegin ()->second.pressure = pres;
         continue;
      }

      if (n > 0 && tokens[0] == "PRES/VMAX")
      {
         const Real pres = stof (tokens[1]);
         const Real vmax = stof (tokens[2]);
         rbegin ()->second.pressure = pres;
         rbegin ()->second.max_wind = vmax;
         continue;
      }

   }

}

Wtpn31_Pgtw::Wtpn31_Pgtw (const Tokens& content,
                          const Dtime& time_stamp)
   : Advisory ("PGTW", "聯")
//   : Advisory ("PGTW", "JTWC")
{

   Integer tau;
   const Reg_Exp initial_re ("^1\\. ");
   const Reg_Exp tc_id_re ("[0-9].[EWC]", true);
   const Reg_Exp tc_name_re ("\\([A-Z\\-]+\\)", true);
   const Reg_Exp tau_re ("^ *[0-9]+ HRS, VALID AT:");
   const Reg_Exp analysis_re ("[0-9].....Z --- NEAR ([0-9]|[0-9].)\\.[0-9](N|S)");
   const Reg_Exp forecast_re ("[0-9].....Z --- ([0-9]|[0-9].)\\.[0-9](N|S)");

   cout << "PGTW " << time_stamp.get_string () << endl;

   for (Tokens::const_iterator iterator = content.begin ();
        iterator != content.end (); iterator++)
   {

      const Dstring& input_string = *(iterator);

      if (input_string.substr (0, 7) == "REMARKS") { break; }

      const Tokens tokens (input_string);
      const Integer n = tokens.size ();

      if (initial_re.match (input_string))
      {
         const Integer tc_id_i = tc_id_re.get_match (input_string).first;
         const Integer tc_id_j = tc_id_re.get_match (input_string).second;
         const Integer tc_name_i = tc_name_re.get_match (input_string).first;
         const Integer tc_name_j = tc_name_re.get_match (input_string).second;
         const Integer tc_id_n = tc_id_j - tc_id_i + 1;
         const Integer tc_name_n = tc_name_j - tc_name_i - 2;
         this->tc_id = input_string.substr (tc_id_i, tc_id_n);
         this->tc_name = input_string.substr (tc_name_i + 1, tc_name_n);
         continue;
      }

      if (analysis_re.match (input_string))
      {

         const Dstring utc = tokens[0].substr (0, 6);
         this->initial_time_string = utc;
         this->initial_time = get_time (initial_time_string, time_stamp);

         cout << "PGTW analysis " << time_stamp.get_string () << " " << initial_time_string << " " << initial_time.get_string () << endl;
         const Lat_Long ll (tokens[3], tokens[4]);
         (*this)[0] = Forecast (ll, GSL_NAN, GSL_NAN);

         continue;

      }

      if (tau_re.match (input_string))
      {
         tau = stoi (tokens[0]);
         continue;
      }

      if (forecast_re.match (input_string))
      {

         const Dstring& utc = tokens[0].substr (0, 6);
         const Lat_Long ll (tokens[2], tokens[3]);

         //const Integer initial_dd = stoi (initial_time_string.substr (0, 2));
         //const Integer initial_hh = stoi (initial_time_string.substr (2, 2));
         //const Integer utc_dd = stoi (utc.substr (0, 2));
         //const Integer utc_hh = stoi (utc.substr (2, 2));
         //const Integer dd = utc_dd - initial_dd;
         //const Integer hh = utc_hh - initial_hh;
         //const Integer tau = (dd < 0 ? -1 : dd) * 24 + hh;

         (*this)[tau] = Forecast (ll, GSL_NAN, GSL_NAN);
         continue;

      }

      if (input_string.substr (0, 25) == "   MAX SUSTAINED WINDS - ")
      {
         const Real max_wind = stof (tokens[4]);
         rbegin ()->second.max_wind = max_wind;
         continue;
      }

   }

}

Real
Wtph20_Rpmm::get_real (const Dstring& digit)
{
   if (digit == "ONE") { return 1; }
   if (digit == "TWO") { return 2; }
   if (digit == "THREE") { return 3; }
   if (digit == "FOUR") { return 4; }
   if (digit == "FIVE") { return 5; }
   if (digit == "SIX") { return 6; }
   if (digit == "SEVEN") { return 7; }
   if (digit == "EIGHT") { return 8; }
   if (digit == "NINE") { return 9; }
   return 0;
}

Lat_Long
Wtph20_Rpmm::get_lat_long (const Dstring& lat_long_string)
{

   const Tokens tokens (lat_long_string);

   const Real latitude_2 = get_real (tokens[0]) * 10;
   const Real latitude_1 = get_real (tokens[1]);
   const Real latitude_0 = get_real (tokens[3]) * 0.1;


   const Real longitude_3 = get_real (tokens[5]) * 100;
   const Real longitude_2 = get_real (tokens[6]) * 10;
   const Real longitude_1 = get_real (tokens[7]);
   const Real longitude_0 = get_real (tokens[9]) * 0.1;

   Real latitude = latitude_2 + latitude_1 + latitude_0;
   if (tokens[4] == "SOUTH") { latitude *= -1; }

   Real longitude = longitude_3 + longitude_2 + longitude_1 + longitude_0;
   if (tokens[10] == "WEST") { longitude *= -1; }

   return Lat_Long (latitude, longitude);

}

Real
Wtph20_Rpmm::get_central_pressure (const Tokens& tokens)
{

   switch (tokens.size ())
   {

      case 6:
      {
         const Real p2 = get_real (tokens[2]) * 100;
         const Real p1 = get_real (tokens[3]) * 10;
         const Real p0 = get_real (tokens[4]);
         return p2 + p1 + p0;
      }

      case 7:
      {
         const Real p3 = get_real (tokens[2]) * 1000;
         const Real p2 = get_real (tokens[3]) * 100;
         const Real p1 = get_real (tokens[4]) * 10;
         const Real p0 = get_real (tokens[5]) * 1;
         return p3 + p2 + p1 + p0;
      }

   }

   return GSL_NAN;

}

Real
Wtph20_Rpmm::get_meters_per_second (const Tokens& tokens)
{

   switch (tokens.size ())
   {

      case 7:
      {
         const Real p1 = get_real (tokens[2]) * 10;
         const Real p0 = get_real (tokens[3]);
         return p1 + p0;
      }

      case 8:
      {
         const Real p2 = get_real (tokens[2]) * 100;
         const Real p1 = get_real (tokens[3]) * 10;
         const Real p0 = get_real (tokens[4]);
         return p2 + p1 + p0;
      }

   }

   return GSL_NAN;

}

Wtph20_Rpmm::Wtph20_Rpmm (const Tokens& content,
                          const Dtime& time_stamp)
   : Advisory ("RPMM", "菲")
//   : Advisory ("RPMM", "PAG")
{

   const Dstring d ("(ZERO|ONE|TWO|THREE|FOUR|FIVE|SIX|SEVEN|EIGHT|NINE)");
   const Dstring m ("(JANUARY|FEBRUARY|MARCH|APRIL|MAY|JUNE|JULY|AUGUST|SEPTEMBER|OCTOBER|NOVEMBER|DECEMBER)");

   const Dstring d2 (d + " +" + d);
   const Dstring d3 (d + " +" + d + " +" + d);
   const Dstring d4 (d + " +" + d + " +" + d + " +" + d);
   const Dstring d2or3 ("(" + d2 + "|" + d3 + ")");
   const Dstring d3or4 ("(" + d3 + "|" + d4 + ")");
   const Dstring lat_str (d2 + " +POINT +" + d + " +(NORTH|SOUTH)");
   const Dstring long_str (d3 + " +POINT +" + d + " +(EAST|WEST)");
   const Dstring cp_str ("CENTRAL +PRESSURE +" + d3or4 + " +HECTOPASCALS?");
   const Dstring mw_str ("MAXIMUM +WINDS? +" + d2or3 + " +METERS? +PER +SECOND");

   Dstring paragraph;
   Tokens paragraphs;

   const Reg_Exp blank ("^ *$");
   const Reg_Exp initial_time_string_re ("^AT +[0-9]... +");

   const Reg_Exp tc_id_re ("(\\(|\\{)[0-9]...(\\)|\\})", true);
   const Reg_Exp tc_name_re ("\\([A-Z\\-]+\\)", true);
   const Reg_Exp ll_re (lat_str + " " + long_str, true);
   const Reg_Exp fcll_re ("AT [0-9]..... " + lat_str + " " + long_str, true);

   const Reg_Exp cp_re (cp_str, true);
   const Reg_Exp mw_re (mw_str, true);

   for (Tokens::const_iterator iterator = content.begin ();
        iterator != content.end (); iterator++)
   {

      if (iterator == content.begin ())
      {
         const Tokens tokens (*(iterator));
         this->initial_time_string = tokens[3];
         this->initial_time = get_time (initial_time_string, time_stamp);
         continue;
      }

      const Dstring& line = *(iterator);
      const bool is_blank = blank.match (line);

      if (is_blank)
      {
         if (paragraph.size () != 0) { paragraphs.push_back (paragraph); }
         paragraph.clear ();
      }
      else
      {
         paragraph += line + " ";
      }

   }

   for (Tokens::const_iterator iterator = paragraphs.begin ();
        iterator != paragraphs.end (); iterator++)
   {

      const Dstring& paragraph = *(iterator);
      if (!initial_time_string_re.match (paragraph)) { continue; }

      const Tokens tokens (paragraph);

      const Integer tc_id_i = tc_id_re.get_match (paragraph).first;
      const Integer tc_id_j = tc_id_re.get_match (paragraph).second;
      const Integer tc_id_n = tc_id_j - tc_id_i - 2;
      if (tc_id_i > 0)
      {
         this->tc_id = paragraph.substr (tc_id_i + 1, tc_id_n);
      }

      const Integer tc_name_i = tc_name_re.get_match (paragraph).first;
      const Integer tc_name_j = tc_name_re.get_match (paragraph).second;
      const Integer tc_name_n = tc_name_j - tc_name_i - 2;
      if (tc_name_i > 0)
      {
         this->tc_name = paragraph.substr (tc_name_i + 1, tc_name_n);
      }

      const Integer ll_i = ll_re.get_match (paragraph).first;
      const Integer ll_j = ll_re.get_match (paragraph).second;
      const Integer ll_n = ll_j - ll_i + 1;
      const Dstring lat_long_str = (paragraph.substr (ll_i, ll_n));
      const Lat_Long ll = get_lat_long (lat_long_str);
      (*this)[0] = Forecast (ll, GSL_NAN, GSL_NAN);

      const Integer cp_i = cp_re.get_match (paragraph).first;
      const Integer cp_j = cp_re.get_match (paragraph).second;
      const Integer cp_n = cp_j - cp_i + 1;
      const Dstring& cp_substr = paragraph.substr (cp_i, cp_n);
      const Tokens tokens_cp (cp_substr);
      (*this)[0].pressure = get_central_pressure (tokens_cp) * 1e2;

      const Integer mw_i = mw_re.get_match (paragraph).first;
      const Integer mw_j = mw_re.get_match (paragraph).second;
      const Integer mw_n = mw_j - mw_i + 1;
      const Dstring& mw_substr = paragraph.substr (mw_i, mw_n);
      const Tokens tokens_mw (mw_substr);
      (*this)[0].max_wind = get_meters_per_second (tokens_mw) * 1.943844492;

      Dstring para = paragraph;

      for (Iduple iduple = fcll_re.get_match (para); iduple.first >= 0;
           iduple = fcll_re.get_match (para))
      {

         const Integer fcll_i = iduple.first;
         const Integer fcll_j = iduple.second;
         const Integer fcll_n = fcll_j - fcll_i + 1;
         const Dstring fc_lat_long_str = (para.substr (fcll_i, fcll_n));
         const Dstring utc = fc_lat_long_str.substr (3, 6);
         const Lat_Long fcll = get_lat_long (fc_lat_long_str.substr (10));

         const Dstring& its = initial_time_string;
         const Integer initial_dd = stoi (its.substr (0, 2));
         const Integer initial_hh = stoi (its.substr (2, 2));
         const Integer utc_dd = stoi (utc.substr (0, 2));
         const Integer utc_hh = stoi (utc.substr (2, 2));
         const Integer dd = utc_dd - initial_dd;
         const Integer hh = utc_hh - initial_hh;

         const Dtime& forecast_time = get_time (utc, time_stamp);
         const Integer tau = Integer (round (forecast_time.t - initial_time.t));
         (*this)[tau] = Forecast (fcll, GSL_NAN, GSL_NAN);

         para = para.substr (fcll_j);

      }

   }

}

Wtnt80_Egrr::Wtnt80_Egrr (const Tokens& content,
                          const Dtime& time_stamp)
   : Advisory ("EGRR", "英")
//   : Advisory ("EGRR", "UK")
{

   const Dstring utc_str ("[0-9].UTC [0-9].\\.[0-9].\\.[0-9]...");
   const Dstring ll_str ("[0-9].\\.[0-9](N|S) [ 0-9][0-9].\\.[0-9](E|W)");

   const Reg_Exp ap_re ("ANALYSED POSITION");
   const Reg_Exp atcf_re ("ATCF IDENTIFIER");
   const Reg_Exp position_re ("^ " + utc_str + "  " + ll_str);

   Dtime initial_dtime;

   for (Tokens::const_iterator iterator = content.begin ();
        iterator != content.end (); iterator++)
   {

      const Dstring& input_string = *(iterator);
      const Tokens tokens (input_string);
      const Integer n = tokens.size ();

      if (ap_re.match (input_string))
      {
         this->tc_name = tokens[n - 6];
         continue;
      }

      if (atcf_re.match (input_string))
      {
         this->tc_id = tokens[3];
         continue;
      }

      if (position_re.match (input_string))
      {
         const Dstring time_string = tokens[0] + " " + tokens[1];
         const Dtime dtime (time_string, "%HUTC %d.%m.%Y");

         if (size () == 0)
         {
            initial_dtime.t = dtime.t;
            this->initial_time_string = dtime.get_string ("%d%H%M");
            this->initial_time = get_time (initial_time_string, time_stamp);
         }

         const Integer tau = Integer (round (dtime.t - initial_dtime.t));
         const Lat_Long ll (tokens[2], tokens[3]);
         (*this)[tau] = Forecast (ll, GSL_NAN, GSL_NAN);

         continue;

      }

   }

}

Wtth20_Vtbb::Wtth20_Vtbb (const Tokens& content,
                          const Dtime& time_stamp)
   : Advisory ("VTBB", "泰")
//   : Advisory ("VTBB", "TH")
{

   const Reg_Exp name_re ("^NAME");
   const Reg_Exp pstn_re ("^PSTN");
   const Reg_Exp hrfc_re ("^[0-9]. HR");
   const Reg_Exp utc_re ("[0-9]..... UTC", true);

   for (Tokens::const_iterator iterator = content.begin ();
        iterator != content.end (); iterator++)
   {

      const Dstring& input_string = *(iterator);
      const Tokens tokens (input_string);
      const Integer n = tokens.size ();

      if (name_re.match (input_string))
      {
         if (tokens.size () > 3)
         {
            this->tc_name = tokens[3];
            this->tc_name = tc_name.substr (1, tc_name.size () - 2);
         }
         if (tokens.size () > 2) { this->tc_id = tokens[2]; }
         continue;
      }

      if (pstn_re.match (input_string))
      {

         const Integer tau = 0;
         const Integer utc_i = utc_re.get_match (input_string).first;

         if (utc_i == -1)
         {
            const Tokens tokens (content[0]);
            this->initial_time_string = tokens[3];
            this->initial_time = get_time (initial_time_string, time_stamp);
            const Lat_Long ll (tokens[1], tokens[3]);
            (*this)[0] = Forecast (ll, GSL_NAN, GSL_NAN);
         }
         else
         {
            this->initial_time_string = input_string.substr (utc_i, 6);
            this->initial_time = get_time (initial_time_string, time_stamp);
            const Lat_Long ll (tokens[3], tokens[5]);
            (*this)[0] = Forecast (ll, GSL_NAN, GSL_NAN);
         }

         continue;
      }

      if (hrfc_re.match (input_string))
      {

         const Dstring& utc = tokens[2];

         if (utc_re.get_match (input_string).first != -1)
         {

            const Dstring& its = initial_time_string;
            const Integer initial_dd = stoi (its.substr (0, 2));
            const Integer initial_hh = stoi (its.substr (2, 2));
            const Integer utc_dd = stoi (utc.substr (0, 2));
            const Integer utc_hh = stoi (utc.substr (2, 2));
            const Integer dd = utc_dd - initial_dd;
            const Integer hh = utc_hh - initial_hh;

            const Dtime& forecast_time = get_time (utc, time_stamp);
            const Real forecast_t = forecast_time.t;
            const Real initial_t = initial_time.t;
            const Integer tau = Integer (round (forecast_t - initial_t));

            const Lat_Long ll (tokens[4], tokens[6]);
            (*this)[tau] = Forecast (ll, GSL_NAN, GSL_NAN);

         }
         else
         {

            const Integer tau = stoi (input_string.substr (0, 2));
            const Lat_Long ll (tokens[4], tokens[6]);
            (*this)[tau] = Forecast (ll, GSL_NAN, GSL_NAN);

         }


         continue;

      }

      if (n > 0 && tokens[0] == "PRES")
      {
         const Real pres = stof (tokens[1]) * 1e2;
         rbegin ()->second.pressure = pres;
         continue;
      }

      if (n > 0 && tokens[0] == "MAXD")
      {
         const Real maxd = stof (tokens[1]);
         rbegin ()->second.max_wind = maxd;
         continue;
      }

   }

}

Knhc_4::Knhc_4 (const Tokens& content,
                const Dtime& time_stamp)
   : Advisory ("KNHC", "颶")
//   : Advisory ("KNHC", "NHC")
{

   Integer tau;
   const Reg_Exp tc_id_re ("[0-9].W", true);
   const Reg_Exp tc_name_re ("\\([A-Z]*\\)", true);
   const Reg_Exp analysis_re ("^INIT  [0-9]./[0-9]...Z");
   const Reg_Exp forecast_re ("^[ 0-9]..H  [0-9]./[0-9]...Z");
   const Reg_Exp advisory_number_re ("(ADVISORY|DISCUSSION) NUMBER");
   const Reg_Exp latitude_re ("([0-9]|[0-9].)\\.[0-9](N|S)");
   const Reg_Exp longitude_re ("([0-9]|[0-9].|[0-9]..)\\.[0-9](E|W)");

   for (Tokens::const_iterator iterator = content.begin ();
        iterator != content.end (); iterator++)
   {

      const Dstring& input_string = *(iterator);
      const Tokens tokens (input_string);
      const Integer n = tokens.size ();

      if (advisory_number_re.match (input_string))
      {
         for (Integer i = 0; i < n; i++)
         {
            if (tokens[i] == "ADVISORY" || tokens[i] == "DISCUSSION")
            {
               this->tc_name = tokens[i - 1];
               break;
            }
         }
         continue;
      }

      if (analysis_re.match (input_string))
      {
         const Dstring& dd_string = input_string.substr (6, 2);
         const Dstring& hhmm_string = input_string.substr (9, 4);
         const Dstring initial_time_string = dd_string + hhmm_string;
         this->initial_time = get_time (initial_time_string, time_stamp);
         const Lat_Long ll (tokens[2], tokens[3]);
         const Real max_wind = stof (tokens[4]);
         (*this)[0] = Forecast (ll, GSL_NAN, max_wind);
         continue;
      }

      if (forecast_re.match (input_string))
      {
         if (n < 7) { continue; }
         if (!latitude_re.match (input_string)) { continue; }
         if (!longitude_re.match (input_string)) { continue; }

         const Dstring& dd_string = input_string.substr (6, 2);
         const Dstring& hhmm_string = input_string.substr (9, 4);
         const Dstring time_string = dd_string + hhmm_string;
         const Dtime dtime = get_time (time_string, time_stamp);
         const Integer tau = Integer (round (dtime.t - initial_time.t));
         const Lat_Long ll (tokens[2], tokens[3]);
         const Real max_wind = stof (tokens[4]);
         (*this)[tau] = Forecast (ll, GSL_NAN, max_wind);
         continue;
      }

   }

}

Advisory_Store::Cluster_Info::Cluster_Info (const Advisory_Store& advisory_store,
                                            const Tokens& key_tokens,
                                            const Dataset& dataset,
                                            const Real cluster_distance)
   : advisory_store (advisory_store),
     key_tokens (key_tokens)
{
   cluster_ptr = new Cluster (dataset);
   cluster_multimap = cluster_ptr->get_multimap (cluster_distance);
}

Advisory_Store::Cluster_Info::~Cluster_Info ()
{
   delete cluster_ptr;
}

set<Integer>
Advisory_Store::Cluster_Info::get_cluster_index_set () const
{
   return cluster_multimap.get_cluster_set ();
}

Dstring
Advisory_Store::Cluster_Info::get_tc_name (const Integer cluster_index) const
{

   Dstring tc_name;
   set<Dstring> tc_name_set;

   typedef denise::Cluster::Multimap::const_iterator Iterator;
   typedef pair<Iterator, Iterator> Range;
   const Range& range = cluster_multimap.equal_range (cluster_index);

   for (Iterator j = range.first; j != range.second; j++)
   {
      const Integer leaf_index = j->second;
      const Dstring& key = key_tokens[leaf_index];
      const Advisory& advisory = advisory_store.at (key);
      const Dstring& tc_name = advisory.tc_name;
      tc_name_set.insert (tc_name);
   }

   for (set<Dstring>::const_iterator iterator = tc_name_set.begin ();
        iterator != tc_name_set.end (); iterator++)
   {
      if (tc_name != "") { tc_name += "/"; }
      tc_name += *(iterator);
   }

   return tc_name;

}

Lat_Long
Advisory_Store::Cluster_Info::get_lat_long (const Integer cluster_index) const
{

   Lat_Long lat_long;

   typedef denise::Cluster::Multimap::const_iterator Iterator;
   typedef pair<Iterator, Iterator> Range;
   const Range& range = cluster_multimap.equal_range (cluster_index);

   for (Iterator j = range.first; j != range.second; j++)
   {
      const Integer leaf_index = j->second;
      const Dstring& key = key_tokens[leaf_index];
      const Advisory& advisory = advisory_store.at (key);
      const Lat_Long& ll = advisory.get_lat_long (0);
      lat_long.latitude += ll.latitude;
      lat_long.longitude += ll.longitude;
   }

   const Integer n = std::distance (range.first, range.second);
   lat_long.latitude /= n;
   lat_long.longitude /= n;

   return lat_long;

}

void
Advisory_Store::insert (const Advisory& advisory)
{

   const Dstring& key = advisory.get_key ();
   Advisory_Store::iterator iterator = find (key);
   if (iterator != end ()) { erase (iterator); } 

   map<Dstring, Advisory>::insert (make_pair (key, advisory));

}

void
Advisory_Store::parse (const Tokens& content,
                       const Dtime& time_stamp)
{

   const Tokens tokens (content[0]);

   if (((tokens[1].substr (0, 5) == "WTSS2") ||
        (tokens[1].substr (0, 5) == "WTPQ2")) &&
       tokens[2] == "VHHH")
   {

      Wtss20_Vhhh advisory (content, time_stamp);
      if (advisory.size () > 0)
      {
         const Dstring& key = advisory.get_key ();
         insert (advisory);
         at (key).make_track ();
      }

      return;

   }

   if (tokens[1].substr (0, 5) == "WTPQ2" && tokens[2] == "BABJ")
   {

      Wtpq20_Babj advisory (content, time_stamp);
      if (advisory.size () > 0)
      {
         const Dstring& key = advisory.get_key ();
         insert (advisory);
         at (key).make_track ();
      }

      return;

   }

   if (tokens[1].substr (0, 5) == "WTPQ2" && tokens[2] == "RJTD")
   {

      Wtpq20_Rjtd advisory (content, time_stamp);
      if (advisory.size () > 0)
      {
         const Dstring& key = advisory.get_key ();
         insert (advisory);
         at (key).make_track ();
      }

      return;

   }

   if (tokens[1].substr (0, 5) == "WTKO2" && tokens[2] == "RKS")
   {

      Wtko20_Rksl advisory (content, time_stamp);
      if (advisory.size () > 0)
      {
         const Dstring& key = advisory.get_key ();
         insert (advisory);
         at (key).make_track ();
      }

      return;

   }

   if (tokens[1].substr (0, 5) == "WTPN3" && tokens[2] == "PGTW")
   {

      Wtpn31_Pgtw advisory (content, time_stamp);
      if (advisory.size () > 0)
      {
         const Dstring& key = advisory.get_key ();
         insert (advisory);
         at (key).make_track ();
      }

      return;

   }

   if (tokens[1].substr (0, 5) == "WTPH2" && tokens[2] == "RPMM")
   {

/*
      Wtph20_Rpmm advisory (content, time_stamp);
      if (advisory.size () > 0)
      {
         const Dstring& key = advisory.get_key ();
         insert (advisory);
         at (key).make_track ();
      }
*/

      return;

   }

   if (tokens[1].substr (0, 5) == "WTNT8" && tokens[2] == "EGRR")
   {

      return;
      Wtnt80_Egrr advisory (content, time_stamp);
      if (advisory.size () > 0)
      {
         const Dstring& key = advisory.get_key ();
         insert (advisory);
         at (key).make_track ();
      }

      return;

   }

   if (tokens[1].substr (0, 5) == "WTTH2" && tokens[2] == "VTBB")
   {

/*
      Wtth20_Vtbb advisory (content, time_stamp);
      if (advisory.size () > 0)
      {
         const Dstring& key = advisory.get_key ();
         insert (advisory);
         at (key).make_track ();
      }
*/

      return;

   }

   if ((tokens[1].substr (0, 5) == "WTNT4" && tokens[2] == "KNHC") ||
       (tokens[1].substr (0, 5) == "WTPZ4" && tokens[2] == "KNHC"))
   {

      Knhc_4 advisory (content, time_stamp);
      if (advisory.size () > 0)
      {
         const Dstring& key = advisory.get_key ();
         insert (advisory);
         at (key).make_track ();
      }

      return;

   }

}

Advisory_Store::Advisory_Store ()
{
}

void
Advisory_Store::ingest_dir (const Dstring& dir_path,
                            const Reg_Exp& file_reg_exp)
{

   const Tokens& dir_listing = get_dir_listing (dir_path, file_reg_exp);

   for (Tokens::const_iterator iterator = dir_listing.begin ();
        iterator != dir_listing.end (); iterator++)
   {
      const Dstring& file_name = *(iterator);
      const Dstring file_path = dir_path + "/" + file_name;
      ingest_file (file_path);
   }

}

void
Advisory_Store::ingest_file (const Dstring& file_path)
{

   string is;
   const string fp (file_path.begin (), file_path.end ());
   ifstream file (fp.c_str ());

   bool first = true;
   const Reg_Exp header ("^\\*. [A-Z]...[0-9]. [A-Z]... [0-9]..... ");

   const Reg_Exp yymmddhh ("[0-9].......");
   const Dstring file_name = Tokens (file_path, "/").back ();
   const bool file_name_is_time = yymmddhh.match (file_name);

   const Dtime time_stamp = (!file_name_is_time ?
      Dtime::last_synoptic () : Dtime (file_name, "%y%m%d%H"));

   Tokens content;

   while (getline (file, is))
   {

      const Dstring input_string (is);
      const bool is_header = header.match (input_string);

      if (is_header && !first)
      {
         parse (content, time_stamp);
         content.clear ();
      }

      content.push_back (input_string);
      first = false;

   }

   parse (content, time_stamp);
   file.close ();

}

set<Dtime>
Advisory_Store::get_initial_time_set (const Area* area_ptr,
                                      const bool synoptic_hours_only,
                                      const bool sub_synoptic_okay) const
{

   set<Dtime> initial_time_set;

   for (Advisory_Store::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Dstring& key = iterator->first;
      const Tokens tokens (key);

      if (area_ptr != NULL)
      {
         const Real latitude = stof (tokens[2]);
         const Real longitude = stof (tokens[3]);
         const Lat_Long lat_long (latitude, longitude);
         if (!area_ptr->contains (lat_long)) { continue; }
      }

      const Dstring& initial_time_string = tokens[1];
      const Dtime& initial_time (initial_time_string);

      if (synoptic_hours_only)
      {
         const Integer hours = initial_time.get_hour ();
         if (hours % 3 != 0) { continue; }
         if (!sub_synoptic_okay && (hours % 6 != 0)) { continue; }
      }

      initial_time_set.insert (initial_time);

   }

   return initial_time_set;

}

Tokens
Advisory_Store::get_key_tokens (const Area* area_ptr) const
{
   set<Dstring> initial_time_string_set;
   return get_key_tokens (initial_time_string_set, area_ptr);
}

Tokens
Advisory_Store::get_key_tokens (const Dstring& initial_time_string,
                                const Area* area_ptr) const
{

   set<Dstring> initial_time_string_set;
   initial_time_string_set.insert (initial_time_string);

   return get_key_tokens (initial_time_string_set, area_ptr);

}

Tokens
Advisory_Store::get_key_tokens (const set<Dstring>& initial_time_string_set,
                                const Area* area_ptr) const
{

   Tokens key_tokens;

   for (Advisory_Store::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Dstring& key = iterator->first;
      const Tokens tokens (key);
      const Dstring& initial_time_string = tokens[1];

      if (initial_time_string_set.size () > 0)
      {
         bool match_its = false;
         typedef set<Dstring>::const_iterator Iterator;
         for (Iterator i = initial_time_string_set.begin ();
              (!match_its && i != initial_time_string_set.end ()); i++)
         {
            const Dstring& it = *(i);
            if (initial_time_string == it) { match_its = true; }
         }
         if (!match_its) { continue; }
      }

      if (area_ptr != NULL)
      {
         const Real latitude = stof (tokens[2]);
         const Real longitude = stof (tokens[3]);
         const Lat_Long lat_long (latitude, longitude);
         if (!area_ptr->contains (lat_long)) { continue; }
      }

      key_tokens.push_back (key);

   }

   return key_tokens;

}

Tokens
Advisory_Store::get_key_tokens (const Dtime& initial_time,
                                const Area* area_ptr) const
{

   set<Dtime> initial_time_set;
   initial_time_set.insert (initial_time);

   return get_key_tokens (initial_time_set, area_ptr);

}

Tokens
Advisory_Store::get_key_tokens (const set<Dtime>& initial_time_set,
                                const Area* area_ptr) const
{

   Tokens key_tokens;

   for (Advisory_Store::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Dstring& key = iterator->first;
      const Tokens tokens (key);
      const Dstring& initial_time_str = tokens[1];
      const Dtime initial_time (initial_time_str);

      if (initial_time_set.size () > 0)
      {
         bool match_initial_time = false;
         for (set<Dtime>::const_iterator i = initial_time_set.begin ();
              (!match_initial_time && i != initial_time_set.end ()); i++)
         {
            const Dtime& it = *(i);
            const Dstring& its = it.get_string ();
            if (initial_time_str == its) { match_initial_time = true; }
         }
         if (!match_initial_time) { continue; }
      }

      if (area_ptr != NULL)
      {
         const Real latitude = stof (tokens[2]);
         const Real longitude = stof (tokens[3]);
         const Lat_Long lat_long (latitude, longitude);
         if (!area_ptr->contains (lat_long)) { continue; }
      }

      key_tokens.push_back (key);

   }

   return key_tokens;

}

Advisory_Store::Cluster_Info
Advisory_Store::get_cluster_info (const Tokens& key_tokens,
                                  const Real criteria_distance) const
{

   const Integer n = key_tokens.size ();
   Dataset dataset (2, n);

   for (Integer i = 0; i < n; i++)
   {

      const Dstring& key = key_tokens[i];
      const Tokens tokens (key);

      const Real latitude = stof (tokens[2]);
      const Real longitude = stof (tokens[3]);
      dataset.set_datum (i, 0, latitude);
      dataset.set_datum (i, 1, longitude);

   }

   return Advisory_Store::Cluster_Info (*this,
      key_tokens, dataset, criteria_distance);

}

Cluster::Multimap
Advisory_Store::get_cluster_multimap (const Tokens& key_tokens,
                                      const Real criteria_distance) const
{

   const Integer n = key_tokens.size ();
   Dataset dataset (2, n);

   for (Integer i = 0; i < n; i++)
   {

      const Dstring& key = key_tokens[i];
      const Tokens tokens (key);

      const Real latitude = stof (tokens[2]);
      const Real longitude = stof (tokens[3]);
      dataset.set_datum (i, 0, latitude);
      dataset.set_datum (i, 1, longitude);

   }

   return Cluster (dataset).get_multimap (criteria_distance);

}

vector<Point_2D>
Advisory_Store::get_position_vector (const Tokens& key_tokens,
                                     const Real tau) const
{

   vector<Point_2D> position_vector;

   for (Tokens::const_iterator iterator = key_tokens.begin ();
        iterator != key_tokens.end (); iterator++)
   {
      const Dstring& key = *(iterator);
      const Advisory& advisory = at (key);
      const Lat_Long ll = advisory.get_lat_long (tau);
      if (!ll.is_nall ()) { position_vector.push_back (ll); }
   }

   return position_vector;

}

Bivariate_Gaussian_Distribution
Advisory_Store::get_bivariate_gaussian_distribution (const Tokens& key_tokens,
                                                     const Real tau) const
{
   const vector<Point_2D>& position_vector =
      get_position_vector (key_tokens, tau);
   return Bivariate_Gaussian_Distribution (position_vector);
}

Ellipse
Advisory_Store::get_ellipse (const Tokens& key_tokens,
                             const Real tau,
                             const Real probability) const
{
   const Bivariate_Gaussian_Distribution& bgd =
      get_bivariate_gaussian_distribution (key_tokens, tau);
   return bgd.get_ellipse (probability);
}

