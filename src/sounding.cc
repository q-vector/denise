//
// sounding.cc
// 
// Copyright (C) 2005-2015 Simon E. Ching
// 
// This file is part of the denise suite.
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

#include "andrea.h"
#include "sounding.h"

using namespace std;
using namespace denise;
using namespace andrea;

Sounding_Package::Sounding_Package (Andrea& andrea)
   : Package (andrea)
{
}

void
Sounding_Package::sounding_load (const Dstring& identifier,
                                 const Dstring& file_path)
{
   Sounding sounding (file_path);
   sounding_map[identifier] = sounding;
}

void
Sounding_Package::sounding_print (const Dstring& identifier,
                                  const Tokens& tokens) const
{

   const Dstring time_fmt (L"%Y%m%d%H%M");
   const Sounding& sounding = sounding_map.at (identifier);
   if (tokens.size () < 1) { return; }

   const Dstring& genre = tokens[0];
   const Tokens& arguments = tokens.subtokens (1);

   if (genre == L"time")
   {
      wcout << sounding.get_time ().get_string (time_fmt) << endl;
      return;
   }
   else
   if (genre == L"basetime")
   {
      wcout << sounding.get_basetime ().get_string (time_fmt) << endl;
      return;
   }
   else
   if (genre == L"location")
   {
      wcout << sounding.get_location_str () << endl;
      return;
   }
   else
   if (genre == L"t")
   {

      if (arguments.size () == 0)
      {
         const Thermo_Line& t_line = sounding.get_t_line ();
         for (auto iterator = t_line.begin ();
              iterator != t_line.end (); iterator++)
         {
            const Real p = iterator->first;
            const Real t = iterator->second;
            wcout << p << L" " << t << endl;
         }
      }
      else
      if (arguments.size () == 1)
      {
         const Real p = stof (arguments[0]);
         const Real t = sounding.get_temperature (tephigram, p);
         wcout << t << endl;
      }

   }
   else
   if (genre == L"td")
   {

      if (arguments.size () == 0)
      {
         const Thermo_Line& t_d_line = sounding.get_t_d_line ();
         for (auto iterator = t_d_line.begin ();
              iterator != t_d_line.end (); iterator++)
         {
            const Real p = iterator->first;
            const Real t_d = iterator->second;
            wcout << p << L" " << t_d << endl;
         }
      }
      else
      if (arguments.size () == 1)
      {
         const Real p = stof (arguments[0]);
         const Real t_d = sounding.get_dew_point (tephigram, p);
         wcout << t_d << endl;
      }

   }
   else
   if (genre == L"wind")
   {

      if (arguments.size () == 0)
      {
         const Wind_Profile& wind_profile = sounding.get_wind_profile ();
         for (auto iterator = wind_profile.begin ();
              iterator != wind_profile.end (); iterator++)
         {
            const Real p = iterator->first;
            const Wind& wind = iterator->second;
            const Real wind_direction = wind.get_direction ();
            const Real wind_speed = wind.get_speed ();
            wcout << p << L" " << wind_direction << L" " << wind_speed << endl;
         }
      }
      else
      if (arguments.size () == 1)
      {
         const Real p = stof (arguments[0]);
         const Wind& wind = sounding.get_wind (p);
         const Real wind_direction = wind.get_direction ();
         const Real wind_speed = wind.get_speed ();
         wcout << wind_direction << L" " << wind_speed << endl;
      }

   }
   else
   if (genre == L"height")
   {

      if (arguments.size () == 0)
      {
         const Height_Profile& height_profile = sounding.get_height_profile ();
         for (auto iterator = height_profile.begin ();
              iterator != height_profile.end (); iterator++)
         {
            const Real p = iterator->first;
            const Real height = iterator->second;
            wcout << p << L" " << height << endl;
         }
      }
      else
      if (arguments.size () == 1)
      {
         const Real p = stof (arguments[0]);
         const Real height = sounding.get_height (p);
         wcout << height << endl;
      }

   }
   else
   if (genre == L"brunt_vaisala")
   {

      const Real_Profile* brunt_vaisala_profile_ptr =
         sounding.get_brunt_vaisala_profile_ptr ();
      for (auto iterator = brunt_vaisala_profile_ptr->begin ();
           iterator != brunt_vaisala_profile_ptr->end (); iterator++)
      {
         const Real p = iterator->first;
         const Real brunt_vaisala = iterator->second;
         wcout << p << L" " << brunt_vaisala << endl;
      }
      delete brunt_vaisala_profile_ptr;

   }

}

void
Sounding_Package::sounding_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == L"load")
   {
      const Dstring& identifier = tokens[1];
      const Dstring& file_path = tokens[2];
      sounding_load (identifier, file_path);
   }
   else
   if (tokens[0] == L"print")
   {
      const Dstring& identifier = tokens[1];
      sounding_print (identifier, tokens.subtokens (2));
   }

}

const map<Dstring, Sounding>&
Sounding_Package::get_sounding_map () const
{
   return sounding_map;
}

const Sounding&
Sounding_Package::get_sounding (const Dstring& identifier) const
{
   Exception e (L"sounding not found: " + identifier);
   try { return sounding_map.at (identifier); }
   catch (const std::out_of_range& oor) { throw e; }
}

void
Sounding_Package::surface_sounding (const Tokens& tokens) const
{

   const Dstring& operation = tokens[0];

   if (tokens[0] == L"tephigram")
   {
      surface_sounding_tephigram (tokens.subtokens (1));
   }
   else
   if (tokens[0] == L"chart")
   {
      surface_sounding_chart (tokens.subtokens (1));
   }

}

void
Sounding_Package::surface_sounding_tephigram (const Tokens& tokens) const
{

   const Dstring& surface_identifier = tokens[0];
   const Dstring& sounding_identifier = tokens[1];

   const RefPtr<Surface>& surface = andrea.get_surface (surface_identifier);
   const RefPtr<Context> cr = andrea.get_cr (surface_identifier);
   const Size_2D& size_2d = andrea.get_size_2d (surface_identifier);
   const Sounding& sounding = andrea.get_sounding (sounding_identifier);
   const Tephigram tephigram (size_2d);

   Color::white ().cairo (cr);
   cr->paint ();

   cr->set_line_width (1);
   tephigram.render (cr);

   cr->set_line_width (2);
   sounding.render (cr, tephigram);

}

void
Sounding_Package::surface_sounding_chart (const Tokens& tokens) const
{

   const Dstring& surface_identifier = tokens[0];
   const Dstring& sounding_identifier = tokens[1];
   const Dstring& x_str = tokens[2];
   const Dstring& y_str = tokens[3];
   const Dstring& genre = tokens[4];

   const RefPtr<Surface>& surface = andrea.get_surface (surface_identifier);
   const RefPtr<Context> cr = andrea.get_cr (surface_identifier);
   const Sounding& sounding = andrea.get_sounding (sounding_identifier);

   const Tokens x_tokens (x_str, L"/");
   const Tokens y_tokens (y_str, L"/");
   const bool is_p = (y_tokens[0][0] == 'p');
   const Domain_1D domain_x (stof (x_tokens[0]), stof (x_tokens[1]));
   const Domain_1D domain_y (stof (y_tokens[1]), stof (y_tokens[2]));
   const Domain_2D domain_2d (domain_x, domain_y);

   const Real title_height = 40;
   const Real margin_l = 75;
   const Real margin_r = 50;
   const Real margin_t = title_height + 50;
   const Real margin_b = 50;

   const Size_2D& size_2d = andrea.get_size_2d (surface_identifier);
   const Real w = size_2d.i - margin_l - margin_r;
   const Real h = size_2d.j - margin_t - margin_b;

   const Color minor_color = Color::black (0.1);
   const Color major_color = Color::black (0.3);

   const Dstring fmt_y (L"%.0f");
   const Dstring unit_y (is_p ? L"Pa" : L"m");
   const Real minor_interval_y = is_p ? 10e2 : 100;
   const Real major_interval_y = is_p ? 100e2 : 1000;

   Mesh_2D mesh_2d (size_2d, domain_2d);
   Dstring fmt_x, unit_x;
   const Real_Profile* real_profile_ptr;

   if (genre == L"height")
   {
      const Real minor_interval_x = 100;
      const Real major_interval_x = 1000;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color);
      fmt_x = Dstring (L"%.0f");
      unit_x = Dstring (L"\u00b0C");
      real_profile_ptr = sounding.get_height_profile_ptr ();
   }
   else
   if (genre == L"theta")
   {
      const Real minor_interval_x = 1;
      const Real major_interval_x = 10;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color);
      fmt_x = Dstring (L"%.0f");
      unit_x = Dstring (L"\u00b0C");
      real_profile_ptr = sounding.get_theta_profile_ptr ();
   }
   else
   if (genre == L"speed")
   {
      const Real minor_interval_x = 1;
      const Real major_interval_x = 10;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color);
      fmt_x = Dstring (L"%.0f");
      unit_x = Dstring (L"ms\u207b\u00b9");
      real_profile_ptr = sounding.get_speed_profile_ptr ();
   }
   else
   if (genre == L"along_speed")
   {
      const Real azimuth = stof (tokens[5]);
      const Real minor_interval_x = 1;
      const Real major_interval_x = 10;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color);
      fmt_x = Dstring (L"%.0f");
      unit_x = Dstring (L"ms\u207b\u00b9");
      real_profile_ptr = sounding.get_along_speed_profile_ptr (azimuth);
   }
   else
   if (genre == L"brunt_vaisala")
   {
      const Real minor_interval_x = 0.001;
      const Real major_interval_x = 0.01;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color);
      fmt_x = Dstring (L"%.2f");
      unit_x = Dstring (L"s\u207b\u00b9");
      real_profile_ptr = sounding.get_brunt_vaisala_profile_ptr ();
   }
   else
   if (genre == L"scorer")
   {
      const Real azimuth = stof (tokens[5]);
      const Real minor_interval_x = 1e-6;
      const Real major_interval_x = 1e-5;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color,
         1e23, 1e23, Color::black (1.0));
      fmt_x = Dstring (L"%g");
      unit_x = Dstring (L"m\u207b\u00b2");
      real_profile_ptr = sounding.get_scorer_profile_ptr (azimuth);
   }

   Affine_Transform_2D transform;
   const Real span_x = domain_x.get_span ();
   const Real span_y = domain_y.get_span ();
   transform.scale (1, -1);
   transform.translate (-domain_x.start, 0);
   transform.scale (w / span_x, h / span_y);
   transform.translate (margin_l, size_2d.j - margin_b);

   surface_sounding_chart (cr, transform, is_p, mesh_2d, fmt_x, fmt_y, unit_x,
      unit_y, sounding, *real_profile_ptr, Ring (4), Color::red (0.4));

   delete real_profile_ptr;

}

void
Sounding_Package::surface_sounding_chart (const RefPtr<Context>& cr,
                                          const Transform_2D& transform,
                                          const bool is_p,
                                          const Mesh_2D& mesh_2d,
                                          const Dstring& fmt_x,
                                          const Dstring& fmt_y,
                                          const Dstring& unit_x,
                                          const Dstring& unit_y,
                                          const Sounding& sounding,
                                          const Real_Profile& real_profile,
                                          const Symbol& symbol,
                                          const Color& color) const
{

   const Domain_2D& domain_2d = mesh_2d.get_domain_2d ();
   const Real start_x = domain_2d.domain_x.start;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_x = domain_2d.domain_x.end;
   const Real end_y = domain_2d.domain_y.end;

   Color::white ().cairo (cr);
   cr->paint ();

   cr->set_line_width (2);
   Color::black ().cairo (cr);
   mesh_2d.render (cr, transform);
   mesh_2d.render_label_x (cr, transform, 0, start_y,
      fmt_x, Label::REAL, 'c', 't', 9);
   mesh_2d.render_label_y (cr, transform, 0, start_x,
      fmt_y, Label::REAL, 'r', 'c', 9);

   Label (unit_x, Point_2D (end_x, start_y), 'l', 'c', 9).cairo (cr, transform);
   Label (unit_y, Point_2D (start_x, end_y), 'r', 'b', 9).cairo (cr, transform);

   color.cairo (cr);

   for (auto iterator = real_profile.begin ();
        iterator != real_profile.end (); iterator++)
   {
      const Real p = iterator->first;
      const Real y = is_p ? p : sounding.get_height (p);
      const Real datum = iterator->second;
      const Point_2D point = transform.transform (Point_2D (datum, y));
      symbol.cairo (cr, point);
      cr->fill ();
   }

}

