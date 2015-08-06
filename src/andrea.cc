#include <denise/dstring.h>
#include <denise/met.h>
#include "andrea.h"

using namespace std;
using namespace denise;

Andrea_Package::Andrea_Package (const Andrea& andrea)
   : andrea (andrea)
{
}

void
Geodesy_Package::geodesy_assign (const string& identifier,
                                 const string& str)
{
   geodesy_map[identifier] = Geodesy (str);
}

void
Geodesy_Package::geodesy_print (const string& identifier) const
{

   auto iterator = geodesy_map.find (identifier);
   const bool is_present = (iterator != geodesy_map.end ());
   if (is_present)
   {
      cout << "geodesy " << identifier << " is present" << endl;
   }

}

void
Geodesy_Package::geodesy_distance (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 2) { throw Exception ("geodesy_distance"); }

   const Lat_Long origin (tokens[0]);
   const Lat_Long destination (tokens[1]);
   const Real distance = Geodesy::get_distance (origin, destination);

   cout << distance << endl;

}

void
Geodesy_Package::geodesy_azimuth (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 2) { throw Exception ("geodesy_azimuth"); }

   const Lat_Long origin (tokens[0]);
   const Lat_Long destination (tokens[1]);
   const Real azimuth_f = Geodesy::get_azimuth_forward (origin, destination);
   const Real azimuth_b = Geodesy::get_azimuth_backward (origin, destination);

   cout << azimuth_f << " " << azimuth_b << endl;

}

void
Geodesy_Package::geodesy_destination (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 3) { throw Exception ("geodesy_destination"); }

   const Lat_Long origin (tokens[0]);
   const Real distance = stof (tokens[1]);
   const Real azimuth_forward = stof (tokens[2]);
   const Lat_Long& destination =
      Geodesy::get_destination (origin, distance, azimuth_forward);

   cout << destination.get_string (lat_long_dp) << endl;

}

Geodesy_Package::Geodesy_Package (const Andrea& andrea)
   : Andrea_Package (andrea),
     lat_long_dp (4)
{
}

void
Geodesy_Package::geodesy_parse (const Tokens& tokens)
{


   const Integer n = tokens.size ();

   if (tokens[0] == "assign")
   {
      const string& identifier = tokens[1];
      geodesy_assign (identifier, tokens[2]);
   }
   else
   if (tokens[0] == "print")
   {
      const string& identifier = tokens[1];
      geodesy_print (identifier);
   }
   else
   if (tokens[0] == "distance")
   {
      geodesy_distance (tokens.subtokens (1));
   }
   else
   if (tokens[0] == "azimuth")
   {
      geodesy_azimuth (tokens.subtokens (1));
   }
   else
   if (tokens[0] == "destination")
   {
      geodesy_destination (tokens.subtokens (1));
   }
   else
   {
      throw Exception ("geodesy_parse");
   }

}

const map<string, Geodesy>&
Geodesy_Package::get_geodesy_map () const
{
   return geodesy_map;
}

const Geodesy&
Geodesy_Package::get_geodesy (const string& identifier) const
{
   return geodesy_map.at (identifier);
}

Journey_Package::Journey_Package (const Andrea& andrea)
   : Andrea_Package (andrea)
{
}

void
Journey_Package::journey_assign (const string& identifier,
                                 const Tokens& arguments)
{

   Journey journey;

   for (auto iterator = arguments.begin ();
        iterator != arguments.end (); iterator++)
   {
      const Lat_Long lat_long (*(iterator));
      journey.push_back (lat_long);
   }

   journey_map[identifier] = journey;

}

void
Journey_Package::journey_print (const string& identifier) const
{
   auto iterator = journey_map.find (identifier);
   const bool is_present = (iterator != journey_map.end ());
   if (is_present)
   {
      cout << "journey " << identifier << " is present" << endl;
   }
}

void
Journey_Package::journey_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "assign")
   {
      const string& identifier = tokens[1];
      journey_assign (identifier, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "print")
   {
      const string& identifier = tokens[1];
      journey_print (identifier);
   }

}

const map<string, Journey>&
Journey_Package::get_journey_map () const
{
   return journey_map;
}

const Journey&
Journey_Package::get_journey (const string& identifier) const
{
   return journey_map.at (identifier);
}


Geodetic_Mesh_Package::Geodetic_Mesh_Package (const Andrea& andrea)
   : Andrea_Package (andrea)
{
}

void
Geodetic_Mesh_Package::geodetic_mesh_assign (const string& identifier,
                                             const Size_2D& size_2d,
                                             const Domain_2D& domain_2d)
{
   const Geodetic_Mesh geodetic_mesh (size_2d, domain_2d);
   geodetic_mesh_map[identifier] = geodetic_mesh;
}

void
Geodetic_Mesh_Package::geodetic_mesh_add (const string& identifier,
                                          const Tokens& arguments)
{

   const Integer n = arguments.size ();
   Geodetic_Mesh& geodetic_mesh = geodetic_mesh_map.at (identifier);

   switch (n)
   {
      case 2:
      {
         const Real interval = stof (arguments[0]); 
         const Color color (arguments[1]);
         geodetic_mesh.add (Simple_Mesh_2D (interval, color));
         break;
      }

      case 3:
      {
         const Real interval_lat = stof (arguments[0]); 
         const Real interval_long = stof (arguments[1]); 
         const Color color (arguments[2]);
         const Simple_Mesh_2D sm (interval_lat, interval_long, color);
         geodetic_mesh.add (sm);
         break;
      }

   } 

}

void
Geodetic_Mesh_Package::geodetic_mesh_print (const string& identifier,
                                            const Tokens& arguments) const
{
   auto iterator = geodetic_mesh_map.find (identifier);
   const bool is_present = (iterator != geodetic_mesh_map.end ());
   if (is_present)
   {
      cout << "geodetic_mesh " << identifier << " is present" << endl;
   }
}

void
Geodetic_Mesh_Package::geodetic_mesh_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "assign")
   {
      const string& identifier = tokens[1];
      const Size_2D size_2d (tokens[2]);
      const Domain_2D domain_2d (tokens[3]);
      geodetic_mesh_assign (identifier, size_2d, domain_2d);
   }
   else
   if (tokens[0] == "add")
   {
      const string& identifier = tokens[1];
      geodetic_mesh_add (identifier, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "print")
   {
      const string& identifier = tokens[1];
      geodetic_mesh_print (identifier, tokens.subtokens (2));
   }

}

const map<string, Geodetic_Mesh>&
Geodetic_Mesh_Package::get_geodetic_mesh_map () const
{
   return geodetic_mesh_map;
}

const Geodetic_Mesh&
Geodetic_Mesh_Package::get_geodetic_mesh (const string& identifier) const
{
   return geodetic_mesh_map.at (identifier);
}

Geodetic_Transform_Package::Geodetic_Transform_Package (const Andrea& andrea)
   : Andrea_Package (andrea)
{
}

void
Geodetic_Transform_Package::geodetic_transform_assign (const string& identifier,
                                                       const string& str)
{
   geodetic_transform_str_map[identifier] = str;
}

void
Geodetic_Transform_Package::geodetic_transform_print (const string& identifier,
                                                      const Tokens& arguments) const
{
   auto iterator = geodetic_transform_str_map.find (identifier);
   const bool is_present = (iterator != geodetic_transform_str_map.end ());
   if (is_present)
   {
      cout << "geodetic_transform " << identifier << " is present" << endl;
   }
}

void
Geodetic_Transform_Package::geodetic_transform_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "assign")
   {
      const string& identifier = tokens[1];
      const string& str = tokens[2];
      geodetic_transform_assign (identifier, str);
   }
   else
   if (tokens[0] == "print")
   {
      const string& identifier = tokens[1];
      geodetic_transform_print (identifier, tokens.subtokens (2));
   }

}

const map<string, string>&
Geodetic_Transform_Package::get_geodetic_transform_str_map () const
{
   return geodetic_transform_str_map;
}

const Geodetic_Transform*
Geodetic_Transform_Package::get_geodetic_transform_ptr (const string& identifier,
                                                        const Point_2D& point) const
{
   typedef Geodetic_Transform Gt;
   const string& str = geodetic_transform_str_map.at (identifier);
   return Gt::get_transform_ptr (str, point);
}

Gshhs_Package::Gshhs_Package (const Andrea& andrea)
   : Andrea_Package (andrea)
{
}

Gshhs_Package::~Gshhs_Package ()
{
   for (auto iterator = gshhs_ptr_map.begin ();
        iterator != gshhs_ptr_map.end (); iterator++)
   {
      Gshhs* gshhs_ptr = iterator->second;
      delete gshhs_ptr;
   }
}

void
Gshhs_Package::gshhs_load (const string& identifier,
                           const string& file_path)
{

   auto iterator = gshhs_ptr_map.find (identifier);
   const bool is_present = (iterator != gshhs_ptr_map.end ());
   if (is_present) { delete iterator->second; }

   Gshhs* gshhs_ptr = new Gshhs (file_path);
   gshhs_ptr_map[identifier] = gshhs_ptr;

}

void
Gshhs_Package::gshhs_print (const string& identifier,
                            const Tokens& arguments) const
{
   auto iterator = gshhs_ptr_map.find (identifier);
   const bool is_present = (iterator != gshhs_ptr_map.end ());
   if (is_present)
   {
      cout << "gshhs " << identifier << " is present" << endl;
   }
}

void
Gshhs_Package::gshhs_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "load")
   {
      const string& identifier = tokens[1];
      const string& file_path = tokens[2];
      gshhs_load (identifier, file_path);
   }
   else
   if (tokens[0] == "print")
   {
      const string& identifier = tokens[1];
      gshhs_print (identifier, tokens.subtokens (2));
   }

}

const map<string, Gshhs*>&
Gshhs_Package::get_gshhs_ptr_map () const
{
   return gshhs_ptr_map;
}

const Gshhs*
Gshhs_Package::get_gshhs_ptr (const string& identifier) const
{
   return gshhs_ptr_map.at (identifier);
}

Sounding_Package::Sounding_Package (const Andrea& andrea)
   : Andrea_Package (andrea)
{
}

void
Sounding_Package::sounding_load (const string& identifier,
                                 const string& file_path)
{
   Sounding sounding (file_path);
   sounding_map[identifier] = sounding;
}

void
Sounding_Package::sounding_print (const string& identifier,
                                  const Tokens& tokens) const
{

   const string time_fmt ("%Y%m%d%H%M");
   const Sounding& sounding = sounding_map.at (identifier);
   if (tokens.size () < 1) { return; }

   const string& genre = tokens[0];
   const Tokens& arguments = tokens.subtokens (1);

   if (genre == "time")
   {
      cout << sounding.get_time ().get_string (time_fmt) << endl;
      return;
   }
   else
   if (genre == "basetime")
   {
      cout << sounding.get_basetime ().get_string (time_fmt) << endl;
      return;
   }
   else
   if (genre == "location")
   {
      cout << sounding.get_location_str () << endl;
      return;
   }
   else
   if (genre == "t")
   {

      if (arguments.size () == 0)
      {
         const Thermo_Line& t_line = sounding.get_t_line ();
         for (auto iterator = t_line.begin ();
              iterator != t_line.end (); iterator++)
         {
            const Real p = iterator->first;
            const Real t = iterator->second;
            cout << p << " " << t << endl;
         }
      }
      else
      if (arguments.size () == 1)
      {
         const Real p = stof (arguments[0]);
         const Real t = sounding.get_temperature (tephigram, p);
         cout << t << endl;
      }

   }
   else
   if (genre == "td")
   {

      if (arguments.size () == 0)
      {
         const Thermo_Line& t_d_line = sounding.get_t_d_line ();
         for (auto iterator = t_d_line.begin ();
              iterator != t_d_line.end (); iterator++)
         {
            const Real p = iterator->first;
            const Real t_d = iterator->second;
            cout << p << " " << t_d << endl;
         }
      }
      else
      if (arguments.size () == 1)
      {
         const Real p = stof (arguments[0]);
         const Real t_d = sounding.get_dew_point (tephigram, p);
         cout << t_d << endl;
      }

   }
   else
   if (genre == "wind")
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
            cout << p << " " << wind_direction << " " << wind_speed << endl;
         }
      }
      else
      if (arguments.size () == 1)
      {
         const Real p = stof (arguments[0]);
         const Wind& wind = sounding.get_wind (p);
         const Real wind_direction = wind.get_direction ();
         const Real wind_speed = wind.get_speed ();
         cout << wind_direction << " " << wind_speed << endl;
      }

   }
   else
   if (genre == "height")
   {

      if (arguments.size () == 0)
      {
         const Height_Profile& height_profile = sounding.get_height_profile ();
         for (auto iterator = height_profile.begin ();
              iterator != height_profile.end (); iterator++)
         {
            const Real p = iterator->first;
            const Real height = iterator->second;
            cout << p << " " << height << endl;
         }
      }
      else
      if (arguments.size () == 1)
      {
         const Real p = stof (arguments[0]);
         const Real height = sounding.get_height (p);
         cout << height << endl;
      }

   }
   else
   if (genre == "brunt_vaisala")
   {

      const Real_Profile* brunt_vaisala_profile_ptr =
         sounding.get_brunt_vaisala_profile_ptr ();
      for (auto iterator = brunt_vaisala_profile_ptr->begin ();
           iterator != brunt_vaisala_profile_ptr->end (); iterator++)
      {
         const Real p = iterator->first;
         const Real brunt_vaisala = iterator->second;
         cout << p << " " << brunt_vaisala << endl;
      }
      delete brunt_vaisala_profile_ptr;

   }

}

void
Sounding_Package::sounding_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "load")
   {
      const string& identifier = tokens[1];
      const string& file_path = tokens[2];
      sounding_load (identifier, file_path);
   }
   else
   if (tokens[0] == "print")
   {
      const string& identifier = tokens[1];
      sounding_print (identifier, tokens.subtokens (2));
   }

}

const map<string, Sounding>&
Sounding_Package::get_sounding_map () const
{
   return sounding_map;
}

const Sounding&
Sounding_Package::get_sounding (const string& identifier) const
{
   return sounding_map.at (identifier);
}

Image_Package::Image_Package (const Andrea& andrea)
   : Andrea_Package (andrea)
{
}

void
Image_Package::image_init (const string& identifier,
                           const string& geometry)
{

   const Tokens tokens (geometry, "x");
   const Size_2D size_2d (stof (tokens[0]), stof (tokens[1]));

   const RefPtr<ImageSurface> image_surface = get_image_surface (size_2d);
   const RefPtr<Context> cr = denise::get_cr (image_surface);

   image_map[identifier] = image_surface;
   cr_map[identifier] = cr;

}

void
Image_Package::image_paint (const string& identifier,
                            const Color& color) const
{
   const RefPtr<Context> cr = cr_map.at (identifier);
   color.cairo (cr);
   cr->paint ();
}

void
Image_Package::image_save (const string& identifier,
                           const string& file_path) const
{
   const RefPtr<ImageSurface>& image_surface = image_map.at (identifier);
   image_surface->write_to_png (file_path);
}

void
Image_Package::image_title (const string& identifier,
                            const Tokens& tokens) const
{

   const RefPtr<ImageSurface> image_surface = image_map.at (identifier);
   const RefPtr<Context> cr = cr_map.at (identifier);
   const Integer w = image_surface->get_width ();
   const Integer h = image_surface->get_height ();
   const Size_2D size_2d (w, h);

   Title title (size_2d);
   title.set (tokens);
   title.cairo (cr);

}

void
Image_Package::image_sounding (const Tokens& tokens) const
{

   const string& operation = tokens[0];

   if (tokens[0] == "tephigram")
   {
      image_sounding_tephigram (tokens.subtokens (1));
   }
   else
   if (tokens[0] == "chart")
   {
      image_sounding_chart (tokens.subtokens (1));
   }

}

void
Image_Package::image_sounding_tephigram (const Tokens& tokens) const
{

   const string& image_identifier = tokens[0];
   const string& sounding_identifier = tokens[1];

   const RefPtr<ImageSurface>& image_surface = image_map.at (image_identifier);
   const RefPtr<Context> cr = cr_map.at (image_identifier);
   const Sounding& sounding =
      andrea.get_sounding_map ().at (sounding_identifier);

   const Integer w = image_surface->get_width ();
   const Integer h = image_surface->get_height ();
   const Size_2D size_2d (w, h);
   const Tephigram tephigram (size_2d);

   Color::white ().cairo (cr);
   cr->paint ();

   cr->set_line_width (1);
   tephigram.render (cr);

   cr->set_line_width (2);
   sounding.render (cr, tephigram);

}

void
Image_Package::image_sounding_chart (const Tokens& tokens) const
{

   const string& image_identifier = tokens[0];
   const string& sounding_identifier = tokens[1];
   const string& x_str = tokens[2];
   const string& y_str = tokens[3];
   const string& genre = tokens[4];

   const RefPtr<ImageSurface>& image_surface = image_map.at (image_identifier);
   const RefPtr<Context> cr = cr_map.at (image_identifier);
   const Sounding& sounding =
      andrea.get_sounding_map ().at (sounding_identifier);

   const Tokens x_tokens (x_str, "/");
   const Tokens y_tokens (y_str, "/");
   const bool is_p = (y_tokens[0][0] == 'p');
   const Domain_1D domain_x (stof (x_tokens[0]), stof (x_tokens[1]));
   const Domain_1D domain_y (stof (y_tokens[1]), stof (y_tokens[2]));
   const Domain_2D domain_2d (domain_x, domain_y);

   const Integer image_width = image_surface->get_width ();
   const Integer image_height = image_surface->get_height ();
   const Size_2D size_2d (image_width, image_height);

   const Real title_height = 40;
   const Real margin_l = 75;
   const Real margin_r = 50;
   const Real margin_t = title_height + 50;
   const Real margin_b = 50;

   const Real w = size_2d.i - margin_l - margin_r;
   const Real h = size_2d.j - margin_t - margin_b;

   const Color minor_color = Color::black (0.1);
   const Color major_color = Color::black (0.3);

   const string fmt_y ("%.0f");
   const string unit_y (is_p ? "Pa" : "m");
   const Real minor_interval_y = is_p ? 10e2 : 100;
   const Real major_interval_y = is_p ? 100e2 : 1000;

   Mesh_2D mesh_2d (size_2d, domain_2d);
   string fmt_x, unit_x;
   const Real_Profile* real_profile_ptr;

   if (genre == "height")
   {
      const Real minor_interval_x = 100;
      const Real major_interval_x = 1000;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color);
      fmt_x = string ("%.0f");
      unit_x = string ("\u00b0C");
      real_profile_ptr = sounding.get_height_profile_ptr ();
   }
   else
   if (genre == "theta")
   {
      const Real minor_interval_x = 1;
      const Real major_interval_x = 10;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color);
      fmt_x = string ("%.0f");
      unit_x = string ("\u00b0C");
      real_profile_ptr = sounding.get_theta_profile_ptr ();
   }
   else
   if (genre == "speed")
   {
      const Real minor_interval_x = 1;
      const Real major_interval_x = 10;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color);
      fmt_x = string ("%.0f");
      unit_x = string ("ms\u207b\u00b9");
      real_profile_ptr = sounding.get_speed_profile_ptr ();
   }
   else
   if (genre == "along_speed")
   {
      const Real azimuth = stof (tokens[5]);
      const Real minor_interval_x = 1;
      const Real major_interval_x = 10;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color);
      fmt_x = string ("%.0f");
      unit_x = string ("ms\u207b\u00b9");
      real_profile_ptr = sounding.get_along_speed_profile_ptr (azimuth);
   }
   else
   if (genre == "brunt_vaisala")
   {
      const Real minor_interval_x = 0.001;
      const Real major_interval_x = 0.01;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color);
      fmt_x = string ("%.2f");
      unit_x = string ("s\u207b\u00b9");
      real_profile_ptr = sounding.get_brunt_vaisala_profile_ptr ();
   }
   else
   if (genre == "scorer")
   {
      const Real azimuth = stof (tokens[5]);
      const Real minor_interval_x = 1e-6;
      const Real major_interval_x = 1e-5;
      mesh_2d = Mesh_2D (Size_2D (2, 2), domain_2d,
         major_interval_x, major_interval_y, major_color,
         minor_interval_x, minor_interval_y, minor_color,
         1e23, 1e23, Color::black (1.0));
      fmt_x = string ("%g");
      unit_x = string ("m\u207b\u00b2");
      real_profile_ptr = sounding.get_scorer_profile_ptr (azimuth);
   }

   Affine_Transform_2D transform;
   const Real span_x = domain_x.get_span ();
   const Real span_y = domain_y.get_span ();
   transform.scale (1, -1);
   transform.translate (-domain_x.start, 0);
   transform.scale (w / span_x, h / span_y);
   transform.translate (margin_l, size_2d.j - margin_b);

   image_sounding_chart (cr, transform, is_p, mesh_2d, fmt_x, fmt_y, unit_x,
      unit_y, sounding, *real_profile_ptr, Ring (4), Color::red (0.4));

   delete real_profile_ptr;

}

void
Image_Package::image_sounding_chart (const RefPtr<Context>& cr,
                                     const Transform_2D& transform,
                                     const bool is_p,
                                     const Mesh_2D& mesh_2d,
                                     const string& fmt_x,
                                     const string& fmt_y,
                                     const string& unit_x,
                                     const string& unit_y,
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
      fmt_x, NUMBER_REAL, 'c', 't', 9);
   mesh_2d.render_label_y (cr, transform, 0, start_x,
      fmt_y, NUMBER_REAL, 'r', 'c', 9);

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

void
Image_Package::image_journey (const string& image_identifier,
                              const string& geodetic_transform_identifier,
                              const string& journey_identifier)
{

   const RefPtr<ImageSurface>& image_surface = image_map.at (image_identifier);
   const RefPtr<Context> cr = cr_map.at (image_identifier);
   const Integer w = image_surface->get_width ();
   const Integer h = image_surface->get_height ();
   const Point_2D centre (Real (w) / 2, Real (h) / 2);

   const Geodetic_Transform* geodetic_transform_ptr =
      andrea.get_geodetic_transform_ptr (geodetic_transform_identifier, centre);
   const Geodetic_Transform& geodetic_transform = *geodetic_transform_ptr;

   const Journey& journey = andrea.get_journey (journey_identifier);

   cr->save ();
   Color::black ().cairo (cr);
   journey.cairo (cr, geodetic_transform);
   cr->restore ();

   delete geodetic_transform_ptr;

}

void
Image_Package::image_geodetic_mesh (const string& image_identifier,
                                    const string& geodetic_transform_identifier,
                                    const string& geodetic_mesh_identifier)
{

   const RefPtr<ImageSurface>& image_surface = image_map.at (image_identifier);
   const RefPtr<Context> cr = cr_map.at (image_identifier);
   const Integer w = image_surface->get_width ();
   const Integer h = image_surface->get_height ();
   const Point_2D centre (Real (w) / 2, Real (h) / 2);

   const Geodetic_Transform* geodetic_transform_ptr =
      andrea.get_geodetic_transform_ptr (geodetic_transform_identifier, centre);
   const Geodetic_Transform& geodetic_transform = *geodetic_transform_ptr;

   const Geodetic_Mesh& geodetic_mesh =
      andrea.get_geodetic_mesh (geodetic_mesh_identifier);

   cr->save ();
   Color::black ().cairo (cr);
   geodetic_mesh.cairo (cr, geodetic_transform);
   cr->stroke ();
   cr->restore ();

   delete geodetic_transform_ptr;

}

void
Image_Package::image_gshhs (const string& image_identifier,
                            const string& geodetic_transform_identifier,
                            const string& gshhs_identifier,
                            const Tokens& arguments)
{

   const RefPtr<ImageSurface>& image_surface = image_map.at (image_identifier);
   const RefPtr<Context> cr = cr_map.at (image_identifier);
   const Integer w = image_surface->get_width ();
   const Integer h = image_surface->get_height ();
   const Point_2D centre (Real (w) / 2, Real (h) / 2);

   const Geodetic_Transform* geodetic_transform_ptr =
      andrea.get_geodetic_transform_ptr (geodetic_transform_identifier, centre);
   const Geodetic_Transform& geodetic_transform = *geodetic_transform_ptr;

   const Gshhs& gshhs = *(andrea.get_gshhs_ptr (gshhs_identifier));

   const Integer n = arguments.size ();
   const Color& color = (n > 0 ? Color (arguments[0]) : Color::black ());
   const bool is_fill = (n > 1 && arguments[1] == "fill");

   cr->save ();
   color.cairo (cr);
   gshhs.cairo (cr, geodetic_transform);
   if (is_fill) { cr->fill (); } else { cr->stroke (); }
   cr->restore ();

   delete geodetic_transform_ptr;
   
}

void
Image_Package::image_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "init")
   {
      const string& identifier = tokens[1];
      const string& geometry = tokens[2];
      image_init (identifier, geometry);
   }
   else
   if (tokens[0] == "paint")
   {
      const string& identifier = tokens[1];
      const Color& color = (n > 2 ? Color (tokens[2]) : Color::white ());
      image_paint (identifier, color);
   }
   else
   if (tokens[0] == "save")
   {
      const string& identifier = tokens[1];
      const string& file_path = tokens[2];
      image_save (identifier, file_path);
   }
   else
   if (tokens[0] == "title")
   {
      const string& image_identifier = tokens[1];
      image_title (image_identifier, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "sounding")
   {
      image_sounding (tokens.subtokens (1));
   }
   else
   if (tokens[0] == "journey")
   {
      //const string& image_identifier = tokens[1];
      //const string& geodetic_transform_identifier = tokens[2];
      //const string& journey_identifier = tokens[3];
      //image_journey (image_identifier, geodetic_transform_identifier,
      //   journey_identifier);
   }
   else
   if (tokens[0] == "geodetic_mesh")
   {
      const string& image_identifier = tokens[1];
      const string& geodetic_transform_identifier = tokens[2];
      const string& geodetic_mesh_identifier = tokens[3];
      image_geodetic_mesh (image_identifier, geodetic_transform_identifier,
         geodetic_mesh_identifier);
   }
   else
   if (tokens[0] == "gshhs")
   {
      const string& image_identifier = tokens[1];
      const string& geodetic_transform_identifier = tokens[2];
      const string& gshhs_identifier = tokens[3];
      image_gshhs (image_identifier, geodetic_transform_identifier,
         gshhs_identifier, tokens.subtokens (4));
   }

}

const map<string, RefPtr<ImageSurface> >&
Image_Package::get_image_map () const
{
   return image_map;
}

const map<string, RefPtr<Context> >&
Image_Package::get_cr_map () const
{
   return cr_map;
}

const RefPtr<ImageSurface>&
Image_Package::get_image (const string& identifier) const
{
   return image_map.at (identifier);
}

const RefPtr<Context>&
Image_Package::get_cr (const string& identifier) const
{
   return cr_map.at (identifier);
}

Entity::Entity (const string& str)
   : string (str)
{
}

Real
Entity::value () const
{

   const bool addition = Reg_Exp::match (*this, "\\+");
   const bool subtraction = Reg_Exp::match (*this, "\\-");
   const bool multiplication = Reg_Exp::match (*this, "\\*");
   const bool division = Reg_Exp::match (*this, "\\/");
   const bool power = Reg_Exp::match (*this, "\\^");

   if (addition)
   {
      const Tokens tokens (*this, "+");
      Real value = 0;
      for (auto iterator = tokens.begin ();
           iterator != tokens.end (); iterator++)
      {
         const Entity entity (*iterator);
         value += entity.value ();
      }
      return value;
   }
   else
   if (subtraction)
   {
      const Tokens tokens (*this, "-");
      Real value = 0;
      for (auto iterator = tokens.begin ();
           iterator != tokens.end (); iterator++)
      {
         const Entity entity (*iterator);
         value -= entity.value ();
      }
      return value;
   }
   else
   if (multiplication)
   {
      const Tokens tokens (*this, "*");
      Real value = 1;
      for (auto iterator = tokens.begin ();
           iterator != tokens.end (); iterator++)
      {
         const Entity entity (*iterator);
         value *= entity.value ();
      }
      return value;
   }
   else
   if (division)
   {
      const Tokens tokens (*this, "/");
      Real value = 1;
      for (auto iterator = tokens.begin ();
           iterator != tokens.end (); iterator++)
      {
         const Entity entity (*iterator);
         value /= entity.value ();
      }
      return value;
   }
   else
   if (power)
   {
      const Integer i = this->find ("^");
      const Entity entity_a (this->substr (0, i));
      const Entity entity_b (this->substr (i + 1));
      return pow (entity_a.value (), entity_b.value ());
   }

   return stof (*this);

}

void
Andrea::wind_shear (const Tokens& arguments) const
{

   const Integer n = arguments.size ();
   if (n != 2) { throw Exception ("wind_shear"); }

   const Tokens tokens_a (arguments[0], ":");
   const Integer n_a = tokens_a.size ();
   if (n_a < 1 || n_a > 2 ) { throw Exception ("wind_shear"); }
   const bool a_with_height = (n_a == 2);
   const Tokens wind_tokens_a (tokens_a[a_with_height ? 1 : 0], "/");
   if (wind_tokens_a.size () != 2) { throw Exception ("wind_shear"); }

   const Tokens tokens_b (arguments[1], ":");
   const Integer n_b = tokens_b.size ();
   if (n_b < 1 || n_b > 2 ) { throw Exception ("wind_shear"); }
   const bool b_with_height = (n_b == 2);
   const Tokens wind_tokens_b (tokens_b[b_with_height ? 1 : 0], "/");
   if (wind_tokens_b.size () != 2) { throw Exception ("wind_shear"); }

   const Real height_a    = a_with_height ? stof (tokens_a[0]) : GSL_NAN;
   const Real direction_a = stof (wind_tokens_a[0]);
   const Real speed_a     = stof (wind_tokens_a[1]);

   const Real height_b    = b_with_height ? stof (tokens_b[0]) : GSL_NAN;
   const Real direction_b = stof (wind_tokens_b[0]);
   const Real speed_b     = stof (wind_tokens_b[1]);

   const Wind& wind_a = Wind::direction_speed (direction_a, speed_a);
   const Wind& wind_b = Wind::direction_speed (direction_b, speed_b);

   const Real du = wind_a.u - wind_b.u;
   const Real dv = wind_a.v - wind_b.v;
   const Wind shear (du, dv);

   const Real d_height = height_a - height_b;
   const bool per_height = gsl_finite (d_height);
   const Real shear_direction = shear.get_direction ();
   const Real shear_speed = shear.get_speed () / (per_height ? d_height : 1);

   cout << shear_direction << "/" << shear_speed << endl;

}

void
Andrea::print (const Entity& entity) const
{
   cout << entity.value () << endl;
}

Andrea::Andrea ()
   : Image_Package (*this),
     Geodesy_Package (*this),
     Gshhs_Package (*this),
     Journey_Package (*this),
     Sounding_Package (*this),
     Geodetic_Mesh_Package (*this),
     Geodetic_Transform_Package (*this)
{
}

void
Andrea::parse (const Tokens& tokens)
{

   if (get_lower_case (tokens[0]) == "image")
   {
      image_parse (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "journey")
   {
      journey_parse (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "geodesy")
   {
      geodesy_parse (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "sounding")
   {
      sounding_parse (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "geodetic_mesh")
   {
      geodetic_mesh_parse (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "geodetic_transform")
   {
      geodetic_transform_parse (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "gshhs")
   {
      gshhs_parse (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "shear" ||
       get_lower_case (tokens[0]) == "wind_shear")
   {
      wind_shear (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "print")
   {
      print (tokens[1]);
      return;
   }

   throw Exception ("parse");

}

