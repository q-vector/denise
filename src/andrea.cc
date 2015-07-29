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
Journey_Package::assign_journey (const string& variable,
                                 const Tokens& arguments)
{

   const Geodesy geodesy;
   const Integer n = arguments.size ();

   switch (n)
   {

      case 2:
      {
         const Lat_Long origin (arguments[0]);
         const Lat_Long destination (arguments[1]);
         Journey journey (origin, destination);
         geodesy.complete (journey);
         journey_map[variable] = journey;
         break;
      }
         
      case 3:
      {
         const Lat_Long origin (arguments[0]);
         const Real distance = stof (arguments[1]);
         const Real azimuth_forward = stof (arguments[2]);
         Journey journey (origin, distance, azimuth_forward);
         geodesy.complete (journey);
         journey_map[variable] = journey;
         break;
      }

      default:
      {
         throw Exception ("assign_journey");
      }

   }

}

void
Journey_Package::print_journey (const string& variable) const
{

   const Journey& journey = journey_map.at (variable);

   cout << journey.get_origin ().get_string (lat_long_dp) << " ";
   cout << journey.get_destination ().get_string (lat_long_dp) << " ";
   cout << journey.get_distance () << " ";
   cout << journey.get_azimuth_forward () << " ";
   cout << journey.get_azimuth_backward () << endl;

}

void
Journey_Package::journey_distance (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 2) { throw Exception ("journey_distance"); }

   const Lat_Long origin (tokens[0]);
   const Lat_Long destination (tokens[1]);
   const Real distance = Geodesy::get_distance (origin, destination);

   cout << distance << endl;

}

void
Journey_Package::journey_azimuth (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 2) { throw Exception ("journey_azimuth"); }

   const Lat_Long origin (tokens[0]);
   const Lat_Long destination (tokens[1]);
   const Real azimuth_f = Geodesy::get_azimuth_forward (origin, destination);
   const Real azimuth_b = Geodesy::get_azimuth_backward (origin, destination);

   cout << azimuth_f << " " << azimuth_b << endl;

}

void
Journey_Package::journey_destination (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 3) { throw Exception ("journey_destination"); }

   const Lat_Long origin (tokens[0]);
   const Real distance = stof (tokens[1]);
   const Real azimuth_forward = stof (tokens[2]);
   const Lat_Long& destination =
      Geodesy::get_destination (origin, distance, azimuth_forward);

   cout << destination.get_string (lat_long_dp) << endl;

}

Journey_Package::Journey_Package (const Andrea& andrea)
   : Andrea_Package (andrea),
     lat_long_dp (4)
{
}

void
Journey_Package::parse_journey (const Tokens& tokens)
{


   const Integer n = tokens.size ();

   if (tokens[0] == "assign")
   {
      const string& variable = tokens[1];
      assign_journey (variable, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "print")
   {
      const string& variable = tokens[1];
      print_journey (variable);
   }
   else
   if (tokens[0] == "distance")
   {
      journey_distance (tokens.subtokens (1));
   }
   else
   if (tokens[0] == "azimuth")
   {
      journey_azimuth (tokens.subtokens (1));
   }
   else
   if (tokens[0] == "destination")
   {
      journey_destination (tokens.subtokens (1));
   }
   else
   {
      throw Exception ("parse_journey");
   }

}

const map<string, Journey>&
Journey_Package::get_journey_map () const
{
   return journey_map;
}

Sounding_Package::Sounding_Package (const Andrea& andrea)
   : Andrea_Package (andrea)
{
}

void
Sounding_Package::load_sounding (const string& variable,
                                 const string& file_path)
{
   Sounding sounding (file_path);
   sounding_map[variable] = sounding;
}

void
Sounding_Package::print_sounding (const string& variable,
                                  const Tokens& tokens) const
{

   const string time_fmt ("%Y%m%d%H%M");
   const Sounding& sounding = sounding_map.at (variable);
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
Sounding_Package::parse_sounding (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "load")
   {
      const string& variable = tokens[1];
      const string& file_path = tokens[2];
      load_sounding (variable, file_path);
   }
   else
   if (tokens[0] == "print")
   {
      const string& variable = tokens[1];
      print_sounding (variable, tokens.subtokens (2));
   }

}

const map<string, Sounding>&
Sounding_Package::get_sounding_map () const
{
   return sounding_map;
}

void
Image_Package::init_image (const string& variable,
                           const string& geometry)
{

   const Tokens tokens (geometry, "x");
   const Size_2D size_2d (stof (tokens[0]), stof (tokens[1]));

   RefPtr<ImageSurface> image_surface = get_image_surface (size_2d);
   image_map[variable] = image_surface;

}

void
Image_Package::save_image (const string& variable,
                           const string& file_path) const
{
   const RefPtr<ImageSurface>& image_surface = image_map.at (variable);
   image_surface->write_to_png (file_path);
}

void
Image_Package::image_title (const string& variable,
                            const Tokens& tokens) const
{

   const RefPtr<ImageSurface>& image_surface = image_map.at (variable);
   const RefPtr<Context> cr = denise::get_cr (image_surface);
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
   const RefPtr<Context> cr = denise::get_cr (image_surface);
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
   const string& genre = tokens[4]; // always speed for now

   const RefPtr<ImageSurface>& image_surface = image_map.at (image_identifier);
   const RefPtr<Context> cr = denise::get_cr (image_surface);
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

   const Real minor_interval_x = 1;
   const Real major_interval_x = 10;
   const Real minor_interval_y = is_p ? 10e2 : 100;
   const Real major_interval_y = is_p ? 100e2 : 1000;
   const Color minor_color = Color::black (0.2);
   const Color major_color = Color::black (0.8);

   const string y_fmt (is_p ? "%.0f Pa" : "%.0f m");
   const string x_fmt ("%.0f ms\u207b\u00b9");

   const Mesh_2D mesh_2d (Size_2D (2, 2), domain_2d,
      major_interval_x, major_interval_y, major_color,
      minor_interval_x, minor_interval_y, minor_color);

   Affine_Transform_2D transform;
   const Real span_x = domain_x.get_span ();
   const Real span_y = domain_y.get_span ();
   transform.scale (1, -1);
   transform.scale (w / span_x, h / span_y);
   transform.translate (margin_l, size_2d.j - margin_b);

   Color::white ().cairo (cr);
   cr->paint ();

   cr->set_line_width (2);
   Color::black ().cairo (cr);
   mesh_2d.render (cr, transform);
   mesh_2d.render_label_x (cr, transform, domain_x.start, domain_y.start,
      x_fmt, NUMBER_REAL, 'c', 't', 5);
   mesh_2d.render_label_y (cr, transform, domain_x.start, domain_y.start,
      y_fmt, NUMBER_REAL, 'r', 'c', 5);

   const Ring ring (4);
   const Wind_Profile& wind_profile = sounding.get_wind_profile ();

   for (auto iterator = wind_profile.begin ();
        iterator != wind_profile.end (); iterator++)
   {

      const Real p = iterator->first;
      const Real z = sounding.get_height (p);
      const Wind& wind = iterator->second;
      const Real speed = wind.get_speed ();
      const Point_2D point = transform.transform (Point_2D (speed, z));

      ring.cairo (cr, point);
      Color::red (0.4).cairo (cr);
      cr->fill ();

   }

}

Image_Package::Image_Package (const Andrea& andrea)
   : Andrea_Package (andrea)
{
}

void
Image_Package::parse_image (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "init")
   {
      const string& variable = tokens[1];
      const string& geometry = tokens[2];
      init_image (variable, geometry);
   }
   else
   if (tokens[0] == "save")
   {
      const string& variable = tokens[1];
      const string& file_path = tokens[2];
      save_image (variable, file_path);
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

}

const map<string, RefPtr<ImageSurface> >&
Image_Package::get_image_map () const
{
   return image_map;
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
   if (n != 2) { throw Exception ("wind_shear error"); }

   const Tokens tokens_a (arguments[0], ":");
   const Integer n_a = tokens_a.size ();
   if (n_a < 1 || n_a > 2 ) { throw Exception ("wind_shear error"); }
   const bool a_with_height = (n_a == 2);
   const Tokens wind_tokens_a (tokens_a[a_with_height ? 1 : 0], "/");
   if (wind_tokens_a.size () != 2) { throw Exception ("wind_shear error"); }

   const Tokens tokens_b (arguments[1], ":");
   const Integer n_b = tokens_b.size ();
   if (n_b < 1 || n_b > 2 ) { throw Exception ("wind_shear error"); }
   const bool b_with_height = (n_b == 2);
   const Tokens wind_tokens_b (tokens_b[b_with_height ? 1 : 0], "/");
   if (wind_tokens_b.size () != 2) { throw Exception ("wind_shear error"); }

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
     Journey_Package (*this),
     Sounding_Package (*this)
{
}

void
Andrea::parse (const Tokens& tokens)
{

   if (get_lower_case (tokens[0]) == "image")
   {
      parse_image (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "journey")
   {
      parse_journey (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "sounding")
   {
      parse_sounding (tokens.subtokens (1));
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

   throw Exception ("cannot process: ");

}

