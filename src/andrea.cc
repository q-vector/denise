#include <denise/dstring.h>
#include <denise/met.h>
#include "andrea.h"

using namespace std;
using namespace denise;

void
Journey_Map::assign (const string& variable,
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
         (*this)[variable] = journey;
         break;
      }
         
      case 3:
      {
         const Lat_Long origin (arguments[0]);
         const Real distance = stof (arguments[1]);
         const Real azimuth_forward = stof (arguments[2]);
         Journey journey (origin, distance, azimuth_forward);
         geodesy.complete (journey);
         (*this)[variable] = journey;
         break;
      }

      default:
      {
         throw Exception ("Journey assign");
      }

   }

}

void
Journey_Map::print (const string& variable) const
{

   const Journey& journey = at (variable);

   cout << journey.get_origin ().get_string (lat_long_dp) << " ";
   cout << journey.get_destination ().get_string (lat_long_dp) << " ";
   cout << journey.get_distance () << " ";
   cout << journey.get_azimuth_forward () << " ";
   cout << journey.get_azimuth_backward () << endl;

}

void
Journey_Map::distance (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 2) { throw Exception ("Journey distance"); }

   const Lat_Long origin (tokens[0]);
   const Lat_Long destination (tokens[1]);
   const Real distance = Geodesy::get_distance (origin, destination);

   cout << distance << endl;

}

void
Journey_Map::azimuth (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 2) { throw Exception ("Journey azimuth"); }

   const Lat_Long origin (tokens[0]);
   const Lat_Long destination (tokens[1]);
   const Real azimuth_f = Geodesy::get_azimuth_forward (origin, destination);
   const Real azimuth_b = Geodesy::get_azimuth_backward (origin, destination);

   cout << azimuth_f << " " << azimuth_b << endl;

}

void
Journey_Map::destination (const Tokens& tokens) const
{

   const Integer n = tokens.size ();
   if (n != 3) { throw Exception ("Journey destination"); }

   const Lat_Long origin (tokens[0]);
   const Real distance = stof (tokens[1]);
   const Real azimuth_forward = stof (tokens[2]);
   const Lat_Long& destination =
      Geodesy::get_destination (origin, distance, azimuth_forward);

   cout << destination.get_string (lat_long_dp) << endl;

}

Journey_Map::Journey_Map ()
   : lat_long_dp (4)
{
}

void
Journey_Map::parse (const Tokens& tokens)
{


   const Integer n = tokens.size ();

   if (tokens[0] == "assign")
   {
      const string& variable = tokens[1];
      assign (variable, tokens.subtokens (2));
   }
   else
   if (tokens[0] == "print")
   {
      const string& variable = tokens[1];
      print (variable);
   }
   else
   if (tokens[0] == "distance")
   {
      distance (tokens.subtokens (1));
   }
   else
   if (tokens[0] == "azimuth")
   {
      azimuth (tokens.subtokens (1));
   }
   else
   if (tokens[0] == "destination")
   {
      destination (tokens.subtokens (1));
   }
   else
   {
      throw Exception ("Journey");
   }

}

void
Sounding_Map::load (const string& variable,
                    const string& file_path)
{
   Sounding sounding (file_path);
   (*this)[variable] = sounding;
}

void
Sounding_Map::print (const string& variable,
                     const Tokens& tokens) const
{

   const string time_fmt ("%Y%m%d%H%M");
   const Sounding& sounding = at (variable);
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

Sounding_Map::Sounding_Map ()
   : tephigram (Size_2D (1000, 1000))
{
}

void
Sounding_Map::parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "load")
   {
      const string& variable = tokens[1];
      const string& file_path = tokens[2];
      load (variable, file_path);
   }
   else
   if (tokens[0] == "print")
   {
      const string& variable = tokens[1];
      print (variable, tokens.subtokens (2));
   }

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
{
}

void
Andrea::parse (const Tokens& tokens)
{

   if (get_lower_case (tokens[0]) == "journey")
   {
      journey_map.parse (tokens.subtokens (1));
      return;
   }
   else
   if (get_lower_case (tokens[0]) == "sounding")
   {
      sounding_map.parse (tokens.subtokens (1));
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

