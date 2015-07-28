#include <denise/dstring.h>
#include <denise/met.h>
#include "andrea.h"

using namespace std;
using namespace denise;

void
Sounding_Map::print (const string& variable,
                     const Tokens& tokens) const
{

   const string time_fmt ("%Y%m%d%H%M");
   const Integer n = tokens.size ();
   const Sounding& sounding = at (variable);

   if (n < 1) { return; }
   if (tokens[0] == "time")
   {
      cout << sounding.get_time ().get_string (time_fmt) << endl;
      return;
   }
   else
   if (tokens[1] == "basetime")
   {
      cout << sounding.get_basetime ().get_string (time_fmt) << endl;
      return;
   }
   else
   if (tokens[1] == "location")
   {
      cout << sounding.get_location_str () << endl;
      return;
   }
   else
   if (tokens[1] == "t")
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
   if (tokens[1] == "td")
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
   if (tokens[1] == "wind")
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
   if (tokens[1] == "height")
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
Andrea::assign_journey (const string& variable,
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
         throw Exception ("assign_journey error");
      }

   }

}

void
Andrea::assign_sounding (const string& variable,
                         const Tokens& arguments)
{

   const Integer n = arguments.size ();

   switch (n)
   {

      case 1:
      {
         const string& file_path = arguments[0];
         Sounding sounding (file_path);
         sounding_map[variable] = sounding;
         break;
      }
         
      default:
      {
         throw Exception ("assign_journey error");
      }

   }

}

void
Andrea::journey (const Tokens& arguments)
{

   const Integer n = arguments.size ();

   const bool is_assignment = (n >= 2) && (arguments[1] == "=");
   const bool is_print = (n == 1);
   const bool is_evaluate = (n >= 2) && (arguments[1] != "=");

   if (is_assignment)
   {
      const string& variable = arguments[0];
      assign_journey (arguments[0], arguments.subtokens (2));
   }
   else
   if (is_evaluate)
   {

      Journey journey;
      const Geodesy geodesy;

      switch (n)
      {

         case 2:
         {
            const Lat_Long origin (arguments[0]);
            const Lat_Long destination (arguments[1]);
            journey = Journey (origin, destination);
            geodesy.complete (journey);
            cout << journey.get_distance () << " " <<
                    journey.get_azimuth_forward () << " " <<
                    journey.get_azimuth_backward () << endl;
            break;
         }
         
         case 3:
         {
            const Lat_Long origin (arguments[0]);
            const Real distance = stof (arguments[1]);
            const Real azimuth_forward = stof (arguments[2]);
            journey = Journey (origin, distance, azimuth_forward);
            geodesy.complete (journey);
            cout << journey.get_destination ().get_string (lat_long_dp) <<
               " " << journey.get_azimuth_backward () << endl;
            break;
         }

         default:
         {
            throw Exception ("journey error");
         }

      }

   }
   else
   if (is_print)
   {
      const string& variable = arguments[0];
      const Journey& journey = journey_map.at (variable);
      cout << journey.get_origin ().get_string (lat_long_dp) << " " <<
              journey.get_destination ().get_string (lat_long_dp) << " " <<
              journey.get_distance () << " " <<
              journey.get_azimuth_forward () << " " <<
              journey.get_azimuth_backward () << endl;
   }

}

void
Andrea::sounding (const Tokens& arguments)
{

   const Integer n = arguments.size ();

   const bool is_assignment = (n >= 2) && (arguments[1] == "=");
   const bool is_print = (n == 1);

   if (is_assignment)
   {
      const string& variable = arguments[0];
      assign_sounding (arguments[0], arguments.subtokens (2));
      return;
   }

   sounding_map.print (arguments[0], arguments[1]);

}

void
Andrea::distance (const Tokens& arguments) const
{

   const Integer n = arguments.size ();
   if (n != 2) { throw Exception ("distance error"); }

   const Lat_Long origin (arguments[0]);
   const Lat_Long destination (arguments[1]);
   cout << Geodesy::get_distance (origin, destination) << endl;

}

void
Andrea::azimuth (const Tokens& arguments) const
{

   const Integer n = arguments.size ();
   if (n != 2) { throw Exception ("distance error"); }

   const Lat_Long origin (arguments[0]);
   const Lat_Long destination (arguments[1]);
   cout << Geodesy::get_azimuth_forward (origin, destination) << " " <<
           Geodesy::get_azimuth_backward (origin, destination) << endl;

}

void
Andrea::destination (const Tokens& arguments) const
{

   const Integer n = arguments.size ();
   if (n != 3) { throw Exception ("distance error"); }

   const Lat_Long origin (arguments[0]);
   const Real distance = stof (arguments[1]);
   const Real azimuth_forward = stof (arguments[2]);
   cout << Geodesy::get_destination (origin, distance, azimuth_forward) << endl;

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
   : lat_long_dp (4)
{
}

void
Andrea::parse (const Tokens& tokens)
{

   if (get_lower_case (tokens[0]) == "journey")
   {
      journey (tokens.subtokens (1));
      return;
   }

   if (get_lower_case (tokens[0]) == "sounding")
   {
      sounding (tokens.subtokens (1));
      return;
   }

   if (get_lower_case (tokens[0]) == "distance")
   {
      distance (tokens.subtokens (1));
      return;
   }

   if (get_lower_case (tokens[0]) == "azimuth")
   {
      azimuth (tokens.subtokens (1));
      return;
   }

   if (get_lower_case (tokens[0]) == "destination")
   {
      destination (tokens.subtokens (1));
      return;
   }

   if (get_lower_case (tokens[0]) == "shear" ||
       get_lower_case (tokens[0]) == "wind_shear")
   {
      wind_shear (tokens.subtokens (1));
      return;
   }

   if (get_lower_case (tokens[0]) == "print")
   {
      print (tokens[1]);
      return;
   }

   throw Exception ("cannot process: ");

}

