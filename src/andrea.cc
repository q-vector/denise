//
// andrea.cc
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

#include <cstdio>
#include <readline/readline.h>
#include <readline/history.h>
#include "andrea.h"

using namespace std;
using namespace denise;
using namespace andrea;

Entity::Entity (const Dstring& str)
   : Dstring (str)
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
Andrea::data_path (const Tokens& arguments)
{

   const Dstring& a = arguments[0].get_lower_case ();
   const Dstring& identifier = arguments[1].get_lower_case ();

   if (a == "assign")
   {
      const Dstring& data_path = arguments[2];
      data_path_map[identifier] = data_path;
   }
   else
   if (a == "print")
   {
      cout << data_path_map[identifier];
   }

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

Andrea::Andrea (const Dstring& prompt)
   : prompt (prompt),
     Best_Tracks_Package (*this),
     Geodesy_Package (*this),
     Geodetic_Mesh_Package (*this),
     Geodetic_Transform_Package (*this),
     Gshhs_Package (*this),
     Journey_Package (*this),
     Sounding_Package (*this),
     Surface_Package (*this)
{
}

const Dstring&
Andrea::get_data_path (const Dstring& identifier) const
{
   return data_path_map.at (identifier);
}

void
Andrea::parse (const Tokens& tokens)
{

   const Dstring& action = tokens[0].get_lower_case ();

   if (action == "data_path")
   {
      data_path (tokens.subtokens (1));
      return;
   }
   else
   if (action == "best_track")
   {
      best_tracks_parse (tokens.subtokens (1));
      return;
   }
   else
   if (action == "journey")
   {
      journey_parse (tokens.subtokens (1));
      return;
   }
   else
   if (action == "geodesy")
   {
      geodesy_parse (tokens.subtokens (1));
      return;
   }
   else
   if (action == "geodetic_mesh")
   {
      geodetic_mesh_parse (tokens.subtokens (1));
      return;
   }
   else
   if (action == "geodetic_transform")
   {
      geodetic_transform_parse (tokens.subtokens (1));
      return;
   }
   else
   if (action == "gshhs")
   {
      gshhs_parse (tokens.subtokens (1));
      return;
   }
   else
   if (action == "sounding")
   {
      sounding_parse (tokens.subtokens (1));
      return;
   }
   else
   if (action == "surface")
   {
      surface_parse (tokens.subtokens (1));
      return;
   }
   else
   if (action == "shear" ||
       action == "wind_shear")
   {
      wind_shear (tokens.subtokens (1));
      return;
   }
   else
   if (action == "print")
   {
      print (tokens[1]);
      return;
   }
   else
   if (action == "set")
   {
      if (tokens[1].get_lower_case () == "prompt")
      {
        prompt = tokens[2];
      }
      return;
   }

   throw Exception ("Unknown action: " + action);

}

void
Andrea::loop ()
{

   char* line;

   for (bool done = false ; !done; )
   {

      line = readline (prompt.c_str ());
      if (!line) { break; }

      Dstring input_line = Dstring (line).get_trimmed ();
      if (input_line[0] == '#') { continue; }

      try
      {

         // This has to come first before assign_tokens
         //    to allow for recursive variable declarations
         apply_variables (input_line);

         // This has to come second after apply_variables
         //    to allow for recursive variable declarations
         const Tokens assign_tokens (input_line, "= ");
         if (assign_tokens.size () == 2)
         {
            const Dstring& variable = assign_tokens[0].get_trimmed ();
            const Dstring& value = assign_tokens[1].get_trimmed ();
            assign_dictionary.insert (make_pair (variable, value));
            continue;
         }

         const Tokens tokens (input_line, " \f\n\t");
         if (tokens.size () == 0) { continue; }

         parse (tokens);

      }
      catch (const Exception& e)
      {
         wcerr << e.what () << endl;
      }
      catch (const std::out_of_range& oor)
      {
         wcerr << "andrea: out_of_range " << oor.what () <<
            " " << line << endl;
      }

      add_history (line);
      free (line);

   }

}

