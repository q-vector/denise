//
// dstring.h
// 
// Copyright (C) 2005-2013 Simon E. Ching
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

#ifndef DENISE_DSTRING_H
#define DENISE_DSTRING_H

#include <ctype.h>
#include <cstdarg>
#include <algorithm>
#include <string>
#include <cstring>
#include <vector>
#include <regex.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <denise/basics.h>
#include <denise/exception.h>

using namespace std;

namespace denise
{

   class Reg_Exp
   {

      private:

         regex_t
         preg;

         Integer
         nmatch;

         regmatch_t*
         pmatch;

      public:

         Reg_Exp (const string& reg_exp_str,
                  const bool match_info = false,
                  const bool posix_extend = true,
                  const bool case_sensitive = true);

         ~Reg_Exp ();

         bool
         match (const string& str) const;

         Iduple
         get_match (const string& str,
                    const bool ignore_no_match = true) const;

         void
         replace (string& str,
                  const string& with,
                  const bool ignore_no_match = true) const;

         Iduple_Vector
         get_match_sub (const string& str,
                        const bool ignore_no_match = true) const;

         static bool
         match (const string& str,
                const string& reg_exp_str,
                const bool posix_extend = true,
                const bool case_sensitive = true);

   };

   class Tokens : public vector<string>
   {

      public:

         Tokens ();

         Tokens (const Tuple& tuple,
                 const string& fmt);

         Tokens (const string& str,
                 const string& delimiters = string (" \f\n\t"));

         Real
         real (const Integer index) const;

         Integer
         integer (const Integer index) const;

         void
         add (const string& str,
              const string& delimiters = string (" \f\n\t"));

         void
         add_prefix (const string& prefix);

         void
         add_suffix (const string& suffix);

   };

   /// Returns a tokenized vector of string
   vector<string>
   tokenize (const string& s,
             const string& delimiters = string (" \f\n\t"));

   /// Returns a tokenized vector of string with seperate delimiters
   vector<string>
   tokenize_s (const string& s,
               const string& delimiters = string (" \f\n\t"));

   /// Converts the given string into lower case.
   void
   to_lower_case (string& s);

   string
   get_lower_case (const string& s);

   /// Converts the given string into upper case.
   void
   to_upper_case (string& s);

   string
   get_upper_case (const string& s);

   /// Converts the given string into capital case.
   void
   to_capital_case (string& s);

   string
   get_captial_case (const string& s);

   /// Chops the last character off the string.
   void
   chop (string& s);

   /// Trims white spaces
   void
   trim (string& s,
         const string& white_string = string (" \f\n\t"));

   /// Trims left white spaces
   void
   left_trim (string& s,
              const string& white_string = string (" \f\n\t"));

   /// Trims right white spaces
   void
   right_trim (string& s,
               const string& white_string = string (" \f\n\t"));

   string
   get_trimmed (const string s,
                const string& white_string = string (" \f\n\t"));

   string
   get_left_trimmed (const string s, 
                     const string& white_string = string (" \f\n\t"));

   string
   get_right_trimmed (const string s,
                      const string& white_string = string (" \f\n\t"));

   /// Returns file extension
   string
   get_file_extension (const string& file_path);

   /// Returns a rendered std::string similar to output of printf
   string
   string_render (const char* format,
                  ...);

   ostream&
   operator << (ostream &out_file,
                const Tokens& tokens);

}


#endif /* DENISE_DSTRING_H */

