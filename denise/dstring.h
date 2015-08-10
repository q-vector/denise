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

         Reg_Exp (const wstring& reg_exp_str,
                  const bool match_info = false,
                  const bool posix_extend = true,
                  const bool case_sensitive = true);

         ~Reg_Exp ();

         bool
         match (const wstring& str) const;

         Iduple
         get_match (const wstring& str,
                    const bool ignore_no_match = true) const;

         void
         replace (wstring& str,
                  const wstring& with,
                  const bool ignore_no_match = true) const;

         Iduple_Vector
         get_match_sub (const wstring& str,
                        const bool ignore_no_match = true) const;

         static bool
         match (const wstring& str,
                const wstring& reg_exp_str,
                const bool posix_extend = true,
                const bool case_sensitive = true);

   };

   class Tokens : public vector<wstring>
   {

      public:

         Tokens ();

         Tokens (const Tuple& tuple,
                 const wstring& fmt);

         Tokens (const wstring& str);

         Tokens (const wstring& str,
                 const wstring& delimiters);

         Tokens
         subtokens (const Integer i,
                    const Integer n = -1) const;

         Real
         real (const Integer index) const;

         Integer
         integer (const Integer index) const;

         void
         add (const wstring& str);

         void
         add (const wstring& str,
              const wstring& delimiters);

         void
         add_prefix (const wstring& prefix);

         void
         add_suffix (const wstring& suffix);

   };

   /// Converts the given string into lower case.
   void
   to_lower_case (wstring& s);

   wstring
   get_lower_case (const wstring& s);

   /// Converts the given string into upper case.
   void
   to_upper_case (wstring& s);

   wstring
   get_upper_case (const wstring& s);

   /// Converts the given string into capital case.
   void
   to_capital_case (wstring& s);

   wstring
   get_captial_case (const wstring& s);

   /// Chops the last character off the string.
   void
   chop (wstring& s);

   /// Trims white spaces
   void
   trim (wstring& s,
         const wstring& white_string = wstring (L" \f\n\t"));

   /// Trims left white spaces
   void
   left_trim (wstring& s,
              const wstring& white_string = wstring (L" \f\n\t"));

   /// Trims right white spaces
   void
   right_trim (wstring& s,
               const wstring& white_string = wstring (L" \f\n\t"));

   wstring
   get_trimmed (const wstring s,
                const wstring& white_string = wstring (L" \f\n\t"));

   wstring
   get_left_trimmed (const wstring s, 
                     const wstring& white_string = wstring (L" \f\n\t"));

   wstring
   get_right_trimmed (const wstring s,
                      const wstring& white_string = wstring (L" \f\n\t"));

   /// Returns file extension
   wstring
   get_file_extension (const wstring& file_path);

   /// Returns a rendered std::string similar to output of printf
   wstring
   string_render (const char* format,
                  ...);

   wostream&
   operator << (wostream &out_file,
                const Tokens& tokens);

}


#endif /* DENISE_DSTRING_H */

