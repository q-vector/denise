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
#include <denise/dstring.h>
#include <denise/exception.h>

using namespace std;

namespace denise
{

   class Dstring : public string
   {

      public:

         Dstring ();

         Dstring (const Dstring& dstring);

         Dstring (const string& str);

         Dstring (const char* buffer);

         string
         get_string () const;

         void
         to_lower_case ();

         void
         to_upper_case ();

         void
         to_capital_case ();

         void
         chop ();

         void
         trim (const Dstring& white_string = " \f\n\t");

         void
         left_trim (const Dstring& white_string = " \f\n\t");

         void
         right_trim (const Dstring& white_string = " \f\n\t");

         Dstring
         get_lower_case () const;

         Dstring
         get_upper_case () const;

         Dstring
         get_captial_case () const;

         Dstring
         get_trimmed (const Dstring& white_string = " \f\n\t") const;

         Dstring
         get_left_trimmed (const Dstring& white_string = " \f\n\t") const;

         Dstring
         get_right_trimmed (const Dstring& white_string = " \f\n\t") const;

         Dstring
         get_file_extension () const;

         static Dstring
         render (const char* format, ...);

         static Dstring
         render (const Dstring& format, ...);

   };

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

         Reg_Exp (const Dstring& reg_exp_str,
                  const bool match_info = false,
                  const bool posix_extend = true,
                  const bool case_sensitive = true);

         ~Reg_Exp ();

         bool
         match (const Dstring& str) const;

         Iduple
         get_match (const Dstring& str,
                    const bool ignore_no_match = true) const;

         void
         replace (Dstring& str,
                  const Dstring& with,
                  const bool ignore_no_match = true) const;

         Iduple_Vector
         get_match_sub (const Dstring& str,
                        const bool ignore_no_match = true) const;

         static bool
         match (const Dstring& str,
                const Dstring& reg_exp_str,
                const bool posix_extend = true,
                const bool case_sensitive = true);

   };

   class Tokens : public vector<Dstring>
   {

      public:

         Tokens ();

         Tokens (const Tuple& tuple,
                 const Dstring& fmt);

         Tokens (const Dstring& str);

         Tokens (const Dstring& str,
                 const Dstring& delimiters);

         Tokens
         subtokens (const Integer i,
                    const Integer n = -1) const;

         Real
         real (const Integer index) const;

         Integer
         integer (const Integer index) const;

         void
         add (const Dstring& str);

         void
         add (const Dstring& str,
              const Dstring& delimiters);

         void
         add_prefix (const Dstring& prefix);

         void
         add_suffix (const Dstring& suffix);

   };

   ostream&
   operator << (ostream &out_file,
                const Tokens& tokens);

}


#endif /* DENISE_DSTRING_H */

