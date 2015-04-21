//
// dstring.cc
// 
// Copyright (C) 2010 Simon E. Ching
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

#include "dstring.h"

using namespace std;
using namespace denise;

Reg_Exp::Reg_Exp (const string& reg_exp_str,
                  const bool match_info,
                  const bool posix_extend,
                  const bool case_sensitive)
{

   int flag = (match_info ? 0 : REG_NOSUB);
   if (posix_extend) { flag |= REG_EXTENDED; }
   if (!case_sensitive) { flag |= REG_ICASE; }

   regcomp (&preg, reg_exp_str.c_str (), flag);

   nmatch = preg.re_nsub + 1;
   pmatch = new regmatch_t[nmatch];

}

Reg_Exp::~Reg_Exp ()
{
   regfree (&preg);
   delete[] pmatch;
}

bool
Reg_Exp::match (const string& str) const
{
   return (regexec (&preg, str.c_str (), nmatch, pmatch, 0) == 0);
}

Iduple
Reg_Exp::get_match (const string& str,
                    const bool ignore_no_match) const
{

   Iduple match;

   if (regexec (&preg, str.c_str (), nmatch, pmatch, 0) == 0)
   {
      match.first = pmatch[0].rm_so;
      match.second = pmatch[0].rm_eo;
      return match;
   }
   else
   {
      if (!ignore_no_match)
      {
         throw Exception ("Reg_Exp: no match ");
      }
      else
      {
         match.first = -1;
         match.second = -1;
         return match;
      }
   }

   return match;

}

void
Reg_Exp::replace (string& str,
                  const string& with,
                  const bool ignore_no_match) const
{
   try
   {
      const Iduple& match = get_match (str, false);
      str.replace (match.first, match.second - match.first, with);
   }
   catch (const Exception& e)
   {
      if (!ignore_no_match) { throw e; }
   }
}

Iduple_Vector
Reg_Exp::get_match_sub (const string& str,
                        const bool ignore_no_match) const
{

   Iduple_Vector iduple_vector;

   if (regexec (&preg, str.c_str (), nmatch, pmatch, 0) == 0)
   {

      for (Integer i = 0; i < nmatch; i++)
      {

         Integer start_index = pmatch[i].rm_so;
         Integer end_index = pmatch[i].rm_eo;

         Iduple iduple = make_pair (start_index, end_index);
         iduple_vector.push_back (iduple);

      }

   }
   else
   {
      if (!ignore_no_match)
      {
         throw Exception ("Reg_Exp: no match ");
      }
   }

   return iduple_vector;


}

bool
Reg_Exp::match (const string& str,
                const string& reg_exp_str,
                const bool posix_extend,
                const bool case_sensitive)
{
   Reg_Exp reg_exp (reg_exp_str, true, posix_extend, case_sensitive);
   return reg_exp.match (str);
}

Tokens::Tokens ()
{
}

Tokens::Tokens (const Tuple& tuple,
                const string& fmt)
{
   for (Tuple::const_iterator iterator = tuple.begin ();
        iterator != tuple.end (); iterator++)
   {
      const Real value = *(iterator);
      const string& token = string_render (fmt.c_str (), value);
      push_back (token);  
   }
}

Tokens::Tokens (const string& str,
                const string& delimiters)
{
   add (str, delimiters);
}

Real
Tokens::real (const Integer index) const
{
   return atof (at (index).c_str ());
}

Integer
Tokens::integer (const Integer index) const
{
   return atoi (at (index).c_str ());
}

void
Tokens::add (const string& str,
             const string& delimiters)
{

   const char* c_string = str.c_str ();
   const char* c_delimiters = delimiters.c_str ();
   int n = strlen (c_string) + 1;

   char* copy = new char[n];
   char* buffer = new char[n];
   strncpy (copy, c_string, n);

   char* token = strtok_r (copy, c_delimiters, (char**)buffer);

   if (token != NULL)
   {

      push_back (string (token));

      while ((token = strtok_r (NULL, c_delimiters, (char**)buffer)) != NULL)
      {
         push_back (string (token));
      }

   }

   delete[] copy;
   delete[] buffer;

}

void
Tokens::add_prefix (const string& prefix)
{
   for (Tokens::iterator iterator = begin (); iterator != end (); iterator++)
   {
      string& token = *(iterator);
      token.assign (prefix + token);
   }
}

void
Tokens::add_suffix (const string& suffix)
{
   for (Tokens::iterator iterator = begin (); iterator != end (); iterator++)
   {
      string& token = *(iterator);
      token.assign (token + suffix);
   }
}

namespace denise
{

   vector<string>
   tokenize (const string& s,
             const string& delimiters)
   {

      vector<string> token_vector;

      const char* c_string = s.c_str ();
      const char* c_delimiters = delimiters.c_str ();
      int n = strlen (c_string) + 1;
 
      char* copy = new char[n];
      char* buffer = new char[n];
      strncpy (copy, c_string, n);

      char* token = strtok_r (copy, c_delimiters, (char**)buffer);

      if (token != NULL)
      {

         token_vector.push_back (string (token));

         while ((token = strtok_r (NULL, c_delimiters, (char**)buffer)) != NULL)
         {
            token_vector.push_back (string (token));
         }

      }

      delete[] copy;
      delete[] buffer;

      return token_vector;

   }

   vector<string>
   tokenize_s (const string& s,
               const string& delimiters)
   {

      vector<string> token_vector;

      const char* c_string = s.c_str ();
      const char* c_delimiters = delimiters.c_str ();
      int n = strlen (c_string) + 1;
 
      char* copy = new char[n];

      strncpy (copy, c_string, n);
      char* token = strsep (&copy, c_delimiters);

      if (token != NULL)
      {

         token_vector.push_back (string (token));

         while ((token = strsep (&copy, c_delimiters)) != NULL)
         {
            token_vector.push_back (string (token));
         }

      }

      delete[] copy;

      return token_vector;

   }

   void
   to_lower_case (string& s)
   {
      transform (s.begin (), s.end (), s.begin (), ::tolower);
   }

   string
   get_lower_case (const string& s)
   {
      string str = s;
      to_lower_case (str);
      return str;
   }

   void
   to_upper_case (string& s)
   {
      transform (s.begin (), s.end (), s.begin (), ::toupper);
   }

   string
   get_upper_case (const string& s)
   {
      string str = s;
      to_upper_case (str);
      return str;
   }

   void
   to_capital_case (string& s)
   {

      bool first_char = true;

      for (string::iterator iterator = s.begin ();
           iterator != s.end (); iterator++)
      {

         char& c = *(iterator);

         if (!first_char) { c = ::tolower (c); }
         else { c = ::toupper (c), first_char = false; }

         if (!isalnum (c)) { first_char = true; }

      }

   }

   string
   get_capital_case (const string& s)
   {
      string str = s;
      to_capital_case (str);
      return str;
   }

   void
   chop (string& s)
   {
      s.erase (s.size () - 1);
   }

   void
   trim (string& s,
         const string& white_string)
   {
      left_trim (s, white_string);
      right_trim (s, white_string);
   }

   void
   left_trim (string& s,
              const string& white_string)
   {
      s.erase (0, s.find_first_not_of (white_string));
   }
   
   void
   right_trim (string& s,
               const string& white_string)
   {
      s.erase (s.find_last_not_of (white_string) + 1); 
   }

   string
   get_file_extension (const string& file_path)
   {

      string::size_type idx = file_path.rfind('.');

      if (idx != string::npos)
      {
         return file_path.substr (idx+1);
      } 
      else
      {
         return "";
      }

   }

   string
   string_render (const char* format,
                  ...)
   {

      va_list ap;

      int n = strlen (format) * 4;
      char* c_string = new char[n];

      va_start (ap, format);
      vsnprintf (c_string, n, format, ap);
      va_end(ap);

      string str (c_string);
      delete[] c_string;
      return str;

   }

}

namespace denise
{

   ostream&
   operator << (ostream& out_file,
                const Tokens& tokens)
   {

      out_file << "(";

      for (Integer i = 0; i < tokens.size (); i++)
      {

         out_file << tokens[i];

         if (i != tokens.size () - 1)
         {
            out_file << ", ";
         }

      }

      out_file << ")";
      return out_file;

   }

}
