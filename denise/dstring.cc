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

Dstring::Dstring ()
{
}

Dstring::Dstring (const Dstring& dstring)
   : wstring (dstring)
{
}

Dstring::Dstring (const wstring& wstr)
   : wstring (wstr)
{
}

Dstring::Dstring (const wchar_t* buffer)
   : wstring (buffer)
{
}

Dstring::Dstring (const string& str)
   : wstring (str.begin (), str.end ())
{
}

Dstring::Dstring (const char* buffer)
{
   const string str (buffer);
   this->assign (wstring (str.begin (), str.end ()));
}

string
Dstring::get_string () const
{
   return string (begin (), end ());
}

void
Dstring::to_lower_case ()
{
   transform (begin (), end (), begin (), ::tolower);
}
   
void
Dstring::to_upper_case ()
{
   transform (begin (), end (), begin (), ::toupper);
}

void
Dstring::to_capital_case ()
{

   bool first_char = true;

   for (auto iterator = begin (); iterator != end (); iterator++)
   {

      wchar_t& c = *(iterator);

      if (!first_char) { c = ::tolower (c); }
      else { c = ::toupper (c), first_char = false; }

      if (!isalnum (c)) { first_char = true; }

   }

}

void
Dstring::chop ()
{
   erase (size () - 1);
}

void
Dstring::trim (const Dstring& white_string)
{
   left_trim (white_string);
   right_trim (white_string);
}

void
Dstring::left_trim (const Dstring& white_string)
{
   erase (0, find_first_not_of (white_string));
}

void
Dstring::right_trim (const Dstring& white_string)
{
   erase (find_last_not_of (white_string) + 1); 
}

Dstring
Dstring::get_lower_case () const
{
   Dstring str (*this);
   str.to_lower_case ();
   return str;
}

Dstring
Dstring::get_upper_case () const
{
   Dstring str (*this);
   str.to_upper_case ();
   return str;
}

Dstring
Dstring::get_captial_case () const
{
   Dstring str (*this);
   str.to_capital_case ();
   return str;
}

Dstring
Dstring::get_trimmed (const Dstring& white_string) const
{
   Dstring trimmed (*this);
   trimmed.trim (white_string);
   return trimmed;
}

Dstring
Dstring::get_left_trimmed (const Dstring& white_string) const
{
   Dstring trimmed (*this);
   trimmed.left_trim (white_string);
   return trimmed;
}

Dstring
Dstring::get_right_trimmed (const Dstring& white_string) const
{
   Dstring trimmed (*this);
   trimmed.right_trim (white_string);
   return trimmed;
}

Dstring
Dstring::get_file_extension () const
{

   Dstring::size_type idx = rfind (L'.');

   if (idx != Dstring::npos)
   {
      return substr (idx + 1);
   } 
   else
   {
      return (*this);
   }

}

Dstring
Dstring::render (const Dstring& format, ...)
{

   va_list ap;

   int n = format.size () * 5;
   char* c_rendered = new char[n];

   va_start (ap, format);
   vsnprintf (c_rendered, n, format.get_string ().c_str (), ap);
   va_end (ap);

   const Dstring str (c_rendered);
   delete[] c_rendered;
   return str;

}

Reg_Exp::Reg_Exp (const Dstring& reg_exp_str,
                  const bool match_info,
                  const bool posix_extend,
                  const bool case_sensitive)
{

   int flag = (match_info ? 0 : REG_NOSUB);
   if (posix_extend) { flag |= REG_EXTENDED; }
   if (!case_sensitive) { flag |= REG_ICASE; }

   const string str (reg_exp_str.begin (), reg_exp_str.end ());
   regcomp (&preg, str.c_str (), flag);

   nmatch = preg.re_nsub + 1;
   pmatch = new regmatch_t[nmatch];

}

Reg_Exp::~Reg_Exp ()
{
   regfree (&preg);
   delete[] pmatch;
}

bool
Reg_Exp::match (const Dstring& str) const
{
   const string s (str.begin (), str.end ());
   return (regexec (&preg, s.c_str (), nmatch, pmatch, 0) == 0);
}

Iduple
Reg_Exp::get_match (const Dstring& str,
                    const bool ignore_no_match) const
{

   Iduple match;
   const string s (str.begin (), str.end ());

   if (regexec (&preg, s.c_str (), nmatch, pmatch, 0) == 0)
   {
      match.first = pmatch[0].rm_so;
      match.second = pmatch[0].rm_eo;
      return match;
   }
   else
   {
      if (!ignore_no_match)
      {
         throw Exception (L"Reg_Exp: no match ");
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
Reg_Exp::replace (Dstring& str,
                  const Dstring& with,
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
Reg_Exp::get_match_sub (const Dstring& str,
                        const bool ignore_no_match) const
{

   Iduple_Vector iduple_vector;
   const string s (str.begin (), str.end ());

   if (regexec (&preg, s.c_str (), nmatch, pmatch, 0) == 0)
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
         throw Exception (L"Reg_Exp: no match ");
      }
   }

   return iduple_vector;


}

bool
Reg_Exp::match (const Dstring& str,
                const Dstring& reg_exp_str,
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
                const Dstring& fmt)
{

   string s_fmt (fmt.begin (), fmt.end ());

   for (Tuple::const_iterator iterator = tuple.begin ();
        iterator != tuple.end (); iterator++)
   {
      const Real value = *(iterator);
      const Dstring& token = string_render (s_fmt.c_str (), value);
      push_back (token);  
   }
}

Tokens::Tokens (const Dstring& str)
{
   add (str);
}

Tokens::Tokens (const Dstring& str,
                const Dstring& delimiters)
{
   add (str, delimiters);
}

Tokens
Tokens::subtokens (const Integer i,
                   const Integer n) const
{

   Tokens subtokens;
   const Integer nn = ((n < 0) ? size () - i : n);

   for (Integer index = i; index < i + nn; index++)
   {
      subtokens.push_back (at (index));
   }

   return subtokens;

}

Real
Tokens::real (const Integer index) const
{
   return stof (at (index));
}

Integer
Tokens::integer (const Integer index) const
{
   return stoi (at (index));
}

void
Tokens::add (const Dstring& str)
{

   Dstring token;
   bool is_number = false;
   bool is_variable = false;
   bool is_equal = false;
   bool in_single_quote = false;
   bool in_double_quote = false;
   bool is_escaped = false;

   for (int i = 0; i < str.size (); i++)
   {

      wchar_t c = str[i];

      if (is_escaped)
      {
         token += c;
         is_escaped = false;
      }
      else
      if (isalpha (c) || c == '_')
      {
         is_variable = true;
         token += c;
      }
      else
      if (isdigit (c))
      {
         if (is_variable) { token += c; }
         else { is_number = true; token += c; }
      }
      else
      if (isspace (c))
      {
         if (in_double_quote) { token += c; }
         else
         {
            is_number = false;
            is_variable = false;
            is_equal = false;
            push_back (token);
            token.clear ();
         }
      }
      else
      if (c == 0x5c) // blackslash
      {
         is_escaped = true;;
      }
      else
      if (c == 0x27) // single quote
      {
         if (in_double_quote) { token += c; }
         else { in_single_quote = !in_single_quote; }
      }
      else
      if (c == 0x22) // double quote
      {
         if (in_single_quote) { token += c; }
         else { in_double_quote = !in_double_quote; }
      }
      else
      {
         if (in_double_quote) { token += c; }
         else
         {
            is_number = false;
            is_variable = false;
            is_equal = false;
            push_back (token);
            token.clear ();
            token += c;
            push_back (token);
            token.clear ();
         }
      }

   }

   push_back (token);

}

void
Tokens::add (const Dstring& str,
             const Dstring& delimiters)
{

   if (delimiters.empty ()) { add (str); }

   int n = str.size ();
   int nd = delimiters.size ();

   char* c_str = new char[n + 1];    c_str[n] = '\0';
   char* c_delim = new char[nd + 1]; c_delim[nd] = '\0';
   char* buffer = new char[n + 1];   buffer[n] = '\0';

   for (Integer i = 0; i < n; i++) { c_str[i] = (char)str[i]; }
   for (Integer i = 0; i < nd; i++) { c_delim[i] = (char)delimiters[i]; }

   char* c_token = strtok_r (c_str, c_delim, (char**)buffer);

   if (c_token != NULL)
   {

      const Dstring token (c_token);
      push_back (token);

      while ((c_token = strtok_r (NULL, c_delim, (char**)buffer)) != NULL)
      {
         const Dstring token (c_token);
         push_back (token);
      }

   }

   delete[] buffer;

}

void
Tokens::add_prefix (const Dstring& prefix)
{
   for (auto iterator = begin (); iterator != end (); iterator++)
   {
      Dstring& token = *(iterator);
      token.assign (prefix + token);
   }
}

void
Tokens::add_suffix (const Dstring& suffix)
{
   for (auto iterator = begin (); iterator != end (); iterator++)
   {
      Dstring& token = *(iterator);
      token.assign (token + suffix);
   }
}

namespace denise
{

   void
   to_lower_case (wstring& s)
   {
      transform (s.begin (), s.end (), s.begin (), ::tolower);
   }

   wstring
   get_lower_case (const wstring& s)
   {
      wstring str = s;
      to_lower_case (str);
      return str;
   }

   void
   to_upper_case (wstring& s)
   {
      transform (s.begin (), s.end (), s.begin (), ::toupper);
   }

   wstring
   get_upper_case (const wstring& s)
   {
      wstring str = s;
      to_upper_case (str);
      return str;
   }

   void
   to_capital_case (wstring& s)
   {

      bool first_char = true;

      for (wstring::iterator iterator = s.begin ();
           iterator != s.end (); iterator++)
      {

         wchar_t& c = *(iterator);

         if (!first_char) { c = ::tolower (c); }
         else { c = ::toupper (c), first_char = false; }

         if (!isalnum (c)) { first_char = true; }

      }

   }

   wstring
   get_capital_case (const wstring& s)
   {
      wstring str = s;
      to_capital_case (str);
      return str;
   }

   void
   chop (wstring& s)
   {
      s.erase (s.size () - 1);
   }

   void
   trim (wstring& s,
         const wstring& white_string)
   {
      left_trim (s, white_string);
      right_trim (s, white_string);
   }

   void
   left_trim (wstring& s,
              const wstring& white_string)
   {
      s.erase (0, s.find_first_not_of (white_string));
   }
   
   void
   right_trim (wstring& s,
               const wstring& white_string)
   {
      s.erase (s.find_last_not_of (white_string) + 1); 
   }

   wstring
   get_trimmed (const wstring s,
                const wstring& white_string)
   {
      wstring trimmed = s;
      trim (trimmed);
      return trimmed;
   }

   wstring
   get_left_trimmed (const wstring s,
                     const wstring& white_string)
   {
      wstring trimmed = s;
      left_trim (trimmed);
      return trimmed;
   }

   wstring
   get_right_trimmed (const wstring s,
                      const wstring& white_string)
   {
      wstring trimmed = s;
      right_trim (trimmed);
      return trimmed;
   }

   wstring
   get_file_extension (const wstring& file_path)
   {

      wstring::size_type idx = file_path.rfind ('.');

      if (idx != wstring::npos)
      {
         return file_path.substr (idx + 1);
      } 
      else
      {
         return wstring ();
      }

   }

   wstring
   string_render (const char* format,
                  ...)
   {

      va_list ap;

      int n = strlen (format) * 4;
      char* c_string = new char[n];

      va_start (ap, format);
      vsnprintf (c_string, n, format, ap);
      va_end(ap);

      const string str (c_string);
      delete[] c_string;
      wstring w_str (str.begin (), str.end ());
      return w_str;

   }

}

namespace denise
{

   wostream&
   operator << (wostream& out_file,
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
