//
// util.cc
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

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <map>
#include <stdint.h>
#include "util.h"

using namespace std;
using namespace denise;

namespace denise
{

#ifndef HAVE_EXP10
   Real
   exp10 (const Real x)
   {
      return pow (10, x);
   }
#endif

   Real
   boolean_sign (const bool b)
   {
      return (b ? 1 : -1);
   }

   bool
   is_even (const Integer i)
   {
      return (i % 2 == 0);
   }


   Real
   modulo (const Real x,
           const Real y)
   {
      Real m = fmod (x, y);
      return (m > 0 ? m : m + y);
   }

   Real
   modulo (const Real x,
           const Real start_x,
           const Real end_x)
   {
      Real span_x = end_x - start_x;
      Real m = denise::modulo ((x - start_x), span_x);
      return m + start_x;
   }

   Real
   modulo (const Real x,
           const Domain_1D domain_x)
   {
      return modulo (x, domain_x.start, domain_x.end);
   }

   Integer
   imodulo (const Integer i,
            const Integer j)
   {
      Integer m = i % j;
      return (m >= 0 ? m : m + j);
   }

   Real
   random (const Real max,
           const Real min)
   {
      Real r = static_cast<Real> (rand ()) / RAND_MAX;
      return r * (max - min) + min;
   }

   Integer
   irandom (const Integer max,
            const Integer min)
   {
      Integer m = max - min + 1;
      return (rand () % m) + min; 
   }

   Real
   determinant (const Real a00,
                const Real a01,
                const Real a10,
                const Real a11)
   {
      return a00*a11 - a01*a10;
   }

   Real
   determinant (const Real a00,
                const Real a01,
                const Real a02,
                const Real a10,
                const Real a11,
                const Real a12,
                const Real a20,
                const Real a21,
                const Real a22)
   {
      return a00*a11*a22 + a10*a21*a02 + a20*a01*a12 -
             a00*a21*a12 - a10*a01*a22 - a20*a11*a02;
   }

   double
   max (double a,
        double b)
   {
      return (a > b ? a : b);
   }

   double
   min (double a,
        double b)
   {
      return (a < b ? a : b);
   }

   FILE*
   get_input_file (const string& file_path)
            throw (IO_Exception)
   {
      FILE* file = fopen (file_path.c_str (), "r");
      if (file == NULL) { throw IO_Exception (file_path); }
      return file;
   }

   FILE*
   get_output_file (const string& file_path,
                    const char* mode)
             throw (IO_Exception)
   {

      FILE* file = NULL;

           if (file_path.empty ()) { file = NULL; }
      else if (file_path == "-") { file = stdout; }
      else
      {
         file = fopen (file_path.c_str (), mode);
         if (file == NULL) { throw IO_Exception (file_path); }
      }

      return file;

   }

   Tokens
   get_dir_listing (const string& dir_path,
                    const Reg_Exp& reg_exp,
                    const bool prepend_dir_path)
             throw (IO_Exception)
   {

      DIR* dp;
      struct dirent* dirp;

      Tokens dir_listing;
      dp = opendir (dir_path.c_str ());
      if (dp == NULL) { throw IO_Exception (dir_path); }

      while ((dirp = readdir (dp)) != NULL)
      {
         string file_name (dirp->d_name);
         if (reg_exp.match (file_name))
         {
            if (prepend_dir_path)
            {
               dir_listing.push_back (dir_path + "/" + file_name);
            }
            else
            {
               dir_listing.push_back (file_name);
            }
         }
      }

      closedir (dp);

      return dir_listing;

   }

   string
   read_line (FILE* file,
              int max_length,
              bool chop)
       throw (IO_Exception)
   {

      char* input_line = new char[max_length];
      char* s = fgets (input_line, max_length, file);

      if (s == 0) { throw IO_Exception ("fgets error or eof"); }
      string str (input_line);
      if (chop) { denise::chop (str); }

      delete[] input_line;
      return str;

   }

   gzFile
   get_gzfile (const string& file_path)
        throw (IO_Exception)
   {

      gzFile file = gzopen (file_path.c_str (), "rb");

      if (file == NULL)
      {

         file = gzopen ((file_path + ".gz").c_str (), "rb");

         if (file == NULL)
         {

            file = gzopen ((file_path + ".Z").c_str (), "rb");

            if (file == NULL)
            {
               throw IO_Exception (file_path);
            }

         }

      }

      return file;

   }

   char*
   gz_readline (char* input_line,
                int input_line_length,
                gzFile file)
   {


      int i = 0;
      char c = gzgetc (file);
      if (gzeof (file)) { return NULL; }

      while (!gzeof (file) && i < input_line_length - 2)
      {

         if (c == 0x0a || c == 0x0d)
         {
            break;
         }

         input_line[i] = c;
         i++;
         c = gzgetc (file);

      }

      input_line[i] = '\0';

      return input_line;

   }

   void
   swap_endian (void* data,
                short size)
   {

      uint8_t* d = (uint8_t*)data;

      switch (size) {
         case 2:
            swap (*(d+1), *(d));
            break;
         case 4:
            swap (*(d+3), *(d)); swap (*(d+2), *(d+1));
            break;
         case 8:
            swap (*(d+7), *(d));   swap (*(d+6), *(d+1));
            swap (*(d+5), *(d+2)); swap (*(d+4), *(d+3));
            break;
         case 16:
            swap (*(d+15), *(d));   swap (*(d+14), *(d+1));
            swap (*(d+13), *(d+2)); swap (*(d+12), *(d+3));
            swap (*(d+11), *(d+4)); swap (*(d+10), *(d+5));
            swap (*(d+9),  *(d+6)); swap (*(d+8),  *(d+7));
            break;
      }

   }

   double
   get_delta (double start_value,
              double end_value,
              int approx_number)
   {
      double mantissa = (end_value - start_value) / approx_number;
      double exponent = floor (log10 (mantissa));
      return rint (mantissa * pow (10, -exponent)) * pow (10, exponent);
   }

   double*
   gen_array (int size_1d,
              double start_value,
              double end_value)
   {

      double* values = new double[size_1d];
      double delta_value = (end_value - start_value) / (size_1d - 1);

      for (int i = 0; i < size_1d; i++) {
         values[i] = start_value + i * delta_value;
      }

      return values;

   }

   unsigned long
   get_file_size (const string& file_path)
           throw (IO_Exception)
   {

      struct stat buffer;

      if (stat (file_path.c_str (), &buffer) != 0)
      {
         throw IO_Exception (file_path);
      }

      return buffer.st_size;

   }

   Dtime
   get_file_last_modify_time (const string& file_path)
           throw (IO_Exception)
   {

      struct stat buffer;

      if (stat (file_path.c_str (), &buffer) != 0)
      {
         throw IO_Exception (file_path);
      }

      return Dtime (Real (buffer.st_mtime) / 3600);

   }

   void
   print_mem_usage (const Integer indent)
   {
      pid_t pid = getpid ();
      string s; for (Integer i = 0; i < indent; i++) { s += " "; }
      const string& str = string_render (
         "echo -n \"%s\"; pmap -x %i | tail -1", s.c_str (), pid);
      system (str.c_str ());
   }

}

bool
Assignable::apply_variables (string& line) const
{

   size_t found;
   for (auto iterator = assign_dictionary.begin ();
        iterator != assign_dictionary.end (); iterator++)
   {
      const string variable ("$" + iterator->first);
      while ((found = line.find (variable)) != string::npos)
      {
         const size_t var_length = variable.length ();
         line.replace (found, var_length, iterator->second);
      }
   }

   return (line.find ("$") != string::npos);

}

void
Config_File::apply_variables ()
{
   for (auto iterator = begin (); iterator != end (); iterator++)
   {
      auto& line = *(iterator);
      while (Assignable::apply_variables (line));
   }
}

Config_File::Config_File (const string& file_path)
{
   ingest (file_path);
}

void
Config_File::ingest (const string& file_path)
{

   ifstream file (file_path);

   for (string input_string; getline (file, input_string); )
   {

      if (input_string.size () == 0) { continue; }
      if (input_string.c_str ()[0] == '#') { continue; }

      const Tokens variable_tokens (input_string, "=");
      if (variable_tokens.size () == 2)
      {
         const string& variable = get_trimmed (variable_tokens[0]);
         const string& value = get_trimmed (variable_tokens[1]);
         assign_dictionary.insert (make_pair (variable, value));
         continue;
      }

      push_back (input_string);

   }

   file.close ();

   apply_variables ();

}

