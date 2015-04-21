//
// util.h
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

#ifndef DENISE_UTIL_H
#define DENISE_UTIL_H

#include <ctype.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <zlib.h>
#include <denise/basics.h>
#include <denise/dstring.h>
#include <denise/dtime.h>
#include <denise/exception.h>

using namespace std;

namespace denise
{

#ifndef HAVE_EXP10
   Real
   exp10 (const Real x);
#endif

   Real
   boolean_sign (const bool b);

   bool
   is_even (const Integer i);

   Real
   modulo (const Real x,
           const Real y);

   Real
   modulo (const Real x,
           const Real start_x,
           const Real end_x);

   Real
   modulo (const Real x,
           const Domain_1D domain_x);

   Integer
   imodulo (const Integer i,
            const Integer j);

   /// Returns a random number between max and min.
   ///
   /// Returns a random number between max and min.
   /// 
   /// \param max  upper bound of random number
   /// \param min  lower bound of random number
   Real
   random (const Real max = 1,
           const Real min = 0);

   /// Returns a random number between max and min.
   ///
   /// Returns a random number between max and min.
   /// 
   /// \param max  upper bound of random number
   /// \param min  lower bound of random number
   Integer
   irandom (const Integer max = 1,
            const Integer min = 0);

   Real
   determinant (const Real a00,
                const Real a01,
                const Real a10,
                const Real a11);

   Real
   determinant (const Real a00,
                const Real a01,
                const Real a02,
                const Real a10,
                const Real a11,
                const Real a12,
                const Real a20,
                const Real a21,
                const Real a22);

   double
   max (double a,
        double b);

   double
   min (double a,
        double b);

   template <class T> T
   bound (const T& t,
          const T& max = 1,
          const T& min = 0)
   {
      if (t > max) { return max; }
      if (t < min) { return min; }
      return t;
   }

   FILE*
   get_input_file (const string& file_path)
            throw (IO_Exception);

   FILE*
   get_output_file (const string& file_path,
                    const char* mode = "w")
             throw (IO_Exception);

   Tokens
   get_dir_listing (const string& dir_path,
                    const Reg_Exp& reg_exp,
                    const bool prepend_dir_path = false)
             throw (IO_Exception);

   string
   read_line (FILE* file,
              int max_length = 2048,
              bool chop = true)
       throw (IO_Exception);

   gzFile
   get_gzfile (const string& file_path)
        throw (IO_Exception);

   char*
   gz_readline (char* input_line,
                int input_line_length,
                gzFile file);

   void
   swap_endian (void* data,
                short size);

   double
   get_delta (double start_value,
              double end_value,
              int approx_number);

   double*
   gen_array (int size_1d,
              double start_value,
              double end_value);

   unsigned long
   get_file_size (const string& file_path)
      throw (IO_Exception);

   Dtime
   get_file_last_modify_time (const string& file_path)
      throw (IO_Exception);

   template <class T> void
   shuffle_vector (vector<T>& v)
   {

      int n = v.size ();

      for (int i = 0; i < v.size (); i++)
      {
         int ii = int (rint (random (double (n))));
         swap (v[i], v[ii]);
      }

   }

   template <class T> void
   shuffle_array (T* array,
                  int n)
   {

      T temp;

      for (int i = 0; i < n; i++)
      {
         int ii = int (rint (random (double (n))));
         temp = array[i];
         array[i] = array[ii];
         array[ii] = temp;
      }

   }

   void
   print_mem_usage (const Integer indent = 0);

}

#endif /* DENISE_UTIL_H */

