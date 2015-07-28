//
// gzstream.h
// 
// Copyright (C) 2005-2015 Simon E. Ching
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
//
// gzstream is derived from the Gzstream Library (LGPL License) written
// by Deepak Bandyopadhyay and Lutz Kettner at the Computational Geometry
// Group at UNC Chapel Hill http://www.cs.unc.edu/Research/compgeom/gzstream/

#ifndef DENISE_GZSTREAM_H
#define DENISE_GZSTREAM_H

#include <fstream>
#include <iostream>
#include <string>
#include <zlib.h>

using namespace std;

namespace denise
{

   class gzstreambuf : public std::streambuf
   {

      private:

         // totals 512 bytes under g++ for igzstream at the end.
         // size of data buff
         static const int
         buffer_size = 47 + 256;

         gzFile
         gz_file;

         char
         buffer[buffer_size];

         bool
         opened;

         int
         mode;

         int
         flush_buffer ();

      public:

         gzstreambuf ();

         bool
         is_open ();

         gzstreambuf*
         open (const string& file_path,
               int open_mode);

         gzstreambuf*
         close ();

         ~gzstreambuf ();

         virtual int
         overflow (int c = EOF);

         virtual int
         underflow ();

         virtual int
         sync ();

   };

   class gzstreambase : virtual public std::ios
   {

      protected:

         gzstreambuf
         buf;

      public:

         gzstreambase ();

         gzstreambase (const string& file_path,
                       int open_mode);

         ~gzstreambase ();

         void
         open (const string& file_path,
               int open_mode);

         void
         close ();

         gzstreambuf*
         rdbuf ();

   };

   class igzstream : public gzstreambase,
                     public std::istream
   {

      public:

         igzstream ();

         igzstream (const string& file_path,
                    int open_mode = std::ios::in);

         gzstreambuf*
         rdbuf ();

         void
         open (const string& file_path,
               int open_mode = std::ios::in);

   };

   class ogzstream : public gzstreambase,
                     public std::ostream
   {

      public:

         ogzstream ();

         ogzstream (const string& file_path,
                    int mode = std::ios::out);

         gzstreambuf*
         rdbuf ();

         void
         open (const string& file_path,
               int open_mode = std::ios::out);

   };

};

#endif /* DENISE_GZSTREAM_H */

