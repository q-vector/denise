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
// gzstream is derived from the Gzstream Library (LGPL license) written
// by Deepak Bandyopadhyay and Lutz Kettner at the Computational Geometry
// Group at UNC Chapel Hill http://www.cs.unc.edu/Research/compgeom/gzstream/

#include <cstring>
#include "exception.h"
#include "gzstream.h"

using namespace std;
using namespace denise;

gzstreambuf::gzstreambuf ()
   : opened (false)
{

   char* read_position = buffer + 4;
   char* end_position = buffer + 4;
   char* beginning_of_putback_area = buffer + 4;

   setp (buffer, buffer + (buffer_size - 1));
   setg (beginning_of_putback_area, read_position, end_position);

}

gzstreambuf::~gzstreambuf ()
{
   close ();
}
       
bool
gzstreambuf::is_open ()
{
   return opened;
}

gzstreambuf*
gzstreambuf::open (const string& file_path,
                   int open_mode)
{

   if (is_open ()) { return (gzstreambuf*)0; }
   this->mode = open_mode;

   const bool is_atend = (mode & std::ios::ate);
   const bool is_append = (mode & std::ios::app);
   const bool is_r = (mode & std::ios::in);
   const bool is_w = (mode & std::ios::out);

   // no atend / append / read-write mode
   if (is_atend) { return (gzstreambuf*)0; } 
   if (is_append) { return (gzstreambuf*)0; } 
   if (is_r && is_w) { return (gzstreambuf*)0; } 

   gz_file = gzopen (file_path.c_str (), is_r ? "rb" : "wb");
   if (gz_file == 0) { return (gzstreambuf*)0; }

   this->opened = true;
   return this;

}

gzstreambuf*
gzstreambuf::close ()
{

   if (is_open ())
   {
      sync ();
      opened = false;
      if (gzclose (gz_file) == Z_OK) { return this; }
   }

   return (gzstreambuf*)0;

}

int
gzstreambuf::underflow ()
{
 
   // used for input buffer only
   if (gptr () && (gptr () < egptr ()))
   {
      return * reinterpret_cast<unsigned char *>(gptr ());
   }

   if (!(mode & std::ios::in) || !opened) { return EOF; }

   const int putback_size = std::min (4, int (gptr () - eback ()));

   char* read_position = buffer + 4;
   char* beginning_of_putback_area = read_position - putback_size;

   memcpy (beginning_of_putback_area, gptr () - putback_size, putback_size);
   const int num = gzread (gz_file, buffer + 4, buffer_size - 4);

   // ERROR or EOF
   if (num <= 0) { return EOF; } 

   // reset buffer pointers
   char* end_position = read_position + num;
   setg (beginning_of_putback_area, read_position, end_position);

   // return next character
   return *reinterpret_cast<unsigned char *>(gptr ());    

}

int
gzstreambuf::flush_buffer ()
{
   // Separate the writing of the buffer from overflow () and
   // sync () operation.
   int w = pptr () - pbase ();
   if (gzwrite (gz_file, pbase (), w) != w) { return EOF; }
   pbump (-w);
   return w;
}

int
gzstreambuf::overflow (int c)
{

   if (!(mode & std::ios::out) || ! opened)
   {
      return EOF;
   }

   if (c != EOF)
   {
      *pptr () = c;
      pbump (1);
   }

   if (flush_buffer () == EOF)
   {
      return EOF;
   }

   return c;

}

int
gzstreambuf::sync ()
{

   if (pptr () && pptr () > pbase ())
   {
      if (flush_buffer () == EOF) { return -1; }
   }

   return 0;

}

gzstreambase::gzstreambase ()
{
   init (&buf);
}

gzstreambase::gzstreambase (const string& file_path,
                            int mode)
{
   init (&buf);
   open (file_path, mode);
}

gzstreambase::~gzstreambase ()
{
   buf.close ();
}

bool
gzstreambase::is_open ()
{
   return buf.is_open ();
}

void
gzstreambase::open (const string& file_path,
                    int open_mode)
{
   if (!buf.open (file_path, open_mode))
   {
      clear (rdstate () | std::ios::badbit);
   }
}

void
gzstreambase::close ()
{
   if (buf.is_open () && !buf.close ())
   {
      clear (rdstate () | std::ios::badbit);
   }
}

gzstreambuf*
gzstreambase::rdbuf ()
{
   return &buf;
}

igzstream::igzstream ()
   : std::istream (&buf)
{
}

igzstream::igzstream (const string& file_path,
                      int open_mode)
   : gzstreambase (file_path, open_mode),
     std::istream (&buf)
{
   if (!is_open ()) { throw IO_Exception ("Can't open " + file_path); }
}

gzstreambuf*
igzstream::rdbuf ()
{
   return gzstreambase::rdbuf ();
}

void
igzstream::open (const string& file_path,
                 int open_mode)
{
   gzstreambase::open (file_path, open_mode);
   if (!is_open ()) { throw IO_Exception ("Can't open " + file_path); }
}

ogzstream::ogzstream ()
   : std::ostream (&buf)
{
}

ogzstream::ogzstream (const string& file_path,
                      int mode)
   : gzstreambase (file_path, mode),
     std::ostream (&buf)
{
   if (!is_open ()) { throw IO_Exception ("Can't open " + file_path); }
}

gzstreambuf*
ogzstream::rdbuf ()
{
   return gzstreambase::rdbuf ();
}

void
ogzstream::open (const string& file_path,
                 int open_mode)
{
   gzstreambase::open (file_path, open_mode);
   if (!is_open ()) { throw IO_Exception ("Can't open " + file_path); }
}

