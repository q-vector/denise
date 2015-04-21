//
// exception.cc
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

#include "exception.h"

using namespace std;
using namespace denise;

Exception::Exception (const string& description) throw ()
   : identifier ("Exception"),
     description (description)
{
}

Exception::Exception (const string& identifier,
                      const string& description) throw ()
   : identifier (identifier),
     description (description)
{
}

Exception::Exception (const Exception& exception) throw ()
   : identifier (exception.identifier),
     description (exception.description)
{
}

Exception::~Exception () throw ()
{
}

const char*
Exception::what () const throw ()
{
   return (identifier + ": " + description).c_str ();
}

IO_Exception::IO_Exception (const string& description) throw ()
   : Exception ("IO_Exception", description)
{
}

No_Match_Exception::No_Match_Exception (const string& description) throw ()
   : Exception ("No_Match_Exception", description)
{
}

Out_Of_Bounds_Exception::Out_Of_Bounds_Exception (const string& description) throw ()
   : Exception ("Out_Of_Bounds_Exception", description)
{
}

namespace denise
{

   ostream&
   operator << (ostream &out_file,
                std::exception& exception)
   {
      out_file << exception.what ();
      return out_file;
   }

   ostream&
   operator << (ostream &out_file,
                const std::exception& exception)
   {
      out_file << exception.what ();
      return out_file;
   }

}

