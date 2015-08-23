//
// exception.h
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

#ifndef DENISE_EXCEPTION_H
#define DENISE_EXCEPTION_H

#include <iostream>
#include <denise/dstring.h>

using namespace std;

namespace denise
{

   class Exception : public std::exception
   {

      public:

         string
         identifier;

         string
         description;

         Exception (const string& description = string ()) throw ();

         Exception (const string& identifier,
                    const string& description) throw ();

         Exception (const Exception& exception) throw ();

         ~Exception () throw ();

         const char*
         what () const throw ();

   };

   class IO_Exception : public Exception
   {

      public:

         IO_Exception (const string& description = string ()) throw ();

   };

   class No_Match_Exception : public Exception
   {

      public:

         No_Match_Exception (const string& description = string ()) throw ();

   };

   class Out_Of_Bounds_Exception : public Exception
   {

      public:

         Out_Of_Bounds_Exception (const string& description = string ()) throw ();

   };

   ostream&
   operator << (ostream &out_file,
                std::exception& exception);

   ostream&
   operator << (ostream &out_file,
                const std::exception& exception);

}

#endif /* DENISE_EXCEPTION_H */

