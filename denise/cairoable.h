//
// cairoable.h
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

#ifndef DENISE_CAIRO_H
#define DENISE_CAIRO_H

#include <cairomm/context.h>
#include <denise/transform.h>

using namespace std;
using namespace Cairo;

namespace denise
{

   class Cairoable
   {

      public:

         virtual void
         cairo (const RefPtr<Context>& cr) const = 0;

         virtual void
         cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform) const = 0;

   };

}

#endif /* DENISE_CAIROABLE_H */

