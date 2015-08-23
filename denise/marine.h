//
// marine.h
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

#ifndef DENISE_MARINE_H
#define DENISE_MARINE_H

#include <cmath>
#include <denise/basics.h>
#include <map>

using namespace std;

namespace denise
{

   class Gd_Sea_Wave_Fetch_Coefficients
   {

      public:

         const Real
         a_0;

         const Real
         a_1;

         const Real
         a_2;

         const Real
         a_3;

         const Real
         a_4;

         const Real
         b_0;

         const Real
         b_1;

         const Real
         b_2;

         const Real
         b_3;

         const Real
         b_4;

         Gd_Sea_Wave_Fetch_Coefficients ();

   };

   class Gd_Sea_Wave_Duration_Coefficients
   {

      public:

         const Real
         a_0;

         const Real
         a_1;

         const Real
         a_2;

         const Real
         b_0;

         const Real
         b_1;

         const Real
         b_2;

         const Real
         c_0;

         const Real
         c_1;

         const Real
         c_2;

         const Real
         c_3;

         const Real
         d_0;

         const Real
         d_1;

         const Real
         d_2;

         const Real
         d_3;

         const Real
         d_4;

         Gd_Sea_Wave_Duration_Coefficients ();

   };

   class Gd_Sea_Wave
   {

      private:

         const Gd_Sea_Wave_Fetch_Coefficients
         gd_sea_wave_fetch_coefficients;

         const Gd_Sea_Wave_Duration_Coefficients
         gd_sea_wave_duration_coefficients;

      public:

         Real
         get_seas_from_fetch (const Real speed,
                              const Real fetch) const;

         Real
         get_seas_from_duration (const Real speed,
                                 const Real duration) const;

   };

}

#endif /* DENISE_MARINE_H */

