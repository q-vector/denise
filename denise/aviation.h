//
// aviation.h
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

#ifndef DENISE_AVIATION_H
#define DENISE_AVIATION_H

#include <cmath>
#include <denise/basics.h>
#include <map>

using namespace std;

namespace denise
{

   class Standard_Atmosphere
   {

       private:

          const Real
          p_0;

          const Real
          rho_0;

          const Real
          T_0;

          const Real
          g_0;

          const Real
          R_0;

          map<Real, Real>
          p_z;

          map<Real, Real>
          z_p;

          void
          init_p_z ();

          void
          init_z_p ();

       public:

          Standard_Atmosphere (const Real p_0 = 1.01325e5,
                               const Real rho_0 = 1.2250,
                               const Real T_0 = 288.16,
                               const Real g_0 = 9.80665,
                               const Real R_0 = 287.0368);

          Real
          get_temperature_gradient (const Real height) const;

          Real
          get_height (const Real pressure) const;

          Real
          get_pressure (const Real height) const;

   };

}

#endif /* DENISE_AVIATION_H */

