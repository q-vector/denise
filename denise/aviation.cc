//
// aviation.cc
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

#include "aviation.h"

using namespace std;
using namespace denise;

void
Standard_Atmosphere::init_p_z ()
{

   p_z.insert (make_pair (1000, 31055));
   p_z.insert (make_pair (2000, 26481));
   p_z.insert (make_pair (3000, 23849));
   p_z.insert (make_pair (4000, 22000));
   p_z.insert (make_pair (5000, 20576));
   p_z.insert (make_pair (7000, 18442));
   p_z.insert (make_pair (10000, 16180));
   p_z.insert (make_pair (15000, 13608));
   p_z.insert (make_pair (20000, 11784));
   p_z.insert (make_pair (22600, 11000));
   p_z.insert (make_pair (25000, 10363));
   p_z.insert (make_pair (30000, 9164));
   p_z.insert (make_pair (40000, 7185));
   p_z.insert (make_pair (50000, 5574));
   p_z.insert (make_pair (60000, 4206));
   p_z.insert (make_pair (70000, 3012));
   p_z.insert (make_pair (80000, 1949));
   p_z.insert (make_pair (85000, 1457));
   p_z.insert (make_pair (90000, 988));
   p_z.insert (make_pair (95000, 540));
   p_z.insert (make_pair (100000, 111));
   p_z.insert (make_pair (101325, 0));
   p_z.insert (make_pair (105000, -302));

}

void
Standard_Atmosphere::init_z_p ()
{

   z_p.insert (make_pair (31055, 1000));
   z_p.insert (make_pair (26481, 2000));
   z_p.insert (make_pair (23849, 3000));
   z_p.insert (make_pair (22000, 4000));
   z_p.insert (make_pair (20576, 5000));
   z_p.insert (make_pair (18442, 7000));
   z_p.insert (make_pair (16180, 10000));
   z_p.insert (make_pair (13608, 15000));
   z_p.insert (make_pair (11784, 20000));
   z_p.insert (make_pair (11000, 22600));
   z_p.insert (make_pair (10363, 25000));
   z_p.insert (make_pair (9164, 30000));
   z_p.insert (make_pair (7185, 40000));
   z_p.insert (make_pair (5574, 50000));
   z_p.insert (make_pair (4206, 60000));
   z_p.insert (make_pair (3012, 70000));
   z_p.insert (make_pair (1949, 80000));
   z_p.insert (make_pair (1457, 85000));
   z_p.insert (make_pair (988, 90000));
   z_p.insert (make_pair (540, 95000));
   z_p.insert (make_pair (111, 100000));
   z_p.insert (make_pair (0, 101325));
   z_p.insert (make_pair (-302, 105000));

}

Standard_Atmosphere::Standard_Atmosphere (const Real p_0,
                                          const Real rho_0,
                                          const Real T_0,
                                          const Real g_0,
                                          const Real R_0)
   : p_0 (p_0),
     T_0 (T_0),
     g_0 (g_0),
     R_0 (R_0),
     rho_0 (rho_0)
{
   init_p_z ();
   init_z_p ();
}

Real
Standard_Atmosphere::get_temperature_gradient (const Real z) const
{

   Real temperature_gradient = GSL_NAN;

   if (z < 11000) { temperature_gradient = -0.0065; }
   else if (z < 20000) { temperature_gradient = 0; }
   else if (z < 32000) { temperature_gradient = 0.001; }
   else if (z < 47000) { temperature_gradient = 0.0028; }
   else if (z < 51000) { temperature_gradient = 0; }
   else if (z < 71000) { temperature_gradient = -0.0028; }
   else if (z < 84852) { temperature_gradient = -0.002; }

   return temperature_gradient;

}

Real
Standard_Atmosphere::get_height (const Real pressure) const
{

   map<Real, Real>::const_iterator this_i = p_z.lower_bound (pressure);
   map<Real, Real>::const_iterator next_i = p_z.upper_bound (pressure);

   if (this_i == next_i)
   {
      if (this_i == p_z.begin ()) { next_i++; }
      else if (next_i == p_z.end ()) { this_i--; this_i--; next_i--; }
      else { this_i--; }
   }
   else if (next_i == p_z.end ()) { this_i--; next_i--; }

   const Real this_p = this_i->first;
   const Real next_p = next_i->first;
   const Real this_z = this_i->second;
   const Real next_z = next_i->second;

   const Real delta_p = next_p - this_p;
   const Real delta_z = next_z - this_z;

   return delta_z * (pressure - this_p) / delta_p + this_z;

}

Real
Standard_Atmosphere::get_pressure (const Real height) const
{

   map<Real, Real>::const_iterator this_i = z_p.lower_bound (height);
   map<Real, Real>::const_iterator next_i = z_p.upper_bound (height);

   if (this_i == next_i)
   {
      if (this_i == z_p.begin ()) { next_i++; }
      else if (next_i == z_p.end ()) { this_i--; this_i--; next_i--; }
      else { this_i--; }
   }
   else if (next_i == z_p.end ()) { this_i--; next_i--; }

   const Real this_z = this_i->first;
   const Real next_z = next_i->first;
   const Real this_p = this_i->second;
   const Real next_p = next_i->second;

   const Real delta_p = next_p - this_p;
   const Real delta_z = next_z - this_z;

   return delta_p * (height - this_z) / delta_z + this_p;

}

