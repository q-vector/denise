//
// marine.cc
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

#include "marine.h"

using namespace std;
using namespace denise;

Gd_Sea_Wave_Fetch_Coefficients::Gd_Sea_Wave_Fetch_Coefficients ()
   : a_0 (-0.07028),
     a_1 (3.32e-3),
     a_2 (-4.358e-5),
     a_3 (-4.801e-7),
     a_4 (9.969e-9),
     b_0 (0.4543),
     b_1 (0.01306),
     b_2 (-8.589e-4),
     b_3 (2.031e-5),
     b_4 (-1.605e-7)
{
}

Gd_Sea_Wave_Duration_Coefficients::Gd_Sea_Wave_Duration_Coefficients ()
   : a_0 (0.059858),
     a_1 (-0.0059951),
     a_2 (0.0001054),
     b_0 (-0.14444),
     b_1 (.0086278),
     b_2 (-0.0001056),
     c_0 (-0.26859),
     c_1 (0.03781),
     c_2 (-0.0018225),
     c_3 (0.000029153),
     d_0 (0.59191),
     d_1 (0.0079206),
     d_2 (-0.0012567),
     d_3 (0.000063899),
     d_4 (0.000001088)
{
}

Real
Gd_Sea_Wave::get_seas_from_fetch (const Real speed,
                                  const Real fetch) const
{

   const Real s = speed / 0.514444;
   const Real log_s = log (s);
   const Real log_f = log (fetch * 1e-3);

   const Real& a_0 = gd_sea_wave_fetch_coefficients.a_0;
   const Real& a_1 = gd_sea_wave_fetch_coefficients.a_1;
   const Real& a_2 = gd_sea_wave_fetch_coefficients.a_2;
   const Real& a_3 = gd_sea_wave_fetch_coefficients.a_3;
   const Real& a_4 = gd_sea_wave_fetch_coefficients.a_4;
   const Real& b_0 = gd_sea_wave_fetch_coefficients.b_0;
   const Real& b_1 = gd_sea_wave_fetch_coefficients.b_1;
   const Real& b_2 = gd_sea_wave_fetch_coefficients.b_2;
   const Real& b_3 = gd_sea_wave_fetch_coefficients.b_3;
   const Real& b_4 = gd_sea_wave_fetch_coefficients.b_4;

   const Real a = (((a_4 * s + a_3) * s + a_2) * s + a_1) * s + a_0;
   const Real b = (((b_4 * s + b_3) * s + b_2) * s + b_1) * s + b_0;

   return exp (((a * log_f + b) * log_f) + 1.1031 * log_s - 4.5285);

}

Real
Gd_Sea_Wave::get_seas_from_duration (const Real speed,
                                     const Real duration) const
{

   const Real s = speed;
   const Real log_s = log (s);
   const Real log_d = log10 (duration);

   const Real& a_0 = gd_sea_wave_duration_coefficients.a_0;
   const Real& a_1 = gd_sea_wave_duration_coefficients.a_1;
   const Real& a_2 = gd_sea_wave_duration_coefficients.a_2;
   const Real& b_0 = gd_sea_wave_duration_coefficients.b_0;
   const Real& b_1 = gd_sea_wave_duration_coefficients.b_1;
   const Real& b_2 = gd_sea_wave_duration_coefficients.b_2;
   const Real& c_0 = gd_sea_wave_duration_coefficients.b_0;
   const Real& c_1 = gd_sea_wave_duration_coefficients.c_1;
   const Real& c_2 = gd_sea_wave_duration_coefficients.c_2;
   const Real& c_3 = gd_sea_wave_duration_coefficients.c_3;
   const Real& d_0 = gd_sea_wave_duration_coefficients.d_0;
   const Real& d_1 = gd_sea_wave_duration_coefficients.d_1;
   const Real& d_2 = gd_sea_wave_duration_coefficients.d_2;
   const Real& d_3 = gd_sea_wave_duration_coefficients.d_3;
   const Real& d_4 = gd_sea_wave_duration_coefficients.d_4;

   const Real a = (a_2 * s + a_1) * s + a_0;
   const Real b = (b_2 * s + b_1) * s + b_0;
   const Real c = ((c_3 * s + c_2) * s + c_1) * s + c_0;
   const Real d = (((d_4 * s + d_3) * s + d_2) * s + d_1) * s + d_0;
   const Real e = 0.628 * log_s - 1.6685;

   const Real z = (((a * log_d + b) * log_d + c) * log_d + d) * log_d + e;
   return pow (10, z);

}

