// thermo.cc
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

#include <fstream>
#include "aviation.h"
#include "thermo.h"

Real
Moisture::get_t (const Real e_s,
                 const Real tolerance_t,
                 const Real start_t,
                 const Real end_t,
                 const Thermo_Medium thermo_medium)
{

   if (end_t - start_t < tolerance_t) { return (start_t + end_t) / 2; }
   else
   {

      Real t = (start_t + end_t) / 2;
      Real middle_e_s = Moisture::get_e_s (t, thermo_medium);

      if (middle_e_s == e_s) { return t; }
      else
      {
         if (middle_e_s > e_s)
         {
            return get_t (e_s, tolerance_t, t, end_t, thermo_medium);
         }
         else
         {
            return get_t (e_s, tolerance_t, start_t, t, thermo_medium);
         }
      }

   }

}

Real
Moisture::get_e_s (const Real t,
                   const Thermo_Medium thermo_medium)
{

   Real e_s;

   switch (thermo_medium)
   {

      default:
      case WATER:
      {

         // Goff-Gratch formula (WMO 1966 & 1979)
         Real a = 1 - 273.16 / (t + K);
         Real b = (t + K) / 273.16 - 1;
         //Real a = 1 - 273.16 / t;
         //Real b = t / 273.16 - 1;
         Real bulk_a = 10.79574 * a - 5.028 * log10 (b + 1);
         Real bulk_b = 1.50475e-4 * (1 - exp10 (-8.2969 * b));
         Real bulk_c = 0.42873e-3 * (exp10 (4.76955 * a) - 1) + 0.78614;
         e_s = exp10 (bulk_a + bulk_b + bulk_c);
         e_s *= 1e2;

         break;

      }

      case ICE:
      {
         Real a = (273.16 / (t + K));
         //Real a = (273.16 / t);
         e_s = exp10 (-9.09685*a - 3.56653*log10 (a) - 0.87682/a + 10.75981);
         e_s *= 1e2;
         break;
      }

   }

   return e_s;

}

Real
Moisture::get_e_s (const Real t,
                   const Real p,
                   const Thermo_Medium thermo_medium)
{

   const Real e_s = Moisture::get_e_s (t, thermo_medium);

   const Real bulk = (t - 12.5 + (75 / p));
   const Real f = 1 + 4.5e-8 * p + 5.6e-7 * bulk * bulk;

   return e_s * f;

}

Real
Moisture::get_r_s (const Real t,
                   const Real p,
                   const Thermo_Medium thermo_medium)
{
   const Real e_s = Moisture::get_e_s (t, p, thermo_medium);
   return epsilon * e_s / (p - e_s);
}

Real
Moisture::get_q_s (const Real t,
                   const Real p,
                   const Thermo_Medium thermo_medium)
{
   const Real e_s = get_e_s (thermo_medium);
   const Real epsilon_e_s = epsilon * e_s;
   return (epsilon_e_s) / (p - e_s + epsilon_e_s);
}

Real
Moisture::get_t (const Real e_s,
                 const Thermo_Medium thermo_medium)
{

   Real t = GSL_NAN;

   // Sargent's formula (Sargent 1980)
   Real v = log (e_s * 1e-2);
   Real a_0 = -22.59529963;
   Real a_1 = 11.33418988;
   Real a_2 = 0.575940348;
   Real a_3 = 0.3025080051e-1;
   Real a_4 = 0.1778276954e-2;
   Real a_5 = 0.7443287646e-4;
   Real a_6 = 0.1129170314e-4;
   Real t_water = (((((a_6 * v + a_5) * v + a_4) *
      v + a_3) * v + a_2) * v + a_1) * v + a_0;

   if (thermo_medium == WATER)
   {
      return t_water;
   }
   else
   {
      return get_t (e_s, 0.01, t_water - 20, t_water, thermo_medium);
   }

}

Real
Moisture::get_rh (const Real t,
                  const Real t_d,
                  const Thermo_Medium thermo_medium)
{
   return get_e_s (t_d, thermo_medium) / get_e_s (t, thermo_medium);
}

Real
Moisture::get_t_d (const Real t,
                   const Real rh,
                   const Thermo_Medium thermo_medium)
{
   return get_t (get_e_s (t - K, thermo_medium) * rh, thermo_medium) + K;
   const Real y = 1 - rh;
   const Real tt = t - K;
   return (tt - (14.55 + 0.114 * tt) * y - pow ((2.5 + 0.007 * tt) * y, 3) -
          (15.9 + 0.117 * tt) * pow (y, 14)) + K;
}

Real
Moisture::get_t_v (const Real t,
                   const Real r)
{
   if (r < 1e-6) { return t; }
   return ((t + K) * (r + epsilon) / (epsilon * (1 + r))) - K;
}

void
Thermo_Point::p_r_s_iterate (const Real p,
                             const Real r_s,
                             const Real tolerance_t,
                             const Real start_t,
                             const Real end_t,
                             const Thermo_Medium thermo_medium,
                             const Real p_0)
{

   const Real t = (start_t + end_t) / 2;
   set_t_p (t, p, p_0);

   if (end_t - start_t >= tolerance_t)
   {

      const Real middle_r_s = get_r_s (thermo_medium);

      if (middle_r_s == r_s) { return; }
      else
      {
         if (middle_r_s > r_s)
         {
            p_r_s_iterate (p, r_s, tolerance_t, start_t, t, thermo_medium, p_0);
         }
         else
         {
            p_r_s_iterate (p, r_s, tolerance_t, t, end_t, thermo_medium, p_0);
         }
      }

   }

}

Thermo_Point::Thermo_Point ()
             : temperature (this->t),
                  pressure (this->p)
{
}

Thermo_Point::Thermo_Point (const Real t,
                            const Real theta,
                            const Real p,
                            const Real p_0)
                     : p_0 (p_0),
                         t (t),
                     theta (theta),
                         p (p),
               temperature (this->t),
                  pressure (this->p)
{
}

Thermo_Point::Thermo_Point (const Thermo_Point& thermo_point)
                     : p_0 (thermo_point.p_0),
                         t (thermo_point.t),
                     theta (thermo_point.theta),
                         p (thermo_point.p),
               temperature (this->t),
                  pressure (this->p)
{
}

void
Thermo_Point::set_t_p (const Real t,
                       const Real p,
                       const Real p_0)
{
   this->t = t;
   this->theta = (t + K) * pow (p_0 / p, kappa) - K;
   this->p = p;
   this->p_0 = p_0;
}
                                                                                
Thermo_Point
Thermo_Point::t_p (const Real t,
                   const Real p,
                   const Real p_0)
{
   Thermo_Point tp;
   tp.set_t_p (t, p, p_0);
   return tp;
}
                                                                                
void
Thermo_Point::set_t_theta (const Real t,
                           const Real theta,
                           const Real p_0)
{
   this->t = t;
   this->theta = theta;
   this->p = p_0 * pow ((t + K) / (theta + K), 1 / kappa);
   this->p_0 = p_0;
}

Thermo_Point
Thermo_Point::t_theta (const Real t,
                       const Real theta,
                       const Real p_0)
{
   Thermo_Point tp;
   tp.set_t_theta (t, theta, p_0);
   return tp;
}

void
Thermo_Point::set_theta_p (const Real theta,
                           const Real p,
                           const Real p_0)
{
   this->t = (theta + K) * pow (p_0 / p, -kappa) - K;
   this->theta = theta;
   this->p = p;
   this->p_0 = p_0;
}

Thermo_Point
Thermo_Point::theta_p (const Real theta,
                       const Real p,
                       const Real p_0)
{
   Thermo_Point tp;
   tp.set_theta_p (theta, p, p_0);
   return tp;
}

void
Thermo_Point::set_t_r_s (const Real t,
                         const Real r_s,
                         const Thermo_Medium thermo_medium,
                         const Real p_0)
{
   const Real e_s = Moisture::get_e_s (t, thermo_medium);
   const Real p = e_s * epsilon / r_s;
   set_t_p (t, p, p_0);
}

Thermo_Point
Thermo_Point::t_r_s (const Real t,
                     const Real r_s,
                     const Thermo_Medium thermo_medium,
                     const Real p_0)
{
   Thermo_Point tp;
   tp.set_t_r_s (t, r_s, thermo_medium, p_0);
   return tp;
}

void
Thermo_Point::set_theta_r_s (const Real theta,
                             const Real r_s,
                             const Real tolerance_t,
                             const Real start_t,
                             const Real end_t,
                             const Thermo_Medium thermo_medium,
                             const Real p_0)
{

   const Real t = (start_t + end_t) / 2;
   set_t_theta (t, theta, p_0);

   if (end_t - start_t >= tolerance_t)
   {

      const Real middle_r_s = get_r_s (thermo_medium);

      if (middle_r_s == r_s) { return; }
      else
      {
         if (middle_r_s > r_s)
         {
            set_theta_r_s (theta, r_s, tolerance_t,
               start_t, t, thermo_medium, p_0);
         }
         else
         {
            set_theta_r_s (theta, r_s, tolerance_t,
               t, end_t, thermo_medium, p_0);
         }
      }

   }

}

Thermo_Point
Thermo_Point::theta_r_s (const Real theta,
                         const Real r_s,
                         const Real tolerance_t,
                         const Real start_t,
                         const Real end_t,
                         const Thermo_Medium thermo_medium,
                         const Real p_0)
{
   Thermo_Point tp;
   tp.set_theta_r_s (theta, r_s, tolerance_t,
      start_t, end_t, thermo_medium, p_0);
   return tp;
}

void
Thermo_Point::set_p_r_s (const Real p,
                         const Real r_s,
                         const Real tolerance_t,
                         const Thermo_Medium thermo_medium,
                         const Real p_0)
{
   const Real e_s = p * r_s / epsilon;
   const Real t = Moisture::get_t (e_s, thermo_medium);
   p_r_s_iterate (p, r_s, tolerance_t, t - 10, t + 10, thermo_medium, p_0);
}

Thermo_Point
Thermo_Point::p_r_s (const Real p,
                     const Real r_s,
                     const Real tolerance_t,
                     const Thermo_Medium thermo_medium,
                     const Real p_0)
{
   Thermo_Point tp;
   tp.set_p_r_s (p, r_s, tolerance_t, thermo_medium, p_0);
   return tp;
}

void
Thermo_Point::set_t_theta_e (const Real t,
                             const Real theta_e,
                             const Real tolerance_p,
                             const Real start_p,
                             const Real end_p,
                             const Thermo_Medium thermo_medium,
                             const Real p_0)
{

   const Real p = (start_p + end_p) / 2;
   set_t_p (t, p, p_0);

   if (end_p - start_p >= tolerance_p)
   {

      const Real middle_theta_e = get_theta_e (thermo_medium);

      if (middle_theta_e == theta_e) { return; }
      else
      {
         if (middle_theta_e > theta_e)
         {
            set_t_theta_e (t, theta_e, tolerance_p, p,
               end_p, thermo_medium, p_0);
         }
         else
         {
            set_t_theta_e (t, theta_e, tolerance_p, start_p,
               p, thermo_medium, p_0);
         }
      }

   }

}

Thermo_Point
Thermo_Point::t_theta_e (const Real t,
                         const Real theta_e,
                         const Real tolerance_p,
                         const Real start_p,
                         const Real end_p,
                         const Thermo_Medium thermo_medium,
                         const Real p_0)
{
   Thermo_Point tp;
   tp.set_t_theta_e (t, theta_e, tolerance_p,
      start_p, end_p, thermo_medium, p_0);
   return tp;
}

void
Thermo_Point::set_theta_theta_e (const Real theta,
                                 const Real theta_e,
                                 const Real tolerance_t,
                                 const Real start_t,
                                 const Real end_t,
                                 const Thermo_Medium thermo_medium,
                                 const Real p_0)
{

   const Real t = (start_t + end_t) / 2;
   set_t_theta (t, theta, p_0);

   if (end_t - start_t >= tolerance_t)
   {

      const Real middle_theta_e = get_theta_e (thermo_medium);

      if (middle_theta_e == theta_e) { return; }
      else
      {
         if (middle_theta_e > theta_e)
         {
            set_theta_theta_e (theta, theta_e, tolerance_t,
               start_t, t, thermo_medium, p_0);
         }
         else
         {
            set_theta_theta_e (theta, theta_e, tolerance_t,
               t, end_t, thermo_medium, p_0);
         }
      }

   }

}

Thermo_Point
Thermo_Point::theta_theta_e (const Real theta,
                             const Real theta_e,
                             const Real tolerance_t,
                             const Real start_t,
                             const Real end_t,
                             const Thermo_Medium thermo_medium,
                             const Real p_0)
{
   Thermo_Point tp;
   tp.set_theta_theta_e (theta, theta_e, tolerance_t,
      start_t, end_t, thermo_medium, p_0);
   return tp;
}

void
Thermo_Point::set_p_theta_e (const Real p,
                             const Real theta_e,
                             const Real tolerance_t,
                             const Real start_t,
                             const Real end_t,
                             const Thermo_Medium thermo_medium,
                             const Real p_0)
{

   const Real t = (start_t + end_t) / 2;
   set_t_p (t, p, p_0);

   if (end_t - start_t >= tolerance_t)
   {

      const Real middle_theta_e = get_theta_e (thermo_medium);

      if (middle_theta_e == theta_e) { return; }
      else
      {
         if (middle_theta_e > theta_e)
         {
            set_p_theta_e (p, theta_e, tolerance_t,
               start_t, t, thermo_medium, p_0);
         }
         else
         {
            set_p_theta_e (p, theta_e, tolerance_t,
               t, end_t, thermo_medium, p_0);
         }
      }

   }

}

Thermo_Point
Thermo_Point::p_theta_e (const Real p,
                         const Real theta_e,
                         const Real tolerance_t,
                         const Real start_t,
                         const Real end_t,
                         const Thermo_Medium thermo_medium,
                         const Real p_0)
{
   Thermo_Point tp;
   tp.set_p_theta_e (p, theta_e, tolerance_t,
      start_t, end_t, thermo_medium, p_0);
   return tp;
}

void
Thermo_Point::set_normand (const Real t,
                           const Real t_d,
                           const Real p,
                           const Real tolerance_t,
                           const Real start_t,
                           const Real end_t,
                           const Thermo_Medium thermo_medium,
                           const Real p_0)
{
   const Real theta = Thermo_Point::t_p (t, p, p_0).get_theta ();
   const Real r_s = Thermo_Point::t_p (t_d, p, p_0).get_r_s ();
   set_theta_r_s (theta, r_s, tolerance_t, start_t, end_t, thermo_medium, p_0);
}

Thermo_Point
Thermo_Point::normand (const Real t,
                       const Real t_d,
                       const Real p,
                       const Real tolerance_t,
                       const Real start_t,
                       const Real end_t,
                       const Thermo_Medium thermo_medium,
                       const Real p_0)
{
   Thermo_Point tp;
   tp.set_normand (t, t_d, p, tolerance_t, start_t, end_t, thermo_medium, p_0);
   return tp;
}

void
Thermo_Point::operator = (const Thermo_Point& thermo_point)
{
   this->p_0 = thermo_point.p_0;
   this->t = thermo_point.t;
   this->theta = thermo_point.theta;
   this->p = thermo_point.p;
}

bool
Thermo_Point::operator == (const Thermo_Point& thermo_point) const
{
   return Integer (round (this->p)) == Integer (round (thermo_point.p));
}

bool
Thermo_Point::operator > (const Thermo_Point& thermo_point) const
{
   return Integer (round (this->p)) > Integer (round (thermo_point.p));
}

bool
Thermo_Point::operator < (const Thermo_Point& thermo_point) const
{
   return Integer (round (this->p)) < Integer (round (thermo_point.p));
}

const Real&
Thermo_Point::get_p_0 () const
{
   return p_0;
}

const Real&
Thermo_Point::get_temperature () const
{
   return t;
}

const Real&
Thermo_Point::get_theta () const
{
   return theta;
}

const Real&
Thermo_Point::get_pressure () const
{
   return p;
}   

const Real&
Thermo_Point::get_t () const
{
   return temperature;
}

const Real&
Thermo_Point::get_p () const
{
   return pressure;
}   

Real
Thermo_Point::get_e_s (const Thermo_Medium thermo_medium) const
{
   return Moisture::get_e_s (t, p, thermo_medium);
}

Real
Thermo_Point::get_r_s (const Thermo_Medium thermo_medium) const
{
   return Moisture::get_r_s (t, p, thermo_medium);
}

Real
Thermo_Point::get_q_s (const Thermo_Medium thermo_medium) const
{
   return Moisture::get_q_s (t, p, thermo_medium);
}

Real
Thermo_Point::get_condense_rate (const Real omega,
                                 const Thermo_Medium thermo_medium) const
{

   if (omega >= 0) { return 0; }

   const Real dp = 10;
   const Real qs_plus = Moisture::get_q_s (t, p + dp, thermo_medium);
   const Real qs_minus = Moisture::get_q_s (t, p - dp, thermo_medium);
   return (qs_plus - qs_minus) / (dp + dp) * omega;

}


Real
Thermo_Point::get_saturated_Q_dot (const Real omega,
                                   const Thermo_Medium thermo_medium) const
{
   const Real condense_rate = get_condense_rate (omega, thermo_medium);
   return L * condense_rate;
}

Real
Thermo_Point::get_saturated_t_dot (const Real omega,
                                   const Thermo_Medium thermo_medium) const
{
   const Real q_dot = get_saturated_Q_dot (omega, thermo_medium);
   const Real alpha = R_d * t / p;
   return (q_dot + alpha * omega) / c_p;
}

Real
Thermo_Point::get_saturated_theta_dot (const Real omega,
                                       const Thermo_Medium thermo_medium) const
{
   const Real q_dot = get_saturated_Q_dot (omega, thermo_medium);
   const Real exner = pow ((p / p_0), kappa);
   return q_dot / (c_p * exner);
}

Real
Thermo_Point::get_theta_e (const Thermo_Medium thermo_medium) const
{

   // Bolton 1980
   const Real r_s = get_r_s (thermo_medium);
   const Real power = 0.2854 * (1 - 0.297 * r_s);
   const Real theta = (t + K) * pow (p_0 / p, power);
   const Real bulk = (3376 / (t + K) - 2.54) * r_s * (1 + 0.81 * r_s);
   return (theta * exp (bulk)) - K;

}

Real
Thermo_Point::get_theta_w (const Real tolerance_t,
                           const Real start_t,
                           const Real end_t,
                           const Thermo_Medium thermo_medium) const
{

   const Real theta_e = get_theta_e (thermo_medium);
   
   Thermo_Point thermo_point = Thermo_Point::p_theta_e (
      p_0, theta_e, tolerance_t, start_t, end_t, thermo_medium, p_0);

   return thermo_point.get_t ();

}


Real
Instability::get_k_index (const Real t_850,
                          const Real t_700,
                          const Real t_500,
                          const Real t_d_850,
                          const Real t_d_700)
{
   return t_850 + t_d_850 - t_700 + t_d_700 - t_500 - K;
}

Real
Instability::get_total_totals (const Real t_850,
                               const Real t_500,
                               const Real t_d_850)
{
   return t_850 + t_d_850 - t_500 - t_500;
}

Real
Instability::get_showalter_index (const Real t_850,
                                  const Real t_500,
                                  const Real t_d_850)
{
   return get_lifted_index (850e2, t_850, t_d_850, 500e2, t_500);
}

Real
Instability::get_lifted_index (const Real start_p,
                               const Real start_t,
                               const Real start_t_d,
                               const Real end_p,
                               const Real end_t)
{

   typedef Thermo_Point Tp;
   const Tp& lc = Tp::normand (start_t-K, start_t_d-K, start_p);
   const Real t_star = Tp::p_theta_e (end_p, lc.get_theta_e ()).get_t () + K;

   return end_t - t_star;

}

Real
International_Standard_Atmosphere::get_p (const Real p_0,
                                          const Real gamma,
                                          const Real dz,
                                          const Real t_0)
{

   Real p;

   if (gamma != 0)
   {
      Real t = gamma * dz + t_0;
      Real a = (t + K) / (t_0 + K);
      Real b = -g / (R * gamma);
      p = p_0 * pow (a, b);
   }
   else
   {
      p = p_0 * exp (-g * dz / (R * (t_0 + K)));
   }

   return p;

}

Real
International_Standard_Atmosphere::get_dz (const Real p,
                                           const Real p_0,
                                           const Real t_0,
                                           const Real gamma)
{
   return (t_0 / gamma) * (pow ((p / p_0), - (R * gamma / g)) - 1);
}

Integer
International_Standard_Atmosphere::get_node (const Real p) const
{

   Integer node;

   Integer size = tuple_p.size ();
   Real span = tuple_p.back () - tuple_p.front ();
   Real start_diff = p - tuple_p.front ();
   Real end_diff   = p - tuple_p.back ();

   if (start_diff * end_diff >= 0)
   {
      if (fabs (start_diff) < fabs (end_diff)) { return 0; }
      else { return size - 2; }
   }

   for (node = 0; node < size - 2; node++)
   {
      if ((p - tuple_p[node]) * (p - tuple_p[node+1]) <= 0) { break; }
   }

   return node;

}

International_Standard_Atmosphere::International_Standard_Atmosphere ()
          : tuple_z (Tuple ("0:11000:20000:32000:47000:51000:71000:84852")),
            tuple_t (Tuple ("15:-56.5:-56.5:-44.5:-2.5:-2.5:-58.5:-86.2"))
{

   Tuple tuple_p ("101325:0:0:0:0:0:0:0");

   for (Integer i = 1; i < 8; i++)
   {
      Real dz = tuple_z[i] - tuple_z[i-1];
      Real gamma = (tuple_t[i] - tuple_t[i-1]) / dz;
      tuple_p[i] = get_p (tuple_p[i-1], gamma, dz, tuple_t[i-1]);
   }

}

const Tuple&
International_Standard_Atmosphere::get_tuple_p () const
{
   return tuple_p;
}

Real
International_Standard_Atmosphere::get_z (const Real p) const
{

   Integer i = get_node (p);

   const Real& p_0 = tuple_p[i];
   const Real& t_0 = tuple_t[i];
   const Real& z_0 = tuple_z[i];

   Real gamma = (tuple_t[i+1] - t_0) / (tuple_z [i+1] - z_0);

   return Isa::get_dz (p, p_0, t_0, gamma) + z_0;

}

Real
International_Standard_Atmosphere::get_t (const Real p) const
{

   Integer i = get_node (p);

   const Real& p_0 = tuple_p[i];
   const Real& t_0 = tuple_t[i];
   const Real& z_0 = tuple_z[i];

   Real gamma = (tuple_t[i+1] - t_0) / (tuple_z [i+1] - z_0);
   Real dz = Isa::get_dz (p, p_0, t_0, gamma);

   return t_0 + gamma * dz;

}

//Real
//Sounding::get_srh (const Real c_x,
//                   const Real c_y,
//                   const Real start_p,
//                   const Real end_p,
//                   const Real delta_p) const
//{
//
//   Real srh = 0;
//   Real span_p = end_p - start_p;
//   Real dp = span_p / (rint (span_p / delta_p));
//
//   for (Real p = start_p; p < end_p - delta_p/2; p += delta_p)
//   {
//
//      Wind wind_l = get_wind_at (p);
//      Wind wind_u = get_wind_at (p + delta_p);
//
//      srh -= ((wind_u.v + wind_l.v) / 2 - c_y) * (wind_u.u - wind_l.u);
//      srh += ((wind_u.u + wind_l.u) / 2 - c_x) * (wind_u.v - wind_l.v);
//
//   }
//
//   return srh;
//
//}

//void
//Sounding::render_z_hko (const RefPtr<Context>& cr,
//                        const Thermo_Diagram& thermo_diagram) const
//{
//   render_z_hko (cr, thermo_diagram, Domain_1D (1001e2, 199e2));
//   render_z_hko (cr, thermo_diagram, Domain_1D (201e2, 69e2));
//   render_z_hko (cr, thermo_diagram, Domain_1D (71e2, 49e2));
//}
//
//void
//Sounding::render_z_hko (const RefPtr<Context>& cr,
//                        const Thermo_Diagram& thermo_diagram,
//                        const Domain_1D& domain_p) const
//{
//
//   cr->set_font_size (12);
//   cr->set_line_width (1);
//   Color (0, 0, 0).cairo (cr);
//
//   Ring ring (3);
 //  Polyline polyline;
//
//   Real offset = 0;
//   if (domain_p.start <= 201e2) { offset = 60; }
//   if (domain_p.start <= 71e2) { offset = 90; }
//
//   bool first_point = true;
//   Tuple z_tuple = get_z_tuple ();
//
//   for (Integer i = 0; i < z_tuple.size (); i++)
//   {
//
//      Real p = z_tuple[i];
//      if (domain_p.is_out_of_bounds (p)) { continue; }
//
//      Real z = get_z_at (p);
//      Real t = -z/100 + 40 + offset;
//
//      Thermo_Point thermo_point = Thermo_Point::t_p (t, p);
//      Point_2D point = thermo_diagram.transform (thermo_point);
//
//      string str = string_render ("%.0f", z);
//      ring.cairo (cr, point);
//
//      Label label (str, point, 'l', 'b', 4);
//      label.cairo (cr);
//
//      polyline.add (point, first_point);
//      if (first_point) { first_point = false; }
//
//   }
//
//   polyline.cairo (cr);
//   cr->stroke ();
//
//}

//void
//Sounding::render_hodograph (const RefPtr<Context>& cr,
//                            const Point_2D& origin,
//                            const Real scale) const
//{
//
//   Cross cross (2, 1);
//   Tuple wind_tuple = get_wind_tuple ();
//   Polar_Transform_2D transform (origin, scale);
//
//   Mesh_Attribute_2D ma_0 (Color (0.9, 0.9, 0.9), 1, M_PI/18);
//   Mesh_Attribute_2D ma_1 (Color (0.7, 0.7, 0.7), 5, M_PI/6);
//   Mesh_2D mesh_2d (Domain_1D (0, 60), Domain_1D (0, 2*M_PI), ma_0, ma_1);
//   mesh_2d.render (cr, transform, Size_2D (2, 90));
//
//   Polyline polyline;
//
//   for (Real p = wind_tuple.front (); p <= wind_tuple.back (); p += 100)
//   {
//
//      Wind wind = get_wind_at (p);
//      Real speed = -wind.get_speed ();
//      //Real direction = (450 - wind.get_direction ()) * DEGREE_TO_RADIAN;
//      Real direction = (270 + wind.get_direction ()) * DEGREE_TO_RADIAN;
//      Point_2D point_2d = transform.transform (speed, direction);
//
//      polyline.add (point_2d);
//
//   }
//
//   cr->set_line_width (2);
//   Color (1, 0.5, 0.5).cairo (cr);
//   polyline.cairo (cr);
//   cr->stroke ();
//
//   Color (0.7, 0.4, 0.4).cairo (cr);
//
//   for (Integer i = 0; i < wind_tuple.size (); i++)
//   {
//
//      Real p = wind_tuple[i];
//
//      if ((p - 1001e2) * (p - 50e2) > 0) { continue; }
//
//      Wind wind = get_wind_at (p);
//
//      Real speed = -wind.get_speed ();
//      //Real direction = (450 - wind.get_direction ()) * DEGREE_TO_RADIAN;
//      Real direction = (270 + wind.get_direction ()) * DEGREE_TO_RADIAN;
//      Point_2D point = transform.transform (speed, direction);
//
//      string text = string_render ("%.0f", p * 1e-2);
//
//      Label label (text, point, 'c', 'c');
//      label.cairo (cr);
//
//      cross.cairo (cr, point);
//      cr->fill ();
//
//   }
//
//   Wind wind = get_mean_wind (1000e2, 400e2);
//   Real speed = -wind.get_speed ();
//   //Real direction = (450 - wind.get_direction ()) * DEGREE_TO_RADIAN;
//   Real direction = (270 + wind.get_direction ()) * DEGREE_TO_RADIAN;
//   Point_2D point_2d = transform.transform (speed, direction);
//
//   Color (0.1, 1.0, 0.3).cairo (cr);
//
//   Label label ("M", point_2d, 'c', 'c');
//   label.cairo (cr);
//
//
//}

Thermo_Conserve::Thermo_Conserve (const P_Layer& p_layer)
   : P_Layer (p_layer)
{
}

const Real
Thermo_Conserve::get_t (const Real p) const
{
   return get_thermo_point (p).get_t ();
}

Isotherm::Isotherm (const Real temperature,
                    const P_Layer& p_layer)
   : Thermo_Conserve (p_layer),
     temperature (temperature)
{
}

const Real&
Isotherm::get_temperature () const
{
   return temperature;
}

void
Isotherm::set_temperature (const Real temperature)
{
   this->temperature = temperature;
}

const Thermo_Point
Isotherm::get_thermo_point (const Real p) const
{
   return Thermo_Point::t_p (temperature, p);
}

Isohume::Isohume (const Real mixing_ratio,
                  const P_Layer& p_layer,
                  const Real tolerance_t,
                  const Thermo_Medium thermo_medium)
   : Thermo_Conserve (p_layer),
     mixing_ratio (mixing_ratio),
     tolerance_t (tolerance_t),
     thermo_medium (thermo_medium)
{
}

const Real&
Isohume::get_mixing_ratio () const
{
   return mixing_ratio;
}

const Thermo_Medium
Isohume::get_thermo_medium () const
{
   return thermo_medium;
}

void
Isohume::set_mixing_ratio (const Real mixing_ratio)
{
   this->mixing_ratio = mixing_ratio;
}

const Thermo_Point
Isohume::get_thermo_point (const Real p) const
{
   return Thermo_Point::p_r_s (p, mixing_ratio, tolerance_t, thermo_medium);
}

Dry_Adiabat::Dry_Adiabat (const Real theta,
                          const P_Layer& p_layer)
   : Thermo_Conserve (p_layer),
     theta (theta)
{
}

const Real&
Dry_Adiabat::get_theta () const
{
   return theta;
}

void
Dry_Adiabat::set_theta (const Real theta)
{
   this->theta = theta;
}

const Thermo_Point
Dry_Adiabat::get_thermo_point (const Real p) const
{
   return Thermo_Point::theta_p (theta, p);
}

const Real
Dry_Adiabat::get_t_v (const Real p,
                      const Real r) const
{
   return Moisture::get_t_v (get_t (p), r);
}

Moist_Adiabat::Moist_Adiabat (const Real theta_e,
                              const P_Layer& p_layer,
                              const Real tolerance_t,
                              const Thermo_Medium thermo_medium,
                              const Real min_p)
   : Thermo_Conserve (P_Layer (max (min_p, p_layer.get_start_p ()),
                               max (min_p, p_layer.get_end_p ()))),
     theta_e (theta_e),
     tolerance_t (tolerance_t),
     thermo_medium (thermo_medium)
{
}

const Real&
Moist_Adiabat::get_theta_e () const
{
   return theta_e;
}

const Thermo_Medium
Moist_Adiabat::get_thermo_medium () const
{
   return thermo_medium;
}

void
Moist_Adiabat::set_theta_e (const Real theta_e)
{
   this->theta_e = theta_e;
}

const Thermo_Point
Moist_Adiabat::get_thermo_point (const Real p) const
{
   return Thermo_Point::p_theta_e (p, theta_e, tolerance_t);
}

const Real
Moist_Adiabat::get_t_v (const Real p) const
{
   const Thermo_Point& tp = get_thermo_point (p);
   const Real t = tp.get_t ();
   const Real r = tp.get_r_s ();
   return Moisture::get_t_v (t, r);
}

Mixed_Layer::Mixed_Layer (const Thermo_Point& normand,
                          const P_Layer& p_layer,
                          const Real tolerance_t,
                          const Thermo_Medium thermo_medium,
                          const Real min_p)
   : Thermo_Conserve (P_Layer (max (min_p, p_layer.get_start_p ()),
                               max (min_p, p_layer.get_end_p ()))),
     dry_adiabat (normand.get_theta (),
                  P_Layer (normand.get_p (),
                           max (min_p, p_layer.get_end_p ()))),
     moist_adiabat (normand.get_theta_e (thermo_medium),
                    P_Layer (normand.get_p (),
                             max (min_p, p_layer.get_end_p ())),
                    tolerance_t, thermo_medium),
     normand (normand),
     thermo_medium (thermo_medium)
{
      const Real theta_e = normand.get_theta_e (thermo_medium);
}

const bool
Mixed_Layer::is_moist (const Real p) const
{
   return (p < normand.get_p ());
}

const Thermo_Point&
Mixed_Layer::get_normand () const
{
   return normand;
}

const Thermo_Medium
Mixed_Layer::get_thermo_medium () const
{
   return moist_adiabat.get_thermo_medium ();
}

void
Mixed_Layer::set_normand (const Thermo_Point& normand)
{
   this->normand = normand;
}

const Real
Mixed_Layer::get_theta_e () const
{
   return normand.get_theta_e (thermo_medium);
}

const Thermo_Point
Mixed_Layer::get_thermo_point (const Real p) const
{
   if (is_moist (p)) { return moist_adiabat.get_thermo_point (p); }
   else              { return dry_adiabat.get_thermo_point (p); }
}

const Real
Mixed_Layer::get_t_v (const Real p) const
{
   if (is_moist (p)) { return moist_adiabat.get_t_v (p); }
   else { return dry_adiabat.get_t_v (p, normand.get_r_s ()); }
}

void
Mixed_Layer::render (const RefPtr<Context>& cr,
                     const Thermo_Diagram& thermo_diagram,
                     const bool render_isohume) const
{

   // Non Virtual Temperature Line
   const Thermo_Line ascent_tl (*this, false);
   ascent_tl.render (cr, thermo_diagram);
   cr->stroke ();

   if (!render_isohume) { return; }

   // Isohume
   const Real mixing_ratio = normand.get_r_s ();
   const P_Layer isohume_p_layer (normand.get_p (), this->get_end_p ());
   const Isohume isohume (mixing_ratio, isohume_p_layer);
   const Thermo_Line isohume_tl (isohume);
   isohume_tl.render (cr, thermo_diagram);
   cr->stroke ();

}


Draft::Draft (const Thermo_Point& normand,
              const P_Layer& p_layer)
   : Mixed_Layer (normand, p_layer),
     thermo_polygon_ptr (NULL)
{
}

Draft::~Draft ()
{
   if (thermo_polygon_ptr != NULL) { delete thermo_polygon_ptr; }
}

void
Draft::update_thermo_polygon (const Thermo_Diagram& thermo_diagram,
                              const Sounding& sounding,
                              const bool use_virtual)
{

   if (thermo_polygon_ptr != NULL) { delete thermo_polygon_ptr; }

   const Thermo_Line draft_tl (*this, use_virtual);
   const Thermo_Line environment_tl (thermo_diagram, sounding, use_virtual);
   const P_Layer& p_layer = *this;

   try
   {
      thermo_polygon_ptr = Thermo_Line::get_thermo_polygon_ptr (
         thermo_diagram, draft_tl, environment_tl, p_layer);
      thermo_polygon_ptr->simplify ();
   }
   catch (const Thermo_Exception& te)
   {
      thermo_polygon_ptr = NULL;
   }

}

Updraft::Updraft (const Thermo_Point& normand,
                  const P_Layer& p_layer)
   : Draft (normand, p_layer)
{
}

Real
Updraft::get_total_cape () const
{
   if (thermo_polygon_ptr == NULL) { return GSL_NAN; }
   return thermo_polygon_ptr->get_positive_energy ();
}

Real
Updraft::get_total_cin () const
{
   if (thermo_polygon_ptr == NULL) { return GSL_NAN; }
   return -thermo_polygon_ptr->get_negative_energy ();
}

void
Downdraft::render_dmape_area (const RefPtr<Context>& cr,
                              const Thermo_Diagram& thermo_diagram,
                              const Real line_width) const
{


   if (thermo_polygon_ptr == NULL) { return; }

   const Color color_pos (0, 1.0, 0, 0.2);
   const Color color_neg (0, 0.6, 0, 0.2);
   thermo_polygon_ptr->render (cr, color_pos, color_neg);

}

Downdraft::Downdraft (const Thermo_Point& normand,
                      const P_Layer& p_layer)
   : Draft (normand, p_layer)
{
}

Real
Downdraft::get_dmape () const
{
   if (thermo_polygon_ptr == NULL) { return GSL_NAN; }
   return -thermo_polygon_ptr->get_energy ();
}

void
Downdraft::render (const RefPtr<Context>& cr,
                   const Thermo_Diagram& thermo_diagram,
                   const Real line_width,
                   const bool use_virtual) const
{

   cr->save ();

   render_dmape_area (cr, thermo_diagram, line_width);

   Color (1, 0, 1, 1).cairo (cr);
   Dashes (Tuple ("8:2:2:2")).cairo (cr);
   cr->set_line_width (line_width);
   cr->set_line_join (LINE_JOIN_ROUND);

   const Thermo_Line thermo_line (*this, false);
   thermo_line.render (cr, thermo_diagram);
   cr->stroke ();

   if (use_virtual)
   {
      Color (1, 0.4, 1, 0.7).cairo (cr);
      const Thermo_Line thermo_line_v (*this, true);
      thermo_line_v.render (cr, thermo_diagram);
      cr->stroke ();
   }

   cr->restore ();

}

Thermo_Polygon::Thermo_Polygon (const Thermo_Diagram& thermo_diagram,
                                const bool positive)
   : thermo_diagram (thermo_diagram),
     positive (positive),
     p_layer (GSL_POSINF, GSL_NEGINF)
{
}

bool
Thermo_Polygon::is_positive () const
{
   return positive;
}

Real
Thermo_Polygon::get_energy () const
{
   return Polygon::get_area () / thermo_diagram.get_jacobian ();
}

Real
Thermo_Polygon::get_energy (const Polygon& polygon,
                            const Thermo_Diagram& thermo_diagram)
{
   return polygon.get_area () / thermo_diagram.get_jacobian ();
}

Real
Thermo_Polygon::get_positive_energy () const
{
   return Polygon::get_positive_area () / thermo_diagram.get_jacobian ();
}

Real
Thermo_Polygon::get_positive_energy (const Polygon& polygon,
                                     const Thermo_Diagram& thermo_diagram)
{
   return polygon.get_positive_area () / thermo_diagram.get_jacobian ();
}

Real
Thermo_Polygon::get_negative_energy () const
{
   return Polygon::get_negative_area () / thermo_diagram.get_jacobian ();
}

Real
Thermo_Polygon::get_negative_energy (const Polygon& polygon,
                                     const Thermo_Diagram& thermo_diagram)
{
   return polygon.get_negative_area () / thermo_diagram.get_jacobian ();
}

void
Thermo_Polygon::add (const Thermo_Point& thermo_point)
{

   const Real p = thermo_point.get_p ();

   const Real& start_p = p_layer.get_start_p ();
   const Real& end_p = p_layer.get_end_p ();
   p_layer.set (min (p, start_p), max (p, end_p));

   const Point_2D& point = thermo_diagram.transform (thermo_point);
   Polygon::add (point, false);

}

void
Thermo_Polygon::cairo (const Cairo::RefPtr<Cairo::Context>& cairo_context) const
{
   Polygon::cairo (cairo_context);
}

void
Thermo_Polygon::render (const Cairo::RefPtr<Cairo::Context>& cr,
                        const Color& color_positive,
                        const Color& color_negative) const
{

   cr->save ();

   const Polygon_Vertex* first_handle_ptr = get_first_handle_ptr ();
   Polygon_Vertex* current_handle_ptr = (Polygon_Vertex*)first_handle_ptr;

   do
   {

      const Real area = Polygon::get_simple_polygon_area (current_handle_ptr);
      const Color& color = (area > 0 ? color_positive : color_negative);

      color.cairo (cr);
      Polygon::cairo (cr, current_handle_ptr);
      cr->fill ();

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);

   cr->restore ();

}

bool
Thermo_Polygon::operator == (const Thermo_Polygon& thermo_polygon) const
{
   const Real& start_p = p_layer.get_start_p ();
   const Real& end_p = p_layer.get_end_p ();
   return (start_p + end_p) == (thermo_polygon.p_layer.get_start_p () +
      thermo_polygon.p_layer.get_end_p ());
}

bool
Thermo_Polygon::operator > (const Thermo_Polygon& thermo_polygon) const
{
   const Real& start_p = p_layer.get_start_p ();
   const Real& end_p = p_layer.get_end_p ();
   return (start_p + end_p) > (thermo_polygon.p_layer.get_start_p () +
      thermo_polygon.p_layer.get_end_p ());
}

bool
Thermo_Polygon::operator < (const Thermo_Polygon& thermo_polygon) const
{
   const Real& start_p = p_layer.get_start_p ();
   const Real& end_p = p_layer.get_end_p ();
   return (start_p + end_p) < (thermo_polygon.p_layer.get_start_p () +
      thermo_polygon.p_layer.get_end_p ());
}

void
Thermo_Line::process_bracketing_iterators (Thermo_Line::const_iterator& lb,
                                           Thermo_Line::const_iterator& ub) const
{

   if (lb == ub) // Not a a node
   {

      if (lb == lower_bound (GSL_NEGINF)) // Before first node
      {
         throw Thermo_Exception ("Out of bounds:");
      }
      else
      if (lb == end ()) // After last node
      {
         throw Thermo_Exception ("Out of bounds:");
      }
      else // Ordinary segment
      {
         lb--;
      }
   }
   else // At a node
   {

      if (ub == end ()) // Last node
      {
         lb--;
         ub--;
      }
      else // Not the last node
      {
         // do nothing
      }

   }

}

Thermo_Line::Thermo_Line ()
{
}

Thermo_Line::Thermo_Line (const Thermo_Conserve& thermo_conserve)
{

   const Tuple& tuple_p = thermo_conserve.get_tuple_p ();

   for (Tuple::const_iterator iterator = tuple_p.begin ();
        iterator != tuple_p.end (); iterator++)
   {
      const Real& p = *(iterator);
      const Real temperature = thermo_conserve.get_t (p);
      add (p, temperature);
   }

}

Thermo_Line::Thermo_Line (const Mixed_Layer& mixed_layer,
                          const bool use_virtual,
                          const Real delta_p)
{

   const bool v = use_virtual;
   const Real p_lcl = mixed_layer.get_normand ().get_p ();
   const Tuple& tuple_p = mixed_layer.get_tuple_p_specify_p (p_lcl, delta_p);

   for (Tuple::const_reverse_iterator iterator = tuple_p.rbegin ();
        iterator != tuple_p.rend (); iterator++)
   {
      const Real p = *(iterator);
      const Real t = (v ? mixed_layer.get_t_v (p) : mixed_layer.get_t (p));
      insert (make_pair (p, t));
   }

}

Thermo_Line::Thermo_Line (const Thermo_Diagram& thermo_diagram,
                          const Sounding& sounding,
                          const bool use_virtual,
                          const Thermo_Medium thermo_medium)
{

   set<Real> p_set;
   const Thermo_Medium& tm = thermo_medium;

   const Thermo_Line& t_line = sounding.get_t_line ();
   const Thermo_Line& t_d_line = sounding.get_t_d_line ();

   for (Thermo_Line::const_iterator iterator = t_line.begin ();
        iterator != t_line.end (); iterator++)
   {
      const Real p = iterator->first;
      p_set.insert (p);
   }

   for (Thermo_Line::const_iterator iterator = t_d_line.begin ();
        iterator != t_d_line.end (); iterator++)
   {
      const Real p = iterator->first;
      p_set.insert (p);
   }

   for (set<Real>::const_iterator iterator = p_set.begin ();
        iterator != p_set.end (); iterator++)
   {

      const Real p = *iterator;

      try
      {
         const Real t = use_virtual ?
            sounding.get_virtual_temperature (thermo_diagram, p, tm):
            sounding.get_temperature (thermo_diagram, p);
         if (gsl_finite (t)) { insert (make_pair (p, t)); }
      }
      catch (const Thermo_Exception& te)
      {
         try
         {
            const Real t = sounding.get_temperature (thermo_diagram, p);
            if (gsl_finite (t)) { insert (make_pair (p, t)); }
         }
         catch (const Thermo_Exception& te)
         {
         }
        
      }

   }

}

Thermo_Line::Thermo_Line (const International_Standard_Atmosphere& isa,
                          const Real delta_p)
{

   const Tuple& isa_tuple_p = isa.get_tuple_p ();

   for (Tuple::const_iterator this_iterator = isa_tuple_p.begin ();
        this_iterator != isa_tuple_p.end (); this_iterator++)
   {

      Tuple::const_iterator next_iterator = this_iterator;
      if ((++next_iterator) == isa_tuple_p.end ()) { break; }

      const bool first_segment = (this_iterator == isa_tuple_p.begin ());
      const Real this_p = *(this_iterator);
      const Real next_p = *(next_iterator);

      const P_Layer p_layer (this_p, next_p);
      const Tuple& tuple_p = p_layer.get_tuple_p (delta_p);

      for (Tuple::const_iterator iterator = tuple_p.begin ();
           iterator != tuple_p.end (); iterator++)
      {
         const Real p = *(iterator);
         const Real t = isa.get_t (p);
         const bool first_point = (iterator == tuple_p.begin ());
         if (first_segment || !first_point) { insert (make_pair (p, t)); }
      }

   }

}

void
Thermo_Line::add (const Real p,
                  const Real value)
{
   insert (make_pair (p, value));
}

void
Thermo_Line::add_node (const Thermo_Diagram& thermo_diagram,
                       const Real p)
{
   add (p, get_thermo_point (thermo_diagram, p).get_t ());
}

void
Thermo_Line::delete_node (const Thermo_Line::iterator iterator)
{

   for (Thermo_Line::iterator i = begin (); i != end (); i++)
   {
      if (i == iterator)
      {
         erase (iterator);
         break;
      }
   }

}

Thermo_Line::iterator
Thermo_Line::get_iterator (const Thermo_Diagram& thermo_diagram,
                           const Point_2D& point,
                           const Real tolerance)
{

   const Real& p_0 = thermo_diagram.get_p_0 ();

   for (Thermo_Line::iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Real& pressure = (iterator->first);
      const Real& temperature = (iterator->second);
      const Thermo_Point& tp = Thermo_Point::t_p (temperature, pressure, p_0);
      const Point_2D& p = thermo_diagram.transform (tp);

      if ((fabs (point.x - p.x) <= tolerance) &&
          (fabs (point.y - p.y) <= tolerance))
      {
         return iterator;
      }

   }

   return end ();

}

Thermo_Line::const_iterator
Thermo_Line::get_iterator (const Thermo_Diagram& thermo_diagram,
                           const Point_2D& point,
                           const Real tolerance) const
{

   const Real& p_0 = thermo_diagram.get_p_0 ();

   for (Thermo_Line::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Real& pressure = (iterator->first);
      const Real& temperature = (iterator->second);
      const Thermo_Point& tp = Thermo_Point::t_p (temperature, pressure, p_0);
      const Point_2D& p = thermo_diagram.transform (tp);

      if ((fabs (point.x - p.x) <= tolerance) &&
          (fabs (point.y - p.y) <= tolerance))
      {
         return iterator;
      }

   }

   return end ();

}

Thermo_Line::iterator
Thermo_Line::get_iterator (const Real p)
{
   return find (p);
}

Thermo_Line::const_iterator
Thermo_Line::get_iterator (const Real p) const
{
   return find (p);
}

Thermo_Point
Thermo_Line::get_thermo_point (const Thermo_Diagram& thermo_diagram,
                               const Real p) const
{

   const Real& p_0 = thermo_diagram.get_p_0 ();

   Thermo_Line::const_iterator lower_iterator = lower_bound (p);
   Thermo_Line::const_iterator upper_iterator = upper_bound (p);

   process_bracketing_iterators (lower_iterator, upper_iterator);

   const Real p_a = lower_iterator->first;
   const Real p_b = upper_iterator->first;
   const Real t_a = lower_iterator->second;
   const Real t_b = upper_iterator->second;

   const Thermo_Point& tp_aa = Thermo_Point::t_p (t_a, p_a, p_0);
   const Thermo_Point& tp_ab = Thermo_Point::t_p (t_b, p_b, p_0);

   if (fabs (p - p_a) < 0.01) { return tp_aa; }
   if (fabs (p - p_b) < 0.01) { return tp_ab; }

   Real start_t = tp_aa.get_temperature ();
   Real end_t = tp_ab.get_temperature ();
   if (fabs (start_t - end_t) < 0.01) { start_t -= 0.1;  end_t += 0.1; }

   const Integer n = Integer (fabs (end_t - start_t)) * 2 + 2;
   const Real delta_t = (end_t - start_t) / (n - 1);

   Thermo_Point tp;

   for (Integer i = 0; i < n - 1; i++)
   {

      const Real t_a = start_t + i * delta_t;
      const Real t_b = t_a + delta_t;

      const Thermo_Point& tp_ba = Thermo_Point::t_p (t_a, p, p_0);
      const Thermo_Point& tp_bb = Thermo_Point::t_p (t_b, p, p_0);

      if (thermo_diagram.intersection (tp, tp_aa, tp_ab, tp_ba, tp_bb))
      {
         return tp;
      }

   }

   const Real surface_p = rbegin ()->first;
   if (surface_p - p < 1e2)
   {
      return Thermo_Point::t_p (rbegin ()->second, surface_p, p_0);
   }

   const string& p_str = string_render ("%f", p);
   throw Thermo_Exception ("Thermo_Line::get_thermo_point ()" + p_str);

}

Real
Thermo_Line::get_nearest_p (const Real p) const
{

   Real nearest_p = GSL_NAN;
   Real min_delta_p = GSL_POSINF;

   for (Thermo_Line::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Real this_p = iterator->first;
      const Real delta_p = fabs (p - this_p);
      if (delta_p < min_delta_p)
      {
         nearest_p = this_p;
         min_delta_p = delta_p;
      }
   }

   return nearest_p;

}

Real
Thermo_Line::get_start_p () const
{
   if (size () == 0) { throw Thermo_Exception ("Empty Thermo_Line"); }
   return begin ()->first;
}

Real
Thermo_Line::get_end_p () const
{
   if (size () == 0) { throw Thermo_Exception ("Empty Thermo_Line"); }
   return rbegin ()->first;
}

void
Thermo_Line::intersection_t (const Thermo_Diagram& thermo_diagram,
                             set<Thermo_Point>& tp_set,
                             const Real t,
                             const Real lower_bound_p,
                             const Real upper_bound_p,
                             const Real tolerance_p) const
{

   Thermo_Point tp, tp_a, tp_b;

   for (Thermo_Line::const_iterator this_iterator = begin ();
        this_iterator != end (); this_iterator++)
   {

      Thermo_Line::const_iterator next_iterator = this_iterator;
      if ((++next_iterator) == end ()) { break; }

      const Real this_p = this_iterator->first;
      const Real next_p = next_iterator->first;

      try
      {
         const Thermo_Point& this_tp = get_thermo_point (thermo_diagram, this_p);
         const Thermo_Point& next_tp = get_thermo_point (thermo_diagram, next_p);

         const Real this_t = this_tp.get_t ();
         const Real next_t = next_tp.get_t ();

         const Real& p_0 = thermo_diagram.get_p_0 ();

         const Real candidate_dp = 100e2;
         const Integer n = Integer ((next_p - this_p) / candidate_dp) + 2;
         const Real delta_p = (next_p - this_p) / Real (n - 1);

         for (Integer i = 0; i < n; i++) 
         {

            const Real p_a = this_p + i * delta_p;
            const Real p_b = p_a + delta_p;
            
            tp_a.set_t_p (t, p_a, p_0);
            tp_b.set_t_p (t, p_b, p_0);

            try
            {

               if (thermo_diagram.intersection (tp, this_tp, next_tp, tp_a, tp_b))
               {
                  if (tp.get_p () < lower_bound_p) { continue; }
                  if (tp.get_p () > upper_bound_p) { continue; }
                  tp_set.insert (tp);
               }

            }
            catch (const Geometry_Exception& ge)
            {
               cerr << "Thermo_Line::intersection_t" << endl;
            }

         }

      }
      catch (const Geometry_Exception& ge)
      {
         cerr << "Thermo_Line::interserction_t 2" << ge << endl;
         continue;
      }
      catch (const Thermo_Exception& e)
      {
         cerr << e << endl;
      }
   }

}

void
Thermo_Line::intersection_theta (const Thermo_Diagram& thermo_diagram,
                                 set<Thermo_Point>& tp_set,
                                 const Real theta,
                                 const Real lower_bound_p,
                                 const Real upper_bound_p,
                                 const Real tolerance_p) const
{

   Thermo_Point tp, tp_a, tp_b;

   for (Thermo_Line::const_iterator this_iterator = begin ();
        this_iterator != end (); this_iterator++)
   {

      Thermo_Line::const_iterator next_iterator = this_iterator;
      if ((++next_iterator) == end ()) { break; }

      const Real this_p = this_iterator->first;
      const Real next_p = next_iterator->first;

      const Thermo_Point& this_tp = get_thermo_point (thermo_diagram, this_p);
      const Thermo_Point& next_tp = get_thermo_point (thermo_diagram, next_p);

      const Real this_theta = this_tp.get_theta ();
      const Real next_theta = next_tp.get_theta ();

      const Real& p_0 = thermo_diagram.get_p_0 ();

      const Real candidate_dp = 10e2;
      const Integer n = Integer ((next_p - this_p) / candidate_dp) + 2;
      const Real delta_p = (next_p - this_p) / Real (n - 1);

      for (Integer i = 0; i < n; i++) 
      {

         const Real p_a = this_p + i * delta_p;
         const Real p_b = p_a + delta_p;
            
         tp_a.set_theta_p (theta, p_a, p_0);
         tp_b.set_theta_p (theta, p_b, p_0);

         if (thermo_diagram.intersection (tp, this_tp, next_tp, tp_a, tp_b))
         {
            if (tp.get_p () < lower_bound_p) { continue; }
            if (tp.get_p () > upper_bound_p) { continue; }
            tp_set.insert (tp);
         }

      }

   }

}

void
Thermo_Line::intersection_theta_e (const Thermo_Diagram& thermo_diagram,
                                   set<Thermo_Point>& tp_set,
                                   const Real theta_e,
                                   const Real lower_bound_p,
                                   const Real upper_bound_p,
                                   const Real tolerance_p,
                                   const Thermo_Medium thermo_medium) const
{

   Thermo_Point tp, tp_a, tp_b;

   for (Thermo_Line::const_iterator this_iterator = begin ();
        this_iterator != end (); this_iterator++)
   {

      Thermo_Line::const_iterator next_iterator = this_iterator;
      if ((++next_iterator) == end ()) { break; }

      const Real this_p = this_iterator->first;
      const Real next_p = next_iterator->first;

      const Thermo_Point this_tp = get_thermo_point (thermo_diagram, this_p);
      const Thermo_Point next_tp = get_thermo_point (thermo_diagram, next_p);

      const Real this_theta_e = this_tp.get_theta_e (thermo_medium);
      const Real next_theta_e = next_tp.get_theta_e (thermo_medium);

      const Real& p_0 = thermo_diagram.get_p_0 ();
      const Real tolerance_t = 0.0005;
      const Real start_t = -120;
      const Real end_t = 80;

      const Real candidate_dp = 10e2;
      const Integer n = Integer ((next_p - this_p) / candidate_dp) + 2;
      const Real delta_p = (next_p - this_p) / Real (n - 1);

      for (Integer i = 0; i < n; i++) 
      {

         const Real p_a = this_p + i * delta_p;
         const Real p_b = p_a + delta_p;
         
         tp_a.set_p_theta_e (p_a, theta_e, tolerance_t,
            start_t, end_t, thermo_medium, p_0);
         tp_b.set_p_theta_e (p_b, theta_e, tolerance_t,
            start_t, end_t, thermo_medium, p_0);

         if (thermo_diagram.intersection (tp, this_tp, next_tp, tp_a, tp_b))
         {
            if (tp.get_p () < lower_bound_p) { continue; }
            if (tp.get_p () > upper_bound_p) { continue; }
            tp_set.insert (tp);
         }

      }

   }

}

Real
Thermo_Line::get_mean_theta (const Thermo_Diagram& thermo_diagram,
                             const P_Layer& p_layer) const
{

   Real sigma_theta = 0;
   Integer theta_n = 0;
   const Tuple& tuple_p = p_layer.get_tuple_p ();

   for (Tuple::const_iterator iterator = tuple_p.begin ();
        iterator != tuple_p.end (); iterator++)
   {

      const Real& p = *(iterator);

      try
      {
         const Thermo_Point& tp = get_thermo_point (thermo_diagram, p);
         sigma_theta += tp.get_theta ();
         theta_n++;
      }
      catch (const Thermo_Exception& te) { }

   }

   return sigma_theta / theta_n;

}

Real
Thermo_Line::get_mean_mixing_ratio (const Thermo_Diagram& thermo_diagram,
                                    const P_Layer& p_layer) const
{

   Real sigma_mixing_ratio = 0;
   Integer mixing_ratio_n = 0;
   const Tuple& tuple_p = p_layer.get_tuple_p ();

   for (Tuple::const_iterator iterator = tuple_p.begin ();
        iterator != tuple_p.end (); iterator++)
   {

      const Real p = *(iterator);

      try
      {
         const Thermo_Point& tp = get_thermo_point (thermo_diagram, p);
         sigma_mixing_ratio += tp.get_r_s ();
         mixing_ratio_n++;
      }
      catch (const Thermo_Exception& te) { }

   }

   return sigma_mixing_ratio / mixing_ratio_n;

}

void
Thermo_Line::render (const RefPtr<Context>& cr,
                     const Thermo_Diagram& thermo_diagram) const
{

   const Real& p_0 = thermo_diagram.get_p_0 ();

   for (Thermo_Line::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Real& p = iterator->first;
      const Real& t = iterator->second;
      const Thermo_Point& tp = Thermo_Point::t_p (t, p, p_0);
      const Point_2D& point = thermo_diagram.transform (tp);

      if (iterator == begin ())
      {
         cr->move_to (point.x, point.y);
      }
      else
      {
         cr->line_to (point.x, point.y);
      }

   }

}

void
Thermo_Line::render (const RefPtr<Context>& cr,
                     const Thermo_Diagram& thermo_diagram,
                     const P_Layer& p_layer,
                     const Real start_x,
                     const Real end_x) const
{

   cr->save ();

   typedef Thermo_Point Tp;

   const Real& start_p = p_layer.get_start_p ();
   const Real& end_p = p_layer.get_end_p ();

   const Tp& tp_a = thermo_diagram.get_thermo_point (start_x, start_p);
   const Tp& tp_b = thermo_diagram.get_thermo_point (end_x, start_p);
   const Tp& tp_c = thermo_diagram.get_thermo_point (start_x, end_p);
   const Tp& tp_d = thermo_diagram.get_thermo_point (end_x, end_p);
   const Real& p_0 = thermo_diagram.get_p_0 ();

   const Integer n_ab = Integer (tp_b.get_t () - tp_a.get_t ()) + 2;
   const Integer n_cd = Integer (tp_d.get_t () - tp_c.get_t ()) + 2;
   const Real delta_t_ab = (tp_b.get_t () - tp_a.get_t ()) / Real (n_ab - 1);
   const Real delta_t_cd = (tp_d.get_t () - tp_c.get_t ()) / Real (n_cd - 1);

   for (Integer i = 0; i < n_ab; i++)
   {

      const Real t = tp_a.get_t () + i * delta_t_ab;
      const Thermo_Point& tp = Tp::t_p (t, start_p, p_0);
      const Point_2D& point = thermo_diagram.transform (tp);

      if (i == 0) { cr->move_to (point.x, point.y); }
      else { cr->line_to (point.x, point.y); }

   }

   for (Integer i = n_cd - 1; i >= 0; i--)
   {

      const Real t = tp_c.get_t () + i * delta_t_cd;
      const Tp& tp = Tp::t_p (t, end_p, p_0);
      const Point_2D& point = thermo_diagram.transform (tp);

      cr->line_to (point.x, point.y);

   }

   cr->clip ();
   render (cr, thermo_diagram);
   cr->stroke ();
   cr->restore ();

}

void
Thermo_Line::render_node (const RefPtr<Context>& cr,
                          const Thermo_Diagram& thermo_diagram,
                          Thermo_Line::const_iterator iterator,
                          const bool fill,
                          const Real node_size) const
{

   const Real& p_0 = thermo_diagram.get_p_0 ();

   const Real temperature = iterator->second;
   const Real pressure = iterator->first;

   const Thermo_Point& tp = Thermo_Point::t_p (temperature, pressure, p_0);
   const Point_2D& point = thermo_diagram.transform (tp);

   Ring (node_size).cairo (cr, point);
   if (fill) { cr->fill (); } else { cr->stroke (); }

}

Real
Thermo_Line::get_top_intersection_p (const Thermo_Diagram& thermo_diagram,
                                     const Thermo_Line& thermo_line_a,
                                     const Thermo_Line& thermo_line_b)
{

   typedef Thermo_Line Tl;
   typedef Thermo_Point Tp;

   Thermo_Point ip;
   const Thermo_Line& tl_a = thermo_line_a;
   const Thermo_Line& tl_b = thermo_line_b;

   for (Tl::const_iterator this_iterator_a = tl_a.begin ();
        this_iterator_a != tl_a.end (); this_iterator_a++)
   {

      Tl::const_iterator next_iterator_a = this_iterator_a;
      if ((++next_iterator_a) == tl_a.end ()) { break; }

      const Real this_p_a = this_iterator_a->first;
      const Real next_p_a = next_iterator_a->first;

      const Tp& this_tp_a = tl_a.get_thermo_point (thermo_diagram, this_p_a);
      const Tp& next_tp_a = tl_a.get_thermo_point (thermo_diagram, next_p_a);

      for (Tl::const_iterator this_iterator_b = tl_b.begin ();
           this_iterator_b != tl_b.end (); this_iterator_b++)
      {

         Tl::const_iterator next_iterator_b = this_iterator_b;
         if ((++next_iterator_b) == tl_b.end ()) { break; }

         const Real this_p_b = this_iterator_b->first;
         const Real next_p_b = next_iterator_b->first;

         const Real a = (this_p_b - this_p_a);
         const Real b = (this_p_b - next_p_a);
         const Real c = (next_p_b - this_p_a);
         const Real d = (next_p_b - next_p_a);

         const bool no_cross =
            ((a * b) >= 0 && (c * d) >= 0 && (a * c) >= 0 && (b * d) >= 0);
         if (no_cross) { continue; }

         const Tp& this_tp_b = tl_b.get_thermo_point (thermo_diagram, this_p_b);
         const Tp& next_tp_b = tl_b.get_thermo_point (thermo_diagram, next_p_b);

         if (thermo_diagram.intersection (ip, this_tp_a,
                next_tp_a, this_tp_b, next_tp_b))
         {
            return ip.get_p ();
         }

      }

   }

   return GSL_NAN;

}

Thermo_Polygon*
Thermo_Line::get_thermo_polygon_ptr (const Thermo_Diagram& thermo_diagram,
                                     const Thermo_Line& thermo_line_up,
                                     const Thermo_Line& thermo_line_down,
                                     const P_Layer& p_layer)
{

   const Real& start_p = p_layer.get_start_p ();
   const Real& end_p = p_layer.get_end_p ();

   const Real base_t_up = thermo_line_up.rbegin ()->second;
   const Real base_t_down = thermo_line_down.rbegin ()->second;
   const bool positive = (base_t_up > base_t_down);
   Thermo_Polygon* thermo_polygon_ptr = new Thermo_Polygon (thermo_diagram, positive);

   for (Thermo_Line::const_reverse_iterator iterator = thermo_line_up.rbegin ();
        iterator != thermo_line_up.rend (); iterator++)
   {
      const Real p = iterator->first;
      if (p < start_p || p > end_p) { continue; }
      const Real t = iterator->second;
      thermo_polygon_ptr->add (Thermo_Point::t_p (t, p));
   }

   {
      const Real t_up = thermo_line_up.get_thermo_point (thermo_diagram, start_p).get_t ();
      const Real t_down = thermo_line_down.get_thermo_point (thermo_diagram, start_p).get_t ();
      thermo_polygon_ptr->add (Thermo_Point::t_p (t_up, start_p));
      thermo_polygon_ptr->add (Thermo_Point::t_p (t_down, start_p));
   }

   for (Thermo_Line::const_iterator iterator = thermo_line_down.begin ();
        iterator != thermo_line_down.end (); iterator++)
   {
      const Real p = iterator->first;
      if (p < start_p || p > end_p) { continue; }
      const Real t = iterator->second;
      thermo_polygon_ptr->add (Thermo_Point::t_p (t, p));
   }

   return thermo_polygon_ptr;

}

set<Thermo_Point>
Thermo_Line::get_tp_set_t (const Thermo_Diagram& thermo_diagram,
                           const Real t,
                           const Real tolerance_p) const
{

   set<Thermo_Point> tp_set;

   try
   {
      if (size () > 0)
      {
         const Real start_p = begin ()->first;
         const Real end_p = rbegin ()->first;
         intersection_t (thermo_diagram, tp_set, t, start_p, end_p, tolerance_p);
      }

   }
   catch (const Geometry_Exception& ge)
   {
      cerr << "Thermo_Line::get_tp-set_t Geometry_Exceptin" << endl;
   }

   return tp_set;

}

set<Thermo_Point>
Thermo_Line::get_intersection_tp_set (const Thermo_Diagram& thermo_diagram,
                                      const Mixed_Layer& mixed_layer,
                                      const Real tolerance_p) const
{

   set<Thermo_Point> tp_set;

   const Real& p_0 = thermo_diagram.get_p_0 ();
   const Thermo_Medium& thermo_medium = mixed_layer.get_thermo_medium ();
   const Thermo_Point lc_tp = mixed_layer.get_normand ();

   const Real theta = lc_tp.get_theta ();
   const Real theta_e = lc_tp.get_theta_e (thermo_medium);

   const Real lcl_p = lc_tp.get_p ();
   const Real start_p = mixed_layer.get_start_p ();
   const Real end_p = mixed_layer.get_end_p ();

   tp_set.insert (Thermo_Point::theta_p (theta, end_p, p_0));

   intersection_theta (thermo_diagram, tp_set, theta,
      lcl_p, end_p + 1e2, tolerance_p);
   intersection_theta_e (thermo_diagram, tp_set, theta_e,
      start_p - 1e2, lcl_p, tolerance_p, thermo_medium);

   return tp_set;

}

Thermo_Polygon*
Thermo_Line::get_thermo_polygon_ptr (const Thermo_Diagram& thermo_diagram,
                                     const Mixed_Layer& mixed_layer,
                                     const P_Layer& p_layer,
                                     const bool use_virtual) const
{

   const bool v = use_virtual;

   const Real p_0 = thermo_diagram.get_p_0 ();
   const Real& start_p = p_layer.get_start_p ();
   const Real& end_p = p_layer.get_end_p ();
   const Real middle_p = p_layer.get_middle_p ();
   const Real parcel_t = mixed_layer.get_thermo_point (middle_p).get_t ();
   const Real environment_t = get_thermo_point (thermo_diagram, middle_p).get_t ();
   const bool positive = (parcel_t > environment_t);

   Thermo_Polygon* thermo_polygon_ptr = new Thermo_Polygon (thermo_diagram, positive);

   // upward with the parcel
   const Real p_lcl = mixed_layer.get_normand ().get_p ();
   const Tuple& tuple_p = p_layer.get_tuple_p_specify_p (p_lcl);
   for (Tuple::const_reverse_iterator iterator = tuple_p.rbegin ();
        iterator != tuple_p.rend (); iterator++)
   {
      const Real p = *(iterator);
      const Real t = (v ? mixed_layer.get_t_v (p) : mixed_layer.get_t (p));
      thermo_polygon_ptr->add (Thermo_Point::t_p (t, p));
   }

   // First point on the environment
   {
      const Real t = get_thermo_point (thermo_diagram, start_p).get_t ();
      const Thermo_Point& tp = Thermo_Point::t_p (t, start_p, p_0);
      thermo_polygon_ptr->add (tp);
   }

   //downward with this thermo_line
   for (Thermo_Line::const_iterator iterator = lower_bound (start_p);
        iterator != upper_bound (end_p); iterator++)
   {
      const Real p = iterator->first;
      const Real t = iterator->second;
      const Thermo_Point& tp = Thermo_Point::t_p (t, p, p_0);
      thermo_polygon_ptr->add (tp);
   }

   return thermo_polygon_ptr;

}

set<Thermo_Polygon*>
Thermo_Line::get_thermo_polygon_ptr_set (const Thermo_Diagram& thermo_diagram,
                                         const Mixed_Layer& mixed_layer,
                                         const set<Thermo_Point>& tp_set) const
{

   set<Thermo_Polygon*> thermo_polygon_ptr_set;

   for (set<Thermo_Point>::const_iterator this_iterator = tp_set.begin ();              this_iterator != tp_set.end (); this_iterator++)
   {

      set<Thermo_Point>::const_iterator next_iterator = this_iterator;
      next_iterator++;
      if (next_iterator == tp_set.end ()) { break; }

      const Real this_p = this_iterator->get_p ();
      const Real next_p = next_iterator->get_p ();
      const P_Layer p_layer (this_p, next_p);

      Thermo_Polygon* thermo_polygon_ptr = get_thermo_polygon_ptr (
         thermo_diagram, mixed_layer, p_layer, false);

      thermo_polygon_ptr_set.insert (thermo_polygon_ptr);

   }

   return thermo_polygon_ptr_set;

}

T_Line::T_Line ()
{
}

T_D_Line::T_D_Line ()
{
}

T_V_Line::T_V_Line (const Thermo_Diagram& thermo_diagram,
                    const Sounding& sounding,
                    const Thermo_Medium thermo_medium)
   : Thermo_Line (thermo_diagram, sounding, true, thermo_medium)
{
}

T_W_Line::T_W_Line (const Thermo_Diagram& thermo_diagram,
                    const Sounding& sounding,
                    const Thermo_Medium& thermo_medium)
{

   set<Real> p_set;

   const Thermo_Line& t_line = sounding.get_t_line ();
   const Thermo_Line& t_d_line = sounding.get_t_d_line ();

   for (Thermo_Line::const_iterator iterator = t_line.begin ();
        iterator != t_line.end (); iterator++)
   {
      const Real p = iterator->first;
      p_set.insert (p);
   }

   for (Thermo_Line::const_iterator iterator = t_d_line.begin ();
        iterator != t_d_line.end (); iterator++)
   {
      const Real p = iterator->first;
      p_set.insert (p);
   }

   for (set<Real>::const_iterator iterator = p_set.begin ();
        iterator != p_set.end (); iterator++)
   {

      const Real p = *iterator;

      try
      {
         const Real t_w = sounding.get_wet_bulb (thermo_diagram, p);
         if (gsl_finite (t_w)) { insert (make_pair (p, t_w)); }
      }
      catch (const Thermo_Exception& te)
      {
      }

   }

}

void
Wind_Profile::process_bracketing_iterators (Wind_Profile::const_iterator& lb,
                                            Wind_Profile::const_iterator& ub) const
{

   if (lb == ub) // Not a a node
   {

      if (lb == lower_bound (GSL_NEGINF)) // Before first node
      {
         throw Thermo_Exception ("Out of bounds:");
      }
      else
      if (lb == end ()) // After last node
      {
         throw Thermo_Exception ("Out of bounds:");
      }
      else // Ordinary segment
      {
         lb--;
      }
   }
   else // At a node
   {

      if (ub == end ()) // Last node
      {
         lb--;
         ub--;
      }
      else // Not the last node
      {
         // do nothing
      }

   }

}

Wind_Profile::Wind_Profile ()
{
}

void
Wind_Profile::add (const Real p,
                   const Wind& wind)
{
   insert (make_pair (p, wind));
}

Wind
Wind_Profile::get_wind (const Real p) const
{

   Wind_Profile::const_iterator lower_iterator = lower_bound (p);
   Wind_Profile::const_iterator upper_iterator = upper_bound (p);

   process_bracketing_iterators (lower_iterator, upper_iterator);

   const Real p_a = lower_iterator->first;
   const Real p_b = upper_iterator->first;

   const Wind& wind_a = lower_iterator->second;
   const Wind& wind_b = upper_iterator->second;
   const Real u_a = wind_a.u;
   const Real v_a = wind_a.v;
   const Real u_b = wind_b.u;
   const Real v_b = wind_b.v;

   const Real u = (p - p_a) * (u_b - u_a) / (p_b - p_a) + u_a;
   const Real v = (p - p_a) * (v_b - v_a) / (p_b - p_a) + v_a;

   return Wind (u, v);

}

Wind
Wind_Profile::get_mean_wind (const P_Layer& p_layer) const
{

   Real sigma_u = 0;
   Real sigma_v = 0;
   Integer u_n = 0;
   Integer v_n = 0;
   const Tuple& tuple_p = p_layer.get_tuple_p ();

   for (Tuple::const_iterator iterator = tuple_p.begin ();
        iterator != tuple_p.end (); iterator++)
   {

      const Real p = *(iterator);

      try
      {
         const Wind& wind = get_wind (p);
         sigma_u += wind.u;
         sigma_v += wind.v;
         u_n++;
         v_n++;
      }
      catch (const Thermo_Exception& te)
      {
      }

   }

   const Real u = sigma_u / u_n;
   const Real v = sigma_v / v_n;
   return Wind (u, v);

}

Wind
Wind_Profile::get_wind_shear (const P_Layer& p_layer) const
{

   const Real& start_p = p_layer.get_start_p ();
   const Real& end_p = p_layer.get_end_p ();
   const Wind& start_wind = get_wind (start_p);
   const Wind& end_wind = get_wind (end_p);

   const Real u = end_wind.u - start_wind.u;
   const Real v = end_wind.v - start_wind.v;

   return Wind (u, v);

}

void
Height_Profile::process_bracketing_iterators (Height_Profile::const_iterator& lb,
                                              Height_Profile::const_iterator& ub) const
{

   if (lb == ub) // Not at a node
   {

      if (lb == lower_bound (GSL_NEGINF)) // Before first node
      {
         ub = lb;
         ub++;
      }
      else
      if (lb == end ()) // After last node
      {
         ub = end ();
         ub--;
         lb = ub;
         lb--;
      }
      else // Ordinary segment
      {
         lb--;
      }
   }
   else // At a node
   {

      if (ub == end ()) // Last node
      {
         lb--;
         ub--;
      }
      else // Not the last node
      {
         // do nothing
      }

   }

}

Real
Height_Profile::get_pressure (const Real this_p,
                              const Real this_z,
                              const Real next_p,
                              const Real next_z,
                              const Real height) const
{
   const Real slope = (height - this_z) / (next_z - this_z);
   return (next_p - this_p) * slope + this_p;
}

Height_Profile::Height_Profile ()
{
}

void
Height_Profile::add (const Real p,
                     const Real height)
{
   insert (make_pair (p, height));
}

Real
Height_Profile::get_height (const Real pressure) const
{

   Height_Profile::const_iterator lower_iterator = lower_bound (pressure);
   Height_Profile::const_iterator upper_iterator = upper_bound (pressure);

   process_bracketing_iterators (lower_iterator, upper_iterator);

   const Real p_a = lower_iterator->first;
   const Real p_b = upper_iterator->first;

   const Real height_a = lower_iterator->second;
   const Real height_b = upper_iterator->second;

   return (pressure - p_a) * (height_b - height_a) / (p_b - p_a) + height_a;

}

Real
Height_Profile::get_pressure (const Real height) const
{

   if (height > begin ()->second)
   {

      Height_Profile::const_iterator this_iterator = begin ();
      Height_Profile::const_iterator next_iterator = this_iterator;
      next_iterator++;

      const Real this_p = this_iterator->first;
      const Real next_p = next_iterator->first;
      const Real this_z = this_iterator->second;
      const Real next_z = next_iterator->second;

      return get_pressure (this_p, this_z, next_p, next_z, height);

   }
   else
   if (height < rbegin ()->second)
   {

      Height_Profile::const_reverse_iterator this_iterator = rbegin ();
      Height_Profile::const_reverse_iterator next_iterator = this_iterator;
      next_iterator++;

      const Real this_p = this_iterator->first;
      const Real next_p = next_iterator->first;
      const Real this_z = this_iterator->second;
      const Real next_z = next_iterator->second;

      return get_pressure (this_p, this_z, next_p, next_z, height);

   }
   else
   {

      for (Height_Profile::const_iterator this_iterator = begin ();
           this_iterator != end (); this_iterator++)
      {

         Height_Profile::const_iterator next_iterator = this_iterator;
         next_iterator++;
         if (next_iterator == end ()) { break; }

         const Real this_z = this_iterator->second;
         const Real next_z = next_iterator->second;

         if ((this_z - height) * (next_z - height) <= 0)
         {
            const Real this_p = this_iterator->first;
            const Real next_p = next_iterator->first;
            return get_pressure (this_p, this_z, next_p, next_z, height);
         }

      }

   }

   throw Thermo_Exception ("Height unknown");

}

Integer
Sounding::get_wmo_id () const
{
   return wmo_id;
}

void
Sounding::render_thermo_line_nodes (const RefPtr<Context>& cr,
                                    const Thermo_Diagram& thermo_diagram,
                                    const Thermo_Line& thermo_line,
                                    const Real node_size) const
{

   for (Thermo_Line::const_iterator iterator = thermo_line.begin ();
        iterator != thermo_line.end (); iterator++)
   {
      thermo_line.render_node (cr, thermo_diagram,
          iterator, false, node_size);
   }

}

Sounding::Sounding ()
   : t_line (),
     t_d_line (),
     wmo_id (-1),
     time (GSL_NAN),
     dry_layer (400e2, 600e2),
     shear_layer (GSL_NAN, GSL_NAN),
     steering_layer (2000, 4500),
     helicity_layer (0, 3000)
{
}

Sounding::Sounding (const Integer wmo_id)
   : t_line (),
     t_d_line (),
     wmo_id (wmo_id),
     time (GSL_NAN),
     dry_layer (400e2, 600e2),
     shear_layer (GSL_NAN, GSL_NAN),
     steering_layer (2000, 4500),
     helicity_layer (0, 3000)
{
}

Sounding::Sounding (const Sounding& sounding)
   : t_line (sounding.t_line),
     t_d_line (sounding.t_d_line),
     wind_profile (sounding.wind_profile),
     height_profile (sounding.height_profile),
     wmo_id (sounding.wmo_id),
     time (sounding.time),
     dry_layer (sounding.dry_layer),
     shear_layer (shear_layer),
     steering_layer (steering_layer),
     helicity_layer (helicity_layer)
{
}

const Dtime&
Sounding::get_time () const
{
   return time;
}

T_Line&
Sounding::get_t_line ()
{
   return t_line;
}

const T_Line&
Sounding::get_t_line () const
{
   return t_line;
}

T_D_Line&
Sounding::get_t_d_line ()
{
   return t_d_line;
}

const T_D_Line&
Sounding::get_t_d_line () const
{
   return t_d_line;
}

Wind_Profile&
Sounding::get_wind_profile ()
{
   return wind_profile;
}

const Wind_Profile&
Sounding::get_wind_profile () const
{
   return wind_profile;
}

Height_Profile&
Sounding::get_height_profile ()
{
   return height_profile;
}

const Height_Profile&
Sounding::get_height_profile () const
{
   return height_profile;
}

const T_W_Line&
Sounding::get_t_w_line () const
{
   return *t_w_line_ptr;
}

const Updraft&
Sounding::get_updraft () const
{
   return *updraft_ptr;
}

const Downdraft&
Sounding::get_downdraft () const
{
   return *downdraft_ptr;
}

void
Sounding::update_t_w_line (const Thermo_Diagram& thermo_diagram,
                           const Thermo_Medium& thermo_medium)
{
   if (t_w_line_ptr == NULL) { delete t_w_line_ptr; }
   t_w_line_ptr = new T_W_Line (thermo_diagram, *this);
}

void
Sounding::update_updraft (const Thermo_Diagram& thermo_diagram,
                          const Real start_pressure)
{

/*
   typedef Thermo_Point Tp;
   if (t_line.size () == 0) { throw Thermo_Exception ("Empty Thermo_Line"); }

   if (updraft_ptr == NULL) { delete updraft_ptr; }
   const Real start_p = get_start_p ();
   const Real end_p = std::min (start_pressure, get_end_p ());
   const P_Layer p_layer (start_p, end_p);
   updraft_ptr = new Updraft (normand, p_layer);

   const Real t = t_line.get_thermo_point (thermo_diagram, p).get_t ();
   const Real t_d = t_d_line.get_thermo_point (thermo_diagram, p).get_t ();
   const Tp& normand = Tp::normand (t, t_d, p);

   const Real start_p = t_line.get_start_p ();
   const Real end_p = t_line.get_end_p ();

   if (p == end_p)
   {
      return get_mixed_layer_surface (thermo_diagram, 50e2);
   }


   return Mixed_Layer (lc_tp, P_Layer (start_p, p));
*/

}

void
Sounding::update_updraft (const Thermo_Diagram& thermo_diagram,
                          const Thermo_Point& normand)
{
}

void
Sounding::update_downdraft (const Thermo_Diagram& thermo_diagram,
                            const Thermo_Point& start_tp)
{
}

Real
Sounding::get_dry_layer_theta_e (const Thermo_Diagram& thermo_diagram,
                                 const Real tolerance_t,
                                 const Thermo_Medium thermo_medium) const
{
   return get_mean_theta_e (thermo_diagram,
      dry_layer, tolerance_t, thermo_medium);
}

Real
Sounding::get_mean_theta_e (const Thermo_Diagram& thermo_diagram,
                            const P_Layer& p_layer,
                            const Real tolerance_t,
                            const Thermo_Medium thermo_medium) const
{

   const Real start_t = -120;
   const Real end_t = 80;

   const Real start_p = p_layer.get_start_p ();
   const Real end_p = p_layer.get_end_p ();
   const Real span_p = end_p - start_p;
   const Integer n = Integer (round (fabs (span_p) / 10e2)) + 2;
   const Real delta_p = span_p / (n - 1);

   Real sigma_theta_e = 0;
   Integer theta_e_n = 0;

   const Thermo_Line& t_line = get_t_line ();
   const Thermo_Line& t_d_line = get_t_d_line ();

   for (Integer i = 0; i < n; i++)
   {

      const Real p = start_p + i * delta_p;

      try
      {
         const Thermo_Point& t_tp = t_line.get_thermo_point (thermo_diagram, p);
         const Thermo_Point& t_d_tp = t_d_line.get_thermo_point (thermo_diagram, p);
         const Real t = t_tp.get_t ();
         const Real t_d = t_d_tp.get_t ();
         const Thermo_Point& tp = Thermo_Point::normand (
            t, t_d, p, tolerance_t, start_t, end_t, thermo_medium);
         sigma_theta_e += tp.get_theta_e (thermo_medium);
         theta_e_n++;
      }
      catch (const Thermo_Exception& te) { }

   }

   return sigma_theta_e / theta_e_n;

}

void
Sounding::add_node (const Thermo_Diagram& thermo_diagram,
                    const Real pressure)
{
   add_node_t (thermo_diagram, pressure);
   add_node_t_d (thermo_diagram, pressure);
}

void
Sounding::add_node_t (const Thermo_Diagram& thermo_diagram,
                      const Real pressure)
{
   t_line.add_node (thermo_diagram, pressure);
}

void
Sounding::add_node_t_d (const Thermo_Diagram& thermo_diagram,
                        const Real pressure)
{
   t_d_line.add_node (thermo_diagram, pressure);
}

void
Sounding::delete_thermo_node (const Thermo_Line::iterator iterator)
{
   t_line.delete_node (iterator);
   t_d_line.delete_node (iterator);
}

void
Sounding::modify_by (const Mixed_Layer& mixed_layer)
{

   typedef Thermo_Point Tp;
   typedef Thermo_Medium Tm;

   const Tp normand = mixed_layer.get_normand ();
   const Tm thermo_medium = mixed_layer.get_thermo_medium ();
   const Real start_p = mixed_layer.get_start_p ();
   const Real end_p = mixed_layer.get_end_p ();

   const Real p_lcl = normand.get_p ();
   const Real t_lcl = normand.get_t ();
   const Real theta = normand.get_theta ();
   const Real mixing_ratio = normand.get_r_s (thermo_medium);

   t_line.erase (t_line.lower_bound (start_p), t_line.upper_bound (end_p));
   t_d_line.erase (t_d_line.lower_bound (start_p), t_d_line.upper_bound (end_p));

   if (p_lcl >= start_p && p_lcl <= end_p)
   {

      t_line.add (p_lcl, t_lcl);
      t_d_line.add (p_lcl, t_lcl);

      t_line.add (end_p, Tp::theta_p (theta, end_p).get_t ());
      t_d_line.add (end_p, Tp::p_r_s (end_p, mixing_ratio).get_t ());

      const Real tolerance_t = 0.02;
      const Integer n = Integer ((p_lcl - start_p) / 10e2) + 2;
      const Real delta_p = (p_lcl - start_p) / Real (n - 1);
      const Real theta_e = normand.get_theta_e ();

      for (Integer i = 0; i < n; i++)
      {
         const Real p = start_p + i * delta_p;
         const Real t = Tp::p_theta_e (p, theta_e, tolerance_t).get_t ();
         t_line.add (p, t);
         t_d_line.add (p, t);
      }

   }
   else
   {
      t_line.add (start_p, Tp::theta_p (theta, start_p).get_t ());
      t_d_line.add (start_p, Tp::p_r_s (start_p, mixing_ratio).get_t ());

      t_line.add (end_p, Tp::theta_p (theta, end_p).get_t ());
      t_d_line.add (end_p, Tp::p_r_s (end_p, mixing_ratio).get_t ());
   }

}

void
Sounding::mix (const Thermo_Diagram& thermo_diagram,
               const P_Layer& p_layer)
{

   typedef Thermo_Point Tp;

   const Real theta = t_line.get_mean_theta (thermo_diagram, p_layer);
   const Real mixing_ratio = t_d_line.get_mean_mixing_ratio (
      thermo_diagram, p_layer);

   const Tp lc_tp = Tp::theta_r_s (theta, mixing_ratio);
   Mixed_Layer mixed_layer (lc_tp, p_layer);
   modify_by (mixed_layer);

}

void
Sounding::surface_heat (const Thermo_Diagram& thermo_diagram,
                        const Real theta)
{

   try
   {

      typedef Thermo_Point Tp;

      const Real tolerance_p = 1e2;
      const Real lb = t_line.begin ()->first;
      const Real ub = t_line.rbegin ()->first;

      typedef Thermo_Point Tp;
      set<Thermo_Point> tp_set;
      t_line.intersection_theta (thermo_diagram,
         tp_set, theta, lb, ub, tolerance_p);
      // intersection_theta throws Thermo_Exception

      const Thermo_Point intersection = *(tp_set.rbegin ());

      typedef Thermo_Point Tp;
      const Real start_p = intersection.get_p ();
      const Real end_p = ub;
      const P_Layer p_layer (start_p, end_p);
      const Real r = t_d_line.get_mean_mixing_ratio (thermo_diagram, p_layer);

      typedef Thermo_Point Tp;
      const Tp lc_tp = Tp::theta_r_s (theta, r);
      Mixed_Layer mixed_layer (lc_tp, p_layer);
      modify_by (mixed_layer);

   }
   catch (const Thermo_Exception& te)
   {
   }

}

Thermo_Line::iterator
Sounding::get_thermo_node (const Thermo_Diagram& thermo_diagram,
                           const Point_2D& point,
                           const Real tolerance)
{

   typedef Thermo_Line::iterator Thermo_Node;

   Thermo_Node node = t_line.get_iterator (thermo_diagram, point, tolerance);
   if (node != t_line.end ()) { return node; }
   else
   {
      Thermo_Node node = t_d_line.get_iterator (thermo_diagram, point, tolerance);
      if (node != t_d_line.end ()) { return node; }
      else { throw Thermo_Exception ("No nearby node."); }
   }

}

Thermo_Line::const_iterator
Sounding::get_thermo_node (const Thermo_Diagram& thermo_diagram,
                           const Point_2D& point,
                           const Real tolerance) const
{

   typedef Thermo_Line::const_iterator Thermo_Node;

   Thermo_Node node = t_line.get_iterator (thermo_diagram, point, tolerance);
   if (node != t_line.end ()) { return node; }
   else
   {
      Thermo_Node node = t_d_line.get_iterator (thermo_diagram, point, tolerance);
      if (node != t_d_line.end ()) { return node; }
      else { throw Thermo_Exception ("No nearby node."); }
   }

}

Thermo_Line::iterator
Sounding::get_thermo_line_iterator (const Thermo_Diagram& thermo_diagram,
                                    bool& is_t_line,
                                    const Point_2D& point,
                                    const Real tolerance)
{

   typedef Thermo_Line Tl;

   Tl::iterator iterator = t_line.get_iterator (
      thermo_diagram, point, tolerance);

   if (iterator != t_line.end ())
   {
      is_t_line = true;
      return iterator;
   }
   else
   {
      Tl::iterator iterator = t_d_line.get_iterator (
         thermo_diagram, point, tolerance);
      if (iterator != t_d_line.end ())
      {
         is_t_line = false;
         return iterator;
      }
      else
      {
         throw Thermo_Exception ("No nearby node.");
      }
   }

}

Thermo_Line::const_iterator
Sounding::get_thermo_line_iterator (const Thermo_Diagram& thermo_diagram,
                                    bool& is_t_line,
                                    const Point_2D& point,
                                    const Real tolerance) const
{

   typedef Thermo_Line Tl;

   Tl::const_iterator iterator = t_line.get_iterator (
      thermo_diagram, point, tolerance);

   if (iterator != t_line.end ())
   {
      is_t_line = true;
      return iterator;
   }
   else
   {
      Tl::const_iterator iterator = t_d_line.get_iterator (
         thermo_diagram, point, tolerance);
      if (iterator != t_d_line.end ())
      {
         is_t_line = false;
         return iterator;
      }
      else
      {
         throw Thermo_Exception ("No nearby node.");
      }
   }

}

Thermo_Point
Sounding::get_prev_thermo_point (const Thermo_Diagram& thermo_diagram,
                                 const Thermo_Line::const_iterator i) const
{

   Thermo_Line::const_iterator iterator = i;
   advance (iterator, -1);

   if (iterator != t_line.end () && iterator != t_d_line.end ())
   {
      const Real& p_0 = thermo_diagram.get_p_0 ();
      const Real p = iterator->first;
      const Real t = iterator->second;
      return Thermo_Point::t_p (t, p, p_0);
   }
   else
   {
      throw Thermo_Exception ("No previous node");
   }

}

Thermo_Point
Sounding::get_next_thermo_point (const Thermo_Diagram& thermo_diagram,
                                 const Thermo_Line::const_iterator i) const
{

   Thermo_Line::const_iterator iterator = i;
   advance (iterator, 1);

   if (iterator != t_line.end () && iterator != t_d_line.end ())
   {
      const Real& p_0 = thermo_diagram.get_p_0 ();
      const Real p = iterator->first;
      const Real t = iterator->second;
      return Thermo_Point::t_p (t, p, p_0);
   }
   else
   {
      throw Thermo_Exception ("No previous node");
   }

}

Real
Sounding::get_temperature (const Thermo_Diagram& thermo_diagram,
                           const Real pressure) const
{
   return t_line.get_thermo_point (thermo_diagram, pressure).get_t ();
}

Real
Sounding::get_dew_point (const Thermo_Diagram& thermo_diagram,
                         const Real pressure) const
{
   return t_d_line.get_thermo_point (thermo_diagram, pressure).get_t ();
}

Real
Sounding::get_mixing_ratio (const Thermo_Diagram& thermo_diagram,
                            const Real pressure,
                            const Thermo_Medium thermo_medium) const
{
   typedef Thermo_Point Tp;
   const Tp& tp = t_d_line.get_thermo_point (thermo_diagram, pressure);
   return tp.get_r_s (thermo_medium);
}

Real
Sounding::get_virtual_temperature (const Thermo_Diagram& thermo_diagram,
                                   const Real pressure,
                                   const Thermo_Medium thermo_medium) const
{

   const Real temperature = get_temperature (thermo_diagram, pressure);
   const Real mixing_ratio = get_mixing_ratio (
      thermo_diagram, pressure, thermo_medium);

   return Moisture::get_t_v (temperature, mixing_ratio);

}

Real
Sounding::get_wet_bulb (const Thermo_Diagram& thermo_diagram,
                        const Real pressure,
                        const Thermo_Medium thermo_medium) const
{

   const Real tolerance_t = 0.02;
   const Real start_t = -120;
   const Real end_t = 80;
   const Real& p_0 = thermo_diagram.get_p_0 ();
   const Real t = get_temperature (thermo_diagram, pressure);
   const Real t_d = get_dew_point (thermo_diagram, pressure);
   const Thermo_Point& lc_tp = Thermo_Point::normand (t, t_d, pressure);
   const Real theta_e = lc_tp.get_theta_e (thermo_medium);

   const Thermo_Point& thermo_point = Thermo_Point::p_theta_e (pressure,
      theta_e, tolerance_t, start_t, end_t, thermo_medium, p_0);

   return thermo_point.get_t ();

}

Real
Sounding::get_height (const Real pressure) const
{
   if (height_profile.size () < 2)
   {
      return Standard_Atmosphere ().get_height (pressure);
   }
   else
   {
      return height_profile.get_height (pressure);
   }
}

Real
Sounding::get_pressure (const Real height) const
{
   if (height_profile.size () < 2)
   {
      return Standard_Atmosphere ().get_pressure (height);
   }
   else
   {
      return height_profile.get_pressure (height);
   }
}

Wind
Sounding::get_wind (const Real pressure) const
{
   return wind_profile.get_wind (pressure);
}

Wind
Sounding::get_mean_wind (const P_Layer& p_layer) const
{
   return wind_profile.get_mean_wind (p_layer);
}

Wind
Sounding::get_wind_z (const Real height) const
{
   return wind_profile.get_wind (get_pressure (height));
}

Wind
Sounding::get_mean_wind (const Z_Layer& z_layer) const
{
   const Real start_p = get_pressure (z_layer.get_start_z ());
   const Real end_p = get_pressure (z_layer.get_end_z ());
   const P_Layer p_layer (start_p, end_p);
   return wind_profile.get_mean_wind (p_layer);
}

Real
Sounding::get_nearest_p (const Real p) const
{
   const Real nearest_p_t = t_line.get_nearest_p (p);
   const Real nearest_p_t_d = t_d_line.get_nearest_p (p);
   const Real delta_p_t = fabs (p - nearest_p_t);
   const Real delta_p_t_d = fabs (p - nearest_p_t_d);
   return ((delta_p_t < delta_p_t_d) ? nearest_p_t : nearest_p_t_d);
}

Real
Sounding::get_start_p () const
{
   const Real start_p_t = t_line.get_start_p ();
   const Real start_p_t_d = t_d_line.get_start_p ();
   return max (start_p_t, start_p_t_d);
}

Real
Sounding::get_end_p () const
{
   const Real end_p_t = t_line.get_end_p ();
   const Real end_p_t_d = t_d_line.get_end_p ();
   return min (end_p_t, end_p_t_d);
}

void
Sounding::render (const RefPtr<Context>& cr,
                  const Thermo_Diagram& thermo_diagram,
                  const Real node_size) const
{

   cr->save ();
   cr->set_line_width (2);

   Color (1, 0, 0).cairo (cr);
   render_t_d (cr, thermo_diagram);
   if (gsl_finite (node_size))
   {
      render_t_d_nodes (cr, thermo_diagram, node_size);
   }

   Color (0, 0, 0).cairo (cr);
   render_t (cr, thermo_diagram);
   if (gsl_finite (node_size))
   {
      render_t_nodes (cr, thermo_diagram, node_size);
   }

   render_winds (cr, thermo_diagram);
   render_heights (cr, thermo_diagram);

   cr->restore ();

}

void
Sounding::render_t (const RefPtr<Context>& cr,
                    const Thermo_Diagram& thermo_diagram) const
{
   t_line.render (cr, thermo_diagram);
   cr->stroke ();
}

void
Sounding::render_t_d (const RefPtr<Context>& cr,
                      const Thermo_Diagram& thermo_diagram) const
{
   t_d_line.render (cr, thermo_diagram);
   cr->stroke ();
}

void
Sounding::render_t_nodes (const RefPtr<Context>& cr,
                          const Thermo_Diagram& thermo_diagram,
                          const Real node_size) const
{
   render_thermo_line_nodes (cr, thermo_diagram, t_line, node_size);
}

void
Sounding::render_t_d_nodes (const RefPtr<Context>& cr,
                            const Thermo_Diagram& thermo_diagram,
                            const Real node_size) const
{
   render_thermo_line_nodes (cr, thermo_diagram, t_d_line, node_size);
}

void
Sounding::render_winds (const RefPtr<Context>& cr,
                        const Thermo_Diagram& thermo_diagram,
                        const Real x) const
{

   const Real width = thermo_diagram.get_size_2d ().i;
   const Real xx = (gsl_isnan (x) ? width - 100 : x);

   Tuple tuple_p ("1000e2:925e2:850e2:800e2:700e2:600e2");
   tuple_p.add_content ("500e2:400e2:300e2:200e2:100e2");

   cr->save ();
   Color::black ().cairo (cr);

   const bool northern_hemisphere = false;
   const Real wind_barb_size = 30;

   for (auto iterator = tuple_p.begin ();
        iterator != tuple_p.end (); iterator++)
   {

      const Real p = *(iterator);

      try
      {

         const Wind& wind = get_wind (p);
         const Thermo_Point& tp = thermo_diagram.get_thermo_point (xx, p);
         const Point_2D& point = thermo_diagram.transform (tp);

         const Wind_Barb wind_barb (wind, wind_barb_size, northern_hemisphere);

         wind_barb.cairo (cr, point);
         cr->fill ();

      }
      catch (const Thermo_Exception& te)
      {
      }

   }

   cr->restore ();

}

void
Sounding::render_heights (const RefPtr<Context>& cr,
                          const Thermo_Diagram& thermo_diagram,
                          const Real x) const
{

   const Real width = thermo_diagram.get_size_2d ().i;
   const Real xx = (gsl_isnan (x) ? width - 200 : x);

   Tuple tuple_p ("1000e2:925e2:850e2:800e2:700e2:600e2");
   tuple_p.add_content ("500e2:400e2:300e2:200e2:100e2");

   cr->save ();
   cr->set_font_size (12);
   Color::black ().cairo (cr);

   const bool show_feet = false;
   const Real multiplier = (show_feet ? 3.2808399 : 1);

   for (auto iterator = tuple_p.begin ();
        iterator != tuple_p.end (); iterator++)
   {
      try
      {
         const Real p = *(iterator);
         const Real& z = get_height (p) * multiplier;
         const Thermo_Point& tp = thermo_diagram.get_thermo_point (xx, p);
         const Point_2D& point = thermo_diagram.transform (tp);
         const string& str = string_render (show_feet ? "%.0fft" : "%.0fm", z);
         Label (str, point, 'l', 'c').cairo (cr);
      }
      catch (const Thermo_Exception& te)
      {
      }
   }

   cr->restore ();

}

Mixed_Layer
Sounding::get_mixed_layer_from (const Thermo_Diagram& thermo_diagram,
                                const Real p) const
{

   typedef Thermo_Point Tp;
   if (t_line.size () == 0) { throw Thermo_Exception ("Empty Thermo_Line"); }

   const Real start_p = t_line.get_start_p ();
   const Real end_p = t_line.get_end_p ();

   if (p == end_p)
   {
      return get_mixed_layer_surface (thermo_diagram, 50e2);
   }

   const Real t = t_line.get_thermo_point (thermo_diagram, p).get_t ();
   const Real t_d = t_d_line.get_thermo_point (thermo_diagram, p).get_t ();
   const Tp lc_tp = Tp::normand (t, t_d, p);

   return Mixed_Layer (lc_tp, P_Layer (start_p, p));

}

Mixed_Layer
Sounding::get_mixed_layer_surface (const Thermo_Diagram& thermo_diagram,
                                   const Real mixing_moisture_depth) const
{

   typedef Thermo_Point Tp;
   if (t_line.size () == 0) { throw Thermo_Exception ("Empty Thermo_Line"); }

   const Real start_p = t_line.get_start_p ();
   const Real end_p = t_line.get_end_p ();

   const Real p = end_p;
   const Real top_p = end_p - mixing_moisture_depth;
   const Real mean_r = t_d_line.get_mean_mixing_ratio (
      thermo_diagram, P_Layer (top_p, p));

   const Real t = t_line.get_thermo_point (thermo_diagram, p).get_t ();
   //const Real t_d = t_d_line.get_thermo_point (thermo_diagram, p).get_t ();
   const Real t_d = Tp::p_r_s (p, mean_r).get_t ();
   const Tp lc_tp = Tp::normand (t, t_d, p);

   return Mixed_Layer (lc_tp, P_Layer (start_p, p));

}

set<Thermo_Point>
Sounding::get_tp_set_t (const Thermo_Diagram& thermo_diagram,
                        const Real t,
                        const Real tolerance_p) const
{
   return t_line.get_tp_set_t (thermo_diagram, t, tolerance_p);
}

set<Thermo_Point>
Sounding::get_intersection_tp_set (const Thermo_Diagram& thermo_diagram,
                                   const Mixed_Layer& mixed_layer,
                                   const Thermo_Line& thermo_line,
                                   const Real tolerance_p) const
{
   return thermo_line.get_intersection_tp_set (
      thermo_diagram, mixed_layer, tolerance_p);
}

Thermo_Polygon*
Sounding::get_thermo_polygon_ptr (const Thermo_Diagram& thermo_diagram,
                                  const Mixed_Layer& mixed_layer,
                                  const Thermo_Line& thermo_line,
                                  const P_Layer& p_layer) const
{
   return thermo_line.get_thermo_polygon_ptr (
      thermo_diagram, mixed_layer, p_layer, false);
}

set<Thermo_Polygon*>
Sounding::get_thermo_polygon_ptr_set (const Thermo_Diagram& thermo_diagram,
                                      const Mixed_Layer& mixed_layer,
                                      const Thermo_Line& thermo_line,
                                      const set<Thermo_Point>& tp_set) const
{
   return thermo_line.get_thermo_polygon_ptr_set (
      thermo_diagram, mixed_layer, tp_set);
}

set<Thermo_Polygon*>
Sounding::get_thermo_polygon_ptr_set (const Thermo_Diagram& thermo_diagram,
                                      const Mixed_Layer& mixed_layer) const
{

   const Mixed_Layer& ml = mixed_layer;

   set<Thermo_Point> intersection_tp_set =
      get_intersection_tp_set (thermo_diagram, ml, t_line);
   set<Thermo_Polygon*> thermo_polygon_ptr_set = get_thermo_polygon_ptr_set (
      thermo_diagram, ml, t_line, intersection_tp_set);

   return thermo_polygon_ptr_set;

}

Thermo_Polygon*
Sounding::get_thermo_polygon_ptr (const Thermo_Diagram& thermo_diagram,
                                  const Mixed_Layer& mixed_layer,
                                  const bool use_virtual) const
{

   const Thermo_Line environment_tl (thermo_diagram, *this, use_virtual);
   const Thermo_Line mixed_layer_tl (mixed_layer, use_virtual);

   const Real start_p = Thermo_Line::get_top_intersection_p (
      thermo_diagram, mixed_layer_tl, environment_tl);
   const Real environment_end_p = environment_tl.get_end_p ();
   const Real mixed_layer_end_p = mixed_layer_tl.get_end_p ();
   const Real end_p = std::min (environment_end_p, mixed_layer_end_p);
   const P_Layer p_layer (start_p, end_p);

   Thermo_Polygon* thermo_polygon_ptr = NULL;

   try
   {
      if (gsl_finite (start_p) && gsl_finite (end_p) && (start_p < end_p))
      {

         thermo_polygon_ptr = Thermo_Line::get_thermo_polygon_ptr (
            thermo_diagram, mixed_layer_tl, environment_tl, p_layer);
         thermo_polygon_ptr->simplify ();
      }
   }
   catch (const Thermo_Exception& te)
   {
   }

   return thermo_polygon_ptr;

}

set<Thermo_Point>
Sounding::get_wbfztp_set (const Thermo_Diagram& thermo_diagram) const
{
   const T_W_Line t_w_line (thermo_diagram, *this);
   return t_w_line.get_tp_set_t (thermo_diagram, 0);
}

Real
Sounding::get_lfs (const Thermo_Diagram& thermo_diagram,
                   const Thermo_Line& downdraft_tl,
                   const Thermo_Line& environment_tl) const
{
   typedef Thermo_Line Tl;
   typedef Thermo_Point Tp;

   Thermo_Point ip;

   for (Tl::const_iterator this_iterator = environment_tl.begin ();
        this_iterator != environment_tl.end (); this_iterator++)
   {

      Tl::const_iterator next_iterator = this_iterator;
      if ((++next_iterator) == environment_tl.end ()) { break; }

      const Real this_p = this_iterator->first;
      const Real next_p = next_iterator->first;

      const Tp& this_tp = environment_tl.get_thermo_point (thermo_diagram, this_p);
      const Tp& next_tp = environment_tl.get_thermo_point (thermo_diagram, next_p);

      for (Tl::const_iterator this_iterator_d = downdraft_tl.begin ();
           this_iterator_d != downdraft_tl.end (); this_iterator_d++)
      {

         Tl::const_iterator next_iterator_d = this_iterator_d;
         if ((++next_iterator_d) == downdraft_tl.end ()) { break; }

         const Real this_p_d = this_iterator_d->first;
         const Real next_p_d = next_iterator_d->first;

         const Real a = (this_p_d - this_p);
         const Real b = (this_p_d - next_p);
         const Real c = (next_p_d - this_p);
         const Real d = (next_p_d - next_p);

         const bool no_cross =
            ((a * b) >= 0 && (c * d) >= 0 && (a * c) >= 0 && (b * d) >= 0);
         if (no_cross) { continue; }

         const Tp& this_tp_d = downdraft_tl.get_thermo_point (thermo_diagram, this_p_d);
         const Tp& next_tp_d = downdraft_tl.get_thermo_point (thermo_diagram, next_p_d);

         const Real this_t = this_tp.get_t ();
         const Real next_t = next_tp.get_t ();
         const Real this_t_d = this_tp_d.get_t ();
         const Real next_t_d = next_tp_d.get_t ();

         const Real gamma = (next_t - this_t) / (next_p - this_p);
         const Real gamma_d = (next_t_d - this_t_d) / (next_p_d - this_p_d);
         if (gamma_d > gamma) { continue; }

         if (thermo_diagram.intersection (ip, this_tp, next_tp, this_tp_d, next_tp_d))
         {
            return ip.get_p ();
         }

      }

   }

   return GSL_NAN;

}

Real
Sounding::get_cross_totals (const Thermo_Diagram& thermo_diagram) const
{

   try
   {
      Real t_d_850 = t_d_line.get_thermo_point (thermo_diagram, 850e2).get_t ();
      Real t_500 = t_line.get_thermo_point (thermo_diagram, 500e2).get_t ();
      return t_d_850 + t_500;
   }
   catch (const Thermo_Exception& te)
   {
      return GSL_NAN;
   }

}

Real
Sounding::get_vertical_totals (const Thermo_Diagram& thermo_diagram) const
{

   try
   {
      Real t_850 = t_line.get_thermo_point (thermo_diagram, 850e2).get_t ();
      Real t_500 = t_line.get_thermo_point (thermo_diagram, 500e2).get_t ();
      return t_850 - t_500;
   }
   catch (const Thermo_Exception& te)
   {
      return GSL_NAN;
   }

}

Real
Sounding::get_total_totals (const Thermo_Diagram& thermo_diagram) const
{

   try
   {
      Real t_d_850 = t_d_line.get_thermo_point (thermo_diagram, 850e2).get_t ();
      Real t_850 = t_line.get_thermo_point (thermo_diagram, 850e2).get_t ();
      Real t_500 = t_line.get_thermo_point (thermo_diagram, 500e2).get_t ();
      return t_850 + t_d_850 - t_500 - t_500;
   }
   catch (const Thermo_Exception& te)
   {
      return GSL_NAN;
   }

}

Real
Sounding::get_k_index (const Thermo_Diagram& thermo_diagram) const
{

   try
   {
      Real t_d_850 = t_d_line.get_thermo_point (thermo_diagram, 850e2).get_t ();
      Real t_d_700 = t_d_line.get_thermo_point (thermo_diagram, 700e2).get_t ();
      Real t_850 = t_line.get_thermo_point (thermo_diagram, 850e2).get_t ();
      Real t_700 = t_line.get_thermo_point (thermo_diagram, 700e2).get_t ();
      Real t_500 = t_line.get_thermo_point (thermo_diagram, 500e2).get_t ();
      return t_850 + t_d_850 - t_700 + t_d_700 - t_500;
   }
   catch (const Thermo_Exception& te)
   {
      return GSL_NAN;
   }

}

Real
Sounding::get_lifted_index (const Thermo_Diagram& thermo_diagram,
                            const Mixed_Layer& mixed_layer) const
{
   try
   {
      Real parcel_t = mixed_layer.get_thermo_point (500e2).get_t ();
      Real environment_t = t_line.get_thermo_point (thermo_diagram, 500e2).get_t ();
      return environment_t - parcel_t;
   }
   catch (const Thermo_Exception& te)
   {
      return GSL_NAN;
   }

}

Real
Sounding::get_lifted_index (const Thermo_Diagram& thermo_diagram,
                            const Real p) const
{
   try
   {
      const Mixed_Layer& mixed_layer = get_mixed_layer_from (thermo_diagram, p);
      return get_lifted_index (thermo_diagram, mixed_layer);
   }
   catch (const Thermo_Exception& te)
   {
      return GSL_NAN;
   }
}

Real
Sounding::get_showalter_index (const Thermo_Diagram& thermo_diagram) const
{
   return get_lifted_index (thermo_diagram, 850e2);
}

Real
Sounding::get_surface_lifted_index (const Thermo_Diagram& thermo_diagram) const
{
   try
   {
      const Real surface_p = get_end_p ();
      return get_lifted_index (thermo_diagram, surface_p);
   }
   catch (const Thermo_Exception& te)
   {
      return GSL_NAN;
   }
}

Real
Sounding::get_total_cin () const
{
   const Updraft& updraft = get_updraft ();
   return updraft.get_total_cin ();
}
   
Real
Sounding::get_total_cape () const
{
   const Updraft& updraft = get_updraft ();
   return updraft.get_total_cape ();
}
   
Real
Sounding::get_dmape () const
{
   const Downdraft& downdraft = get_downdraft ();
   return downdraft.get_dmape ();
}

Real
Sounding::get_convective_gust (const Real dmape_factor) const
{
/*
   const Downdraft& downdraft = get_downdraft ();
   const Real dmape = downdraft.get_dmape () * dmape_factor;
   const Real steering_speed = get_steering ().get_speed ();
   const Real ke = 0.5 * (steering_speed * steering_speed);
   return sqrt (2 * (dmape + ke));
*/
}

Real
Sounding::get_warm_cloud_depth () const
{

/*
   if (ascent_ptr == NULL) { return GSL_NAN; }

   const Thermo_Point& normand = ascent_ptr->get_normand ();
   const Real theta = normand.get_theta ();
   const Real theta_e = normand.get_theta_e ();

   const Real p_lcl = normand.get_pressure ();
   const Real p_frozen_parcel = Thermo_Point::t_theta_e (0, theta_e).get_p();

   Real warm_cloud_depth = GSL_NAN;

   if (p_frozen_parcel <= p_lcl)
   {
      const Real z_frozen_parcel = this->get_height (p_frozen_parcel);
      const Real z_lcl = this->get_height (p_lcl);
      warm_cloud_depth = z_frozen_parcel - z_lcl;
   }

   return warm_cloud_depth;
*/

}
               
Real
Sounding::get_bulk_richardson_number () const
{
}
               
Real
Sounding::get_helicity () const
{

/*
   const Real start_z = helicity_layer.get_start_z ();
   const Real end_z = helicity_layer.get_start_z ();
   const Real start_p = get_pressure (start_z);
   const Real end_p = get_pressure (end_z);
   const P_Layer p_layer (start_p, end_p);

   list<Wind> wind_list;

   try
   {
      const Wind& start_wind = wind_profile.get_wind (start_p);
      wind_list.push_back (start_wind);
   }
   catch (const Thermo_Exception& te)
   {
   }

   Wind_Profile::const_iterator begin_i = wind_profile.lower_bound (start_p);
   Wind_Profile::const_iterator end_i = wind_profile.upper_bound (end_p);

   for (Wind_Profile::const_iterator iterator = begin_i;
        iterator != end_i; iterator++)
   {
      try
      {
         const Wind& wind = iterator->second;
         wind_list.push_back (wind);
      }
      catch (const Thermo_Exception& te)
      {
      }
   }

   try
   {
      const Wind& end_wind = wind_profile.get_wind (end_p);
      wind_list.push_back (end_wind);
   }
   catch (const Thermo_Exception& te)
   {
   }

   if (wind_list.size () == 0)
   {
      return GSL_NAN;
   }

   Real helicity = 0;

   for (list<Wind>::const_iterator this_i = wind_list.begin ();
        this_i != wind_list.end (); this_i++)
   {

      list<Wind>::const_iterator next_i = this_i;
      if ((++next_i) == wind_list.end ()) { break; }

      const Wind& this_wind = *(this_i);
      const Wind& next_wind = *(next_i);

      const Real u_0 = this_wind.u - storm_motion.u;
      const Real v_0 = this_wind.v - storm_motion.v;
      const Real u_1 = next_wind.u - this_wind.u;
      const Real v_1 = next_wind.v - this_wind.v;

      const Real d_helicity = (u_0 * v_1 - u_1 * v_0);
      helicity += d_helicity;

   }

   return helicity;
*/

}

Real
Thermo_Diagram::find_temperature (const Real x,
                                  const Real p,
                                  const Real start_t,
                                  const Real end_t,
                                  const Real tolerance_t) const
{

   const Real t = (start_t + end_t) / 2;
   if (end_t - start_t < tolerance_t) { return t; }

   const Thermo_Point& thermo_point = Thermo_Point::t_p (t, p);
   const Real middle_x = transform (thermo_point).x;
   if (middle_x == x) { return t; }

   const Real lower_t = (middle_x > x ? start_t : t);
   const Real upper_t = (middle_x > x ? t : end_t);

   return find_temperature (x, p, lower_t, upper_t, tolerance_t);
   
}

Thermo_Point
Thermo_Diagram::get_thermo_point (const Real x,
                                  const Real p) const
{
   const Real temperature = find_temperature (x, p);
   return Thermo_Point::t_p (temperature, p);
}

void
Thermo_Diagram::render_isobar (const RefPtr<Context> cr,
                               const Real p,
                               const Real start_t,
                               const Real delta_t,
                               const Integer n) const
{

   for (Integer i = 0; i < n; i++)
   {

      const Real t = start_t + i * delta_t;
      const Thermo_Point& tp = Thermo_Point::t_p (t, p, p_0);
      const Point_2D& point = transform (tp);

      if (i == 0)
      {
         cr->move_to (point.x, point.y);
      }
      else
      {
         cr->line_to (point.x, point.y);
      }

   }

}

void
Thermo_Diagram::render_isobars (const RefPtr<Context> cr,
                                const Color& color_0,
                                const Color& color_1,
                                const Color& color_2,
                                const bool fine) const
{

   const Real start_t = -90;
   const Real end_t = 60;
   const Integer n = Integer (fabs (end_t - start_t)) + 2;
   const Real delta_t = (end_t - start_t) / Real (n - 1);

   const Real delta_p_s = (fine ? 10e2 : 100e2);
   const Real delta_p_u = (fine ? 1e2 : 10e2);
   const Color& color = (fine ? color_2 : color_0);

   for (Real p = 100e2; p <= 1050e2; p += delta_p_s)
   {
      render_isobar (cr, p, start_t, delta_t, n);
   }

   for (Real p = 50e2; p < 100e2; p += delta_p_u)
   {
      render_isobar (cr, p, start_t, delta_t, n);
   }

   color.cairo (cr);
   cr->stroke ();

}

void
Thermo_Diagram::render_isotherms (const RefPtr<Context> cr,
                                  const Color& color_0,
                                  const Color& color_1,
                                  const Color& color_2,
                                  const bool fine) const
{


   const Real delta_p = 10e2;
   const Real start_p = 50e2;
   const Real end_p = 1050e2;

   const Real start_t = -90;
   const Real end_t = 60;
   const Real delta_t = (fine ? 1 : 10);

   const P_Layer p_layer (start_p, end_p);

   for (Real t = start_t; t <= end_t + delta_t/2; t += delta_t)
   {
      const Isotherm isotherm (t, p_layer);
      const Thermo_Line thermo_line (isotherm);
      thermo_line.render (cr, *this);
   }

   const Color& color = (fine ? color_2 : color_0);
   color.cairo (cr);
   cr->stroke ();

}

void
Thermo_Diagram::render_dry_adiabats (const RefPtr<Context>& cr,
                                     const Color& color_0,
                                     const Color& color_1,
                                     const Color& color_2) const
{

   color_1.cairo (cr);

   const Real delta_p = 10e2;

   for (Real theta = -40; theta <= 60.1; theta += 10)
   {
      const Real start_p = 1050e2;
      const Real end_p = Thermo_Point::t_theta (-90, theta).get_p ();
      const P_Layer p_layer (start_p, end_p);
      const Dry_Adiabat dry_adiabat (theta, p_layer);
      const Thermo_Line thermo_line (dry_adiabat);
      thermo_line.render (cr, *this);
   }

   for (Real theta = 70; theta <= 150.1; theta += 10)
   {
      const Real start_p = Thermo_Point::t_theta (-90, theta).get_p ();
      const Real end_p = Thermo_Point::t_theta (70, theta).get_p ();
      const P_Layer p_layer (start_p, end_p);
      const Dry_Adiabat dry_adiabat (theta, p_layer);
      const Thermo_Line thermo_line (dry_adiabat);
      thermo_line.render (cr, *this);
   }

   for (Real theta = 160; theta <= 250.1; theta += 10)
   {
      const Real start_p = 50e2;
      const Real end_p = Thermo_Point::t_theta (70, theta).get_p ();
      const P_Layer p_layer (start_p, end_p);
      const Dry_Adiabat dry_adiabat (theta, p_layer);
      const Thermo_Line thermo_line (dry_adiabat);
      thermo_line.render (cr, *this);
   }

   cr->stroke ();

}

void
Thermo_Diagram::render_saturated_adiabats (const RefPtr<Context>& cr,
                                           const Color& color_0,
                                           const Color& color_1,
                                           const Color& color_2,
                                           const bool fine) const
{

   const Real end_t = -50;
   const Real start_p = 1000e2;
   const Real delta_t = (fine ? 2 : 10);
   const Color& color = (fine ? color_2 : color_1);

   color.cairo (cr);

   for (Real t = -40; t <= 60.1; t += delta_t)
   {

      const Integer up_n = Integer (fabs (t - end_t));
      const Real& theta_e = Thermo_Point::t_p (t, p_0).get_theta_e (WATER);
      const Real end_pp = Thermo_Point::t_theta_e (end_t,
         theta_e, 2, 1e2, 1200e2, WATER, p_0).get_p ();
      const Real end_p = (end_pp < 50e2 ? 50e2 : end_pp);

      const P_Layer p_layer (start_p, end_p);
      const Moist_Adiabat moist_adiabat (theta_e, p_layer, 0.02, WATER);
      const Thermo_Line thermo_line (moist_adiabat);
      thermo_line.render (cr, *this);

   }

   cr->stroke (); 

}

void
Thermo_Diagram::render_isohumes (const RefPtr<Context>& cr,
                                 const Color& color_0,
                                 const Color& color_1,
                                 const Color& color_2) const
{

   cr->save ();
   color_1.cairo (cr);
   Dashes ("4:2").cairo (cr);

   Tuple r_s_tuple ("0.001e-3:0.002e-3:0.005e-3:0.01e-3:0.02e-3:0.03e-3");
   r_s_tuple.add_content ("0.05e-3:0.1e-3:0.15e-3:0.2e-3:0.3e-3:0.4e-3");
   r_s_tuple.add_content ("0.5e-3:0.6e-3:0.8e-3:1e-3:1.5e-3:2e-3:2.5e-3");
   r_s_tuple.add_content ("3e-3:4e-3:5e-3:6e-3:7e-3:8e-3:9e-3:10e-3");
   r_s_tuple.add_content ("12e-3:14e-3:16e-3:18e-3:20e-3:24e-3:28e-3");
   r_s_tuple.add_content ("32e-3:36e-3:40e-3:44e-3:48e-3:52e-3:56e-3");
   r_s_tuple.add_content ("60e-3:68e-3:80e-3");

   for (Tuple::iterator iterator = r_s_tuple.begin ();
        iterator != r_s_tuple.end (); iterator++)
   {

      const Real& r_s = *(iterator);
      render_isohume (cr, r_s);
   }

   cr->restore ();

}

void
Thermo_Diagram::render_labels (const RefPtr<Context>& cr,
                               const Real label_x,
                               const Color& color_0,
                               const Color& color_1,
                               const Color& color_2) const
{

   color_0.cairo (cr);
   cr->set_font_size (label_size);

   for (Real t = -40; t <= 70.1; t += 10)
   {
    
      const Thermo_Point& thermo_point = Thermo_Point::t_p (t, 1000e2);
      const Point_2D& point = transform (thermo_point);
      const string& text = string_render ("%.0f", t);
      Label (text, point, 'r', 't', 2).cairo (cr);
   }

   for (Real t = -80; t <= 0.1; t += 10)
   {
      const Thermo_Point& thermo_point = Thermo_Point::t_p (t, 190e2);
      const Point_2D& point = transform (thermo_point);
      const string& text = string_render ("%.0f", t);
      Label (text, point, 'c', 'c').cairo (cr);
   }

   Tuple tuple_p (5, Real (50e2), Real (90e2));
   tuple_p.add_content (2, Real (100e2), Real (150e2));
   tuple_p.add_content (9, Real (200e2), Real (1000e2));

   //const Real width = Real (transform_ptr->size_2d.i);
   //const Real x = width * 0.92;

   typedef Thermo_Point Tp;

   for (Tuple::const_iterator iterator = tuple_p.begin ();
        iterator != tuple_p.end (); iterator++)
   {

      const Real& p = *(iterator);
      const string& text = string_render ("%.0fhPa", p * 1e-2);

      const Thermo_Point& tp = get_thermo_point (label_x, p);
      const Point_2D& point = transform (tp);
      Label (text, point, 'c', 'c', label_size / 7).cairo (cr);

//      Thermo_Point thermo_point = Thermo_Point::t_p (-90, p);
//      Point_2D point = transform.transform (thermo_point);
//
//      if (point.x < 0)
//      {
//         thermo_point = transform.get_thermo_point (0, p);
//         point = transform.transform (thermo_point);
//         Label (text, point, 'l', 'b', 2).cairo (cr);
//      }
//      else
//      {
//         Label (text, point, 'r', 'b', 2).cairo (cr);
//      }

   }

   for (Real theta = -40; theta <= 220.1; theta += 20)
   {
      const Tp& thermo_point = Tp::t_theta (-45, theta);
      const Point_2D& point = transform (thermo_point);
      const string& text = string_render ("%.0f", theta);
      Label (text, point, 'c', 'c').cairo (cr);
   }

   color_1.cairo (cr);
   cr->set_font_size (label_size * 0.75);

   const Real p_a = 1050e2;
   const Real p_b = 500e2;
   const Real padding = label_size * 0.4;

   Tuple r_s_tuple ("0.001e-3:0.002e-3:0.005e-3:0.01e-3:0.02e-3:0.03e-3");
   r_s_tuple.add_content ("0.05e-3:0.1e-3:0.15e-3:0.2e-3:0.3e-3:0.4e-3");
   r_s_tuple.add_content ("0.5e-3:0.6e-3:0.8e-3:1e-3:1.5e-3:2e-3:2.5e-3");
   r_s_tuple.add_content ("3e-3:4e-3:5e-3:6e-3:7e-3:8e-3:9e-3:10e-3");
   r_s_tuple.add_content ("12e-3:14e-3:16e-3:18e-3:20e-3:24e-3:28e-3");
   r_s_tuple.add_content ("32e-3:36e-3:40e-3:44e-3:48e-3:52e-3:56e-3");
   r_s_tuple.add_content ("60e-3:68e-3:80e-3");

   for (Tuple::iterator iterator = r_s_tuple.begin ();
        iterator != r_s_tuple.end (); iterator++)
   {

      const Real& r_s = *(iterator);
      const string& text = string_render ("%g", r_s * 1e3);

      const Tp& tp_a = Tp::p_r_s (p_a, r_s, 0.02, WATER, p_0);
      const Real t_a = tp_a.get_t ();

      const Point_2D& point_a = transform (Tp::t_p (t_a, p_a));
      Label label_a (text, point_a, 'r', 'c', padding);
      label_a.set_text_angle (-M_PI/3);
      label_a.cairo (cr);

      if (r_s < 0.1e-3) { continue; }

      const Tp& tp_b = Tp::p_r_s (p_b, r_s, 0.02, WATER, p_0);
      const Real t_b = tp_b.get_t ();

      const Point_2D& point_b = transform (Tp::t_p (t_b, p_b));
      Label label_b (text, point_b, 'r', 'c', padding);
      label_b.set_text_angle (-M_PI/3);
      label_b.cairo (cr);

   }

}

Thermo_Diagram::Thermo_Diagram (const Size_2D& size_2d,
                                const Real p_0,
                                const Thermo_Point& ref_thermo_point)
   : size_2d (size_2d),
     ref_thermo_point (ref_thermo_point),
     label_size (size_2d.i * 0.014),
     p_0 (p_0)
{
}

Thermo_Diagram::~Thermo_Diagram ()
{
}

const Size_2D&
Thermo_Diagram::get_size_2d () const
{
   return size_2d;
}

void
Thermo_Diagram::set_anchor (const Point_2D& anchor)
{
   this->anchor = anchor;
}

void
Thermo_Diagram::zoom (const Point_2D& point, 
                      const Real scale)
{
   Affine_Transform_2D::translate (-point.x, -point.y);
   Affine_Transform_2D::scale (scale, scale);
   Affine_Transform_2D::translate (point.x, point.y);
}

const Real
Thermo_Diagram::get_p_0 () const
{
   return p_0;
}

void
Thermo_Diagram::render (const RefPtr<Context>& cr,
                        const Real label_x,
                        const Color& color_0,
                        const Color& color_1,
                        const Color& color_2,
                        const bool fine) const
{

   cr->save ();

   render_isohumes (cr, color_0, color_1, color_2);

   if (fine)
   {
      render_isobars (cr, color_0, color_1, color_2, true);
      render_isotherms (cr, color_0, color_1, color_2, true);
      render_saturated_adiabats (cr, color_0, color_1, color_2, true);
   }

   render_dry_adiabats (cr, color_0, color_1, color_2);
   render_saturated_adiabats (cr, color_0, color_1, color_2, false);
   render_isobars (cr, color_0, color_1, color_2, false);
   render_isotherms (cr, color_0, color_1, color_2, false);
   render_labels (cr, label_x, color_0, color_1, color_2);

   cr->restore ();

}

Point_2D
Thermo_Diagram::transform (const Thermo_Point& thermo_point) const
{
   Point_2D point_2d;
   transform (point_2d, thermo_point);
   return point_2d;
}

Thermo_Point
Thermo_Diagram::get_thermo_point (const Point_2D& point_2d) const
{
   Thermo_Point thermo_point;
   reverse_tp (thermo_point, point_2d);
   return thermo_point;
}

bool
Thermo_Diagram::intersection (Thermo_Point& intersection,
                              const Thermo_Point& tp_aa,
                              const Thermo_Point& tp_ab,
                              const Thermo_Point& tp_ba,
                              const Thermo_Point& tp_bb) const
{

   Point_2D ip;
   const Point_2D& p_aa = transform (tp_aa);
   const Point_2D& p_ab = transform (tp_ab);
   const Point_2D& p_ba = transform (tp_ba);
   const Point_2D& p_bb = transform (tp_bb);

   if (!Edge::intersection (ip, p_aa, p_ab, p_ba, p_bb)) { return false; }
   else
   {
      reverse_tp (intersection, ip);
      return true;
   }

}

void
Thermo_Diagram::render_isohume (const RefPtr<Context>& cr,
                                const Real r_s,
                                const Real start_p,
                                const Real end_p) const
{

   const Real end_t = -49.5;
   const Real end_pp = std::min (50e2, (gsl_finite (end_p) ? end_p:
      Thermo_Point::t_r_s (end_t, r_s, WATER, p_0).get_p ()));

   const Real delta_p = -50;
   const Real tolerance_t = 0.02;
   const Integer n = Integer ((end_pp - start_p) / delta_p) + 2;

   const P_Layer p_layer (start_p, end_pp);
   const Isohume isohume (r_s, p_layer, tolerance_t, WATER);
   const Thermo_Line thermo_line (isohume);
   thermo_line.render (cr, *this);

   cr->stroke ();

}

Real
Tephigram::get_jacobian () const
{
   return Affine_Transform_2D::get_jacobian () / -c_p;
}

void
Tephigram::transform (Point_2D& point_2d,
                      const Thermo_Point& thermo_point) const
{
   const Real t = thermo_point.get_t ();
   const Real log_theta = log (thermo_point.get_theta () + K);
   Affine_Transform_2D::transform (point_2d.x, point_2d.y, t, log_theta);
}

void
Tephigram::reverse_tp (Thermo_Point& thermo_point,
                       const Point_2D& point) const
{
   Real t, log_theta;
   Affine_Transform_2D::reverse (t, log_theta, point.x, point.y);
   thermo_point.set_t_theta (t, exp (log_theta) - K, p_0);
}

Tephigram::Tephigram (const Size_2D& size_2d,
                      const Real p_0,
                      const Thermo_Point& ref_thermo_point)
   : Thermo_Diagram (size_2d, p_0, ref_thermo_point)
{
   reset (size_2d);
}

void
Tephigram::reset (const Size_2D& size_2d)
{

   this->size_2d = size_2d;

   const Real scale = size_2d.i / 150.0;
   const Real x = ref_thermo_point.get_t ();
   const Real y = log (ref_thermo_point.get_theta () + K);

   const Real anchor_x = (size_2d.i / 33.0);
   const Real anchor_y = size_2d.j - anchor_x * 2;
   set_anchor (Point_2D (anchor_x, anchor_y));

   Affine_Transform_2D::set_identity ();
   Affine_Transform_2D::translate (-x, -y);
   Affine_Transform_2D::scale (scale, -K * scale);
   Affine_Transform_2D::rotate (M_PI/4);
   Affine_Transform_2D::translate (anchor.x, anchor.y);

}

Real
Emagram::get_jacobian () const
{
   return Affine_Transform_2D::get_jacobian () / R_d;
}

void
Emagram::transform (Point_2D& point_2d,
                    const Thermo_Point& thermo_point) const
{
   const Real t = thermo_point.get_t ();
   const Real log_p = log (thermo_point.get_p ());
   Affine_Transform_2D::transform (point_2d.x, point_2d.y, t, log_p);
}

void
Emagram::reverse_tp (Thermo_Point& thermo_point,
                     const Point_2D& point) const
{
   Real t, log_p;
   Affine_Transform_2D::reverse (t, log_p, point.x, point.y);
   thermo_point.set_t_p (t, exp (log_p), p_0);
}

Emagram::Emagram (const Size_2D& size_2d,
                  const Real magic_ratio,
                  const Real p_0,
                  const Thermo_Point& ref_thermo_point)
   : Thermo_Diagram (size_2d, p_0, ref_thermo_point),
     magic_ratio (magic_ratio)
{
   reset (size_2d);
}

void
Emagram::reset (const Size_2D& size_2d)
{

   this->size_2d = size_2d;


   const Real scale = size_2d.i / 200.0;
   const Real x = ref_thermo_point.get_t ();
   const Real y = log (ref_thermo_point.get_p ());

   const Real anchor_x = (size_2d.i / 33.0);
   const Real anchor_y = size_2d.j - anchor_x * 2;
   set_anchor (Point_2D (anchor_x, anchor_y));

   Affine_Transform_2D::set_identity ();
   Affine_Transform_2D::translate (-x, -y);
   Affine_Transform_2D::scale (scale, magic_ratio * scale);
   Affine_Transform_2D::translate (anchor.x, anchor.y);

}

Skew_T::Skew_T (const Size_2D& size_2d,
                const Real magic_ratio,
                const Real p_0,
                const Thermo_Point& ref_thermo_point)
   : Emagram (size_2d, magic_ratio, p_0, ref_thermo_point)
{
   reset (size_2d);
}

void
Skew_T::reset (const Size_2D& size_2d)
{

   this->size_2d = size_2d;

   const Real scale = size_2d.i / 100.0;
   const Real x = ref_thermo_point.get_t ();
   const Real y = log (ref_thermo_point.get_p ());

   const Real anchor_x = (size_2d.i / 33.0);
   const Real anchor_y = size_2d.j - anchor_x * 2;
   set_anchor (Point_2D (anchor_x, anchor_y));

   Affine_Transform_2D::set_identity ();
   Affine_Transform_2D::translate (-x, -y);
   Affine_Transform_2D::scale (scale, magic_ratio * scale);
   Affine_Transform_2D::shear_x (-1);
   Affine_Transform_2D::translate (anchor.x, anchor.y);

}

Ttxx_Sounding::Ttxx_Sounding (const Integer wmo_id,
                              const Integer yygg)
   : Sounding (wmo_id),
     yygg (yygg)
{
}

const string
Ttxx_Sounding::get_key () const
{
   return string_render ("%05d:%04d", wmo_id, yygg);
}

Thermo_Exception::Thermo_Exception (const string& description)
   : Exception ("Thermo_Exception", description)
{
}

namespace denise
{

   ostream&
   operator << (ostream &out_file,
                const Thermo_Point& thermo_point)
   {
      out_file << "(" << thermo_point.get_t () << ", " <<
         thermo_point.get_theta () << ", " <<
         thermo_point.get_p () << ", " <<
         thermo_point.get_p_0 () << ")";
      return out_file;
   }

}

