//
// nwp.cc
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

#include "marine.h"
#include "nwp.h"

using namespace std;
using namespace denise;

//Level_Tuple::Level_Vector ()
//{
//}

Level_Tuple::Level_Tuple (const Level_Type type,
                          const string& str,
                          const string& delimiter)
   : Tuple (str, delimiter),
     type (type)
{
   add (str, delimiter);
}

void
Level_Tuple::add (const string& str,
                  const string& delimiter,
                  const bool clear_first)
{
   if (clear_first) { clear (); }
   Tuple::add_content (str, delimiter);
}

const Real
Level_Tuple::get_next_up (const Real value) const
{

   const Real v = value;
   typedef Level_Tuple::const_reverse_iterator I;

   switch (type)
   {

      case PRESSURE_LEVEL:
      {
         for (I i = rbegin (); i != rend (); i++) { if (*i < v) { return *i; } }
         break;
      }

      case THETA_LEVEL:
      case SIGMA_LEVEL:
      {
         for (I i = rbegin (); i != rend (); i++) { if (*i > v) { return *i; } }
      }

   }

   throw Nwp_Exception ("Level_Vector::get_next_up confused");

}

const Real
Level_Tuple::get_next_down (const Real value) const
{

   const Real v = value;
   typedef Level_Tuple::const_iterator I;

   switch (type)
   {

      case PRESSURE_LEVEL:
      {
         for (I i = begin (); i != end (); i++) { if (*i > v) { return *i; } }
         break;
      }

      case SIGMA_LEVEL:
      case THETA_LEVEL:
      {
         for (I i = begin (); i != end (); i++) { if (*i < v) { return *i; } }
         break;
      }

   }

   throw Nwp_Exception ("Level_Vector::get_next_up confused");

}

Level_Element::Level_Element (const Level& level,
                              const Nwp_Element nwp_element)
   : level (level),
     nwp_element (nwp_element)
{
}

Nwp::Key::Key (const Dtime& base_time,
               const Integer forecast_hour)
   : base_time (base_time),
     forecast_hour (forecast_hour)
{
}

Nwp::Key::Key (const Key& key)
   : base_time (key.base_time),
     forecast_hour (key.forecast_hour)
{
}

bool
Nwp::Key::operator == (const Key& key) const
{
   const bool a = fabs (base_time.t - key.base_time.t) <= TIME_TOLERANCE;
   const bool b = (forecast_hour == key.forecast_hour);
   return (a && b);
}

bool
Nwp::Key::operator > (const Key& key) const
{
   if (fabs (base_time.t - key.base_time.t) > TIME_TOLERANCE)
   {
      return (base_time > key.base_time);
   }
   else
   {
      return forecast_hour > key.forecast_hour;
   }
}

bool
Nwp::Key::operator < (const Key& key) const
{
   if (fabs (base_time.t - key.base_time.t) > TIME_TOLERANCE)
   {
      return (base_time < key.base_time);
   }
   else
   {
      return forecast_hour < key.forecast_hour;
   }
}

Nwp::Sounding::Sounding (const Key& key)
   : key (key)
{
   time.t = key.base_time.t + key.forecast_hour;
}

Nwp::Element_Vector::Element_Vector (const vector<Nwp_Element>& nwp_element_vector)
   : vector<Nwp_Element> (nwp_element_vector)
{

   for (Integer i = 0; i < nwp_element_vector.size (); i++)
   {
      const Nwp_Element& nwp_element = nwp_element_vector[i];
      reverse_map.insert (make_pair (nwp_element, i));
   }

}

Integer
Nwp::Element_Vector::get_index (const Nwp_Element nwp_element) const
{

   return reverse_map.at (nwp_element);

   for (Integer i = 0; i < size (); i++)
   {
      if (at (i) == nwp_element) { return i; }
   }

   throw Nwp_Exception ("nwp_element not available.");

}

Nwp_Element
Nwp::Element_Vector::get_nwp_element (const Integer element_index) const
{
   return at (element_index);
}

Nwp::Data_3D::Data_3D (const vector<Nwp_Element>& nwp_element_vector,
                       const Key& key)
  : key (key),
    available (false)
{

   typedef vector<Nwp_Element>::const_iterator Iterator;

   for (Iterator iterator = nwp_element_vector.begin ();
        iterator != nwp_element_vector.end (); iterator++)
   {
      const Nwp_Element& nwp_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = NULL;
      insert (make_pair (nwp_element, gvd_3d_ptr));
   }

}

Nwp::Data_3D::~Data_3D ()
{
   for (Nwp::Data_3D::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Geodetic_Vector_Data_3D* gvd_3d_ptr = iterator->second;
      if (gvd_3d_ptr != NULL) { delete gvd_3d_ptr; }
   }
}

void
Nwp::Data_3D::set_gvd_3d_ptr (const Nwp_Element nwp_element,
                              Geodetic_Vector_Data_3D* gvd_3d_ptr)
{
   at (nwp_element) = gvd_3d_ptr;
}

const Geodetic_Vector_Data_3D&
Nwp::Data_3D::get_gvd_3d (const Nwp_Element nwp_element) const
{
   try
   {
      const Geodetic_Vector_Data_3D* gvd_3d_ptr = at (nwp_element);
      if (gvd_3d_ptr == NULL)
      {
         throw Nwp_Exception ("Nwp::Data_3D gvd_3d_ptr is NULL");
      }
      return *gvd_3d_ptr;
   }
   catch (const std::exception& se)
   {
      throw Nwp_Exception ("Nwp::Data_3D has no such nwp_element");
   }
}

Geodetic_Vector_Data_3D&
Nwp::Data_3D::get_gvd_3d (const Nwp_Element nwp_element)
{
   Geodetic_Vector_Data_3D* gvd_3d_ptr = at (nwp_element);
   try
   {
      Geodetic_Vector_Data_3D* gvd_3d_ptr = at (nwp_element);
      if (gvd_3d_ptr == NULL)
      {
         throw Nwp_Exception ("Nwp::Data_3D gvd_3d_ptr is NULL");
      }
      return *gvd_3d_ptr;
   }
   catch (const std::exception& se)
   {
      throw Nwp_Exception ("Nwp::Data_3D has no such nwp_element");
   }
}

const Tuple&
Nwp::Data_3D::get_tuple_p (const Nwp_Element nwp_element) const
{
   try
   {
      const Geodetic_Vector_Data_3D& gvd_3d = get_gvd_3d (nwp_element);
      return gvd_3d.get_coordinate_tuple (0);
   }
   catch (const Nwp_Exception& ne)
   {
      const Geodetic_Vector_Data_3D& gvd_3d = *(begin ()->second);
      return gvd_3d.get_coordinate_tuple (0);
   }
}

Real
Nwp::Data_3D::get_p (const Nwp_Element nwp_element,
                     const Integer k) const
{
   const Tuple& tuple_p = get_tuple_p (nwp_element);
   return tuple_p[k];
}

Lat_Long
Nwp::Data_3D::get_lat_long (const Nwp_Element nwp_element,
                            const Integer i,
                            const Integer j) const
{
   const Geodetic_Vector_Data_3D& gvd_3d = get_gvd_3d (nwp_element);
   const Real latitude = gvd_3d.get_coordinate (1, i);
   const Real longitude = gvd_3d.get_coordinate (2, j);
   return Lat_Long (latitude, longitude);
}

bool
Nwp::Data_3D::is_available () const
{
   return available;
}

void
Nwp::Data_3D::set_available ()
{
   available = true;
}

void
Nwp::Data_3D::initialize (const Real datum)
{
   for (Nwp::Data_3D::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Geodetic_Vector_Data_3D* gvd_3d_ptr = iterator->second;
      if (gvd_3d_ptr == NULL) { continue; }
      gvd_3d_ptr->initialize_all (datum);
   }
}

Real
Nwp::Data_3D::get_p_from_element (const Nwp_Element nwp_element,
                                  const Integer i,
                                  const Integer j,
                                  const Real element_value) const
{

   const Geodetic_Vector_Data_3D* gvd_3d_ptr = at (nwp_element);
   if (gvd_3d_ptr == NULL) { return GSL_NAN; }
   const Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   Real p = GSL_NAN;
   const Real x = element_value;
   const Tuple& tuple_p = get_tuple_p (nwp_element);

   for (Integer k = 0; k < tuple_p.size () - 1; k++)
   {

      const Real lower_x = gvd_3d.get_datum (0, k, i, j); 
      const Real upper_x = gvd_3d.get_datum (0, k+1, i, j); 
      const Real delta_x = x - lower_x;
      const bool match = (delta_x * (x - upper_x) <= 0);

      if (match)
      {
         const Real lower_p = tuple_p[k];
         const Real upper_p = tuple_p[k + 1];
         const Real dp = upper_p - lower_p;
         const Real dx = upper_x - lower_x;
         p = lower_p + (delta_x  * dp / dx);
         break;
      }

   }

   return p;

}

Real
Nwp::Data_3D::get_p_from_element (const Nwp_Element nwp_element,
                                  const Real latitude,
                                  const Real longitude,
                                  const Real element_value) const
{

   const Geodetic_Vector_Data_3D* gvd_3d_ptr = at (nwp_element);
   if (gvd_3d_ptr == NULL) { return GSL_NAN; }
   const Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   Real p = GSL_NAN;
   const Real x = element_value;
   const Tuple& tuple_p = get_tuple_p (nwp_element);

   for (Integer k = 0; k < tuple_p.size () - 1; k++)
   {

      const Real lower_x = evaluate (nwp_element, k, latitude, longitude); 
      const Real upper_x = evaluate (nwp_element, k+1, latitude, longitude); 
      const Real delta_x = x - lower_x;
      const bool match = (delta_x * (x - upper_x) <= 0);

      if (match)
      {
         const Real lower_p = tuple_p[k];
         const Real upper_p = tuple_p[k + 1];
         const Real dp = upper_p - lower_p;
         const Real dx = upper_x - lower_x;
         p = lower_p + (delta_x  * dp / dx);
         break;
      }

   }

   return p;

}


Real
Nwp::Data_3D::evaluate (const Nwp_Element element,
                        const Integer k,
                        const Integer i,
                        const Integer j,
                        const Evaluate_Op evaluate_op) const
{

   switch (element)
   {

      case WIND_SPEED:
      {
         const Real u = evaluate (ZONAL_WIND, k, i, j);
         const Real v = evaluate (MERIDIONAL_WIND, k, i, j);
         return sqrt (u*u + v*v);
      }

      case WIND_DIRECTION:
      {
         const Real u = evaluate (ZONAL_WIND, k, i, j);
         const Real v = evaluate (MERIDIONAL_WIND, k, i, j);
         return Wind (u, v).get_direction ();
      }

      case DEW_POINT_DEPRESSION:
      {
         const Real t = evaluate (TEMPERATURE, k, i, j);
         const Real t_d = evaluate (DEW_POINT, k, i, j);
         return t - t_d;
      }

      case THETA:
      {
         const Real p = get_p (TEMPERATURE, k);
         const Real t = evaluate (TEMPERATURE, k, i, j);
         return Thermo_Point::t_p (t - K, p).get_theta () + K;
      }

      case THETA_E:
      {
         const Real p = get_p (TEMPERATURE, k);
         const Real t = evaluate (TEMPERATURE, k, i, j);
         const Real t_d = evaluate (DEW_POINT, k, i, j);
         return Thermo_Point::normand (t - K, t_d - K, p).get_theta_e () + K;
      }

      case THETA_W:
      {
         const Real p = get_p (TEMPERATURE, k);
         const Real t = evaluate (TEMPERATURE, k, i, j);
         const Real t_d = evaluate (DEW_POINT, k, i, j);
         return Thermo_Point::normand (t - K, t_d - K, p).get_theta_w () + K;
      }

      case TEMPERATURE_ADVECTION:
      {
         const Real t_x = evaluate (TEMPERATURE, k, i, j, DX);
         const Real t_y = evaluate (TEMPERATURE, k, i, j, DY);
         const Real u = evaluate (ZONAL_WIND, k, i, j);
         const Real v = evaluate (MERIDIONAL_WIND, k, i, j);
         return -(t_x * u + t_y * v);
      }

      case ADIABATIC_HEATING:
      {
         const Real p = get_p (TEMPERATURE, k);
         const Real t = evaluate (TEMPERATURE, k, i, j);
         const Real omega = evaluate (OMEGA, k, i, j);
         const Real alpha = (R_d * t) / p;
         return alpha / c_p * omega;
      }

      case LATENT_HEATING:
      {
         const Real rh = evaluate (RELATIVE_HUMIDITY, k, i, j);
         if (rh < 0.9) { return 0; }
         const Real p = get_p (TEMPERATURE, k);
         const Real t = evaluate (TEMPERATURE, k, i, j);
         const Real omega = evaluate (OMEGA, k, i, j);
         const Thermo_Point& tp = Thermo_Point::t_p (t, p);
         return tp.get_saturated_Q_dot (omega) * rh;
      }

      case MONTGOMERY:
      {
         const Real t = evaluate (TEMPERATURE, k, i, j);
         const Real z = evaluate (GEOPOTENTIAL_HEIGHT, k, i, j);
         return g * z + c_p * t;
      }

      case ABSOLUTE_VORTICITY:
      {

         const Nwp_Element U = ZONAL_WIND;
         const Nwp_Element V = MERIDIONAL_WIND;

         const Lat_Long& lat_long = get_lat_long (U, i, j);
         const Real& latitude = lat_long.latitude;
         const Real& longitude = lat_long.longitude;
         const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

         const Real dv_dx = evaluate (V, k, latitude, longitude, DX);
         const Real du_dy = evaluate (U, k, latitude, longitude, DY);

         const Real zeta = dv_dx - du_dy + f;
         return zeta;

      };

      case POTENTIAL_VORTICITY:
      {

         const Real p = get_p (TEMPERATURE, k);
         const Real exner = pow (Real (p / 1000e2), Real (kappa));

         const Lat_Long& lat_long = get_lat_long (TEMPERATURE, i, j);
         const Real& latitude = lat_long.latitude;
         const Real& longitude = lat_long.longitude;
         const Real f = Geodetic_Vector_Data_2D::get_f (latitude);
         const Real f_hat = Geodetic_Vector_Data_2D::get_f_hat (latitude);

         const Nwp_Element T = TEMPERATURE;
         const Nwp_Element U = ZONAL_WIND;
         const Nwp_Element V = MERIDIONAL_WIND;
         const Nwp_Element W = VERTICAL_VELOCITY;
         const Nwp_Element Z = GEOPOTENTIAL_HEIGHT;

         const Real t = evaluate (T, k, i, j);
         const Real dt_dx = evaluate (T, k, latitude, longitude, DX);
         const Real dt_dy = evaluate (T, k, latitude, longitude, DY);
         const Real dt_dp = evaluate (T, k, latitude, longitude, DZ);
         const Real dv_dx = evaluate (V, k, latitude, longitude, DX);
         const Real dv_dp = evaluate (V, k, latitude, longitude, DZ);
         const Real du_dy = evaluate (U, k, latitude, longitude, DY);
         const Real du_dp = evaluate (U, k, latitude, longitude, DZ);
         const Real dz_dp = evaluate (Z, k, latitude, longitude, DZ);

         const Real rho = p / (R_d * t);
         const Real rho_g = rho * g;

         const Real dw_dx = evaluate (W, k, latitude, longitude, DX);
         const Real dw_dy = evaluate (W, k, latitude, longitude, DY);

         const Real dt_dz = (dt_dp - (kappa * t / p)) / dz_dp;

         const Real xi = dw_dy - dv_dp / dz_dp;
         const Real eta = du_dp / dz_dp - dw_dx + f_hat;
         const Real zeta = dv_dx - du_dy + f;

         const Real pv_x = xi * dt_dx;
         const Real pv_y = eta * dt_dy;
         const Real pv_z = zeta * dt_dz;
         const Real pv = (pv_x + pv_y + pv_z) / (exner * rho);

         return pv;

      };

   }

}

Real
Nwp::Data_3D::evaluate (const Nwp_Element element,
                        const Real p,
                        const Integer i,
                        const Integer j,
                        const Evaluate_Op evaluate_op) const
{

   switch (element)
   {

      case WIND_SPEED:
      {
         const Real u = evaluate (ZONAL_WIND, p, i, j);
         const Real v = evaluate (MERIDIONAL_WIND, p, i, j);
         return sqrt (u*u + v*v);
      }

      case WIND_DIRECTION:
      {
         const Real u = evaluate (ZONAL_WIND, p, i, j);
         const Real v = evaluate (MERIDIONAL_WIND, p, i, j);
         return Wind (u, v).get_direction ();
      }

      case DEW_POINT_DEPRESSION:
      {
         const Real t = evaluate (TEMPERATURE, p, i, j);
         const Real t_d = evaluate (DEW_POINT, p, i, j);
         return t - t_d;
      }

      case THETA:
      {
         const Real t = evaluate (TEMPERATURE, p, i, j);
         return Thermo_Point::t_p (t - K, p).get_theta () + K;
      }

      case THETA_E:
      {
         const Real t = evaluate (TEMPERATURE, p, i, j);
         const Real t_d = evaluate (DEW_POINT, p, i, j);
         return Thermo_Point::normand (t - K, t_d - K, p).get_theta_e () + K;
      }

      case THETA_W:
      {
         const Real t = evaluate (TEMPERATURE, p, i, j);
         const Real t_d = evaluate (DEW_POINT, p, i, j);
         return Thermo_Point::normand (t - K, t_d - K, p).get_theta_w () + K;
      }

      case TEMPERATURE_ADVECTION:
      {
         const Real t_x = evaluate (TEMPERATURE, p, i, j, DX);
         const Real t_y = evaluate (TEMPERATURE, p, i, j, DY);
         const Real u = evaluate (ZONAL_WIND, p, i, j);
         const Real v = evaluate (MERIDIONAL_WIND, p, i, j);
         return -(t_x * u + t_y * v);
      }

      case ADIABATIC_HEATING:
      {
         const Real t = evaluate (TEMPERATURE, p, i, j);
         const Real omega = evaluate (OMEGA, p, i, j);
         const Real alpha = (R_d * t) / p;
         return alpha / c_p * omega;
      }

      case LATENT_HEATING:
      {
         const Real rh = evaluate (RELATIVE_HUMIDITY, p, i, j);
         if (rh < 0.9) { return 0; }
         const Real t = evaluate (TEMPERATURE, p, i, j);
         const Real omega = evaluate (OMEGA, p, i, j);
         const Thermo_Point& tp = Thermo_Point::t_p (t, p);
         return tp.get_saturated_Q_dot (omega) * rh;
      }

      case MONTGOMERY:
      {
         const Real t = evaluate (TEMPERATURE, p, i, j, evaluate_op);
         const Real z = evaluate (GEOPOTENTIAL_HEIGHT, p, i, j, evaluate_op);
         return g * z + c_p * t;
      }

      case ABSOLUTE_VORTICITY:
      {

         const Lat_Long& lat_long = get_lat_long (MERIDIONAL_WIND, i, j);
         const Real& latitude = lat_long.latitude;
         const Real& longitude = lat_long.longitude;
         const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

         const Real dv_dx = evaluate (MERIDIONAL_WIND, p, i, j, DX);
         const Real du_dy = evaluate (ZONAL_WIND, p, i, j, DY);

         const Real zeta = dv_dx - du_dy + f;
         return zeta;

      }

      case POTENTIAL_VORTICITY:
      {

         const Real exner = pow (Real (p / 1000e2), Real (kappa));
         const Lat_Long& lat_long = get_lat_long (TEMPERATURE, i, j);
         const Real& latitude = lat_long.latitude;
         const Real& longitude = lat_long.longitude;
         const Real f = Geodetic_Vector_Data_2D::get_f (latitude);
         const Real f_hat = Geodetic_Vector_Data_2D::get_f_hat (latitude);

         const Real t = evaluate (TEMPERATURE, p, i, j);
         const Real dt_dx = evaluate (TEMPERATURE, p, i, j, DX);
         const Real dt_dy = evaluate (TEMPERATURE, p, i, j, DY);
         const Real dt_dp = evaluate (TEMPERATURE, p, i, j, DZ);
         const Real dv_dx = evaluate (MERIDIONAL_WIND, p, i, j, DX);
         const Real dv_dp = evaluate (MERIDIONAL_WIND, p, i, j, DZ);
         const Real du_dy = evaluate (ZONAL_WIND, p, i, j, DY);
         const Real du_dp = evaluate (ZONAL_WIND, p, i, j, DZ);
         const Real dz_dp = evaluate (GEOPOTENTIAL_HEIGHT, p, i, j, DZ);
         const Real dw_dx = evaluate (VERTICAL_VELOCITY, p, i, j, DX);
         const Real dw_dy = evaluate (VERTICAL_VELOCITY, p, i, j, DY);

         const Real rho = p / (R_d * t);
         const Real rho_g = rho * g;
         const Real dt_dz = (dt_dp - (kappa * t / p)) / dz_dp;

         const Real xi = dw_dy - dv_dp / dz_dp;
         const Real eta = du_dp / dz_dp - dw_dx + f_hat;
         const Real zeta = dv_dx - du_dy + f;

         const Real pv_x = xi * dt_dx;
         const Real pv_y = eta * dt_dy;
         const Real pv_z = zeta * dt_dz;
         const Real pv = (pv_x + pv_y + pv_z) / (exner * rho);

         return pv;

      };

   }

   Data_3D::const_iterator iterator = find (element);
   if (iterator == end ()) { return GSL_NAN; }

   const Geodetic_Vector_Data_3D* gvd_3d_ptr = iterator->second;
   if (gvd_3d_ptr == NULL) { return GSL_NAN; }

   const Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;
   return gvd_3d.evaluate (0, p, i, j, evaluate_op);
}

Real
Nwp::Data_3D::evaluate (const Nwp_Element element,
                        const Integer k,
                        const Real latitude,
                        const Real longitude,
                        const Evaluate_Op evaluate_op) const
{

   switch (element)
   {

      case WIND_SPEED:
      {
         const Real u = evaluate (ZONAL_WIND, k, latitude, longitude);
         const Real v = evaluate (MERIDIONAL_WIND, k, latitude, longitude);
         return sqrt (u*u + v*v);
      }

      case WIND_DIRECTION:
      {
         const Real u = evaluate (ZONAL_WIND, k, latitude, longitude);
         const Real v = evaluate (MERIDIONAL_WIND, k, latitude, longitude);
         return Wind (u, v).get_direction ();
      }

      case DEW_POINT_DEPRESSION:
      {
         const Real t = evaluate (TEMPERATURE, k, latitude, longitude);
         const Real t_d = evaluate (DEW_POINT, k, latitude, longitude);
         return t - t_d;
      }

      case THETA:
      {
         const Real p = get_p (TEMPERATURE, k);
         const Real t = evaluate (TEMPERATURE, k, latitude, longitude);
         return Thermo_Point::t_p (t - K, p).get_theta () + K;;
      }

      case THETA_E:
      {
         const Real p = get_p (TEMPERATURE, k);
         const Real t = evaluate (TEMPERATURE, k, latitude, longitude);
         const Real t_d = evaluate (DEW_POINT, k, latitude, longitude);
         return Thermo_Point::normand (t - K, t_d - K, p).get_theta_e () + K;
      }

      case THETA_W:
      {
         const Real p = get_p (TEMPERATURE, k);
         const Real t = evaluate (TEMPERATURE, k, latitude, longitude);
         const Real t_d = evaluate (DEW_POINT, k, latitude, longitude);
         return Thermo_Point::normand (t - K, t_d - K, p).get_theta_w () + K;
      }

      case TEMPERATURE_ADVECTION:
      {
         const Real t_x = evaluate (TEMPERATURE, k, latitude, longitude, DX);
         const Real t_y = evaluate (TEMPERATURE, k, latitude, longitude, DY);
         const Real u = evaluate (ZONAL_WIND, k, latitude, longitude);
         const Real v = evaluate (MERIDIONAL_WIND, k, latitude, longitude);
         return -(t_x * u + t_y * v);
      }

      case ADIABATIC_HEATING:
      {
         const Real p = get_p (TEMPERATURE, k);
         const Real t = evaluate (TEMPERATURE, k, latitude, longitude);
         const Real omega = evaluate (OMEGA, k, latitude, longitude);
         const Real alpha = (R_d * t) / p;
         return alpha / c_p * omega;
      }

      case LATENT_HEATING:
      {
         const Real rh = evaluate (RELATIVE_HUMIDITY, k, latitude, longitude);
         if (rh < 0.9) { return 0; }
         const Real p = get_p (TEMPERATURE, k);
         const Real t = evaluate (TEMPERATURE, k, latitude, longitude);
         const Real omega = evaluate (OMEGA, k, latitude, longitude);
         const Thermo_Point& tp = Thermo_Point::t_p (t, p);
         return tp.get_saturated_Q_dot (omega) * rh;
      }

      case MONTGOMERY:
      {
         const Nwp_Element T = TEMPERATURE;
         const Nwp_Element Z = GEOPOTENTIAL_HEIGHT;
         const Real t = evaluate (T, k, latitude, longitude, evaluate_op);
         const Real z = evaluate (Z, k, latitude, longitude, evaluate_op);
         return g * z + c_p * t;
      }

      case ABSOLUTE_VORTICITY:
      {

         const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

         const Nwp_Element U = ZONAL_WIND;
         const Nwp_Element V = MERIDIONAL_WIND;

         const Real dv_dx = evaluate (V, k, latitude, longitude, DX);
         const Real du_dy = evaluate (U, k, latitude, longitude, DY);
         const Real zeta = dv_dx - du_dy + f;
         return zeta;

      }

      case POTENTIAL_VORTICITY:
      {

         const Real p = get_p (TEMPERATURE, k);
         const Real exner = pow (Real (p / 1000e2), Real (kappa));
         const Real f = Geodetic_Vector_Data_2D::get_f (latitude);
         const Real f_hat = Geodetic_Vector_Data_2D::get_f_hat (latitude);

         const Nwp_Element T = TEMPERATURE;
         const Nwp_Element U = ZONAL_WIND;
         const Nwp_Element V = MERIDIONAL_WIND;
         const Nwp_Element W = VERTICAL_VELOCITY;
         const Nwp_Element Z = GEOPOTENTIAL_HEIGHT;

         const Real t = evaluate (T, k, latitude, longitude);
         const Real dt_dx = evaluate (T, k, latitude, longitude, DX);
         const Real dt_dy = evaluate (T, k, latitude, longitude, DY);
         const Real dt_dp = evaluate (T, k, latitude, longitude, DZ);
         const Real dv_dx = evaluate (V, k, latitude, longitude, DX);
         const Real dv_dp = evaluate (V, k, latitude, longitude, DZ);
         const Real du_dy = evaluate (U, k, latitude, longitude, DY);
         const Real du_dp = evaluate (U, k, latitude, longitude, DZ);
         const Real dz_dp = evaluate (Z, k, latitude, longitude, DZ);
         const Real dw_dx = evaluate (W, k, latitude, longitude, DX);
         const Real dw_dy = evaluate (W, k, latitude, longitude, DY);

         const Real rho = p / (R_d * t);
         const Real rho_g = rho * g;
         const Real dt_dz = (dt_dp - (kappa * t / p)) / dz_dp;

         const Real xi = dw_dy - dv_dp / dz_dp;
         const Real eta = du_dp / dz_dp - dw_dx + f_hat;
         const Real zeta = dv_dx - du_dy + f;

         const Real pv_x = xi * dt_dx;
         const Real pv_y = eta * dt_dy;
         const Real pv_z = zeta * dt_dz;
         const Real pv = (pv_x + pv_y + pv_z) / (exner * rho);

         return pv;

      }

   }

   Data_3D::const_iterator iterator = find (element);
   if (iterator == end ()) { return GSL_NAN; }

   const Geodetic_Vector_Data_3D* gvd_3d_ptr = iterator->second;
   if (gvd_3d_ptr == NULL) { return GSL_NAN; }

   const Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;
   return gvd_3d.evaluate (0, k, latitude, longitude, evaluate_op);

}

Real
Nwp::Data_3D::evaluate (const Nwp_Element element,
                        const Integer k,
                        const Lat_Long& lat_long,
                        const Evaluate_Op evaluate_op) const
{
   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;
   return evaluate (element, k, latitude, longitude, evaluate_op);
}

Real
Nwp::Data_3D::evaluate (const Nwp_Element element,
                        const Real p,
                        const Real latitude,
                        const Real longitude,
                        const Evaluate_Op evaluate_op) const
{

   switch (element)
   {

      case WIND_SPEED:
      {
         const Real u = evaluate (ZONAL_WIND, p, latitude, longitude);
         const Real v = evaluate (MERIDIONAL_WIND, p, latitude, longitude);
         return sqrt (u*u + v*v);
      }

      case WIND_DIRECTION:
      {
         const Real u = evaluate (ZONAL_WIND, p, latitude, longitude);
         const Real v = evaluate (MERIDIONAL_WIND, p, latitude, longitude);
         return Wind (u, v).get_direction ();
      }

      case DEW_POINT_DEPRESSION:
      {
         const Evaluate_Op& eo = evaluate_op;
         const Real t = evaluate (TEMPERATURE, p, latitude, longitude);
         const Real t_d = evaluate (DEW_POINT, p, latitude, longitude);
         return t - t_d;
      }

      case THETA:
      {
         const Evaluate_Op& eo = evaluate_op;
         const Real t = evaluate (TEMPERATURE, p, latitude, longitude);
         return Thermo_Point::t_p (t - K, p).get_theta () + K;
      }

      case THETA_E:
      {
         const Evaluate_Op& eo = evaluate_op;
         const Real t = evaluate (TEMPERATURE, p, latitude, longitude);
         const Real t_d = evaluate (DEW_POINT, p, latitude, longitude);
         return Thermo_Point::normand (t - K, t_d - K, p).get_theta_e () + K;
      }

      case THETA_W:
      {
         const Evaluate_Op& eo = evaluate_op;
         const Real t = evaluate (TEMPERATURE, p, latitude, longitude);
         const Real t_d = evaluate (DEW_POINT, p, latitude, longitude);
         return Thermo_Point::normand (t - K, t_d - K, p).get_theta_w () + K;
      }

      case TEMPERATURE_ADVECTION:
      {
         const Real t_x = evaluate (TEMPERATURE, p, latitude, longitude, DX);
         const Real t_y = evaluate (TEMPERATURE, p, latitude, longitude, DY);
         const Real u = evaluate (ZONAL_WIND, p, latitude, longitude);
         const Real v = evaluate (MERIDIONAL_WIND, p, latitude, longitude);
         return -(t_x * u + t_y * v);
      }

      case ADIABATIC_HEATING:
      {
         const Real t = evaluate (TEMPERATURE, p, latitude, longitude);
         const Real omega = evaluate (OMEGA, p, latitude, longitude);
         const Real alpha = (R_d * t) / p;
         return alpha / c_p * omega;
      }

      case LATENT_HEATING:
      {
         const Real rh = evaluate (RELATIVE_HUMIDITY, p, latitude, longitude);
         if (rh < 0.9) { return 0; }
         const Real t = evaluate (TEMPERATURE, p, latitude, longitude);
         const Real omega = evaluate (OMEGA, p, latitude, longitude);
         const Thermo_Point& tp = Thermo_Point::t_p (t, p);
         return tp.get_saturated_Q_dot (omega) * rh;
      }

      case MONTGOMERY:
      {
         const Nwp_Element T = TEMPERATURE;
         const Nwp_Element Z = GEOPOTENTIAL_HEIGHT;
         const Real t = evaluate (T, p, latitude, longitude);
         const Real z = evaluate (Z, p, latitude, longitude);
         return g * z + c_p * t;
      }

      case ABSOLUTE_VORTICITY:
      {

         const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

         const Nwp_Element U = ZONAL_WIND;
         const Nwp_Element V = MERIDIONAL_WIND;

         const Real dv_dx = evaluate (V, p, latitude, longitude, DX);
         const Real du_dy = evaluate (U, p, latitude, longitude, DY);
         const Real zeta = dv_dx - du_dy + f;
         return zeta;

      }

      case POTENTIAL_VORTICITY:
      {

         const Real exner = pow (Real (p / 1000e2), Real (kappa));
         const Real f = Geodetic_Vector_Data_2D::get_f (latitude);
         const Real f_hat = Geodetic_Vector_Data_2D::get_f_hat (latitude);

         const Nwp_Element T = TEMPERATURE;
         const Nwp_Element U = ZONAL_WIND;
         const Nwp_Element V = MERIDIONAL_WIND;
         const Nwp_Element W = VERTICAL_VELOCITY;
         const Nwp_Element Z = GEOPOTENTIAL_HEIGHT;

         const Real t = evaluate (T, p, latitude, longitude);
         const Real dt_dx = evaluate (T, p, latitude, longitude, DX);
         const Real dt_dy = evaluate (T, p, latitude, longitude, DY);
         const Real dt_dp = evaluate (T, p, latitude, longitude, DZ);
         const Real dv_dx = evaluate (V, p, latitude, longitude, DX);
         const Real dv_dp = evaluate (V, p, latitude, longitude, DZ);
         const Real du_dy = evaluate (U, p, latitude, longitude, DY);
         const Real du_dp = evaluate (U, p, latitude, longitude, DZ);
         const Real dz_dp = evaluate (Z, p, latitude, longitude, DZ);

         const Real rho = p / (R_d * t);
         const Real rho_g = rho * g;

         const Real dw_dx = evaluate (W, p, latitude, longitude, DX);
         const Real dw_dy = evaluate (W, p, latitude, longitude, DY);

         const Real dt_dz = (dt_dp - (kappa * t / p)) / dz_dp;

         const Real xi = dw_dy - dv_dp / dz_dp;
         const Real eta = du_dp / dz_dp - dw_dx + f_hat;
         const Real zeta = dv_dx - du_dy + f;

         const Real pv_x = xi * dt_dx;
         const Real pv_y = eta * dt_dy;
         const Real pv_z = zeta * dt_dz;
         const Real pv = (pv_x + pv_y + pv_z) / (exner * rho);

         return pv;

      }

   }

   
   Data_3D::const_iterator iterator = find (element);
   if (iterator == end ()) { return GSL_NAN; }

   const Geodetic_Vector_Data_3D* gvd_3d_ptr = iterator->second;
   if (gvd_3d_ptr == NULL) { return GSL_NAN; }

   const Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;
   return gvd_3d.evaluate (0, p, latitude, longitude, evaluate_op);

}

Real
Nwp::Data_3D::evaluate (const Nwp_Element element,
                        const Real p,
                        const Lat_Long& lat_long,
                        const Evaluate_Op evaluate_op) const
{
   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;
   return evaluate (element, p, latitude, longitude, evaluate_op);
}

Real
Nwp::Data_3D::get_li_thunder (const Integer k,
                              const Integer i,
                              const Integer j,
                              const Real thunder_p,
                              const Real thunder_t) const
{
   const Real p = get_p (TEMPERATURE, k);
   if (p < thunder_p) { return GSL_NAN; }
   const Real t = evaluate (TEMPERATURE, k, i, j);
   const Real td = evaluate (DEW_POINT, k, i, j);
   return Instability::get_lifted_index (p, t, td, thunder_p, thunder_t-K) + K;
}

Real
Nwp::Data_3D::get_li_thunder (const Real p,
                              const Integer i,
                              const Integer j,
                              const Real thunder_p,
                              const Real thunder_t) const
{
   if (p < thunder_p) { return GSL_NAN; }
   const Real t = evaluate (TEMPERATURE, p, i, j);
   const Real td = evaluate (DEW_POINT, p, i, j);
   return Instability::get_lifted_index (p, t, td, thunder_p, thunder_t-K) + K;
}

Real
Nwp::Data_3D::get_li_thunder (const Integer k,
                              const Real latitude,
                              const Real longitude,
                              const Real thunder_p,
                              const Real thunder_t) const
{
   const Real p = get_p (TEMPERATURE, k);
   if (p < thunder_p) { return GSL_NAN; }
   const Real t = evaluate (TEMPERATURE, k, latitude, longitude);
   const Real td = evaluate (DEW_POINT, k, latitude, longitude);
   return Instability::get_lifted_index (p, t, td, thunder_p, thunder_t-K) + K;
}

Real
Nwp::Data_3D::get_li_thunder (const Integer k,
                              const Lat_Long& lat_long,
                              const Real thunder_p,
                              const Real thunder_t) const
{
   const Real p = get_p (TEMPERATURE, k);
   if (p < thunder_p) { return GSL_NAN; }
   const Real t = evaluate (TEMPERATURE, k, lat_long);
   const Real td = evaluate (DEW_POINT, k, lat_long);
   return Instability::get_lifted_index (p, t, td, thunder_p, thunder_t-K) + K;
}

Real
Nwp::Data_3D::get_li_thunder (const Real p,
                              const Real latitude,
                              const Real longitude,
                              const Real thunder_p,
                              const Real thunder_t) const
{
   if (p < thunder_p) { return GSL_NAN; }
   const Real t = evaluate (TEMPERATURE, p, latitude, longitude);
   const Real td = evaluate (DEW_POINT, p, latitude, longitude);
   return Instability::get_lifted_index (p, t, td, thunder_p, thunder_t-K) + K;
}

Real
Nwp::Data_3D::get_li_thunder (const Real p,
                              const Lat_Long& lat_long,
                              const Real thunder_p,
                              const Real thunder_t) const
{
   if (p < thunder_p) { return GSL_NAN; }
   const Real t = evaluate (TEMPERATURE, p, lat_long);
   const Real td = evaluate (DEW_POINT, p, lat_long);
   return Instability::get_lifted_index (p, t, td, thunder_p, thunder_t-K) + K;
}

Nwp::Cross_Section::Cross_Section ()
   : terrain_profile_ptr (NULL),
     rainfall_profile_ptr (NULL)
{
}

Nwp::Cross_Section::~Cross_Section ()
{
   for (Nwp::Cross_Section::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Scalar_Data_2D* sd_2d_ptr = iterator->second;
      if (sd_2d_ptr != NULL) { delete sd_2d_ptr; }
   }
   if (terrain_profile_ptr != NULL) { delete terrain_profile_ptr; }
   if (rainfall_profile_ptr != NULL) { delete rainfall_profile_ptr; }
}

void
Nwp::Cross_Section::set_profile_ptrs (Scalar_Data_1D* terrain_profile_ptr,
                                      Scalar_Data_1D* rainfall_profile_ptr)
{
   this->terrain_profile_ptr = terrain_profile_ptr;
   this->rainfall_profile_ptr = rainfall_profile_ptr;
}

void
Nwp::Cross_Section::insert_nwp_element_if_needed (const Nwp_Element nwp_element,
                                                  const Tuple& tuple_x,
                                                  const Tuple& tuple_p)
{
   if (find (nwp_element) == end ())
   {
      Scalar_Data_2D* sd_2d_ptr = new Scalar_Data_2D (tuple_x, tuple_p);
      insert (make_pair (nwp_element, sd_2d_ptr));
   }
}

const Scalar_Data_1D&
Nwp::Cross_Section::get_terrain_profile () const
{
   return *terrain_profile_ptr;
}

const Scalar_Data_1D&
Nwp::Cross_Section::get_rainfall_profile () const
{
   return *rainfall_profile_ptr;
}

const Scalar_Data_2D&
Nwp::Cross_Section::get_sd_2d (const Nwp_Element nwp_element) const
{
   return *at (nwp_element);
}

Scalar_Data_2D&
Nwp::Cross_Section::get_sd_2d (const Nwp_Element nwp_element)
{
   return *at (nwp_element);
}

const Tuple&
Nwp::Cross_Section::get_tuple_p (const Nwp_Element nwp_element) const
{
   return at (nwp_element)->get_coordinate_tuple (1);
}

Nwp::Time_Cross::Time_Cross ()
   : terrain_profile_ptr (NULL),
     rainfall_profile_ptr (NULL)
{
}

Nwp::Time_Cross::~Time_Cross ()
{
   for (Nwp::Cross_Section::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Scalar_Data_2D* sd_2d_ptr = iterator->second;
      if (sd_2d_ptr != NULL) { delete sd_2d_ptr; }
   }
   if (terrain_profile_ptr != NULL) { delete terrain_profile_ptr; }
   if (rainfall_profile_ptr != NULL) { delete rainfall_profile_ptr; }
}

void
Nwp::Time_Cross::set_profile_ptrs (Scalar_Data_1D* terrain_profile_ptr,
                                   Scalar_Data_1D* rainfall_profile_ptr)
{
   this->terrain_profile_ptr = terrain_profile_ptr;
   this->rainfall_profile_ptr = rainfall_profile_ptr;
}

void
Nwp::Time_Cross::insert_nwp_element_if_needed (const Nwp_Element nwp_element,
                                               const Tuple& tuple_t,
                                               const Tuple& tuple_p)
{
   if (find (nwp_element) == end ())
   {
      Scalar_Data_2D* sd_2d_ptr = new Scalar_Data_2D (tuple_t, tuple_p);
      insert (make_pair (nwp_element, sd_2d_ptr));
   }
}

const Scalar_Data_1D&
Nwp::Time_Cross::get_terrain_profile () const
{
   return *terrain_profile_ptr;
}

const Scalar_Data_1D&
Nwp::Time_Cross::get_rainfall_profile () const
{
   return *rainfall_profile_ptr;
}

const Scalar_Data_2D&
Nwp::Time_Cross::get_sd_2d (const Nwp_Element nwp_element) const
{
   return *at (nwp_element);
}

Scalar_Data_2D&
Nwp::Time_Cross::get_sd_2d (const Nwp_Element nwp_element)
{
   return *at (nwp_element);
}

const Tuple&
Nwp::Time_Cross::get_tuple_p (const Nwp_Element nwp_element) const
{
   return at (nwp_element)->get_coordinate_tuple (1);
}

void
Nwp::Key_Multimap::clear ()
{
   multimap<Dtime, Integer>::clear ();
   key_set.clear ();
}

void
Nwp::Key_Multimap::add (const Key& key)
{
   if (key_set.find (key) == key_set.end ())
   {
      key_set.insert (key);
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      insert (make_pair (bt, fh));
   }
}

bool
Nwp::Key_Multimap::is_no_match (const Key& key) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   try
   {

      typedef Key_Multimap::const_iterator Iterator;
      pair<Iterator, Iterator> range = equal_range (base_time);

      for (Iterator iterator = range.first;
           iterator != range.second; iterator++)
      {
         const Integer fh = iterator->second;
         if (forecast_hour == fh) { return false; }
      }

   }
   catch (const std::exception& se)
   {
   }

   return true;

}

bool
Nwp::Key_Multimap::is_first_step (const Key& key) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   try
   {

      typedef Key_Multimap::const_iterator Iterator;
      pair<Iterator, Iterator> range = equal_range (base_time);

      const Integer first_fh = range.first->second;
      return (forecast_hour == first_fh);

   }
   catch (const std::exception& se)
   {
   }

   return false;

}

Nwp::Key
Nwp::Key_Multimap::get_previous_key (const Key& key) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   typedef Key_Multimap::const_iterator Iterator;
   pair<Iterator, Iterator> range = equal_range (base_time);

   for (Iterator iterator = range.first;
        iterator != range.second; iterator++)
   {
      const Integer fh = iterator->second;
      if (forecast_hour == fh)
      {
         const Integer prev_fh = (--iterator)->second;
         const Key prev_key (base_time, prev_fh);
         return prev_key;
      }
   }

   throw Nwp_Exception ("There is no previous key");

}

set<Dtime>
Nwp::Key_Multimap::get_base_time_set () const
{

   set<Dtime> base_time_set;
   typedef Key_Multimap::const_iterator Iterator;

   for (Iterator iterator = begin (); iterator != end (); iterator++)
   {
      const Dtime& base_time = iterator->first;
      base_time_set.insert (base_time);
   }

   return base_time_set;

}

set<Dtime>
Nwp::Key_Multimap::get_valid_time_set (const Dtime& base_time) const
{

   set<Dtime> valid_time_set;

   if (base_time.is_nat ())
   {
      for (Key_Multimap::const_reverse_iterator iterator = rbegin ();
           iterator != rend (); iterator++)
      {
         const Dtime& bt = iterator->first;
         const Integer fh = iterator->second;
         const Dtime dtime (bt.t + fh);
         valid_time_set.insert (dtime);
      }

   }
   else
   {

      typedef Key_Multimap::const_iterator Iterator;
      const pair<Iterator, Iterator> range = equal_range (base_time);

      for (Iterator iterator = range.first;
           iterator != range.second; iterator++)
      {
         const Dtime& bt = iterator->first;
         const Integer fh = iterator->second;
         const Dtime dtime (bt.t + fh);
         valid_time_set.insert (dtime);
      }

   }

   return valid_time_set;

}

Nwp::Key
Nwp::Key_Multimap::get_key (const Dtime& dtime,
                            const Dtime& base_time) const
{

   if (size () == 0) { throw Nwp_Exception ("Nwp is empty."); }

   if (!base_time.is_nat ())
   {
      const Dtime& bt = base_time;
      const Integer fh = Integer (round (dtime.t - bt.t));
      return Key (bt, fh);
   }

   typedef Key_Multimap::const_iterator Iterator;
   set<Nwp::Key> key_set;
   for (Iterator iterator = begin (); iterator != end (); iterator++)
   {
      const Dtime& bt = iterator->first;
      const Integer fh = iterator->second;
      if (Dtime (bt.t + fh) == dtime) { key_set.insert (Key (bt, fh)); }
   }

   // Return the key with latest base_time, except when it's fh = 0
   const Integer n = key_set.size ();
   if (n == 0) { throw Nwp_Exception ("No Match Key"); }
   else
   {
      const Nwp::Key& last_key = *(key_set.rbegin ());
      if (last_key.forecast_hour != 0 || n == 1) { return last_key; }
      else
      {
         const Nwp::Key& second_last_key = *(++key_set.rbegin ());
         return second_last_key;
      }
   }

}
 
void
Nwp::fill_lapse_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Level& level)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* next_p_data_ptr = get_next_p_up_data_ptr (key, level);
   Gvd_2d* this_t_data_ptr = get_temperature_data_ptr (key, level);
   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   const Data_3D& data_3d = get_3d_data (key);
   const Tuple& tuple_p = data_3d.get_tuple_p (TEMPERATURE);

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real next_p = next_p_data_ptr->get_datum (0, i, j);
         const Real next_t = data_3d.evaluate (TEMPERATURE, next_p, i, j);
         const Real this_t = this_t_data_ptr->get_datum (0, i, j);
         gvd_2d.set_datum (vector_index, i, j, next_t - this_t);
      }
   }

   delete next_p_data_ptr;
   delete this_t_data_ptr;

}

void
Nwp::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                     const Integer vector_index,
                     const Key& key,
                     const Integer hours)
{

   const Level& nil = Level::nil_level ();

   if (hours == key.forecast_hour)
   {
      fill_data (gvd_2d, vector_index, key, nil, RAINFALL_CUMULATIVE);
      return;
   }

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
   const Integer prev_forecast_hour = key.forecast_hour - hours;
   const Key prev_key (key.base_time, prev_forecast_hour);

   try
   {

      fill_data (*data_ptr, 0, key, nil, RAINFALL_CUMULATIVE);
      fill_data (*data_ptr, 1, prev_key, nil, RAINFALL_CUMULATIVE);

      const Size_2D& size_2d = gvd_2d.get_size_2d ();

      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real c_datum = data_ptr->get_datum (0, i, j);
            const Real p_datum = data_ptr->get_datum (1, i, j);
            const Real precipitation = c_datum - p_datum;
            gvd_2d.set_datum (vector_index, i, j, precipitation);
         }
      }

   }
   catch (const Nwp_Exception& ne)
   {
      delete data_ptr;
      throw Nwp_Exception ("Time Span Not Valid");
   }

   delete data_ptr;

}

void
Nwp::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                     const Integer vector_index,
                     const Key& key,
                     const Nwp_Element nwp_element)
{

   switch (nwp_element)
   {

      case RAINFALL_STEP:
      {

         if (key_multimap.is_no_match (key))
         {
            gvd_2d.initialize (vector_index, GSL_NAN);
            return;
         }

         if (key_multimap.is_first_step (key))
         {
            gvd_2d.initialize (vector_index, 0);
            return;
         }

         const Key& prev_key = key_multimap.get_previous_key (key);

         Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
         fill_rain_data (*data_ptr, 0, prev_key, RAINFALL_CUMULATIVE);
         fill_rain_data (*data_ptr, 1, key, RAINFALL_CUMULATIVE);

         const Size_2D& size_2d = gvd_2d.get_size_2d ();

         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real& prev_rf = data_ptr->get_datum (0, i, j);
               const Real& this_rf = data_ptr->get_datum (1, i, j);
               const Real& step_rf = this_rf - prev_rf;
               gvd_2d.set_datum (vector_index, i, j, step_rf);
            }
         }

         delete data_ptr;
         return;

      }

   }

   throw Nwp_Exception ("Nwp::fill_rain_data N/A");

}

void
Nwp::fill_p_potential_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                      const Integer vector_index,
                                      const Key& key,
                                      const Real p)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   Nwp::Data_3D& data_3d = get_3d_data (key);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real pv = data_3d.evaluate (POTENTIAL_VORTICITY, p, i, j);
         gvd_2d.set_datum (vector_index, i, j, pv);
      }
   }

}

void
Nwp::fill_theta_potential_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                          const Integer vector_index,
                                          const Key& key,
                                          const Real theta)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   vector<Nwp_Element> ev;
   Gvd_2d* data_ptr = get_theta_level_data_ptr (key, theta, ev, true);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real datum = data_ptr->get_datum (0, i, j);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete data_ptr;

}

void
Nwp::fill_potential_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                    const Integer vector_index,
                                    const Key& key,
                                    const Level& level)
{

   const Integer vi = vector_index;

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         const Real p = level.value;
         fill_p_potential_vorticity_data (gvd_2d, vi, key, p);
         return;
      }

      case THETA_LEVEL:
      {
         const Real theta = level.value;
         fill_theta_potential_vorticity_data (gvd_2d, vi, key, theta);
         return;
      }

   }

   throw Nwp_Exception ("Nwp::fill_potential_vorticity_data N/A");

}

void
Nwp::fill_p_temperature_advection_data (Geodetic_Vector_Data_2D& gvd_2d,
                                        const Integer vector_index,
                                        const Key& key,
                                        const Real p)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   Nwp::Data_3D& data_3d = get_3d_data (key);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real pv = data_3d.evaluate (TEMPERATURE_ADVECTION, p, i, j);
         gvd_2d.set_datum (vector_index, i, j, pv);
      }
   }

}

void
Nwp::fill_theta_temperature_advection_data (Geodetic_Vector_Data_2D& gvd_2d,
                                            const Integer vector_index,
                                            const Key& key,
                                            const Real theta)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   vector<Nwp_Element> ev;
   ev.push_back (TEMPERATURE_ADVECTION);
   Gvd_2d* data_ptr = get_theta_level_data_ptr (key, theta, ev, false);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real datum = data_ptr->get_datum (0, i, j);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete data_ptr;

}

void
Nwp::fill_temperature_advection_data (Geodetic_Vector_Data_2D& gvd_2d,
                                      const Integer vector_index,
                                      const Key& key,
                                      const Level& level)
{

   const Integer vi = vector_index;

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         const Real p = level.value;
         fill_p_temperature_advection_data (gvd_2d, vi, key, p);
         return;
      }

      case THETA_LEVEL:
      {
         const Real theta = level.value;
         fill_theta_temperature_advection_data (gvd_2d, vi, key, theta);
         return;
      }

   }

   throw Nwp_Exception ("Nwp::fill_potential_vorticity_data N/A");

}

void
Nwp::fill_p_adiabatic_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                    const Integer vector_index,
                                    const Key& key,
                                    const Real p)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   Nwp::Data_3D& data_3d = get_3d_data (key);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real pv = data_3d.evaluate (ADIABATIC_HEATING, p, i, j);
         gvd_2d.set_datum (vector_index, i, j, pv);
      }
   }

}

void
Nwp::fill_theta_adiabatic_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                        const Integer vector_index,
                                        const Key& key,
                                        const Real theta)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   vector<Nwp_Element> ev;
   ev.push_back (ADIABATIC_HEATING);
   Gvd_2d* data_ptr = get_theta_level_data_ptr (key, theta, ev, false);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real datum = data_ptr->get_datum (0, i, j);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete data_ptr;

}

void
Nwp::fill_adiabatic_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                  const Integer vector_index,
                                  const Key& key,
                                  const Level& level)
{

   const Integer vi = vector_index;

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         const Real p = level.value;
         fill_p_adiabatic_heating_data (gvd_2d, vi, key, p);
         return;
      }

      case THETA_LEVEL:
      {
         const Real theta = level.value;
         fill_theta_adiabatic_heating_data (gvd_2d, vi, key, theta);
         return;
      }

   }

   throw Nwp_Exception ("Nwp::fill_adiabatic_heating_data N/A");

}

void
Nwp::fill_p_latent_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Real p)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   Nwp::Data_3D& data_3d = get_3d_data (key);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real pv = data_3d.evaluate (LATENT_HEATING, p, i, j);
         gvd_2d.set_datum (vector_index, i, j, pv);
      }
   }

}

void
Nwp::fill_theta_latent_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                     const Integer vector_index,
                                     const Key& key,
                                     const Real theta)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   vector<Nwp_Element> ev;
   ev.push_back (LATENT_HEATING);
   Gvd_2d* data_ptr = get_theta_level_data_ptr (key, theta, ev, false);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real datum = data_ptr->get_datum (0, i, j);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete data_ptr;

}

void
Nwp::fill_latent_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Level& level)
{

   const Integer vi = vector_index;

   try
   {

      switch (level.type)
      {

         case PRESSURE_LEVEL:
         {
            const Real p = level.value;
            fill_p_latent_heating_data (gvd_2d, vi, key, p);
            return;
         }

         case THETA_LEVEL:
         {
            const Real theta = level.value;
            fill_theta_latent_heating_data (gvd_2d, vi, key, theta);
            return;
         }

      }

   }
   catch (const std::exception& se)
   {
   }

   throw Nwp_Exception ("Nwp::latent_heating_data N/A");

}

void
Nwp::fill_pv_p_data (Geodetic_Vector_Data_2D& gvd_2d,
                     const Integer vector_index,
                     const Key& key,
                     const Real pv_threshold)
{


   const Nwp_Element PV = POTENTIAL_VORTICITY;

   const Integer nk = tuple_p.size ();
   const Data_3D& data_3d = get_3d_data (key);

   const Level& surface = Level::surface_level ();
   fill_data (gvd_2d, vector_index, key, surface, PRESSURE);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Real lowest_p = 500e2;

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real latitude = gvd_2d.get_coordinate (0, i);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         Real p = tuple_p.front ();
         const Real longitude = gvd_2d.get_coordinate (1, j);

         for (Integer k = 0; k < tuple_p.size () - 1; k++)
         {

            const Real lower_p = tuple_p[k]; 
            const Real upper_p = tuple_p[k + 1];
            const Real lower_pv = data_3d.evaluate (PV, lower_p, i, j); 
            const Real upper_pv = data_3d.evaluate (PV, upper_p, i, j);

            if (lower_p > lowest_p)
            {
               if (fabs (lower_pv) > fabs (pv_threshold))
               {
                  p = lowest_p;
               }
               break;
            }

            const Real lower_d_pv = fabs (pv_threshold) - fabs (lower_pv);
            const Real upper_d_pv = fabs (pv_threshold) - fabs (upper_pv);

            if (lower_d_pv * upper_d_pv <= 0)
            {

               const Real d_p = upper_p - lower_p;
               const Real d_pv = fabs (upper_pv) - fabs (lower_pv);

               p = lower_p + (lower_d_pv * d_p / d_pv);
               break;

            }

         }

         gvd_2d.set_datum (vector_index, i, j, p);

      }
   }

}

void
Nwp::fill_absolute_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                   const Integer vector_index,
                                   const Key& key,
                                   const Level& level)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* wind_ptr = get_wind_data_ptr (key, level);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real latitude = gvd_2d.get_coordinate (0, i);
      const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real zeta = wind_ptr->evaluate_2d (0, 1, i, j, VORTICITY_OP);
         gvd_2d.set_datum (vector_index, i, j, zeta + f);
      }

   }

   delete wind_ptr;

}

void
Nwp::fill_shear_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Level& level)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* wind_ptr = get_wind_data_ptr (key, level);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real shear = wind_ptr->evaluate_2d (0, 1, i, j, SHEAR_OP);
         gvd_2d.set_datum (vector_index, i, j, shear);
      }
   }

   delete wind_ptr;

}

void
Nwp::fill_curvature_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                    const Integer vector_index,
                                    const Key& key,
                                    const Level& level)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* wind_ptr = get_wind_data_ptr (key, level);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real datum = wind_ptr->evaluate_2d (0, 1, i, j, CURVATURE_OP);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete wind_ptr;

}

void
Nwp::fill_wind_data (Geodetic_Vector_Data_2D& gvd_2d,
                     const Integer vector_index_u,
                     const Integer vector_index_v,
                     const Key& key,
                     const Level& level)
{
   fill_data (gvd_2d, vector_index_u, key, level, ZONAL_WIND);
   fill_data (gvd_2d, vector_index_v, key, level, MERIDIONAL_WIND);
}

void
Nwp::fill_mix_down_temperature_data (Geodetic_Vector_Data_2D& gvd_2d,
                                     const Integer vector_index,
                                     const Key& key,
                                     const Level& level)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = get_initialized_vd_2d (2);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface = Level::surface_level ();

   fill_data (*data_ptr, 0, key, level, THETA);
   fill_data (*data_ptr, 1, key, surface, PRESSURE);

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real theta = data_ptr->get_datum (0, i, j);
         const Real surface_p = data_ptr->get_datum (1, i, j);
         const Real t = Thermo_Point::theta_p (theta, surface_p).get_t ();
         gvd_2d.set_datum (vector_index, i, j, t);
      }
   }

   delete data_ptr;

}

void
Nwp::fill_snow_probability_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = get_initialized_vd_2d (2);

   Nwp::Data_3D& data_3d = get_3d_data (key);
   const Level& surface = Level::surface_level ();
   fill_nil_level_data (*data_ptr, 0, key, THICKNESS);
   fill_data (*data_ptr, 1, key, surface, PRESSURE);

   const Real b0 = -134.0066;
   const Real b1 = -0.0079;
   const Real b2 = 0.57193;
   const Real b3 = 0.0042365;
   const Real b4 = 0.02552;
   const Nwp_Element T = TEMPERATURE;
   const Nwp_Element Z = GEOPOTENTIAL_HEIGHT;

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real t = data_3d.evaluate (T, Real (850e2), i, j) - K;
         const Real z = data_3d.evaluate (Z, Real (850e2), i, j);
         const Real thick = data_ptr->get_datum (0, i, j);
         const Real surface_p = data_ptr->get_datum (1, i, j);
         const Real e = data_3d.evaluate (Z, surface_p, i, j);
         const Real p = 1 / (1 + exp (b0 + b1*e + b2*t + b3*z + b4*thick));
         gvd_2d.set_datum (vector_index, i, j, p);
      }
   }

   delete data_ptr;

}


void
Nwp::fill_snow_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key,
                           const Real probability)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = get_initialized_vd_2d (3);

   fill_pressure_level_data (*data_ptr, 0, key, 850e2, TEMPERATURE);
   fill_pressure_level_data (*data_ptr, 1, key, 850e2, GEOPOTENTIAL_HEIGHT);
   fill_nil_level_data (*data_ptr, 2, key, THICKNESS);

   const Real b0 = -134.0066;
   const Real b1 = -0.0079;
   const Real b2 = 0.57193;
   const Real b3 = 0.0042365;
   const Real b4 = 0.02552;
   const Real k = log ((1 - probability) / probability);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real t = data_ptr->get_datum (0, i, j) - K;
         const Real z = data_ptr->get_datum (1, i, j);
         const Real thick = data_ptr->get_datum (2, i, j);
         const Real snow_level = (k - (b0 + b2*t + b3*z + b4*thick)) / b1;
         gvd_2d.set_datum (vector_index, i, j, snow_level);
      }
   }

   delete data_ptr;

}

void
Nwp::fill_fdi_data (Geodetic_Vector_Data_2D& gvd_2d,
                    const Integer vector_index,
                    const Key& key,
                    const Level& level,
                    const Nwp_Element nwp_element)
{

   if (nwp_element != FFDI && nwp_element != GFDI)
   {
      throw Nwp_Exception ("Nwp::fill_fdi_data N/A");
   }

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   const Level& screen = Level::screen_level ();
   Gvd_2d* data_ptr = get_initialized_vd_2d (4);

   fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
   fill_data (*data_ptr, 1, key, screen, RELATIVE_HUMIDITY);
   fill_data (*data_ptr, 2, key, level, ZONAL_WIND);
   fill_data (*data_ptr, 3, key, level, MERIDIONAL_WIND);

   const Real df = 10;
   const Real curing = 100;
   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real t = data_ptr->get_datum (0, i, j) - K;
         const Real rh = data_ptr->get_datum (1, i, j) * 100;
         const Real u = data_ptr->get_datum (2, i, j);
         const Real v = data_ptr->get_datum (3, i, j);
         const Real kph = sqrt (u*u + v*v) * 3.6;

         Real fdi;

         switch (nwp_element)
         {

            case FFDI:
            {
               fdi = Fire::get_ffdi (t, rh, kph, df);
               break;
            }

            case GFDI:
            {
               fdi = Fire::get_gfdi (t, rh, kph, curing);
               break;
            }

         }

         gvd_2d.set_datum (vector_index, i, j, fdi);
      
      }
   }

   delete data_ptr;

}

void
Nwp::fill_total_totals_data (Geodetic_Vector_Data_2D& gvd_2d,
                             const Integer vector_index,
                             const Key& key)
{

   const Nwp_Element T = TEMPERATURE;
   const Nwp_Element TD = DEW_POINT;

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (3);
   fill_pressure_level_data (*temp_data_ptr, 0, key, 850e2, TD);
   fill_pressure_level_data (*temp_data_ptr, 1, key, 850e2, T);
   fill_pressure_level_data (*temp_data_ptr, 2, key, 500e2, T);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real td_850 = temp_data_ptr->get_datum (0, i, j);
         const Real t_850 = temp_data_ptr->get_datum (1, i, j);
         const Real t_500 = temp_data_ptr->get_datum (2, i, j);

         const Real datum = Instability::get_total_totals (
            t_850, t_500, td_850);
         gvd_2d.set_datum (vector_index, i, j, datum);

      }
   }

   delete temp_data_ptr;

}

void
Nwp::fill_continuous_haines_data (Geodetic_Vector_Data_2D& gvd_2d,
                                  const Integer vector_index,
                                  const Key& key)
{

   const Nwp_Element T = TEMPERATURE;
   const Nwp_Element TD = DEW_POINT;

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (3);
   fill_pressure_level_data (*temp_data_ptr, 0, key, 850e2, TD);
   fill_pressure_level_data (*temp_data_ptr, 1, key, 850e2, T);
   fill_pressure_level_data (*temp_data_ptr, 2, key, 700e2, T);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real td_850 = temp_data_ptr->get_datum (0, i, j);
         const Real t_850 = temp_data_ptr->get_datum (1, i, j);
         const Real t_700 = temp_data_ptr->get_datum (2, i, j);

         const Real ca = (t_850 - t_700) / 2 - 2;
         const Real cb = std::min ((t_850 - td_850), Real (30.0)) / 3 - 1;
         const Real datum = ca + (cb > 5 ? (cb - 5) / 2 + 5 : cb);

         gvd_2d.set_datum (vector_index, i, j, datum);

      }
   }

   delete temp_data_ptr;

}

void
Nwp::fill_k_index_data (Geodetic_Vector_Data_2D& gvd_2d,
                        const Integer vector_index,
                        const Key& key)
{

   const Nwp_Element T = TEMPERATURE;
   const Nwp_Element TD = DEW_POINT;

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (5);
   fill_pressure_level_data (*temp_data_ptr, 0, key, 850e2, TD);
   fill_pressure_level_data (*temp_data_ptr, 1, key, 700e2, TD);
   fill_pressure_level_data (*temp_data_ptr, 2, key, 850e2, T);
   fill_pressure_level_data (*temp_data_ptr, 3, key, 700e2, T);
   fill_pressure_level_data (*temp_data_ptr, 4, key, 500e2, T);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real td_850 = temp_data_ptr->get_datum (0, i, j);
         const Real td_700 = temp_data_ptr->get_datum (1, i, j);
         const Real t_850 = temp_data_ptr->get_datum (2, i, j);
         const Real t_700 = temp_data_ptr->get_datum (3, i, j);
         const Real t_500 = temp_data_ptr->get_datum (4, i, j);

         const Real datum = Instability::get_k_index (
            t_850, t_700, t_500, td_850, td_700);
         gvd_2d.set_datum (vector_index, i, j, datum);

      }
   }

   delete temp_data_ptr;

}

void
Nwp::fill_lifted_index_data (Geodetic_Vector_Data_2D& gvd_2d,
                             const Integer vector_index,
                             const Key& key,
                             const Real start_p,
                             const Real end_p)
{

   const Nwp_Element T = TEMPERATURE;
   const Nwp_Element TD = DEW_POINT;

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (3);
   fill_pressure_level_data (*temp_data_ptr, 0, key, start_p, T);
   fill_pressure_level_data (*temp_data_ptr, 1, key, start_p, TD);
   fill_pressure_level_data (*temp_data_ptr, 2, key, end_p, T);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real start_t = temp_data_ptr->evaluate (0, i, j);
         const Real start_td = temp_data_ptr->evaluate (1, i, j);
         const Real end_t = temp_data_ptr->evaluate (2, i, j);

         const Real datum = Instability::get_lifted_index (
            start_p, start_t, start_td, end_p, end_t);

         gvd_2d.set_datum (vector_index, i, j, datum);

      }
   }

   delete temp_data_ptr;

}

void
Nwp::fill_surface_lifted_index_data (Geodetic_Vector_Data_2D& gvd_2d,
                                     const Integer vector_index,
                                     const Key& key,
                                     const Real top_p)
{

   const Level& screen = Level::screen_level ();
   const Level& surface = Level::surface_level ();
   const Level& top_p_level = Level::pressure_level (top_p);

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (4);
   fill_data (*temp_data_ptr, 0, key, screen, TEMPERATURE);
   fill_data (*temp_data_ptr, 1, key, screen, DEW_POINT);
   fill_data (*temp_data_ptr, 2, key, surface, PRESSURE);
   fill_data (*temp_data_ptr, 3, key, top_p_level, TEMPERATURE);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real screen_t = temp_data_ptr->get_datum (0, i, j);
         const Real screen_t_d = temp_data_ptr->get_datum (1, i, j);
         const Real surface_p = temp_data_ptr->get_datum (2, i, j);
         const Real top_t = temp_data_ptr->get_datum (3, i, j);

         const Real datum = Instability::get_lifted_index (
            surface_p, screen_t, screen_t_d, top_p, top_t);
         gvd_2d.set_datum (vector_index, i, j, datum);

      }
   }

   delete temp_data_ptr;

}

void
Nwp::fill_li_thunder_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key,
                           const Level& level,
                           const Real thunder_t)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   const Level& surface = Level::surface_level ();

   Gvd_2d* t_data_ptr = get_temperature_data_ptr (key, level);
   Gvd_2d* td_data_ptr = get_dew_point_data_ptr (key, level);
   Gvd_2d* p_data_ptr = get_temperature_p_data_ptr (key, thunder_t);

   Gvd_2d* surface_p_data_ptr = get_initialized_vd_2d (1);
   fill_data (*surface_p_data_ptr, 0, key, surface, PRESSURE);

   Real start_p;
   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {

         switch (level.type)
         {

            case PRESSURE_LEVEL:
               start_p = level.value;
               break;

            case MEAN_SEA_LEVEL:
            case TEN_METRE_LEVEL:
            case FIFTY_METRE_LEVEL:
            case SCREEN_LEVEL:
               start_p = surface_p_data_ptr->evaluate (0, i, j);
               break;
 
            default:
               string error_str = "Nwp::fill_thunder_li_data ";
               error_str += level.get_string () + " invalid";
               throw Nwp_Exception (error_str);
               break;

         }

         const Real end_p = p_data_ptr->evaluate (0, i, j);
         const bool invalid_p = (start_p < end_p || gsl_isnan (end_p));

         const Real start_t = t_data_ptr->evaluate (0, i, j);
         const Real start_td = td_data_ptr->evaluate (0, i, j);

         const Real datum = (invalid_p ? GSL_NAN :
            Instability::get_lifted_index (start_p, start_t,
               start_td, end_p, thunder_t - K) + K);

         gvd_2d.set_datum (vector_index, i, j, datum);

      }
   }

   delete t_data_ptr;
   delete td_data_ptr;
   delete p_data_ptr;
   delete surface_p_data_ptr;

}

void
Nwp::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                             const Integer vector_index,
                             const Key& key,
                             const Level& level,
                             const Nwp_Element nwp_element)
{

   const Integer vi = vector_index;

   switch (nwp_element)
   {

      case SLI:
      {
         fill_surface_lifted_index_data (gvd_2d, vi, key, 500e2);
         return;
      }

      case SHOWALTER:
      {
         fill_lifted_index_data (gvd_2d, vi, key, 850e2, 500e2);
         return;
      }

      case LI_700:
      {
         fill_lifted_index_data (gvd_2d, vi, key, 700e2, 400e2);
         return;
      }

      case K_INDEX:
      {
         fill_k_index_data (gvd_2d, vi, key);
         return;
      }

      case TOTAL_TOTALS:
      {
         fill_total_totals_data (gvd_2d, vi, key);
         return;
      }

      case LI_THUNDER:
      {
         // thunder_t is always -20°C
         fill_li_thunder_data (gvd_2d, vi, key, level/*,-20 + K*/);
         return;
      }

   }

   throw Nwp_Exception ("Nwp::fill_ts_diagnosis_data N/A");

}

void
Nwp::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Nwp_Element nwp_element)
{
   throw Nwp_Exception ("Nwp::fill_cloud_data N/A");
}

void
Nwp::fill_pressure_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Real p,
                               const Nwp_Element nwp_element)
{

   switch (nwp_element)
   {

      case denise::MIX_DOWN_TEMPERATURE:
      {
         const Level& level = Level::pressure_level (p);
         fill_mix_down_temperature_data (gvd_2d, vector_index, key, level);
         return;
      }

      case denise::POTENTIAL_VORTICITY:
      {
         fill_p_potential_vorticity_data (gvd_2d, vector_index, key, p);
         return;
      }

      case denise::ADIABATIC_HEATING:
      {
         fill_p_adiabatic_heating_data (gvd_2d, vector_index, key, p);
         return;
      }

   }

   Nwp::Data_3D& data_3d = get_3d_data (key);
   const Nwp_Element element = nwp_element;

   if ((p - tuple_p.front ()) * (p - tuple_p.back ()) > 0)
   {
      string error_str = "Nwp::fill_pressure_level_data_ptr ";
      error_str += string_render (" %f", p);
      throw Nwp_Exception (error_str);
   }

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_coordinate (0, i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_coordinate (1, j);
         const Real datum = data_3d.evaluate (element, p, latitude, longitude);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

}

void
Nwp::fill_theta_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                            const Integer vector_index,
                            const Key& key,
                            const Real theta,
                            const Nwp_Element nwp_element)
{

   const Integer vi = vector_index;

   switch (nwp_element)
   {

      case denise::POTENTIAL_VORTICITY:
      {
         fill_theta_potential_vorticity_data (gvd_2d, vi, key, theta);
         return;
      }

      case denise::ADIABATIC_HEATING:
      {
         fill_theta_adiabatic_heating_data (gvd_2d, vi, key, theta);
         return;
      }

   }

   throw Nwp_Exception ("Nwp::fill_theta_level_data N/A");

}

void
Nwp::fill_sigma_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                            const Integer vector_index,
                            const Key& key,
                            const Real sigma,
                            const Nwp_Element nwp_element,
                            const Geodetic_Vector_Data_2D& surface_p_data)
{

   const Nwp_Element& ne = nwp_element;
   const Real& start_p = tuple_p.front ();
   const Real& end_p = tuple_p.back ();

   const Data_3D& data_3d = get_3d_data (key);
   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   { 

      const Real latitude = gvd_2d.get_coordinate (0, i);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real longitude = gvd_2d.get_coordinate (1, j);
         const Real surface_p = surface_p_data.get_datum (0, i, j);
         const Real p = sigma * surface_p;
         const bool nan = ((p - start_p) * (p - end_p)) > 0;

         if (nan)
         {
            gvd_2d.set_datum (vector_index, i, j, GSL_NAN);
         }
         else
         {
            const Real datum = data_3d.evaluate (ne, p, latitude, longitude);
            gvd_2d.set_datum (vector_index, i, j, datum);
         }

      }

   }

}

void
Nwp::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                             const Integer vector_index,
                             const Key& key,
                             const Nwp_Element nwp_element)
{

   switch (nwp_element)
   {

      case DEW_POINT_DEPRESSION:
      {

         const Level& screen = Level::screen_level ();

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, DEW_POINT);
         const Size_2D& size_2d = gvd_2d.get_size_2d ();

         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real t_d = data_ptr->get_datum (1, i, j);
               const Real dew_point_depression = t - t_d;
               gvd_2d.set_datum (vector_index, i, j, dew_point_depression);
            }
         }

         delete data_ptr;
         return;

      }

      case THETA_E:
      {

         const Level& screen = Level::screen_level ();
         const Level& surface = Level::surface_level ();

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (3);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, DEW_POINT);
         fill_data (*data_ptr, 2, key, surface, PRESSURE);
         const Size_2D& size_2d = gvd_2d.get_size_2d ();

         for (Integer i = 0; i < size_2d.i; i++)
         {
            Thermo_Point thermo_point;
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real t_d = data_ptr->get_datum (1, i, j);
               const Real p = data_ptr->get_datum (2, i, j);

               thermo_point.set_normand (t - K, t_d - K, p);
               const Real theta_e = thermo_point.get_theta_e () + K;
               gvd_2d.set_datum (vector_index, i, j, theta_e);

            }
         }

         delete data_ptr;
         return;

      }

      case THETA_W:
      {

         const Level& screen = Level::screen_level ();
         const Level& surface = Level::surface_level ();

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (3);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, DEW_POINT);
         fill_data (*data_ptr, 2, key, surface, PRESSURE);
         const Size_2D& size_2d = gvd_2d.get_size_2d ();

         for (Integer i = 0; i < size_2d.i; i++)
         {
            Thermo_Point thermo_point;
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real t_d = data_ptr->get_datum (1, i, j);
               const Real p = data_ptr->get_datum (2, i, j);

               thermo_point.set_normand (t - K, t_d - K, p);
               const Real theta_w = thermo_point.get_theta_w () + K;
               gvd_2d.set_datum (vector_index, i, j, theta_w);

            }
         }

         delete data_ptr;
         return;

      }

   }

   throw Nwp_Exception ("Nwp::fill_screen_data N/A");

}

void
Nwp::fill_50m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Nwp_Element nwp_element)
{  
   throw Nwp_Exception ("Nwp::fill_50m_level_data N/A");
}

void
Nwp::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Nwp_Element nwp_element)
{  
   throw Nwp_Exception ("Nwp::fill_10m_level_data N/A");
}

void
Nwp::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                    const Integer vector_index,
                    const Key& key,
                    const Nwp_Element nwp_element)
{
   throw Nwp_Exception ("Nwp::fill_msl_data N/A");
}

void
Nwp::fill_nil_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Nwp_Element nwp_element)
{

   const Integer vi = vector_index;

   switch (nwp_element)
   {

      case CONTINUOUS_HAINES:
      {
         fill_continuous_haines_data (gvd_2d, vector_index, key);
         return;
      };

      case SLI:
      case SHOWALTER:
      case LI_700:
      case K_INDEX:
      case TOTAL_TOTALS:
      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& level = Level ();
         fill_ts_diagnosis_data (gvd_2d, vi, key, level, nwp_element);
         return;
      }

      case THICKNESS:
      {

         const Size_2D& size_2d = gvd_2d.get_size_2d ();
         const Level& p500 = Level::pressure_level (500e2);
         const Level& p1000 = Level::pressure_level (1000e2);

         Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (2);
         Geodetic_Vector_Data_2D& temp_data = *temp_data_ptr;
         fill_data (temp_data, 0, key, p500, GEOPOTENTIAL_HEIGHT);
         fill_data (temp_data, 1, key, p1000, GEOPOTENTIAL_HEIGHT);

         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real z500 = temp_data.get_datum (0, i, j);
               const Real z1000 = temp_data.get_datum (1, i, j);
               gvd_2d.set_datum (vector_index, i, j, z500 - z1000);
            }
         }

         delete temp_data_ptr;
         return;

      }

      case HIGH_CLOUD:
      case MIDDLE_CLOUD:
      case LOW_CLOUD:
      case TOTAL_CLOUD:
      {
         fill_cloud_data (gvd_2d, vi, key, nwp_element);
         return;
      }

      case PV1_5_PRESSURE:
      {
         fill_pv_p_data (gvd_2d, vi, key);
         return;
      }

   }

   throw Nwp_Exception ("Nwp::fill_nil_data N/A");

}

void
Nwp::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Nwp_Element nwp_element)
{
   throw Nwp_Exception ("Nwp::fill_surface_data N/A");
}

void
Nwp::fill_data (Geodetic_Vector_Data_2D& gvd_2d,
                const Integer vector_index,
                const Key& key,
                const Level& level,
                const Nwp_Element nwp_element)
{

   switch (nwp_element)
   {

      case SNOW_LEVEL:
      {
         fill_snow_level_data (gvd_2d, vector_index, key, 0.5);
         return;
      };

      case LAPSE:
      {
         fill_lapse_data (gvd_2d, vector_index, key, level);
         return;
      };

      case RAINFALL_STEP:
      case RAINFALL_CUMULATIVE:
      {
         fill_rain_data (gvd_2d, vector_index, key, nwp_element);
         return;
      };

      case GFDI:
      case FFDI:
      {
         fill_fdi_data (gvd_2d, vector_index, key, level, nwp_element);
         return;
      };

      case ADIABATIC_HEATING:
      {
         fill_adiabatic_heating_data (gvd_2d, vector_index, key, level);
         return;
      };

      case LATENT_HEATING:
      {
         fill_latent_heating_data (gvd_2d, vector_index, key, level);
         return;
      };

      case POTENTIAL_VORTICITY:
      {
         fill_potential_vorticity_data (gvd_2d, vector_index, key, level);
         return;
      };

      case ABSOLUTE_VORTICITY:
      {
         fill_absolute_vorticity_data (gvd_2d, vector_index, key, level);
         return;
      };

      case SHEAR_VORTICITY:
      {
         fill_shear_vorticity_data (gvd_2d, vector_index, key, level);
         return;
      };

      case CURVATURE_VORTICITY:
      {
         fill_curvature_vorticity_data (gvd_2d, vector_index, key, level);
         return;
      };

      case LI_THUNDER:
      {
         fill_li_thunder_data (gvd_2d, vector_index, key, level);
         return;
      };

      case WIND_SPEED:
      {

         Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
         fill_wind_data (*data_ptr, 0, 1, key, level);
         const Size_2D& size_2d = data_ptr->get_size_2d ();

         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real u = data_ptr->get_datum (0, i, j);
               const Real v = data_ptr->get_datum (1, i, j);
               const Real speed = sqrt (u*u + v*v);
               gvd_2d.set_datum (vector_index, i, j, speed);
            }
         }

         delete data_ptr;
         return;

      };

      case WIND_DIRECTION:
      {

         Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
         fill_wind_data (*data_ptr, 0, 1, key, level);
         const Size_2D& size_2d = data_ptr->get_size_2d ();

         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real u = data_ptr->get_datum (0, i, j);
               const Real v = data_ptr->get_datum (1, i, j);
               const Real direction = Wind (u, v).get_direction ();
               gvd_2d.set_datum (vector_index, i, j, direction);
            }
         }

         delete data_ptr;
         return;

      };

   }

   switch (level.type)
   {
   
      case PRESSURE_LEVEL:
      {
         const Real v = level.value;
         const Integer vi = vector_index;
         fill_pressure_level_data (gvd_2d, vi, key, v, nwp_element);
         return;
      }

      case THETA_LEVEL:
      {
         const Real v = level.value;
         const Integer vi = vector_index;
         fill_theta_level_data (gvd_2d, vi, key, v, nwp_element);
         return;
      }

      case SIGMA_LEVEL:
      {
         const Real v = level.value;
         const Integer vi = vector_index;
         const Level& surface = Level::surface_level ();
         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (1);
         fill_data (*data_ptr, 0, key, surface, PRESSURE);
         fill_sigma_level_data (gvd_2d, vi, key, v, nwp_element, *data_ptr);
         delete data_ptr;
         return;
      }

      case SCREEN_LEVEL:
      {
         fill_screen_level_data (gvd_2d, vector_index, key, nwp_element);
         return;
      }

      case FIFTY_METRE_LEVEL:
      {
         fill_50m_level_data (gvd_2d, vector_index, key, nwp_element);
         return;
      }

      case TEN_METRE_LEVEL:
      {
         fill_10m_level_data (gvd_2d, vector_index, key, nwp_element);
         return;
      }

      case MEAN_SEA_LEVEL:
      {
         fill_msl_data (gvd_2d, vector_index, key, nwp_element);
         return;
      }

      case NIL_LEVEL:
      {
         fill_nil_level_data (gvd_2d, vector_index, key, nwp_element);
         return;
      }

      case SURFACE_LEVEL:
      {
         fill_surface_level_data (gvd_2d, vector_index, key, nwp_element);
         return;
      }

   }

   throw Nwp_Exception ("Nwp::fill_data N/A");

}

Nwp::Nwp (const string& description,
          const string& path)
   : description (description),
     path (path)
{
}

Nwp::~Nwp ()
{
   clear_data_3d_ptr_map ();
}

void
Nwp::clear_data_3d_ptr_map ()
{

   typedef map<Key, Data_3D*>::iterator Iterator;

   for (Iterator iterator = data_3d_ptr_map.begin ();
        iterator != data_3d_ptr_map.end (); iterator++)
   {
      Data_3D* data_3d_ptr = iterator->second;
      delete data_3d_ptr;
   }

   data_3d_ptr_map.clear ();

}

void
Nwp::set_domain_2d (const Domain_2D& domain_2d)
{
}

const string&
Nwp::get_description () const
{
   return description;
}

const string&
Nwp::get_status () const
{
   return status;
}

set<Dtime>
Nwp::get_valid_time_set (const Dtime& base_time) const
{
   return key_multimap.get_valid_time_set (base_time);
}

Tuple
Nwp::get_valid_t_tuple (const Dtime& base_time,
                        const Dtime& start_time,
                        const Dtime& end_time) const
{

   Tuple t_tuple;
   const set<Dtime>& full_vts = get_valid_time_set (base_time);

   for (set<Dtime>::const_iterator this_iterator = full_vts.begin ();
        this_iterator != full_vts.end (); this_iterator++)
   {

      set<Dtime>::const_iterator next_iterator = this_iterator;
      const Dtime& this_dtime = *(this_iterator);
      const Dtime& next_dtime = *(next_iterator);

      bool this_t_tick, next_t_tick;
      this_t_tick = (this_dtime.t > (start_time.t - TIME_TOLERANCE) &&
                     this_dtime.t < (end_time.t + TIME_TOLERANCE));
      next_t_tick = (next_dtime.t > (start_time.t - TIME_TOLERANCE) &&
                     next_dtime.t < (end_time.t + TIME_TOLERANCE));

      if (this_t_tick || next_t_tick)
      {
         t_tuple.push_back (this_dtime.t);
      }

   }

   return t_tuple;

}

set<Dtime>
Nwp::get_base_time_set () const
{
   return key_multimap.get_base_time_set ();
}

Nwp::Key
Nwp::get_key (const Dtime& dtime,
              const Dtime& base_time) const
{
   return key_multimap.get_key (dtime, base_time);
}

Nwp::Data_3D&
Nwp::get_3d_data (const Key& key)
{

   typedef Data_3D Nd_3d;
   typedef map<Key, Nd_3d*>::const_iterator Iterator;

   Iterator iterator = data_3d_ptr_map.find (key);

   if (iterator == data_3d_ptr_map.end ())
   {
      string error_str = "Nwp::get_3d_data get_3d_data failed  ";
      const Dtime& base_time = key.base_time;
      const Integer forecast_hour = key.forecast_hour;
      error_str += base_time.get_string ();
      error_str += string_render (" +%d hr", forecast_hour);
      throw Nwp_Exception (error_str);
   }

   Data_3D& data_3d = *(iterator->second);
   if (!data_3d.is_available ()) { load_3d_data (data_3d); }
   return data_3d;

}

Geodetic_Vector_Data_2D*
Nwp::get_ts_steering_data_ptr (const Key& key)
{
   const Layer layer (PRESSURE_LEVEL, 600e2, 800e2);
   return get_steering_data_ptr (key, layer);
}

Geodetic_Vector_Data_2D*
Nwp::get_steering_data_ptr (const Key& key,
                            const Layer& layer)
{

   if (layer.type != PRESSURE_LEVEL)
   {
      throw Nwp_Exception ("Can't do non-P Layer");
   }

   const Real p_a = layer.value;
   const Real p_b = layer.value_;

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (2);

   const Integer nk = tuple_p.size ();
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   Integer n = 0;
   data_ptr->initialize (0, 0);
   data_ptr->initialize (1, 0);

   for (Integer k = 0; k < nk; k++)
   {

      const Real p = tuple_p[k];
      if ((p - p_a) * (p - p_b) > 0) { continue; }
      const Level& level = Level::pressure_level (p);

      try
      {

         fill_wind_data (*temp_data_ptr, 0, 1, key, level);

         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {

               const Real u = temp_data_ptr->get_datum (0, i, j);
               const Real v = temp_data_ptr->get_datum (1, i, j);

               data_ptr->get_datum (0, i, j) += u;
               data_ptr->get_datum (1, i, j) += v;

            }
         }

         n++;

      }
      catch (const Exception& e)
      {
      }

   }

   delete temp_data_ptr;

   const Real reciprocal = Real (1) / n;
   data_ptr->scale_offset (0, reciprocal, 0);
   data_ptr->scale_offset (1, reciprocal, 0);

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_vertical_shear_data_ptr (const Key& key,
                                  const Layer& layer)
{

   const Level& level_a = layer.get_level_0 ();
   const Level& level_b = layer.get_level_1 ();

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (2);

   fill_wind_data (*data_ptr, 0, 1, key, level_a);
   fill_wind_data (*temp_data_ptr, 0, 1, key, level_b);

   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real u_b = temp_data_ptr->get_datum (0, i, j);
         const Real v_b = temp_data_ptr->get_datum (1, i, j);

         data_ptr->get_datum (0, i, j) -= u_b;
         data_ptr->get_datum (1, i, j) -= v_b;

      }
   }

   delete temp_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_freezing_level_data_ptr (const Key& key)

{

   const Nwp_Element T = TEMPERATURE;
   const Nwp_Element Z = GEOPOTENTIAL_HEIGHT;

   const Data_3D& data_3d = get_3d_data (key);
   const Geodetic_Vector_Data_3D& t_data_3d = data_3d.get_gvd_3d (T);
   const Geodetic_Vector_Data_3D& z_data_3d = data_3d.get_gvd_3d (Z);
   const Tuple& tuple_p = data_3d.get_tuple_p (T);
   const Integer nk = tuple_p.size ();

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {

      Real* array = new Real[3];

      for (Integer j = 0; j < size_2d.j; j++)
      {

         Integer filled = 0;
         array[0] = GSL_NAN;
         array[1] = GSL_NAN;
         array[2] = GSL_NAN;

         for (Integer k = 0; k < nk - 1; k++)
         {

            const Real lower_p = tuple_p[k];
            const Real upper_p = tuple_p[k + 1];
            const Real lower_t = t_data_3d.evaluate (0, lower_p, i, j) - K; 
            const Real upper_t = t_data_3d.evaluate (0, upper_p, i, j) - K;

            if (lower_t * upper_t <= 0)
            {

               const Real lower_z = z_data_3d.evaluate (0, lower_p, i, j);
               const Real upper_z = z_data_3d.evaluate (0, upper_p, i, j);
               const Real dz = upper_z - lower_z;
               const Real dt = upper_t - lower_t;

               const Real z = lower_z - (lower_t * dz / dt);

               array[filled] = z;
               if (filled < 2) { filled++; }

            }

         }

         data_ptr->set_datum (0, i, j, array[0]);
         data_ptr->set_datum (1, i, j, array[1]);
         data_ptr->set_datum (2, i, j, array[2]);

      }

      delete[] array;

   }

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_temperature_p_data_ptr (const Key& key,
                                 const Real temperature)

{

   const Data_3D& data_3d = get_3d_data (key);

   const Real t = temperature;
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
   const Size_2D& size_2d = data_ptr->get_size_2d ();
   const Geodetic_Vector_Data_3D& t_data_3d = data_3d.get_gvd_3d (TEMPERATURE);

   const Tuple& tuple_p = data_3d.get_tuple_p (TEMPERATURE);
   const Integer nk = tuple_p.size ();

   for (Integer i = 0; i < size_2d.i; i++)
   {

      for (Integer j = 0; j < size_2d.j; j++)
      {

         Real p = GSL_NAN;

         for (Integer k = 0; k < nk - 1; k++)
         {

            const Real lower_p = tuple_p[k];
            const Real upper_p = tuple_p[k + 1];

            const Real lower_t = t_data_3d.evaluate (0, lower_p, i, j); 
            const bool invalid_lower_t =
               gsl_isnan (lower_t) || (lower_t > 350) || (lower_t < -150);
            if (invalid_lower_t) { continue; }

            const Real upper_t = t_data_3d.evaluate (0, upper_p, i, j);
            const bool invalid_upper_t =
               gsl_isnan (upper_t) || (upper_t > 350) || (upper_t < -150);
            if (invalid_upper_t) { continue; }

            const Real delta_t = t - lower_t;
            const bool match = (delta_t * (t - upper_t) <= 0);

            if (match)
            {
               const Real dp = upper_p - lower_p;
               const Real dt = upper_t - lower_t;
               p = lower_p + (delta_t  * dp / dt);
               break;
            }

         }

         data_ptr->set_datum (0, i, j, p);

      }
   }

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_p_data_ptr (const Key& key,
                     const Level& level)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
         const Size_2D& size_2d = data_ptr->get_size_2d ();
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               data_ptr->set_datum (0, i, j, level.value);
            }
         }
         return data_ptr;
      }

      case THETA_LEVEL:
      {

         Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
         const Size_2D& size_2d = data_ptr->get_size_2d ();
         const Real theta = level.value;
         const Data_3D& data_3d = get_3d_data (key);
         const Geodetic_Vector_Data_3D& t_data_3d =
            data_3d.get_gvd_3d (TEMPERATURE);
         const Tuple& tuple_p = data_3d.get_tuple_p (TEMPERATURE);

         for (Integer i = 0; i < size_2d.i; i++)
         {

            Thermo_Point thermo_point;

            for (Integer j = 0; j < size_2d.j; j++)
            {

               Real p = GSL_NAN;

               for (Integer k = 0; k < tuple_p.size () - 1; k++)
               {

                  const Real lower_p = tuple_p[k]; 
                  const Real upper_p = tuple_p[k + 1];
                  const Real lower_t = t_data_3d.get_datum (0, k, i, j); 
                  const Real upper_t = t_data_3d.get_datum (0, k + 1, i, j);

                  thermo_point.set_t_p (lower_t - K, lower_p);
                  const Real lower_theta = thermo_point.get_theta () + K;
                  thermo_point.set_t_p (upper_t - K, upper_p);
                  const Real upper_theta = thermo_point.get_theta () + K;

                  const Real lower_d_theta = theta - lower_theta;
                  const Real upper_d_theta = theta - upper_theta;

                  if (lower_d_theta * upper_d_theta <= 0)
                  {

                     const Real d_p = upper_p - lower_p;
                     const Real d_theta = upper_theta - lower_theta;

                     p = lower_p + (lower_d_theta * d_p / d_theta);
                     break;

                  }

               }

               data_ptr->set_datum (0, i, j, p);

            }
         }

      }

      //case SIGMA_LEVEL:
      //{
      //   break;
      //}

      case SCREEN_LEVEL:
      case FIFTY_METRE_LEVEL:
      case TEN_METRE_LEVEL:
      case MEAN_SEA_LEVEL:
      case SURFACE_LEVEL:
      {
         return get_data_ptr (key, Level::surface_level (), PRESSURE);
      }
   }

   throw Exception ("Not Available for this level type");

}

Geodetic_Vector_Data_2D*
Nwp::get_next_p_up_data_ptr (const Key& key,
                             const Level& level)
{

   Geodetic_Vector_Data_2D* this_p_data_ptr = get_p_data_ptr (key, level);
   Geodetic_Vector_Data_2D* next_p_data_ptr = get_initialized_vd_2d (1);
   const Size_2D& size_2d = next_p_data_ptr->get_size_2d ();

   const Data_3D& data_3d = get_3d_data (key);
   const Geodetic_Vector_Data_3D& t_data_3d = data_3d.get_gvd_3d (TEMPERATURE);
   const Tuple& tuple_p = data_3d.get_tuple_p (TEMPERATURE);

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         Real next_p = GSL_NAN;
         const Real this_p = this_p_data_ptr->get_datum (0, i, j);
         for (Integer k = tuple_p.size () - 1; k >= 0; k--)
         {
            const Real p = tuple_p[k];
            if (p < this_p) { next_p = p; break;}
         }
         next_p_data_ptr->set_datum (0, i, j, next_p);
      }
   }

   delete this_p_data_ptr;
   return next_p_data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_cloud_data_ptr (const Key& key)
{

   const Level& nil_level = Level::nil_level ();

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
   fill_data (*data_ptr, 0, key, nil_level, denise::HIGH_CLOUD);
   fill_data (*data_ptr, 1, key, nil_level, denise::MIDDLE_CLOUD);
   fill_data (*data_ptr, 2, key, nil_level, denise::LOW_CLOUD);

   return data_ptr;


}

Geodetic_Vector_Data_2D*
Nwp::get_cloud_base_data_ptr (const Key& key,
                              const Level& level)
{

   if (level.type != PRESSURE_LEVEL) { throw Exception ("Only for P Levels"); }

   const Nwp_Element T = TEMPERATURE;
   const Nwp_Element Z = GEOPOTENTIAL_HEIGHT;

   const Data_3D& data_3d = get_3d_data (key);

   const Geodetic_Vector_Data_3D& z_data_3d = data_3d.get_gvd_3d (Z);
   const Tuple& tuple_p = data_3d.get_tuple_p (Z);
   const Integer nk = tuple_p.size ();

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   const Level& surface = Level::surface_level ();
   fill_data (*data_ptr, 0, key, surface, PRESSURE);

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real latitude = data_ptr->get_coordinate (0, i);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real longitude = data_ptr->get_coordinate (1, j);

         const Real surface_p = data_ptr->get_datum (0, i, j);
         data_ptr->set_datum (0, i, j, GSL_NAN);

         for (Integer k = nk - 1; k >= 0; k--)
         {

            const Real p = tuple_p[k];
            if (p > surface_p || p > level.value) { continue; }

            const Real rh = data_3d.evaluate (
               RELATIVE_HUMIDITY, k, latitude, longitude);

            if ((rh >= 0.85) ||
                (p < 850e2 && rh >= 0.8) ||
                (p < 600e2 && rh >= 0.75))
            {
               const Real datum = z_data_3d.evaluate (0, p, i, j);
               data_ptr->set_datum (0, i, j, datum);
               break;
            }

         }

      }
   }

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_cloud_top_temp_data_ptr (const Key& key)
{

   const Data_3D& data_3d = get_3d_data (key);

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   const Geodetic_Vector_Data_3D& t_data_3d = data_3d.get_gvd_3d (TEMPERATURE);
   const Tuple& tuple_p = data_3d.get_tuple_p (TEMPERATURE);
   const Integer nk = tuple_p.size ();

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real latitude = data_ptr->get_coordinate (0, i);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real longitude = data_ptr->get_coordinate (1, j);

         data_ptr->set_datum (0, i, j, 30 + K);

         for (Integer k = 0; k < nk; k++)
         {

            const Real p = tuple_p[k];
            const Real rh = data_3d.evaluate (
               RELATIVE_HUMIDITY, p, latitude, longitude);

            if ((rh >= 0.85) ||
                (p < 850e2 && rh >= 0.8) ||
                (p < 600e2 && rh >= 0.75))
            {
               const Real t = t_data_3d.evaluate (0, p, i, j);
               data_ptr->set_datum (0, i, j, t);
               break;
            }

         }

      }
   }

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_ageostrophic_wind_data_ptr (const Key& key,
                                     const Real p)
{

   const Level& level = Level::pressure_level (p);

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (3);
   fill_data (*temp_data_ptr, 0, key, level, GEOPOTENTIAL_HEIGHT);
   fill_data (*temp_data_ptr, 1, key, level, ZONAL_WIND);
   fill_data (*temp_data_ptr, 2, key, level, MERIDIONAL_WIND);
   const Geodetic_Vector_Data_2D& temp_data = *temp_data_ptr;

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real latitude = data_ptr->get_coordinate (0, i);
      const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real longitude = data_ptr->get_coordinate (1, j);

         const Real z_x = temp_data.evaluate (0, i, j, DX);
         const Real z_y = temp_data.evaluate (0, i, j, DY);
         const Real u = temp_data.get_datum (1, i, j);
         const Real v = temp_data.get_datum (2, i, j);

         const Real u_g = -g * z_y / f;
         const Real v_g = g * z_x / f;
         const Real u_a = u - u_g;
         const Real v_a = v - v_g;

         data_ptr->set_datum (0, i, j, u_a);
         data_ptr->set_datum (1, i, j, v_a);

      }
   }

   delete temp_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_surface_ageostrophic_wind_data_ptr (const Key& key,
                                             const Level& level)
{

   const Level& msl = Level::mean_sea_level ();
   const Level& screen = Level::screen_level ();

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (4);
   fill_data (*temp_data_ptr, 0, key, msl, PRESSURE);
   fill_data (*temp_data_ptr, 1, key, screen, TEMPERATURE);
   fill_data (*temp_data_ptr, 2, key, level, ZONAL_WIND);
   fill_data (*temp_data_ptr, 3, key, level, MERIDIONAL_WIND);
   const Geodetic_Vector_Data_2D& temp_data = *temp_data_ptr;

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real latitude = data_ptr->get_coordinate (0, i);
      const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real longitude = data_ptr->get_coordinate (1, j);

         const Real p = temp_data.get_datum (0, i, j);
         const Real t = temp_data.get_datum (1, i, j);
         const Real u = temp_data.get_datum (2, i, j);
         const Real v = temp_data.get_datum (3, i, j);

         const Real p_x = temp_data.evaluate (0, latitude, longitude, DX);
         const Real p_y = temp_data.evaluate (0, latitude, longitude, DY);
         const Real rho = p / (R_d * t);

         const Real u_g = -p_y / (f * rho);
         const Real v_g = p_x / (f * rho);
         const Real u_a = u - u_g;
         const Real v_a = v - v_g;

         data_ptr->set_datum (0, i, j, u_a);
         data_ptr->set_datum (1, i, j, v_a);

      }
   }

   delete temp_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_ageostrophic_wind_data_ptr (const Key& key,
                                     const Level& level)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         const Real p = level.value;
         return get_ageostrophic_wind_data_ptr (key, p);
      }

      case TEN_METRE_LEVEL:
      case FIFTY_METRE_LEVEL:
      {
         return get_surface_ageostrophic_wind_data_ptr (key, level);
      }

   }

   throw Nwp_Exception ("Ageostrophic Wind N/A for this level");

}

Geodetic_Vector_Data_2D*
Nwp::get_geostrophic_wind_data_ptr (const Key& key,
                                    const Real p)
{

   const Level& level = Level::pressure_level (p);

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (1);
   fill_data (*temp_data_ptr, 0, key, level, GEOPOTENTIAL_HEIGHT);
   const Geodetic_Vector_Data_2D& temp_data = *temp_data_ptr;

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real latitude = data_ptr->get_coordinate (0, i);
      const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real longitude = data_ptr->get_coordinate (1, j);

         const Real z_x = temp_data_ptr->evaluate (0, i, j, DX);
         const Real z_y = temp_data_ptr->evaluate (0, i, j, DY);

         const Real u_g = -g * z_y / f;
         const Real v_g = g * z_x / f;

         data_ptr->set_datum (0, i, j, u_g);
         data_ptr->set_datum (1, i, j, v_g);

      }
   }

   delete temp_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_surface_geostrophic_wind_data_ptr (const Key& key)
{

   const Level& msl = Level::mean_sea_level ();
   const Level& screen = Level::screen_level ();

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (4);
   fill_data (*temp_data_ptr, 0, key, msl, PRESSURE);
   fill_data (*temp_data_ptr, 1, key, screen, TEMPERATURE);
   const Geodetic_Vector_Data_2D& temp_data = *temp_data_ptr;

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real latitude = data_ptr->get_coordinate (0, i);
      const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real longitude = data_ptr->get_coordinate (1, j);

         const Real p = temp_data.get_datum (0, i, j);
         const Real p_x = temp_data.evaluate (0, latitude, longitude, DX);
         const Real p_y = temp_data.evaluate (0, latitude, longitude, DY);
         const Real temperature = temp_data.get_datum (1, i, j);
         const Real rho = p / (R_d * temperature);

         const Real u_g = -p_y / (f * rho);
         const Real v_g = p_x / (f * rho);

         data_ptr->set_datum (0, i, j, u_g);
         data_ptr->set_datum (1, i, j, v_g);

      }
   }

   delete temp_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_geostrophic_wind_data_ptr (const Key& key,
                                    const Level& level)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         const Real p = level.value;
         return get_geostrophic_wind_data_ptr (key, p);
      }

      case FIFTY_METRE_LEVEL:
      case TEN_METRE_LEVEL:
      {
         return get_surface_geostrophic_wind_data_ptr (key);
      }

   }

   throw Nwp_Exception ("Geostrophic Wind N/A for this level");

}

Geodetic_Vector_Data_2D*
Nwp::get_wind_data_ptr (const Key& key,
                        const Level& level)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         vector<Nwp_Element> ev;
         ev.push_back (ZONAL_WIND);
         ev.push_back (MERIDIONAL_WIND);
         return get_pressure_level_data_ptr (key, level.value, ev);
      }

      case THETA_LEVEL:
      {
         vector<Nwp_Element> ev;
         ev.push_back (ZONAL_WIND);
         ev.push_back (MERIDIONAL_WIND);
         return get_theta_level_data_ptr (key, level.value, ev, false);
      }

      case SIGMA_LEVEL:
      {
         vector<Nwp_Element> ev;
         ev.push_back (ZONAL_WIND);
         ev.push_back (MERIDIONAL_WIND);
         return get_sigma_level_data_ptr (key, level.value, ev);
      }

      case FIFTY_METRE_LEVEL:
      case TEN_METRE_LEVEL:
      {
         Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, level, ZONAL_WIND);
         fill_data (*data_ptr, 1, key, level, MERIDIONAL_WIND);
         return data_ptr;
      }

   }

   throw Nwp_Exception ("Wind N/A for this level");

}

Geodetic_Vector_Data_2D*
Nwp::get_horizontal_shear_data_ptr (const Key& key,
                                    const Level& level,
                                    const bool with_wind)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* wind_data_ptr = get_wind_data_ptr (key, level);

   const Integer n = (with_wind ? 3 : 1);
   Gvd_2d* data_ptr = get_initialized_vd_2d (n);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real u = wind_data_ptr->get_datum (0, i, j);
         const Real v = wind_data_ptr->get_datum (1, i, j);
         const Real u_x = wind_data_ptr->evaluate (0, i, j, DX);
         const Real u_y = wind_data_ptr->evaluate (0, i, j, DY);
         const Real v_x = wind_data_ptr->evaluate (1, i, j, DX);
         const Real v_y = wind_data_ptr->evaluate (1, i, j, DY);
         const Real V_2 = u*u + v*v;

         const Real datum = (v*v*v_x - u*u*u_y - u*v*(v_y-u_x)) / V_2;
         data_ptr->set_datum (0, i, j, datum);

         if (with_wind)
         {
            data_ptr->set_datum (1, i, j, u);
            data_ptr->set_datum (2, i, j, v);
         }

      }
   }

   delete wind_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_omega_data_ptr (const Key& key,
                         const Level& level,
                         const bool with_wind)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;

         Integer n = 1;
         if (with_wind) { n += 2; }
         Gvd_2d* data_ptr = get_initialized_vd_2d (n);
         const Real p = level.value;
         fill_pressure_level_data (*data_ptr, 0, key, p, OMEGA);

         if (with_wind)
         {
            fill_wind_data (*data_ptr, 1, 2, key, level);
         }

         return data_ptr;

      }

   }

   throw Nwp_Exception ("Omega N/A for this level");

}

Geodetic_Vector_Data_2D*
Nwp::get_temperature_24hr_data_ptr (const Key& key,
                                    const Level& level)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;
   const Key prev_key (base_time, forecast_hour - 24);

   Gvd_2d* this_data_ptr = get_temperature_data_ptr (key, level);
   Gvd_2d* prev_data_ptr = get_temperature_data_ptr (prev_key, level);

   const Size_2D& size_2d = this_data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real this_datum = this_data_ptr->get_datum (0, i, j);
         const Real prev_datum = prev_data_ptr->get_datum (0, i, j);
         const Real datum = this_datum - prev_datum;
         this_data_ptr->set_datum (0, i, j, datum);
      }
   }

   delete prev_data_ptr;
   return this_data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_mix_down_temperature_data_ptr (const Key& key,
                                        const Level& level,
                                        const bool with_wind)
{


   Integer n = 1;
   if (with_wind) { n += 2; }

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = get_initialized_vd_2d (n);
   fill_mix_down_temperature_data (*data_ptr, 0, key, level);

   if (with_wind)
   {
      fill_wind_data (*data_ptr, 1, 2, key, level);
   }

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_stratus_data_ptr (const Key& key,
                           const Level& level,
                           const bool with_wind)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   Integer n = 3;
   if (with_wind) { n += 2; }

   Gvd_2d* data_ptr = get_initialized_vd_2d (n);
   const Level& screen = Level::screen_level ();
   Gvd_2d* next_p_data_ptr = get_next_p_up_data_ptr (key, screen);
   Gvd_2d* this_t_data_ptr = get_temperature_data_ptr (key, screen);
   Gvd_2d* this_td_data_ptr = get_dew_point_data_ptr (key, screen);

   Gvd_2d& gvd_2d = *data_ptr;
   const Data_3D& data_3d = get_3d_data (key);
   const Tuple& tuple_p = data_3d.get_tuple_p (TEMPERATURE);
   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real next_p = next_p_data_ptr->get_datum (0, i, j);
         const Real next_t = data_3d.evaluate (TEMPERATURE, next_p, i, j);
         const Real next_td = data_3d.evaluate (DEW_POINT, next_p, i, j);
         const Real this_t = this_t_data_ptr->get_datum (0, i, j);
         const Real this_td = this_td_data_ptr->get_datum (0, i, j);
         const Real next_ttd = next_t - next_td;
         const Real this_ttd = this_t - this_td;
         const Real stability = next_t - this_t;
         gvd_2d.set_datum (0, i, j, this_ttd);
         gvd_2d.set_datum (1, i, j, next_ttd);
         gvd_2d.set_datum (2, i, j, stability);
      }
   }

   delete next_p_data_ptr;
   delete this_t_data_ptr;
   delete this_td_data_ptr;

   if (with_wind)
   {
      fill_wind_data (*data_ptr, 2, 3, key, level);
   }

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_lapse_data_ptr (const Key& key,
                         const Level& level,
                         const bool with_wind)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   Integer n = 1;
   if (with_wind) { n += 2; }
   Gvd_2d* data_ptr = get_initialized_vd_2d (n);
   fill_lapse_data (*data_ptr, 0, key, level);

   if (with_wind)
   {
      fill_wind_data (*data_ptr, 1, 2, key, level);
   }

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_temperature_data_ptr (const Key& key,
                               const Level& level,
                               const bool with_wind)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {

         vector<Nwp_Element> ev;
         ev.push_back (TEMPERATURE);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case THETA_LEVEL:
      {

         vector<Nwp_Element> ev;
         ev.push_back (TEMPERATURE);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case SCREEN_LEVEL:
      {
         if (with_wind)
         {
            const Level& screen = Level::screen_level ();
            const Level& ten = Level::ten_metre_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
            fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
            fill_data (*data_ptr, 1, key, ten, ZONAL_WIND);
            fill_data (*data_ptr, 2, key, ten, MERIDIONAL_WIND);
            return data_ptr;
         }
         else
         {
            const Level& screen = Level::screen_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
            fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
            return data_ptr;
         }
      }

   }

   throw Nwp_Exception ("Temperature N/A for this level");

}

Geodetic_Vector_Data_2D*
Nwp::get_dew_point_data_ptr (const Key& key,
                             const Level& level,
                             const bool with_wind)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         vector<Nwp_Element> ev;
         ev.push_back (DEW_POINT);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case THETA_LEVEL:
      {

         vector<Nwp_Element> ev;
         ev.push_back (DEW_POINT);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case SCREEN_LEVEL:
      {
         if (with_wind)
         {
            const Level& screen = Level::screen_level ();
            const Level& ten = Level::ten_metre_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
            fill_data (*data_ptr, 0, key, screen, DEW_POINT);
            fill_data (*data_ptr, 1, key, ten, ZONAL_WIND);
            fill_data (*data_ptr, 2, key, ten, MERIDIONAL_WIND);
            return data_ptr;
         }
         else
         {
            const Level& screen = Level::screen_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
            fill_data (*data_ptr, 0, key, screen, DEW_POINT);
            return data_ptr;
         }
      }

   }

   throw Nwp_Exception ("Dew Point N/A for this level");

}

Geodetic_Vector_Data_2D*
Nwp::get_rh_data_ptr (const Key& key,
                      const Level& level,
                      const bool with_wind)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         vector<Nwp_Element> ev;
         ev.push_back (RELATIVE_HUMIDITY);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case THETA_LEVEL:
      {

         vector<Nwp_Element> ev;
         ev.push_back (RELATIVE_HUMIDITY);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case SCREEN_LEVEL:
      {
         if (with_wind)
         {
            const Level& screen = Level::screen_level ();
            const Level& ten = Level::ten_metre_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
            fill_data (*data_ptr, 0, key, screen, RELATIVE_HUMIDITY);
            fill_data (*data_ptr, 1, key, ten, ZONAL_WIND);
            fill_data (*data_ptr, 2, key, ten, MERIDIONAL_WIND);
            return data_ptr;
         }
         else
         {
            const Level& screen = Level::screen_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
            fill_data (*data_ptr, 0, key, screen, RELATIVE_HUMIDITY);
            return data_ptr;
         }
      }

   }

   throw Nwp_Exception ("Relative Humidity N/A for this level");

}

Geodetic_Vector_Data_2D*
Nwp::get_dew_point_depression_data_ptr (const Key& key,
                                        const Level& level,
                                        const bool with_wind)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         vector<Nwp_Element> ev;
         ev.push_back (DEW_POINT_DEPRESSION);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case THETA_LEVEL:
      {

         vector<Nwp_Element> ev;
         ev.push_back (DEW_POINT_DEPRESSION);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case SCREEN_LEVEL:
      {
         if (with_wind)
         {
            const Level& screen = Level::screen_level ();
            const Level& ten = Level::ten_metre_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
            fill_data (*data_ptr, 0, key, screen, DEW_POINT_DEPRESSION);
            fill_data (*data_ptr, 1, key, ten, ZONAL_WIND);
            fill_data (*data_ptr, 2, key, ten, MERIDIONAL_WIND);
            return data_ptr;
         }
         else
         {
            const Level& screen = Level::screen_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
            fill_data (*data_ptr, 0, key, screen, DEW_POINT_DEPRESSION);
            return data_ptr;
         }
      }

   }

   throw Nwp_Exception ("Dew Point Depression N/A for this level");

}

Geodetic_Vector_Data_2D*
Nwp::get_theta_e_data_ptr (const Key& key,
                           const Level& level,
                           const bool with_wind)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         vector<Nwp_Element> ev;
         ev.push_back (THETA_E);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case THETA_LEVEL:
      {

         vector<Nwp_Element> ev;
         ev.push_back (THETA_E);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case SCREEN_LEVEL:
      {
         if (with_wind)
         {
            const Level& screen = Level::screen_level ();
            const Level& ten = Level::ten_metre_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
            fill_data (*data_ptr, 0, key, screen, THETA_E);
            fill_data (*data_ptr, 1, key, ten, ZONAL_WIND);
            fill_data (*data_ptr, 2, key, ten, MERIDIONAL_WIND);
            return data_ptr;
         }
         else
         {
            const Level& screen = Level::screen_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
            fill_data (*data_ptr, 0, key, screen, THETA_E);
            return data_ptr;
         }
      }

   }

   throw Nwp_Exception ("Theta_e N/A for this level");

}

Geodetic_Vector_Data_2D*
Nwp::get_theta_w_data_ptr (const Key& key,
                           const Level& level,
                           const bool with_wind)
{

   switch (level.type)
   {

      case PRESSURE_LEVEL:
      {
         vector<Nwp_Element> ev;
         ev.push_back (THETA_W);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case THETA_LEVEL:
      {

         vector<Nwp_Element> ev;
         ev.push_back (THETA_W);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case SCREEN_LEVEL:
      {
         if (with_wind)
         {
            const Level& screen = Level::screen_level ();
            const Level& ten = Level::ten_metre_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
            fill_data (*data_ptr, 0, key, screen, THETA_W);
            fill_data (*data_ptr, 1, key, ten, ZONAL_WIND);
            fill_data (*data_ptr, 2, key, ten, MERIDIONAL_WIND);
            return data_ptr;
         }
         else
         {
            const Level& screen = Level::screen_level ();
            Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
            fill_data (*data_ptr, 0, key, screen, THETA_W);
            return data_ptr;
         }
      }

   }

   throw Nwp_Exception ("Theta_w N/A for this level");

}

Geodetic_Vector_Data_2D*
Nwp::get_q_vector_data_ptr (const Key& key,
                            const Level& level)
{

   if (level.type != PRESSURE_LEVEL)
   {
      throw Nwp_Exception ("Q Vector only available for pressure levels");
   }

   const Real p = level.value;

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (2);
   fill_data (*temp_data_ptr, 0, key, level, TEMPERATURE);
   fill_data (*temp_data_ptr, 1, key, level, GEOPOTENTIAL_HEIGHT);

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
   const Size_2D& size_2d = data_ptr->get_size_2d ();
   const Real Rg_p = R * g / p;

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real latitude = data_ptr->get_coordinate (0, i);
      const Real f = Geodetic_Vector_Data_2D::get_f (latitude);
      const Real Rg_pf = Rg_p / f;

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real t_x = temp_data_ptr->evaluate (0, i, j, DX);
         const Real t_y = temp_data_ptr->evaluate (0, i, j, DY);
         const Real z_xx = temp_data_ptr->evaluate (1, i, j, DX2);
         const Real z_xy = temp_data_ptr->evaluate (1, i, j, DXY);
         const Real z_yy = temp_data_ptr->evaluate (1, i, j, DY2);

         const Real q_x = Rg_pf * (z_xy * t_x - z_xx * t_y);
         const Real q_y = Rg_pf * (z_yy * t_x - z_xy * t_y);

         data_ptr->set_datum (0, i, j, q_x);
         data_ptr->set_datum (1, i, j, q_y);

      }
   }

   delete temp_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_temperature_advection_data_ptr (const Key& key,
                                         const Level& level)
{

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (1);
   fill_data (*temp_data_ptr, 0, key, level, TEMPERATURE);

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
   fill_data (*data_ptr, 1, key, level, ZONAL_WIND);
   fill_data (*data_ptr, 2, key, level, MERIDIONAL_WIND);

   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real t_x = temp_data_ptr->evaluate (0, i, j, DX);
         const Real t_y = temp_data_ptr->evaluate (0, i, j, DY);
         const Real u = data_ptr->evaluate (1, i, j);
         const Real v = data_ptr->evaluate (2, i, j);

         const Real temperature_advection = - (t_x * u + t_y * v);
         data_ptr->set_datum (0, i, j, temperature_advection);

      }
   }

   delete temp_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_gd_sea_wave_data_ptr (const Key& key,
                               const Level& level,
                               const Fetch& fetch,
                               const Real max_fetch)
{


   const Real bearing = 90;

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   Geodetic_Vector_Data_2D* wind_data_ptr = get_wind_data_ptr (key, level);

   Lat_Long lat_long;
   Real& latitude = lat_long.latitude;
   Real& longitude = lat_long.longitude;

   Gd_Sea_Wave gsw;

   for (Integer i = 0; i < size_2d.i; i++)
   {
      latitude = data_ptr->get_coordinate (0, i);
      for (Integer j = 0; j < size_2d.j; j++)
      {

         longitude = data_ptr->get_coordinate (1, j);

         const Real u = wind_data_ptr->get_datum (0, i, j);
         const Real v = wind_data_ptr->get_datum (1, i, j);
         const Wind wind (u, v);
         const Real bearing = wind.get_direction ();

         const Real raw_fetch = fetch.get_fetch (bearing, lat_long);

         if (raw_fetch < 1) { data_ptr->set_datum (0, i, j, 0); }
         else
         {

            const Real speed = wind.get_speed ();
            const Real f = min (raw_fetch, max_fetch);
            const Real d = (speed < 15 ? 7 : 6);

            const Real gd_sea_wave_f = gsw.get_seas_from_fetch (speed, f);
            const Real gd_sea_wave_d = gsw.get_seas_from_duration (speed, d);

            const Real gd_sea_wave = std::min (gd_sea_wave_f, gd_sea_wave_d);
            data_ptr->set_datum (0, i, j, gd_sea_wave);

         }

         data_ptr->set_datum (1, i, j, u);
         data_ptr->set_datum (2, i, j, v);
         
      }
   }

   delete wind_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_step_rainfall_data_ptr (const Key& key)
{
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
   fill_rain_data (*data_ptr, 0, key, RAINFALL_STEP);
   return data_ptr;
}

Geodetic_Vector_Data_2D*
Nwp::get_rainfall_data_ptr (const Key& key,
                            const Integer hours)
{

   if (hours < 0)
   {
      return get_step_rainfall_data_ptr (key);
   }
   else
   {
      Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
      fill_rain_data (*data_ptr, 0, key, hours);
      return data_ptr;
   }


}

Geodetic_Vector_Data_2D*
Nwp::get_sutcliffe_data_ptr (const Key& key)
{

   const Level& msl = Level::mean_sea_level ();
   const Level& nil = Level::nil_level ();
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);

   fill_data (*data_ptr, 0, key, msl, PRESSURE);
   fill_data (*data_ptr, 1, key, nil, THICKNESS);
   fill_data (*data_ptr, 2, key, nil, RAINFALL_STEP);

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_pv_p_data_ptr (const Key& key,
                        const Real pv_threshold)
{
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
   fill_pv_p_data (*data_ptr, 0, key, pv_threshold);
   return data_ptr;
}

Geodetic_Vector_Data_2D*
Nwp::get_theta_level_data_ptr (const Key& key,
                               const Real theta,
                               const vector<Nwp_Element> element_vector,
                               const bool with_pv)
{


   const Integer nn = element_vector.size ();
   const Integer n = (with_pv ? nn + 4 : nn);

   const Data_3D& data_3d = get_3d_data (key);

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (n);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   const Geodetic_Vector_Data_3D& t_data_3d = data_3d.get_gvd_3d (TEMPERATURE);
   const Tuple& tuple_p = data_3d.get_tuple_p (TEMPERATURE);

   for (Integer i = 0; i < size_2d.i; i++)
   {

      Thermo_Point thermo_point;
      const Real latitude = data_ptr->get_coordinate (0, i);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         Real p = GSL_NAN;
         const Real longitude = data_ptr->get_coordinate (1, j);

         for (Integer k = 0; k < tuple_p.size () - 1; k++)
         {

            const Real lower_p = tuple_p[k]; 
            const Real upper_p = tuple_p[k + 1];
            const Real lower_t = t_data_3d.get_datum (0, k, i, j); 
            const Real upper_t = t_data_3d.get_datum (0, k + 1, i, j);

            thermo_point.set_t_p (lower_t - K, lower_p);
            const Real lower_theta = thermo_point.get_theta () + K;
            thermo_point.set_t_p (upper_t - K, upper_p);
            const Real upper_theta = thermo_point.get_theta () + K;

            const Real lower_d_theta = theta - lower_theta;
            const Real upper_d_theta = theta - upper_theta;

            if (lower_d_theta * upper_d_theta <= 0)
            {

               const Real d_p = upper_p - lower_p;
               const Real d_theta = upper_theta - lower_theta;

               p = lower_p + (lower_d_theta * d_p / d_theta);
               break;

            }

         }

         if (with_pv)
         {

            const Nwp_Element& T = TEMPERATURE;
            const Nwp_Element& U = ZONAL_WIND;
            const Nwp_Element& V = MERIDIONAL_WIND;

            const Real u = data_3d.evaluate (U, p, latitude, longitude);
            const Real v = data_3d.evaluate (V, p, latitude, longitude);
            const Real dt_dp = data_3d.evaluate (T, p, latitude, longitude, DZ);

            data_ptr->set_datum (0, i, j, p);
            data_ptr->set_datum (1, i, j, u);
            data_ptr->set_datum (2, i, j, v);
            data_ptr->set_datum (3, i, j, dt_dp);

         }

         for (Integer ee = 0; ee < nn; ee++)
         {

            const Nwp_Element& ne = element_vector[ee];
            const Integer e = (with_pv ? ee + 3 : ee);

            Real datum = p;

            if (ne != denise::PRESSURE)
            {
               datum = data_3d.evaluate (ne, p, latitude, longitude);
            }

            data_ptr->set_datum (e, i, j, datum);

         }

      }
   }

   // Calculates pv, before this point slot 3 is dt_dp instead of pv
   if (with_pv)
   {

      for (Integer i = 0; i < size_2d.i; i++)
      {

         const Real latitude = data_ptr->get_coordinate (0, i);
         const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

         for (Integer j = 0; j < size_2d.j; j++)
         {

            const Real longitude = data_ptr->get_coordinate (1, j);

            const Real p = data_ptr->get_datum (0, i, j);
            const Real dv_dx = data_ptr->evaluate (2, latitude, longitude, DX);
            const Real du_dy = data_ptr->evaluate (1, latitude, longitude, DY);
            const Real dt_dp = data_ptr->get_datum (3, i, j);

            const Real exner = pow (Real (p / 1000e2), Real (kappa));
            const Real dtheta_dp = dt_dp / exner - kappa * theta / p;
            const Real pv = -g * (f + dv_dx - du_dy) * dtheta_dp;

            data_ptr->set_datum (3, i, j, pv);

         }

      }

   }

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_sigma_level_data_ptr (const Key& key,
                               const Real sigma,
                               const vector<Nwp_Element> element_vector)
{

   const Integer n = element_vector.size ();
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (n);

   const Level& surface = Level::surface_level ();
   Geodetic_Vector_Data_2D* surface_p_data_ptr = get_initialized_vd_2d (1);
   fill_data (*surface_p_data_ptr, 0, key, surface, PRESSURE);
   const Geodetic_Vector_Data_2D& surf_p_data = *surface_p_data_ptr;

   for (Integer e = 0; e < n; e++)
   {
      const Nwp_Element& ne = element_vector[e];
      fill_sigma_level_data (*data_ptr, e, key, sigma, ne, surf_p_data);
   }

   delete surface_p_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_pressure_level_data_ptr (const Key& key,
                                  const Real p,
                                  const vector<Nwp_Element> element_vector)
{

   const Integer n = element_vector.size ();

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (n);

   for (Integer e = 0; e < n; e++)
   {
      const Nwp_Element& nwp_element = element_vector[e];
      fill_pressure_level_data (*data_ptr, e, key, p, nwp_element);
   }

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_data_ptr (const Key& key,
                   const Level& level,
                   const Nwp_Element nwp_element)
{
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
   fill_data (*data_ptr, 0, key, level, nwp_element);
   return data_ptr;
}

Geodetic_Vector_Data_2D*
Nwp::get_data_ptr (const Key& key,
                   const Level& level,
                   const vector<Nwp_Element>& nwp_element_vector)
{

   const Integer n = nwp_element_vector.size ();
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (n);

   for (Integer i = 0; i < n; i++)
   {
      const Nwp_Element& nwp_element = nwp_element_vector[i];
      fill_data (*data_ptr, i, key, level, nwp_element);
   }

   return data_ptr;

}

Scalar_Data_1D*
Nwp::get_terrain_profile_ptr (const Key& key,
                              const Multi_Journey& multi_journey)
{

   const Geodesy geodesy;
   const Tuple& tuple_x = multi_journey.get_tuple_x (geodesy);
   const Real distance = tuple_x.back ();
   const Level& surface = Level::surface_level ();

   if (tuple_x.size () < 2) { throw Nwp_Exception ("Invalid multi_journey"); }
   if (gsl_isnan (distance)) { throw Nwp_Exception ("Invalid multi_journey"); }
   if (distance < 1) { throw Nwp_Exception ("multi_journey too short"); }

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (1);
   fill_data (*temp_data_ptr, 0, key, surface, PRESSURE);

   Scalar_Data_1D* data_ptr = new Scalar_Data_1D (tuple_x);

   for (Multi_Journey::const_iterator iterator = multi_journey.begin ();
        iterator != multi_journey.end (); iterator++)
   {

      const Integer i = std::distance (multi_journey.begin (), iterator);
      const Real latitude = iterator->x;
      const Real longitude = iterator->y;

      const Real datum = temp_data_ptr->evaluate (0, latitude, longitude);
      data_ptr->set_datum (i, datum);

   }

   delete temp_data_ptr;

   data_ptr->set_interpolation ();
   return data_ptr;

}

Scalar_Data_1D*
Nwp::get_rainfall_profile_ptr (const Key& key,
                               const Multi_Journey& multi_journey)
{

   const Geodesy geodesy;
   const Tuple& tuple_x = multi_journey.get_tuple_x (geodesy);
   const Real distance = tuple_x.back ();

   if (tuple_x.size () < 2) { throw Nwp_Exception ("Invalid multi_journey"); }
   if (gsl_isnan (distance)) { throw Nwp_Exception ("Invalid multi_journey"); }
   if (distance < 1) { throw Nwp_Exception ("multi_journey too short"); }

   const Level& nil = Level::nil_level ();
   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (1);
   fill_data (*temp_data_ptr, 0, key, nil, RAINFALL_STEP);

   Scalar_Data_1D* data_ptr = new Scalar_Data_1D (tuple_x);

   for (Multi_Journey::const_iterator iterator = multi_journey.begin ();
        iterator != multi_journey.end (); iterator++)
   {

      const Integer i = std::distance (multi_journey.begin (), iterator);
      const Real latitude = iterator->x;
      const Real longitude = iterator->y;

      const Real datum = temp_data_ptr->evaluate (0, latitude, longitude);
      data_ptr->set_datum (i, datum);

   }

   delete temp_data_ptr;

   data_ptr->set_interpolation ();
   return data_ptr;

}

Nwp::Cross_Section*
Nwp::get_cross_section_ptr (const Key& key,
                            const Multi_Journey& multi_journey,
                            const Nwp_Element nwp_element,
                            const bool with_wind)
{
   vector<Nwp_Element> nwp_element_vector;
   nwp_element_vector.push_back (nwp_element);
   return get_cross_section_ptr (key, multi_journey,
      nwp_element_vector, with_wind);
}

Nwp::Cross_Section*
Nwp::get_cross_section_ptr (const Key& key,
                            const Multi_Journey& multi_journey,
                            const vector<Nwp_Element>& nwp_element_vector,
                            const bool with_wind)
{

   const Geodesy geodesy;
   const Tuple& tuple_x = multi_journey.get_tuple_x (geodesy);
   const Real distance = tuple_x.back ();

   if (tuple_x.size () < 2) { throw Nwp_Exception ("Invalid multi_journey"); }
   if (gsl_isnan (distance)) { throw Nwp_Exception ("Invalid multi_journey"); }
   if (distance < 1) { throw Nwp_Exception ("multi_journey too short"); }

   const Data_3D& data_3d = get_3d_data (key);

   const Integer nn = nwp_element_vector.size ();
   const Integer nn_1 = nn + 1;
   const Integer nn_2 = nn + 2;

   const Integer n = (with_wind ? nn + 3 : nn);
   Cross_Section* cs_ptr = new Cross_Section ();

   typedef vector<Nwp_Element>::const_iterator Iterator;
   for (Iterator iterator = nwp_element_vector.begin ();
        iterator != nwp_element_vector.end (); iterator++)
   {

      const Nwp_Element& nwp_element = *(iterator);
      const Tuple& tuple_p = data_3d.get_tuple_p (nwp_element);

      for (Multi_Journey::const_iterator i = multi_journey.begin ();
           i != multi_journey.end (); i++)
      {

         const Integer ii = std::distance (multi_journey.begin (), i);
         const Lat_Long& lat_long = *(i);
         const Real latitude = lat_long.latitude;
         const Real longitude = lat_long.longitude;
         const Real angle = multi_journey.get_azimuth_forward (i, geodesy);
         const Real theta = angle * DEGREE_TO_RADIAN;
         const Real c = cos (theta);
         const Real s = sin (theta);

         // This is only computed for ne = LI_THUNDER
         const Real thunder_t = -20 + K;
         const Real thunder_p = data_3d.get_p_from_element (
            TEMPERATURE, latitude, longitude, thunder_t);

         cs_ptr->insert_nwp_element_if_needed (nwp_element, tuple_x, tuple_p);
         Scalar_Data_2D& sd_2d = cs_ptr->get_sd_2d (nwp_element);

         for (Integer k = 0; k < tuple_p.size (); k++)
         {

            const Real p = tuple_p[k];

            try
            {

               const Real datum = (nwp_element == LI_THUNDER) ?
                  data_3d.get_li_thunder (p, lat_long, thunder_p, thunder_t) :
                  data_3d.evaluate (nwp_element, p, lat_long);
               sd_2d.set_datum (ii, k, datum);
            }
            catch (const std::exception& se)
            {
               sd_2d.set_datum (ii, k, GSL_NAN);
            }

         }

      }

   }

   if (with_wind)
   {

      const Nwp_Element U = ZONAL_WIND;
      const Nwp_Element V = MERIDIONAL_WIND;
      const Nwp_Element O = OMEGA;

      const Tuple& tuple_p = data_3d.get_tuple_p (ZONAL_WIND);
      cs_ptr->insert_nwp_element_if_needed (STREAMLINE_WIND, tuple_x, tuple_p);
      cs_ptr->insert_nwp_element_if_needed (NORMAL_WIND, tuple_x, tuple_p);
      cs_ptr->insert_nwp_element_if_needed (OMEGA, tuple_x, tuple_p);
      Scalar_Data_2D& s_sd_2d = cs_ptr->get_sd_2d (STREAMLINE_WIND);
      Scalar_Data_2D& n_sd_2d = cs_ptr->get_sd_2d (NORMAL_WIND);
      Scalar_Data_2D& o_sd_2d = cs_ptr->get_sd_2d (OMEGA);

      for (Multi_Journey::const_iterator i = multi_journey.begin ();
           i != multi_journey.end (); i++)
      {

         const Integer ii = std::distance (multi_journey.begin (), i);
         const Lat_Long& lat_long = *(i);
         const Real latitude = lat_long.latitude;
         const Real longitude = lat_long.longitude;
         const Real angle = multi_journey.get_azimuth_forward (i, geodesy);
         const Real theta = angle * DEGREE_TO_RADIAN;
         const Real c = cos (theta);
         const Real s = sin (theta);

         for (Integer k = 0; k < tuple_p.size (); k++)
         {

            try
            {
               const Real u = data_3d.evaluate (U, k, latitude, longitude);
               const Real v = data_3d.evaluate (V, k, latitude, longitude);
               const Real uu = u * s + v * c;
               const Real vv = v * s - u * c;
               s_sd_2d.set_datum (ii, k, uu);
               n_sd_2d.set_datum (ii, k, vv);
            }
            catch (const std::exception& se)
            {
               s_sd_2d.set_datum (ii, k, GSL_NAN);
               n_sd_2d.set_datum (ii, k, GSL_NAN);
            }

            try
            {
               const Real o = data_3d.evaluate (O, k, latitude, longitude);
               o_sd_2d.set_datum (ii, k, o);
            }
            catch (const std::exception& se)
            {
               o_sd_2d.set_datum (ii, k, 0);
            }

         }

      }

   }

   typedef Scalar_Data_1D Sd_1d;
   Sd_1d* terrain_profile_ptr = get_terrain_profile_ptr (key, multi_journey);
   Sd_1d* rainfall_profile_ptr = get_rainfall_profile_ptr (key, multi_journey);
   cs_ptr->set_profile_ptrs (terrain_profile_ptr, rainfall_profile_ptr);
   return cs_ptr;

}


Nwp::Sounding*
Nwp::get_sounding_ptr (const Lat_Long& lat_long,
                       const Key& key)
{

   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;
   const Level screen = Level::screen_level ();
   const Level surface = Level::surface_level ();
   const Level ten = Level::ten_metre_level ();

   Nwp::Sounding* sounding_ptr = new Nwp::Sounding (key);

   Real surface_p = GSL_POSINF;
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (5);

   try
   {

      fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
      fill_data (*data_ptr, 1, key, screen, DEW_POINT);
      fill_data (*data_ptr, 2, key, surface, PRESSURE);
      fill_data (*data_ptr, 3, key, ten, ZONAL_WIND);
      fill_data (*data_ptr, 4, key, ten, MERIDIONAL_WIND);

      surface_p = data_ptr->evaluate (2, latitude, longitude);
      const Real screen_t = data_ptr->evaluate (0, latitude, longitude);
      const Real screen_t_d = data_ptr->evaluate (1, latitude, longitude);
      const Real u = data_ptr->evaluate (3, latitude, longitude);
      const Real v = data_ptr->evaluate (4, latitude, longitude);
      const Wind wind (u, v);

      sounding_ptr->get_t_line ().add (surface_p, screen_t - K);
      sounding_ptr->get_t_d_line ().add (surface_p, screen_t_d - K);
      sounding_ptr->get_height_profile ().add (surface_p, 0);
      sounding_ptr->get_wind_profile ().add (surface_p, wind);

   }
   catch (const Nwp_Exception& be)
   {
      cerr << "Nwp::get_sounding_ptr " << be << endl;
   }

   delete data_ptr;

   const Data_3D& data_3d = get_3d_data (key);

   typedef Geodetic_Vector_Data_3D Gvd_3d;
   const Tuple& t_tuple_p = data_3d.get_tuple_p (TEMPERATURE);
   const Tuple& td_tuple_p = data_3d.get_tuple_p (DEW_POINT);
   const Tuple& z_tuple_p = data_3d.get_tuple_p (GEOPOTENTIAL_HEIGHT);
   const Tuple& u_tuple_p = data_3d.get_tuple_p (ZONAL_WIND);
   const Tuple& v_tuple_p = data_3d.get_tuple_p (MERIDIONAL_WIND);

   for (Integer k = 0; k < t_tuple_p.size (); k++)
   {
      const Real p = t_tuple_p[k];
      if (p > surface_p) { continue; }
      const Real t = data_3d.evaluate (TEMPERATURE, k, latitude, longitude);
      const bool invalid = gsl_isnan (t) || (t > 350) || (t < -150);
      if (!invalid) { sounding_ptr->get_t_line ().add (p, t - K); }
   }

   for (Integer k = 0; k < td_tuple_p.size (); k++)
   {
      const Real p = td_tuple_p[k];
      if (p > surface_p) { continue; }
      const Real t_d = data_3d.evaluate (DEW_POINT, k, latitude, longitude);
      const bool invalid = gsl_isnan (t_d) || (t_d > 350) || (t_d < -150);
      if (!invalid) { sounding_ptr->get_t_d_line ().add (p, t_d - K); }
   }

   for (Integer k = 0; k < z_tuple_p.size (); k++)
   {
      const Real p = z_tuple_p[k];
      if (p > surface_p) { continue; }
      const Real z = data_3d.evaluate (GEOPOTENTIAL_HEIGHT, k, latitude, longitude);
      const bool invalid = gsl_isnan (z) || (z > 30000) || (z < -500);
      if (!invalid) { sounding_ptr->get_height_profile ().add (p, z); }
   }

   for (Integer k = 0; k < u_tuple_p.size (); k++)
   {
      const Real p = u_tuple_p[k];
      if (p > surface_p) { continue; }
      const Real u = data_3d.evaluate (ZONAL_WIND, k, latitude, longitude);
      const Real v = data_3d.evaluate (MERIDIONAL_WIND, k, latitude, longitude);
      const bool invalid_u = gsl_isnan (u) || (u > 300) || (u < -300);
      const bool invalid_v = gsl_isnan (v) || (v > 300) || (v < -300);
      const bool valid = !invalid_u && !invalid_v;
      if (valid) { sounding_ptr->get_wind_profile ().add (p, Wind (u, v)); }

   }

   return sounding_ptr;

}

Tuple
Nwp::get_point_tuple (const Lat_Long& lat_long,
                      const Key& key)
{

   Tuple tuple;

   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;
   const Level& msl = Level::mean_sea_level ();

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (4);

   const Level screen = Level::screen_level ();
   const Level hPa850 = Level::pressure_level (850e2);

   fill_data (*data_ptr, 0, key, msl, PRESSURE);
   fill_data (*data_ptr, 1, key, screen, TEMPERATURE);
   fill_data (*data_ptr, 2, key, hPa850, TEMPERATURE);
   fill_data (*data_ptr, 3, key, screen, DEW_POINT);

   const Real mslp = data_ptr->evaluate (0, latitude, longitude);
   const Real screen_t = data_ptr->evaluate (1, latitude, longitude);
   const Real t_850 = data_ptr->evaluate (2, latitude, longitude);
   const Real screen_t_d = data_ptr->evaluate (3, latitude, longitude);

   tuple.push_back (mslp);
   tuple.push_back (screen_t);
   tuple.push_back (t_850);
   tuple.push_back (screen_t_d);

   delete data_ptr;

   return tuple;

}

Vector_Data_1D*
Nwp::get_time_series_ptr (const Lat_Long& lat_long,
                          const Dtime& base_time,
                          const Dtime& start_time,
                          const Dtime& end_time,
                          const Level_Element level_element)
{

   vector<Level_Element> level_element_vector;
   level_element_vector.push_back (level_element);

   return get_time_series_ptr (lat_long, base_time,
      start_time, end_time, level_element_vector);

}

Vector_Data_1D*
Nwp::get_time_series_ptr (const Lat_Long& lat_long,
                          const Dtime& base_time,
                          const Dtime& start_time,
                          const Dtime& end_time,
                          const vector<Level_Element>& level_element_vector)
{

   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;
   const Integer n = level_element_vector.size ();
   const Tuple& t_tuple = get_valid_t_tuple (base_time, start_time, end_time);

   if (t_tuple.size () == 0) { return NULL; }

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
   Vector_Data_1D* time_series_ptr = new Vector_Data_1D (n, t_tuple);

   for (Integer i = 0; i < t_tuple.size (); i++)
   {

      const Dtime& dtime = t_tuple[i];
      const Key& key = get_key (dtime, base_time);

      for (Integer e = 0; e < n; e++)
      {

         const Level_Element& level_element = level_element_vector[e];
         const Level& level = level_element.level;
         const Nwp_Element& nwp_element = level_element.nwp_element;

         try
         {
            fill_data (*data_ptr, 0, key, level, nwp_element);
            const Real datum = data_ptr->evaluate (0, latitude, longitude);
            time_series_ptr->set_datum (e, i, datum);
         }
         catch (const Nwp_Exception& ne)
         {
            time_series_ptr->set_datum (e, i, GSL_NAN);
         }


      }

   }

   delete data_ptr;
   const Integer nn = time_series_ptr->size ();

   if (nn < 2)
   {
      delete time_series_ptr;
      return NULL;
   }

   const gsl_interp_type* interp_type = (nn > 3 ?
      gsl_interp_cspline : gsl_interp_linear);
   time_series_ptr->set_interpolation (interp_type);
   return time_series_ptr;

}

Nwp::Time_Cross*
Nwp::get_time_cross_data_ptr (const Lat_Long& lat_long,
                              const Dtime& base_time,
                              const Dtime& start_time,
                              const Dtime& end_time,
                              const vector<Nwp_Element>& nwp_element_vector)
{

   const Tuple& tuple_t = get_valid_t_tuple (base_time, start_time, end_time);
   if (tuple_t.size () < 2) { return NULL; }

   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;

   typedef Scalar_Data_2D Sd_2d;
   typedef vector<Nwp_Element>::const_iterator Iterator;
   Nwp::Time_Cross* tc_ptr = new Time_Cross ();

   for (Integer i = 0; i < tuple_t.size (); i++)
   {

      const Dtime& dtime = tuple_t[i];
      const Key& key = get_key (dtime, base_time);

      for (Iterator iterator = nwp_element_vector.begin ();
           iterator != nwp_element_vector.end (); iterator++)
      {

         const Data_3D& data_3d = get_3d_data (key);
         const Nwp_Element& nwp_element = *(iterator);

         const Real thunder_t = -20 + K;
         const Real thunder_p = data_3d.get_p_from_element (
            TEMPERATURE, latitude, longitude, thunder_t);

         const Tuple& tuple_p = data_3d.get_tuple_p (nwp_element);
         tc_ptr->insert_nwp_element_if_needed (nwp_element, tuple_t, tuple_p);
         Scalar_Data_2D& sd_2d = tc_ptr->get_sd_2d (nwp_element);

         for (Integer k = 0; k < tuple_p.size (); k++)
         {

            try
            {

               const Real datum = (nwp_element == LI_THUNDER) ?
                  data_3d.get_li_thunder (k, lat_long, thunder_p, thunder_t) :
                  data_3d.evaluate (nwp_element, k, lat_long);
               sd_2d.set_datum (i, k, datum);
            }
            catch (const std::exception& se)
            {
               sd_2d.set_datum (i, k, 0);
            }

         }

      }
   }

   typedef Scalar_Data_1D Sd_1d;
   //Sd_1d* terrain_profile_ptr = get_terrain_profile_ptr (key, multi_journey);
   //Sd_1d* rainfall_profile_ptr = get_rainfall_profile_ptr (key, multi_journey);
   //cs_ptr->set_profile_ptrs (terrain_profile_ptr, rainfall_profile_ptr);
   return tc_ptr;

}

Nwp_Exception::Nwp_Exception (const string& str)
   : Exception ("Nwp_Exception", str)
{
}

Access::Data_3D::Data_3D (const vector<Nwp_Element>& nwp_element_vector,
                          const Key& key)
   : Nwp::Data_3D (nwp_element_vector, key)
{
}

Real
Access::Data_3D::evaluate (const Nwp_Element element,
                           const Real p,
                           const Real latitude,
                           const Real longitude,
                           const Evaluate_Op evaluate_op) const
{

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& RH = denise::RELATIVE_HUMIDITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real rh = Nwp::Data_3D::evaluate (RH, p, latitude, longitude);
         const Thermo_Medium thermo_medium = (t < 0 ? ICE : WATER);
         return Moisture::get_t_d (t, rh, WATER);
      }

      case denise::OMEGA:
      {
         const Nwp_Element& T = TEMPERATURE;
         const Nwp_Element& W = VERTICAL_VELOCITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real w = Nwp::Data_3D::evaluate (W, p, latitude, longitude);
         const Real rho = p / (R_d * t);
         return -rho * g * w;
      }

   }

   return Nwp::Data_3D::evaluate (element, p,
      latitude, longitude, evaluate_op);

}

Grib::Key
Access::get_grib_key (const Key& key,
                      const Nwp_Element nwp_element,
                      const Level& level) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   Grib::Key grib_key;
   set_grib_key (grib_key, nwp_element, base_time, forecast_hour);
   set_grib_key (grib_key, nwp_element, level);

   return grib_key;

}

void
Access::set_grib_key (Grib::Key& grib_key,
                      const Nwp_Element nwp_element,
                      const Dtime& base_time,
                      const Integer forecast_hour) const
{

   uint8_t* buffer = grib_key.buffer;

   const Integer yyyy = base_time.get_year ();
   const Integer mm = base_time.get_month ();
   const Integer dd = base_time.get_day ();
   const Integer HH = base_time.get_hour ();
   const Integer MM = base_time.get_minute ();

   buffer[0] = uint8_t (yyyy / 100) + 1;
   buffer[1] = uint8_t (yyyy % 100);
   buffer[2] = uint8_t (mm);
   buffer[3] = uint8_t (dd);
   buffer[4] = uint8_t (HH);
   buffer[5] = uint8_t (MM);
   buffer[6] = 1; // hour
   buffer[7] = uint8_t (forecast_hour);
   buffer[8] = 0;
   buffer[9] = 0;
   buffer[10] = 0;
   buffer[11] = 0;

   if (nwp_element == RAINFALL_CUMULATIVE)
   {
      buffer[7] = 0;
      buffer[8] = uint8_t (forecast_hour);
      buffer[9] = 4;
      buffer[10] = 0;
      buffer[11] = 1;
   }

}

void
Access::set_grib_key (Grib::Key& grib_key,
                      const Nwp_Element nwp_element,
                      const denise::Level& level) const
{

   uint8_t* buffer = grib_key.buffer;
   const uint8_t n = (omega_as_w ? 135 : 134);

   if (level.type == PRESSURE_LEVEL)
   {
      switch (nwp_element)
      {
         case denise::ZONAL_WIND:          { buffer[12] = 131; break; }
         case denise::MERIDIONAL_WIND:     { buffer[12] = 132; break; }
         case denise::TEMPERATURE:         { buffer[12] = 130; break; }
         case denise::RELATIVE_HUMIDITY:   { buffer[12] = 157; break; }
         case denise::VERTICAL_VELOCITY:   { buffer[12] = n; break; }
         case denise::GEOPOTENTIAL_HEIGHT: { buffer[12] = 156; break; }
      }
      buffer[13] = 100;
      uint16_t hpa = uint16_t (round (level.get_value () * 1e-2)); 
#ifndef WORDS_BIGENDIAN
      swap_endian (&hpa, sizeof (uint16_t));
#endif
      memcpy (buffer + 14, &hpa, sizeof (uint16_t));
   }
   else
   {

      switch (nwp_element)
      {
         case MEAN_SEA_LEVEL_PRESSURE:
            buffer[12] = 151;
            buffer[13] = 102;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case PRESSURE:
            buffer[12] = 134;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case PPTN:
         case RAINFALL_CUMULATIVE:
            buffer[12] = 61;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case HIGH_CLOUD:
            buffer[12] = 186;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case MIDDLE_CLOUD:
            buffer[12] = 187;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case LOW_CLOUD:
            buffer[12] = 188;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case ZONAL_WIND:
            buffer[12] = 165;
            buffer[13] = 105;
            buffer[14] = 0;
            buffer[15] = 10;
            break;
         case MERIDIONAL_WIND:
            buffer[12] = 166;
            buffer[13] = 105;
            buffer[14] = 0;
            buffer[15] = 10;
            break;
         case TEMPERATURE:
            buffer[12] = 167;
            buffer[13] = 105;
            buffer[14] = 0;
            buffer[15] = 2;
            break;
         case DEW_POINT:
            buffer[12] = 168;
            buffer[13] = 105;
            buffer[14] = 0;
            buffer[15] = 2;
            break;
      }

   }

}

void
Access::initialize_3d_data (const Key& key)
{
   typedef Access::Data_3D Ad_3d;
   Ad_3d* ad_3d_ptr = new Ad_3d (nwp_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, ad_3d_ptr));
}

void          
Access::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Nwp_Element>::const_iterator Iterator;
   const Key& key = data_3d.key;

   for (Iterator iterator = nwp_element_vector.begin ();
        iterator != nwp_element_vector.end (); iterator++)
   {
      const Nwp_Element& nwp_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = get_gvd_3d_ptr (nwp_element, key);
      data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
   }

   data_3d.set_available ();

}

Geodetic_Vector_Data_2D*
Access::get_initialized_vd_2d (const Integer vector_size) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   Gvd_2d* data_ptr = new Gvd_2d (vector_size, size_2d, domain_2d, false);
   return data_ptr;

}

void 
Access::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Level& level,
                                const Nwp_Element nwp_element)
{

   const Integer vi = vector_index;

   switch (nwp_element)
   {

      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& surface = Level::surface_level ();
         fill_grib_data (gvd_2d, vi, nwp_element, key, surface);
         return;
      }

   }

   Nwp::fill_ts_diagnosis_data (gvd_2d, vi, key, level, nwp_element);

}

void
Access::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                        const Integer vector_index,
                        const Key& key,
                        const Nwp_Element nwp_element)
{

   switch (nwp_element)
   {

      case RAINFALL_CUMULATIVE:
      {

         if (key.forecast_hour == 0)
         {
            gvd_2d.initialize (vector_index, 0);
            return;
         }

         const Level& nil_level = Level::nil_level ();
         fill_grib_data (gvd_2d, vector_index,
            RAINFALL_CUMULATIVE, key, nil_level);
         return;

      }

   }

   Nwp::fill_rain_data (gvd_2d, vector_index, key, nwp_element);

}

void
Access::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Key& key,
                         const Nwp_Element nwp_element)
{
   const Level& nil_level = Level::nil_level ();
   fill_grib_data (gvd_2d, vector_index, nwp_element, key, nil_level);
}

void
Access::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Nwp_Element nwp_element)
{

   const Level& screen = Level::screen_level ();

   switch (nwp_element)
   {

      case TEMPERATURE:
      case DEW_POINT:
      {
         fill_grib_data (gvd_2d, vector_index, nwp_element, key, screen);
         return;
      }

      case RELATIVE_HUMIDITY:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, DEW_POINT);

         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real t_d = data_ptr->get_datum (1, i, j);
               const Real rh = Moisture::get_rh (t - K, t_d - K);
               gvd_2d.set_datum (vector_index, i, j, rh);
            }
         }

         delete data_ptr;
         return;

      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Access::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                             const Integer vector_index,
                             const Key& key,
                             const Nwp_Element nwp_element)
{

   const Level& ten = Level::ten_metre_level ();

   switch (nwp_element)
   {

      case ZONAL_WIND:
      case MERIDIONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, nwp_element, key, ten);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Access::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Key& key,
                       const Nwp_Element nwp_element)
{

   const Level& msl = Level::mean_sea_level ();

   switch (nwp_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         const Nwp_Element mslp = MEAN_SEA_LEVEL_PRESSURE;
         fill_grib_data (gvd_2d, vector_index, mslp, key, msl);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, nwp_element);

}

void
Access::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Nwp_Element nwp_element)
{

   const Level& surface = Level::surface_level ();

   switch (nwp_element)
   {

      case PRESSURE:
      {
         fill_grib_data (gvd_2d, vector_index, PRESSURE, key, surface);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, nwp_element);

}

Geodetic_Vector_Data_3D*
Access::get_gvd_3d_ptr (const Nwp_Element nwp_element,
                        const Key& key) const
{

   Geodetic_Vector_Data_3D* gvd_3d_ptr =
      new Geodetic_Vector_Data_3D (1, tuple_p, size_2d, domain_2d);
   Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   const bool is_rh = (nwp_element == RELATIVE_HUMIDITY);
   const bool is_zonal_wind = (nwp_element == ZONAL_WIND);
   const bool is_meridional_wind = (nwp_element == MERIDIONAL_WIND);
   const bool is_wind = (is_zonal_wind || is_meridional_wind);

   for (Integer k = 0; k < tuple_p.size (); k++)
   {

      const Real p = tuple_p[k];
      const denise::Level level (PRESSURE_LEVEL, p);
      const Grib::Key& grib_key = get_grib_key (key, nwp_element, level);

      Access::const_iterator iterator = find (grib_key);
      if (iterator == end ())
      {
         cout << key.base_time << " " << key.forecast_hour << " " <<
                 nwp_element << " " << level.get_string () <<
                 " throw exception " << endl;
         throw Nwp_Exception ("Access::gvd_3d_ptr Not Available");
      }
      const Grib& grib = *(iterator->second);

      Geodetic_Vector_Data_2D* grib_data_ptr =
         get_grib_data_ptr (grib, grib_key);

      if (is_rh) { grib_data_ptr->scale_offset (0, 0.01, 0); }

      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Integer jj = (j == size_2d.j - 1 ? j - 1 : j);
            Real datum = grib_data_ptr->get_datum (0, i, jj);
            if (is_wind && fabs (datum) > 500) { datum = 0; }
            gvd_3d.set_datum (0, k, i, j, datum);
         }
      }

      delete grib_data_ptr;

   }

   return gvd_3d_ptr;

}

Geodetic_Vector_Data_2D*
Access::get_grib_data_ptr (const Grib& grib,
                           const Grib::Key& grib_key) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = new Gvd_2d (1, size_2d, domain_2d);
   grib.fill_data (*data_ptr, 0, grib_key);

   const Integer last_j = size_2d.j - 1;

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real datum = data_ptr->get_datum (0, i, 0);
      data_ptr->set_datum (0, i, last_j, datum);
   }

   return data_ptr;

}

void
Access::fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                        const Integer vector_index,
                        const Nwp_Element nwp_element,
                        const Key& key,
                        const Level& level) const
{

   const Grib::Key& grib_key = get_grib_key (key, nwp_element, level);
   const Nwp_Exception exception ("Access::fill_grib_data Not Available");

   Access::const_iterator iterator = find (grib_key);
   if (iterator == end ()) { throw exception; }
   const Grib& grib = *(iterator->second);

   Geodetic_Vector_Data_2D* grib_data_ptr = get_grib_data_ptr (grib, grib_key);

   const bool is_zonal_wind = (nwp_element == ZONAL_WIND);
   const bool is_meridional_wind = (nwp_element == MERIDIONAL_WIND);
   const bool is_wind = (is_zonal_wind || is_meridional_wind);

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Integer jj = (j == size_2d.j - 1 ? j - 1 : j);
         Real datum = grib_data_ptr->get_datum (0, i, jj);
         if (is_wind && fabs (datum) > 500) { datum = 0; }
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete grib_data_ptr;

}

Access::Access (const string& description,
                const string& data_path,
                const string& search_string,
                const bool omega_as_w)
   : Nwp (description, data_path),
     data_path (data_path),
     search_string (search_string),
     omega_as_w (omega_as_w)
{

   this->status = "Unloaded";

   nwp_element_vector.push_back (denise::TEMPERATURE);
   nwp_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   nwp_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::ZONAL_WIND);
   nwp_element_vector.push_back (denise::MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::VERTICAL_VELOCITY);

}

Access::~Access ()
{
   clean_up ();
}

void
Access::survey ()
{

   typedef map<Grib::Key, Grib::Header*> Header_Ptr_Map;

   const vector<string>& dir_listing = get_dir_listing (path, search_string);
   for (vector<string>::const_iterator iterator = dir_listing.begin ();
        iterator != dir_listing.end (); iterator++)
   {

      const string& file_name = *(iterator);
      const string& file_path = path + "/" + file_name;

cout << "Grib file_path " << file_path << endl;
      Grib* grib_ptr = new Grib (file_path);
      const Header_Ptr_Map& header_ptr_map = grib_ptr->get_header_ptr_map ();

      for (Header_Ptr_Map::const_iterator i = header_ptr_map.begin ();
           i != header_ptr_map.end (); i++)
      {
         const Grib::Header& header = *(i->second);
         const Grib::Pds& pds = header.get_pds ();
         const Dtime base_time = pds.get_base_time ();
         const Integer forecast_hour = pds.get_forecast_time ().get_p1 ();

         const Grib::Key grib_key (pds);
         const Nwp::Key key (base_time, forecast_hour);

         key_multimap.add (key);
         insert (make_pair (grib_key, grib_ptr));

      }

      grib_ptr_set.insert (grib_ptr);

   }

   if (size () > 0)
   {

      set<uint16_t> set_p;
      const Grib& grib = *(begin ()->second);
      const Header_Ptr_Map& header_ptr_map = grib.get_header_ptr_map ();

      const Grib::Header& first_header = *(header_ptr_map.begin ()->second);
      const Grib::Gds& gds = first_header.get_gds ();
      size_2d = gds.get_size_2d ();
      const Real latitude_0 = gds.get_int (10, 3) * 1e-3;
      const Real longitude_0 = gds.get_int (13, 3) * 1e-3;
      const Real latitude_1 = gds.get_int (17, 3) * 1e-3;
      const Real longitude_1 = gds.get_int (20, 3) * 1e-3;

      Domain_1D& domain_latitude = domain_2d.domain_x;
      Domain_1D& domain_longitude = domain_2d.domain_y;
      domain_latitude.start = std::min (latitude_0, latitude_1);
      domain_latitude.end = std::max (latitude_0, latitude_1);
      domain_longitude.start = std::min (longitude_0, longitude_1);
      domain_longitude.end = std::max (longitude_0, longitude_1);

      for (Header_Ptr_Map::const_iterator iterator = header_ptr_map.begin ();
           iterator != header_ptr_map.end (); iterator++)
      {

         const Grib::Header& header = *(iterator->second);
         const Grib::Pds& pds = header.get_pds ();
         const Grib::Pds::Level& level = pds.get_level ();

         if (level.get_uint (0, 1) == 100)
         {
            const uint16_t p = level.get_uint (1, 2);
            if (p < 100) { continue; }
            set_p.insert (p);
         }

      }

      for (set<uint16_t>::const_iterator iterator = set_p.begin ();
           iterator != set_p.end (); iterator++)
      {
         const Real p = Real (*(iterator)) * 100;
         tuple_p.push_back (p);
      }

   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

   for (Key_Multimap::const_iterator iterator = key_multimap.begin ();
        iterator != key_multimap.end (); iterator++)
   {
      const Key key (iterator->first, iterator->second);
      initialize_3d_data (key);
   }

}

void
Access::clean_up ()
{

   for (set<Grib*>::iterator iterator = grib_ptr_set.begin ();
        iterator != grib_ptr_set.end (); iterator++)
   {
      Grib* grib_ptr = *(iterator);
      delete grib_ptr;
   }

   key_multimap.clear ();
   tuple_p.clear ();
   clear_data_3d_ptr_map ();

   clear ();
   grib_ptr_set.clear ();

}

Ecmwf::Data_3D::Data_3D (const vector<Nwp_Element>& nwp_element_vector,
                         const Key& key)
   : Nwp::Data_3D (nwp_element_vector, key)
{
}

Real
Ecmwf::Data_3D::evaluate (const Nwp_Element element,
                          const Real p,
                          const Real latitude,
                          const Real longitude,
                          const Evaluate_Op evaluate_op) const
{

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& R = denise::MIXING_RATIO;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real q = Nwp::Data_3D::evaluate (R, p, latitude, longitude);
         const Real rh = 0.00263 * p*q / (exp ((17.67 * (t-K) / (t-29.65))));
         return Moisture::get_t_d (t, rh);
      }

      case denise::RELATIVE_HUMIDITY:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& R = denise::MIXING_RATIO;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real q = Nwp::Data_3D::evaluate (R, p, latitude, longitude);
         return 0.00263 * p * q / (exp ((17.67 * (t - K) / (t - 29.65))));
      }

      case denise::OMEGA:
      {
         const Nwp_Element& T = TEMPERATURE;
         const Nwp_Element& W = VERTICAL_VELOCITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real w = Nwp::Data_3D::evaluate (W, p, latitude, longitude);
         const Real rho = p / (R_d * t);
         return -rho * g * w;
      }

   }

   return Nwp::Data_3D::evaluate (element, p,
      latitude, longitude, evaluate_op);

}

Grib::Key
Ecmwf::get_grib_key (const Key& key,
                     const Nwp_Element nwp_element,
                     const Level& level) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   Grib::Key grib_key;
   set_grib_key (grib_key, nwp_element, base_time, forecast_hour);
   set_grib_key (grib_key, nwp_element, level);

   return grib_key;

}

void
Ecmwf::set_grib_key (Grib::Key& grib_key,
                     const Nwp_Element nwp_element,
                     const Dtime& base_time,
                     const Integer forecast_hour) const
{

   uint8_t* buffer = grib_key.buffer;

   const Integer yyyy = base_time.get_year ();
   const Integer mm = base_time.get_month ();
   const Integer dd = base_time.get_day ();
   const Integer HH = base_time.get_hour ();
   const Integer MM = base_time.get_minute ();

   buffer[0] = uint8_t (yyyy / 100) + 1;
   buffer[1] = uint8_t (yyyy % 100);
   buffer[2] = uint8_t (mm);
   buffer[3] = uint8_t (dd);
   buffer[4] = uint8_t (HH);
   buffer[5] = uint8_t (MM);
   buffer[6] = 1; // hour
   buffer[7] = uint8_t (forecast_hour);
   buffer[8] = 0;
   buffer[9] = 0;
   buffer[10] = 0;
   buffer[11] = 0;

   if (nwp_element == RAINFALL_CUMULATIVE)
   {
      buffer[7] = 0;
      buffer[8] = uint8_t (forecast_hour);
      buffer[9] = 4;
      buffer[10] = 0;
      buffer[11] = 1;
   }

//cout << "Ecmwf::set_grib_key " << grib_key << endl;

}

void
Ecmwf::set_grib_key (Grib::Key& grib_key,
                     const Nwp_Element nwp_element,
                     const denise::Level& level) const
{

   uint8_t* buffer = grib_key.buffer;
   //const uint8_t n = (omega_as_w ? 135 : 134);

   if (level.type == PRESSURE_LEVEL)
   {
      switch (nwp_element)
      {
         case denise::ZONAL_WIND:          { buffer[12] = 131; break; }
         case denise::MERIDIONAL_WIND:     { buffer[12] = 132; break; }
         case denise::TEMPERATURE:         { buffer[12] = 130; break; }
         case denise::MIXING_RATIO:        { buffer[12] = 133; break; }
         case denise::VERTICAL_VELOCITY:   { buffer[12] = 135; break; }
         case denise::GEOPOTENTIAL_HEIGHT: { buffer[12] = 156; break; }
      }
      buffer[13] = 100;
      uint16_t hpa = uint16_t (round (level.get_value () * 1e-2)); 
#ifndef WORDS_BIGENDIAN
      swap_endian (&hpa, sizeof (uint16_t));
#endif
      memcpy (buffer + 14, &hpa, sizeof (uint16_t));
   }
   else
   {

      switch (nwp_element)
      {
         case MEAN_SEA_LEVEL_PRESSURE:
            buffer[12] = 151;
            buffer[13] = 102;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case PRESSURE:
            buffer[12] = 134;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case PPTN:
         case RAINFALL_CUMULATIVE:
            buffer[12] = 61;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case HIGH_CLOUD:
            buffer[12] = 186;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case MIDDLE_CLOUD:
            buffer[12] = 187;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case LOW_CLOUD:
            buffer[12] = 188;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case ZONAL_WIND:
            buffer[12] = 165;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case MERIDIONAL_WIND:
            buffer[12] = 166;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case TEMPERATURE:
            buffer[12] = 167;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case DEW_POINT:
            buffer[12] = 168;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
      }

//cout << "!!Ecmwf::set_grib_key " << grib_key << endl;
   }

}

void
Ecmwf::initialize_3d_data (const Key& key)
{
   typedef Ecmwf::Data_3D Ed_3d;
   Ed_3d* ed_3d_ptr = new Ed_3d (nwp_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, ed_3d_ptr));
}

void          
Ecmwf::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Nwp_Element>::const_iterator Iterator;
   const Key& key = data_3d.key;

   for (Iterator iterator = nwp_element_vector.begin ();
        iterator != nwp_element_vector.end (); iterator++)
   {
      const Nwp_Element& nwp_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = get_gvd_3d_ptr (nwp_element, key);
      data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
   }

   data_3d.set_available ();

}

Geodetic_Vector_Data_2D*
Ecmwf::get_initialized_vd_2d (const Integer vector_size) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   const Size_2D& size_2d = size_2d_map.at (1000);
   Gvd_2d* data_ptr = new Gvd_2d (vector_size, size_2d, domain_2d, false);
   return data_ptr;

}

void 
Ecmwf::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Level& level,
                               const Nwp_Element nwp_element)
{

   const Integer vi = vector_index;

   switch (nwp_element)
   {

      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& surface = Level::surface_level ();
         fill_grib_data (gvd_2d, vi, nwp_element, key, surface);
         return;
      }

   }

   Nwp::fill_ts_diagnosis_data (gvd_2d, vi, key, level, nwp_element);

}

void
Ecmwf::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Key& key,
                       const Nwp_Element nwp_element)
{

   switch (nwp_element)
   {

      case RAINFALL_CUMULATIVE:
      {

         if (key.forecast_hour == 0)
         {
            gvd_2d.initialize (vector_index, 0);
            return;
         }

         const Level& nil_level = Level::nil_level ();
         fill_grib_data (gvd_2d, vector_index,
            RAINFALL_CUMULATIVE, key, nil_level);
         return;

      }

   }

   Nwp::fill_rain_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ecmwf::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                        const Integer vector_index,
                        const Key& key,
                        const Nwp_Element nwp_element)
{
   const Level& nil_level = Level::nil_level ();
   fill_grib_data (gvd_2d, vector_index, nwp_element, key, nil_level);
}

void
Ecmwf::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Nwp_Element nwp_element)
{

   const Level& screen = Level::screen_level ();

   switch (nwp_element)
   {

      case TEMPERATURE:
      case DEW_POINT:
      {
         fill_grib_data (gvd_2d, vector_index, nwp_element, key, screen);
         return;
      }

      case RELATIVE_HUMIDITY:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, DEW_POINT);

         const Size_2D& size_2d = gvd_2d.get_size_2d ();
         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            const Real latitude = gvd_2d.get_coordinate (0, i);
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real longitude = gvd_2d.get_coordinate (1, j);
               const Real t = data_ptr->evaluate (0, latitude, longitude);
               const Real t_d = data_ptr->evaluate (1, latitude, longitude);
               const Real rh = Moisture::get_rh (t - K, t_d - K);
               gvd_2d.set_datum (vector_index, i, j, rh);
            }
         }

         delete data_ptr;
         return;

      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ecmwf::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                            const Integer vector_index,
                            const Key& key,
                            const Nwp_Element nwp_element)
{

   const Level& ten = Level::ten_metre_level ();

   switch (nwp_element)
   {

      case ZONAL_WIND:
      case MERIDIONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, nwp_element, key, ten);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ecmwf::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Nwp_Element nwp_element)
{

   const Level& msl = Level::mean_sea_level ();

   switch (nwp_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         const Nwp_Element mslp = MEAN_SEA_LEVEL_PRESSURE;
         fill_grib_data (gvd_2d, vector_index, mslp, key, msl);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ecmwf::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Nwp_Element nwp_element)
{

   const Level& surface = Level::surface_level ();

   switch (nwp_element)
   {

      case PRESSURE:
      {
         fill_grib_data (gvd_2d, vector_index, PRESSURE, key, surface);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, nwp_element);

}

Geodetic_Vector_Data_3D*
Ecmwf::get_gvd_3d_ptr (const Nwp_Element nwp_element,
                       const Key& key) const
{

   const Size_2D& size_2d = size_2d_map.at (1000);
cout << "GVD 3D " << size_2d << " " << tuple_p << " " << nwp_element << endl;
   Geodetic_Vector_Data_3D* gvd_3d_ptr =
      new Geodetic_Vector_Data_3D (1, tuple_p, size_2d, domain_2d);
   Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   const bool is_rh = (nwp_element == RELATIVE_HUMIDITY);
   const bool is_zonal_wind = (nwp_element == ZONAL_WIND);
   const bool is_meridional_wind = (nwp_element == MERIDIONAL_WIND);
   const bool is_wind = (is_zonal_wind || is_meridional_wind);

   for (Integer k = 0; k < tuple_p.size (); k++)
   {

      const Real p = tuple_p[k];
      const denise::Level level (PRESSURE_LEVEL, p);
      const Grib::Key& grib_key = get_grib_key (key, nwp_element, level);

      Ecmwf::const_iterator iterator = find (grib_key);
      if (iterator == end ())
      {
         //cout << key.base_time << " " << key.forecast_hour << " " <<
         //        nwp_element << " " << level.get_string () <<
         //        " throw exception " << endl;
         throw Nwp_Exception ("Ecmwf::gvd_3d_ptr Not Available");
      }
      const Grib& grib = *(iterator->second);

      Geodetic_Vector_Data_2D* grib_data_ptr =
         get_grib_data_ptr (grib, grib_key);

      if (is_rh) { grib_data_ptr->scale_offset (0, 0.01, 0); }

      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         const Real latitude = gvd_3d_ptr->get_coordinate (1, i);
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real longitude = gvd_3d_ptr->get_coordinate (2, j);
            const Integer jj = (j == size_2d.j - 1 ? j - 1 : j);
            Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
            if (is_wind && fabs (datum) > 500) { datum = 0; }
            gvd_3d.set_datum (0, k, i, j, datum);
         }
      }

      delete grib_data_ptr;

   }

   return gvd_3d_ptr;

}

Geodetic_Vector_Data_2D*
Ecmwf::get_grib_data_ptr (const Grib& grib,
                          const Grib::Key& grib_key) const
{

   //const Size_2D& size_2d = size_2d_map.at (1000);

   const Grib::Gds& gds = grib.get_header_ptr_map ().at (grib_key)->get_gds ();
   const Size_2D& size_2d = gds.get_size_2d ();

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = new Gvd_2d (1, size_2d, domain_2d);
   grib.fill_data (*data_ptr, 0, grib_key);

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   const Integer last_j = size_2d.j - 1;

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real datum = data_ptr->get_datum (0, i, 0);
      data_ptr->set_datum (0, i, last_j, datum);
   }

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   return data_ptr;

}

void
Ecmwf::fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Nwp_Element nwp_element,
                       const Key& key,
                       const Level& level) const
{

   const Grib::Key& grib_key = get_grib_key (key, nwp_element, level);
   const Nwp_Exception exception ("Ecmwf::fill_grib_data Not Available");

   Ecmwf::const_iterator iterator = find (grib_key);
   if (iterator == end ()) { throw exception; }
   const Grib& grib = *(iterator->second);

   Geodetic_Vector_Data_2D* grib_data_ptr = get_grib_data_ptr (grib, grib_key);

   const bool is_zonal_wind = (nwp_element == ZONAL_WIND);
   const bool is_meridional_wind = (nwp_element == MERIDIONAL_WIND);
   const bool is_wind = (is_zonal_wind || is_meridional_wind);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_coordinate (0, i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_coordinate (1, j);
         Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
         if (is_wind && fabs (datum) > 500) { datum = 0; }
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete grib_data_ptr;

}

Ecmwf::Ecmwf (const string& description,
              const string& data_path,
              const string& search_string,
              const bool omega_as_w)
   : Nwp (description, data_path),
     data_path (data_path),
     search_string (search_string),
     omega_as_w (omega_as_w)
{

   this->status = "Unloaded";

   nwp_element_vector.push_back (denise::TEMPERATURE);
   nwp_element_vector.push_back (denise::MIXING_RATIO);
   nwp_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::ZONAL_WIND);
   nwp_element_vector.push_back (denise::MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::VERTICAL_VELOCITY);

}

Ecmwf::~Ecmwf ()
{
   clean_up ();
}

void
Ecmwf::survey ()
{

   typedef map<Grib::Key, Grib::Header*> Header_Ptr_Map;

   const vector<string>& dir_listing = get_dir_listing (path, search_string);
   for (vector<string>::const_iterator iterator = dir_listing.begin ();
        iterator != dir_listing.end (); iterator++)
   {

      const string& file_name = *(iterator);
      const string& file_path = path + "/" + file_name;

      Grib* grib_ptr = new Grib (file_path);
      const Header_Ptr_Map& header_ptr_map = grib_ptr->get_header_ptr_map ();

      for (Header_Ptr_Map::const_iterator i = header_ptr_map.begin ();
           i != header_ptr_map.end (); i++)
      {
         const Grib::Header& header = *(i->second);
         const Grib::Pds& pds = header.get_pds ();
         const Dtime base_time = pds.get_base_time ();
         const Integer forecast_hour = pds.get_forecast_time ().get_p1 ();

         const Grib::Key grib_key (pds);
         const Nwp::Key key (base_time, forecast_hour);

         key_multimap.add (key);
         insert (make_pair (grib_key, grib_ptr));

      }

      grib_ptr_set.insert (grib_ptr);

   }

   if (size () > 0)
   {

      set<uint16_t> set_p;
      const Grib& grib = *(begin ()->second);
      const Header_Ptr_Map& header_ptr_map = grib.get_header_ptr_map ();

      const Grib::Header& first_header = *(header_ptr_map.begin ()->second);
      const Grib::Gds& gds = first_header.get_gds ();
      const Real latitude_0 = gds.get_int (10, 3) * 1e-3;
      const Real longitude_0 = gds.get_int (13, 3) * 1e-3;
      const Real latitude_1 = gds.get_int (17, 3) * 1e-3;
      const Real longitude_1 = gds.get_int (20, 3) * 1e-3;

      Domain_1D& domain_latitude = domain_2d.domain_x;
      Domain_1D& domain_longitude = domain_2d.domain_y;
      domain_latitude.start = std::min (latitude_0, latitude_1);
      domain_latitude.end = std::max (latitude_0, latitude_1);
      domain_longitude.start = std::min (longitude_0, longitude_1);
      domain_longitude.end = std::max (longitude_0, longitude_1);

      for (Header_Ptr_Map::const_iterator iterator = header_ptr_map.begin ();
           iterator != header_ptr_map.end (); iterator++)
      {

         const Grib::Header& header = *(iterator->second);
         const Grib::Pds& pds = header.get_pds ();
         const Grib::Pds::Level& level = pds.get_level ();

         const Grib::Gds& gds = header.get_gds ();

         if (level.get_uint (0, 1) == 100)
         {
            const uint16_t p = level.get_uint (1, 2);
            if (p < 100) { continue; }
            size_2d_map[p] = gds.get_size_2d ();
            set_p.insert (p);
         }

      }

      for (set<uint16_t>::const_iterator iterator = set_p.begin ();
           iterator != set_p.end (); iterator++)
      {
         const Real p = Real (*(iterator)) * 100;
         tuple_p.push_back (p);
      }

   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

   for (Key_Multimap::const_iterator iterator = key_multimap.begin ();
        iterator != key_multimap.end (); iterator++)
   {
      const Key key (iterator->first, iterator->second);
      initialize_3d_data (key);
   }

}

void
Ecmwf::clean_up ()
{

   for (set<Grib*>::iterator iterator = grib_ptr_set.begin ();
        iterator != grib_ptr_set.end (); iterator++)
   {
      Grib* grib_ptr = *(iterator);
      delete grib_ptr;
   }

   key_multimap.clear ();
   tuple_p.clear ();
   clear_data_3d_ptr_map ();

   clear ();
   grib_ptr_set.clear ();

}

Gfs3::Data_3D::Data_3D (const vector<Nwp_Element>& nwp_element_vector,
                        const Key& key)
   : Nwp::Data_3D (nwp_element_vector, key)
{
}

Real
Gfs3::Data_3D::evaluate (const Nwp_Element element,
                         const Real p,
                         const Real latitude,
                         const Real longitude,
                         const Evaluate_Op evaluate_op) const
{

   typedef Nwp::Data_3D Nd_3d;
   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& RH = denise::RELATIVE_HUMIDITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real rh = Nwp::Data_3D::evaluate (RH, p, latitude, longitude);
         return Moisture::get_t_d (t, rh);
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& O = denise::OMEGA;

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Evaluate_Op& eo = evaluate_op;
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real ts = Nd_3d::evaluate (T, p, latitude, longitude, eo);
            const Real os = Nd_3d::evaluate (O, p, latitude, longitude, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nd_3d::evaluate (element, p, latitude, longitude, eo);

}

Grib::Key
Gfs3::get_grib_key (const Key& key,
                    const Nwp_Element nwp_element,
                    const Level& level) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   Grib::Key grib_key;
   set_grib_key (grib_key, base_time);
   set_grib_key (grib_key, nwp_element, forecast_hour);
   set_grib_key (grib_key, nwp_element);
   set_grib_key (grib_key, nwp_element, level);

   return grib_key;

}

void
Gfs3::set_grib_key (Grib::Key& grib_key,
                    const Dtime& base_time) const
{

   uint8_t* buffer = grib_key.buffer;

   const Integer yyyy = base_time.get_year ();
   const Integer mm = base_time.get_month ();
   const Integer dd = base_time.get_day ();
   const Integer HH = base_time.get_hour ();
   const Integer MM = base_time.get_minute ();

   buffer[0] = uint8_t (yyyy / 100) + 1;
   buffer[1] = uint8_t (yyyy % 100);
   buffer[2] = uint8_t (mm);
   buffer[3] = uint8_t (dd);
   buffer[4] = uint8_t (HH);
   buffer[5] = uint8_t (MM);

}

void
Gfs3::set_grib_key (Grib::Key& grib_key,
                    const Nwp_Element nwp_element,
                    const Integer forecast_hour) const
{

   uint8_t* buffer = grib_key.buffer;

   buffer[6] = 1;
   buffer[10] = 0;
   buffer[11] = 0;

   switch (nwp_element)
   {

      default:
         buffer[7] = uint8_t (forecast_hour >> 8);
         buffer[8] = uint8_t (forecast_hour % 256);
         buffer[9] = 10;
         buffer[10] = 0;
         buffer[11] = 0;
         break;

      case PPT3:
         buffer[7] = uint8_t (forecast_hour - 3);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;

      case PPT6:
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;

      case PPTN:
         buffer[7] = 0;
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;

      case TOTAL_CLOUD:
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;

      case HIGH_CLOUD:
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;

      case MIDDLE_CLOUD:
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;

      case LOW_CLOUD:
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;

   }

}

void
Gfs3::set_grib_key (Grib::Key& grib_key,
                    const Nwp_Element nwp_element) const
{

   uint8_t* buffer = grib_key.buffer;
   uint8_t& octet = buffer[12];

   switch (nwp_element)
   {
      case PRESSURE:
         octet = 1;
         break;
      case TEMPERATURE:
         octet = 11;
         break;
      case TEMPERATURE_CELCIUS:
         octet = 0;
         break;
      case MIX_DOWN_TEMPERATURE:
         octet = 0;
         break;
      case ZONAL_WIND:
         octet = 33;
         break;
      case MERIDIONAL_WIND:
         octet = 34;
         break;
      case WIND_SPEED:
         octet = 32;
         break;
      case VERTICAL_VELOCITY:
         octet = 40;
         break;
      case GEOPOTENTIAL_HEIGHT:
         octet = 7;
         break;
      case DEW_POINT:
         octet = 17;
         break;
      case MEAN_SEA_LEVEL_PRESSURE:
         octet = 2;
         break;
      case HIGH_CLOUD:
         octet = 71;
         break;
      case MIDDLE_CLOUD:
         octet = 71;
         break;
      case LOW_CLOUD:
         octet = 71;
         break;
      case TOTAL_CLOUD:
         octet = 71;
         break;
      case OMEGA:
         octet = 39;
         break;
      case PPT3:
         octet = 61;
         break;
      case PPT6:
         octet = 61;
         break;
      case PPTN:
         octet = 61;
         break;
      case RAINFALL_STEP:
         octet = 0;
         break;
      case RELATIVE_HUMIDITY:
         octet = 52;
         break;
      case DEW_POINT_DEPRESSION:
         octet = 0;
         break;
      case POTENTIAL_TEMPERATURE:
         octet = 13;
         break;
      case THETA:
         octet = 13;
         break;
      case THETA_E:
         octet = 14;
         break;
      case MIXING_RATIO:
         octet = 53;
         break;
      case MONTGOMERY:
         octet = 37;
         break;
      case SLI:
         octet = 131;
         break;
      case SHOWALTER:
         octet = 0;
         break;
      case LI_700:
         octet = 0;
         break;
      case LI_THUNDER:
         octet = 0;
         break;
      case K_INDEX:
         octet = 133;
         break;
      case TOTAL_TOTALS:
         octet = 0;
         break;
      case CAPE:
         octet = 157;
         break;
      case PRECIPITABLE_WATER:
         octet = 54;
         break;
      case FOG_FRACTION:
         octet = 0;
         break;
      case THICKNESS:
         octet = 0;
         break;
      case POTENTIAL_VORTICITY:
         octet = 149;
         break;
      case ABSOLUTE_VORTICITY:
         octet = 41;
         break;
      case SHEAR_VORTICITY:
         octet = 0;
         break;
      case CURVATURE_VORTICITY:
         octet = 0;
         break;
   }

}

void
Gfs3::set_grib_key (Grib::Key& grib_key,
                    const Nwp_Element nwp_element,
                    const denise::Level& level) const
{

   uint8_t* buffer = grib_key.buffer;

   switch (nwp_element)
   {
      case TOTAL_CLOUD:
         buffer[13] = 200;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
      case HIGH_CLOUD:
         buffer[13] = 214;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
      case MIDDLE_CLOUD:
         buffer[13] = 224;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
      case LOW_CLOUD:
         buffer[13] = 234;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
      case PRECIPITABLE_WATER:
         buffer[13] = 200;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
   }

   switch (level.type)
   {
      case PRESSURE_LEVEL:
      {
         const uint16_t p = uint16_t (round (level.value * 1e-2));
         buffer[13] = 100;
         buffer[14] = uint8_t (p >> 8);
         buffer[15] = uint8_t (p % 256);
         break;
      }
      case THETA_LEVEL:
      {
         const uint16_t theta = uint16_t (round (level.value));
         buffer[13] = 113;
         buffer[14] = uint8_t (theta >> 8);
         buffer[15] = uint8_t (theta % 256);
         break;
      }
      case SIGMA_LEVEL:
      {
         const uint16_t sigma = uint16_t (round (level.value * 10000));
         buffer[13] = 107;
         buffer[14] = uint8_t (sigma >> 8);
         buffer[15] = uint8_t (sigma % 256);
         break;
      }
      case SCREEN_LEVEL:
      {
         const uint16_t z = uint16_t (2);
         buffer[13] = 105;
         buffer[14] = uint8_t (z >> 8);
         buffer[15] = uint8_t (z % 256);
         break;
      }
      case FIFTY_METRE_LEVEL:
      {
         const uint16_t z = uint16_t (50);
         buffer[13] = 105;
         buffer[14] = uint8_t (z >> 8);
         buffer[15] = uint8_t (z % 256);
         break;
      }
      case TEN_METRE_LEVEL:
      {
         const uint16_t z = uint16_t (10);
         buffer[13] = 105;
         buffer[14] = uint8_t (z >> 8);
         buffer[15] = uint8_t (z % 256);
         break;
      }
      case MEAN_SEA_LEVEL:
      {
         buffer[13] = 102;
         buffer[14] = 0;
         buffer[15] = 0;
         break;
      }
      case SURFACE_LEVEL:
      {
         buffer[13] = 1;
         buffer[14] = 0;
         buffer[15] = 0;
         break;
      }
      case NIL_LEVEL:
      case NOT_A_LEVEL:
      {
         buffer[13] = 0;
         buffer[14] = 0;
         buffer[15] = 0;
         break;
      }
   }

}

void
Gfs3::initialize_3d_data (const Key& key)
{
   typedef Gfs3::Data_3D G3d_3d;
   G3d_3d* g3d_3d_ptr = new G3d_3d (nwp_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, g3d_3d_ptr));
cout << "data_3d_ptr_map.size () = " << &data_3d_ptr_map << ": " << data_3d_ptr_map.size () << "  key = " << key.base_time.get_string () << " " << key.forecast_hour << endl;
}

void          
Gfs3::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Nwp_Element>::const_iterator Iterator;
   const Key& key = data_3d.key;

   for (Iterator iterator = nwp_element_vector.begin ();
        iterator != nwp_element_vector.end (); iterator++)
   {
      const Nwp_Element& nwp_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = get_gvd_3d_ptr (nwp_element, key);
      data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
   }

   data_3d.set_available ();

}

Geodetic_Vector_Data_2D*
Gfs3::get_initialized_vd_2d (const Integer vector_size) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   Gvd_2d* data_ptr = new Gvd_2d (vector_size, size_2d, domain_2d, false);
   return data_ptr;

}

void 
Gfs3::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Level& level,
                              const Nwp_Element nwp_element)
{

   const Integer vi = vector_index;

   switch (nwp_element)
   {

      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& surface = Level::surface_level ();
         fill_grib_data (gvd_2d, vi, nwp_element, key, surface);
         return;
      }

   }

   Nwp::fill_ts_diagnosis_data (gvd_2d, vi, key, level, nwp_element);

}

void
Gfs3::fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key)
{

   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour < 0)
   {
      throw Nwp_Exception ("Forecast Hour < 0");
      return;
   }

   gvd_2d.initialize (vector_index, 0);
   if (forecast_hour == 0) { return; }

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();

   Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
   Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

   for (Integer fh = 0; fh < forecast_hour; fh += 6)
   {
      fill_grib_data (precip_data, 0, PPT6, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (precip_data, 0, PPT3, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   delete precip_data_ptr;

}

void
Gfs3::fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();
   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour == 0)
   {
      gvd_2d.initialize (vector_index, 0);
      return;
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (gvd_2d, vector_index, PPT3, key, surface_level);
   }
   else
   {

      fill_grib_data (gvd_2d, vector_index, PPT6, key, surface_level);

      Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
      Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

      try
      {
         const Dtime& base_time = key.base_time;
         const Integer forecast_hour = key.forecast_hour;
         const Key prev_key (base_time, forecast_hour - 3);
         fill_grib_data (precip_data, 0, PPT3, prev_key, surface_level);
         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real ppt3 = precip_data.get_datum (0, i, j);
               gvd_2d.get_datum (vector_index, i, j) -= ppt3;
            }
         }
      }
      catch (const Nwp_Exception& ne)
      {
      }

      delete precip_data_ptr;

   }

}

void
Gfs3::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Nwp_Element nwp_element)
{

   switch (nwp_element)
   {

      case RAINFALL_CUMULATIVE:
      {
         fill_cumulative_rain_data (gvd_2d, vector_index, key);
         return;
      }

      case RAINFALL_STEP:
      {
         fill_step_rain_data (gvd_2d, vector_index, key);
         return;
      }

   }

   Nwp::fill_rain_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs3::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Key& key,
                       const Nwp_Element nwp_element)
{
   const Level& nil_level = Level::nil_level ();
   fill_grib_data (gvd_2d, vector_index, nwp_element, key, nil_level);
}

void
Gfs3::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Nwp_Element nwp_element)
{

   const Level& screen = Level::screen_level ();
   const Level& surface = Level::surface_level ();

   switch (nwp_element)
   {

      case TEMPERATURE:
      {
         fill_grib_data (gvd_2d, vector_index, TEMPERATURE, key, screen);
         return;
      }

      case DEW_POINT:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, RELATIVE_HUMIDITY);

         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real rh = data_ptr->get_datum (1, i, j);
               const Real t_d = Moisture::get_t_d (t, rh);
               gvd_2d.set_datum (vector_index, i, j, t_d);
            }
         }

         delete data_ptr;
         return;

      }

      case RELATIVE_HUMIDITY:
      {
         fill_grib_data (gvd_2d, vector_index, RELATIVE_HUMIDITY, key, screen);
         return;
      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs3::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key,
                           const Nwp_Element nwp_element)
{

   const Level& ten = Level::ten_metre_level ();

   switch (nwp_element)
   {

      case ZONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, ZONAL_WIND, key, ten);
         return;
      }

      case MERIDIONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, MERIDIONAL_WIND, key, ten);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs3::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                     const Integer vector_index,
                     const Key& key,
                     const Nwp_Element nwp_element)
{

   const Level& msl = Level::mean_sea_level ();

   switch (nwp_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         const Nwp_Element mslp = MEAN_SEA_LEVEL_PRESSURE;
         fill_grib_data (gvd_2d, vector_index, mslp, key, msl);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs3::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Nwp_Element nwp_element)
{

   const Level& surface = Level::surface_level ();

   switch (nwp_element)
   {

      case PRESSURE:
      {
         fill_grib_data (gvd_2d, vector_index, PRESSURE, key, surface);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, nwp_element);

}

/*
void
Gfs3::fill_grib_data (Geodetic_Vector_Data_3D& gvd_3d,
                      const Integer vector_index,
                      const Nwp_Element nwp_element,
                      const Key& key) const
{

   Gfs3::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib& grib = *(iterator->second);

   for (Integer k = 0; k < tuple_p.size (); k++)
   {
      const Real p = tuple_p[k];
      const Level level (PRESSURE_LEVEL, p);
      const Grib::Key& key = get_grib_key (key, nwp_element, level);
      grib.fill_data (gvd_3d, vector_index, k, key);
   }

}
*/

Geodetic_Vector_Data_3D*
Gfs3::get_gvd_3d_ptr (const Nwp_Element nwp_element,
                      const Key& key) const
{

   Geodetic_Vector_Data_3D* gvd_3d_ptr =
      new Geodetic_Vector_Data_3D (1, tuple_p, size_2d, domain_2d);
   Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   for (Integer k = 0; k < tuple_p.size (); k++)
   {

      const Real p = tuple_p[k];
      const denise::Level level (PRESSURE_LEVEL, p);
      const Grib::Key& grib_key = get_grib_key (key, nwp_element, level);

      Gfs3::const_iterator iterator = find (key);
      if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
      const Grib& grib = *(iterator->second);

      Geodetic_Vector_Data_2D* grib_data_ptr =
         get_grib_data_ptr (grib, grib_key);

      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         const Real latitude = grib_data_ptr->get_latitude (i);
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real longitude = grib_data_ptr->get_longitude (j);
            const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
            gvd_3d.set_datum (0, k, i, j, datum);
         }
      }

      delete grib_data_ptr;

   }

   return gvd_3d_ptr;

}

Geodetic_Vector_Data_2D*
Gfs3::get_grib_data_ptr (const Grib& grib,
                         const Grib::Key& grib_key) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = new Gvd_2d (1, size_2d, domain_2d, true);
   grib.fill_data (*data_ptr, 0, grib_key);

   const Integer last_j = size_2d.j - 1;

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real datum = data_ptr->get_datum (0, i, 0);
      data_ptr->set_datum (0, i, last_j, datum);
   }

   return data_ptr;

}

void
Gfs3::fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Nwp_Element nwp_element,
                      const Key& key,
                      const Level& level) const
{

   Gfs3::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib& grib = *(iterator->second);

   const Grib::Key& grib_key = get_grib_key (key, nwp_element, level);
   Geodetic_Vector_Data_2D* grib_data_ptr = get_grib_data_ptr (grib, grib_key);

   Geodetic_Vector_Data_2D grib_data (1, size_2d, domain_2d, true);

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_latitude (i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_longitude (j);
         const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete grib_data_ptr;

}

Gfs3::Gfs3 (const string& data_path)
   : Nwp ("Gfs3", data_path),
     data_path (data_path),
     size_2d (181, 361),
     domain_2d (Domain_1D (-90, 90), Domain_1D (0, 360))
{

   this->status = "Unloaded";

   nwp_element_vector.push_back (denise::TEMPERATURE);
   nwp_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   nwp_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::ZONAL_WIND);
   nwp_element_vector.push_back (denise::MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::OMEGA);

}

Gfs3::~Gfs3 ()
{
   clean_up ();
}

void
Gfs3::survey ()
{

   const string re_ym ("[0-9].....");
   const string re_ymd ("[0-9].......");
   const string fmt ("gfs_3_[0-9]......._[0-9]..._[0-9]..\\.grb");

   typedef map<Grib::Key, Grib::Header*> Header_Ptr_Map;

   const vector<string>& dir_ym = get_dir_listing (path, re_ym);
   for (vector<string>::const_iterator i = dir_ym.begin ();
        i != dir_ym.end (); i++)
   {

      const string& ym = *(i);
      const string& ym_path = path + "/" + ym;

      const vector<string>& dir_ymd = get_dir_listing (ym_path, re_ymd);
      for (vector<string>::const_iterator j = dir_ymd.begin ();
           j != dir_ymd.end (); j++)
      {

         const string& ymd = *(j);
         const string& ymd_path = ym_path + "/" + ymd;

         const vector<string>& dir_listing = get_dir_listing (ymd_path, fmt);
         for (vector<string>::const_iterator iterator = dir_listing.begin ();
              iterator != dir_listing.end (); iterator++)
         {

            // fn = filename
            const string& fn = *(iterator);
            const string& file_path = ymd_path + "/" + fn;
            const string& bt_str = fn.substr (6, 8) + fn.substr (15, 2);
            const string& fh_str = fn.substr (20, 3);

            const Dtime base_time (bt_str, string ("%Y%m%d%H"));
            const Integer forecast_hour = atoi (fh_str.c_str ());
            const Key key (base_time, forecast_hour);

            key_multimap.add (key);
            Grib* grib_ptr = new Grib (file_path);
            insert (make_pair (key, grib_ptr));

         }

      }

   }

   if (size () > 0)
   {

      set<uint16_t> set_p;
      const Grib& grib = *(begin ()->second);
      const Header_Ptr_Map& header_ptr_map = grib.get_header_ptr_map ();

      for (Header_Ptr_Map::const_iterator iterator = header_ptr_map.begin ();
           iterator != header_ptr_map.end (); iterator++)
      {

         const Grib::Header& header = *(iterator->second);
         const Grib::Pds& pds = header.get_pds ();
         const Grib::Pds::Level& level = pds.get_level ();

         if (level.get_uint (0, 1) == 100)
         {
            const uint16_t p = level.get_uint (1, 2);
            if (p < 200) { continue; }
            set_p.insert (p);
         }

      }

      for (set<uint16_t>::const_iterator iterator = set_p.begin ();
           iterator != set_p.end (); iterator++)
      {
         const Real p = Real (*(iterator)) * 100;
         tuple_p.push_back (p);
      }

   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

   for (Key_Multimap::const_iterator iterator = key_multimap.begin ();
        iterator != key_multimap.end (); iterator++)
   {
      const Key key (iterator->first, iterator->second);
      initialize_3d_data (key);
   }

}

void
Gfs3::clear_3d_data ()
{
}

void
Gfs3::clean_up ()
{
   for (Gfs3::iterator iterator = begin (); iterator != end (); iterator++)
   {
      Grib* grib_ptr = iterator->second;
      delete grib_ptr;
   }
   clear ();
}

vector<Dtime>
Gfs3::get_valid_time_vector () const
{

   vector<Dtime> valid_time_vector;

   for (Gfs3::const_iterator iterator = begin ();
       iterator != end (); iterator++)
   {
      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime dtime (bt.t + fh);
      valid_time_vector.push_back (dtime);
   }

   return valid_time_vector;

/*
   vector<Dtime> valid_time_vector;
   Dtime start_time ("2012011300");
   Dtime end_time ("2012011400");
   for (Real t = start_time.t; t < end_time.t; t += 3)
   {
      Dtime dtime (t);
      valid_time_vector.push_back (dtime);
   }
   return valid_time_vector;
*/
}

Nwp::Key
Gfs3::get_key (const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs3::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());

   const Dtime& base_time = base_time_vector.back ();
   const Integer forecast_hour = Integer (round (dtime.t - base_time.t));

   return Key (base_time, forecast_hour);

}

void
Gfs3::acquire_base_time_forecast_hour (Dtime& base_time,
                                       Integer& forecast_hour,
                                       const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs3::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());
   base_time.t = base_time_vector.back ().t;
   forecast_hour = Integer (round (dtime.t - base_time.t));

}

Gfs4::Data_3D::Data_3D (const vector<Nwp_Element>& nwp_element_vector,
                        const Key& key)
   : Nwp::Data_3D (nwp_element_vector, key)
{
}

Real
Gfs4::Data_3D::evaluate (const Nwp_Element element,
                         const Real p,
                         const Real latitude,
                         const Real longitude,
                         const Evaluate_Op evaluate_op) const
{

   typedef Nwp::Data_3D Nd_3d;
   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& RH = denise::RELATIVE_HUMIDITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real rh = Nwp::Data_3D::evaluate (RH, p, latitude, longitude);
         return Moisture::get_t_d (t, rh);
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& O = denise::OMEGA;

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Evaluate_Op& eo = evaluate_op;
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real ts = Nd_3d::evaluate (T, p, latitude, longitude, eo);
            const Real os = Nd_3d::evaluate (O, p, latitude, longitude, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nd_3d::evaluate (element, p, latitude, longitude, eo);

}

Grib2::Key
Gfs4::get_grib_key (const Key& key,
                    const Nwp_Element nwp_element,
                    const Level& level) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   Grib2::Key grib_key;
   set_grib_key (grib_key, base_time);
   set_grib_key (grib_key, nwp_element, forecast_hour);
   set_grib_key (grib_key, nwp_element);
   set_grib_key (grib_key, nwp_element, level);
   return grib_key;

}

void
Gfs4::set_grib_key (Grib2::Key& grib_key,
                    const Dtime& base_time) const
{

   uint8_t* buffer = grib_key.buffer;

   const Integer yyyy = base_time.get_year ();
   const Integer mm = base_time.get_month ();
   const Integer dd = base_time.get_day ();
   const Integer HH = base_time.get_hour ();
   const Integer MM = base_time.get_minute ();
   const Integer SS = base_time.get_second ();

   buffer[0] = uint8_t (yyyy / 256);
   buffer[1] = uint8_t (yyyy % 256);
   buffer[2] = uint8_t (mm);
   buffer[3] = uint8_t (dd);
   buffer[4] = uint8_t (HH);
   buffer[5] = uint8_t (MM);
   buffer[6] = uint8_t (SS);

}

void
Gfs4::set_grib_key (Grib2::Key& grib_key,
                    const Nwp_Element nwp_element,
                    const Integer forecast_hour) const
{

   uint8_t* buffer = grib_key.buffer;
   buffer[7] = 1;

   switch (nwp_element)
   {

      default:
      {
         uint32_t value = forecast_hour;
#ifndef WORDS_BIGENDIAN
         swap_endian (&value, sizeof (uint32_t));
#endif
         memcpy (buffer + 8, &value, sizeof (uint32_t));
         break;
      }

      case PPT3:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 3);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case PPT6:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case PPTN:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = 0;
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case TOTAL_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case HIGH_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case MIDDLE_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case LOW_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

   }

}

void
Gfs4::set_grib_key (Grib2::Key& grib_key,
                    const Nwp_Element nwp_element) const
{

   uint8_t* buffer = grib_key.buffer;

   switch (nwp_element)
   {
      case PRESSURE:
         buffer[12] = 3;
         buffer[13] = 0;
         break;
      case TEMPERATURE:
         buffer[12] = 0;
         buffer[13] = 0;
         break;
      case TEMPERATURE_CELCIUS:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case MIX_DOWN_TEMPERATURE:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case ZONAL_WIND:
         buffer[12] = 2;
         buffer[13] = 2;
         break;
      case MERIDIONAL_WIND:
         buffer[12] = 2;
         buffer[13] = 3;
         break;
      case WIND_SPEED:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case VERTICAL_VELOCITY:
         buffer[12] = 2;
         buffer[13] = 9;
         break;
      case GEOPOTENTIAL_HEIGHT:
         buffer[12] = 3;
         buffer[13] = 5;
         break;
      case DEW_POINT:
         buffer[12] = 0;
         buffer[13] = 6;
         break;
      case MEAN_SEA_LEVEL_PRESSURE:
         buffer[12] = 3;
         buffer[13] = 1;
         break;
      case HIGH_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case MIDDLE_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case LOW_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case TOTAL_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case OMEGA:
         buffer[12] = 2;
         buffer[13] = 8;
         break;
      case PPT3:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case PPT6:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case PPTN:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case RAINFALL_STEP:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case RELATIVE_HUMIDITY:
         buffer[12] = 1;
         buffer[13] = 1;
         break;
      case DEW_POINT_DEPRESSION:
         buffer[12] = 0;
         buffer[13] = 7;
         break;
      case POTENTIAL_TEMPERATURE:
         buffer[12] = 0;
         buffer[13] = 2;
         break;
      case THETA:
         buffer[12] = 0;
         buffer[13] = 2;
         break;
      case THETA_E:
         buffer[12] = 0;
         buffer[13] = 3;
         break;
      case MIXING_RATIO:
         buffer[12] = 1;
         buffer[13] = 2;
         break;
      case MONTGOMERY:
         buffer[12] = 2;
         buffer[13] = 6;
         break;
      case SLI:
         buffer[12] = 7;
         buffer[13] = 10;
         break;
      case SHOWALTER:
         buffer[12] = 7;
         buffer[13] = 13;
         break;
      case LI_700:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case LI_THUNDER:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case K_INDEX:
         buffer[12] = 7;
         buffer[13] = 2;
         break;
      case TOTAL_TOTALS:
         buffer[12] = 7;
         buffer[13] = 4;
         break;
      case CAPE:
         buffer[12] = 7;
         buffer[13] = 6;
         break;
      case PRECIPITABLE_WATER:
         buffer[12] = 1;
         buffer[13] = 3;
         break;
      case FOG_FRACTION:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case THICKNESS:
         buffer[12] = 3;
         buffer[13] = 12;
         break;
      case POTENTIAL_VORTICITY:
         buffer[12] = 2;
         buffer[13] = 14;
         break;
      case ABSOLUTE_VORTICITY:
         buffer[12] = 2;
         buffer[13] = 10;
         break;
      case SHEAR_VORTICITY:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case CURVATURE_VORTICITY:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
   }

}

void
Gfs4::set_grib_key (Grib2::Key& grib_key,
                    const Nwp_Element nwp_element,
                    const denise::Level& level) const
{

   uint8_t* buffer = grib_key.buffer;

   switch (nwp_element)
   {
      case TOTAL_CLOUD:
         buffer[14] = 200;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case HIGH_CLOUD:
         buffer[14] = 214;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case MIDDLE_CLOUD:
         buffer[14] = 224;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case LOW_CLOUD:
         buffer[14] = 234;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case PRECIPITABLE_WATER:
         buffer[13] = 200;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
   }

   switch (level.type)
   {
      case PRESSURE_LEVEL:
      {
         uint32_t p = uint32_t (round (level.value));
         buffer[14] = 100;
         buffer[15] = 0;
#ifndef WORDS_BIGENDIAN
         swap_endian (&p, sizeof (uint32_t));
#endif
         memcpy (buffer + 16, &p, sizeof (uint32_t));
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case THETA_LEVEL:
      {
         uint32_t theta = uint32_t (round (level.value));
         buffer[14] = 107;
         buffer[15] = 0;
#ifndef WORDS_BIGENDIAN
         swap_endian (&theta, sizeof (uint32_t));
#endif
         memcpy (buffer + 16, &theta, sizeof (uint32_t));
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case SIGMA_LEVEL:
      {
         uint32_t sigma = uint32_t (round (level.value * 10000));
         buffer[14] = 113;
         buffer[15] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case SCREEN_LEVEL:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 2;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case FIFTY_METRE_LEVEL:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 50;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case TEN_METRE_LEVEL:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 10;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case MEAN_SEA_LEVEL:
      {
         buffer[14] = 101;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case SURFACE_LEVEL:
      {
         buffer[14] = 1;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case NIL_LEVEL:
      case NOT_A_LEVEL:
      {
         buffer[14] = 1;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
   }

}

void
Gfs4::initialize_3d_data (const Key& key)
{
   typedef Gfs4::Data_3D G3d_3d;
   G3d_3d* g3d_3d_ptr = new G3d_3d (nwp_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, g3d_3d_ptr));
}

void          
Gfs4::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Nwp_Element>::const_iterator Iterator;
   const Key& key = data_3d.key;

   for (Iterator iterator = nwp_element_vector.begin ();
        iterator != nwp_element_vector.end (); iterator++)
   {
      const Nwp_Element& nwp_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = get_gvd_3d_ptr (nwp_element, key);
      data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
   }

   data_3d.set_available ();

}

Geodetic_Vector_Data_2D*
Gfs4::get_initialized_vd_2d (const Integer vector_size) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   Gvd_2d* data_ptr = new Gvd_2d (vector_size, size_2d, domain_2d, false);
   return data_ptr;

}

void 
Gfs4::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Level& level,
                              const Nwp_Element nwp_element)
{

   const Integer vi = vector_index;

   switch (nwp_element)
   {

      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& surface = Level::surface_level ();
         fill_grib_data (gvd_2d, vi, nwp_element, key, surface);
         return;
      }

   }

   Nwp::fill_ts_diagnosis_data (gvd_2d, vi, key, level, nwp_element);

}

void
Gfs4::fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key)
{

   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour < 0)
   {
      throw Nwp_Exception ("Forecast Hour < 0");
      return;
   }

   gvd_2d.initialize (vector_index, 0);
   if (forecast_hour == 0) { return; }

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();

   Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
   Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

   for (Integer fh = 0; fh < forecast_hour; fh += 6)
   {
      fill_grib_data (precip_data, 0, PPT6, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (precip_data, 0, PPT3, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   delete precip_data_ptr;

}

void
Gfs4::fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();
   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour == 0)
   {
      gvd_2d.initialize (vector_index, 0);
      return;
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (gvd_2d, vector_index, PPT3, key, surface_level);
   }
   else
   {

      fill_grib_data (gvd_2d, vector_index, PPT6, key, surface_level);

      Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
      Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

      try
      {

         const Dtime& base_time = key.base_time;
         const Integer forecast_hour = key.forecast_hour;
         const Key prev_key (base_time, forecast_hour - 3);

         fill_grib_data (precip_data, 0, PPT3, prev_key, surface_level);
         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real ppt3 = precip_data.get_datum (0, i, j);
               gvd_2d.get_datum (vector_index, i, j) -= ppt3;
            }
         }
      }
      catch (const Nwp_Exception& ne)
      {
      }

      delete precip_data_ptr;

   }

}

void
Gfs4::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Nwp_Element nwp_element)
{

   switch (nwp_element)
   {

      case RAINFALL_CUMULATIVE:
      {
         fill_cumulative_rain_data (gvd_2d, vector_index, key);
         return;
      }

      case RAINFALL_STEP:
      {
         fill_step_rain_data (gvd_2d, vector_index, key);
         return;
      }

   }

   Nwp::fill_rain_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs4::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Key& key,
                       const Nwp_Element nwp_element)
{
   const Level& nil_level = Level::nil_level ();
   fill_grib_data (gvd_2d, vector_index, nwp_element, key, nil_level);
}

void
Gfs4::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Nwp_Element nwp_element)
{

   const Level& screen = Level::screen_level ();
   const Level& surface = Level::surface_level ();

   switch (nwp_element)
   {

      case TEMPERATURE:
      {
         fill_grib_data (gvd_2d, vector_index, TEMPERATURE, key, screen);
         return;
      }

      case DEW_POINT:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, RELATIVE_HUMIDITY);

         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real rh = data_ptr->get_datum (1, i, j);
               const Real t_d = Moisture::get_t_d (t, rh);
               gvd_2d.set_datum (vector_index, i, j, t_d);
            }
         }

         delete data_ptr;
         return;

      }

      case RELATIVE_HUMIDITY:
      {
         fill_grib_data (gvd_2d, vector_index, RELATIVE_HUMIDITY, key, screen);
         return;
      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs4::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key,
                           const Nwp_Element nwp_element)
{

   const Level& ten = Level::ten_metre_level ();

   switch (nwp_element)
   {

      case ZONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, ZONAL_WIND, key, ten);
         return;
      }

      case MERIDIONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, MERIDIONAL_WIND, key, ten);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs4::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                     const Integer vector_index,
                     const Key& key,
                     const Nwp_Element nwp_element)
{

   const Level& msl = Level::mean_sea_level ();

   switch (nwp_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         const Nwp_Element mslp = MEAN_SEA_LEVEL_PRESSURE;
         fill_grib_data (gvd_2d, vector_index, mslp, key, msl);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs4::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Nwp_Element nwp_element)
{

   const Level& surface = Level::surface_level ();

   switch (nwp_element)
   {

      case PRESSURE:
      {
         fill_grib_data (gvd_2d, vector_index, PRESSURE, key, surface);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, nwp_element);

}

/*
void
Gfs4::fill_grib_data (Geodetic_Vector_Data_3D& gvd_3d,
                      const Integer vector_index,
                      const Nwp_Element nwp_element,
                      const Key& key) const
{

   const Key key (key);

   Gfs4::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib2& grib = *(iterator->second);

   for (Integer k = 0; k < tuple_p.size (); k++)
   {
      const Real p = tuple_p[k];
      const Level level (PRESSURE_LEVEL, p);
      const Grib2::Key& key = get_grib_key (key, nwp_element, level);
      grib.fill_data (gvd_3d, vector_index, k, key);
   }

}
*/

Geodetic_Vector_Data_3D*
Gfs4::get_gvd_3d_ptr (const Nwp_Element nwp_element,
                      const Key& key) const
{

   Geodetic_Vector_Data_3D* gvd_3d_ptr =
      new Geodetic_Vector_Data_3D (1, tuple_p, size_2d, domain_2d);
   Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   for (Integer k = 0; k < tuple_p.size (); k++)
   {

      const Real p = tuple_p[k];
      const denise::Level level (PRESSURE_LEVEL, p);
      const Grib2::Key& grib_key = get_grib_key (key, nwp_element, level);

      Gfs4::const_iterator iterator = find (key);
      if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
      const Grib2& grib = *(iterator->second);

      Geodetic_Vector_Data_2D* grib_data_ptr =
         get_grib_data_ptr (grib, grib_key);

      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         const Real latitude = grib_data_ptr->get_latitude (i);
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real longitude = grib_data_ptr->get_longitude (j);
            const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
            gvd_3d.set_datum (0, k, i, j, datum);
         }
      }

      delete grib_data_ptr;

   }

   return gvd_3d_ptr;

}

Geodetic_Vector_Data_2D*
Gfs4::get_grib_data_ptr (const Grib2& grib,
                         const Grib2::Key& grib_key) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = new Gvd_2d (1, grib_size_2d, grib_domain_2d, true);
   grib.fill_data (*data_ptr, 0, grib_key);

   const Integer last_j = grib_size_2d.j - 1;

   #pragma omp parallel for
   for (Integer i = 0; i < grib_size_2d.i; i++)
   {
      const Real datum = data_ptr->get_datum (0, i, 0);
      data_ptr->set_datum (0, i, last_j, datum);
   }

   return data_ptr;

}

void
Gfs4::fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Nwp_Element nwp_element,
                      const Key& key,
                      const Level& level) const
{

   Gfs4::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib2& grib = *(iterator->second);

   const Grib2::Key& grib_key = get_grib_key (key, nwp_element, level);
   Geodetic_Vector_Data_2D* grib_data_ptr = get_grib_data_ptr (grib, grib_key);

   Geodetic_Vector_Data_2D grib_data (1, grib_size_2d, grib_domain_2d, true);

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_latitude (i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_longitude (j);
         const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete grib_data_ptr;

}

Gfs4::Gfs4 (const string& data_path)
   : Nwp ("Gfs4", data_path),
     data_path (data_path),
     grib_size_2d (361, 720),
     grib_domain_2d (Domain_1D (-90, 90), Domain_1D (0, 360)),
     size_2d (51, 71),
     domain_2d (Domain_1D (-60, 0), Domain_1D (80, 160))
{

   this->status = "Unloaded";

   nwp_element_vector.push_back (denise::TEMPERATURE);
   nwp_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   nwp_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::ZONAL_WIND);
   nwp_element_vector.push_back (denise::MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::OMEGA);

}

Gfs4::~Gfs4 ()
{
   clean_up ();
}

void
Gfs4::survey ()
{

   const string re_ym ("[0-9].....");
   const string re_ymd ("[0-9].......");
   const string fmt ("gfs_4_[0-9]......._[0-9]..._[0-9]..\\.grb");

   typedef map<Grib2::Key, Grib2::Header*> Header_Ptr_Map;

   const vector<string>& dir_ym = get_dir_listing (path, re_ym);
   for (vector<string>::const_iterator i = dir_ym.begin ();
        i != dir_ym.end (); i++)
   {

      const string& ym = *(i);
      const string& ym_path = path + "/" + ym;

      const vector<string>& dir_ymd = get_dir_listing (ym_path, re_ymd);
      for (vector<string>::const_iterator j = dir_ymd.begin ();
           j != dir_ymd.end (); j++)
      {

         const string& ymd = *(j);
         const string& ymd_path = ym_path + "/" + ymd;

         const vector<string>& dir_listing = get_dir_listing (ymd_path, fmt);
         for (vector<string>::const_iterator iterator = dir_listing.begin ();
              iterator != dir_listing.end (); iterator++)
         {

            // fn = filename
            const string& fn = *(iterator);
            const string& file_path = ymd_path + "/" + fn;
            const string& bt_str = fn.substr (6, 8) + fn.substr (15, 2);
            const string& fh_str = fn.substr (20, 3);

            const Dtime base_time (bt_str, string ("%Y%m%d%H"));
            const Integer forecast_hour = atoi (fh_str.c_str ());

            const Key key (base_time, forecast_hour);
            key_multimap.add (key);
            Grib2* grib_ptr = new Grib2 (file_path);
            insert (make_pair (key, grib_ptr));
            //valid_time_set.insert (valid_time);
            //key_set.insert (Key (base_time, fh));

         }

      }

   }

   if (size () > 0)
   {

      set<uint32_t> set_p;
      const Grib2& grib = *(begin ()->second);
      const Header_Ptr_Map& header_ptr_map = grib.get_header_ptr_map ();

      for (Header_Ptr_Map::const_iterator iterator = header_ptr_map.begin ();
           iterator != header_ptr_map.end (); iterator++)
      {

         const Grib2::Header& header = *(iterator->second);
         const Grib2::Block_4& block_4 = header.block_4;
         const Grib2::Block_4::Level& level = block_4.get_level ();

         if (level.get_first_level_type () == 100)
         {
            const uint32_t p = level.get_first_level ();
            if (p < 200) { continue; }
            set_p.insert (p);
         }

      }

      for (set<uint32_t>::const_iterator iterator = set_p.begin ();
           iterator != set_p.end (); iterator++)
      {
         const Real p = Real (*(iterator));
         tuple_p.push_back (p);
      }

   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

}

void
Gfs4::clean_up ()
{
   for (Gfs4::iterator iterator = begin (); iterator != end (); iterator++)
   {
      Grib2* grib_ptr = iterator->second;
      delete grib_ptr;
   }
}

void
Gfs4::set_domain_2d (const Domain_2D& domain_2d)
{

   this->domain_2d.domain_x.start = ceil (domain_2d.domain_x.start);
   this->domain_2d.domain_x.end = floor (domain_2d.domain_x.end);
   this->domain_2d.domain_y.start = ceil (domain_2d.domain_y.start);
   this->domain_2d.domain_y.end = floor (domain_2d.domain_y.end);

   const Real& start_latitude = this->domain_2d.domain_x.start;
   const Real& end_latitude = this->domain_2d.domain_x.end;
   const Real& start_longitude = this->domain_2d.domain_y.start;
   const Real& end_longitude = this->domain_2d.domain_y.end;
 
   size_2d.i = Integer (round (end_latitude - start_latitude));
   size_2d.j = Integer (round (end_longitude - start_longitude));

}

vector<Dtime>
Gfs4::get_valid_time_vector () const
{

   vector<Dtime> valid_time_vector;

   for (Gfs4::const_iterator iterator = begin ();
       iterator != end (); iterator++)
   {
      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime dtime (bt.t + fh);
      valid_time_vector.push_back (dtime);
   }

   return valid_time_vector;

/*
   vector<Dtime> valid_time_vector;
   Dtime start_time ("2012011300");
   Dtime end_time ("2012011400");
   for (Real t = start_time.t; t < end_time.t; t += 3)
   {
      Dtime dtime (t);
      valid_time_vector.push_back (dtime);
   }
   return valid_time_vector;
*/
}

Nwp::Key
Gfs4::get_key (const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs4::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());

   const Dtime& base_time = base_time_vector.back ();
   const Integer forecast_hour = Integer (round (dtime.t - base_time.t));

   return Key (base_time, forecast_hour);

}

void
Gfs4::acquire_base_time_forecast_hour (Dtime& base_time,
                                       Integer& forecast_hour,
                                       const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs4::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());
   base_time.t = base_time_vector.back ().t;
   forecast_hour = Integer (round (dtime.t - base_time.t));

}

Gfs::Data_3D::Data_3D (const vector<Nwp_Element>& nwp_element_vector,
                       const Key& key)
   : Nwp::Data_3D (nwp_element_vector, key)
{
}

Real
Gfs::Data_3D::evaluate (const Nwp_Element element,
                        const Real p,
                        const Real latitude,
                        const Real longitude,
                        const Evaluate_Op evaluate_op) const
{

   typedef Nwp::Data_3D Nd_3d;
   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& RH = denise::RELATIVE_HUMIDITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real rh = Nwp::Data_3D::evaluate (RH, p, latitude, longitude);
         return Moisture::get_t_d (t, rh);
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& O = denise::OMEGA;

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Evaluate_Op& eo = evaluate_op;
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real ts = Nd_3d::evaluate (T, p, latitude, longitude, eo);
            const Real os = Nd_3d::evaluate (O, p, latitude, longitude, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nd_3d::evaluate (element, p, latitude, longitude, eo);

}

Grib2::Key
Gfs::get_grib_key (const Key& key,
                   const Nwp_Element nwp_element,
                   const Level& level) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   Grib2::Key grib_key;
   set_grib_key (grib_key, base_time);
   set_grib_key (grib_key, nwp_element, forecast_hour);
   set_grib_key (grib_key, nwp_element);
   set_grib_key (grib_key, nwp_element, level);
   return grib_key;

}

void
Gfs::set_grib_key (Grib2::Key& grib_key,
                   const Dtime& base_time) const
{

   uint8_t* buffer = grib_key.buffer;

   const Integer yyyy = base_time.get_year ();
   const Integer mm = base_time.get_month ();
   const Integer dd = base_time.get_day ();
   const Integer HH = base_time.get_hour ();
   const Integer MM = base_time.get_minute ();
   const Integer SS = base_time.get_second ();

   buffer[0] = uint8_t (yyyy / 256);
   buffer[1] = uint8_t (yyyy % 256);
   buffer[2] = uint8_t (mm);
   buffer[3] = uint8_t (dd);
   buffer[4] = uint8_t (HH);
   buffer[5] = uint8_t (MM);
   buffer[6] = uint8_t (SS);

}

void
Gfs::set_grib_key (Grib2::Key& grib_key,
                   const Nwp_Element nwp_element,
                   const Integer forecast_hour) const
{

   uint8_t* buffer = grib_key.buffer;
   buffer[7] = 1;

   switch (nwp_element)
   {

      default:
      {
         uint32_t value = forecast_hour;
#ifndef WORDS_BIGENDIAN
         swap_endian (&value, sizeof (uint32_t));
#endif
         memcpy (buffer + 8, &value, sizeof (uint32_t));
         break;
      }

      case PPT3:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 3);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case PPT6:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case PPTN:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = 0;
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case TOTAL_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case HIGH_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case MIDDLE_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case LOW_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

   }

}

void
Gfs::set_grib_key (Grib2::Key& grib_key,
                   const Nwp_Element nwp_element) const
{

   uint8_t* buffer = grib_key.buffer;

   switch (nwp_element)
   {
      case PRESSURE:
         buffer[12] = 3;
         buffer[13] = 0;
         break;
      case TEMPERATURE:
         buffer[12] = 0;
         buffer[13] = 0;
         break;
      case TEMPERATURE_CELCIUS:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case MIX_DOWN_TEMPERATURE:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case ZONAL_WIND:
         buffer[12] = 2;
         buffer[13] = 2;
         break;
      case MERIDIONAL_WIND:
         buffer[12] = 2;
         buffer[13] = 3;
         break;
      case WIND_SPEED:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case VERTICAL_VELOCITY:
         buffer[12] = 2;
         buffer[13] = 9;
         break;
      case GEOPOTENTIAL_HEIGHT:
         buffer[12] = 3;
         buffer[13] = 5;
         break;
      case DEW_POINT:
         buffer[12] = 0;
         buffer[13] = 6;
         break;
      case MEAN_SEA_LEVEL_PRESSURE:
         buffer[12] = 3;
         buffer[13] = 1;
         break;
      case HIGH_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case MIDDLE_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case LOW_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case TOTAL_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case OMEGA:
         buffer[12] = 2;
         buffer[13] = 8;
         break;
      case PPT3:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case PPT6:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case PPTN:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case RAINFALL_STEP:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case RELATIVE_HUMIDITY:
         buffer[12] = 1;
         buffer[13] = 1;
         break;
      case DEW_POINT_DEPRESSION:
         buffer[12] = 0;
         buffer[13] = 7;
         break;
      case POTENTIAL_TEMPERATURE:
         buffer[12] = 0;
         buffer[13] = 2;
         break;
      case THETA:
         buffer[12] = 0;
         buffer[13] = 2;
         break;
      case THETA_E:
         buffer[12] = 0;
         buffer[13] = 3;
         break;
      case MIXING_RATIO:
         buffer[12] = 1;
         buffer[13] = 2;
         break;
      case MONTGOMERY:
         buffer[12] = 2;
         buffer[13] = 6;
         break;
      case SLI:
         buffer[12] = 7;
         buffer[13] = 10;
         break;
      case SHOWALTER:
         buffer[12] = 7;
         buffer[13] = 13;
         break;
      case LI_700:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case LI_THUNDER:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case K_INDEX:
         buffer[12] = 7;
         buffer[13] = 2;
         break;
      case TOTAL_TOTALS:
         buffer[12] = 7;
         buffer[13] = 4;
         break;
      case CAPE:
         buffer[12] = 7;
         buffer[13] = 6;
         break;
      case PRECIPITABLE_WATER:
         buffer[12] = 1;
         buffer[13] = 3;
         break;
      case FOG_FRACTION:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case THICKNESS:
         buffer[12] = 3;
         buffer[13] = 12;
         break;
      case POTENTIAL_VORTICITY:
         buffer[12] = 2;
         buffer[13] = 14;
         break;
      case ABSOLUTE_VORTICITY:
         buffer[12] = 2;
         buffer[13] = 10;
         break;
      case SHEAR_VORTICITY:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case CURVATURE_VORTICITY:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
   }

}

void
Gfs::set_grib_key (Grib2::Key& grib_key,
                   const Nwp_Element nwp_element,
                   const denise::Level& level) const
{

   uint8_t* buffer = grib_key.buffer;

   switch (nwp_element)
   {
      case TOTAL_CLOUD:
         buffer[14] = 200;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case HIGH_CLOUD:
         buffer[14] = 214;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case MIDDLE_CLOUD:
         buffer[14] = 224;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case LOW_CLOUD:
         buffer[14] = 234;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case PRECIPITABLE_WATER:
         buffer[13] = 200;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
   }

   switch (level.type)
   {
      case PRESSURE_LEVEL:
      {
         uint32_t p = uint32_t (round (level.value));
         buffer[14] = 100;
         buffer[15] = 0;
#ifndef WORDS_BIGENDIAN
         swap_endian (&p, sizeof (uint32_t));
#endif
         memcpy (buffer + 16, &p, sizeof (uint32_t));
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case THETA_LEVEL:
      {
         uint32_t theta = uint32_t (round (level.value));
         buffer[14] = 107;
         buffer[15] = 0;
#ifndef WORDS_BIGENDIAN
         swap_endian (&theta, sizeof (uint32_t));
#endif
         memcpy (buffer + 16, &theta, sizeof (uint32_t));
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case SIGMA_LEVEL:
      {
         uint32_t sigma = uint32_t (round (level.value * 10000));
         buffer[14] = 113;
         buffer[15] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case SCREEN_LEVEL:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 2;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case FIFTY_METRE_LEVEL:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 50;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case TEN_METRE_LEVEL:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 10;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case MEAN_SEA_LEVEL:
      {
         buffer[14] = 101;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case SURFACE_LEVEL:
      {
         buffer[14] = 1;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case NIL_LEVEL:
      case NOT_A_LEVEL:
      {
         buffer[14] = 1;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
   }

}

void
Gfs::initialize_3d_data (const Key& key)
{
   typedef Gfs::Data_3D G3d_3d;
   G3d_3d* g3d_3d_ptr = new G3d_3d (nwp_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, g3d_3d_ptr));
}

void          
Gfs::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Nwp_Element>::const_iterator Iterator;
   const Key& key = data_3d.key;

   for (Iterator iterator = nwp_element_vector.begin ();
        iterator != nwp_element_vector.end (); iterator++)
   {
      const Nwp_Element& nwp_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = get_gvd_3d_ptr (nwp_element, key);
      data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
   }

   data_3d.set_available ();

}

Geodetic_Vector_Data_2D*
Gfs::get_initialized_vd_2d (const Integer vector_size) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   Gvd_2d* data_ptr = new Gvd_2d (vector_size, size_2d, domain_2d, false);
   return data_ptr;

}

void 
Gfs::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Level& level,
                              const Nwp_Element nwp_element)
{

   const Integer vi = vector_index;

   switch (nwp_element)
   {

      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& surface = Level::surface_level ();
         fill_grib_data (gvd_2d, vi, nwp_element, key, surface);
         return;
      }

   }

   Nwp::fill_ts_diagnosis_data (gvd_2d, vi, key, level, nwp_element);

}

void
Gfs::fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key)
{

   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour < 0)
   {
      throw Nwp_Exception ("Forecast Hour < 0");
      return;
   }

   gvd_2d.initialize (vector_index, 0);
   if (forecast_hour == 0) { return; }

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();

   Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
   Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

   for (Integer fh = 0; fh < forecast_hour; fh += 6)
   {
      fill_grib_data (precip_data, 0, PPT6, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (precip_data, 0, PPT3, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   delete precip_data_ptr;

}

void
Gfs::fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();
   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour == 0)
   {
      gvd_2d.initialize (vector_index, 0);
      return;
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (gvd_2d, vector_index, PPT3, key, surface_level);
   }
   else
   {

      fill_grib_data (gvd_2d, vector_index, PPT6, key, surface_level);

      Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
      Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

      try
      {

         const Dtime& base_time = key.base_time;
         const Integer forecast_hour = key.forecast_hour;
         const Key prev_key (base_time, forecast_hour - 3);

         fill_grib_data (precip_data, 0, PPT3, prev_key, surface_level);
         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real ppt3 = precip_data.get_datum (0, i, j);
               gvd_2d.get_datum (vector_index, i, j) -= ppt3;
            }
         }
      }
      catch (const Nwp_Exception& ne)
      {
      }

      delete precip_data_ptr;

   }

}

void
Gfs::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Nwp_Element nwp_element)
{

   switch (nwp_element)
   {

      case RAINFALL_CUMULATIVE:
      {
         fill_cumulative_rain_data (gvd_2d, vector_index, key);
         return;
      }

      case RAINFALL_STEP:
      {
         fill_step_rain_data (gvd_2d, vector_index, key);
         return;
      }

   }

   Nwp::fill_rain_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Nwp_Element nwp_element)
{
   const Level& nil_level = Level::nil_level ();
   fill_grib_data (gvd_2d, vector_index, nwp_element, key, nil_level);
}

void
Gfs::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                             const Integer vector_index,
                             const Key& key,
                             const Nwp_Element nwp_element)
{

   const Level& screen = Level::screen_level ();
   const Level& surface = Level::surface_level ();

   switch (nwp_element)
   {

      case TEMPERATURE:
      {
         fill_grib_data (gvd_2d, vector_index, TEMPERATURE, key, screen);
         return;
      }

      case DEW_POINT:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, RELATIVE_HUMIDITY);

         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real rh = data_ptr->get_datum (1, i, j);
               const Real t_d = Moisture::get_t_d (t, rh);
               gvd_2d.set_datum (vector_index, i, j, t_d);
            }
         }

         delete data_ptr;
         return;

      }

      case RELATIVE_HUMIDITY:
      {
         fill_grib_data (gvd_2d, vector_index, RELATIVE_HUMIDITY, key, screen);
         return;
      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Nwp_Element nwp_element)
{

   const Level& ten = Level::ten_metre_level ();

   switch (nwp_element)
   {

      case ZONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, ZONAL_WIND, key, ten);
         return;
      }

      case MERIDIONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, MERIDIONAL_WIND, key, ten);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                    const Integer vector_index,
                    const Key& key,
                    const Nwp_Element nwp_element)
{

   const Level& msl = Level::mean_sea_level ();

   switch (nwp_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         const Nwp_Element mslp = MEAN_SEA_LEVEL_PRESSURE;
         fill_grib_data (gvd_2d, vector_index, mslp, key, msl);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, nwp_element);

}

void
Gfs::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Nwp_Element nwp_element)
{

   const Level& surface = Level::surface_level ();

   switch (nwp_element)
   {

      case PRESSURE:
      {
         fill_grib_data (gvd_2d, vector_index, PRESSURE, key, surface);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, nwp_element);

}

/*
void
Gfs::fill_grib_data (Geodetic_Vector_Data_3D& gvd_3d,
                     const Integer vector_index,
                     const Nwp_Element nwp_element,
                     const Key& key) const
{

   const Key key (key);

   Gfs::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib2& grib = *(iterator->second);

   for (Integer k = 0; k < tuple_p.size (); k++)
   {
      const Real p = tuple_p[k];
      const Level level (PRESSURE_LEVEL, p);
      const Grib2::Key& key = get_grib_key (key, nwp_element, level);
      grib.fill_data (gvd_3d, vector_index, k, key);
   }

}
*/

Geodetic_Vector_Data_3D*
Gfs::get_gvd_3d_ptr (const Nwp_Element nwp_element,
                     const Key& key) const
{

   Geodetic_Vector_Data_3D* gvd_3d_ptr =
      new Geodetic_Vector_Data_3D (1, tuple_p, size_2d, domain_2d);
   Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   for (Integer k = 0; k < tuple_p.size (); k++)
   {

      const Real p = tuple_p[k];
      const denise::Level level (PRESSURE_LEVEL, p);
      const Grib2::Key& grib_key = get_grib_key (key, nwp_element, level);

      Gfs::const_iterator iterator = find (key);
      if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
      const Grib2& grib = *(iterator->second);

      Geodetic_Vector_Data_2D* grib_data_ptr =
         get_grib_data_ptr (grib, grib_key);

      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         const Real latitude = grib_data_ptr->get_latitude (i);
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real longitude = grib_data_ptr->get_longitude (j);
            const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
            gvd_3d.set_datum (0, k, i, j, datum);
         }
      }

      delete grib_data_ptr;

   }

   return gvd_3d_ptr;

}

Geodetic_Vector_Data_2D*
Gfs::get_grib_data_ptr (const Grib2& grib,
                        const Grib2::Key& grib_key) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = new Gvd_2d (1, grib_size_2d, grib_domain_2d, true);
   grib.fill_data (*data_ptr, 0, grib_key);

   const Integer last_j = grib_size_2d.j - 1;

   #pragma omp parallel for
   for (Integer i = 0; i < grib_size_2d.i; i++)
   {
      const Real datum = data_ptr->get_datum (0, i, 0);
      data_ptr->set_datum (0, i, last_j, datum);
   }

   return data_ptr;

}

void
Gfs::fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                     const Integer vector_index,
                     const Nwp_Element nwp_element,
                     const Key& key,
                     const Level& level) const
{

   Gfs::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib2& grib = *(iterator->second);

   const Grib2::Key& grib_key = get_grib_key (key, nwp_element, level);
   Geodetic_Vector_Data_2D* grib_data_ptr = get_grib_data_ptr (grib, grib_key);

   Geodetic_Vector_Data_2D grib_data (1, grib_size_2d, grib_domain_2d, true);

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_latitude (i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_longitude (j);
         const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete grib_data_ptr;

}

Gfs::Gfs (const string& data_path)
   : Nwp ("Gfs", data_path),
     data_path (data_path),
     grib_size_2d (361, 720),
     grib_domain_2d (Domain_1D (-90, 90), Domain_1D (0, 360)),
     size_2d (51, 71),
     domain_2d (Domain_1D (-60, 0), Domain_1D (80, 160))
{

   this->status = "Unloaded";

   nwp_element_vector.push_back (denise::TEMPERATURE);
   nwp_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   nwp_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::ZONAL_WIND);
   nwp_element_vector.push_back (denise::MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::OMEGA);

}

Gfs::~Gfs ()
{
   clean_up ();
}

void
Gfs::survey ()
{

   const string re_ym ("[0-9].....");
   const string re_ymd ("[0-9].......");
   const string fmt ("gfs_4_[0-9]......._[0-9]..._[0-9]..\\.grb");

   typedef map<Grib2::Key, Grib2::Header*> Header_Ptr_Map;

   const vector<string>& dir_ym = get_dir_listing (path, re_ym);
   for (vector<string>::const_iterator i = dir_ym.begin ();
        i != dir_ym.end (); i++)
   {

      const string& ym = *(i);
      const string& ym_path = path + "/" + ym;

      const vector<string>& dir_ymd = get_dir_listing (ym_path, re_ymd);
      for (vector<string>::const_iterator j = dir_ymd.begin ();
           j != dir_ymd.end (); j++)
      {

         const string& ymd = *(j);
         const string& ymd_path = ym_path + "/" + ymd;

         const vector<string>& dir_listing = get_dir_listing (ymd_path, fmt);
         for (vector<string>::const_iterator iterator = dir_listing.begin ();
              iterator != dir_listing.end (); iterator++)
         {

            // fn = filename
            const string& fn = *(iterator);
            const string& file_path = ymd_path + "/" + fn;
            const string& bt_str = fn.substr (6, 8) + fn.substr (15, 2);
            const string& fh_str = fn.substr (20, 3);

            const Dtime base_time (bt_str, string ("%Y%m%d%H"));
            const Integer forecast_hour = atoi (fh_str.c_str ());

            const Key key (base_time, forecast_hour);
            key_multimap.add (key);
            Grib2* grib_ptr = new Grib2 (file_path);
            insert (make_pair (key, grib_ptr));
            //valid_time_set.insert (valid_time);
            //key_set.insert (Key (base_time, fh));

         }

      }

   }

   if (size () > 0)
   {

      set<uint32_t> set_p;
      const Grib2& grib = *(begin ()->second);
      const Header_Ptr_Map& header_ptr_map = grib.get_header_ptr_map ();

      for (Header_Ptr_Map::const_iterator iterator = header_ptr_map.begin ();
           iterator != header_ptr_map.end (); iterator++)
      {

         const Grib2::Header& header = *(iterator->second);
         const Grib2::Block_4& block_4 = header.block_4;
         const Grib2::Block_4::Level& level = block_4.get_level ();

         if (level.get_first_level_type () == 100)
         {
            const uint32_t p = level.get_first_level ();
            if (p < 200) { continue; }
            set_p.insert (p);
         }

      }

      for (set<uint32_t>::const_iterator iterator = set_p.begin ();
           iterator != set_p.end (); iterator++)
      {
         const Real p = Real (*(iterator));
         tuple_p.push_back (p);
      }

   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

}

void
Gfs::clean_up ()
{
   for (Gfs::iterator iterator = begin (); iterator != end (); iterator++)
   {
      Grib2* grib_ptr = iterator->second;
      delete grib_ptr;
   }
}

void
Gfs::set_domain_2d (const Domain_2D& domain_2d)
{

   this->domain_2d.domain_x.start = ceil (domain_2d.domain_x.start);
   this->domain_2d.domain_x.end = floor (domain_2d.domain_x.end);
   this->domain_2d.domain_y.start = ceil (domain_2d.domain_y.start);
   this->domain_2d.domain_y.end = floor (domain_2d.domain_y.end);

   const Real& start_latitude = this->domain_2d.domain_x.start;
   const Real& end_latitude = this->domain_2d.domain_x.end;
   const Real& start_longitude = this->domain_2d.domain_y.start;
   const Real& end_longitude = this->domain_2d.domain_y.end;
 
   size_2d.i = Integer (round (end_latitude - start_latitude));
   size_2d.j = Integer (round (end_longitude - start_longitude));

}

vector<Dtime>
Gfs::get_valid_time_vector () const
{

   vector<Dtime> valid_time_vector;

   for (Gfs::const_iterator iterator = begin ();
       iterator != end (); iterator++)
   {
      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime dtime (bt.t + fh);
      valid_time_vector.push_back (dtime);
   }

   return valid_time_vector;

/*
   vector<Dtime> valid_time_vector;
   Dtime start_time ("2012011300");
   Dtime end_time ("2012011400");
   for (Real t = start_time.t; t < end_time.t; t += 3)
   {
      Dtime dtime (t);
      valid_time_vector.push_back (dtime);
   }
   return valid_time_vector;
*/
}

Nwp::Key
Gfs::get_key (const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());

   const Dtime& base_time = base_time_vector.back ();
   const Integer forecast_hour = Integer (round (dtime.t - base_time.t));

   return Key (base_time, forecast_hour);

}

void
Gfs::acquire_base_time_forecast_hour (Dtime& base_time,
                                      Integer& forecast_hour,
                                      const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());
   base_time.t = base_time_vector.back ().t;
   forecast_hour = Integer (round (dtime.t - base_time.t));

}

Integer
Ncep_Ncar_File::get_dimension_size (const string& dimension_string) const
{

   int ret, dim_id;
   size_t dim_length;
   const char* str = dimension_string.c_str ();

   ret = nc_inq_dimid (nc_id, dimension_string.c_str (), &dim_id);
   if (ret != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_inq_dimid " + dimension_string); }

   ret = nc_inq_dimlen (nc_id, dim_id, &dim_length);
   if (ret != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_inq_dimlen " + dimension_string); }

   return Integer (dim_length);

}

Ncep_Ncar_File::Ncep_Ncar_File (const string& file_path)
                   : file_path (file_path)

{

   int ret, dim_id;
   size_t latitude_dim_length, longitude_dim_length;

   ret = nc_open (file_path.c_str (), NC_NOWRITE, &nc_id);
   if (ret != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_open " + file_path); }

}

Ncep_Ncar_File::~Ncep_Ncar_File ()
{
   int ret = nc_close (nc_id);
   if (ret != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_close " + file_path); }
}

Ncep_Ncar_3d_File::Ncep_Ncar_3d_File (const string& file_path)
   : Ncep_Ncar_File (file_path),
     actual_number_of_levels (get_dimension_size ("level"))
{
}

void
Ncep_Ncar_3d_File::fill_data (Geodetic_Vector_Data_3D& data_3d,
                              const Integer vector_index,
                              const string& variable_string,
                              const Integer time_index) const
{

   int r, var_id;
   float offset, scale;

   const Tuple tuple_level = data_3d.get_coordinate_tuple (0);
   const Tuple tuple_latitude = data_3d.get_coordinate_tuple (1);
   const Tuple tuple_longitude = data_3d.get_coordinate_tuple (2);

   const Integer nk = tuple_level.size ();
   const Integer ni = tuple_latitude.size ();
   const Integer nj = tuple_longitude.size () - 1;
   size_t n = nk * ni * nj;

   short* array = new short[n];
   size_t start[] = { time_index, 0, 0, 0 };
   size_t count[] = { 1, actual_number_of_levels, ni, nj };

   const string& vs = variable_string;
   typedef Ncep_Ncar_Exception Nne;

   r = nc_inq_varid (nc_id, vs.c_str (), &var_id);
   if (r != NC_NOERR) { throw Nne ("nc_inq_varid " + vs); }

   r = nc_get_vara_short (nc_id, var_id, start, count, array);
   if (r != NC_NOERR) { throw Nne ("nc_get_vara_short " + vs); }

   r = nc_get_att_float (nc_id, var_id, "add_offset", &offset);
   if (r != NC_NOERR) { throw Nne ("nc_get_att_float offset " + vs); }

   r = nc_get_att_float (nc_id, var_id, "scale_factor", &scale);
   if (r != NC_NOERR) { throw Nne ("nc_get_att_float scale " + vs); }

   //cout << "actual number of levels = " << actual_number_of_levels << endl;

   for (Integer kk = 0; kk < actual_number_of_levels; kk++)
   {

      const Integer k = nk - kk - 1;

      for (Integer i = 0; i < ni; i++)
      {

         const Integer ii = ni - i - 1;

         for (Integer j = 0; j < nj; j++)
         {
            const Integer index = (kk * ni + ii) * nj + j;
            const Real datum = array[index] * scale + offset;
            data_3d.set_datum (vector_index, k, i, j, datum);
         }

         const Real datum = data_3d.get_datum (vector_index, k, i, 0);
         data_3d.set_datum (vector_index, k, i, nj, datum);

      }
   }

   if (variable_string == "rhum")
   {
      data_3d.scale_offset (vector_index, 0.01, 0);
   }

   if (variable_string == "omega")
   {
      data_3d.scale_offset (vector_index, -10, 0);
   }

   delete[] array;

}

void
Ncep_Ncar_3d_File::fill_data (Geodetic_Vector_Data_2D& data_2d,
                              const Integer vector_index,
                              const string& variable_string,
                              const Integer time_index,
                              const Integer level_index) const
{

   int r, var_id;
   float offset, scale;

   const Tuple tuple_latitude = data_2d.get_coordinate_tuple (0);
   const Tuple tuple_longitude = data_2d.get_coordinate_tuple (1);

   const Integer ni = tuple_latitude.size ();
   const Integer nj = tuple_longitude.size () - 1;
   size_t n = ni * nj;

   short* array = new short[n];
   size_t start[] = { time_index, level_index, 0, 0 };
   size_t count[] = { 1, 1, ni, nj };

   const string& vs = variable_string;

   r = nc_inq_varid (nc_id, vs.c_str (), &var_id);
   if (r != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_inq_varid " + vs); }

   r = nc_get_vara_short (nc_id, var_id, start, count, array);
   if (r != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_get_vara_short " + vs); }

   r = nc_get_att_float (nc_id, var_id, "add_offset", &offset);
   if (r != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_get_att_float offset " + vs); }

   r = nc_get_att_float (nc_id, var_id, "scale_factor", &scale);
   if (r != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_get_att_float scale " + vs); }

   #pragma omp parallel for
   for (Integer i = 0; i < ni; i++)
   {

      const Integer ii = ni - i - 1;

      for (Integer j = 0; j < nj; j++)
      {
         const Integer index = ii * nj + j;
         const Real datum = array[index] * scale + offset;
         data_2d.set_datum (vector_index, i, j, datum);
      }

      const Real datum = data_2d.get_datum (vector_index, i, 0);
      data_2d.set_datum (vector_index, i, nj, datum);

   }

   delete[] array;

}

Ncep_Ncar_2d_File::Ncep_Ncar_2d_File (const string& file_path)
                              : Ncep_Ncar_File (file_path)
{
}

void
Ncep_Ncar_2d_File::fill_data (Geodetic_Vector_Data_2D& data_2d,
                              const Integer vector_index,
                              const string& variable_string,
                              const Integer time_index) const
{

   int r, var_id;
   float offset, scale;

   const Tuple tuple_latitude = data_2d.get_coordinate_tuple (0);
   const Tuple tuple_longitude = data_2d.get_coordinate_tuple (1);

   const Integer ni = tuple_latitude.size ();
   const Integer nj = tuple_longitude.size () - 1;
   size_t n = ni * nj;

   short* array = new short[n];
   size_t start[] = { time_index, 0, 0 };
   size_t count[] = { 1, ni, nj };

   const string& vs = variable_string;

   r = nc_inq_varid (nc_id, vs.c_str (), &var_id);
   if (r != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_inq_varid " + vs); }

   r = nc_get_vara_short (nc_id, var_id, start, count, array);
   if (r != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_get_vara_short " + vs); }

   r = nc_get_att_float (nc_id, var_id, "add_offset", &offset);
   if (r != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_get_att_float offset " + vs); }

   r = nc_get_att_float (nc_id, var_id, "scale_factor", &scale);
   if (r != NC_NOERR) { throw Ncep_Ncar_Exception ("nc_get_att_float scale " + vs); }

   #pragma omp parallel for
   for (Integer i = 0; i < ni; i++)
   {

      const Integer ii = ni - i - 1;

      for (Integer j = 0; j < nj; j++)
      {
         const Integer index = ii * nj + j;
         const Real datum = array[index] * scale + offset;
         data_2d.set_datum (vector_index, i, j, datum);
      }

      const Real datum = data_2d.get_datum (vector_index, i, 0);
      data_2d.set_datum (vector_index, i, nj, datum);

   }

   if (variable_string == "rhum")
   {
      data_2d.scale_offset (vector_index, 0.01, 0);
   }

   delete[] array;

}

Ncep_Ncar_Exception::Ncep_Ncar_Exception (const string& description)
   : Exception ("Ncep_Ncar_Exception", description)
{
}

Ncep_Ncar::Data_3D::Data_3D (const vector<Nwp_Element>& nwp_element_vector,
                             const Key& key)
   : Nwp::Data_3D (nwp_element_vector, key)
{
}

Real
Ncep_Ncar::Data_3D::evaluate (const Nwp_Element element,
                              const Real p,
                              const Real latitude,
                              const Real longitude,
                              const Evaluate_Op evaluate_op) const
{

   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& RH = denise::RELATIVE_HUMIDITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real rh = Nwp::Data_3D::evaluate (RH, p, latitude, longitude);
         return Moisture::get_t_d (t, std::max (0.01, rh));
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& O = denise::OMEGA;

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
            const Real o = Nwp::Data_3D::evaluate (O, p, latitude, longitude);
            Real ts = Nwp::Data_3D::evaluate (T, p, latitude, longitude, eo);
            Real os = Nwp::Data_3D::evaluate (O, p, latitude, longitude, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
            const Real o = Nwp::Data_3D::evaluate (O, p, latitude, longitude);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nwp::Data_3D::evaluate (element, p, latitude, longitude, eo);

}

string
Ncep_Ncar::get_file_path (const Element element,
                          const Dtime& time) const
{

   string identifier;

   switch (element)
   {
      case ncepncar::U:             identifier = "uwnd.";          break;
      case ncepncar::V:             identifier = "vwnd.";          break;
      case ncepncar::OMG:           identifier = "omega.";         break;
      case ncepncar::T:             identifier = "air.";           break;
      case ncepncar::RH:            identifier = "rhum.";          break;
      case ncepncar::SH:            identifier = "shum.";          break;
      case ncepncar::Z:             identifier = "hgt.";           break;
      case ncepncar::SURFACE_U:     identifier = "uwnd.sig995.";   break;
      case ncepncar::SURFACE_V:     identifier = "vwnd.sig995.";   break;
      case ncepncar::SURFACE_OMG:   identifier = "omega.sig995.";  break;
      case ncepncar::SURFACE_P:     identifier = "pres.sfc.";      break;
      case ncepncar::SURFACE_T:     identifier = "air.sig995.";    break;
      case ncepncar::SURFACE_RH:    identifier = "rhum.sig995.";   break;
      case ncepncar::SURFACE_THETA: identifier = "pottmp.sig995."; break;
      case ncepncar::PRECIP_WATER:  identifier = "pr_wtr.eatm.";   break;
      case ncepncar::SLP:           identifier = "slp.";           break;
      case ncepncar::U_10M:         identifier = "uwnd.10m.gauss."; break;
      case ncepncar::V_10M:         identifier = "vwnd.10m.gauss."; break;
      case ncepncar::T_2M:          identifier = "air.2m.";        break;
   }

   const Dtime quantized_time ((rint (time.t / 6) * 6));
   const string year_string = quantized_time.get_string ("%Y");
   const string path_base (data_path + "/" + year_string + "/");

   const string file_path = path_base + identifier + year_string + ".nc";
   return file_path;

}

string
Ncep_Ncar::get_element_string (const Element element) const
{

   string str;

   switch (element)
   {
      case ncepncar::U:             str = "uwnd";   break;
      case ncepncar::V:             str = "vwnd";   break;
      case ncepncar::OMG:           str = "omega";  break;
      case ncepncar::T:             str = "air";    break;
      case ncepncar::RH:            str = "rhum";   break;
      case ncepncar::SH:            str = "shum";   break;
      case ncepncar::Z:             str = "hgt";    break;
      case ncepncar::SURFACE_U:     str = "uwnd";   break;
      case ncepncar::SURFACE_V:     str = "vwnd";   break;
      case ncepncar::SURFACE_OMG:   str = "omega";  break;
      case ncepncar::SURFACE_P:     str = "pres";   break;
      case ncepncar::SURFACE_T:     str = "air";    break;
      case ncepncar::SURFACE_RH:    str = "rhum";   break;
      case ncepncar::SURFACE_THETA: str = "pottmp"; break;
      case ncepncar::PRECIP_WATER:  str = "pr_wtr"; break;
      case ncepncar::SLP:           str = "slp";    break;
      case ncepncar::U_10M:         str = "uwnd";   break;
      case ncepncar::V_10M:         str = "vwnd";   break;
      case ncepncar::T_2M:          str = "air";    break;
   }

   return str;

}

Integer
Ncep_Ncar::get_time_index (const Dtime& time) const
{

   const Integer day_of_year = time.get_day_of_year ();
   const Integer hour = time.get_hour ();

   return Integer (rint (Real ((day_of_year - 1) * 24 + hour) / 6));

}

Integer
Ncep_Ncar::get_p_index (const Element element,
                        const Real p) const
{

   const Tuple tuple_p = get_tuple_p (element);
   const Integer n = tuple_p.size ();

   for (Integer i = 0; i < n; i++)
   {
      const Real delta_p = fabs (tuple_p[i] - p);
      if (delta_p < 0.1) { return n - i - 1; }
   }

   throw Ncep_Ncar_Exception ("Level not available");

}

const Tuple&
Ncep_Ncar::get_tuple_p (const Element element) const
{

   switch (element)
   {

      case ncepncar::U:
      case ncepncar::V:
      case ncepncar::T:
      case ncepncar::Z:
      {
         return tuple_p;
      }

      case ncepncar::OMG:
      {
         return tuple_p_omega;
      }

      case ncepncar::RH:
      case ncepncar::SH:
      {
         return tuple_p_humid;
      }

      default:
      {
         throw Ncep_Ncar_Exception ("Not 3D Element");
      }

   }

}

Geodetic_Vector_Data_2D*
Ncep_Ncar::get_initialized_vd_2d (const Integer vector_size) const
{
   typedef Geodetic_Vector_Data_2D Gvd_2d;
   const Tuple tuple_latitude = Ncep_Ncar::tuple_latitude ();
   const Tuple tuple_longitude = Ncep_Ncar::tuple_longitude ();
   return new Gvd_2d (vector_size, tuple_latitude, tuple_longitude, true);
}

void
Ncep_Ncar::initialize_3d_data (const Key& key)
{

   typedef Ncep_Ncar::Data_3D Nnd_3d;

   Nnd_3d* nnd_3d_ptr = new Nnd_3d (nwp_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, nnd_3d_ptr));

}

void          
Ncep_Ncar::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Nwp_Element>::const_iterator Iterator;
   typedef Geodetic_Vector_Data_3D Gvd_3d;
   const Key& key = data_3d.key;
   const Dtime& dtime = key.base_time;
   const Tuple tuple_latitude = Ncep_Ncar::tuple_latitude ();
   const Tuple tuple_longitude = Ncep_Ncar::tuple_longitude ();

   for (Iterator iterator = nwp_element_vector.begin ();
        iterator != nwp_element_vector.end (); iterator++)
   {

      const Nwp_Element& nwp_element = *(iterator);

      switch (nwp_element)
      {

         case TEMPERATURE:
         {
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, ncepncar::T, dtime);
            data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
            break;
         }

         case RELATIVE_HUMIDITY:
         {
            Tuple tuple_p;
            tuple_p.add_content ("300e2:400e2:500e2:600e2");
            tuple_p.add_content ("700e2:850e2:925e2:1000e2");
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, ncepncar::RH, dtime);
            data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
            break;
         }

         case GEOPOTENTIAL_HEIGHT:
         {
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, ncepncar::Z, dtime);
            data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
            break;
         }

         case ZONAL_WIND:
         {
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, ncepncar::U, dtime);
            data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
            break;
         }

         case MERIDIONAL_WIND:
         {
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, ncepncar::V, dtime);
            data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
            break;
         }

         case OMEGA:
         {
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, ncepncar::OMG, dtime);
            gvd_3d.scale_offset (0, -1, 0);
            data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
            break;
         }

      }

   }

   data_3d.set_available ();

}

void
Ncep_Ncar::fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                      const Integer vector_index,
                                      const Key& key)
{
   gvd_2d.initialize (vector_index, 0);
}

void
Ncep_Ncar::fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key)
{
   gvd_2d.initialize (vector_index, 0);
}

void
Ncep_Ncar::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                   const Integer vector_index,
                                   const Key& key,
                                   const Nwp_Element nwp_element)
{

   const Level& screen = Level::screen_level ();
   const Level& surface = Level::surface_level ();
   const Dtime& dtime = key.base_time;

   switch (nwp_element)
   {

      case TEMPERATURE:
      {
         fill_nc_data (gvd_2d, vector_index, SURFACE_T, dtime);
         return;
      }

      case DEW_POINT:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);

         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, RELATIVE_HUMIDITY);

         const Size_2D& size_2d = data_ptr->get_size_2d ();

         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.i; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real rh = data_ptr->get_datum (1, i, j);
               const Real t_d = Moisture::get_t_d (t, std::max (0.01, rh));
               gvd_2d.set_datum (vector_index, i, j, t_d);
            }
         }

         delete data_ptr;
         return;

      }

      case RELATIVE_HUMIDITY:
      {
         fill_nc_data (gvd_2d, vector_index, SURFACE_RH, dtime);
         return;
      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ncep_Ncar::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Nwp_Element nwp_element)
{

   const Dtime& dtime = key.base_time;

   switch (nwp_element)
   {

      case ZONAL_WIND:
      {
         fill_nc_data (gvd_2d, vector_index, SURFACE_U, dtime);
         return;
      }

      case MERIDIONAL_WIND:
      {
         fill_nc_data (gvd_2d, vector_index, SURFACE_V, dtime);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ncep_Ncar::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Nwp_Element nwp_element)
{

   const Dtime& dtime = key.base_time;

   switch (nwp_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         fill_nc_data (gvd_2d, vector_index, SLP, dtime);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ncep_Ncar::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                    const Integer vector_index,
                                    const Key& key,
                                    const Nwp_Element nwp_element)
{

   const Dtime& dtime = key.base_time;

   switch (nwp_element)
   {

      case PRESSURE:
      {
         fill_nc_data (gvd_2d, vector_index, SURFACE_P, dtime);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ncep_Ncar::fill_nc_data (Geodetic_Vector_Data_3D& gvd_3d,
                         const Integer vector_index,
                         const Element element,
                         const Dtime& time) const
{

   const string file_path = get_file_path (element, time);

   const string element_string = get_element_string (element);
   const Integer time_index = get_time_index (time);

   switch (element)
   {

      case ncepncar::U:
      case ncepncar::V:
      case ncepncar::OMG:
      case ncepncar::T:
      case ncepncar::RH:
      case ncepncar::SH:
      case ncepncar::Z:
      {
         const Ncep_Ncar_3d_File file (file_path);
         file.fill_data (gvd_3d, vector_index, element_string, time_index);
         break;
      }

   }

}

void
Ncep_Ncar::fill_nc_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Element element,
                         const Dtime& time,
                         const Real p) const
{

   const string file_path = get_file_path (element, time);

   const string element_string = get_element_string (element);
   const Integer time_index = get_time_index (time);

   switch (element)
   {

      case ncepncar::U:
      case ncepncar::V:
      case ncepncar::OMG:
      case ncepncar::T:
      case ncepncar::RH:
      case ncepncar::SH:
      case ncepncar::Z:
      {
         const Ncep_Ncar_3d_File file (file_path);
         const Integer pi = get_p_index (element, p);
         file.fill_data (gvd_2d, vector_index, element_string, time_index, pi);
         break;
      }

      case ncepncar::SURFACE_U:
      case ncepncar::SURFACE_V:
      case ncepncar::SURFACE_OMG:
      case ncepncar::SURFACE_P:
      case ncepncar::SURFACE_T:
      case ncepncar::SURFACE_RH:
      case ncepncar::SURFACE_THETA:
      case ncepncar::PRECIP_WATER:
      case ncepncar::SLP:
      case ncepncar::U_10M:
      case ncepncar::V_10M:
      case ncepncar::T_2M:
      {
         const Ncep_Ncar_2d_File file (file_path);
         file.fill_data (gvd_2d, vector_index, element_string, time_index);
         break;
      }

   }

}

void
Ncep_Ncar::init ()
{

   tuple_longitude_gaussian.add_content (193, Real (0.0), Real (360.0));

   Tuple& tuple_lat_g = tuple_latitude_gaussian;
   tuple_lat_g.add_content ("-88.542:-86.6531:-84.7532:-82.8508:-80.9473");
   tuple_lat_g.add_content ("-79.0435:-77.1394:-75.2351:-73.3307:-71.4262");
   tuple_lat_g.add_content ("-69.5217:-67.6171:-65.7125:-63.8079:-61.9033");
   tuple_lat_g.add_content ("-59.9986:-58.0939:-56.1893:-54.2846:-52.3799");
   tuple_lat_g.add_content ("-50.4752:-48.5705:-46.6658:-44.7611:-42.8564");
   tuple_lat_g.add_content ("-40.9517:-39.047:-37.1422:-35.2375:-33.3328");
   tuple_lat_g.add_content ("-31.4281:-29.5234:-27.6186:-25.7139:-23.8092");
   tuple_lat_g.add_content ("-21.9044:-19.9997:-18.095:-16.1902:-14.2855");
   tuple_lat_g.add_content ("-12.3808:-10.47604:-8.57131:-6.66657:-4.76184");
   tuple_lat_g.add_content ("-2.8571:-0.952368:0.952368:2.8571:4.76184");
   tuple_lat_g.add_content ("6.66657:8.57131:10.47604:12.3808:14.2855");
   tuple_lat_g.add_content ("16.1902:18.095:19.9997:21.9044:23.8092");
   tuple_lat_g.add_content ("25.7139:27.6186:29.5234:31.4281:33.3328");
   tuple_lat_g.add_content ("35.2375:37.1422:39.047:40.9517:42.8564");
   tuple_lat_g.add_content ("44.7611:46.6658:48.5705:50.4752:52.3799");
   tuple_lat_g.add_content ("54.2846:56.1893:58.0939:59.9986:61.9033");
   tuple_lat_g.add_content ("63.8079:65.7125:67.6171:69.5217:71.4262");
   tuple_lat_g.add_content ("73.3307:75.2351:77.1394:79.0435:80.9473");
   tuple_lat_g.add_content ("82.8508:84.7532:86.6531:88.542");

   tuple_p.add_content ("10e2:20e2:30e2:50e2:70e2:100e2:150e2");
   tuple_p.add_content ("200e2:250e2:300e2:400e2:500e2:600e2");
   tuple_p.add_content ("700e2:850e2:925e2:1000e2");

   tuple_p_omega.add_content ("100e2:150e2:200e2:250e2:300e2:400e2");
   tuple_p_omega.add_content ("500e2:600e2:700e2:850e2:925e2:1000e2");

   tuple_p_humid.add_content ("300e2:400e2:500e2:600e2:700e2");
   tuple_p_humid.add_content ("850e2:925e2:1000e2");

}

Ncep_Ncar::Ncep_Ncar (const string& data_path,
                      const Dtime& start_time,
                      const Real time_span,
                      const Real time_interval)
   : Nwp ("Ncep_Ncar", data_path),
     data_path (data_path),
     start_time (start_time),
     time_span (time_span),
     time_interval (time_interval)
{

   init ();

   nwp_element_vector.push_back (denise::TEMPERATURE);
   nwp_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   nwp_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::ZONAL_WIND);
   nwp_element_vector.push_back (denise::MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::OMEGA);

   survey ();

}

Ncep_Ncar::Ncep_Ncar (const string& data_path,
                      const Dtime& start_time,
                      const Dtime& end_time,
                      const Real time_interval)
   : Nwp ("Ncep_Ncar", data_path),
     data_path (data_path),
     start_time (start_time),
     time_span (end_time.t - start_time.t),
     time_interval (time_interval)
{

   init ();

   nwp_element_vector.push_back (denise::TEMPERATURE);
   nwp_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   nwp_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::ZONAL_WIND);
   nwp_element_vector.push_back (denise::MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::OMEGA);

   survey ();

}

void
Ncep_Ncar::survey ()
{

   for (Real t = 0; t <= time_span; t += time_interval)
   {
      const Dtime dtime (start_time.t + t);
      const Key key (dtime, 0);
      key_multimap.add (key);
   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

   for (Key_Multimap::const_iterator iterator = key_multimap.begin ();
        iterator != key_multimap.end (); iterator++)
   {
      const Key key (iterator->first, iterator->second);
      initialize_3d_data (key);
   }

}

void
Ncep_Ncar::clean_up ()
{
   clear_data_3d_ptr_map ();
}

set<Dtime>
Ncep_Ncar::get_valid_time_set () const
{

   set<Dtime> valid_time_set;

   for (Real t = start_time.t;
        t < start_time.t + time_span;
        t += time_interval)
   {
      Dtime dtime (t);
      valid_time_set.insert (dtime);
   }

   return valid_time_set;

}

Tuple
Ncep_Ncar::tuple_latitude ()
{
   Tuple tuple_latitude;
   tuple_latitude.add_content (1, Real (-89.99));
   tuple_latitude.add_content (71, Real (-87.5), Real (87.5));
   tuple_latitude.add_content (1, Real (89.999));
   return tuple_latitude;
}

Tuple
Ncep_Ncar::tuple_longitude ()
{
   return Tuple (145, Real (0.0), Real (360.0));
}

const Tuple
Ncep_Ncar::get_tuple_latitude () const
{
   return Ncep_Ncar::tuple_latitude ();
}

const Tuple
Ncep_Ncar::get_tuple_longitude () const
{
   return Ncep_Ncar::tuple_longitude ();
}

Nwp::Key
Ncep_Ncar::get_key (const Dtime& dtime) const
{
   return Key (dtime, 0);
}

void
Ncep_Ncar::acquire_base_time_forecast_hour (Dtime& base_time,
                                            Integer& forecast_hour,
                                            const Dtime& dtime) const
{
   base_time.t = dtime.t;
   forecast_hour = 0;
}

Geodetic_Vector_Data_2D*
Ncep_Ncar::get_surface_data_2d_ptr (const Key& key)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* nn_data_2d_ptr = get_initialized_vd_2d (6);
   Gvd_2d& nn_data_2d = *nn_data_2d_ptr;

   const Level& msl = Level::mean_sea_level ();
   const Level& screen = Level::screen_level ();
   const Level& surface = Level::surface_level ();
   const Level& ten = Level::ten_metre_level ();

   fill_data (nn_data_2d, 0, key, screen, TEMPERATURE);
   fill_data (nn_data_2d, 1, key, screen, DEW_POINT);
   fill_data (nn_data_2d, 2, key, ten, ZONAL_WIND);
   fill_data (nn_data_2d, 3, key, ten, MERIDIONAL_WIND);
   fill_data (nn_data_2d, 4, key, msl, PRESSURE);
   fill_data (nn_data_2d, 5, key, surface, PRESSURE);

   return nn_data_2d_ptr;

}

Ncep_Ncar_Climate::Data_3D::Data_3D (const vector<Nwp_Element>& nwp_element_vector,
                                     const Key& key)
   : Nwp::Data_3D (nwp_element_vector, key)
{
}

Real
Ncep_Ncar_Climate::Data_3D::evaluate (const Nwp_Element element,
                                      const Integer k,
                                      const Integer i,
                                      const Integer j,
                                      const Evaluate_Op evaluate_op) const
{

   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::RELATIVE_HUMIDITY:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& TD = denise::DEW_POINT;
         const Real t = Nwp::Data_3D::evaluate (T, k, i, j);
         const Real td = Nwp::Data_3D::evaluate (TD, k, i, j);
         return Moisture::get_rh (t - K, td - K);
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& O = denise::OMEGA;
         const Real p = get_p (T, k);

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Real t = Nwp::Data_3D::evaluate (T, k, i, j);
            const Real o = Nwp::Data_3D::evaluate (O, k, i, j);
            Real ts = Nwp::Data_3D::evaluate (T, k, i, j, eo);
            Real os = Nwp::Data_3D::evaluate (O, k, i, j, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nwp::Data_3D::evaluate (T, k, i, j);
            const Real o = Nwp::Data_3D::evaluate (O, k, i, j);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nwp::Data_3D::evaluate (element, k, i, j, eo);

}

Real
Ncep_Ncar_Climate::Data_3D::evaluate (const Nwp_Element element,
                                      const Real p,
                                      const Integer i,
                                      const Integer j,
                                      const Evaluate_Op evaluate_op) const
{

   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::RELATIVE_HUMIDITY:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& TD = denise::DEW_POINT;
         const Real t = Nwp::Data_3D::evaluate (T, p, i, j);
         const Real td = Nwp::Data_3D::evaluate (TD, p, i, j);
         return Moisture::get_rh (t - K, td - K);
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& O = denise::OMEGA;

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Evaluate_Op& eo = evaluate_op;
            const Real t = Nwp::Data_3D::evaluate (T, p, i, j);
            const Real o = Nwp::Data_3D::evaluate (O, p, i, j);
            Real ts = Nwp::Data_3D::evaluate (T, p, i, j, eo);
            Real os = Nwp::Data_3D::evaluate (O, p, i, j, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nwp::Data_3D::evaluate (T, p, i, j);
            const Real o = Nwp::Data_3D::evaluate (O, p, i, j);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nwp::Data_3D::evaluate (element, p, i, j, eo);

}

Real
Ncep_Ncar_Climate::Data_3D::evaluate (const Nwp_Element element,
                                      const Integer k,
                                      const Real latitude,
                                      const Real longitude,
                                      const Evaluate_Op evaluate_op) const
{

   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::RELATIVE_HUMIDITY:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& TD = denise::DEW_POINT;
         const Real t = Nwp::Data_3D::evaluate (T, k, latitude, longitude);
         const Real td = Nwp::Data_3D::evaluate (TD, k, latitude, longitude);
         return Moisture::get_rh (t - K, td - K);
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& O = denise::OMEGA;
         const Real p = get_p (T, k);

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Real t = Nwp::Data_3D::evaluate (T, k, latitude, longitude);
            const Real o = Nwp::Data_3D::evaluate (O, k, latitude, longitude);
            Real ts = Nwp::Data_3D::evaluate (T, k, latitude, longitude, eo);
            Real os = Nwp::Data_3D::evaluate (O, k, latitude, longitude, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nwp::Data_3D::evaluate (T, k, latitude, longitude);
            const Real o = Nwp::Data_3D::evaluate (O, k, latitude, longitude);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nwp::Data_3D::evaluate (element, k, latitude, longitude, eo);

}

Real
Ncep_Ncar_Climate::Data_3D::evaluate (const Nwp_Element element,
                                      const Real p,
                                      const Real latitude,
                                      const Real longitude,
                                      const Evaluate_Op evaluate_op) const
{

   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::RELATIVE_HUMIDITY:
      {
         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& TD = denise::DEW_POINT;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real td = Nwp::Data_3D::evaluate (TD, p, latitude, longitude);
         return Moisture::get_rh (t - K, td - K);
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Nwp_Element& T = denise::TEMPERATURE;
         const Nwp_Element& O = denise::OMEGA;

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
            const Real o = Nwp::Data_3D::evaluate (O, p, latitude, longitude);
            Real ts = Nwp::Data_3D::evaluate (T, p, latitude, longitude, eo);
            Real os = Nwp::Data_3D::evaluate (O, p, latitude, longitude, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
            const Real o = Nwp::Data_3D::evaluate (O, p, latitude, longitude);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nwp::Data_3D::evaluate (element, p, latitude, longitude, eo);

}

void
Ncep_Ncar_Climate::Work_Data::sum_up_2d (const Geodetic_Vector_Data_2D& nn_data)
{

   const Integer in = Ncep_Ncar::tuple_latitude ().size ();
   const Integer jn = Ncep_Ncar::tuple_longitude ().size ();

   #pragma omp parallel for
   for (Integer i = 0; i < in; i++)
   {

      for (Integer j = 0; j < jn; j++)
      {

         const Real t = nn_data.get_datum (0, i, j);
         const Real td = nn_data.get_datum (1, i, j);
         const Real u = nn_data.get_datum (2, i, j);
         const Real v = nn_data.get_datum (3, i, j);
         const Real mslp = nn_data.get_datum (4, i, j);
         const Real surfp = nn_data.get_datum (5, i, j);

         data_2d_ptr->get_datum (0, i, j) += t;
         data_2d_ptr->get_datum (1, i, j) += td;
         data_2d_ptr->get_datum (2, i, j) += u;
         data_2d_ptr->get_datum (3, i, j) += v;
         data_2d_ptr->get_datum (4, i, j) += mslp;
         data_2d_ptr->get_datum (5, i, j) += surfp;

         data_2d_ptr->get_datum (0+6, i, j) += t * t;
         data_2d_ptr->get_datum (1+6, i, j) += td * td;
         data_2d_ptr->get_datum (2+6, i, j) += u * u;
         data_2d_ptr->get_datum (3+6, i, j) += v * v;
         data_2d_ptr->get_datum (4+6, i, j) += mslp * mslp;
         data_2d_ptr->get_datum (5+6, i, j) += surfp * surfp;

      }

   }

}

void
Ncep_Ncar_Climate::Work_Data::sum_up_3d (Nwp::Data_3D& nn_data,
                                         const Tuple& tuple_p)
{

   const Integer pn = tuple_p.size ();
   const Integer in = Ncep_Ncar::tuple_latitude ().size ();
   const Integer jn = Ncep_Ncar::tuple_longitude ().size ();

   vector<Nwp_Element> element_vector;
   element_vector.push_back (TEMPERATURE);
   element_vector.push_back (DEW_POINT);
   element_vector.push_back (GEOPOTENTIAL_HEIGHT);
   element_vector.push_back (ZONAL_WIND);
   element_vector.push_back (MERIDIONAL_WIND);
   element_vector.push_back (OMEGA);
   const Integer ne = element_vector.size ();

   #pragma omp parallel for
   for (Integer k = 0; k < pn; k++)
   {
      const Real p = tuple_p[k];
      for (Integer i = 0; i < in; i++)
      {
         for (Integer j = 0; j < jn; j++)
         {
            for (Integer e = 0; e < ne; e++)
            {
               const Nwp_Element& nwp_element = element_vector[e];
               const Real datum = nn_data.evaluate (nwp_element, p, i, j);
               data_3d_ptr->get_datum (e, k, i, j) += datum;
               data_3d_ptr->get_datum (e+ne, k, i, j) += datum * datum;
            }
         }
      }
   }

}

void
Ncep_Ncar_Climate::Work_Data::divide_2d ()
{

   const Integer ne = 6;

   #pragma omp parallel for
   for (Integer i = 0; i < size_3d.i; i++)
   {
      for (Integer j = 0; j < size_3d.j; j++)
      {
         for (Integer e = 0; e < ne; e++)
         {
            Real& mean = data_2d_ptr->get_datum (e, i, j);
            Real& variance = data_2d_ptr->get_datum (e+ne, i, j);
            mean /= Real (n);
            variance = variance / Real (n) - mean;
         }
      }
   }

}

void
Ncep_Ncar_Climate::Work_Data::divide_3d ()
{

   const Integer ne = 6;

   for (Integer k = 0; k < size_3d.k; k++)
   {
      #pragma omp parallel for
      for (Integer i = 0; i < size_3d.i; i++)
      {
         for (Integer j = 0; j < size_3d.j; j++)
         {
            for (Integer e = 0; e < ne; e++)
            {
               Real& mean = data_3d_ptr->get_datum (e, k, i, j);
               Real& variance = data_3d_ptr->get_datum (e+ne, k, i, j);
               mean /= Real (n);
               variance = variance / Real (n) - mean;
            }
         }
      }
   }

}

Ncep_Ncar_Climate::Work_Data::Work_Data (const Tuple& tuple_p,
                                         const Tuple& tuple_latitude,
                                         const Tuple& tuple_longitude)
   : n (0),
     size_3d (tuple_p.size (), tuple_latitude.size (), tuple_longitude.size ())
{
   data_2d_ptr = new Geodetic_Vector_Data_2D (12,
      tuple_latitude, tuple_longitude, true);
   data_3d_ptr = new Geodetic_Vector_Data_3D (12, tuple_p,
      tuple_latitude, tuple_longitude, true);
}

Ncep_Ncar_Climate::Work_Data::~Work_Data ()
{
   delete data_2d_ptr;
   delete data_3d_ptr;
}

void
Ncep_Ncar_Climate::Work_Data::write (const string& path,
                                     const Integer jjj) const
{

   const string& file_name = string_render ("nnc_%03d", jjj);
   const string& file_path = path + "/" + file_name;
   ofstream file (file_path.c_str ());

   data_2d_ptr->write (file);
   data_3d_ptr->write (file);

   file.close ();

}

void
Ncep_Ncar_Climate::Work_Data::sum_up (const Geodetic_Vector_Data_2D& nn_data_2d,
                                      Nwp::Data_3D& nn_data_3d,
                                      const Tuple& tuple_p)
{
   n++;
   sum_up_2d (nn_data_2d);
   sum_up_3d (nn_data_3d, tuple_p);
}

void
Ncep_Ncar_Climate::Work_Data::divide ()
{
   divide_2d ();
   divide_3d ();
}

void
Ncep_Ncar_Climate::Work_Data_Map::sum_up (Ncep_Ncar& ncep_ncar,
                                          const Integer tolerance)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   typedef Geodetic_Vector_Data_3D Gvd_3d;
   const Tuple& tuple_latitude = Ncep_Ncar::tuple_latitude ();
   const Tuple& tuple_longitude = Ncep_Ncar::tuple_longitude ();
   const set<Dtime> time_set = ncep_ncar.get_valid_time_set ();
   const Integer in = tuple_latitude.size ();
   const Integer jn = tuple_longitude.size ();

   vector<Nwp_Element> element_vector;
   element_vector.push_back (TEMPERATURE);
   element_vector.push_back (DEW_POINT);
   element_vector.push_back (PRESSURE);
   element_vector.push_back (ZONAL_WIND);
   element_vector.push_back (MERIDIONAL_WIND);
   element_vector.push_back (MEAN_SEA_LEVEL_PRESSURE);
   const Integer ne = element_vector.size ();

   for (set<Dtime>::const_iterator iterator_t = time_set.begin ();
        iterator_t != time_set.end (); iterator_t++)
   {

      const Dtime& dtime = *(iterator_t);
      const Key key (dtime, 0);

      Nwp::Data_3D& nn_data_3d = ncep_ncar.get_3d_data (key);
      Gvd_2d* nn_data_2d_ptr = ncep_ncar.get_surface_data_2d_ptr (key);
      Gvd_2d& nn_data_2d = *nn_data_2d_ptr;

      const Integer day_of_year = dtime.get_day_of_year ();

      for (Ituple::const_iterator iterator_j = ituple_jjj.begin ();
           iterator_j != ituple_jjj.end (); iterator_j++)
      {

         const Integer jjj = *(iterator_j);
         const Integer delta_j = abs (jjj - day_of_year);
         const bool criteria_a = delta_j < (365 - tolerance);
         const bool criteria_b = delta_j > (tolerance);

         if (criteria_a && criteria_b) { continue; }

         Work_Data& work_data = *(at (jjj));
         work_data.sum_up (nn_data_2d, nn_data_3d, tuple_p);

      }

      delete nn_data_2d_ptr;
      nn_data_3d.unload ();

   }

}

void
Ncep_Ncar_Climate::Work_Data_Map::divide ()
{

   for (Ituple::const_iterator iterator_j = ituple_jjj.begin ();
        iterator_j != ituple_jjj.end (); iterator_j++)
   {
      const Integer jjj = *(iterator_j);
      Work_Data& work_data = *(at (jjj));
      work_data.divide ();
   }

}

Ncep_Ncar_Climate::Work_Data_Map::Work_Data_Map (const Ituple& ituple_jjj,
                                                 const Tuple& tuple_p,
                                                 const Tuple& tuple_latitude,
                                                 const Tuple& tuple_longitude)
   : ituple_jjj (ituple_jjj),
     tuple_p (tuple_p)
{

   for (Ituple::const_iterator iterator = ituple_jjj.begin ();
        iterator != ituple_jjj.end (); iterator++)
   {

      const Integer jjj = *(iterator);

      Work_Data* work_data_ptr = new Work_Data (tuple_p,
         tuple_latitude, tuple_longitude);
      insert (make_pair (jjj, work_data_ptr));

   }

}

Ncep_Ncar_Climate::Work_Data_Map::~Work_Data_Map ()
{

   for (Work_Data_Map::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Work_Data* work_data_ptr = iterator->second;
      delete work_data_ptr;
   }

}

void
Ncep_Ncar_Climate::Work_Data_Map::calculate (Ncep_Ncar& ncep_ncar,
                                             const Integer tolerance)
{
   sum_up (ncep_ncar, tolerance);
   divide ();
}

void
Ncep_Ncar_Climate::Work_Data_Map::write (const string& path) const
{

   for (Ituple::const_iterator iterator = ituple_jjj.begin ();
        iterator != ituple_jjj.end (); iterator++)
   {
      const Integer jjj = *(iterator);
      const Work_Data& work_data = *(at (jjj));
      work_data.write (path, jjj);
   }

}

Integer
Ncep_Ncar_Climate::get_nearest_jjj (const Integer jjj) const
{

   Integer min_jjj;
   Integer min_delta = 9999;

   for (Ituple::const_iterator iterator = ituple_jjj.begin ();
        iterator != ituple_jjj.end (); iterator++)
   {
      const Integer j = *(iterator);
      const Integer delta = abs (j - jjj);
      if (delta < min_delta)
      {
         min_delta = delta;
         min_jjj = j;
      }
   } 

   return min_jjj;

}

string
Ncep_Ncar_Climate::get_file_path (const Integer jjj) const
{
   const string& file_name = string_render ("nnc_%03d", jjj);
   const string& file_path = path + "/" + file_name;
   return file_path;
}

void
Ncep_Ncar_Climate::read_3d (Nwp::Data_3D& data_3d,
                            const Integer jjj) const
{

   const Nwp_Element nwp_elements[] = { TEMPERATURE, DEW_POINT,
      GEOPOTENTIAL_HEIGHT, ZONAL_WIND, MERIDIONAL_WIND, OMEGA,
      VARIANCE_TEMPERATURE, VARIANCE_DEW_POINT,
      VARIANCE_GEOPOTENTIAL_HEIGHT, VARIANCE_ZONAL_WIND,
      VARIANCE_MERIDIONAL_WIND, VARIANCE_OMEGA };

   typedef Geodetic_Vector_Data_3D Gvd_3d;
   const Tuple& tuple_latitude = Ncep_Ncar::tuple_latitude ();
   const Tuple& tuple_longitude = Ncep_Ncar::tuple_longitude ();

   for (Integer i = 0; i < 12; i++)
   {
      const Nwp_Element& nwp_element = nwp_elements[i];
      Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
         tuple_latitude, tuple_longitude, true);
      data_3d.set_gvd_3d_ptr (nwp_element, gvd_3d_ptr);
   }

   const Integer ni = tuple_latitude.size ();
   const Integer nj = tuple_longitude.size ();
   const Integer offset_3d = 12 * ni * nj * sizeof (float);

   ifstream file (get_file_path (jjj).c_str ());
   file.seekg (offset_3d);
   data_3d.read (file);
   file.close ();

}

void
Ncep_Ncar_Climate::write (const Work_Data_Map& work_data_map) const
{

   const string file_path = path + "/nnc_header";
   ofstream header (file_path.c_str ());
   header << tuple_p << endl;
   header << tolerance << endl;
   header.close ();

   work_data_map.write (path);

}

Geodetic_Vector_Data_2D*
Ncep_Ncar_Climate::get_initialized_vd_2d (const Integer vector_size) const
{
   typedef Geodetic_Vector_Data_2D Gvd_2d;
   const Tuple tuple_latitude = Ncep_Ncar::tuple_latitude ();
   const Tuple tuple_longitude = Ncep_Ncar::tuple_longitude ();
   return new Gvd_2d (vector_size, tuple_latitude, tuple_longitude, true);
}

void
Ncep_Ncar_Climate::initialize_3d_data (const Key& key)
{

   typedef Ncep_Ncar_Climate::Data_3D Nncd_3d;

   Nncd_3d* nncd_3d_ptr = new Nncd_3d (nwp_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, nncd_3d_ptr));

}

void
Ncep_Ncar_Climate::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef Ncep_Ncar_Climate::Data_3D Nncd_3d;
   Nncd_3d& nncd_3d = dynamic_cast<Nncd_3d&>(data_3d);

   const Dtime& dtime = data_3d.key.base_time;
   const Integer jjj = dtime.get_day_of_year ();

   read_3d (data_3d, jjj);
   data_3d.set_available ();

}

void
Ncep_Ncar_Climate::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                           const Integer vector_index,
                                           const Key& key,
                                           const Nwp_Element nwp_element)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Integer plate_size = size_2d.i * size_2d.j * sizeof (float);

   const Dtime& dtime = key.base_time;
   const Integer jjj = dtime.get_day_of_year ();

   switch (nwp_element)
   {

      case TEMPERATURE:
      {
         ifstream file (get_file_path (jjj).c_str ());
         //gvd_2d.read_chunk (file, vector_index);
         file.close ();
         return;
      }

      case DEW_POINT:
      {
         ifstream file (get_file_path (jjj).c_str ());
         file.seekg (plate_size);
         //gvd_2d.read_chunk (file, vector_index);
         file.close ();
         return;
      }

      case RELATIVE_HUMIDITY:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);

         fill_screen_level_data (*data_ptr, 0, key, TEMPERATURE);
         fill_screen_level_data (*data_ptr, 0, key, DEW_POINT);

         const Size_2D& size_2d = data_ptr->get_size_2d ();

         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.i; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real td = data_ptr->get_datum (1, i, j);
               const Real rh = Moisture::get_rh (t - K, td - K);
               gvd_2d.set_datum (vector_index, i, j, rh);
            }
         }

         delete data_ptr;
         return;

      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ncep_Ncar_Climate::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                        const Integer vector_index,
                                        const Key& key,
                                        const Nwp_Element nwp_element)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Integer plate_size = size_2d.i * size_2d.j * sizeof (float);

   const Dtime& dtime = key.base_time;
   const Integer jjj = dtime.get_day_of_year ();

   switch (nwp_element)
   {

      case ZONAL_WIND:
      {
         ifstream file (get_file_path (jjj).c_str ());
         file.seekg (2 * plate_size);
         //gvd_2d.read_chunk (file, vector_index);
         file.close ();
         return;
      }

      case MERIDIONAL_WIND:
      {
         ifstream file (get_file_path (jjj).c_str ());
         file.seekg (3 * plate_size);
         //gvd_2d.read_chunk (file, vector_index);
         file.close ();
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ncep_Ncar_Climate::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                                  const Integer vector_index,
                                  const Key& key,
                                  const Nwp_Element nwp_element)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Integer plate_size = size_2d.i * size_2d.j * sizeof (float);

   const Dtime& dtime = key.base_time;
   const Integer jjj = dtime.get_day_of_year ();

   switch (nwp_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         ifstream file (get_file_path (jjj).c_str ());
         file.seekg (4 * plate_size);
         //gvd_2d.read_chunk (file, vector_index);
         file.close ();
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, nwp_element);

}

void
Ncep_Ncar_Climate::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                            const Integer vector_index,
                                            const Key& key,
                                            const Nwp_Element nwp_element)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Integer plate_size = size_2d.i * size_2d.j * sizeof (float);

   const Dtime& dtime = key.base_time;
   const Integer jjj = dtime.get_day_of_year ();

   switch (nwp_element)
   {

      case PRESSURE:
      {
         ifstream file (get_file_path (jjj).c_str ());
         file.seekg (5 * plate_size);
         //gvd_2d.read_chunk (file, vector_index);
         file.close ();
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, nwp_element);

}

Ncep_Ncar_Climate::Ncep_Ncar_Climate (const string& description,
                                      const string& path,
                                      const Integer nominal_year)
   : Nwp (description, path)
{

   nwp_element_vector.push_back (denise::TEMPERATURE);
   nwp_element_vector.push_back (denise::DEW_POINT);
   nwp_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::ZONAL_WIND);
   nwp_element_vector.push_back (denise::MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::OMEGA);
   nwp_element_vector.push_back (denise::VARIANCE_TEMPERATURE);
   nwp_element_vector.push_back (denise::VARIANCE_DEW_POINT);
   nwp_element_vector.push_back (denise::VARIANCE_GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::VARIANCE_ZONAL_WIND);
   nwp_element_vector.push_back (denise::VARIANCE_MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::VARIANCE_OMEGA);

   set_nominal_year (nominal_year);
   survey ();

   string input_line;
   const string header_file_path = path + "/" + "nnc_header";
   ifstream header_file (header_file_path.c_str ());
   std::getline (header_file, input_line);
   header_file.close ();
   this->tuple_p = Tuple (input_line);

}

Ncep_Ncar_Climate::Ncep_Ncar_Climate (Ncep_Ncar& ncep_ncar,
                                      const string& path,
                                      const Ituple ituple_jjj,
                                      const Integer tolerance,
                                      const Tuple& tuple_p,
                                      const Integer nominal_year)
   : Nwp ("NN_Climate", path),
     ituple_jjj (ituple_jjj),
     tolerance (tolerance)
{

   this->tuple_p = tuple_p;
   set_nominal_year (nominal_year);

   nwp_element_vector.push_back (denise::TEMPERATURE);
   nwp_element_vector.push_back (denise::DEW_POINT);
   nwp_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::ZONAL_WIND);
   nwp_element_vector.push_back (denise::MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::OMEGA);
   nwp_element_vector.push_back (denise::VARIANCE_TEMPERATURE);
   nwp_element_vector.push_back (denise::VARIANCE_DEW_POINT);
   nwp_element_vector.push_back (denise::VARIANCE_GEOPOTENTIAL_HEIGHT);
   nwp_element_vector.push_back (denise::VARIANCE_ZONAL_WIND);
   nwp_element_vector.push_back (denise::VARIANCE_MERIDIONAL_WIND);
   nwp_element_vector.push_back (denise::VARIANCE_OMEGA);

   map<Integer, Integer> n_map_ptr;
   map<Integer, Geodetic_Vector_Data_3D*> climate_map;

   const Tuple& tuple_latitude = Ncep_Ncar::tuple_latitude ();
   const Tuple& tuple_longitude = Ncep_Ncar::tuple_longitude ();

   Work_Data_Map work_data_map (ituple_jjj, tuple_p,
      tuple_latitude, tuple_longitude);
   work_data_map.calculate (ncep_ncar, tolerance);
   write (work_data_map);

};

Ncep_Ncar_Climate::~Ncep_Ncar_Climate ()
{
}

void
Ncep_Ncar_Climate::set_nominal_year (const Integer nominal_year)
{
   if (nominal_year > -998) { this->nominal_year = nominal_year; }
   else { this->nominal_year = Dtime ().get_year (); }
}

void
Ncep_Ncar_Climate::survey ()
{

   string tuple_p_str, tolerance_str;

   const string& file_path = path + "/nnc_header";
   ifstream header (file_path.c_str ());
   header >> tuple_p_str;
   header >> tolerance_str;
   header.close ();

   tuple_p = Tuple (tuple_p_str);
   tolerance = Integer (atoi (tolerance_str.c_str ()));

   const string str ("nnc_[0-9]..");
   const vector<string>& dir_listing = get_dir_listing (path, str);

   ituple_jjj.clear ();

   for (vector<string>::const_iterator iterator = dir_listing.begin ();
        iterator != dir_listing.end (); iterator++)
   {
      const string& file_name = *(iterator);
      const string& file_path = path + "/" + file_name;
      const Integer jjj = Integer (atoi (file_name.substr (4, 3).c_str ()));
      ituple_jjj.push_back (jjj);
   }

   std::sort (ituple_jjj.begin (), ituple_jjj.end ());


   set<Dtime> valid_time_set = get_valid_time_set ();

   for (set<Dtime>::const_iterator iterator = valid_time_set.begin ();
        iterator != valid_time_set.end (); iterator++)
   {
      const Dtime& dtime = *(iterator);
      const Key key (dtime, 0);
      key_multimap.add (key);
   }

   status = "Loaded";

   for (Key_Multimap::const_iterator iterator = key_multimap.begin ();
        iterator != key_multimap.end (); iterator++)
   {
      const Key key (iterator->first, iterator->second);
      initialize_3d_data (key);
   }


}

void
Ncep_Ncar_Climate::clean_up ()
{
   clear_data_3d_ptr_map ();
}

set<Dtime>
Ncep_Ncar_Climate::get_valid_time_set () const
{

   set<Dtime> valid_time_set;
   const Integer yyyy = nominal_year;

   for (Ituple::const_iterator iterator = ituple_jjj.begin ();
        iterator != ituple_jjj.end (); iterator++)
   {
      const Integer jjj = *(iterator);
      const string& time_str = string_render ("%04d%03d", yyyy, jjj);
      const Dtime dtime (time_str, string ("%Y%j"));
      valid_time_set.insert (dtime);
   }

   return valid_time_set;

}
         
Nwp::Key
Ncep_Ncar_Climate::get_key (const Dtime& dtime) const
{
   const Integer yyyy = nominal_year;
   const Integer jjj = get_nearest_jjj (dtime.get_day_of_year ());
   const string fmt ("%Y%j");
   return Key (Dtime (string_render ("%04d%03d", yyyy, jjj), fmt), 0);
}

void
Ncep_Ncar_Climate::acquire_base_time_forecast_hour (Dtime& base_time,
                                                    Integer& forecast_hour,
                                                    const Dtime& dtime) const
{
   const Integer yyyy = nominal_year;
   const Integer jjj = get_nearest_jjj (dtime.get_day_of_year ());
   const string fmt ("%Y%j");
   base_time = Dtime (string_render ("%04d%03d", yyyy, jjj), fmt);
   forecast_hour = 0;
}

