//
// nwp.cc
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
#include "nwp.h"

using namespace std;
using namespace denise;

//Level_Tuple::Level_Vector ()
//{
//}

Level_Tuple::Level_Tuple (const Level::Type type,
                          const Dstring& str,
                          const Dstring& delimiter)
   : Tuple (str, delimiter),
     type (type)
{
   add (str, delimiter);
}

void
Level_Tuple::add (const Dstring& str,
                  const Dstring& delimiter,
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

      case Level::PRESSURE:
      {
         for (I i = rbegin (); i != rend (); i++) { if (*i < v) { return *i; } }
         break;
      }

      case Level::THETA:
      case Level::SIGMA:
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

      case Level::PRESSURE:
      {
         for (I i = begin (); i != end (); i++) { if (*i > v) { return *i; } }
         break;
      }

      case Level::SIGMA:
      case Level::THETA:
      {
         for (I i = begin (); i != end (); i++) { if (*i < v) { return *i; } }
         break;
      }

   }

   throw Nwp_Exception ("Level_Vector::get_next_up confused");

}

Level_Element::Level_Element (const Level& level,
                              const Met_Element met_element)
   : level (level),
     met_element (met_element)
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

Nwp::Element_Vector::Element_Vector (const vector<Met_Element>& met_element_vector)
   : vector<Met_Element> (met_element_vector)
{

   for (Integer i = 0; i < met_element_vector.size (); i++)
   {
      const Met_Element& met_element = met_element_vector[i];
      reverse_map.insert (make_pair (met_element, i));
   }

}

Integer
Nwp::Element_Vector::get_index (const Met_Element met_element) const
{

   return reverse_map.at (met_element);

   for (Integer i = 0; i < size (); i++)
   {
      if (at (i) == met_element) { return i; }
   }

   throw Nwp_Exception ("met_element not available.");

}

Met_Element
Nwp::Element_Vector::get_met_element (const Integer element_index) const
{
   return at (element_index);
}

Nwp::Data_3D::Data_3D (const vector<Met_Element>& met_element_vector,
                       const Key& key)
  : key (key),
    available (false),
    met_element_vector (met_element_vector)
{

   typedef vector<Met_Element>::const_iterator Iterator;

   for (Iterator iterator = met_element_vector.begin ();
        iterator != met_element_vector.end (); iterator++)
   {
      const Met_Element& met_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = NULL;
      insert (make_pair (met_element, gvd_3d_ptr));
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
Nwp::Data_3D::read (ifstream& file,
                    const bool float_length)
{

   typedef vector<Met_Element>::const_iterator Iterator;
   for (Iterator iterator = met_element_vector.begin ();
        iterator != met_element_vector.end (); iterator++)
   {
      const Met_Element met_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = at (met_element);
      if (gvd_3d_ptr == NULL) { continue; }
      gvd_3d_ptr->read (file, float_length);
   }
}

void
Nwp::Data_3D::write (ofstream& file,
                     const bool float_length) const
{
   typedef vector<Met_Element>::const_iterator Iterator;
   for (Iterator iterator = met_element_vector.begin ();
        iterator != met_element_vector.end (); iterator++)
   {
      const Met_Element met_element = *(iterator);
      const Geodetic_Vector_Data_3D* gvd_3d_ptr = at (met_element);
      if (gvd_3d_ptr == NULL) { continue; }
      gvd_3d_ptr->write (file, float_length);
   }
}

void
Nwp::Data_3D::unload ()
{
   for (Nwp::Data_3D::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Geodetic_Vector_Data_3D* gvd_3d_ptr = iterator->second;
      if (gvd_3d_ptr != NULL) { delete gvd_3d_ptr; iterator->second = NULL; }
   }
}

void
Nwp::Data_3D::unload (const Met_Element met_element)
{
   Geodetic_Vector_Data_3D* gvd_3d_ptr = at (met_element);
   if (gvd_3d_ptr != NULL) { delete gvd_3d_ptr; at (met_element) = NULL; }
}

void
Nwp::Data_3D::set_gvd_3d_ptr (const Met_Element met_element,
                              Geodetic_Vector_Data_3D* gvd_3d_ptr)
{
   at (met_element) = gvd_3d_ptr;
}

const Geodetic_Vector_Data_3D&
Nwp::Data_3D::get_gvd_3d (const Met_Element met_element) const
{
   try
   {
      const Geodetic_Vector_Data_3D* gvd_3d_ptr = at (met_element);
      if (gvd_3d_ptr == NULL)
      {
         throw Nwp_Exception ("Nwp::Data_3D gvd_3d_ptr is NUL");
      }
      return *gvd_3d_ptr;
   }
   catch (const std::exception& se)
   {
      throw Nwp_Exception ("Nwp::Data_3D has no such met_element");
   }
}

Geodetic_Vector_Data_3D&
Nwp::Data_3D::get_gvd_3d (const Met_Element met_element)
{
   Geodetic_Vector_Data_3D* gvd_3d_ptr = at (met_element);
   try
   {
      Geodetic_Vector_Data_3D* gvd_3d_ptr = at (met_element);
      if (gvd_3d_ptr == NULL)
      {
         throw Nwp_Exception ("Nwp::Data_3D gvd_3d_ptr is NUL");
      }
      return *gvd_3d_ptr;
   }
   catch (const std::exception& se)
   {
      throw Nwp_Exception ("Nwp::Data_3D has no such met_element");
   }
}

const Tuple&
Nwp::Data_3D::get_tuple_p (const Met_Element met_element) const
{
   try
   {
      const Geodetic_Vector_Data_3D& gvd_3d = get_gvd_3d (met_element);
      return gvd_3d.get_coordinate_tuple (0);
   }
   catch (const Nwp_Exception& ne)
   {
      const Geodetic_Vector_Data_3D& gvd_3d = *(begin ()->second);
      return gvd_3d.get_coordinate_tuple (0);
   }
}

Real
Nwp::Data_3D::get_p (const Met_Element met_element,
                     const Integer k) const
{
   const Tuple& tuple_p = get_tuple_p (met_element);
   return tuple_p[k];
}

Lat_Long
Nwp::Data_3D::get_lat_long (const Met_Element met_element,
                            const Integer i,
                            const Integer j) const
{
   const Geodetic_Vector_Data_3D& gvd_3d = get_gvd_3d (met_element);
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
Nwp::Data_3D::get_p_from_element (const Met_Element met_element,
                                  const Real latitude,
                                  const Real longitude,
                                  const Real element_value) const
{

   const Geodetic_Vector_Data_3D* gvd_3d_ptr = at (met_element);
   if (gvd_3d_ptr == NULL) { return GSL_NAN; }
   const Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   Real p = GSL_NAN;
   const Real x = element_value;
   const Tuple& tuple_p = get_tuple_p (met_element);

   for (Integer k = 0; k < tuple_p.size () - 1; k++)
   {

      const Real lower_p = tuple_p[k];
      const Real upper_p = tuple_p[k + 1];
      const Real lower_x = evaluate (met_element, lower_p, latitude, longitude);
      const Real upper_x = evaluate (met_element, upper_p, latitude, longitude);
      const Real delta_x = x - lower_x;
      const bool match = (delta_x * (x - upper_x) <= 0);

      if (match)
      {
         const Real dp = upper_p - lower_p;
         const Real dx = upper_x - lower_x;
         p = lower_p + (delta_x  * dp / dx);
         break;
      }

   }

   return p;

}

Real
Nwp::Data_3D::evaluate (const Met_Element element,
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
         const Met_Element T = TEMPERATURE;
         const Met_Element Z = GEOPOTENTIAL_HEIGHT;
         const Real t = evaluate (T, p, latitude, longitude);
         const Real z = evaluate (Z, p, latitude, longitude);
         return g * z + c_p * t;
      }

      case ABSOLUTE_VORTICITY:
      {

         const Real f = Geodetic_Vector_Data_2D::get_f (latitude);

         const Met_Element U = ZONAL_WIND;
         const Met_Element V = MERIDIONAL_WIND;

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

         const Met_Element T = TEMPERATURE;
         const Met_Element U = ZONAL_WIND;
         const Met_Element V = MERIDIONAL_WIND;
         const Met_Element W = VERTICAL_VELOCITY;
         const Met_Element Z = GEOPOTENTIAL_HEIGHT;

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
Nwp::Data_3D::evaluate (const Met_Element element,
                        const Real p,
                        const Lat_Long& lat_long,
                        const Evaluate_Op evaluate_op) const
{
   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;
   return evaluate (element, p, latitude, longitude, evaluate_op);
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
Nwp::Cross_Section::insert_met_element_if_needed (const Met_Element met_element,
                                                  const Tuple& tuple_x,
                                                  const Tuple& tuple_p)
{
   if (find (met_element) == end ())
   {
      Scalar_Data_2D* sd_2d_ptr = new Scalar_Data_2D (tuple_x, tuple_p);
      insert (make_pair (met_element, sd_2d_ptr));
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
Nwp::Cross_Section::get_sd_2d (const Met_Element met_element) const
{
   return *at (met_element);
}

Scalar_Data_2D&
Nwp::Cross_Section::get_sd_2d (const Met_Element met_element)
{
   return *at (met_element);
}

const Tuple&
Nwp::Cross_Section::get_tuple_p (const Met_Element met_element) const
{
   return at (met_element)->get_coordinate_tuple (1);
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
Nwp::Time_Cross::insert_met_element_if_needed (const Met_Element met_element,
                                               const Tuple& tuple_t,
                                               const Tuple& tuple_p)
{
   if (find (met_element) == end ())
   {
      Scalar_Data_2D* sd_2d_ptr = new Scalar_Data_2D (tuple_t, tuple_p);
      insert (make_pair (met_element, sd_2d_ptr));
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
Nwp::Time_Cross::get_sd_2d (const Met_Element met_element) const
{
   return *at (met_element);
}

Scalar_Data_2D&
Nwp::Time_Cross::get_sd_2d (const Met_Element met_element)
{
   return *at (met_element);
}

const Tuple&
Nwp::Time_Cross::get_tuple_p (const Met_Element met_element) const
{
   return at (met_element)->get_coordinate_tuple (1);
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
   const Met_Element T = TEMPERATURE;

   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_coordinate (0, i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_coordinate (1, j);
         Real next_p = next_p_data_ptr->get_datum (0, latitude, longitude);
         Real next_t = data_3d.evaluate (T, next_p, latitude, longitude);
         Real this_t = this_t_data_ptr->get_datum (0, latitude, longitude);
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
                     const Met_Element met_element)
{

   switch (met_element)
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
   const Met_Element PV = POTENTIAL_VORTICITY;

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      const Real latitude = gvd_2d.get_coordinate (0, i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_coordinate (1, j);
         const Real pv = data_3d.evaluate (PV, p, latitude, longitude);
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

   vector<Met_Element> ev;
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

      case Level::PRESSURE:
      {
         const Real p = level.value;
         fill_p_potential_vorticity_data (gvd_2d, vi, key, p);
         return;
      }

      case Level::THETA:
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

   const Met_Element TA = TEMPERATURE_ADVECTION;
   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   Nwp::Data_3D& data_3d = get_3d_data (key);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      const Real latitude = gvd_2d.get_coordinate (0, i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_coordinate (1, j);
         const Real ta = data_3d.evaluate (TA, p, latitude, longitude);
         gvd_2d.set_datum (vector_index, i, j, ta);
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

   vector<Met_Element> ev;
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

      case Level::PRESSURE:
      {
         const Real p = level.value;
         fill_p_temperature_advection_data (gvd_2d, vi, key, p);
         return;
      }

      case Level::THETA:
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

   const Met_Element AH = ADIABATIC_HEATING;
   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   Nwp::Data_3D& data_3d = get_3d_data (key);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      const Real latitude = gvd_2d.get_coordinate (0, i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_coordinate (1, j);
         const Real pv = data_3d.evaluate (AH, p, latitude, longitude);
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

   vector<Met_Element> ev;
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

      case Level::PRESSURE:
      {
         const Real p = level.value;
         fill_p_adiabatic_heating_data (gvd_2d, vi, key, p);
         return;
      }

      case Level::THETA:
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

   const Met_Element LH = LATENT_HEATING;
   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   Nwp::Data_3D& data_3d = get_3d_data (key);

   for (Integer i = 0; i < size_2d.i; i++)
   { 
      const Real latitude = gvd_2d.get_coordinate (0, i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_coordinate (1, j);
         const Real pv = data_3d.evaluate (LH, p, latitude, longitude);
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

   vector<Met_Element> ev;
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

         case Level::PRESSURE:
         {
            const Real p = level.value;
            fill_p_latent_heating_data (gvd_2d, vi, key, p);
            return;
         }

         case Level::THETA:
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


   const Met_Element PV = POTENTIAL_VORTICITY;

   const Data_3D& data_3d = get_3d_data (key);
   const Tuple& tuple_p = data_3d.get_tuple_p (TEMPERATURE);
   const Integer nk = tuple_p.size ();

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
            Real lower_pv = data_3d.evaluate (PV, lower_p, latitude, longitude);
            Real upper_pv = data_3d.evaluate (PV, upper_p, latitude, longitude);

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
   const Real p_850 = 850e2;
   const Met_Element T = TEMPERATURE;
   const Met_Element Z = GEOPOTENTIAL_HEIGHT;

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_coordinate (0, i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_coordinate (1, j);
         const Real t = data_3d.evaluate (T, p_850, latitude, longitude) - K;
         const Real z = data_3d.evaluate (Z, p_850, latitude, longitude);
         const Real thick = data_ptr->get_datum (0, i, j);
         const Real surface_p = data_ptr->get_datum (1, i, j);
         const Real e = data_3d.evaluate (Z, surface_p, latitude, longitude);
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
                    const Met_Element met_element)
{

   if (met_element != FFDI && met_element != GFDI)
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

         switch (met_element)
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

   const Met_Element T = TEMPERATURE;
   const Met_Element TD = DEW_POINT;

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

   const Met_Element T = TEMPERATURE;
   const Met_Element TD = DEW_POINT;

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

   const Met_Element T = TEMPERATURE;
   const Met_Element TD = DEW_POINT;

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

   const Met_Element T = TEMPERATURE;
   const Met_Element TD = DEW_POINT;

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

            case Level::PRESSURE:
               start_p = level.value;
               break;

            case Level::MEAN_SEA:
            case Level::TEN_METRE:
            case Level::FIFTY_METRE:
            case Level::SCREEN:
               start_p = surface_p_data_ptr->evaluate (0, i, j);
               break;
 
            default:
               Dstring error_str = "Nwp::fill_thunder_li_data ";
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
                             const Met_Element met_element)
{

   const Integer vi = vector_index;

   switch (met_element)
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
                      const Met_Element met_element)
{
   throw Nwp_Exception ("Nwp::fill_cloud_data N/A");
}

void
Nwp::fill_pressure_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Real p,
                               const Met_Element met_element)
{

   switch (met_element)
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

   //if ((p - tuple_p.front ()) * (p - tuple_p.back ()) > 0)
   //{
   //   Dstring error_str = "Nwp::fill_pressure_level_data_ptr ";
   //   error_str += string_render (" %f", p);
   //   throw Nwp_Exception (error_str);
   //}

   Nwp::Data_3D& data_3d = get_3d_data (key);
   const Met_Element element = met_element;
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
                            const Met_Element met_element)
{

   const Integer vi = vector_index;

   switch (met_element)
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
                            const Met_Element met_element,
                            const Geodetic_Vector_Data_2D& surface_p_data)
{

   const Data_3D& data_3d = get_3d_data (key);
   const Tuple& tuple_p = data_3d.get_tuple_p (TEMPERATURE);
   const Real& start_p = tuple_p.front ();
   const Real& end_p = tuple_p.back ();

   const Met_Element& ne = met_element;
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
                             const Met_Element met_element)
{

   switch (met_element)
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
                          const Met_Element met_element)
{  
   throw Nwp_Exception ("Nwp::fill_50m_level_data N/A");
}

void
Nwp::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Met_Element met_element)
{  
   throw Nwp_Exception ("Nwp::fill_10m_level_data N/A");
}

void
Nwp::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                    const Integer vector_index,
                    const Key& key,
                    const Met_Element met_element)
{
   throw Nwp_Exception ("Nwp::fill_msl_data N/A");
}

void
Nwp::fill_nil_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Met_Element met_element)
{

   const Integer vi = vector_index;

   switch (met_element)
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
         fill_ts_diagnosis_data (gvd_2d, vi, key, level, met_element);
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
         fill_cloud_data (gvd_2d, vi, key, met_element);
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
                              const Met_Element met_element)
{
   throw Nwp_Exception ("Nwp::fill_surface_data N/A");
}

void
Nwp::fill_data (Geodetic_Vector_Data_2D& gvd_2d,
                const Integer vector_index,
                const Key& key,
                const Level& level,
                const Met_Element met_element)
{

   switch (met_element)
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
         fill_rain_data (gvd_2d, vector_index, key, met_element);
         return;
      };

      case GFDI:
      case FFDI:
      {
         fill_fdi_data (gvd_2d, vector_index, key, level, met_element);
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
   
      case Level::PRESSURE:
      {
         const Real v = level.value;
         const Integer vi = vector_index;
         fill_pressure_level_data (gvd_2d, vi, key, v, met_element);
         return;
      }

      case Level::THETA:
      {
         const Real v = level.value;
         const Integer vi = vector_index;
         fill_theta_level_data (gvd_2d, vi, key, v, met_element);
         return;
      }

      case Level::SIGMA:
      {
         const Real v = level.value;
         const Integer vi = vector_index;
         const Level& surface = Level::surface_level ();
         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (1);
         fill_data (*data_ptr, 0, key, surface, PRESSURE);
         fill_sigma_level_data (gvd_2d, vi, key, v, met_element, *data_ptr);
         delete data_ptr;
         return;
      }

      case Level::SCREEN:
      {
         fill_screen_level_data (gvd_2d, vector_index, key, met_element);
         return;
      }

      case Level::FIFTY_METRE:
      {
         fill_50m_level_data (gvd_2d, vector_index, key, met_element);
         return;
      }

      case Level::TEN_METRE:
      {
         fill_10m_level_data (gvd_2d, vector_index, key, met_element);
         return;
      }

      case Level::MEAN_SEA:
      {
         fill_msl_data (gvd_2d, vector_index, key, met_element);
         return;
      }

      case Level::NIL:
      {
         fill_nil_level_data (gvd_2d, vector_index, key, met_element);
         return;
      }

      case Level::SURFACE:
      {
         fill_surface_level_data (gvd_2d, vector_index, key, met_element);
         return;
      }

   }

   throw Nwp_Exception ("Nwp::fill_data N/A");

}

Nwp::Nwp (const Dstring& description,
          const Dstring& path)
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

const Dstring&
Nwp::get_description () const
{
   return description;
}

const Dstring&
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
      Dstring error_str = "Nwp::get_3d_data get_3d_data failed  ";
      const Dtime& base_time = key.base_time;
      const Integer forecast_hour = key.forecast_hour;
      error_str += base_time.get_string ();
      error_str += Dstring::render (" +%d hr", forecast_hour);
      throw Nwp_Exception (error_str);
   }

   Data_3D& data_3d = *(iterator->second);
   if (!data_3d.is_available ()) { load_3d_data (data_3d); }
   return data_3d;

}

Geodetic_Vector_Data_2D*
Nwp::get_ts_steering_data_ptr (const Key& key)
{
   const Level level (Level::PRESSURE, 600e2, 800e2);
   return get_steering_data_ptr (key, level);
}

Geodetic_Vector_Data_2D*
Nwp::get_steering_data_ptr (const Key& key,
                            const Level& level)
{

   if (level.type != Level::PRESSURE)
   {
      throw Nwp_Exception ("Can't do non-P Level");
   }

   const Met_Element& U = ZONAL_WIND;
   const Met_Element& V = MERIDIONAL_WIND;
   const Data_3D& data_3d = get_3d_data (key);

   Integer n = 0;
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (2);
   data_ptr->initialize (0, 0);
   data_ptr->initialize (1, 0);

   const Tuple& tuple_p = data_3d.get_tuple_p (ZONAL_WIND);
   const Integer nk = tuple_p.size ();
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer k = 0; k < nk; k++)
   {

      // Iterate each level
      // If this level is out of boundaries p_a, p_b, then continue
      const Real p = tuple_p[k];
      if ((p - level.value) * (p - level.value_) > 0) { continue; }
      const Level& level = Level::pressure_level (p);

      try
      {

         for (Integer i = 0; i < size_2d.i; i++)
         {
            const Real latitude = data_ptr->get_coordinate (0, i);
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real longitude = data_ptr->get_coordinate (1, j);

               const Real u = data_3d.evaluate (U, p, latitude, longitude);
               const Real v = data_3d.evaluate (V, p, latitude, longitude);

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

   const Real reciprocal = Real (1) / n;
   data_ptr->scale_offset (0, reciprocal, 0);
   data_ptr->scale_offset (1, reciprocal, 0);

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_vertical_shear_data_ptr (const Key& key,
                                  const Level& level)
{

   const Level& level_a = level.get_level_0 ();
   const Level& level_b = level.get_level_1 ();

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
Nwp::get_min_li_thunder_data_ptr (const Key& key,
                                  const Real thunder_t)
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   const Met_Element T = TEMPERATURE;
   const Met_Element TD = DEW_POINT;

   const Data_3D& data_3d = get_3d_data (key);
   const Tuple& tuple_p = data_3d.get_tuple_p (T);

   Gvd_2d* data_ptr = get_initialized_vd_2d (2);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   Gvd_2d* p_data_ptr = get_temperature_p_data_ptr (key, thunder_t);

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = data_ptr->get_coordinate (0, i);

      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = data_ptr->get_coordinate (1, j);

         Real min_li_thunder_p = GSL_NAN;
         Real min_li_thunder = GSL_POSINF;

         for (Integer k = 0; k < tuple_p.size (); k++)
         {

            const Real p = tuple_p[k];
            const Real end_p = p_data_ptr->evaluate (0, i, j);
            if (p < end_p || gsl_isnan (end_p)) { continue; }

            const Real start_t = data_3d.evaluate (T, p, latitude, longitude);
            const Real start_td = data_3d.evaluate (TD, p, latitude, longitude);
            const Real li_thunder = Instability::get_lifted_index (
               p, start_t, start_td, end_p, thunder_t - K) + K;

            if (gsl_finite (li_thunder))
            {
               if (li_thunder < min_li_thunder)
               {
                  min_li_thunder = li_thunder;
                  min_li_thunder_p = p;
               }
            }

         }

         data_ptr->set_datum (0, i, j, min_li_thunder);
         data_ptr->set_datum (1, i, j, min_li_thunder_p);

      }
   }

   delete p_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_freezing_level_data_ptr (const Key& key)

{

   const Met_Element T = TEMPERATURE;
   const Met_Element Z = GEOPOTENTIAL_HEIGHT;

   const Data_3D& data_3d = get_3d_data (key);
   const Tuple& tuple_p = data_3d.get_tuple_p (T);
   const Integer nk = tuple_p.size ();

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {

      Real* array = new Real[3];
      const Real latitude = data_ptr->get_coordinate (0, i);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         Integer filled = 0;
         array[0] = GSL_NAN;
         array[1] = GSL_NAN;
         array[2] = GSL_NAN;
         const Real longitude = data_ptr->get_coordinate (1, j);

         for (Integer k = 0; k < nk - 1; k++)
         {

            const Real lower_p = tuple_p[k];
            const Real upper_p = tuple_p[k + 1];
            Real lt = data_3d.evaluate (T, lower_p, latitude, longitude) - K; 
            Real ut = data_3d.evaluate (T, upper_p, latitude, longitude) - K;

            if (lt * ut <= 0)
            {

               Real lz = data_3d.evaluate (Z, lower_p, latitude, longitude);
               Real uz = data_3d.evaluate (Z, upper_p, latitude, longitude);
               const Real dz = uz - lz;
               const Real dt = ut - lt;

               const Real z = lz - (lt * dz / dt);

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

      const Real latitude = data_ptr->get_coordinate (0, i);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         Real p = GSL_NAN;
         const Real longitude = data_ptr->get_coordinate (1, j);

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

      case Level::PRESSURE:
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

      case Level::THETA:
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

      //case Level::SIGMA:
      //{
      //   break;
      //}

      case Level::SCREEN:
      case Level::FIFTY_METRE:
      case Level::TEN_METRE:
      case Level::MEAN_SEA:
      case Level::SURFACE:
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

   if (level.type != Level::PRESSURE) { throw Exception ("Only for P Levels"); }

   const Met_Element T = TEMPERATURE;
   const Met_Element Z = GEOPOTENTIAL_HEIGHT;
   const Met_Element U = ZONAL_WIND;
   const Met_Element V = MERIDIONAL_WIND;

   const Data_3D& data_3d = get_3d_data (key);

   const Tuple& tuple_p = data_3d.get_tuple_p (Z);
   const Integer nk = tuple_p.size ();

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (3);
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
         data_ptr->set_datum (1, i, j, GSL_NAN);
         data_ptr->set_datum (2, i, j, GSL_NAN);

         for (Integer k = nk - 1; k >= 0; k--)
         {

            const Real p = tuple_p[k];
            if (p > surface_p || p > level.value) { continue; }

            const Real rh = data_3d.evaluate (
               RELATIVE_HUMIDITY, p, latitude, longitude);

            if ((rh >= 0.85) ||
                (p < 850e2 && rh >= 0.8) ||
                (p < 600e2 && rh >= 0.75))
            {
               const Real datum = data_3d.evaluate (Z, p, latitude, longitude);
               const Real u = data_3d.evaluate (U, p, latitude, longitude);
               const Real v = data_3d.evaluate (V, p, latitude, longitude);
               data_ptr->set_datum (0, i, j, datum);
               data_ptr->set_datum (1, i, j, u);
               data_ptr->set_datum (2, i, j, v);
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

      case Level::PRESSURE:
      {
         const Real p = level.value;
         return get_ageostrophic_wind_data_ptr (key, p);
      }

      case Level::TEN_METRE:
      case Level::FIFTY_METRE:
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

      case Level::PRESSURE:
      {
         const Real p = level.value;
         return get_geostrophic_wind_data_ptr (key, p);
      }

      case Level::FIFTY_METRE:
      case Level::TEN_METRE:
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

      case Level::PRESSURE:
      {
         vector<Met_Element> ev;
         ev.push_back (ZONAL_WIND);
         ev.push_back (MERIDIONAL_WIND);
         return get_pressure_level_data_ptr (key, level.value, ev);
      }

      case Level::THETA:
      {
         vector<Met_Element> ev;
         ev.push_back (ZONAL_WIND);
         ev.push_back (MERIDIONAL_WIND);
         return get_theta_level_data_ptr (key, level.value, ev, false);
      }

      case Level::SIGMA:
      {
         vector<Met_Element> ev;
         ev.push_back (ZONAL_WIND);
         ev.push_back (MERIDIONAL_WIND);
         return get_sigma_level_data_ptr (key, level.value, ev);
      }

      case Level::FIFTY_METRE:
      case Level::TEN_METRE:
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

      case Level::PRESSURE:
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
   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real latitude = data_ptr->get_coordinate (0, i);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real longitude = data_ptr->get_coordinate (0, longitude);

         const Real next_p = next_p_data_ptr->get_datum (0, i, j);
         const Real next_t = data_3d.evaluate (TEMPERATURE,
            next_p, latitude, longitude);
         const Real next_td = data_3d.evaluate (DEW_POINT,
            next_p, latitude, longitude);
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

      case Level::PRESSURE:
      {

         vector<Met_Element> ev;
         ev.push_back (TEMPERATURE);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case Level::THETA:
      {

         vector<Met_Element> ev;
         ev.push_back (TEMPERATURE);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case Level::SCREEN:
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

      case Level::PRESSURE:
      {
         vector<Met_Element> ev;
         ev.push_back (DEW_POINT);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case Level::THETA:
      {

         vector<Met_Element> ev;
         ev.push_back (DEW_POINT);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case Level::SCREEN:
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

      case Level::PRESSURE:
      {
         vector<Met_Element> ev;
         ev.push_back (RELATIVE_HUMIDITY);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case Level::THETA:
      {

         vector<Met_Element> ev;
         ev.push_back (RELATIVE_HUMIDITY);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case Level::SCREEN:
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

      case Level::PRESSURE:
      {
         vector<Met_Element> ev;
         ev.push_back (DEW_POINT_DEPRESSION);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case Level::THETA:
      {

         vector<Met_Element> ev;
         ev.push_back (DEW_POINT_DEPRESSION);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case Level::SCREEN:
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

      case Level::PRESSURE:
      {
         vector<Met_Element> ev;
         ev.push_back (THETA_E);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case Level::THETA:
      {

         vector<Met_Element> ev;
         ev.push_back (THETA_E);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case Level::SCREEN:
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

      case Level::PRESSURE:
      {
         vector<Met_Element> ev;
         ev.push_back (THETA_W);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_pressure_level_data_ptr (key, level.value, ev);

      }

      case Level::THETA:
      {

         vector<Met_Element> ev;
         ev.push_back (THETA_W);

         if (with_wind)
         {
            ev.push_back (ZONAL_WIND);
            ev.push_back (MERIDIONAL_WIND);
         }

         return get_theta_level_data_ptr (key, level.value, ev, false);

      }

      case Level::SCREEN:
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

   if (level.type != Level::PRESSURE)
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
Nwp::get_isallobar_data_ptr (const Key& key)
{

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);

   Geodetic_Vector_Data_2D* mslp_data_ptr = get_initialized_vd_2d (2);
   const Level& msl = Level::mean_sea_level ();
   const Key& prev_key = key_multimap.get_previous_key (key);

   fill_data (*mslp_data_ptr, 0, prev_key, msl, PRESSURE);
   fill_data (*mslp_data_ptr, 1, key, msl, PRESSURE);

   const Size_2D& size_2d = data_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real this_mslp = mslp_data_ptr->get_datum (1, i, j); 
         const Real prev_mslp = mslp_data_ptr->get_datum (0, i, j); 
         data_ptr->set_datum (0, i, j, this_mslp - prev_mslp);
      }
   }

   delete mslp_data_ptr;
   return data_ptr;

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
                               const vector<Met_Element> element_vector,
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

            const Met_Element& T = TEMPERATURE;
            const Met_Element& U = ZONAL_WIND;
            const Met_Element& V = MERIDIONAL_WIND;

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

            const Met_Element& ne = element_vector[ee];
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
                               const vector<Met_Element> element_vector)
{

   const Integer n = element_vector.size ();
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (n);

   const Level& surface = Level::surface_level ();
   Geodetic_Vector_Data_2D* surface_p_data_ptr = get_initialized_vd_2d (1);
   fill_data (*surface_p_data_ptr, 0, key, surface, PRESSURE);
   const Geodetic_Vector_Data_2D& surf_p_data = *surface_p_data_ptr;

   for (Integer e = 0; e < n; e++)
   {
      const Met_Element& ne = element_vector[e];
      fill_sigma_level_data (*data_ptr, e, key, sigma, ne, surf_p_data);
   }

   delete surface_p_data_ptr;
   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_pressure_level_data_ptr (const Key& key,
                                  const Real p,
                                  const vector<Met_Element> element_vector)
{

   const Integer n = element_vector.size ();

   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (n);

   for (Integer e = 0; e < n; e++)
   {
      const Met_Element& met_element = element_vector[e];
      fill_pressure_level_data (*data_ptr, e, key, p, met_element);
   }

   return data_ptr;

}

Geodetic_Vector_Data_2D*
Nwp::get_data_ptr (const Key& key,
                   const Level& level,
                   const Met_Element met_element)
{
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (1);
   fill_data (*data_ptr, 0, key, level, met_element);
   return data_ptr;
}

Geodetic_Vector_Data_2D*
Nwp::get_data_ptr (const Key& key,
                   const Level& level,
                   const vector<Met_Element>& met_element_vector)
{

   const Integer n = met_element_vector.size ();
   Geodetic_Vector_Data_2D* data_ptr = get_initialized_vd_2d (n);

   for (Integer i = 0; i < n; i++)
   {
      const Met_Element& met_element = met_element_vector[i];
      fill_data (*data_ptr, i, key, level, met_element);
   }

   return data_ptr;

}

Scalar_Data_1D*
Nwp::get_terrain_profile_ptr (const Key& key,
                              const Journey& journey)
{

   const Geodesy geodesy;
   const Tuple& tuple_x = journey.get_tuple_x (geodesy);
   const Real distance = tuple_x.back ();
   const Level& surface = Level::surface_level ();

   if (tuple_x.size () < 2) { throw Nwp_Exception ("Invalid journey"); }
   if (gsl_isnan (distance)) { throw Nwp_Exception ("Invalid journey"); }
   if (distance < 1) { throw Nwp_Exception ("journey too short"); }

   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (1);
   fill_data (*temp_data_ptr, 0, key, surface, PRESSURE);

   Scalar_Data_1D* data_ptr = new Scalar_Data_1D (tuple_x);

   for (Journey::const_iterator iterator = journey.begin ();
        iterator != journey.end (); iterator++)
   {

      const Integer i = std::distance (journey.begin (), iterator);
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
                               const Journey& journey)
{

   const Geodesy geodesy;
   const Tuple& tuple_x = journey.get_tuple_x (geodesy);
   const Real distance = tuple_x.back ();

   if (tuple_x.size () < 2) { throw Nwp_Exception ("Invalid journey"); }
   if (gsl_isnan (distance)) { throw Nwp_Exception ("Invalid journey"); }
   if (distance < 1) { throw Nwp_Exception ("journey too short"); }

   const Level& nil = Level::nil_level ();
   Geodetic_Vector_Data_2D* temp_data_ptr = get_initialized_vd_2d (1);
   fill_data (*temp_data_ptr, 0, key, nil, RAINFALL_STEP);

   Scalar_Data_1D* data_ptr = new Scalar_Data_1D (tuple_x);

   for (Journey::const_iterator iterator = journey.begin ();
        iterator != journey.end (); iterator++)
   {

      const Integer i = std::distance (journey.begin (), iterator);
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
                            const Journey& journey,
                            const Met_Element met_element,
                            const bool with_wind)
{
   vector<Met_Element> met_element_vector;
   met_element_vector.push_back (met_element);
   return get_cross_section_ptr (key, journey,
      met_element_vector, with_wind);
}

Nwp::Cross_Section*
Nwp::get_cross_section_ptr (const Key& key,
                            const Journey& journey,
                            const vector<Met_Element>& met_element_vector,
                            const bool with_wind)
{

   const Geodesy geodesy;
   const Tuple& tuple_x = journey.get_tuple_x (geodesy);
   const Real distance = tuple_x.back ();

   if (tuple_x.size () < 2) { throw Nwp_Exception ("Invalid journey"); }
   if (gsl_isnan (distance)) { throw Nwp_Exception ("Invalid journey"); }
   if (distance < 1) { throw Nwp_Exception ("journey too short"); }

   const Data_3D& data_3d = get_3d_data (key);

   const Integer nn = met_element_vector.size ();
   const Integer nn_1 = nn + 1;
   const Integer nn_2 = nn + 2;

   const Integer n = (with_wind ? nn + 3 : nn);
   Cross_Section* cs_ptr = new Cross_Section ();

   typedef vector<Met_Element>::const_iterator Iterator;
   for (Iterator iterator = met_element_vector.begin ();
        iterator != met_element_vector.end (); iterator++)
   {

      const Met_Element& met_element = *(iterator);
      const Tuple& tuple_p = data_3d.get_tuple_p (met_element);

      for (Journey::const_iterator i = journey.begin ();
           i != journey.end (); i++)
      {

         const Integer ii = std::distance (journey.begin (), i);
         const Lat_Long& lat_long = *(i);
         const Real latitude = lat_long.latitude;
         const Real longitude = lat_long.longitude;
         const Real angle = journey.get_azimuth_forward (i, geodesy);
         const Real theta = angle * DEGREE_TO_RADIAN;
         const Real c = cos (theta);
         const Real s = sin (theta);

         // This is only computed for ne = LI_THUNDER
         const Real thunder_t = -20 + K;
         const Real thunder_p = data_3d.get_p_from_element (
            TEMPERATURE, latitude, longitude, thunder_t);

         cs_ptr->insert_met_element_if_needed (met_element, tuple_x, tuple_p);
         Scalar_Data_2D& sd_2d = cs_ptr->get_sd_2d (met_element);

         for (Integer k = 0; k < tuple_p.size (); k++)
         {

            const Real p = tuple_p[k];

            try
            {
               const Real datum = (met_element == LI_THUNDER) ?
                  data_3d.get_li_thunder (p, lat_long, thunder_p, thunder_t) :
                  data_3d.evaluate (met_element, p, lat_long);
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

      const Met_Element U = ZONAL_WIND;
      const Met_Element V = MERIDIONAL_WIND;
      const Met_Element O = OMEGA;

      const Tuple& tuple_p = data_3d.get_tuple_p (ZONAL_WIND);
      cs_ptr->insert_met_element_if_needed (STREAMLINE_WIND, tuple_x, tuple_p);
      cs_ptr->insert_met_element_if_needed (NORMAL_WIND, tuple_x, tuple_p);
      cs_ptr->insert_met_element_if_needed (OMEGA, tuple_x, tuple_p);
      Scalar_Data_2D& s_sd_2d = cs_ptr->get_sd_2d (STREAMLINE_WIND);
      Scalar_Data_2D& n_sd_2d = cs_ptr->get_sd_2d (NORMAL_WIND);
      Scalar_Data_2D& o_sd_2d = cs_ptr->get_sd_2d (OMEGA);

      for (Journey::const_iterator i = journey.begin ();
           i != journey.end (); i++)
      {

         const Integer ii = std::distance (journey.begin (), i);
         const Lat_Long& lat_long = *(i);
         const Real latitude = lat_long.latitude;
         const Real longitude = lat_long.longitude;
         const Real angle = journey.get_azimuth_forward (i, geodesy);
         const Real theta = angle * DEGREE_TO_RADIAN;
         const Real c = cos (theta);
         const Real s = sin (theta);

         for (Integer k = 0; k < tuple_p.size (); k++)
         {

            const Real p = tuple_p[k];

            try
            {
               const Real u = data_3d.evaluate (U, p, latitude, longitude);
               const Real v = data_3d.evaluate (V, p, latitude, longitude);
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
               const Real o = data_3d.evaluate (O, p, latitude, longitude);
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
   Sd_1d* terrain_profile_ptr = get_terrain_profile_ptr (key, journey);
   Sd_1d* rainfall_profile_ptr = get_rainfall_profile_ptr (key, journey);
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
      const Real t = data_3d.evaluate (TEMPERATURE, p, latitude, longitude);
      const bool invalid = gsl_isnan (t) || (t > 350) || (t < -150);
      if (!invalid) { sounding_ptr->get_t_line ().add (p, t - K); }
   }

   for (Integer k = 0; k < td_tuple_p.size (); k++)
   {
      const Real p = td_tuple_p[k];
      if (p > surface_p) { continue; }
      const Real t_d = data_3d.evaluate (DEW_POINT, p, latitude, longitude);
      const bool invalid = gsl_isnan (t_d) || (t_d > 350) || (t_d < -150);
      if (!invalid) { sounding_ptr->get_t_d_line ().add (p, t_d - K); }
   }

   for (Integer k = 0; k < z_tuple_p.size (); k++)
   {
      const Real p = z_tuple_p[k];
      if (p > surface_p) { continue; }
      const Real z = data_3d.evaluate (GEOPOTENTIAL_HEIGHT, p, latitude, longitude);
      const bool invalid = gsl_isnan (z) || (z > 30000) || (z < -500);
      if (!invalid) { sounding_ptr->get_height_profile ().add (p, z); }
   }

   for (Integer k = 0; k < u_tuple_p.size (); k++)
   {
      const Real p = u_tuple_p[k];
      if (p > surface_p) { continue; }
      const Real u = data_3d.evaluate (ZONAL_WIND, p, latitude, longitude);
      const Real v = data_3d.evaluate (MERIDIONAL_WIND, p, latitude, longitude);
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
         const Met_Element& met_element = level_element.met_element;

         try
         {
            fill_data (*data_ptr, 0, key, level, met_element);
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
                              const vector<Met_Element>& met_element_vector)
{

   const Tuple& tuple_t = get_valid_t_tuple (base_time, start_time, end_time);
   if (tuple_t.size () < 2) { return NULL; }

   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;

   typedef Scalar_Data_2D Sd_2d;
   typedef vector<Met_Element>::const_iterator Iterator;
   Nwp::Time_Cross* tc_ptr = new Time_Cross ();

   for (Integer i = 0; i < tuple_t.size (); i++)
   {

      const Dtime& dtime = tuple_t[i];
      const Key& key = get_key (dtime, base_time);

      for (Iterator iterator = met_element_vector.begin ();
           iterator != met_element_vector.end (); iterator++)
      {

         const Data_3D& data_3d = get_3d_data (key);
         const Met_Element& met_element = *(iterator);

         const Real thunder_t = -20 + K;
         const Real thunder_p = data_3d.get_p_from_element (
            TEMPERATURE, latitude, longitude, thunder_t);

         const Tuple& tuple_p = data_3d.get_tuple_p (met_element);
         tc_ptr->insert_met_element_if_needed (met_element, tuple_t, tuple_p);
         Scalar_Data_2D& sd_2d = tc_ptr->get_sd_2d (met_element);

         for (Integer k = 0; k < tuple_p.size (); k++)
         {

            const Real p = tuple_p[k];

            try
            {

               const Real datum = (met_element == LI_THUNDER) ?
                  data_3d.get_li_thunder (p, lat_long, thunder_p, thunder_t) :
                  data_3d.evaluate (met_element, p, lat_long);
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
   //Sd_1d* terrain_profile_ptr = get_terrain_profile_ptr (key, journey);
   //Sd_1d* rainfall_profile_ptr = get_rainfall_profile_ptr (key, journey);
   //cs_ptr->set_profile_ptrs (terrain_profile_ptr, rainfall_profile_ptr);
   return tc_ptr;

}

Nwp_Exception::Nwp_Exception (const Dstring& str)
   : Exception ("Nwp_Exception", str)
{
}

