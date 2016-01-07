//
// nc.cc
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
// but WITHOUT ANY WARRANTY), without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with libdenise.  If not, see <http://www.gnu.org/licenses/>.

#include "nc.h"
#include <denise/thermo.h>

using namespace std;
using namespace denise;

Integer
Ncep_Ncar::File::get_dim_size (const Dstring& dim_str) const
{

   int ret, dim_id;
   size_t dim_length;
   const char* str = dim_str.c_str ();

   ret = nc_inq_dimid (nc_id, dim_str.c_str (), &dim_id);
   if (ret != NC_NOERR) { throw Exception ("nc_inq_dimid " + dim_str); }

   ret = nc_inq_dimlen (nc_id, dim_id, &dim_length);
   if (ret != NC_NOERR) { throw Exception ("nc_inq_dimlen " + dim_str); }

   return Integer (dim_length);

}

Ncep_Ncar::File::File (const Dstring& file_path)
   : file_path (file_path)

{

   int ret, dim_id;
   size_t latitude_dim_length, longitude_dim_length;

   ret = nc_open (file_path.c_str (), NC_NOWRITE, &nc_id);
   if (ret != NC_NOERR) { throw Exception ("nc_open " + file_path); }

}

Ncep_Ncar::File::~File ()
{
   int ret = nc_close (nc_id);
   if (ret != NC_NOERR) { throw Exception ("nc_close " + file_path); }
}

Ncep_Ncar::File_3D::File_3D (const Dstring& file_path)
   : Ncep_Ncar::File (file_path),
     actual_number_of_levels (get_dim_size ("level"))
{
}

void
Ncep_Ncar::File_3D::fill_data (Geodetic_Vector_Data_3D& data_3d,
                               const Integer vector_index,
                               const Dstring& variable_string,
                               const Integer time_index) const
{

   int r, var_id;
   float offset, scale;

   const Tuple tuple_level = data_3d.get_coordinate_tuple (0);
   const Tuple tuple_latitude = data_3d.get_coordinate_tuple (1);
   const Tuple tuple_longitude = data_3d.get_coordinate_tuple (2);

   const size_t nk = tuple_level.size ();
   const size_t ni = tuple_latitude.size ();
   const size_t nj = tuple_longitude.size () - 1;
   size_t n = nk * ni * nj;

   short* array = new short[n];
   size_t start[] = { size_t (time_index), 0, 0, 0 };
   size_t count[] = { 1, size_t (actual_number_of_levels), ni, nj };

   const string& vs = variable_string;

   r = nc_inq_varid (nc_id, vs.c_str (), &var_id);
   if (r != NC_NOERR) { throw Exception ("nc_inq_varid " + vs); }

   r = nc_get_vara_short (nc_id, var_id, start, count, array);
   if (r != NC_NOERR) { throw Exception ("nc_get_vara_short " + vs); }

   r = nc_get_att_float (nc_id, var_id, "add_offset", &offset);
   if (r != NC_NOERR) { throw Exception ("nc_get_att_float offset " + vs); }

   r = nc_get_att_float (nc_id, var_id, "scale_factor", &scale);
   if (r != NC_NOERR) { throw Exception ("nc_get_att_float scale " + vs); }

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
Ncep_Ncar::File_3D::fill_data (Geodetic_Vector_Data_2D& data_2d,
                               const Integer vector_index,
                               const Dstring& variable_string,
                               const Integer time_index,
                               const Integer level_index) const
{

   int r, var_id;
   float offset, scale;

   const Tuple tuple_latitude = data_2d.get_coordinate_tuple (0);
   const Tuple tuple_longitude = data_2d.get_coordinate_tuple (1);

   const size_t ni = tuple_latitude.size ();
   const size_t nj = tuple_longitude.size () - 1;
   size_t n = ni * nj;

   short* array = new short[n];
   size_t start[] = { size_t (time_index), size_t (level_index), 0, 0 };
   size_t count[] = { 1, 1, ni, nj };

   const string& vs = variable_string;

   r = nc_inq_varid (nc_id, vs.c_str (), &var_id);
   if (r != NC_NOERR) { throw Exception ("nc_inq_varid " + vs); }

   r = nc_get_vara_short (nc_id, var_id, start, count, array);
   if (r != NC_NOERR) { throw Exception ("nc_get_vara_short " + vs); }

   r = nc_get_att_float (nc_id, var_id, "add_offset", &offset);
   if (r != NC_NOERR) { throw Exception ("nc_get_att_float offset " + vs); }

   r = nc_get_att_float (nc_id, var_id, "scale_factor", &scale);
   if (r != NC_NOERR) { throw Exception ("nc_get_att_float scale " + vs); }

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

Ncep_Ncar::File_2D::File_2D (const Dstring& file_path)
   : Ncep_Ncar::File (file_path)
{
}

void
Ncep_Ncar::File_2D::fill_data (Geodetic_Vector_Data_2D& data_2d,
                               const Integer vector_index,
                               const Dstring& variable_string,
                               const Integer time_index) const
{

   int r, var_id;
   float offset, scale;

   const Tuple tuple_latitude = data_2d.get_coordinate_tuple (0);
   const Tuple tuple_longitude = data_2d.get_coordinate_tuple (1);

   const size_t ni = tuple_latitude.size ();
   const size_t nj = tuple_longitude.size () - 1;
   size_t n = ni * nj;

   short* array = new short[n];
   size_t start[] = { size_t (time_index), 0, 0 };
   size_t count[] = { 1, ni, nj };

   const string& vs = variable_string;

   r = nc_inq_varid (nc_id, vs.c_str (), &var_id);
   if (r != NC_NOERR) { throw Exception ("nc_inq_varid " + vs); }

   r = nc_get_vara_short (nc_id, var_id, start, count, array);
   if (r != NC_NOERR) { throw Exception ("nc_get_vara_short " + vs); }

   r = nc_get_att_float (nc_id, var_id, "add_offset", &offset);
   if (r != NC_NOERR) { throw Exception ("nc_get_att_float offset " + vs); }

   r = nc_get_att_float (nc_id, var_id, "scale_factor", &scale);
   if (r != NC_NOERR) { throw Exception ("nc_get_att_float scale " + vs); }

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

Ncep_Ncar::Data_3D::Data_3D (const vector<Met_Element>& met_element_vector,
                             const Key& key)
   : Nwp::Data_3D (met_element_vector, key)
{
}

Real
Ncep_Ncar::Data_3D::evaluate (const Met_Element met_element,
                              const Real p,
                              const Real latitude,
                              const Real longitude,
                              const Evaluate_Op evaluate_op) const
{

   const Evaluate_Op& eo = evaluate_op;

   switch (met_element)
   {

      case denise::DEW_POINT:
      {
         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& RH = denise::RELATIVE_HUMIDITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real rh = Nwp::Data_3D::evaluate (RH, p, latitude, longitude);
         return Moisture::get_t_d (t, std::max (0.01, rh));
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& O = denise::OMEGA;

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

   return Nwp::Data_3D::evaluate (met_element, p, latitude, longitude, eo);

}

Dstring
Ncep_Ncar::get_file_path (const Met_Element met_element,
                          const Dtime& time) const
{

   string identifier;

   switch (met_element)
   {
      case denise::U:             identifier = "uwnd.";           break;
      case denise::V:             identifier = "vwnd.";           break;
      case denise::OMEGA:         identifier = "omega.";          break;
      case denise::TEMPERATURE:   identifier = "air.";            break;
      case denise::RH:            identifier = "rhum.";           break;
      case denise::Q:             identifier = "shum.";           break;
      case denise::Z:             identifier = "hgt.";            break;
      //case denise::SURFACE_U:     identifier = "uwnd.sig995.";    break;
      //case denise::SURFACE_V:     identifier = "vwnd.sig995.";    break;
      //case denise::SURFACE_OMEGA:   identifier = "omega.sig995.";   break;
      //case denise::SURFACE_P:     identifier = "pres.sfc.";       break;
      //case denise::SURFACE_T:     identifier = "air.sig995.";     break;
      //case denise::SURFACE_RH:    identifier = "rhum.sig995.";    break;
      //case denise::SURFACE_THETA: identifier = "pottmp.sig995.";  break;
      //case denise::PRECIP_WATER:  identifier = "pr_wtr.eatm.";    break;
      //case denise::SLP:           identifier = "slp.";            break;
      //case denise::U_10M:         identifier = "uwnd.10m.gauss."; break;
      //case denise::V_10M:         identifier = "vwnd.10m.gauss."; break;
      //case denise::T_2M:          identifier = "air.2m.";         break;
   }

   const Dtime quantized_time ((rint (time.t / 6) * 6));
   const string year_string = quantized_time.get_string ("%Y");
   const string path_base (data_path + "/" + year_string + "/");

   const Dstring file_path = path_base + identifier + year_string + ".nc";
   return file_path;

}

Dstring
Ncep_Ncar::get_element_string (const Met_Element met_element) const
{

   Dstring str;

   switch (met_element)
   {
      case denise::U:             str = "uwnd";   break;
      case denise::V:             str = "vwnd";   break;
      case denise::OMEGA:         str = "omega";  break;
      case denise::T:             str = "air";    break;
      case denise::RH:            str = "rhum";   break;
      case denise::Q:             str = "shum";   break;
      case denise::Z:             str = "hgt";    break;
      //case denise::SURFACE_U:     str = "uwnd";   break;
      //case denise::SURFACE_V:     str = "vwnd";   break;
      //case denise::SURFACE_OMEGA:   str = "omega";  break;
      //case denise::SURFACE_P:     str = "pres";   break;
      //case denise::SURFACE_T:     str = "air";    break;
      //case denise::SURFACE_RH:    str = "rhum";   break;
      //case denise::SURFACE_THETA: str = "pottmp"; break;
      //case denise::PRECIP_WATER:  str = "pr_wtr"; break;
      //case denise::SLP:           str = "slp";    break;
      //case denise::U_10M:         str = "uwnd";   break;
      //case denise::V_10M:         str = "vwnd";   break;
      //case denise::T_2M:          str = "air";    break;
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
Ncep_Ncar::get_p_index (const Met_Element met_element,
                        const Real p) const
{

   const Tuple tuple_p = get_tuple_p (met_element);
   const Integer n = tuple_p.size ();

   for (Integer i = 0; i < n; i++)
   {
      const Real delta_p = fabs (tuple_p[i] - p);
      if (delta_p < 0.1) { return n - i - 1; }
   }

   throw Exception ("Level not available");

}

const Tuple&
Ncep_Ncar::get_tuple_p (const Met_Element met_element) const
{

   switch (met_element)
   {

      case denise::U:
      case denise::V:
      case denise::T:
      case denise::Z:
      {
         return tuple_p;
      }

      case denise::OMEGA:
      {
         return tuple_p_omega;
      }

      case denise::RH:
      case denise::Q:
      {
         return tuple_p_humid;
      }

      default:
      {
         throw Exception ("Not 3D Element");
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

   Nnd_3d* nnd_3d_ptr = new Nnd_3d (met_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, nnd_3d_ptr));

}

void          
Ncep_Ncar::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Met_Element>::const_iterator Iterator;
   typedef Geodetic_Vector_Data_3D Gvd_3d;
   const Key& key = data_3d.key;
   const Dtime& dtime = key.base_time;
   const Tuple tuple_latitude = Ncep_Ncar::tuple_latitude ();
   const Tuple tuple_longitude = Ncep_Ncar::tuple_longitude ();

   for (Iterator iterator = met_element_vector.begin ();
        iterator != met_element_vector.end (); iterator++)
   {

      const Met_Element& met_element = *(iterator);

      switch (met_element)
      {

         case TEMPERATURE:
         {
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, denise::T, dtime);
            data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
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
            fill_nc_data (gvd_3d, 0, denise::RH, dtime);
            data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
            break;
         }

         case GEOPOTENTIAL_HEIGHT:
         {
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, denise::Z, dtime);
            data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
            break;
         }

         case ZONAL_WIND:
         {
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, denise::U, dtime);
            data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
            break;
         }

         case MERIDIONAL_WIND:
         {
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, denise::V, dtime);
            data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
            break;
         }

         case OMEGA:
         {
            Gvd_3d* gvd_3d_ptr = new Gvd_3d (1, tuple_p,
               tuple_latitude, tuple_longitude, true);
            Gvd_3d& gvd_3d = *gvd_3d_ptr;
            fill_nc_data (gvd_3d, 0, denise::OMEGA, dtime);
            gvd_3d.scale_offset (0, -1, 0);
            data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
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
                                   const Met_Element met_element)
{

   const Level& screen = Level::screen_level ();
   const Level& surface = Level::surface_level ();
   const Dtime& dtime = key.base_time;

   switch (met_element)
   {

      case TEMPERATURE:
      {
         //fill_nc_data (gvd_2d, vector_index, SURFACE_T, dtime);
         return;
      }

      case DEW_POINT:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);

         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, RELATIVE_HUMIDITY);

         const Size_2D& size_2d = data_ptr->get_size_2d ();

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
         //fill_nc_data (gvd_2d, vector_index, SURFACE_RH, dtime);
         return;
      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, met_element);

}

void
Ncep_Ncar::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Met_Element met_element)
{

   const Dtime& dtime = key.base_time;

   switch (met_element)
   {

      case ZONAL_WIND:
      {
         //fill_nc_data (gvd_2d, vector_index, SURFACE_U, dtime);
         return;
      }

      case MERIDIONAL_WIND:
      {
         //fill_nc_data (gvd_2d, vector_index, SURFACE_V, dtime);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, met_element);

}

void
Ncep_Ncar::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Met_Element met_element)
{

   const Dtime& dtime = key.base_time;

   switch (met_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         //fill_nc_data (gvd_2d, vector_index, SLP, dtime);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, met_element);

}

void
Ncep_Ncar::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                    const Integer vector_index,
                                    const Key& key,
                                    const Met_Element met_element)
{

   const Dtime& dtime = key.base_time;

   switch (met_element)
   {

      case PRESSURE:
      {
         //fill_nc_data (gvd_2d, vector_index, SURFACE_P, dtime);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, met_element);

}

void
Ncep_Ncar::fill_nc_data (Geodetic_Vector_Data_3D& gvd_3d,
                         const Integer vector_index,
                         const Met_Element met_element,
                         const Dtime& time) const
{

   const Dstring& file_path = get_file_path (met_element, time);

   const Dstring& element_string = get_element_string (met_element);
   const Integer time_index = get_time_index (time);

   switch (met_element)
   {

      case denise::U:
      case denise::V:
      case denise::OMEGA:
      case denise::T:
      case denise::RH:
      case denise::Q:
      case denise::Z:
      {
         const Ncep_Ncar::File_3D file (file_path);
         file.fill_data (gvd_3d, vector_index, element_string, time_index);
         break;
      }

   }

}

void
Ncep_Ncar::fill_nc_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Met_Element met_element,
                         const Dtime& time,
                         const Real p) const
{

   const Dstring& file_path = get_file_path (met_element, time);

   const Dstring& element_string = get_element_string (met_element);
   const Integer time_index = get_time_index (time);

   switch (met_element)
   {

      case denise::U:
      case denise::V:
      case denise::OMEGA:
      case denise::T:
      case denise::RH:
      case denise::Q:
      case denise::Z:
      {
         const Ncep_Ncar::File_3D file (file_path);
         const Integer pi = get_p_index (met_element, p);
         file.fill_data (gvd_2d, vector_index, element_string, time_index, pi);
         break;
      }

      //case denise::SURFACE_U:
      //case denise::SURFACE_V:
      //case denise::SURFACE_OMEGA:
      //case denise::SURFACE_P:
      //case denise::SURFACE_T:
      //case denise::SURFACE_RH:
      //case denise::SURFACE_THETA:
      //case denise::PRECIP_WATER:
      //case denise::SLP:
      //case denise::U_10M:
      //case denise::V_10M:
      //case denise::T_2M:
      //{
      //   const Ncep_Ncar_2d_File file (file_path);
      //   file.fill_data (gvd_2d, vector_index, element_string, time_index);
      //   break;
      //}

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

Ncep_Ncar::Ncep_Ncar (const Dstring& data_path,
                      const Dtime& start_time,
                      const Real time_span,
                      const Integer time_interval)
   : Nwp ("Ncep_Ncar", data_path),
     data_path (data_path),
     start_time (start_time),
     time_span (time_span),
     time_interval (time_interval)
{

   init ();

   met_element_vector.push_back (denise::TEMPERATURE);
   met_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   met_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   met_element_vector.push_back (denise::ZONAL_WIND);
   met_element_vector.push_back (denise::MERIDIONAL_WIND);
   met_element_vector.push_back (denise::OMEGA);

   survey ();

}

Ncep_Ncar::Ncep_Ncar (const Dstring& data_path,
                      const Dtime& start_time,
                      const Dtime& end_time,
                      const Integer time_interval)
   : Nwp ("Ncep_Ncar", data_path),
     data_path (data_path),
     start_time (start_time),
     time_span (end_time.t - start_time.t),
     time_interval (time_interval)
{

   init ();

   met_element_vector.push_back (denise::TEMPERATURE);
   met_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   met_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   met_element_vector.push_back (denise::ZONAL_WIND);
   met_element_vector.push_back (denise::MERIDIONAL_WIND);
   met_element_vector.push_back (denise::OMEGA);

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
Ncep_Ncar::acquire_base_time_forecast_second (Dtime& base_time,
                                              Integer& forecast_second,
                                              const Dtime& dtime) const
{
   base_time.t = dtime.t;
   forecast_second = 0;
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

