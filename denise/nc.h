//
// nc.h
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

#ifndef DENISE_NC_H
#define DENISE_NC_H

#include <map>
#include <netcdf.h>
#include <denise/gis.h>
#include <denise/geodesy.h>
#include <denise/nwp.h>
#include <denise/thermo.h>

using namespace std;

namespace denise
{

   class Ncep_Ncar : public Nwp
   {

      private:

         class File
         {

            protected:

               int
               nc_id;

               const Dstring&
               file_path;

               Integer
               get_dim_size (const Dstring& dim_str) const;

            public:

               File (const Dstring& file_path);

               ~File ();

               const Tuple&
               get_tuple_latitude () const;

               const Tuple&
               get_tuple_longitude () const;

         };

         class File_3D : public File
         {

            private:

               Integer
               actual_number_of_levels;

            public:

               File_3D (const Dstring& file_path);

               void
               fill_data (Geodetic_Vector_Data_3D& data_3d,
                          const Integer vector_index,
                          const Dstring& variable_string,
                          const Integer time_index) const;

               void
               fill_data (Geodetic_Vector_Data_2D& data_2d,
                          const Integer vector_index,
                          const Dstring& variable_string,
                          const Integer time_index,
                          const Integer level_index) const;

         };

         class File_2D : public File
         {

            public:

               File_2D (const Dstring& file_path);

               void
               fill_data (Geodetic_Vector_Data_2D& data_2d,
                          const Integer vector_index,
                          const Dstring& variable_string,
                          const Integer time_index) const;

         };


         const Dtime
         start_time;

         const Real
         time_span;

         const Integer
         time_interval;

         class Data_3D : public Nwp::Data_3D
         {

            public:

               Data_3D (const vector<Met_Element>& met_element_vector,
                        const Key& key);

               virtual Real
               evaluate (const Met_Element element,
                         const Real p,
                         const Real latitude,
                         const Real longitude,
                         const Evaluate_Op evaluate_op = VALUE) const;


         };

         const Dstring
         data_path;

         Tuple
         tuple_latitude_gaussian;

         Tuple
         tuple_longitude_gaussian;

         Tuple
         tuple_p_omega;

         Tuple
         tuple_p_humid;

         Dstring
         get_file_path (const Met_Element met_element,
                        const Dtime& time) const;

         Dstring
         get_element_string (const Met_Element met_element) const;

         Integer
         get_time_index (const Dtime& time) const;

         Integer
         get_p_index (const Met_Element met_element,
                      const Real p) const;

         const Tuple&
         get_tuple_p (const Met_Element met_element) const;

         Geodetic_Vector_Data_2D*
         get_initialized_vd_2d (const Integer vector_size) const;

         void
         initialize_3d_data (const Key& key);

         void
         load_3d_data (Nwp::Data_3D& data_3d);

         void
         fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                    const Integer vector_index,
                                    const Key& key);

         void
         fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key);

         void
         fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Met_Element Met_element);

         void
         fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Met_Element met_element);

         void
         fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                        const Integer vector_index,
                        const Key& key,
                        const Met_Element met_element);
         void
         fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                  const Integer vector_index,
                                  const Key& key,
                                  const Met_Element met_element);

         void
         fill_nc_data (Geodetic_Vector_Data_3D& gvd_3d,
                       const Integer vector_index,
                       const Met_Element met_element,
                       const Dtime& time) const;

         void
         fill_nc_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Met_Element met_element,
                       const Dtime& time,
                       const Real p = GSL_NAN) const;

         void
         init ();

      public:

         Ncep_Ncar (const Dstring& data_path,
                    const Dtime& start_time,
                    const Real time_span,
                    const Integer time_interval = 21600);

         Ncep_Ncar (const Dstring& data_path,
                    const Dtime& start_time,
                    const Dtime& end_time,
                    const Integer time_interval = 21600);

         void
         survey ();

         void
         clean_up ();

         static Tuple
         tuple_latitude ();

         static Tuple
         tuple_longitude ();

         const Tuple
         get_tuple_latitude () const;

         const Tuple
         get_tuple_longitude () const;

         set<Dtime>
         get_valid_time_set () const;

         Key
         get_key (const Dtime& dtime) const;

         void
         acquire_base_time_forecast_second (Dtime& base_time,
                                            Integer& forecast_second,
                                            const Dtime& dtime) const;

         Geodetic_Vector_Data_2D*
         get_surface_data_2d_ptr (const Key& key);

   };

};

#endif /* DENISE_NC_H */ 
