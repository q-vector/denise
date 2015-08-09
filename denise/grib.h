//
// grib.h
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

#ifndef DENISE_GRIB_H
#define DENISE_GRIB_H

#include <map>
#include <denise/gis.h>
#include <denise/geodesy.h>
#include <denise/nwp.h>
#include <denise/thermo.h>

using namespace std;

namespace denise
{

   class Grib
   {

      public:

         class Block
         {

            public:

               uint8_t*
               buffer;

               uint32_t
               n;

               Block ();

               Block (const uint32_t n,
                      const bool set_zero = false);

               Block (FILE* file,
                      const uint32_t address);

               Block (const Block& block);

               ~Block ();

               static uint32_t
               get_size (FILE* file,
                         const uint32_t address);

               bool
               get_bit (const uint32_t position,
                        const uint8_t bit_number) const;

               void
               set_bit (const uint32_t position,
                        const uint8_t bit_number);

               void
               unset_bit (const uint32_t position,
                          const uint8_t bit_number);

               static uint32_t
               get_uint (FILE* file,
                         const uint32_t address,
                         const uint32_t n);

               static int32_t
               get_int (FILE* file,
                        const uint32_t address,
                        const uint32_t n);

               uint32_t
               get_uint (const uint32_t position,
                         const uint32_t n) const;

               int32_t
               get_int (const uint32_t position,
                        const uint32_t n) const;

               uint32_t
               get_uint (const uint32_t position,
                         const uint8_t start_bit,
                         const uint8_t number_of_bits) const;

               int32_t
               get_int (const uint32_t position,
                        const uint8_t start_bit,
                        const uint8_t number_of_bits) const;

               float
               get_float (const uint32_t position) const;

               bool
               operator== (const Block& block) const;

               bool
               operator> (const Block& block) const;

               bool
               operator< (const Block& block) const;

         };

         class Section
         {

            private:

               const Block&
               block;

               uint32_t
               offset;

               uint32_t
               length;

            public:

               Section (const Block& block,
                        const uint32_t offset,
                        const uint32_t length);

               uint32_t
               get_uint (const uint32_t position,
                         const uint32_t n) const;

               int32_t
               get_int (const uint32_t position,
                        const uint32_t n) const;

               uint32_t
               get_uint (const uint32_t position,
                         const uint8_t start_bit,
                         const uint8_t number_of_bits) const;

               int32_t
               get_int (const uint32_t position,
                        const uint8_t start_bit,
                        const uint8_t number_of_bits) const;

               float
               get_float (const uint32_t position) const;

         };

         class Pds : public Block
         {

            public:

               class Level : public Section
               {

                  public:

                     Level (const Pds& pds);

               };

               class Forecast_Time : public Section
               {

                  public:

                     Forecast_Time (const Pds& pds);

                     uint8_t
                     get_p1 () const;

                     uint8_t
                     get_p2 () const;

               };

               Pds (FILE* file,
                    const uint32_t address);

               uint8_t
               get_center_id () const;

               uint8_t
               get_grid_id () const;

               uint8_t
               get_element_id () const;

               bool
               is_gds_present () const;

               bool
               is_bms_present () const;

               Dtime
               get_base_time () const;

               Grib::Pds::Level
               get_level () const;

               Grib::Pds::Forecast_Time
               get_forecast_time () const;

               int16_t
               get_D () const;

         };

         class Gds : public Block
         {

            public:

               Gds (FILE* file,
                    const uint32_t address);

               const Size_2D
               get_size_2d () const;

         };

         class Bms : public Block
         {

            public:

               Bms (FILE* file,
                    const uint32_t address);

         };

         class Bds : public Block
         {

            public:

               Bds (FILE* file,
                    const uint32_t address);

               float
               get_float (const uint32_t position,
                          const uint8_t start_bit,
                          const uint8_t number_of_bits,
                          const int16_t E,
                          const float R,
                          const int16_t D) const;

         };

         class Key : public Block
         {

            public:

               Key ();

               Key (const Pds& pds);

               Key (const Key& key);

               ostream&
               operator << (ostream& out);

         };

         class Header
         {

            public:

               uint32_t
               address;

               uint32_t
               record_size;

               uint32_t
               gds_position;

               uint32_t
               bms_position;

               uint32_t
               bds_position;

               Pds*
               pds_ptr;

               Gds*
               gds_ptr;

               Header (FILE* file,
                       uint32_t& offset);

               ~Header ();

               const Pds&
               get_pds () const;

               const Gds&
               get_gds () const;

         };

      private:

         const string
         file_path;

         map<Grib::Key, Grib::Header*>
         header_ptr_map;

      public:

         Grib (const string& file_path);

         ~Grib ();

         const map<Grib::Key, Grib::Header*>&
         get_header_ptr_map () const;

         Gds*
         get_gds_ptr (FILE* file,
                      const Key& key) const;

         Bms*
         get_bms_ptr (FILE* file,
                      const Key& key) const;

         Bds*
         get_bds_ptr (FILE* file,
                      const Key& key) const;

         void
         fill_data (Geodetic_Vector_Data_3D& gvd_3d,
                    const Integer element_index,
                    const Integer k,
                    const Key& key) const;

         void
         fill_data (Geodetic_Vector_Data_2D& gvd_2d,
                    const Integer element_index,
                    const Key& key) const;

   };
   
   class Grib2
   {

      public:

         class Block
         {

            public:

               uint8_t*
               buffer;

               uint32_t
               n;

               Block ();

               Block (const uint32_t n,
                      const bool set_zero = false);

               Block (FILE* file,
                      uint32_t& address);

               Block (const Block& block);

               ~Block ();

               static uint32_t
               get_size (FILE* file,
                         const uint32_t address);

               bool
               get_bit (const uint32_t position,
                        const uint8_t bit_number) const;

               void
               set_bit (const uint32_t position,
                        const uint8_t bit_number);

               void
               unset_bit (const uint32_t position,
                          const uint8_t bit_number);

               static uint64_t
               get_uint (FILE* file,
                         const uint32_t address,
                         const uint32_t n);

               static int64_t
               get_int (FILE* file,
                        const uint32_t address,
                        const uint32_t n);

               uint64_t
               get_uint (const uint32_t position,
                         const uint32_t n) const;

               int64_t
               get_int (const uint32_t position,
                        const uint32_t n) const;

               uint32_t
               get_uint (const uint32_t position,
                         const uint8_t start_bit,
                         const uint8_t number_of_bits) const;

               int32_t
               get_int (const uint32_t position,
                        const uint8_t start_bit,
                        const uint8_t number_of_bits) const;

               float
               get_float (const uint32_t position) const;

               bool
               operator== (const Block& block) const;

               bool
               operator> (const Block& block) const;

               bool
               operator< (const Block& block) const;

         };

         class Section
         {

            public:

               const Block&
               block;

               uint32_t
               offset;

               uint32_t
               length;

               Section (const Block& block,
                        const uint32_t offset,
                        const uint32_t length);

               uint32_t
               get_uint (const uint32_t position,
                         const uint32_t n) const;

               int32_t
               get_int (const uint32_t position,
                        const uint32_t n) const;

               uint32_t
               get_uint (const uint32_t position,
                         const uint8_t start_bit,
                         const uint8_t number_of_bits) const;

               int32_t
               get_int (const uint32_t position,
                        const uint8_t start_bit,
                        const uint8_t number_of_bits) const;

               float
               get_float (const uint32_t position) const;

         };

         class Block_1 : public Block
         {

            public:

               class Base_Time : public Section
               {

                  public:

                     Base_Time (const Block_1& block_1);

               };

               Block_1 (FILE* file,
                        uint32_t& address);

               Dtime
               get_reference_time () const;

               Grib2::Block_1::Base_Time
               get_base_time () const;

         };

         class Block_2 : public Block
         {

            public:

               Block_2 (FILE* file,
                        uint32_t& address);

         };

         class Block_3 : public Block
         {

            private:

               uint16_t
               get_template_number () const;

               uint8_t
               get_scanning_mode () const;

            public:

               Block_3 (FILE* file,
                        uint32_t& address);

               uint32_t
               get_number_of_points () const;

               const Size_2D
               get_size_2d () const;

               const uint8_t
               get_scan_flag () const;

         };

         class Block_4 : public Block
         {

            public:

               class Parameter : public Section
               {

                  public:

                     Parameter (const Block_4& block_4);

               };

               class Level : public Section
               {

                  public:

                     Level (const Block_4& block_4);

                     uint8_t
                     get_first_level_type () const;

                     uint32_t
                     get_first_level () const;

                     uint8_t
                     get_second_level_type () const;

                     uint32_t
                     get_second_level () const;

               };

               class Forecast_Time : public Section
               {

                  public:

                     Forecast_Time (const Block_4& block_4);

               };

            private:

               uint16_t
               get_template_number () const;

            public:

               Block_4 (FILE* file,
                        uint32_t& address);

               Grib2::Block_4::Parameter
               get_parameter () const;

               Grib2::Block_4::Level
               get_level () const;

               Grib2::Block_4::Forecast_Time
               get_forecast_time () const;

               bool
               parameter_is_percentage () const;

         };

         class Block_5 : public Block
         {

            public:

               Block_5 (FILE* file,
                        uint32_t& address);

               uint16_t
               get_template_number () const;

               uint32_t
               get_number_of_points () const;

               float
               get_reference_value () const;

               int16_t
               get_binary_scale_factor () const;

               int16_t
               get_decimal_scale_factor () const;

               uint8_t
               get_bits_per_datum () const;

         };

         class Block_6 : public Block
         {

            public:

               Block_6 (FILE* file,
                        uint32_t& address);

               uint8_t
               get_bitmap_indicator () const;

               bool
               bitmap_present () const;

         };

         class Block_7 : public Block
         {

            public:

               Block_7 (FILE* file,
                        uint32_t& address);

         };

         class Data
         {

            protected:

               const uint32_t
               n;

               float*
               buffer;

            public:

               Data (const uint32_t n);

               ~Data ();

               const uint32_t
               get_n () const;

               float*
               get_buffer () const;

               void
               simple (const Block_5& block_5,
                       const Block_7& block_7);

               void
               jpeg2000 (const Block_5& block_5,
                         const Block_7& block_7);

         };

         class Header
         {

            public:

               uint32_t
               block_7_address;

               const Block_1
               block_1;

               const Block_3
               block_3;

               const Block_4
               block_4;

               const Block_5
               block_5;

               const Block_6
               block_6;

               Header (const Block_1& block_1,
                       const Block_3& block_3,
                       const Block_4& block_4,
                       const Block_5& block_5,
                       const Block_6& block_6,
                       const uint32_t block_7_address);

         };

         class Key : public Block
         {

            public:

               Key ();

               Key (const Header& header);

               Key (const Key& key);

         };

      private:

         const string
         file_path;

         map<Grib2::Key, Grib2::Header*>
         header_ptr_map;

         const Header&
         get_header (const Key& key) const;

         Data*
         data_ptr;

      public:

         Grib2 (const string& file_path);

         ~Grib2 ();

         const map<Grib2::Key, Grib2::Header*>&
         get_header_ptr_map () const;

         void
         fill_data (Geodetic_Vector_Data_3D& gvd_3d,
                    const Integer element_index,
                    const Integer k,
                    const Key& key) const;

         void
         fill_data (Geodetic_Vector_Data_2D& gvd_2d,
                    const Integer element_index,
                    const Key& key) const;

   };

   ostream&
   operator << (ostream &out_file,
                const Grib::Key& key);

   ostream&
   operator << (ostream &out_file,
                const Grib2::Key& key);

   class Access : public Nwp,
                  public map<Grib::Key, Grib*>
   {

      private:

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

         set<Grib*>
         grib_ptr_set;

         const bool
         omega_as_w;

         const string
         data_path;

         const string
         search_string;

         Size_2D
         size_2d;

         Domain_2D
         domain_2d;

         Grib::Key
         get_grib_key (const Key& key,
                       const Met_Element met_element,
                       const Level& level) const;

         void
         set_grib_key (Grib::Key& grib_key,
                       const Met_Element met_element,
                       const Dtime& base_time,
                       const Integer forecast_hour) const;

         void
         set_grib_key (Grib::Key& grib_key,
                       const Met_Element met_element,
                       const denise::Level& level) const;

         void
         initialize_3d_data (const Key& key);

         void
         load_3d_data (Nwp::Data_3D& data_3d);

         Geodetic_Vector_Data_2D*
         get_initialized_vd_2d (const Integer vector_size) const;

         void
         fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Level& level,
                                 const Met_Element met_element);

         void
         fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Key& key,
                         const Met_Element met_element);

         void
         fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Met_Element met_element);

         void
         fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Met_Element met_element);

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

         Geodetic_Vector_Data_3D*
         get_gvd_3d_ptr (const Met_Element met_element,
                         const Key& key) const;

         Geodetic_Vector_Data_2D*
         get_grib_data_ptr (const Grib& grib,
                            const Grib::Key& grib_key) const;

         void
         fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Met_Element met_element,
                         const Key& key,
                         const Level& level) const;

      public:

         Access (const string& description,
                 const string& data_path,
                 const string& search_string,
                 const bool omega_as_w = false);

         ~Access ();

         void
         survey ();

         void
         clean_up ();

   };

   class Ecmwf : public Nwp,
                 public map<Grib::Key, Grib*>
   {

      private:

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

         set<Grib*>
         grib_ptr_set;

         const bool
         omega_as_w;

         const string
         data_path;

         const string
         search_string;

         map<Real, Size_2D>
         size_2d_map;

         Domain_2D
         domain_2d;

         Grib::Key
         get_grib_key (const Key& key,
                       const Met_Element met_element,
                       const Level& level) const;

         void
         set_grib_key (Grib::Key& grib_key,
                       const Met_Element met_element,
                       const Dtime& base_time,
                       const Integer forecast_hour) const;

         void
         set_grib_key (Grib::Key& grib_key,
                       const Met_Element met_element,
                       const denise::Level& level) const;

         void
         initialize_3d_data (const Key& key);

         void
         load_3d_data (Nwp::Data_3D& data_3d);

         Geodetic_Vector_Data_2D*
         get_initialized_vd_2d (const Integer vector_size) const;

         void
         fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Level& level,
                                 const Met_Element met_element);

         void
         fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Key& key,
                         const Met_Element met_element);

         void
         fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Met_Element met_element);

         void
         fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Met_Element met_element);

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

         Geodetic_Vector_Data_3D*
         get_gvd_3d_ptr (const Met_Element met_element,
                         const Key& key) const;

         Geodetic_Vector_Data_2D*
         get_grib_data_ptr (const Grib& grib,
                            const Grib::Key& grib_key) const;

         void
         fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Met_Element met_element,
                         const Key& key,
                         const Level& level) const;

      public:

         Ecmwf (const string& description,
                const string& data_path,
                const string& search_string,
                const bool omega_as_w = false);

         ~Ecmwf ();

         void
         survey ();

         void
         clean_up ();

   };

   class Gfs3 : public Nwp,
                public map<Nwp::Key, Grib*>
   {

      private:

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

         const string
         data_path;

         Size_2D
         size_2d;

         Domain_2D
         domain_2d;

         Grib::Key
         get_grib_key (const Key& key,
                       const Met_Element met_element,
                       const Level& level) const;

         void
         set_grib_key (Grib::Key& grib_key,
                       const Dtime& base_time) const;

         void
         set_grib_key (Grib::Key& grib_key,
                       const Met_Element met_element,
                       const Integer forecast_hour) const;

         void
         set_grib_key (Grib::Key& grib_key,
                       const Met_Element met_element) const;

         void
         set_grib_key (Grib::Key& grib_key,
                       const Met_Element met_element,
                       const denise::Level& level) const;

         void
         initialize_3d_data (const Key& key);

         void
         load_3d_data (Nwp::Data_3D& data_3d);

         Geodetic_Vector_Data_2D*
         get_initialized_vd_2d (const Integer vector_size) const;

         void
         fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Level& level,
                                 const Met_Element met_element);

         void
         fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                    const Integer vector_index,
                                    const Key& key);

         void
         fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key);

         void
         fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Key& key,
                         const Met_Element met_element);

         void
         fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Met_Element met_element);

         void
         fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Met_Element met_element);

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

/*
         void
         fill_grib_data (Geodetic_Vector_Data_3D& gvd_3d,
                         const Integer vector_index,
                         const Met_Element met_element,
                         const Key& key) const;
*/
         Geodetic_Vector_Data_3D*
         get_gvd_3d_ptr (const Met_Element met_element,
                         const Key& key) const;

         Geodetic_Vector_Data_2D*
         get_grib_data_ptr (const Grib& grib,
                            const Grib::Key& grib_key) const;

         void
         fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Met_Element met_element,
                         const Key& key,
                         const Level& level) const;

      public:

         Gfs3 (const string& data_path);

         ~Gfs3 ();

         void
         survey ();

         void
         clear_3d_data ();

         void
         clean_up ();

         vector<Dtime>
         get_valid_time_vector () const;

         Key
         get_key (const Dtime& dtime) const;

         void
         acquire_base_time_forecast_hour (Dtime& base_time,
                                          Integer& forecast_hour,
                                          const Dtime& dtime) const;

   };

   class Gfs4 : public Nwp,
                public map<Nwp::Key, Grib2*>
   {

      private:

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

         const string
         data_path;

         Size_2D
         grib_size_2d;

         Domain_2D
         grib_domain_2d;

         Size_2D
         size_2d;

         Domain_2D
         domain_2d;

         Grib2::Key
         get_grib_key (const Key& key,
                       const Met_Element met_element,
                       const Level& level) const;

         void
         set_grib_key (Grib2::Key& grib_key,
                       const Dtime& base_time) const;

         void
         set_grib_key (Grib2::Key& grib_key,
                       const Met_Element met_element,
                       const Integer forecast_hour) const;

         void
         set_grib_key (Grib2::Key& grib_key,
                       const Met_Element met_element) const;

         void
         set_grib_key (Grib2::Key& grib_key,
                       const Met_Element met_element,
                       const denise::Level& level) const;

         void
         initialize_3d_data (const Key& key);

         void
         load_3d_data (Nwp::Data_3D& data_3d);

         Geodetic_Vector_Data_2D*
         get_initialized_vd_2d (const Integer vector_size) const;

         void
         fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Level& level,
                                 const Met_Element met_element);

         void
         fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                    const Integer vector_index,
                                    const Key& key);

         void
         fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key);

         void
         fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Key& key,
                         const Met_Element met_element);

         void
         fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Met_Element met_element);

         void
         fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Met_Element met_element);

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

/*
         void
         fill_grib_data (Geodetic_Vector_Data_3D& gvd_3d,
                         const Integer vector_index,
                         const Met_Element met_element,
                         const Key& key) const;
*/

         Geodetic_Vector_Data_3D*
         get_gvd_3d_ptr (const Met_Element met_element,
                         const Key& key) const;

         Geodetic_Vector_Data_2D*
         get_grib_data_ptr (const Grib2& grib,
                            const Grib2::Key& grib_key) const;

         void
         fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Met_Element met_element,
                         const Key& key,
                         const Level& level) const;

      public:

         Gfs4 (const string& data_path);

         ~Gfs4 ();

         void
         survey ();

         void
         clean_up ();

         void
         set_domain_2d (const Domain_2D& domain_2d);

         vector<Dtime>
         get_valid_time_vector () const;

         Key
         get_key (const Dtime& dtime) const;

         void
         acquire_base_time_forecast_hour (Dtime& base_time,
                                          Integer& forecast_hour,
                                          const Dtime& dtime) const;

   };

   class Gfs : public Nwp,
               public map<Nwp::Key, Grib2*>
   {

      private:

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

         const string
         data_path;

         Size_2D
         grib_size_2d;

         Domain_2D
         grib_domain_2d;

         Size_2D
         size_2d;

         Domain_2D
         domain_2d;

         Grib2::Key
         get_grib_key (const Key& key,
                       const Met_Element met_element,
                       const Level& level) const;

         void
         set_grib_key (Grib2::Key& grib_key,
                       const Dtime& base_time) const;

         void
         set_grib_key (Grib2::Key& grib_key,
                       const Met_Element met_element,
                       const Integer forecast_hour) const;

         void
         set_grib_key (Grib2::Key& grib_key,
                       const Met_Element met_element) const;

         void
         set_grib_key (Grib2::Key& grib_key,
                       const Met_Element met_element,
                       const denise::Level& level) const;

         void
         initialize_3d_data (const Key& key);

         void
         load_3d_data (Nwp::Data_3D& data_3d);

         Geodetic_Vector_Data_2D*
         get_initialized_vd_2d (const Integer vector_size) const;

         void
         fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Level& level,
                                 const Met_Element met_element);

         void
         fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                    const Integer vector_index,
                                    const Key& key);

         void
         fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key);

         void
         fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Key& key,
                         const Met_Element met_element);

         void
         fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Met_Element met_element);

         void
         fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Met_Element met_element);

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

/*
         void
         fill_grib_data (Geodetic_Vector_Data_3D& gvd_3d,
                         const Integer vector_index,
                         const Met_Element met_element,
                         const Key& key) const;
*/

         Geodetic_Vector_Data_3D*
         get_gvd_3d_ptr (const Met_Element met_element,
                         const Key& key) const;

         Geodetic_Vector_Data_2D*
         get_grib_data_ptr (const Grib2& grib,
                            const Grib2::Key& grib_key) const;

         void
         fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Met_Element met_element,
                         const Key& key,
                         const Level& level) const;

      public:

         Gfs (const string& data_path);

         ~Gfs ();

         void
         survey ();

         void
         clean_up ();

         void
         set_domain_2d (const Domain_2D& domain_2d);

         vector<Dtime>
         get_valid_time_vector () const;

         Key
         get_key (const Dtime& dtime) const;

         void
         acquire_base_time_forecast_hour (Dtime& base_time,
                                          Integer& forecast_hour,
                                          const Dtime& dtime) const;

   };

};

#endif /* DENISE_GRIB_H */ 
