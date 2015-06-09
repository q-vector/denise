//
// nwp.h
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

#ifndef DENISE_NWP_H
#define DENISE_NWP_H

#include <denise/gis.h>
#include <denise/geodesy.h>
#include <denise/thermo.h>

using namespace std;

#define TIME_TOLERANCE 0.1

namespace denise
{

   enum Nwp_Element
   {
      PRESSURE = 0,          P = 0,
      TEMPERATURE = 1,       T = 1,
      TEMPERATURE_CELCIUS,
      LAPSE,
      MIX_DOWN_TEMPERATURE,
      ZONAL_WIND = 5,        U = 5,
      MERIDIONAL_WIND = 6,   V = 6,
      WIND_SPEED,
      WIND_DIRECTION,
      VERTICAL_VELOCITY = 9, W = 9,
      STREAMLINE_WIND,
      NORMAL_WIND,
      GEOPOTENTIAL_HEIGHT = 12, Z = 12,
      DEW_POINT = 13,         TD = 13,
      MEAN_SEA_LEVEL_PRESSURE = 14, MSLP = 14,
      HIGH_CLOUD,
      MIDDLE_CLOUD,
      LOW_CLOUD,
      TOTAL_CLOUD,
      OMEGA,
      SNOW_LEVEL,
      FFDI,
      GFDI,
      CONTINUOUS_HAINES,
      PPT3,
      PPT6,
      PPTN,
      RAINFALL_STEP,
      RAINFALL_CUMULATIVE,
      RELATIVE_HUMIDITY = 29, RH = 29,
      SPECIFIC_HUMIDITY = 30, Q = 30,
      DEW_POINT_DEPRESSION,
      POTENTIAL_TEMPERATURE = 32, THETA = 32,
      THETA_E,
      THETA_W,
      MIXING_RATIO,
      MONTGOMERY,
      SLI,
      SHOWALTER,
      LI_700,
      LI_THUNDER,
      K_INDEX,
      TOTAL_TOTALS,
      CAPE,
      PRECIPITABLE_WATER,
      FOG_FRACTION,
      THICKNESS,
      POTENTIAL_VORTICITY,
      PV1_5_PRESSURE,
      ABSOLUTE_VORTICITY = 49, ZETA_A = 49,
      RELATIVE_VORTICITY = 50, ZETA = 50,
      SHEAR_VORTICITY,
      CURVATURE_VORTICITY,
      TEMPERATURE_ADVECTION,
      ADIABATIC_HEATING,
      LATENT_HEATING,
      DENSITY = 56, RHO = 56,
      P_THETA,
      P_RHO,
      VARIANCE_TEMPERATURE,
      VARIANCE_DEW_POINT,
      VARIANCE_GEOPOTENTIAL_HEIGHT,
      VARIANCE_ZONAL_WIND,
      VARIANCE_MERIDIONAL_WIND,
      VARIANCE_OMEGA
   };

   class Level_Tuple : public Tuple
   {

      public:

         const Level_Type
         type;

         //Level_Tuple ();

         Level_Tuple (const Level_Type type,
                      const string& str,
                      const string& delimiter = ":");

         void
         add (const string& str,
              const string& delimiter = ":",
              const bool clear_first = false);

         const Real
         get_next_up (const Real value) const;

         const Real
         get_next_down (const Real value) const;

   };

   class Level_Element
   {

      public:

         Level
         level;

         Nwp_Element
         nwp_element;

         Level_Element (const Level& level,
                        const Nwp_Element nwp_element);

   };

   class Nwp
   {

      public:

         class Key
         {

            public:

               Dtime
               base_time;

               Integer
               forecast_hour;

               Key () { }

               Key (const Dtime& base_time,
                    const Integer forecast_hour = 0);

               Key (const Key& key);

               bool
               operator == (const Key& key) const;

               bool
               operator > (const Key& key) const;

               bool
               operator < (const Key& key) const;

         };

         class Sounding : public denise::Sounding
         {

            protected:

               const Key&
               key;

            public:

               Sounding (const Key& key);

         };

         class Element_Vector : public vector<Nwp_Element>
         {

             private:

                map<Nwp_Element, Integer>
                reverse_map;

             public:

                Element_Vector (const vector<Nwp_Element>& nwp_element_vector);

                Integer
                get_index (const Nwp_Element nwp_element) const;

                Nwp_Element
                get_nwp_element (const Integer element_index) const;

         };

         // Assumes vertical coordinate is p
         class Data_3D : public map<Nwp_Element, Geodetic_Vector_Data_3D*>
         {

            private:

               bool
               available;

               const vector<Nwp_Element>
               nwp_element_vector;

            public:

               const Key
               key;

               Data_3D (const vector<Nwp_Element>& nwp_element_vector,
                        const Key& key);

               ~Data_3D ();

               void
               read (FILE* file,
                     const bool float_length = true);

               void
               read (ifstream& file,
                     const bool float_length = true);

               void
               write (FILE* file,
                      const bool float_length = true) const;

               void
               write (ofstream& file,
                      const bool float_length = true) const;

               void
               unload ();

               void
               unload (const Nwp_Element nwp_element);

               void
               set_gvd_3d_ptr (const Nwp_Element nwp_element,
                               Geodetic_Vector_Data_3D* gvd_3d_ptr);

               const Geodetic_Vector_Data_3D&
               get_gvd_3d (const Nwp_Element nwp_element) const;

               Geodetic_Vector_Data_3D&
               get_gvd_3d (const Nwp_Element nwp_element);

               virtual const Tuple&
               get_tuple_p (const Nwp_Element nwp_element) const;

               Real
               get_p (const Nwp_Element nwp_element,
                      const Integer k) const;

               Lat_Long
               get_lat_long (const Nwp_Element nwp_element,
                             const Integer i,
                             const Integer j) const;

               bool
               is_available () const;

               void
               set_available ();

               void
               initialize (const Real datum);

               virtual Real
               get_p_from_element (const Nwp_Element nwp_element,
                                   const Real latitude,
                                   const Real longitude,
                                   const Real element_value) const;

               virtual Real
               evaluate (const Nwp_Element element,
                         const Real p,
                         const Real latitude,
                         const Real longitude,
                         const Evaluate_Op evaluate_op = VALUE) const;

               virtual Real
               evaluate (const Nwp_Element element,
                         const Real p,
                         const Lat_Long& lat_long,
                         const Evaluate_Op evaluate_op = VALUE) const;

               Real
               get_li_thunder (const Real p,
                               const Real latitude,
                               const Real longitude,
                               const Real thunder_p,
                               const Real thunder_t) const;

               Real
               get_li_thunder (const Real p,
                               const Lat_Long& lat_long,
                               const Real thunder_p,
                               const Real thunder_t) const;

         };

         class Time_Cross : public map<Nwp_Element, Scalar_Data_2D*>
         {

            private:

               Scalar_Data_1D*
               terrain_profile_ptr;

               Scalar_Data_1D*
               rainfall_profile_ptr;

            public:

               Time_Cross ();

               ~Time_Cross ();

               void
               set_profile_ptrs (Scalar_Data_1D* terrain_profile_ptr,
                                 Scalar_Data_1D* rainfall_profile_ptr);

               void
               insert_nwp_element_if_needed (const Nwp_Element nwp_element,
                                             const Tuple& tuple_t,
                                             const Tuple& tuple_p);

               const Scalar_Data_1D&
               get_terrain_profile () const;

               const Scalar_Data_1D&
               get_rainfall_profile () const;

               const Scalar_Data_2D&
               get_sd_2d (const Nwp_Element nwp_element) const;

               Scalar_Data_2D&
               get_sd_2d (const Nwp_Element nwp_element);

               virtual const Tuple&
               get_tuple_p (const Nwp_Element nwp_element) const;

         };

         class Cross_Section : public map<Nwp_Element, Scalar_Data_2D*>
         {

            private:

               Scalar_Data_1D*
               terrain_profile_ptr;

               Scalar_Data_1D*
               rainfall_profile_ptr;

            public:

               Cross_Section ();

               ~Cross_Section ();

               void
               set_profile_ptrs (Scalar_Data_1D* terrain_profile_ptr,
                                 Scalar_Data_1D* rainfall_profile_ptr);

               void
               insert_nwp_element_if_needed (const Nwp_Element nwp_element,
                                             const Tuple& tuple_x,
                                             const Tuple& tuple_p);

               const Scalar_Data_1D&
               get_terrain_profile () const;

               const Scalar_Data_1D&
               get_rainfall_profile () const;

               const Scalar_Data_2D&
               get_sd_2d (const Nwp_Element nwp_element) const;

               Scalar_Data_2D&
               get_sd_2d (const Nwp_Element nwp_element);

               virtual const Tuple&
               get_tuple_p (const Nwp_Element nwp_element) const;

         };

         class Key_Multimap : public multimap<Dtime, Integer>
         {

            private:

               set<Key>
               key_set;

            public:

               void
               clear ();

               void
               add (const Key& key);

               bool
               is_no_match (const Key& key) const;

               bool
               is_first_step (const Key& key) const;

               Key
               get_previous_key (const Key& key) const;

               set<Dtime>
               get_base_time_set () const;

               set<Dtime>
               get_valid_time_set (const Dtime& base_time) const;

               Key
               get_key (const Dtime& dtime) const;

               Key
               get_key (const Dtime& dtime,
                        const Dtime& base_time) const;

         };

      protected:

         typedef map<Key, Data_3D*>
         Data_3d_Ptr_Map;

         vector<Nwp_Element>
         nwp_element_vector;

         Tuple
         tuple_p;

         const string
         path;

         const string
         description;

         string
         status;

         Key_Multimap
         key_multimap;

         Data_3d_Ptr_Map
         data_3d_ptr_map;

         virtual Geodetic_Vector_Data_2D*
         get_initialized_vd_2d (const Integer vector_size) const = 0;

         virtual void
         initialize_3d_data (const Key& key) = 0;

         virtual void
         load_3d_data (Data_3D& data_3d) = 0;

         virtual void
         fill_lapse_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Level& level);

         virtual void
         fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Key& key,
                         const Nwp_Element nwp_element);

         virtual void
         fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Key& key,
                         const Integer hours = -1);

         void
         fill_p_potential_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                          const Integer vector_index,
                                          const Key& key,
                                          const Real p);

         void
         fill_theta_potential_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                              const Integer vector_index,
                                              const Key& key,
                                              const Real theta);

         void
         fill_potential_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                        const Integer vector_index,
                                        const Key& key,
                                        const Level& level);

         void
         fill_p_temperature_advection_data (Geodetic_Vector_Data_2D& gvd_2d,
                                            const Integer vector_index,
                                            const Key& key,
                                            const Real p);

         void
         fill_theta_temperature_advection_data (Geodetic_Vector_Data_2D& gvd_2d,
                                                const Integer vector_index,
                                                const Key& key,
                                                const Real theta);

         void
         fill_temperature_advection_data (Geodetic_Vector_Data_2D& gvd_2d,
                                          const Integer vector_index,
                                          const Key& key,
                                          const Level& level);

         void
         fill_p_adiabatic_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                        const Integer vector_index,
                                        const Key& key,
                                        const Real p);

         void
         fill_theta_adiabatic_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                            const Integer vector_index,
                                            const Key& key,
                                            const Real theta);

         void
         fill_adiabatic_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                      const Integer vector_index,
                                      const Key& key,
                                      const Level& level);

         void
         fill_p_latent_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                     const Integer vector_index,
                                     const Key& key,
                                     const Real p);

         void
         fill_theta_latent_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                         const Integer vector_index,
                                         const Key& key,
                                         const Real theta);

         void
         fill_latent_heating_data (Geodetic_Vector_Data_2D& gvd_2d,
                                   const Integer vector_index,
                                   const Key& key,
                                   const Level& level);

         virtual void
         fill_pv_p_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Key& key,
                         const Real pv_threshold = 1.5e-6);

         void
         fill_absolute_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                       const Integer vector_index,
                                       const Key& key,
                                       const Level& level);

         void
         fill_shear_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                    const Integer vector_index,
                                    const Key& key,
                                    const Level& level);

         void
         fill_curvature_vorticity_data (Geodetic_Vector_Data_2D& gvd_2d,
                                        const Integer vector_index,
                                        const Key& key,
                                        const Level& level);

         void
         fill_wind_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index_u,
                         const Integer vector_index_v,
                         const Key& key,
                         const Level& level);

         void
         fill_mix_down_temperature_data (Geodetic_Vector_Data_2D& gvd_2d,
                                         const Integer vector_index,
                                         const Key& key,
                                         const Level& level);

         virtual void
         fill_snow_probability_data (Geodetic_Vector_Data_2D& gvd_2d,
                                     const Integer vector_index,
                                     const Key& key);

         virtual void
         fill_snow_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Real probability);

         virtual void
         fill_fdi_data (Geodetic_Vector_Data_2D& gvd_2d,
                        const Integer vector_index,
                        const Key& key,
                        const Level& level,
                        const Nwp_Element nwp_element);

         virtual void
         fill_continuous_haines_data (Geodetic_Vector_Data_2D& gvd_2d,
                                      const Integer vector_index,
                                      const Key& key);

         virtual void
         fill_total_totals_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key);

         virtual void
         fill_k_index_data (Geodetic_Vector_Data_2D& gvd_2d,
                            const Integer vector_index,
                            const Key& key);

         virtual void
         fill_lifted_index_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Real bottom_p,
                                 const Real top_p);

         virtual void
         fill_surface_lifted_index_data (Geodetic_Vector_Data_2D& gvd_2d,
                                         const Integer vector_index,
                                         const Key& key,
                                         const Real top_p);

         virtual void
         fill_li_thunder_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Level& level,
                               const Real thunder_t = -20 + K);

         virtual void
         fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Level& level,
                                 const Nwp_Element nwp_element);

         virtual void
         fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Nwp_Element nwp_element);

         virtual void
         fill_pressure_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                   const Integer vector_index,
                                   const Key& key,
                                   const Real p,
                                   const Nwp_Element nwp_element);

         virtual void
         fill_theta_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Real theta,
                                const Nwp_Element nwp_element);

         virtual void
         fill_sigma_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Real sigma,
                                const Nwp_Element nwp_element,
                                const Geodetic_Vector_Data_2D& surface_p_data);

         virtual void
         fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Nwp_Element nwp_element);

         virtual void
         fill_50m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Nwp_Element nwp_element);

         virtual void
         fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Nwp_Element nwp_element);

         virtual void
         fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                        const Integer vector_index,
                        const Key& key,
                        const Nwp_Element nwp_element);

         virtual void
         fill_nil_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Nwp_Element nwp_element);

         virtual void
         fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                  const Integer vector_index,
                                  const Key& key,
                                  const Nwp_Element nwp_element);

         virtual void
         fill_data (Geodetic_Vector_Data_2D& gvd_2d,
                    const Integer vector_index,
                    const Key& key,
                    const Level& level,
                    const Nwp_Element nwp_element);

      public:

         Nwp (const string& description,
              const string& path);

         ~Nwp ();

         virtual void
         survey () = 0;

         virtual void
         clean_up () = 0;

         void
         clear_data_3d_ptr_map ();

         virtual void
         set_domain_2d (const Domain_2D& domain_2d);

         const string&
         get_description () const;

         const string&
         get_status () const;

         virtual set<Dtime>
         get_valid_time_set (const Dtime& base_time) const;

         virtual Tuple
         get_valid_t_tuple (const Dtime& base_time,
                            const Dtime& start_time,
                            const Dtime& end_time) const;

         virtual set<Dtime>
         get_base_time_set () const;

         virtual Key
         get_key (const Dtime& dtime,
                  const Dtime& base_time) const;

         Nwp::Data_3D&
         get_3d_data (const Key& key);

         virtual Geodetic_Vector_Data_2D*
         get_ts_steering_data_ptr (const Key& key);

         virtual Geodetic_Vector_Data_2D*
         get_steering_data_ptr (const Key& key,
                                const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_vertical_shear_data_ptr (const Key& key,
                                      const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_min_li_thunder_data_ptr (const Key& key,
                                      const Real thunder_t = -20 + K);

         virtual Geodetic_Vector_Data_2D*
         get_freezing_level_data_ptr (const Key& key);

         virtual Geodetic_Vector_Data_2D*
         get_temperature_p_data_ptr (const Key& key,
                                     const Real temperature);

         virtual Geodetic_Vector_Data_2D*
         get_p_data_ptr (const Key& key,
                         const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_next_p_up_data_ptr (const Key& key,
                                 const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_cloud_data_ptr (const Key& key);

         virtual Geodetic_Vector_Data_2D*
         get_cloud_base_data_ptr (const Key& key,
                                  const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_cloud_top_temp_data_ptr (const Key& key);

         virtual Geodetic_Vector_Data_2D*
         get_ageostrophic_wind_data_ptr (const Key& key,
                                         const Real p);

         virtual Geodetic_Vector_Data_2D*
         get_surface_ageostrophic_wind_data_ptr (const Key& key,
                                                 const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_ageostrophic_wind_data_ptr (const Key& key,
                                         const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_geostrophic_wind_data_ptr (const Key& key,
                                        const Real p);

         virtual Geodetic_Vector_Data_2D*
         get_surface_geostrophic_wind_data_ptr (const Key& key);

         virtual Geodetic_Vector_Data_2D*
         get_geostrophic_wind_data_ptr (const Key& key,
                                        const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_wind_data_ptr (const Key& key,
                            const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_horizontal_shear_data_ptr (const Key& key,
                                        const Level& level,
                                        const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_omega_data_ptr (const Key& key,
                             const Level& level,
                             const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_temperature_24hr_data_ptr (const Key& key,
                                        const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_mix_down_temperature_data_ptr (const Key& key,
                                            const Level& level,
                                            const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_stratus_data_ptr (const Key& key,
                               const Level& level,
                               const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_lapse_data_ptr (const Key& key,
                             const Level& level,
                             const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_temperature_data_ptr (const Key& key,
                                   const Level& level,
                                   const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_dew_point_data_ptr (const Key& key,
                                 const Level& level,
                                 const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_rh_data_ptr (const Key& key,
                          const Level& level,
                          const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_dew_point_depression_data_ptr (const Key& key,
                                            const Level& level,
                                            const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_theta_e_data_ptr (const Key& key,
                               const Level& level,
                               const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_theta_w_data_ptr (const Key& key,
                               const Level& level,
                               const bool with_wind = false);

         virtual Geodetic_Vector_Data_2D*
         get_q_vector_data_ptr (const Key& key,
                                const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_temperature_advection_data_ptr (const Key& key,
                                             const Level& level);

         virtual Geodetic_Vector_Data_2D*
         get_gd_sea_wave_data_ptr (const Key& key,
                                   const Level& level,
                                   const Fetch& fetch,
                                   const Real max_fetch = 150e3);

         virtual Geodetic_Vector_Data_2D*
         get_step_rainfall_data_ptr (const Key& key);

         virtual Geodetic_Vector_Data_2D*
         get_rainfall_data_ptr (const Key& key,
                                const Integer hours = -1);

         virtual Geodetic_Vector_Data_2D*
         get_isallobar_data_ptr (const Key& key);

         virtual Geodetic_Vector_Data_2D*
         get_sutcliffe_data_ptr (const Key& key);

         virtual Geodetic_Vector_Data_2D*
         get_pv_p_data_ptr (const Key& key,
                            const Real pv_threshold = 1.5e-6);

         virtual Geodetic_Vector_Data_2D*
         get_theta_level_data_ptr (const Key& key,
                                   const Real theta,
                                   const vector<Nwp_Element> element_vector,
                                   const bool with_pv);

         virtual Geodetic_Vector_Data_2D*
         get_sigma_level_data_ptr (const Key& key,
                                   const Real sigma,
                                   const vector<Nwp_Element> element_vector);

         virtual Geodetic_Vector_Data_2D*
         get_pressure_level_data_ptr (const Key& key,
                                      const Real p,
                                      const vector<Nwp_Element> element_vector);

         virtual Geodetic_Vector_Data_2D*
         get_data_ptr (const Key& key,
                       const Level& level,
                       const Nwp_Element nwp_element);

         virtual Geodetic_Vector_Data_2D*
         get_data_ptr (const Key& key,
                       const Level& level,
                       const vector<Nwp_Element>& nwp_element_vector);

         virtual Scalar_Data_1D*
         get_terrain_profile_ptr (const Key& key,
                                  const Multi_Journey& multi_journey);

         virtual Scalar_Data_1D*
         get_rainfall_profile_ptr (const Key& key,
                                   const Multi_Journey& multi_journey);

         virtual Cross_Section*
         get_cross_section_ptr (const Key& key,
                                const Multi_Journey& multi_journey,
                                const Nwp_Element nwp_element,
                                const bool with_wind);

         virtual Cross_Section*
         get_cross_section_ptr (const Key& key,
                                const Multi_Journey& multi_journey,
                                const vector<Nwp_Element>& element_vector,
                                const bool with_wind);

         virtual Nwp::Sounding*
         get_sounding_ptr (const Lat_Long& lat_long,
                           const Key& key);

         virtual Tuple
         get_point_tuple (const Lat_Long& lat_long,
                          const Key& key);

         virtual Vector_Data_1D*
         get_time_series_ptr (const Lat_Long& lat_long,
                              const Dtime& base_time,
                              const Dtime& start_time,
                              const Dtime& end_time,
                              const Level_Element level_element);

         virtual Vector_Data_1D*
         get_time_series_ptr (const Lat_Long& lat_long,
                              const Dtime& base_time,
                              const Dtime& start_time,
                              const Dtime& end_time,
                              const vector<Level_Element>& level_element_vector);

         virtual Time_Cross*
         get_time_cross_data_ptr (const Lat_Long& lat_long,
                                  const Dtime& base_time,
                                  const Dtime& start_time,
                                  const Dtime& end_time,
                                  const vector<Nwp_Element>& element_vector);

   };

   class Nwp_Exception : public Exception
   {

      public:

         Nwp_Exception (const string& str = string (""));

   };

};

#endif /* DENISE_NWP_H */ 
