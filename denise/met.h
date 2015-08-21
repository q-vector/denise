//
// met.h
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

#ifndef DENISE_MET_H
#define DENISE_MET_H

#include <cmath>
#include <cstdio>
#include <iostream>
#include <cairomm/context.h>
#include <denise/analysis.h>
#include <denise/basics.h>
#include <denise/geodesy.h>
#include <denise/graphics.h>
#include <denise/geometry.h>

using namespace std;

namespace denise
{

   enum Met_Element
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
      THETA_V,
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
      ABSOLUTE_VORTICITY = 50, ZETA_A = 50,
      RELATIVE_VORTICITY = 51, ZETA = 51,
      SHEAR_VORTICITY,
      CURVATURE_VORTICITY,
      TEMPERATURE_ADVECTION,
      ADIABATIC_HEATING,
      LATENT_HEATING,
      DENSITY = 57, RHO = 57,
      PRECIP_RATE,
      BRUNT_VAISALA,
      SCORER_PARAMETER,
      Q_H_ADVECTION,
      Q_V_ADVECTION,
      Q_S_ADVECTION,
      Q_N_ADVECTION,
      P_THETA,
      P_RHO,
      VARIANCE_TEMPERATURE,
      VARIANCE_DEW_POINT,
      VARIANCE_GEOPOTENTIAL_HEIGHT,
      VARIANCE_ZONAL_WIND,
      VARIANCE_MERIDIONAL_WIND,
      VARIANCE_OMEGA
   };

   /// Tropical Cyclone Category
   enum Tc_Category
   {
      UNKNOWN_TC_CATEGORY   = -2,
      ALL_TC_CATEGORIES     = -1,
      TROPICAL_DEPRESSION   = 0, TC_TD  = 0,
      TROPICAL_STORM        = 1, TC_TS  = 1,
      SEVERE_TROPICAL_STORM = 2, TC_STS = 2,
      TYPHOON               = 3, TC_T   = 3, TC_TY = 3,
      LOW                   = 4,
      NUMBER_OF_TC_CATEGORIES = 5
   };

   /// Returns TC_Category given wind speed
   ///
   /// \param max_wind Wind speed
   Tc_Category
   get_tc_category (const Real max_wind);

   /// Returns Tc_Category given a char* represention of a Tc_Category
   ///
   /// \param tc_category_string character string
   Tc_Category
   get_tc_category (const Dstring& tc_category_string);

   /// Returns maximal intensity of a Tc_Category
   ///
   /// \param tc_category TC Category
   Real
   get_maximal_intensity (const Tc_Category tc_category);

   /// Returns minimal intensity of a Tc_Category
   ///
   /// \param tc_category TC Category
   Real
   get_minimal_intensity (const Tc_Category tc_category);

   /// Returns a string representation of a Tc_Category
   ///
   /// \param tc_category TC Category
   Dstring
   get_tc_category_string (const Tc_Category tc_category);

   class Tc_Symbol : public Symbol
   {

      private:

         void
         init_north (const Point_2D& point,
                     const Point_2D& left,
                     const Point_2D& right,
                     const Point_2D& up,
                     const Point_2D& down,
                     const Real r,
                     const Real r_o,
                     const Real r_i,
                     const Real s_o,
                     const Real s_i,
                     const Real w,
                     const Real w_2,
                     const Real theta,
                     const Real phi,
                     const Real lambda);

         void
         init_south (const Point_2D& point,
                     const Point_2D& left,
                     const Point_2D& right,
                     const Point_2D& up,
                     const Point_2D& down,
                     const Real r,
                     const Real r_o,
                     const Real r_i,
                     const Real s_o,
                     const Real s_i,
                     const Real w,
                     const Real w_2,
                     const Real theta,
                     const Real phi,
                     const Real lambda);

      public:

         Tc_Symbol (const bool southern_hemisphere,
                    const Real size,
                    const Real width = GSL_NAN);

   };

   class Octa : public Symbol
   {

      public:

         enum Number
         {
            OCTA_ZERO,
            OCTA_ONE,
            OCTA_TWO,
            OCTA_THREE,
            OCTA_FOUR,
            OCTA_FIVE,
            OCTA_SIX,
            OCTA_SEVEN,
            OCTA_EIGHT,
            OCTA_OBSCURED,
            OCTA_GARBLED
         };

         Octa (const Octa::Number octa_number,
               const Real size,
               const Real width = 1);

   };

   /// Wind with m/s as unit
   class Wind
   {

      public:

         /// Wind Category
         enum Category
         {
            CALM,
            LIGHT,
            MODERATE,
            FRESH,
            STRONG,
            GALE,
            STORM,
            HURRICANE
         };

         /// u-component
         Real
         u;

         /// v-component
         Real
         v;

         /// Constructor
         Wind (const Real u = 0,
               const Real v = 0);

         /// Copy Constuctor
         Wind (const Wind& wind);

         static Wind
         direction_speed (const Real direction,
                          const Real speed);

         void
         set_from_direction_speed (const Real direction,
                                   const Real speed);

         /// Normalizes the Wind to magnitude
         void
         normalize (const Real magnitude = 1);

         /// veers the wind by degree degrees
         void
         veer (const Real degree);

         /// backs the wind by degree degrees
         void
         back (const Real degree);

         /// Returns a Wind to normalized magnitude
         Wind
         get_normalized_wind (const Real magnitude = 1) const;

         /// Returns a Wind veered by degree
         Wind
         get_veered_wind (const Real degree) const;

         /// Returns a Wind backed by degree
         Wind
         get_backed_wind (const Real degree) const;

         /// Returns direction of Wind
         Real
         get_direction () const;

         /// Returns string_representation of direction of Wind
         Dstring
         get_direction_string (const Integer n = 16,
                               const Dstring& format = "") const;

         /// Returns speed of Wind
         Real
         get_speed () const;

         bool
         is_naw () const;

         Dstring
         get_string (const Real speed_multiplier = 3.6/1.852,
                     const Dstring& format = "%03.0f/%02.0f") const;

         static Category
         get_category (const Real speed);

         static Dstring
         get_category_string (const Category category);

         Wind
         operator + (const Wind& wind) const;

         void
         operator += (const Wind& wind);

         Wind
         operator - (const Wind& wind) const;

         void
         operator -= (const Wind& wind);

         Wind
         operator * (const Real a) const;

         void
         operator *= (const Real a);

         Wind
         operator / (const Real a) const;

         void
         operator /= (const Real a);

   };

   class Wind_Barb : public Symbol
   {

      private:

         const bool
         north;

         void
         analysis (Integer& pennants,
                   Integer& barbs,
                   Real& residue,
                   Real speed);

         void
         add_pennant (const Real theta,
                      const Point_2D& point,
                      const Real length);

         void
         add_barb (const Real theta,
                   const Point_2D& point,
                   const Real length,
                   const Real width,
                   const bool invisible_barb);

         void
         barb (const Wind& wind,
               const Real size,
               const Real width);

      public:

         Wind_Barb (const Wind& wind,
                    const Real size,
                    const bool north,
                    const Real width = gsl_nan (),
                    const Real calm_threshold = 1);

   };


   class Wind_Transform : public Polar_Transform_2D
   {

      public:

         Wind_Transform ();

         Wind_Transform (const Point_2D& origin,
                         const Real scale);

         void
         transform (Real& x,
                    Real& y,
                    const Real direction, 
                    const Real speed) const;

         void
         reverse (Real& direction,
                  Real& speed,
                  const Real x, 
                  const Real y) const;

   };

   class Wind_Rose_Threshold
   {

      public:

         Real
         value;

         Dstring
         label_str;

         Wind_Rose_Threshold (const Real threshold = GSL_NAN,
                              const Dstring& label_str = "");

   };

   class Wind_Rose_Record
   {

      public:

         Wind_Rose_Threshold
         threshold;

         Real
         percentage;

         Wind_Rose_Record (const Real threshold,
                           const Real percentage,
                           const Dstring& label_str = "");

         Wind_Rose_Record (const Wind_Rose_Threshold& threshold,
                           const Real percentage);

   };

   class Wind_Rose_Arm : public vector<Wind_Rose_Record>
   {

      public:

         void
         render (const RefPtr<Context> cr,
                 const Color& pen_color,
                 const Color_Chooser& color_chooser,
                 const Point_2D& point,
                 const Real direction,
                 const Real percentage_size,
                 const Real delta_width,
                 const bool show_label) const;

   };

   /// Wind_Rose
   class Wind_Rose
   {

      protected:

         Dstring
         unit_string;

         Real
         multiplier;

         Integer
         total_count;

         Integer
         calm_count;

         Real
         delta_direction;

         Integer*
         count_data;

         Integer
         number_of_directions;

         vector<Wind_Rose_Threshold>
         threshold_vector;

         void
         init (const Integer number_of_directions,
               const vector<Wind_Rose_Threshold>& threshold_vector,
               const Dstring& unit_string,
               const Real multiplier);

         void
         init (const Integer number_of_directions,
               const Tuple& threshold_tuple,
               const Dstring& unit_string,
               const Real multiplier);

         void
         init (const Dstring& unit_string,
               const Real multiplier);

         static Real
         get_multiplier (const Dstring& unit_string);

      public:

         Wind_Rose (const Integer number_of_directions,
                    const vector<Wind_Rose_Threshold>& threshold_vector,
                    const Dstring& unit_string = "kt",
                    const Real multiplier = GSL_NAN);

         /// Constructor
         ///
         /// \param number_of_directions number_of_directions of wind_rose
         /// \param threshold_vector     vector holding speed thresholds
         Wind_Rose (const Integer number_of_directions,
                    const Tuple& threshold_tuple,
                    const Dstring& unit_string = "kt",
                    const Real multiplier = GSL_NAN);

         /// Constructor
         ///
         /// \param number_of_directions    number_of_directions of wind_rose
         /// \param threshold_vector_string string representing threshold vector
         Wind_Rose (const Integer number_of_directions,
                    const Dstring& threshold_vector_string,
                    const Dstring& unit_string = "kt",
                    const Real multiplier = GSL_NAN);

         ~Wind_Rose ();

         const Integer&
         get_number_of_directions () const;

         const vector<Wind_Rose_Threshold>&
         get_threshold_vector () const;

         const Wind_Rose_Threshold&
         get_max_threshold () const;

         Tuple
         get_threshold_tuple () const;

         /// Returns total count
         const Integer&
         get_total_count () const;

         /// Returns count
         Integer&
         get_count (const Integer direction_index,
                    const Integer speed_index);

         /// Returns count
         const Integer&
         get_count (const Integer direction_index,
                    const Integer speed_index) const;

         /// Returns percentage at or stronger at this speed index
         Integer
         get_count_stronger (const Integer direction_index,
                             const Integer speed_index) const;

         /// Returns count
         Integer
         get_count_d (const Integer direction_index) const;

         /// Returns percentage
         Real
         get_percentage (const Integer direction_index,
                         const Integer speed_index) const;

         /// Returns percentage at or stronger at this speed index
         Real
         get_percentage_stronger (const Integer direction_index,
                                  const Integer speed_index) const;

         /// Returns percentage
         Real
         get_percentage_d (const Integer direction_index) const;

         /// Returns max percentage
         Real
         get_max_percentage () const;

         /// Returns calm count
         Integer
         get_calm_count () const;

         /// Returns calm percentage
         Real
         get_calm_percentage () const;

         Index_2D
         get_index (const Wind& wind) const;

         void
         clear ();

         /// Add one wind observation.
         virtual void
         add_wind (const Wind& wind);

         /// Add a vector of wind observations
         void
         add_wind (const vector<Wind>& wind_vector);

         Wind_Rose_Arm
         get_wind_rose_arm (const Integer direction_index) const;

         Wind_Rose_Arm
         get_unit_wind_rose_arm () const;

         void
         render (const RefPtr<Context>& cr,
                 const Point_2D& point_2d,
                 const Color& pen_color,
                 const Color& calm_color,
                 const Color& ring_color,
                 const Color_Chooser& color_chooser,
                 const Real calm_ring_size,
                 const Real percentage_size,
                 const Real delta_width,
                 const Tuple& ring_tuple = Tuple ()) const;

         void
         render (const RefPtr<Context>& cr,
                 const Transform_2D& transform_2d,
                 const Point_2D& point,
                 const Color& pen_color,
                 const Color& calm_color,
                 const Color& ring_color,
                 const Color_Chooser& color_chooser,
                 const Real calm_ring_size,
                 const Real percentage_size,
                 const Real delta_width,
                 const Tuple& ring_tuple = Tuple ()) const;

   };

   class Wind_Disc : public Wind_Rose
   {

      private:

         class Transform : public Transform_2D
         {

            private:

               const Wind_Disc&
               wind_disc;

               Point_2D
               origin;

               Real
               max_speed;

               Real
               max_radius;

               Real
               scale;

               Real
               calm_radius;

               Real
               label_height;

            public:

               Transform (const Wind_Disc& wind_disc,
                          const Point_2D& origin,
                          const Real max_speed,
                          const Real max_radius,
                          const Real calm_radius,
                          const Real label_height);

               bool
               out_of_domain (const Real x,
                              const Real y) const;

               void
               transform (Real& x,
                          Real& y,
                          const Real direction, 
                          const Real speed) const;

               void
               reverse (Real& direction,
                        Real& speed,
                        const Real x, 
                        const Real y) const;

               const Point_2D&
               get_origin () const;

               const Real
               get_max_speed () const;

               const Real
               get_max_radius () const;

               const Real
               get_calm_radius () const;

               const Real
               get_label_height () const;

               Real
               get_radius (const Real speed) const;

               Real
               get_speed (const Real radius) const;

         };

         Wind_Disc::Transform*
         transform_ptr;

         vector<Wind>
         wind_vector;

         const Tuple
         speed_label_tuple;
        
         const Color
         major_color;
        
         const Color
         middle_color;
        
         const Color
         minor_color;
        
         const Color
         shade_color;
        
         const Real
         font_size;
        
         const Real
         scatter_ring_size;
        
         void
         render_background (const RefPtr<Context> cr) const;

         void
         render_directions (const RefPtr<Context> cr) const;

         void
         render_direction_labels (const RefPtr<Context> cr) const;

         void
         render_directions_even (const RefPtr<Context> cr) const;

         void
         render_directions_odd (const RefPtr<Context> cr) const;

         Polygon*
         render_ring_label (const RefPtr<Context> cr) const;

         void
         render_mesh (const RefPtr<Context> cr) const;

         void
         render_mesh (const RefPtr<Context> cr,
                      const Color& color,
                      const Tuple& speed_tuple,
                      const Real line_width) const;

         void
         render_scatter_plot (const RefPtr<Context> cr,
                              const Real hue,
                              const Real dir_scatter = 0) const;

         void
         render_scatter_ring (const RefPtr<Context> cr,
                              const Real radius,
                              const Real line_width) const;

         void
         render_percentages (const RefPtr<Context> cr) const;

         void
         render_percentage (const RefPtr<Context> cr,
                            const Point_2D& point,
                            const Real percentage) const;

         void
         render_percentage_d (const RefPtr<Context> cr,
                              const Real hue) const;

      public:

         Wind_Disc (const Integer number_of_directions,
                    const Tuple& threshold_tuple,
                    const Point_2D& origin,
                    const Real max_radius,
                    const Tuple& speed_label_tuple = Tuple (),
                    const Color& major_color = Color (0.0, 0.0, 0.0, 0.8),
                    const Color& middle_color = Color (0.2, 0.2, 0.2, 0.5),
                    const Color& minor_color = Color (0.4, 0.4, 0.4, 0.2),
                    const Color& shade_color = Color (0.6, 0.6, 0.6, 0.2),
                    const Real max_speed = GSL_NAN,
                    const Real font_size = 12,
                    const Real scatter_ring_size = 10,
                    const Real calm_radius = 40,
                    const Real label_height = 24,
                    const Dstring& unit_string = "kt",
                    const Real multiplier = GSL_NAN);

         ~Wind_Disc ();

         void
         set_position (const Point_2D& origin,
                       const Real max_radius);

         void
         clear ();

         void
         add_wind (const Wind& wind);

         Wind
         get_wind (const Point_2D& point) const;

         void
         render (const RefPtr<Context> cr,
                 const Real hue,
                 const bool outline,
                 const Real dir_scatter = 0) const;

         void
         render_bg (const RefPtr<Context> cr) const;

         void
         render_index (const RefPtr<Context> cr,
                       const Index_2D& index) const;

         const Point_2D&
         get_origin () const;

         const Real
         get_max_radius () const;

   };

   /// Axisymmetric vortex with constant inflow angle
   /// 
   /// This is an abstract base class.
   class Axisymmetric_Vortex
   {

      protected:

         const Point_2D
         center;

         const Real
         inflow;

         virtual Real
         get_V (const Real r) const = 0;

      public:

         Axisymmetric_Vortex (const Real inflow = 0);

         Axisymmetric_Vortex (const Point_2D& center,
                              const Real inflow = 0);

         /// Returns the center of the vortex
         const Point_2D&
         get_center () const;

         /// Returns the inflow angle in degrees
         const Real&
         get_inflow () const;

         /// Returns the wind at point
         Wind
         get_wind (const Real x,
                   const Real y) const;

         /// Returns the wind at point
         Wind
         get_wind (const Point_2D& point) const;

         /// Returns the radial wind at point wrt to radar position
         Real
         get_radial_wind (const Point_2D& radar,
                          const Real x,
                          const Real y) const;

         /// Returns the radial wind at point wrt to radar position
         Real
         get_radial_wind (const Point_2D& radar,
                          const Point_2D& point) const;

         /// Returns the isodop when viewed from
         ///        a doppler radar at radar_center
         ///
         /// The isodop is the locus of zero tangential wind and
         /// can be proved to be a Circle.
         Circle
         get_isodop (const Point_2D& radar_center) const;

   };

   /// Axisymmetric Vortex with wind speed decays exponentially from center
   class Exponential_Vortex : public Axisymmetric_Vortex
   {

      private:

         const Real
         max_wind;

         const Real
         decay;

      public:

         /// Constructor
         ///
         /// Center is at (0, 0)
         ///
         /// \param max_wind Wind speed at center
         /// \param decay    Exponential decay factor
         /// \param inflow   Inflow angle in degrees.
         Exponential_Vortex (const Real max_wind,
                             const Real decay,
                             const Real inflow = 0);

         /// Constructor
         ///
         /// \param center   Vortex center
         /// \param max_wind Wind speed at center
         /// \param decay    Exponential decay factor
         /// \param inflow   Inflow angle in degrees.
         Exponential_Vortex (const Point_2D& center,
                             const Real max_wind,
                             const Real decay,
                             const Real inflow = 0);

         /// Returns wind speed at radius r
         Real
         get_V (const Real r) const;

         /// Returns wind speed at center
         const Real&
         get_max_wind () const;

         /// Returns exponential decay factor
         const Real&
         get_decay () const;

   };

   /// Axisymmetric Vortex with Rankine wind speed profile
   class Rankine_Vortex : public Axisymmetric_Vortex
   {

      private:

         const Real
         max_radius;

         const Real
         k_in;

         const Real
         k_out;

      public:

         /// Constructor
         ///
         /// Center is at (0, 0)
         ///
         /// \param max_wind   Wind speed at radius max_radius
         /// \param max_radius Radius with max_wind
         /// \param inflow     Inflow angle in degrees.
         Rankine_Vortex (const Real max_wind,
                         const Real max_radius,
                         const Real inflow = 0);

         /// Constructor
         ///
         /// \param center     Vortex center
         /// \param max_wind   Wind speed at radius max_radius
         /// \param max_radius Radius with max_wind
         /// \param inflow     Inflow angle in degrees.
         Rankine_Vortex (const Point_2D& center,
                         const Real max_wind,
                         const Real max_radius,
                         const Real inflow = 0);

         /// Returns wind speed at radius r
         Real
         get_V (const Real r) const;

         /// Returns wind speed at radius max_radius
         Real
         get_max_wind () const;

         /// Returns radius with max_wind
         const Real&
         get_max_radius () const;

   };

   class Fire
   {

      public:

         static Real
         get_gfdi (const Real t,
                   const Real rh,
                   const Real kph,
                   const Real curing = 100);

         static Real
         get_ffdi (const Real t,
                   const Real rh,
                   const Real kph,
                   const Real df = 10);

   };

   class Level
   {

      public:

         enum Type
         {
            PRESSURE,
            HEIGHT,
            THETA,
            SIGMA,
            SCREEN,
            FIFTY_METRE,
            TEN_METRE,
            MEAN_SEA,
            NIL,
            SURFACE,
            MODEL,
            NAL
         };

         Level::Type
         type;

         Real
         value;

         Real
         value_;

         Level ();

         Level (const Level& level);

         Level (const Dstring& str);

         Level (const Level::Type type,
                const Real value = GSL_NAN);

         Level (const Level::Type type,
                const Real value_0,
                const Real value_1);

         bool
         is_layer () const;

         void
         order ();

         static Level
         theta_level (const Real theta);

         static Level
         sigma_level (const Real sigma);

         static Level
         pressure_level (const Real pressure);

         static Level
         z_level (const Real z);

         static Level
         model_level (const Real m);

         static Level
         screen_level ();

         static Level
         fifty_metre_level ();

         static Level
         ten_metre_level ();

         static Level
         mean_sea_level ();

         static Level
         nil_level ();

         static Level
         surface_level ();

         void
         set_height (const Real z);

         void
         set_pressure (const Real pressure);

         void
         set_theta (const Real theta);

         void
         set_sigma (const Real sigma);

         void
         set_model (const Real m);

         void
         set_screen ();

         void
         set_fifty_metre ();

         void
         set_ten_metre ();

         void
         set_mean_sea ();

         void
         set_nil ();

         void
         set_surface ();

         void
         set_value (const Real value);

         Real
         get_value () const;

         Dstring
         get_string () const;

         Level
         get_level_0 () const;

         Level
         get_level_1 () const;

         void
         set_nal ();

         bool
         is_nal () const;

   };

   class Layer : public Level
   {

      public:

         Layer (const Level::Type type,
                const Real value_0,
                const Real value_1);

         void
         set (const Real value_0,
              const Real value_1);

   };

   class Z_Layer : public Level
   {

      public:

         Z_Layer (const Z_Layer& z_layer);

         Z_Layer (const Real start_z,
                  const Real end_z);

         void
         set (const Real start_z,
              const Real end_z);

         void
         set_start_z (const Real start_z);

         void
         set_end_z (const Real end_z);

         Real&
         get_start_z ();

         const Real&
         get_start_z () const;

         Real&
         get_end_z ();

         const Real&
         get_end_z () const;

         Real
         get_span_z () const;

   };

   class P_Layer : public Level
   {

      public:

         P_Layer (const P_Layer& p_layer);

         P_Layer (const Real start_p,
                  const Real end_p);

         void
         set (const Real start_p,
              const Real end_p);

         void
         set_start_p (const Real start_p);

         void
         set_end_p (const Real end_p);

         Real&
         get_start_p ();

         const Real&
         get_start_p () const;

         Real&
         get_end_p ();

         const Real&
         get_end_p () const;

         Real
         get_span_p () const;

         Real
         get_middle_p () const;

         Tuple
         get_tuple_p_specify_p (const Real specified_p,
                                const Real approx_delta_p = 10e2) const;

         Tuple
         get_tuple_p (const Real approx_delta_p = 10e2) const;

   };

}

#endif /* DENISE_MET_H */

