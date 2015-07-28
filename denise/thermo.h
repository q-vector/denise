//
// thermo.h
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

#ifndef DENISE_THERMO_H
#define DENISE_THERMO_H

#include <cmath>
#include <cairomm/context.h>
#include <denise/basics.h>
#include <denise/analysis.h>
#include <denise/geometry.h>
#include <denise/graphics.h>
#include <denise/dstring.h>
#include <denise/dtime.h>
#include <denise/met.h>

using namespace std;
using namespace denise;
using namespace Cairo;

/// denise
namespace denise
{

   /// Thermo Medium
   enum Thermo_Medium
   {
      WATER,
      ICE
   };

   //enum Thermo_Element
  // {
  //    MIXING_RATIO,
  //    TEMPERATURE,
  //    THETA,
  //    THETA_E
  // };

   /// Absolute zero
   const Real
   K = 273.15;

   /// Universal Gas Constant
   const Real
   R = 8314.32;

   /// Atomic mass of dry air
   const Real
   M_d = 28.9644;

   /// Atomic mass of water vapor
   const Real
   M_w = 18.0153;

   /// Constant pressure heat capacitity
   const Real
   c_p = 1005;

   /// Constant volume heat capacitity
   const Real
   c_v = 717;

   /// Latent heat
   const Real
   L = 2.501e6;

   /// Gravitational constant
   const Real
   g = 9.80665;

   /// Universal Gas Constant of dry air
   const Real
   R_d = 287.053072; // R / M_d

   /// Universal Gas Constant of water vapor
   const Real
   R_v = 461.5143794; // R / M_w

   /// \f$epsilon\f$
   const Real
   epsilon = 0.621980776; // R_d / R_v

   /// \f$kappa\f$
   const Real
   kappa = 0.285624947; // R_d / c_p

   class Sounding;
   class Thermo_Point;
   class Thermo_Polygon;
   class Thermo_Conserve;

   class Moisture
   {

      private:

         static Real
         get_t (const Real e_s,
                const Real tolerance_t,
                const Real start_t,
                const Real end_t,
                const Thermo_Medium thermo_medium);

      public:

         static Real
         get_e_s (const Real t,
                  const Thermo_Medium thermo_medium = WATER);

         static Real
         get_e_s (const Real t,
                  const Real p,
                  const Thermo_Medium thermo_medium = WATER);

         static Real
         get_r_s (const Real t,
                  const Real p,
                  const Thermo_Medium thermo_medium = WATER);

         static Real
         get_q_s (const Real t,
                  const Real p,
                  const Thermo_Medium thermo_medium = WATER);

         static Real
         get_t (const Real e_s,
                const Thermo_Medium thermo_medium = WATER);

         static Real
         get_rh (const Real t,
                 const Real t_d,
                 const Thermo_Medium thermo_medium = WATER);

         static Real
         get_t_d (const Real t,
                  const Real rh,
                  const Thermo_Medium thermo_medium = WATER);

         static Real
         get_t_v (const Real t,
                  const Real r);

   };

   /// Represents a point on the thermodynamic space
   ///
   /// The thermodynamic space is a 2-D space which can be
   /// represent by 3 variables: temperature, potential temperature,
   /// and pressure.  Knowledge of any 2 of the 3 variables determines
   /// the remaining variable through the Poission equation. p_0
   /// is the reference pressure in the Poisson equation and is
   /// set to 1000e2 by default.
   ///
   /// Thermo_Point extends Point_3D; the variables are kept
   /// in the inherited z, x and y fields of Point_3D.  Temperature
   /// is kept in the z field, theta in x field and pressure in y field.
   class Thermo_Point
   {

      protected:

         Real
         p_0;

         Real
         p;

         Real
         theta;

         Real
         t;

         Real&
         pressure;

         Real&
         temperature;

         void
         p_r_s_iterate (const Real p,
                        const Real r_s,
                        const Real tolerance_t = 0.02,
                        const Real start_t = -120,
                        const Real end_t = 80,
                        const Thermo_Medium thermo_medium = WATER,
                        const Real p_0 = 1000e2);

      public:

         Thermo_Point ();

         Thermo_Point (const Real t,
                       const Real theta,
                       const Real p,
                       const Real p_0 = 1000e2);

         Thermo_Point (const Thermo_Point& thermo_point);

         void
         set_t_p (const Real t,
                  const Real p,
                  const Real p_0 = 1000e2);
                                                                                
         static Thermo_Point
         t_p (const Real t,
              const Real p,
              const Real p_0 = 1000e2);
                                                                                
         void
         set_t_theta (const Real t,
                      const Real p,
                      const Real p_0 = 1000e2);
                                                                                
         static Thermo_Point
         t_theta (const Real t,
                  const Real theta,
                  const Real p_0 = 1000e2);
                                                                                
         void
         set_theta_p (const Real theta,
                      const Real p,
                      const Real p_0 = 1000e2);

         static Thermo_Point
         theta_p (const Real theta,
                  const Real p,
                  const Real p_0 = 1000e2);

         void
         set_t_r_s (const Real t,
                    const Real r_s,
                    const Thermo_Medium thermo_medium = WATER,
                    const Real p_0 = 1000e2);

         static Thermo_Point
         t_r_s (const Real t,
                const Real r_s,
                const Thermo_Medium thermo_medium = WATER,
                const Real p_0 = 1000e2);

         void
         set_theta_r_s (const Real theta,
                        const Real r_s,
                        const Real tolerance_t = 0.02,
                        const Real start_t = -120,
                        const Real end_t = 80,
                        const Thermo_Medium thermo_medium = WATER,
                        const Real p_0 = 1000e2);

         static Thermo_Point
         theta_r_s (const Real theta,
                    const Real r_s,
                    const Real tolerance_t = 0.02,
                    const Real start_t = -120,
                    const Real end_t = 80,
                    const Thermo_Medium thermo_medium = WATER,
                    const Real p_0 = 1000e2);

         void
         set_p_r_s (const Real p,
                    const Real r_s,
                    const Real tolerance_t = 0.02,
                    const Thermo_Medium thermo_medium = WATER,
                    const Real p_0 = 1000e2);

         static Thermo_Point
         p_r_s (const Real p,
                const Real r_s,
                const Real tolerance_t = 0.02,
                const Thermo_Medium thermo_medium = WATER,
                const Real p_0 = 1000e2);

         void
         set_t_theta_e (const Real t,
                        const Real theta_e,
                        const Real tolerance_p = 1,
                        const Real start_p = 1e2,
                        const Real end_p = 1200e2,
                        const Thermo_Medium thermo_medium = WATER,
                        const Real p_0 = 1000e2);

         static Thermo_Point
         t_theta_e (const Real t,
                    const Real theta_e,
                    const Real tolerance_p = 1,
                    const Real start_p = 1e2,
                    const Real end_p = 1200e2,
                    const Thermo_Medium thermo_medium = WATER,
                    const Real p_0 = 1000e2);

         void
         set_theta_theta_e (const Real theta,
                            const Real theta_e,
                            const Real tolerance_t = 0.02,
                            const Real start_t = -120,
                            const Real end_t = 80,
                            const Thermo_Medium thermo_medium = WATER,
                            const Real p_0 = 1000e2);

         static Thermo_Point
         theta_theta_e (const Real theta,
                        const Real theta_e,
                        const Real tolerance_t = 0.02,
                        const Real start_t = -120,
                        const Real end_t = 80,
                        const Thermo_Medium thermo_medium = WATER,
                        const Real p_0 = 1000e2);

         void
         set_p_theta_e (const Real p,
                        const Real theta_e,
                        const Real tolerance_t = 0.02,
                        const Real start_t = -120,
                        const Real end_t = 80,
                        const Thermo_Medium thermo_medium = WATER,
                        const Real p_0 = 1000e2);

         static Thermo_Point
         p_theta_e (const Real p,
                    const Real theta_e,
                    const Real tolerance_t = 0.02,
                    const Real start_t = -120,
                    const Real end_t = 80,
                    const Thermo_Medium thermo_medium = WATER,
                    const Real p_0 = 1000e2);

         void
         set_normand (const Real t,
                      const Real t_d,
                      const Real p,
                      const Real tolerance_t = 0.02,
                      const Real start_t = -120,
                      const Real end_t = 80,
                      const Thermo_Medium thermo_medium = WATER,
                      const Real p_0 = 1000e2);

         static Thermo_Point
         normand (const Real t,
                  const Real t_d,
                  const Real p,
                  const Real tolerance_t = 0.02,
                  const Real start_t = -120,
                  const Real end_t = 80,
                  const Thermo_Medium thermo_medium = WATER,
                  const Real p_0 = 1000e2);

         void
         operator = (const Thermo_Point& thermo_point);

         bool
         operator == (const Thermo_Point& thermo_point) const;

         bool
         operator > (const Thermo_Point& thermo_point) const;

         bool
         operator < (const Thermo_Point& thermo_point) const;

         const Real&
         get_p_0 () const;

         const Real&
         get_temperature () const;

         const Real&
         get_theta () const;

         const Real&
         get_pressure () const;     

         const Real&
         get_t () const;

         const Real&
         get_p () const;     

         /// Returns saturated vapour pressure (Pa)
         ///
         /// \param thermo_medium    thermo_medium
         Real
         get_e_s (const Thermo_Medium thermo_medium = WATER) const;

         /// Returns saturated mixing ratio (unit)
         ///
         /// \param thermo_medium    thermo_medium
         Real
         get_r_s (const Thermo_Medium thermo_medium = WATER) const;

         /// Returns saturated specific humidity (unit)
         ///
         /// \param thermo_medium    thermo_medium
         Real
         get_q_s (const Thermo_Medium thermo_medium = WATER) const;

         /// Returns Condensation Rate in an ascent
         ///
         /// \param thermo_medium    thermo_medium
         Real
         get_condense_rate (const Real omega,
                            const Thermo_Medium thermo_medium = WATER) const;

         /// Returns Latent Heating Rate from saturated ascent (W)
         ///
         /// \param thermo_medium    thermo_medium
         Real
         get_saturated_Q_dot (const Real omega,
                              const Thermo_Medium thermo_medium = WATER) const;

         /// Returns Temperature Dot from saturated ascent (W)
         ///
         /// \param thermo_medium    thermo_medium
         Real
         get_saturated_t_dot (const Real omega,
                              const Thermo_Medium thermo_medium = WATER) const;

         /// Returns Theta Dot from saturated ascent (W)
         ///
         /// \param thermo_medium    thermo_medium
         Real
         get_saturated_theta_dot (const Real omega,
                                  const Thermo_Medium thermo_medium = WATER) const;

         /// Returns Equivalent Potential Temperature
         ///
         /// \param thermo_medium    thermo_medium
         Real
         get_theta_e (const Thermo_Medium thermo_medium = WATER) const;

         /// Returns Wet Bulb Potential Temperature
         ///
         /// \param thermo_medium    thermo_medium
         Real
         get_theta_w (const Real tolerance_t = 0.02,
                      const Real start_t = -120,
                      const Real end_t = 80,
                      const Thermo_Medium thermo_medium = WATER) const;

   };

   class Instability
   {

      public:

         static Real
         get_k_index (const Real t_850,
                      const Real t_700,
                      const Real t_500,
                      const Real t_d_850,
                      const Real t_d_700);

         static Real
         get_total_totals (const Real t_850,
                           const Real t_500,
                           const Real t_d_850);

         static Real
         get_showalter_index (const Real t_850,
                              const Real t_500,
                              const Real t_d_850);

         static Real
         get_lifted_index (const Real start_p,
                           const Real start_t,
                           const Real start_t_d,
                           const Real end_p,
                           const Real end_t);

   };

   class Thermo_Diagram : public Affine_Transform_2D
   {

      protected:

         Size_2D
         size_2d;

         const Thermo_Point
         ref_thermo_point;

         const Real
         label_size;

         Point_2D
         anchor;

         const Real
         p_0;

         Real
         find_temperature (const Real x,
                           const Real p,
                           const Real start_t = -150,
                           const Real end_t = 70,
                           const Real epsilon = 0.01) const;

         virtual void
         transform (Point_2D& point_2d,
                    const Thermo_Point& thermo_point) const = 0;

         virtual void
         reverse_tp (Thermo_Point& thermo_point,
                     const Point_2D& point_2d) const = 0;

         void
         render_isobar (const RefPtr<Context> cr,
                        const Real p,
                        const Real start_t,
                        const Real delta_t,
                        const Integer n) const;

         void
         render_isobars (const RefPtr<Context> cr,
                         const Color& color_0,
                         const Color& color_1,
                         const Color& color_2,
                         const bool fine) const;

         void
         render_isotherms (const RefPtr<Context> cr,
                           const Color& color_0,
                           const Color& color_1,
                           const Color& color_2,
                           const bool fine) const;

         void
         render_dry_adiabats (const RefPtr<Context>& cr,
                              const Color& color_0,
                              const Color& color_1,
                              const Color& color_2) const;

         void
         render_saturated_adiabats (const RefPtr<Context>& cr,
                                    const Color& color_0,
                                    const Color& color_1,
                                    const Color& color_2,
                                    const bool fine) const;

         void
         render_isohumes (const RefPtr<Context>& cr,
                          const Color& color_0,
                          const Color& color_1,
                          const Color& color_2) const;

         void
         render_labels (const RefPtr<Context>& cr,
                        const Real label_x,
                        const Color& color_0,
                        const Color& color_1,
                        const Color& color_2) const;

      public:

         Thermo_Diagram (const Size_2D& size_2d,
                         const Real p_0 = 1000e2,
                         const Thermo_Point& ref_thermo_point = Thermo_Point::t_p (-40, 1000e3));

         ~Thermo_Diagram ();

         virtual void
         reset (const Size_2D& size_2d) = 0;

         const Size_2D&
         get_size_2d () const;

         virtual void
         set_anchor (const Point_2D& anchor);

         void
         zoom (const Point_2D& point,
               const Real scale);

         const Real
         get_p_0 () const;

         void
         render (const RefPtr<Context>& cr,
                 const Real label_x,
                 const Color& color_0 = Color (0.5, 0.5, 0.5),
                 const Color& color_1 = Color (0.7, 0.7, 0.7),
                 const Color& color_2 = Color (0.9, 0.9, 0.9),
                 const bool fine = true) const;

         Point_2D
         transform (const Thermo_Point& thermo_point) const;

         Thermo_Point
         get_thermo_point (const Point_2D& point) const;

         Thermo_Point
         get_thermo_point (const Real x,
                           const Real p) const;

         bool
         intersection (Thermo_Point& intersection,
                       const Thermo_Point& tp_aa,
                       const Thermo_Point& tp_ab,
                       const Thermo_Point& tp_ba,
                       const Thermo_Point& tp_bb) const;

         void
         render_isohume (const RefPtr<Context>& cr,
                         const Real r_s,
                         const Real start_p = 1050e2,
                         const Real end_p = GSL_NAN) const;

   };

   class International_Standard_Atmosphere
   {

      private:

         Tuple
         tuple_p;

         const Tuple
         tuple_t;

         const Tuple
         tuple_z;

         static Real
         get_p (const Real p_0,
                const Real gamma,
                const Real dz,
                const Real t_0);


         static Real
         get_dz (const Real p,
                 const Real p_0,
                 const Real t_0,
                 const Real gamma);

         Integer
         get_node (const Real p) const;

      public:

         International_Standard_Atmosphere ();

         const Tuple&
         get_tuple_p () const;

         Real
         get_z (const Real p) const;

         Real
         get_t (const Real p) const;

   };

   typedef International_Standard_Atmosphere Isa;

   class Thermo_Conserve : public P_Layer
   {

      public:

         Thermo_Conserve (const P_Layer& p_layer);

         virtual const Thermo_Point
         get_thermo_point (const Real p) const = 0;

         const Real
         get_t (const Real p) const;

   };

   class Isotherm : public Thermo_Conserve
   {

      protected:

         Real
         temperature;

      public:

         Isotherm (const Real temperature,
                   const P_Layer& p_layer);

         const Real&
         get_temperature () const;

         void
         set_temperature (const Real temperature);

         const Thermo_Point
         get_thermo_point (const Real p) const;

   };

   class Isohume : public Thermo_Conserve
   {

      protected:

         Real
         mixing_ratio;

         const Thermo_Medium
         thermo_medium;

         const Real
         tolerance_t;

      public:

         Isohume (const Real mixing_ratio,
                  const P_Layer& p_layer,
                  const Real tolerance_t = 0.02,
                  const Thermo_Medium thermo_medium = WATER);

         const Real&
         get_mixing_ratio () const;

         const Thermo_Medium
         get_thermo_medium () const;

         void
         set_mixing_ratio (const Real mixing_ratio);

         const Thermo_Point
         get_thermo_point (const Real p) const;

   };

   class Dry_Adiabat : public Thermo_Conserve
   {

      private:

         Real
         theta;

      public:

         Dry_Adiabat (const Real theta,
                      const P_Layer& p_layer);

         const Real&
         get_theta () const;

         void
         set_theta (const Real theta);

         const Thermo_Point
         get_thermo_point (const Real p) const;

         const Real
         get_t_v (const Real p,
                  const Real r) const;

   };

   class Moist_Adiabat : public Thermo_Conserve
   {

      private:

         Real
         theta_e;

         const Thermo_Medium
         thermo_medium;

         const Real
         tolerance_t;

      public:

         Moist_Adiabat (const Real theta_e,
                        const P_Layer& p_layer,
                        const Real tolerance_t = 0.02,
                        const Thermo_Medium thermo_medium = WATER,
                        const Real min_p = 80e2);

         const Real&
         get_theta_e () const;

         const Thermo_Medium
         get_thermo_medium () const;

         void
         set_theta_e (const Real theta);

         const Thermo_Point
         get_thermo_point (const Real p) const;

         const Real
         get_t_v (const Real p) const;

   };

   class Mixed_Layer : public Thermo_Conserve
   {

      private:

         Dry_Adiabat
         dry_adiabat;

         Moist_Adiabat
         moist_adiabat;

         Thermo_Point
         normand;

         Thermo_Medium
         thermo_medium;

      public:

         Mixed_Layer (const Thermo_Point& normand,
                      const P_Layer& p_layer,
                      const Real tolerance_t = 0.02,
                      const Thermo_Medium thermo_medium = WATER,
                      const Real min_p = 80e2);

         const bool
         is_moist (const Real p) const;

         const Thermo_Point&
         get_normand () const;

         const Thermo_Medium
         get_thermo_medium () const;

         void
         set_normand (const Thermo_Point& normand);

         const Real
         get_theta_e () const;

         const Thermo_Point
         get_thermo_point (const Real p) const;

         const Real
         get_t_v (const Real p) const;

         void
         render (const RefPtr<Context>& cr,
                 const Thermo_Diagram& thermo_diagram,
                 const bool render_isohume = true) const;

   };

   class Draft : public Mixed_Layer
   {

      protected:

         Thermo_Polygon*
         thermo_polygon_ptr;

      public:

         Draft (const Thermo_Point& normand,
                const P_Layer& p_layer);

         ~Draft ();

         void
         update_thermo_polygon (const Thermo_Diagram& thermo_diagram,
                                const Sounding& sounding,
                                const bool use_virtual);

   };

   class Updraft : public Draft
   {

      public:

         Updraft (const Real start_pressure,
                  const P_Layer& p_layer);

         Updraft (const Thermo_Point& normand,
                  const P_Layer& p_layer);

         Real
         get_total_cape () const;

         Real
         get_total_cin () const;

   };

   class Downdraft : public Draft
   {

      private:

         void
         render_dmape_area (const RefPtr<Context>& cr,
                            const Thermo_Diagram& thermo_diagram,
                            const Real line_width) const;

      public:

         Downdraft (const Thermo_Point& normand,
                    const P_Layer& p_layer);

         Real
         get_dmape () const;

         void
         render (const RefPtr<Context>& cr,
                 const Thermo_Diagram& thermo_diagram,
                 const Real line_width,
                 const bool use_virtual) const;

   };

   class Thermo_Polygon : public Polygon
   {

      private:

         const Thermo_Diagram& 
         thermo_diagram;

         const bool
         positive;

         P_Layer
         p_layer;

      public:

         Thermo_Polygon (const Thermo_Diagram& thermo_diagram,
                         const bool positive);

         bool
         is_positive () const;

         Real
         get_energy () const;

         static Real
         get_energy (const Polygon& polygon,
                     const Thermo_Diagram& thermo_diagram);

         Real
         get_positive_energy () const;

         static Real
         get_positive_energy (const Polygon& polygon,
                              const Thermo_Diagram& thermo_diagram);

         Real
         get_negative_energy () const;

         static Real
         get_negative_energy (const Polygon& polygon,
                              const Thermo_Diagram& thermo_diagram);

         void
         add (const Thermo_Point& thermo_point);

         void
         cairo (const Cairo::RefPtr<Cairo::Context>& cr) const;

         void
         render (const Cairo::RefPtr<Cairo::Context>& cr,
                 const Color& color_positive,
                 const Color& color_negative) const;

         bool
         operator == (const Thermo_Polygon& thermo_polygon) const;

         bool
         operator > (const Thermo_Polygon& thermo_polygon) const;

         bool
         operator < (const Thermo_Polygon& thermo_polygon) const;

   };

   class Thermo_Line : public Real_Profile
   {

      private:

         void
         process_bracketing_iterators (Thermo_Line::const_iterator& lb,
                                       Thermo_Line::const_iterator& ub) const;

      public:

         Thermo_Line ();

         Thermo_Line (const Thermo_Conserve& thermo_conserve);

         Thermo_Line (const Mixed_Layer& mixed_layer,
                      const bool use_virtual,
                      const Real delta_p = 10e2);

         Thermo_Line (const Thermo_Diagram& thermo_diagram,
                      const Sounding& sounding,
                      const bool use_virtual,
                      const Thermo_Medium thermo_medium = WATER);

         Thermo_Line (const International_Standard_Atmosphere& isa,
                      const Real delta_p = 10e2);

         void
         write (ofstream& file,
                const string& identifier) const;

         set<Real>
         get_p_set () const;

         void
         add (const Real p,
              const Real value);

         void
         add_node (const Thermo_Diagram& thermo_diagram,
                   const Real p);

         void
         delete_node (const Thermo_Line::iterator iterator);

         Thermo_Line::iterator
         get_iterator (const Thermo_Diagram& thermo_diagram,
                       const Point_2D& point,
                       const Real tolerance);

         Thermo_Line::const_iterator
         get_iterator (const Thermo_Diagram& thermo_diagram,
                       const Point_2D& point,
                       const Real tolerance) const;

         Thermo_Line::iterator
         get_iterator (const Real p);

         Thermo_Line::const_iterator
         get_iterator (const Real p) const;

         Thermo_Point
         get_thermo_point (const Thermo_Diagram& thermo_diagram,
                           const Real p) const;

         Real
         get_nearest_p (const Real p) const;

         Real
         get_start_p () const;

         Real
         get_end_p () const;

         void
         intersection_t (const Thermo_Diagram& thermo_diagram,
                         set<Thermo_Point>& tp_set,
                         const Real t,
                         const Real lower_bound_p,
                         const Real upper_bound_p,
                         const Real tolerance_p = 10e2) const;

         void
         intersection_theta (const Thermo_Diagram& thermo_diagram,
                             set<Thermo_Point>& tp_set,
                             const Real theta,
                             const Real lower_bound_p,
                             const Real upper_bound_p,
                             const Real tolerance_p = 10e2) const;

         void
         intersection_theta_e (const Thermo_Diagram& thermo_diagram,
                               set<Thermo_Point>& tp_set,
                               const Real theta_e,
                               const Real lower_bound_p,
                               const Real upper_bound_p,
                               const Real tolerance_p = 10e2,
                               const Thermo_Medium thermo_medium = WATER) const;

         Real
         get_mean_theta (const Thermo_Diagram& thermo_diagram,
                         const P_Layer& p_layer) const;

         Real
         get_mean_mixing_ratio (const Thermo_Diagram& thermo_diagram,
                                const P_Layer& p_layer) const;

         void
         render (const RefPtr<Context>& cr,
                 const Thermo_Diagram& thermo_diagram) const;

         void
         render (const RefPtr<Context>& cr,
                 const Thermo_Diagram& thermo_diagram,
                 const P_Layer& p_layer,
                 const Real start_x,
                 const Real end_x) const;

         void
         render_node (const RefPtr<Context>& cr,
                      const Thermo_Diagram& thermo_diagram,
                      Thermo_Line::const_iterator iterator,
                      const bool fill,
                      const Real node_size) const;

         static Real
         get_top_intersection_p (const Thermo_Diagram& thermo_diagram,
                                 const Thermo_Line& thermo_line_a,
                                 const Thermo_Line& thermo_line_b);

         static Thermo_Polygon*
         get_thermo_polygon_ptr (const Thermo_Diagram& thermo_diagram,
                                 const Thermo_Line& thermo_line_up,
                                 const Thermo_Line& thermo_line_down,
                                 const P_Layer& p_layer);

         set<Thermo_Point>
         get_tp_set_t (const Thermo_Diagram& thermo_diagram,
                       const Real t,
                       const Real tolerance_p = 1) const;


         set<Thermo_Point>
         get_intersection_tp_set (const Thermo_Diagram& thermo_diagram,
                                  const Mixed_Layer& mixed_layer,
                                  const Real tolerance_p = 1) const;

         Thermo_Polygon*
         get_thermo_polygon_ptr (const Thermo_Diagram& thermo_diagram,
                                 const Mixed_Layer& mixed_layer,
                                 const P_Layer& p_layer,
                                 const bool use_virtual) const;

         set<Thermo_Polygon*>
         get_thermo_polygon_ptr_set (const Thermo_Diagram& thermo_diagram,
                                     const Mixed_Layer& mixed_layer,
                                     const set<Thermo_Point>& tp_set) const;

   };

   class T_Line : public Thermo_Line
   {

      public:

         T_Line ();

   };

   class T_D_Line : public Thermo_Line
   {

      public:

         T_D_Line ();

   };

   class T_V_Line : public Thermo_Line
   {

      public:

         T_V_Line (const Thermo_Diagram& thermo_diagram,
                   const Sounding& sounding,
                   const Thermo_Medium thermo_medium = WATER);
   };

   class T_W_Line : public Thermo_Line
   {

      public:

         T_W_Line (const Thermo_Diagram& thermo_diagram,
                   const Sounding& sounding,
                   const Thermo_Medium& thermo_medium = WATER);

   };

   class Wind_Profile : public map<Real, Wind>
   {

      private:

         void
         process_bracketing_iterators (Wind_Profile::const_iterator& lb,
                                       Wind_Profile::const_iterator& ub) const;

      public:

         Wind_Profile ();

         void
         write (ofstream& file) const;

         set<Real>
         get_p_set () const;

         void
         add (const Real p,
              const Wind& wind);

         Wind
         get_wind (const Real p) const;

         Wind
         get_mean_wind (const P_Layer& p_layer) const; 

         Wind
         get_wind_shear (const P_Layer& p_layer) const;

   };

   class Height_Profile : public Real_Profile
   {

      private:

         void
         process_bracketing_iterators (Height_Profile::const_iterator& lb,
                                       Height_Profile::const_iterator& ub) const;

         Real
         get_pressure (const Real this_p,
                       const Real this_z,
                       const Real next_p,
                       const Real next_z,
                       const Real height) const;

      public:

         Height_Profile ();

         void
         write (ofstream& file) const;

         set<Real>
         get_p_set () const;

         void
         add (const Real pressure,
              const Real height);

         Real
         get_height (const Real pressure) const;

         Real
         get_pressure (const Real height) const;

   };

   class Sounding
   {

      protected:

         Dtime
         time;

         Dtime
         basetime;

         string
         location_str;

         T_Line
         t_line;

         T_D_Line
         t_d_line;

         Wind_Profile
         wind_profile;

         Height_Profile
         height_profile;

         P_Layer
         dry_layer;

         Z_Layer
         shear_layer;

         Z_Layer
         steering_layer;

         Z_Layer
         helicity_layer;

         T_W_Line*
         t_w_line_ptr;

         Updraft*
         updraft_ptr;

         Downdraft*
         downdraft_ptr;

         void
         render_thermo_line_nodes (const RefPtr<Context>& cr,
                                   const Thermo_Diagram& thermo_diagram,
                                   const Thermo_Line& thermo_line,
                                   const Real node_size) const;

      public:

         Sounding ();

         Sounding (const string& file_path);

         Sounding (const Sounding& sounding);

         void
         load (const string& file_path);

         void
         save (const string& file_path) const;

         void
         write (ofstream& file) const;

         void
         set_time (const Dtime& dtime);

         void
         set_basetime (const Dtime& basetime);

         void
         set_location_str (const string& location_str);

         static Sounding*
         get_mean_sounding_ptr (const list<const Sounding*>& sounding_ptr_list,
                                const Thermo_Diagram& thermo_diagram);

         Real_Profile*
         get_scorer_profile_ptr (const Real azimuth,
                                 const Thermo_Diagram& thermo_diagram) const;

         Real_Profile*
         get_brunt_vaisala_profile_ptr () const;

         set<Real>
         get_p_set () const;

         const Dtime&
         get_time () const;

         T_Line&
         get_t_line ();

         const T_Line&
         get_t_line () const;

         T_D_Line&
         get_t_d_line ();

         const T_D_Line&
         get_t_d_line () const;

         Wind_Profile&
         get_wind_profile ();

         const Wind_Profile&
         get_wind_profile () const;

         Height_Profile&
         get_height_profile ();

         const Height_Profile&
         get_height_profile () const;

         const T_W_Line&
         get_t_w_line () const;

         const Updraft&
         get_updraft () const;

         const Downdraft&
         get_downdraft () const;

         void
         update_t_w_line (const Thermo_Diagram& thermo_diagram,
                          const Thermo_Medium& thermo_medium);

         void
         update_updraft (const Thermo_Diagram& thermo_diagram,
                         const Real start_pressure);

         void
         update_updraft (const Thermo_Diagram& thermo_diagram,
                         const Thermo_Point& normand);

         void
         update_downdraft (const Thermo_Diagram& thermo_diagram,
                           const Thermo_Point& start_tp);

         Real
         get_dry_layer_theta_e (const Thermo_Diagram& thermo_diagram,
                                const Real tolerance_t = 0.2,
                                const Thermo_Medium thermo_medium = WATER) const;

         Real
         get_mean_theta_e (const Thermo_Diagram& thermo_diagram,
                           const P_Layer& p_layer,
                           const Real tolerance_t = 0.2,
                           const Thermo_Medium thermo_medium = WATER) const;

         void
         add_node (const Thermo_Diagram& thermo_diagram,
                   const Real pressure);

         void
         add_node_t (const Thermo_Diagram& thermo_diagram,
                     const Real pressure);

         void
         add_node_t_d (const Thermo_Diagram& thermo_diagram,
                       const Real pressure);

         void
         delete_thermo_node (const Thermo_Line::iterator iterator);

         void
         modify_by (const Mixed_Layer& mixed_layer);

         void
         mix (const Thermo_Diagram& thermo_diagram,
              const P_Layer& p_layer);

         void
         surface_heat (const Thermo_Diagram& thermo_diagram,
                       const Real theta);

         Thermo_Line::iterator
         get_thermo_node (const Thermo_Diagram& thermo_diagram,
                          const Point_2D& point,
                          const Real tolerance);

         Thermo_Line::const_iterator
         get_thermo_node (const Thermo_Diagram& thermo_diagram,
                          const Point_2D& point,
                          const Real tolerance) const;

         Thermo_Line::iterator
         get_thermo_line_iterator (const Thermo_Diagram& thermo_diagram,
                                   bool& is_t_line,
                                   const Point_2D& point,
                                   const Real tolerance);

         Thermo_Line::const_iterator
         get_thermo_line_iterator (const Thermo_Diagram& thermo_diagram,
                                   bool& is_t_line,
                                   const Point_2D& point,
                                   const Real tolerance) const;

         Thermo_Point
         get_prev_thermo_point (const Thermo_Diagram& thermo_diagram,
                                const Thermo_Line::const_iterator i) const;

         Thermo_Point
         get_next_thermo_point (const Thermo_Diagram& thermo_diagram,
                                const Thermo_Line::const_iterator i) const;

         Real
         get_temperature (const Thermo_Diagram& thermo_diagram,
                          const Real pressure) const;

         Real
         get_dew_point (const Thermo_Diagram& thermo_diagram,
                        const Real pressure) const;

         Real
         get_mixing_ratio (const Thermo_Diagram& thermo_diagram,
                           const Real pressure,
                           const Thermo_Medium thermo_medium = WATER) const;

         Real
         get_virtual_temperature (const Thermo_Diagram& thermo_diagram,
                                  const Real pressure,
                                  const Thermo_Medium thermo_medium = WATER) const;

         Real
         get_wet_bulb (const Thermo_Diagram& thermo_diagram,
                       const Real pressure,
                       const Thermo_Medium thermo_medium = WATER) const;

         Real
         get_height (const Real pressure) const;

         Real
         get_pressure (const Real height) const;

         Wind
         get_wind (const Real pressure) const;

         Wind
         get_mean_wind (const P_Layer& p_layer) const;

         Wind
         get_wind_z (const Real height) const;

         Wind
         get_mean_wind (const Z_Layer& z_layer) const;

         Real
         get_nearest_p (const Real p) const;

         Real
         get_start_p () const;

         Real
         get_end_p () const;

         void
         render (const RefPtr<Context>& cr,
                 const Thermo_Diagram& thermo_diagram,
                 const Real node_size = GSL_NAN) const;

         void
         render_t (const RefPtr<Context>& cr,
                   const Thermo_Diagram& thermo_diagram) const;

         void
         render_t_d (const RefPtr<Context>& cr,
                     const Thermo_Diagram& thermo_diagram) const;

         void
         render_t_nodes (const RefPtr<Context>& cr,
                         const Thermo_Diagram& thermo_diagram,
                         const Real node_size) const;

         void
         render_t_d_nodes (const RefPtr<Context>& cr,
                           const Thermo_Diagram& thermo_diagram,
                           const Real node_size) const;

         void
         render_winds (const RefPtr<Context>& cr,
                       const Thermo_Diagram& thermo_diagram,
                       const Real x = GSL_NAN) const;

         void
         render_heights (const RefPtr<Context>& cr,
                         const Thermo_Diagram& thermo_diagram,
                         const Real x = GSL_NAN) const;

         Mixed_Layer
         get_mixed_layer_from (const Thermo_Diagram& thermo_diagram,
                               const Real p) const;

         Mixed_Layer
         get_mixed_layer_surface (const Thermo_Diagram& thermo_diagram,
                                  const Real mixing_dew_point_depth = 50e2) const;

         set<Thermo_Point>
         get_tp_set_t (const Thermo_Diagram& thermo_diagram,
                       const Real t,
                       const Real tolerance_p = 1e2) const;
 
         set<Thermo_Point>
         get_intersection_tp_set (const Thermo_Diagram& thermo_diagram,
                                  const Mixed_Layer& mixed_layer,
                                  const Thermo_Line& thermo_line,
                                  const Real tolerance_p = 1e2) const;

         Thermo_Polygon*
         get_thermo_polygon_ptr (const Thermo_Diagram& thermo_diagram,
                                 const Mixed_Layer& mixed_layer,
                                 const Thermo_Line& thermo_line,
                                 const P_Layer& p_layer) const;

         set<Thermo_Polygon*>
         get_thermo_polygon_ptr_set (const Thermo_Diagram& thermo_diagram,
                                     const Mixed_Layer& mixed_layer,
                                     const Thermo_Line& thermo_line,
                                     const set<Thermo_Point>& tp_set) const;

         set<Thermo_Polygon*>
         get_thermo_polygon_ptr_set (const Thermo_Diagram& thermo_diagram,
                                     const Mixed_Layer& mixed_layer) const;

         Thermo_Polygon*
         get_thermo_polygon_ptr (const Thermo_Diagram& thermo_diagram,
                                 const Mixed_Layer& mixed_layer,
                                 const bool use_virtual) const;

         set<Thermo_Point>
         get_wbfztp_set (const Thermo_Diagram& thermo_diagram) const;

         Real
         get_lfs (const Thermo_Diagram& thermo_diagram,
                  const Thermo_Line& downdraft_tl,
                  const Thermo_Line& environment_tl) const;

         Real
         get_cross_totals (const Thermo_Diagram& thermo_diagram) const;

         Real
         get_vertical_totals (const Thermo_Diagram& thermo_diagram) const;

         Real
         get_total_totals (const Thermo_Diagram& thermo_diagram) const;

         Real
         get_k_index (const Thermo_Diagram& thermo_diagram) const;

         Real
         get_lifted_index (const Thermo_Diagram& thermo_diagram,
                           const Mixed_Layer& mixed_layer) const;

         Real
         get_lifted_index (const Thermo_Diagram& thermo_diagram,
                           const Real p) const;

         Real
         get_showalter_index (const Thermo_Diagram& thermo_diagram) const;

         Real
         get_surface_lifted_index (const Thermo_Diagram& thermo_diagram) const;

         Real
         get_total_cin () const;

         Real
         get_total_cape () const;

         Real
         get_dmape () const;

         Real
         get_convective_gust (const Real dmape_factor) const;

         Real
         get_warm_cloud_depth () const;

         Real
         get_bulk_richardson_number () const;

         Real
         get_helicity () const;

   };

   class Ttxx_Sounding : public Sounding
   {

      private:

         class Ttxx : public Tokens
         {

            protected:

               string
               yygg;

               Integer
               wmo_id;

               Ttxx (const string& str);

               bool
               use_knots () const;

               string
               get_yygg () const;

               const string&
               get_wmo_id () const;

               string
               get_key () const;

               pair<Real, Real>
               parse_tttdd (const string& tttdd) const;

               Wind
               parse_ddfff (const string& ddfff) const;

            public:

               virtual void
               parse_to (Sounding& sounding) const = 0;

         };

         class Ttac : public Ttxx
         {

            protected:

               const char
               get_highest_wind_code () const;

               Integer
               parse_special_level (Sounding& sounding,
                                    Integer& index) const;

               Integer
               parse_standard_levels (Sounding& sounding,
                                      Integer& index) const;

               pair<Real, Real>
               parse_pphhh (const string& pphhh) const;

               virtual void
               interpret_p_z (const Integer pp,
                              Real& pressure,
                              Real& geopotential_height) const = 0;

               virtual Integer
               get_number_of_standard_levels () const = 0;

            public:

               Ttac (const string& str);

         };

         class Ttbd : public Ttxx
         {

            protected:

               virtual Real
               parse_nnppp (const string& nnppp) const = 0;

               bool
               parse_significant_temperature_levels (Sounding& sounding,
                                                     Integer& index) const;

               bool
               parse_significant_wind_levels (Sounding& sounding,
                                              Integer& index) const;

            public:

               Ttbd (const string& str);

         };

         class Ttaa : public Ttac
         {

            protected:

               void
               interpret_p_z (const Integer pp,
                              Real& pressure,
                              Real& geopotential_height) const;

               Integer
               get_number_of_standard_levels () const;

            public:

               Ttaa (const string& str);

               void
               parse_to (Sounding& sounding) const;

         };

         class Ttbb : public Ttbd
         {

            protected:

               Real
               parse_nnppp (const string& nnppp) const;

            public:

               Ttbb (const string& str);

               void
               parse_to (Sounding& sounding) const;

         };

         class Ttcc : public Ttac
         {

            protected:

               void
               interpret_p_z (const Integer pp,
                              Real& pressure,
                              Real& geopotential_height) const;

               Integer
               get_number_of_standard_levels () const;

            public:

               Ttcc (const string& str);

               void
               parse_to (Sounding& sounding) const;

         };

         class Ttdd : public Ttbd
         {

            protected:

               Real
               parse_nnppp (const string& nnppp) const;

            public:

               Ttdd (const string& str);

               void
               parse_to (Sounding& sounding) const;

         };

         Integer
         yygg;

      public:

         Ttxx_Sounding (const Integer wmo_id,
                        const Integer yygg);

         const string
         get_key () const;

   };

   class Tephigram : public Thermo_Diagram
   {

      private:

         void
         transform (Point_2D& point_2d,
                    const Thermo_Point& thermo_point) const;

         void
         reverse_tp (Thermo_Point& thermo_point,
                     const Point_2D& point_2d) const;

      public:

         Tephigram (const Size_2D& size_2d,
                    const Real p_0 = 1000e2,
                    const Thermo_Point& ref_thermo_point = Thermo_Point::t_p (-40, 1000e2));

         void
         reset (const Size_2D& size_2d);

         Real
         get_jacobian () const;

   };

   class Emagram : public Thermo_Diagram
   {

      protected:

         const Real
         magic_ratio;

         void
         transform (Point_2D& point_2d,
                    const Thermo_Point& thermo_point) const;

         void
         reverse_tp (Thermo_Point& thermo_point,
                     const Point_2D& point_2d) const;

      public:

         Emagram (const Size_2D& size_2d,
                  const Real magic_ratio = 45,
                  const Real p_0 = 1000e2,
                  const Thermo_Point& ref_thermo_point = Thermo_Point::t_p (-120, 1000e2));

         virtual void
         reset (const Size_2D& size_2d);

         Real
         get_jacobian () const;

   };

   class Skew_T : public Emagram
   {

      public:

         Skew_T (const Size_2D& size_2d,
                 const Real magic_ratio = 36.6,
                 const Real p_0 = 1000e2,
                 const Thermo_Point& ref_thermo_point = Thermo_Point::t_p (-40, 1000e2));

         void
         reset (const Size_2D& size_2d);

   };

   class Thermo_Exception : public Exception
   {

      public:

         Thermo_Exception (const string& description = string (""));

   };

   ostream&
   operator << (ostream &out_file,
                const Thermo_Point& thermo_point);

}

#endif /* DENISE_THERMO_H */
