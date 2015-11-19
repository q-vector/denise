//
// streamline.h
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

#ifndef DENISE_STREAMLINE_H
#define DENISE_STREAMLINE_H

#include <gsl/gsl_poly.h>
#include <cairomm/context.h>
#include <denise/analysis.h>
#include <denise/geometry.h>
#include <denise/graphics.h>
#include <denise/util.h>

using namespace std;
using namespace denise;

namespace denise
{

   enum Critical_Point_Type_2D
   {
      SADDLE,
      VORTEX,
      SOURCE_SINK
   };

/*
   class Duple_Field_2D
   {

      public:

         const Tuple_Field_2D&
         tuple_field_2d;

         const Integer
         tuple_u_index;

         const Integer
         tuple_v_index;

         const Real
         u_multiplier;

         const Real
         v_multiplier;

         Duple_Field_2D (const Tuple_Field_2D& tuple_field_2d,
                         const Integer tuple_u_index,
                         const Integer tuple_v_index,
                         const Real u_multiplier,
                         const Real v_multiplier);

         Real
         evaluate_u_at (const Real x,
                        const Real y,
                        const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate_u_at (const Point_2D& point,
                        const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate_v_at (const Real x,
                        const Real y,
                        const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate_v_at (const Point_2D& point,
                        const Evaluate_Op evaluate_op = VALUE) const;

   };

*/

   class Intensity_Data
   {

      private:

         Real*
         data;

         const Index_2D
         anchor;

         const Size_2D
         size_2d;

      public:

         Intensity_Data (const Box_2D& box_2d);

         ~Intensity_Data ();

         void
         set_noise ();

         void
         set (const Real intensity,
              const Integer i,
              const Integer j);

         void
         set_raw (const Real intensity,
                  const Integer i,
                  const Integer j);

         Real&
         get (const Integer i,
              const Integer j) const;

         Real&
         get_raw (const Integer i,
                  const Integer j) const;

         void
         enhance (const Real min_intensity,
                  const Real max_intensity);

   };

   class Hits
   {

      private:

         uint8_t*
         data;

         const Index_2D
         anchor;

         const Size_2D
         size_2d;

      public:

         Hits (const Box_2D& box_2d);

         ~Hits ();

         uint8_t&
         get (const Index_2D& index_2d) const;

   };

   class Streamliner : public Box_2D
   {

      protected:

         const Integration_Scheme
         integration_scheme;

         const Transform_2D&
         transform;

         const Vector_Data_2D&
         vector_data_2d;

         const Integer
         u_index;

         const Integer
         v_index;

         const Real
         aspect;

         void
         step (Point_2D& delta,
               const Real x,
               const Real y,
               const Real h) const;

         void
         integrate (Point_2D& next_point,
                    const Point_2D& point_2d,
                    const Real h) const;

         Streamliner (const Box_2D& box_2d,
                      const Transform_2D& transform,
                      const Vector_Data_2D& vector_data_2d,
                      const Integer u_index,
                      const Integer v_index,
                      const Integration_Scheme integration_scheme,
                      const Real aspect = 1);

   };

   class Lic_Raster : public Raster,
                      public Streamliner
   {

      private:

         const Integer
         L;

         const Integer
         M;

         const Integer
         L_plus_M;

         const Real
         h;

         void
         grow (const Point_2D& point_2d,
               Point_2D* point_array,
               Real* noise_array,
               const Intensity_Data& white_noise) const;

         bool
         out_of_bounds (const Point_2D& point_2d) const;

         void
         get_sequence (Index_2D* sequence_array,
                       const Box_2D& box_2d) const;

         Intensity_Data*
         get_slow_lic_intensity_data_ptr (Real& min_intensity,
                                          Real& max_intensity,
                                          const Box_2D& box_2d) const;

         Intensity_Data*
         get_lic_intensity_data_ptr (Real& min_intensity,
                                     Real& max_intensity,
                                     const Box_2D& box_2d) const;

      public:

         Lic_Raster (const Box_2D& box_2d,
                     const Transform_2D& transform,
                     const Vector_Data_2D& vector_data_2d,
                     const Integer u_index,
                     const Integer v_index,
                     const bool enhance = false,
                     const Real aspect = 1,
                     const Integration_Scheme integration_scheme = EULER,
                     const Integer L = 40,
                     const Integer M = 100,
                     const Real h = 0.8);

   };

   class Streamlines
   {

      public:

         class Uv_Data
         {               
                         
            public:      
                         
               const Vector_Data_2D&
               data;
                         
               const Integer
               u_index;  
                         
               const Integer
               v_index;
         
               const Real
               u_multiplier;
         
               const Real
               v_multiplier;
         
               Uv_Data (const Vector_Data_2D& data,
                        const Integer u_index,
                        const Integer v_index,
                        const Real u_multiplier = 1,
                        const Real v_multiplier = 1);
         
         };

         class Singular_Point : public Point_2D
         {

            private:

               Jacobian_2D
               jacobian;

            public:

               Singular_Point (const Point_2D& point,
                               const Jacobian_2D& jacobian);

               Singular_Point (const Singular_Point& singular_point);

               bool
               is_saddle () const;

               pair<Real, Real>
               get_eigenvalue_pair () const;

               pair<Real, Real>
               get_eigentheta_pair () const;

               Real
               get_smaller_eigentheta () const;

         };

         class Singular_Points : public vector<Singular_Point>
         {

            private:

               Real
               residual;

               Real
               epsilon;

               const Uv_Data&
               uv_data;

               static int
               gsl_get_value (const gsl_vector* x,
                              void* p,
                              gsl_vector* f);

               void
               potentially_add (const Uv_Data& uv_data,
                                const Real a,
                                const Real c,
                                const Real delta_x,
                                const Real delta_y,
                                const Real s,
                                const Real t);

            public:

               Singular_Points (const Uv_Data& uv_data,
                                const Domain_2D& domain_2d,
                                const Real residual = 1e-9,
                                const Real epsilon = 1e-9);

               list<Singular_Point>
               get_saddle_point_list () const;

         };

         class Sand_Box : public Grid_nD
         {

            private:

               list<Point_2D>*
               point_lists;

               const Real
               separation;

            public:

               Integer
               get_offset (const Integer i,
                           const Integer j) const;

            public:

               Sand_Box (const Domain_2D& domain_2d,
                         const Real separation);

               ~Sand_Box ();

               bool
               same_box (const Point_2D& point_0,
                         const Point_2D& point_1) const;

               void
               add (const Point_2D& point);

               bool
               too_close (const Point_2D& point,
                          const Real theta = GSL_NAN) const;

         };


         class Streamline : public Simple_Polyline
         {

            private:

               void
               set_nan (Real& dx,
                        Real& dy) const;

            protected:

               bool
               forward_stopped;

               bool
               backward_stopped;

               void
               euler (Real& dx,
                      Real& dy,
                      const Uv_Data& uv_data,
                      const Real x,
                      const Real y,
                      const Real h) const;

               void
               runga_kutta (Real& dx,
                            Real& dy,
                            const Uv_Data& uv_data,
                            const Real x,
                            const Real y,
                            const Real h) const;

               Point_2D
               integrate (const Uv_Data& uv_data,
                          const Real x,
                          const Real y,
                          const Real h,
                          const Integration_Scheme integration_scheme) const;

            public:

               const bool
               forward;

               const bool
               backward;

               Streamline (const bool forward,
                           const bool backward,
                           const bool forward_stopped,
                           const bool backward_stopped);

               Streamline (const Point_2D& point,
                           const bool forward,
                           const bool backward);

               Streamline (const Point_2D& point_0,
                           const Point_2D& point_1,
                           const bool forward,
                           const bool backward);

               void
               stop ();

               void
               forward_stop ();

               void
               backward_stop ();

               bool
               is_stopped () const;

               bool
               is_forward_stopped () const;

               bool
               is_backward_stopped () const;

               const Point_2D&
               tip () const;

               void
               add_to (Sand_Box& sand_box) const;

               void
               grow_step (const Uv_Data& uv_data,
                          const Domain_2D& domain_2d,
                          Sand_Box& sand_box,
                          const Real h,
                          const Integration_Scheme integration_scheme);

               void
               grow (const Uv_Data& uv_data,
                     const Domain_2D& domain_2d,
                     Sand_Box& sand_box,
                     const Real h,
                     const Integration_Scheme integration_scheme);

               void
               render (const RefPtr<Context>& cr,
                       const Transform_2D& transform) const;

         };

         class Separatrix : public Streamline
         {

            private:

               void
               grow_step (const Uv_Data& uv_data,
                          const Domain_2D& domain_2d,
                          Sand_Box& sand_box,
                          const Real h,
                          const Integration_Scheme integration_scheme);

            public:

               Separatrix (const bool forward,
                           const bool backward,
                           const bool forward_stopped,
                           const bool backward_stopped);

               Separatrix (const Point_2D& point,
                           const bool outgoing);

               const Point_2D&
               get_origin () const;

               const Point_2D&
               get_tip () const;

               void
               grow (const Uv_Data& uv_data,
                     const Domain_2D& domain_2d,
                     Sand_Box& sand_box,
                     const Real step_size,
                     const Integration_Scheme integration_scheme);

               void
               render (const RefPtr<Context>& cr,
                       const Transform_2D& transform) const;

         };

         class Web_Stick : public Separatrix
         {

            private:

               void
               grow_step (const Uv_Data& uv_data,
                          const Domain_2D& domain_2d,
                          const Real step_size,
                          const Integration_Scheme integration_scheme);

            public:

               Web_Stick (const Singular_Point& saddle_point,
                          const Point_2D& point,
                          const bool outgoing);

               void
               grow (const Uv_Data& uv_data,
                     const Domain_2D& domain_2d,
                     const Integer steps,
                     const Real h,
                     const Integration_Scheme integration_scheme);

         };

         class Saddle : public Singular_Point,
                        public list<Web_Stick>
         {

            public:

               Saddle (const Singular_Point& saddle_point,
                       const Uv_Data& uv_data,
                       const Domain_2D& domain_2d,
                       const Integer size,
                       const Real step_size,
                       const Integration_Scheme integration_scheme);

               void
               add_to (Sand_Box& sand_box) const;

               void
               render (const RefPtr<Context>& cr,
                       const Transform_2D& transform) const;

         };

         class Saddle_List : public list<Saddle>
         {

            private:

               list<Separatrix>
               separatrix_list;

               void
               grow_separatrices (const Uv_Data& uv_data,
                                  const Domain_2D& domain_2d,
                                  Sand_Box& sand_box,
                                  const Real step_size,
                                  const Integration_Scheme integration_scheme);

               void
               add_to (Sand_Box& sand_box) const;

            public:

               Saddle_List (const list<Singular_Point>& saddle_point_list,
                            const Uv_Data& uv_data,
                            const Domain_2D& domain_2d,
                            Sand_Box& sand_box,
                            const Integer size,
                            const Real h,
                            const Integration_Scheme integration_scheme);

               void
               render (const RefPtr<Context>& cr,
                       const Transform_2D& transform) const;

         };

         Saddle_List*
         saddle_list_ptr;

         list<Streamline>
         vortex_line_list;
 
         list<Streamline>
         side_w_list;
 
         list<Streamline>
         side_n_list;
 
         list<Streamline>
         side_e_list;
 
         list<Streamline>
         side_s_list;
 
         list<Streamline>
         streamline_list;
 
      public:

         const Uv_Data&
         uv_data;

         const Domain_2D&
         domain_2d;

         const Integration_Scheme
         integration_scheme;

         void
         construct_sides (Sand_Box& sand_box,
                          const Real separation,
                          const Real step_size);

         void
         construct_streamlines (Sand_Box& sand_box,
                                const Real separation,
                                const Real step_size);

         Streamlines (const Uv_Data& uv_data,
                      const Real separation = 0.5,
                      const Real step_size = 0.05,
                      const Integer saddle_size = 20,
                      const Integration_Scheme scheme = RUNGA_KUTTA);

         Streamlines (const Uv_Data& uv_data,
                      const Domain_2D& domain_2d,
                      const Real separation = 0.5,
                      const Real step_size = 0.05,
                      const Integer saddle_size = 20,
                      const Integration_Scheme scheme = RUNGA_KUTTA);

         ~Streamlines ();

         void
         render (const RefPtr<Context>& cr,
                 const Transform_2D& transform) const;

   };

}

#endif /* DENISE_STREAMLINE_H */

