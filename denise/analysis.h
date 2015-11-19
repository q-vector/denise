//
// analysis.h
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
//
// Differentiation code partially adapted from
//    Singh and Bhadauria, 2009: Finite Difference Formulae for Unequal
//       Sub-Intervals Using Lagrange's Interpolation Formula.  Int.
//       Hournal of Math. Analysis, Vol 3, 2009, 17, 815-827.

#ifndef DENISE_ANALYSIS_H
#define DENISE_ANALYSIS_H

#if HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

#include <fstream>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multiroots.h>
#include <denise/basics.h>
#include <denise/exception.h>
#include <denise/geometry.h>
#include <denise/graphics.h>
#include <denise/gzstream.h>
#include <denise/transform.h>
#include <denise/util.h>

using namespace std;

namespace denise
{

   enum Evaluate_Op_2D
   {
      MAGNITUDE_OP,
      VORTICITY_OP,
      DIVERGENCE_OP,
      SHEAR_OP,
      CURVATURE_OP
   };

   enum Integration_Scheme
   {
      EULER,
      RUNGA_KUTTA
   };

   class Differentiation
   {

      public:

         static Real
         d (const Real y_0,
            const Real y_1,
            const Real h);

         static Real
         d_0 (const Real y_0,
              const Real y_1,
              const Real y_2,
              const Real h);

         static Real
         d_1 (const Real y_0,
              const Real y_2,
              const Real h);

         static Real
         d_2 (const Real y_0,
              const Real y_1,
              const Real y_2,
              const Real h);

         static Real
         d_0 (const Real y_0,
              const Real y_1,
              const Real y_2,
              const Real x_0,
              const Real x_1,
              const Real x_2);

         static Real
         d_1 (const Real y_0,
              const Real y_1,
              const Real y_2,
              const Real x_0,
              const Real x_1,
              const Real x_2);

         static Real
         d_2 (const Real y_0,
              const Real y_1,
              const Real y_2,
              const Real x_0,
              const Real x_1,
              const Real x_2);

         static Real
         d2 (const Real y_0,
             const Real y_1,
             const Real y_2,
             const Real h);

         static Real
         d2 (const Real y_0,
             const Real y_1,
             const Real y_2,
             const Real x_0,
             const Real x_1,
             const Real x_2);

   };

   class Jacobian_2D
   {

      private:

         Real
         u_x;

         Real
         u_y;

         Real
         v_x;

         Real
         v_y;

      public:

         Jacobian_2D (const Real u_x,
                      const Real u_y,
                      const Real v_x,
                      const Real v_y);

         Jacobian_2D (const Jacobian_2D& jacobian_2d);

         Real
         get_u_x () const;

         Real
         get_u_y () const;

         Real
         get_v_x () const;

         Real
         get_v_y () const;

         const Real
         get_divergence () const;

         const Real
         get_determinant () const;

   };

   class Vector_Field_1D
   {

      public:

         virtual Real
         evaluate (const Integer vector_element,
                   const Real coordinate,
                   const Evaluate_Op evaluate_op = VALUE) const = 0;

   };

   class Vector_Field_2D
   {

      public:

         virtual Real
         evaluate (const Integer vector_element,
                   const Real coordinate_x,
                   const Real coordinate_y,
                   const Evaluate_Op evaluate_op = VALUE) const = 0;

   };

   class Vector_Field_3D
   {

      public:

         virtual Real
         evaluate (const Integer vector_element,
                   const Real coordinate_z,
                   const Real coordinate_x,
                   const Real coordinate_y,
                   const Evaluate_Op evaluate_op = VALUE) const = 0;

   };

   class Scalar_Field_1D
   {

      public:

         virtual Real
         evaluate (const Real coordinate,
                   const Evaluate_Op evaluate_op = VALUE) const = 0;

   };

   class Scalar_Field_2D
   {

      public:

         virtual Real
         evaluate (const Real coordinate_x,
                   const Real coordinate_y,
                   const Evaluate_Op evaluate_op = VALUE) const = 0;

   };

   class Scalar_Field_3D
   {

      public:

         virtual Real
         evaluate (const Real coordinate_z,
                   const Real coordinate_x,
                   const Real coordinate_y,
                   const Evaluate_Op evaluate_op = VALUE) const = 0;

   };

   class Grid_nD
   {

      protected:

         const Integer
         n;

         Tuple*
         coordinate_tuples;

         Real*
         spacings;

         bool*
         periodics;

         Size_nD
         size_nd;

         void
         init () throw (std::bad_alloc);

         void
         init (const Tuple* coordinate_tuples,
               const Real* spacings,
               const bool* periodics);

      public:

         Grid_nD (const Integer n);

         Grid_nD (const Grid_nD& grid_nd);

         ~Grid_nD ();

         Integer
         size () const;

         Tuple&
         get_coordinate_tuple (const Integer dimension) const;

         const Real&
         get_coordinate (const Integer dimension,
                         const Integer node) const;

         static Integer
         get_node (const Tuple& coordinate_tuple,
                   const Real spacing,
                   const Real coordinate);

         Integer
         get_node (const Integer dimension,
                   const Real coordinate) const;

         static Integer
         get_nearest_node (const Tuple& coordinate_tuple,
                           const Real spacing,
                           const Real coordinate);

         Integer
         get_nearest_node (const Integer dimension,
                           const Real coordinate) const;

         const Real&
         get_spacing (const Integer dimension) const;

         Real
         get_spacing (const Integer dimension,
                      const Integer node) const;

         bool
         node_out_of_bounds (const Integer dimension,
                             const Integer node) const;

         bool
         out_of_bounds (const Integer dimension,
                        const Real coordinate) const;

         void
         translate (const Integer dimension,
                    const Real delta);
         
         void
         standardize_node (const Integer dimension,
                           Integer& node) const;

         void
         standardize_coordinate (const Integer dimension,
                                 Real& coordinate) const;

   };

   class Chunk
   {

      protected:

         size_t
         chunk_size;

         Real*
         buffer;

      public:

         Chunk ();

         Chunk (const Chunk& chunk);

         ~Chunk ();

         void
         read (igzstream& file,
               const bool float_length = true);

         void
         write (ogzstream& file,
                const bool float_length = true) const;

         void
         init (const Integer chunk_size);

         void
         copy (const Chunk& chunk,
               const Integer address);

         void
         set (const Integer i,
              const Real datum);

         const Real&
         get (const Integer i) const;

         Real&
         get (const Integer i);

         void
         initialize (const Real datum,
                     const size_t address = 0,
                     const size_t size = -1);

         void
         scale_offset (const Real scale,
                       const Real offset,
                       const size_t address = 0,
                       const size_t size = -1);

         Domain_1D
         get_max_min (const size_t address = 0,
                      const size_t size = -1) const;

         Real
         get_mean (const size_t address = 0,
                   const size_t size = -1) const;

         Real
         subtract_mean (const size_t address = 0,
                        const size_t size = -1);

   };

   class Vector_Data_nD : public Grid_nD,
                          public Chunk
   {

      protected:

         const Integer
         vector_size;

         void
         init (const Tuple* coordinate_tuples,
               const Real* spacings,
               const bool* periodics);

         Vector_Data_nD (const Integer vector_size,
                         const Integer n);

         Vector_Data_nD (const Vector_Data_nD& vector_data_nd);

         ~Vector_Data_nD ();

         Chunk*
         get_chunk_ptr (const Integer vector_element) const;

      public:

         const Integer&
         get_vector_size () const;

         void
         initialize_all (const Real datum);

         virtual void
         initialize (const Integer vector_element,
                     const Real datum) = 0;

         virtual void
         scale_offset (const Integer vector_element,
                       const Real scale,
                       const Real offset) = 0;

         virtual Domain_1D
         get_max_min (const Integer vector_element) const = 0;

         Real
         get_epsilon (const Integer vector_element) const;

         void
         subtract_mean (const Integer vector_element);

   };

   class Vector_Data_1D : public Vector_Data_nD,
                          public virtual Vector_Field_1D
   {

      private:

         gsl_interp_accel**
         gia_ptrs;

         gsl_spline**
         gs_ptrs;

         void
         init (const Tuple& coordinate_tuple,
               const Real spacing,
               const bool periodic);

      public:

         Vector_Data_1D (const Vector_Data_1D& vector_data_1d);

         Vector_Data_1D (const Integer vector_size,
                         const Integer size_1d,
                         const Domain_1D& domain_1d = Domain_1D (0, 1),
                         const bool periodic = false);

         Vector_Data_1D (const Integer vector_size,
                         const Tuple coordinate_tuple,
                         const bool periodic = false);

         ~Vector_Data_1D ();

         void
         modify_coordinate_tuple (const Integer node,
                                  const Real coordinate);

         void
         set_interpolation (const gsl_interp_type* interp_type = gsl_interp_linear);

         void
         set_interpolation (const Integer vector_element,
                            const gsl_interp_type* interp_type = gsl_interp_linear);

         void
         initialize (const Integer vector_element,
                     const Real datum);

         void
         scale_offset (const Integer vector_element,
                       const Real scale,
                       const Real offset);

         Domain_1D
         get_max_min (const Integer vector_element) const;

         void
         set_datum (const Integer vector_element,
                    const Integer node,
                    const Real datum);

         const Real&
         get_datum (const Integer vector_element,
                    const Integer node) const;

         Real&
         get_datum (const Integer vector_element,
                    const Integer node);

         Real
         evaluate (const Integer vector_element,
                   const Real coordinate,
                   const Evaluate_Op evaluate_op = VALUE) const;

   };

   class Bicubic_Coefficients
   {

      private:

         Real***
         data;

         const Size_2D
         size_2d;

      public:

         Bicubic_Coefficients (const Size_2D& size_2d);

         ~Bicubic_Coefficients ();

         Real*
         get_array (const Integer node_x,
                    const Integer node_y) const;

   };

   class Tricubic_Coefficients
   {

      private:

         Real****
         data;

         const Size_3D
         size_3d;

      public:

         Tricubic_Coefficients (const Size_3D& size_3d);

         ~Tricubic_Coefficients ();

         Real*
         get_array (const Integer node_z,
                    const Integer node_x,
                    const Integer node_y) const;

   };

   class Vector_Data_2D : public Vector_Data_nD,
                          public virtual Vector_Field_2D
   {

      private:

         Bicubic_Coefficients**
         bicubic_coefficients_ptrs;

         void
         init_bicubic_coefficients ();

         Bicubic_Coefficients*
         get_bicubic_coefficients_ptr (const Chunk& chunk,
                                       const Chunk& chunk_x,
                                       const Chunk& chunk_y,
                                       const Chunk& chunk_xy) const;

         static int
         gsl_multiroot_f (const gsl_vector* x,
                          void* params,
                          gsl_vector* f);

         static int
         gsl_multiroot_df (const gsl_vector* x,
                           void* params,
                           gsl_matrix* J);

         static int
         gsl_multiroot_fdf (const gsl_vector* x,
                            void* params,
                            gsl_vector* f,
                            gsl_matrix* J);

      protected:

         void
         init (const Tuple& coordinate_tuple_x,
               const Tuple& coordinate_tuple_y,
               const Real spacing_x,
               const Real spacing_y,
               const bool periodic_x,
               const bool periodic_y);

         Chunk*
         get_chunk_x_ptr (const Chunk& chunk);

         Chunk*
         get_chunk_y_ptr (const Chunk& chunk);

         Real
         evaluate_bilinear (const Integer vector_element,
                            const Real coordinate_x,
                            const Real coordinate_y,
                            const Evaluate_Op evaluate_op) const;

         Real
         evaluate_bicubic (const Integer vector_element,
                           const Real coordinate_x,
                           const Real coordinate_y,
                           const Evaluate_Op evaluate_op) const;

         Real
         evaluate_bicubic (const Real* a,
                           const Real u,
                           const Real v,
                           const Real dx,
                           const Real dy,
                           const Evaluate_Op evaluate_op) const;

         virtual Real
         get_dmagnitude_dx (const Integer vector_element_u,
                            const Integer vector_element_v,
                            const Integer node_x,
                            const Integer node_y,
                            const Real magnitude) const;

         virtual Real
         get_dmagnitude_dy (const Integer vector_element_u,
                            const Integer vector_element_v,
                            const Integer node_x,
                            const Integer node_y,
                            const Real magnitude) const;

         virtual Real
         get_dmagnitude_dx (const Integer vector_element_u,
                            const Integer vector_element_v,
                            const Real x,
                            const Real y,
                            const Real magnitude) const;

         virtual Real
         get_dmagnitude_dy (const Integer vector_element_u,
                            const Integer vector_element_v,
                            const Real x,
                            const Real y,
                            const Real magnitude) const;

      public:

         Vector_Data_2D (const Vector_Data_2D& vector_data_2d,
                         const bool copy_data);

         Vector_Data_2D (const Integer vector_size,
                         const Size_2D& size_2d,
                         const Domain_2D& domain_2d,
                         const bool periodic_x = false,
                         const bool periodic_y = false);

         Vector_Data_2D (const Integer vector_size,
                         const Tuple coordinate_tuple_x,
                         const Tuple coordinate_tuple_y,
                         const bool periodic_x = false,
                         const bool periodic_y = false);

         virtual
         ~Vector_Data_2D ();

         void
         set_bicubic_interpolation ();

         void
         set_bicubic_interpolation (const Integer vector_element);

         bool
         out_of_bounds (const Integer node_x,
                        const Integer node_y) const;

         bool
         out_of_bounds (const Real coordinate_x,
                        const Real coordinate_y) const;

         Size_2D
         get_size_2d () const;

         Domain_2D
         get_domain_2d () const;

         void
         initialize (const Integer vector_element,
                     const Real datum);

         void
         scale_offset (const Integer vector_element,
                       const Real scale,
                       const Real offset);

         Domain_1D
         get_max_min (const Integer vector_element) const;

         void
         set_datum (const Integer vector_element,
                    const Integer node_x,
                    const Integer node_y,
                    const Real datum);

         const Real&
         get_datum (const Integer vector_element,
                    const Integer node_x,
                    const Integer node_y) const;

         Real&
         get_datum (const Integer vector_element,
                    const Integer node_x,
                    const Integer node_y);

         virtual Real
         evaluate (const Integer vector_element,
                   const Integer node_x,
                   const Integer node_y,
                   const Evaluate_Op evaluate_op = VALUE) const;

         virtual Real
         evaluate (const Integer vector_element,
                   const Real coordinate_x,
                   const Real coordinate_y,
                   const Evaluate_Op evaluate_op = VALUE) const;

         virtual Real
         evaluate_2d (const Integer vector_element_u,
                      const Integer vector_element_v,
                      const Integer node_x,
                      const Integer node_y,
                      const Evaluate_Op_2D evaluate_op_2d) const;

         virtual Real
         evaluate_2d (const Integer vector_element_u,
                      const Integer vector_element_v,
                      const Real coordinate_x,
                      const Real coordinate_y,
                      const Evaluate_Op_2D evaluate_op_2d) const;

         Real
         evaluate_nocheck (const Integer vector_element,
                           const Integer node_x,
                           const Integer node_y,
                           const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate_nocheck (const Integer vector_element,
                           const Real coordinate_x,
                           const Real coordinate_y,
                           const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate_nocheck (const Integer vector_element,
                           const Point_2D& point_2d,
                           const Evaluate_Op evaluate_op = VALUE) const;

         void
         acquire_root (Real& root_coordinate_x,
                       Real& root_coordinate_y,
                       const Real residual = 1e-9) const;

         Jacobian_2D
         get_jacobian (const Integer u_index,
                       const Integer v_index,
                       const Point_2D& point_2d) const;

         Jacobian_2D
         get_jacobian_nocheck (const Integer u_index,
                               const Integer v_index,
                               const Point_2D& point_2d) const;

         void
         subtract_x_mean (const Integer vector_element);

         void
         subtract_y_mean (const Integer vector_element);

   };

   class Vector_Data_3D : public Vector_Data_nD,
                          public virtual Vector_Field_3D
   {

      private:

         Tricubic_Coefficients**
         tricubic_coefficients_ptrs;

         Tricubic_Coefficients*
         get_tricubic_coefficients_ptr (const Chunk& chunk,
                                        const Chunk& chunk_y,
                                        const Chunk& chunk_x,
                                        const Chunk& chunk_xy,
                                        const Chunk& chunk_z,
                                        const Chunk& chunk_zy,
                                        const Chunk& chunk_zx,
                                        const Chunk& chunk_zxy) const;

      protected:

         Integer
         get_offset (const Integer vector_element,
                     const Integer node_z,
                     const Integer node_x,
                     const Integer node_y) const;

         void
         init (const Tuple& coordinate_tuple_z,
               const Tuple& coordinate_tuple_x,
               const Tuple& coordinate_tuple_y,
               const Real spacing_z,
               const Real spacing_x,
               const Real spacing_y,
               const bool periodic_z,
               const bool periodic_x,
               const bool periodic_y);

         Chunk*
         get_chunk_z_ptr (const Chunk& chunk);

         Chunk*
         get_chunk_x_ptr (const Chunk& chunk);

         Chunk*
         get_chunk_y_ptr (const Chunk& chunk);

         Real
         evaluate_linear_z (const Integer vector_element,
                            const Real coordinate_z,
                            const Integer node_x,
                            const Integer node_y) const;

         Real
         evaluate_bilinear (const Integer vector_element,
                            const Integer node_z,
                            const Real coordinate_x,
                            const Real coordinate_y,
                            const Evaluate_Op evaluate_op) const;

         Real
         evaluate_trilinear (const Integer vector_element,
                             const Real coordinate_z,
                             const Real coordinate_x,
                             const Real coordinate_y,
                             const Evaluate_Op evaluate_op) const;

         Real
         evaluate_tricubic (const Integer vector_element,
                            const Real coordinate_z,
                            const Real coordinate_x,
                            const Real coordinate_y,
                            const Evaluate_Op evaluate_op) const;

         Real
         evaluate_tricubic (const Real* a,
                            const Real w,
                            const Real u,
                            const Real v,
                            const Real dz,
                            const Real dx,
                            const Real dy,
                            const Evaluate_Op evaluate_op) const;

      public:

         Vector_Data_3D (const Integer vector_size,
                         const Size_3D& size_3d,
                         const Domain_3D& domain_3d,
                         const bool periodic_z = false,
                         const bool periodic_x = false,
                         const bool periodic_y = false);

         Vector_Data_3D (const Integer vector_size,
                         const Tuple coordinate_tuple_z,
                         const Tuple coordinate_tuple_x,
                         const Tuple coordinate_tuple_y,
                         const bool periodic_z = false,
                         const bool periodic_x = false,
                         const bool periodic_y = false);

         Vector_Data_3D (const Integer vector_size,
                         const Tuple coordinate_tuple_z,
                         const Size_2D& size_2d,
                         const Domain_2D& domain_2d,
                         const bool periodic_z = false,
                         const bool periodic_x = false,
                         const bool periodic_y = false);

         virtual
         ~Vector_Data_3D ();

         void
         set_tricubic_interpolation ();

         void
         set_tricubic_interpolation (const Integer vector_element);

         bool
         out_of_bounds (const Integer node_z,
                        const Integer node_x,
                        const Integer node_y) const;

         bool
         out_of_bounds (const Integer node_z,
                        const Real coordinate_x,
                        const Real coordinate_y) const;

         bool
         out_of_bounds (const Real coordinate_z,
                        const Integer node_x,
                        const Integer node_y) const;

         bool
         out_of_bounds (const Real coordinate_z,
                        const Real coordinate_x,
                        const Real coordinate_y) const;

         Size_2D
         get_size_2d () const;

         Domain_2D
         get_domain_2d () const;

         Size_3D
         get_size_3d () const;

         void
         initialize (const Integer vector_element,
                     const Real datum);

         void
         scale_offset (const Integer vector_element,
                       const Real scale,
                       const Real offset);

         Domain_1D
         get_max_min (const Integer vector_element) const;

         void
         set_datum (const Integer vector_element,
                    const Integer k,
                    const Integer i,
                    const Integer j,
                    const Real datum);

         const Real&
         get_datum (const Integer vector_element,
                    const Integer k,
                    const Integer i,
                    const Integer j) const;

         Real&
         get_datum (const Integer vector_element,
                    const Integer k,
                    const Integer i,
                    const Integer j);

         virtual Real
         evaluate (const Integer vector_element,
                   const Integer k,
                   const Integer i,
                   const Integer j,
                   const Evaluate_Op evaluate_op = VALUE) const;

         virtual Real
         evaluate (const Integer vector_element,
                   const Integer k,
                   const Real x,
                   const Real y,
                   const Evaluate_Op evaluate_op = VALUE) const;

         virtual Real
         evaluate (const Integer vector_element,
                   const Real z,
                   const Integer i,
                   const Integer j,
                   const Evaluate_Op evaluate_op = VALUE) const;

         virtual Real
         evaluate (const Integer vector_element,
                   const Real z,
                   const Real x,
                   const Real y,
                   const Evaluate_Op evaluate_op = VALUE) const;

         virtual Real
         evaluate_uv (const Integer vector_element_u,
                      const Integer vector_element_v,
                      const Integer k,
                      const Integer i,
                      const Integer j,
                      const Evaluate_Op_2D evaluate_op_2d) const;

         virtual Real
         evaluate_uv (const Integer vector_element_u,
                      const Integer vector_element_v,
                      const Integer k,
                      const Real x,
                      const Real y,
                      const Evaluate_Op_2D evaluate_op_2d) const;

         virtual Real
         evaluate_uv (const Integer vector_element_u,
                      const Integer vector_element_v,
                      const Real z,
                      const Integer i,
                      const Integer j,
                      const Evaluate_Op_2D evaluate_op_2d) const;

         virtual Real
         evaluate_uv (const Integer vector_element_u,
                      const Integer vector_element_v,
                      const Real z,
                      const Real x,
                      const Real y,
                      const Evaluate_Op_2D evaluate_op_2d) const;

         Real
         evaluate_nocheck (const Integer vector_element,
                           const Integer z,
                           const Integer i,
                           const Integer j,
                           const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate_nocheck (const Integer vector_element,
                           const Integer z,
                           const Real x,
                           const Real y,
                           const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate_nocheck (const Integer vector_element,
                           const Real z,
                           const Integer i,
                           const Integer j,
                           const Evaluate_Op evaluate_op = VALUE) const;

         Real
         evaluate_nocheck (const Integer vector_element,
                           const Real z,
                           const Real x,
                           const Real y,
                           const Evaluate_Op evaluate_op = VALUE) const;

         void
         subtract_z_mean (const Integer vector_element);

         void
         subtract_x_mean (const Integer vector_element);

         void
         subtract_y_mean (const Integer vector_element);

         void
         subtract_zx_mean (const Integer vector_element);

         void
         subtract_zy_mean (const Integer vector_element);

         void
         subtract_xy_mean (const Integer vector_element);

   };

   class Scalar_Data_1D : public virtual Vector_Data_1D,
                          public Scalar_Field_1D
   {

      public:

         Scalar_Data_1D (const Integer size_1d,
                         const Domain_1D& domain_1d,
                         const bool periodic = false);

         Scalar_Data_1D (const Tuple coordinate_tuple,
                         const bool periodic = false);

         void
         set_datum (const Integer node,
                    const Real datum);

         const Real&
         get_datum (const Integer node) const;

         Real&
         get_datum (const Integer node);

         Real
         evaluate (const Real coordinate,
                   const Evaluate_Op evaluate_op = VALUE) const;

   };

   class Scalar_Data_2D : public virtual Vector_Data_2D,
                          public Scalar_Field_2D
   {

      public:

         Scalar_Data_2D (const Size_2D& size_2d,
                         const Domain_2D& domain_2d,
                         const bool periodic_x = false,
                         const bool periodic_y = false);

         Scalar_Data_2D (const Tuple coordinate_tuple_x,
                         const Tuple coordinate_tuple_y,
                         const bool periodic_x = false,
                         const bool periodic_y = false);

         Scalar_Data_2D (const Vector_Data_2D& vector_data_2d,
                         const Integer vector_element = -1);

         Domain_1D
         get_max_min () const;

         Domain_1D
         get_max_min (const Domain_2D& domain_2d) const;

         void
         set_datum (const Integer node_x,
                    const Integer node_y,
                    const Real datum);

         const Real&
         get_datum (const Integer node_x,
                    const Integer node_y) const;

         Real&
         get_datum (const Integer node_x,
                    const Integer node_y);

         virtual Real
         evaluate (const Real coordinate_x,
                   const Real coordinate_y,
                   const Evaluate_Op evaluate_op = VALUE) const;

   };

   class Scalar_Data_3D : public virtual Vector_Data_3D,
                          public Scalar_Field_3D
   {

      public:

         Scalar_Data_3D (const Size_3D& size_3d,
                         const Domain_3D& domain_3d,
                         const bool periodic_z = false,
                         const bool periodic_x = false,
                         const bool periodic_y = false);

         Scalar_Data_3D (const Tuple coordinate_tuple_z,
                         const Tuple coordinate_tuple_x,
                         const Tuple coordinate_tuple_y,
                         const bool periodic_z = false,
                         const bool periodic_x = false,
                         const bool periodic_y = false);

         Scalar_Data_3D (const Tuple coordinate_tuple_z,
                         const Size_2D& size_2d,
                         const Domain_2D& domain_2d,
                         const bool periodic_z = false,
                         const bool periodic_x = false,
                         const bool periodic_y = false);

         void
         set_datum (const Integer node_z,
                    const Integer node_x,
                    const Integer node_y,
                    const Real datum);

         const Real&
         get_datum (const Integer node_z,
                    const Integer node_x,
                    const Integer node_y) const;

         Real&
         get_datum (const Integer node_z,
                    const Integer node_x,
                    const Integer node_y);

         virtual Real
         evaluate (const Real coordinate_z,
                   const Real coordinate_x,
                   const Real coordinate_y,
                   const Evaluate_Op evaluate_op = VALUE) const;

         virtual Real
         evaluate (const Integer node_z,
                   const Real coordinate_x,
                   const Real coordinate_y,
                   const Evaluate_Op evaluate_op = VALUE) const;

   };

   class Uv_Field
   {

      private:

         const Vector_Data_2D*
         vector_data_2d_ptr;

         const Integer
         u_index;

         const Integer
         v_index;

         void
         step (Point_2D& rk,
               const Real x,
               const Real y,
               const Real h) const;

      public:

         Uv_Field ();

         Uv_Field (const Vector_Data_2D& vector_data_2d,
                   const Integer u_index,
                   const Integer v_index);

         virtual Vector_2D
         evaluate (const Real x,
                   const Real y) const;

         void
         integrate (Point_2D& next_point_2d,
                    const Point_2D& point_2d,
                    const Real h,
                    const Integration_Scheme integration_scheme) const;

   };

   class Bezier_Curve : public vector<Tuple>
   {

      private:

         Integer
         tuple_size;

         Real
         bernstein (const Real t,
                    const Integer k,
                    const Integer n) const;

      public:

         Bezier_Curve (const vector<Tuple>& tuple_vector);

         Real
         get_value_at (const Real t,
                       const Integer tuple_index) const;

         Tuple
         get_value_at (const Real t) const;

   };

   class B_Spline : protected vector<Tuple>
   {
   
      private:

         Integer
         n;
   
         Integer
         d;
   
         Tuple
         knot_tuple;
   
         Real
         cox_deboor (const Real t, 
                     const Integer k, 
                     const Integer d) const;
   
      public:
   
         B_Spline (const vector<Tuple>& tuple_vector,
                   const Integer degree = 2,
                   const bool open = false,
                   const bool repeat = false);
   
         B_Spline (const vector<Tuple>& tuple_vector,
                   const Tuple& knot_tuple);
   
         Real
         get_value_at (const Real t,
                       const Integer tuple_index) const;
        
         Tuple
         get_value_at (const Real t) const;
      
   }; 

   class Cartesian_Transform_1D : public Affine_Transform_1D
   {

      private:

         const bool
         log;

         void
         init (const Domain_1D& domain,
               const Domain_1D& range);

      public:

         Cartesian_Transform_1D (const Domain_1D& domain,
                                 const Domain_1D& range,
                                 const bool log = false);

         void
         transform (Real& transformed,
                    const Real x) const;

         Real
         transform (const Real x) const;

         Real
         reverse (const Real x) const;

         void
         reverse (Real& reversed,
                  const Real x) const;

   };

   class Cartesian_Transform_2D : public Affine_Transform_2D
   {

      private:

         void
         init (const Domain_1D& domain_x,
               const Domain_1D& domain_y,
               const Real width,
               const Real height,
               const Point_2D& origin,
               const bool transpose);

      public:

         Cartesian_Transform_2D (const Domain_1D& domain_x,
                                 const Domain_1D& domain_y,
                                 const Box_2D& box_2d,
                                 const bool transpose = false);

         Cartesian_Transform_2D (const Domain_1D& domain_x,
                                 const Domain_1D& domain_y,
                                 const Real width,
                                 const Real height,
                                 const Point_2D& origin = Point_2D (0, 0),
                                 const bool transpose = false);

   };

   class Log_Transform_2D : public Cartesian_Transform_2D
   {

      private:

         bool
         log_x;

         bool
         log_y;

      public:

         Log_Transform_2D (const Domain_1D& domain_x,
                           const Domain_1D& domain_y,
                           const Real width,
                           const Real height,
                           const bool log_x,
                           const bool log_y,
                           const Point_2D& origin = Point_2D (0, 0));

         bool
         out_of_domain (const Real x,
                        const Real y) const;

         void
         transform (Real& transformed_x,
                    Real& transformed_y,
                    const Real x,
                    const Real y) const;

         void
         reverse (Real& reversed_x,
                  Real& reversed_y,
                  const Real x,
                  const Real y) const;

         void
         transform_uv (Real& u,
                       Real& v,
                       const Real x,
                       const Real y) const;

         Real
         get_theta (const Real u,
                    const Real v,
                    const Real x,
                    const Real y) const;

   };

   class Polar_Transform_2D : public Transform_2D
   {
   
      private:

         Point_2D
         origin;

         Real
         scale;

      public:

         Polar_Transform_2D ();

         Polar_Transform_2D (const Point_2D& origin,
                             const Real scale);

         void
         set (const Point_2D& origin,
              const Real scale);

         virtual bool     
         out_of_domain (const Real r,
                        const Real theta) const;

         virtual void
         transform (Real& x,
                    Real& y,
                    const Real r,
                    const Real theta) const;

         virtual void
         reverse (Real& r,
                  Real& theta,
                  const Real x,
                  const Real y) const;

   };

   class Parabolic_Transform_2D : public Transform_2D
   {

      private:

         const Point_2D
         origin;

      public:

         Parabolic_Transform_2D (const Point_2D& origin);

         bool
         out_of_domain (const Real u,
                        const Real v) const;

         void
         transform (Real& x,
                    Real& y,
                    const Real u,
                    const Real v) const;

         void
         reverse (Real& u,
                  Real& v,
                  const Real x,
                  const Real y) const;

   };

   class Elliptic_Transform_2D : public Transform_2D
   {

      private:

         const Point_2D
         origin;

         const Real
         scale;

      public:

         Elliptic_Transform_2D (const Point_2D& origin,
                                const Real scale);

         bool
         out_of_domain (const Real u,
                        const Real v) const;

         void
         transform (Real& x,
                    Real& y,
                    const Real u,
                    const Real v) const;

         void
         reverse (Real& u,
                  Real& v,
                  const Real x,
                  const Real y) const;

   };

   class Bipolar_Transform_2D : public Transform_2D
   {

      private:

         const Point_2D
         origin;

         const Real
         scale;

      public:

         Bipolar_Transform_2D (const Point_2D& origin,
                               const Real scale);

         bool
         out_of_domain (const Real u,
                        const Real v) const;

         void
         transform (Real& x,
                    Real& y,
                    const Real u,
                    const Real v) const;

         void
         reverse (Real& u,
                  Real& v,
                  const Real x,
                  const Real y) const;

   };

   class Poisson_Disk : public vector<Point_2D>
   {

      private:

         class Process_Vector : public vector<Point_2D>
         {

            public:

               void
               push (const Point_2D& point);

               Point_2D
               pop ();

         };

         class Grid : public Vector_Data_2D
         {

            private:

               const Real
               r;

               const Real
               h;

            public:

               Grid (const Size_2D& size_2d,
                     const Domain_2D& domain_2d,
                     const Real r,
                     const Real h);

               void
               set_point (const Point_2D& point);

               bool
               has_neighbour (const Point_2D& point);

         };

         Point_2D
         get_random_nearby_point (const Point_2D& p,
                                  const Real r) const;

      public:

         Poisson_Disk (const Domain_2D& domain_2d,
                       const Real r);

   };

   class Graph_Raster : public Raster
   {

      protected:

         Real
         f (const Real x,
            const Real y) const;

      public:

         Graph_Raster (const Transform_2D& transform,
                       const Box_2D& box_2d,
                       const Color& color,
                       const Real width);

   };

}

#endif /* DENISE_ANALYSIS_H */

