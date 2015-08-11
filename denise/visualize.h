//
// visualize.h
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

#ifndef DENISE_VISUALIZE_H
#define DENISE_VISUALIZE_H

#include <gsl/gsl_poly.h>
#include <cairomm/context.h>
#include <denise/analysis.h>
#include <denise/exception.h>
#include <denise/geodesy.h>
#include <denise/geometry.h>
#include <denise/graphics.h>

using namespace std;

namespace denise
{

   enum Scalarization_2d
   {
      MAGNITUDE,
      DIVERGENCE,
      VORTICITY
   };

   class Contour
   {

      public:

         class Isoline;

         enum Cell_State
         {
            UNPROCESSED,
            SEMI_PROCESSED_NE,
            SEMI_PROCESSED_NW,
            SEMI_PROCESSED_SE,
            SEMI_PROCESSED_SW,
            PROCESSED
         };

         enum Step_Type
         {
            TURN_LEFT,
            TURN_RIGHT,
            GO_FORWARD
         };

         enum Side
         {
            SIDE_N,
            SIDE_S,
            SIDE_E,
            SIDE_W
         };

         enum Side_End
         {
            SE_LEFT,
            SE_RIGHT,
            SE_LOWER,
            SE_UPPER
         };

         class Label_Point_Set : set<Point_2D>
         {

            private:

               const Real
               tolerance;

            public:

               Label_Point_Set (const Real tolerance);

               void
               add (const Point_2D& point);

               bool
               interfere_with (const Point_2D& point);

         };

         class Boundary
         {

            public:

               class Segment
               {

                  //private:
                  public:

                     Integer
                     index;

                     Real
                     value_a;

                     Real
                     value_b;

                     Point_2D
                     point_a;

                     Point_2D
                     point_b;

                     Isoline**
                     isoline_ptrs;

                     bool*
                     heads;

                  public:

                     Segment (const Integer index,
                              const Real value_a,
                              const Real value_b,
                              const Point_2D& point_a,
                              const Point_2D& point_b,
                              const Integer number_of_levels);

                     ~Segment ();

                     Isoline*
                     get_isoline_ptr (const Integer level_index) const;

                     bool
                     connected_to_head (const Integer level_index) const;

                     bool
                     connected_to_head (const Isoline& isoline) const;

                     void
                     register_isoline (const Isoline& isoline,
                                       const bool head) const;

                     void
                     detach_isoline (const Integer level_index);

                     const Integer&
                     get_index () const;

                     const Real&
                     get_value_a () const;

                     const Real&
                     get_value_b () const;

                     const Point_2D&
                     get_point_a () const;

                     const Point_2D&
                     get_point_b () const;

               };

            private:

               vector<Segment*>
               segment_ptr_vector;

            public:

               Boundary (const Scalar_Data_2D& scalar_data_2d,
                         const Tuple& level_tuple);

               ~Boundary ();

               Integer
               size () const;

               void
               set_segment_ptr (const Segment* segment_ptr,
                                const Integer index);

               Boundary::Segment*
               get_segment_ptr (const Integer index) const;

               void
               advance_segment_ptr (Segment*& segment_ptr,
                                    const bool backward) const;

         };

         class Edge : public denise::Edge
         {

            public:

               const Integer
               level_index;

               const Side
               side_a;

               const Side
               side_b;

               Edge (const Integer level_index,
                     const Point_2D& point_a,
                     const Point_2D& point_b,
                     const Side side_a,
                     const Side side_b);

               bool
               operator== (const Edge& edge) const;

         };

         class Side_Edge
         {

            public:

               const Integer
               lower_level_index;

               const Side
               side;

               const Real
               position_a;

               const Real
               position_b;

               const Side_End
               side_end_a;

               const Side_End
               side_end_b;

               Side_Edge (const Integer lower_level_index,
                          const Real position_a,
                          const Real position_b,
                          const Side side,
                          const Side_End side_end_a,
                          const Side_End side_end_b);

               bool
               operator== (const Side_Edge& side_edge) const;

         };

         class Cell : public Index_2D
         {

            private:

               list<Edge>*
               edge_lists;

               list<Side_Edge>*
               side_edge_lists;

               list<Isoline*>*
               start_isoline_ptr_lists;

               list<Isoline*>*
               end_isoline_ptr_lists;

            public:

               const Tuple&
               level_tuple;

               Cell (const Tuple& level_tuple,
                     const Integer i,
                     const Integer j);

               ~Cell ();

               bool
               has_edge (const Integer level_index) const;

               bool
               has_edge (const Integer level_index,
                         const Side side) const;

               Edge
               get_edge (const Integer level_index) const;

               Edge
               get_edge (const Integer level_index,
                         const Side side) const;

               bool
               has_side_edge (const Integer lower_level_index) const;

               Side_Edge
               get_side_edge (const Integer lower_level_index) const;

               void
               insert (const Edge& edge);

               void
               insert (const Side_Edge& side_edge);

               void
               register_start_isoline (Isoline& isoline);

               void
               register_end_isoline (Isoline& isoline);

         };

         class Isoline : public Simple_Polyline
         {

            //private:
            public:

               Boundary::Segment*
               segment_ptr_a;

               Boundary::Segment*
               segment_ptr_b;

               bool
               soiled;

            public:

               const Integer
               level_index;

               Isoline (const Integer level_index);

               Boundary::Segment*
               get_segment_ptr_a () const;

               Boundary::Segment*
               get_segment_ptr_b () const;

               void
               set_segment_ptr_a (Boundary::Segment* segment_ptr);

               void
               set_segment_ptr_b (Boundary::Segment* segment_ptr);

               bool
               is_closed () const;

               const bool
               is_soiled () const;

               void
               mark_soiled ();

               void
               reset_soiled ();

         };

      private:

         class Valid_Edge : public Index_2D
         {

            public:

               bool
               along_x;

               Valid_Edge () {};

               Valid_Edge (const Integer i,
                           const Integer j,
                           const bool along_x);

               bool
               operator== (const Valid_Edge& valid_edge) const;

               bool
               operator> (const Valid_Edge& valid_edge) const;

               bool
               operator< (const Valid_Edge& valid_edge) const;

         };

         class Valid_Edge_Set : public set<Valid_Edge>
         {

            public:

               bool
               contains (const Valid_Edge& valid_edge) const;

               bool
               contains (const Integer i,
                         const Integer j,
                         const bool along_x) const;

         };

         class Valid_Data
         {

            private:

               Size_2D
               size_2d;

               bool*
               buffer;

            public:

               Valid_Data (const Scalar_Data_2D& scalar_data_2d);

               ~Valid_Data ();

               const Size_2D&
               get_size_2d () const;

               bool
               get (const Integer i,
                    const Integer j) const;

               void
               set (const Integer i,
                    const Integer j,
                    const bool b) const;

               bool
               is_valid_edge_x (const Integer i,
                                const Integer j) const;

               bool
               is_valid_edge_y (const Integer i,
                                const Integer j) const;

               bool
               is_valid (const Valid_Edge& valid_edge) const;

         };

         class Valid_Path : public Valid_Edge
         {

            public:

               bool
               forward;

               Valid_Path () {};

               Valid_Path (const Integer i,
                           const Integer j,
                           const bool along_x,
                           const bool forward);

               Valid_Path (const Valid_Edge& valid_edge,
                           const bool forward);

               void
               set (const Integer i,
                    const Integer j,
                    const bool along_x,
                    const bool forward);

               bool
               should_start_trace (const Valid_Data& valid_data,
                                   const Valid_Edge_Set& valid_edge_set) const;

         };

         class Valid_Polygon : public Polygon
         {

            private:

               Valid_Edge_Set
               valid_edge_set;

               Valid_Path
               get_valid_path_a (const Valid_Path& valid_path) const;

               Valid_Path
               get_valid_path_b (const Valid_Path& valid_path) const;

               Valid_Path
               get_valid_path_c (const Valid_Path& valid_path) const;

               Point_2D
               get_point_head (const Scalar_Data_2D& scalar_data_2d,
                               const Valid_Path& valid_path) const;

               Point_2D
               get_point_tail (const Scalar_Data_2D& scalar_data_2d,
                               const Valid_Path& valid_path) const;

               void
               trace (const Scalar_Data_2D& scalar_data_2d,
                      const Valid_Data& valid_data,
                      const Valid_Path& valid_path,
                      const bool new_handle = false);

            public:

               Valid_Polygon (const Scalar_Data_2D& scalar_data_2d);

         };

         Valid_Polygon*
         valid_polygon_ptr;

      //private:
      public:

         Tuple
         level_tuple;

         Real
         max_value;

         Real
         min_value;

         Scalar_Data_2D*
         scalar_data_2d_ptr;

         Cell**
         cell_ptrs;

         Boundary*
         boundary_ptr;

         list<Isoline*>*
         isoline_ptr_lists;

         Polygon**
         polygon_ptrs;

         static Side
         opposite (const Side side);

         Cell&
         get_cell (const Integer i,
                   const Integer j) const; 

         Cell**
         get_cell_ptrs () const;

         void
         substitude_nan (const Real substitude);

         void
         work_out_edges ();

         void
         work_out_side_edges ();

         void
         trace_isolines ();

         void
         work_out_edges (const Integer level_index);

         void
         work_out_side_edges (const Integer lower_level_index,
                              const Integer upper_level_index);

         void
         work_out_side_edge (Cell& cell,
                             const Side side,
                             const Real left,
                             const Real right,
                             const Real l_value,
                             const Real r_value,
                             const Integer lower_level_index,
                             const Integer upper_level_index);

         void
         work_out_polygon (const Integer lower_level_index);

         Boundary::Segment*
         follow_isoline (Polygon& polygon,
                         const Isoline& isoline,
                         const bool forward,
                         const bool new_handle) const;

         void
         trace_isolines (const Integer level_index);

         bool
         trace_isoline (Isoline& isoline,
                        const Edge& edge,
                        const Integer i,
                        const Integer j,
                        const Side side,
                        const bool forward);

         Polygon*
         render_label (const RefPtr<Context>& cr,
                       const Transform_2D& transform,
                       const Dstring& format,
                       const Real label_multiplier,
                       const Real label_offset,
                       const Real label_distance,
                       const Integer label_stride) const;

         void
         init_a (const Vector_Field_2D& vector_field_2d,
                 const Integer vector_element,
                 const Tuple& coordinate_tuple_x,
                 const Tuple& coordinate_tuple_y);

         void
         init_a (const Vector_Field_2D& vector_field_2d,
                 const Scalarization_2d scalarization_2d,
                 const Integer vector_element_0,
                 const Integer vector_element_1,
                 const Tuple& coordinate_tuple_x,
                 const Tuple& coordinate_tuple_y);

         void
         init_b (const Tuple& level_tuple,
                 const Real epsilon);

         void
         init (const Vector_Field_2D& vector_field_2d,
               const Integer vector_element,
               const Real step,
               const Tuple& coordinate_tuple_x,
               const Tuple& coordinate_tuple_y,
               const Real epsilon);

         void
         init (const Vector_Field_2D& vector_field_2d,
               const Integer vector_element,
               const Tuple& level_tuple,
               const Tuple& coordinate_tuple_x,
               const Tuple& coordinate_tuple_y,
               const Real epsilon);

         void
         init (const Vector_Field_2D& vector_field_2d,
               const Scalarization_2d scalarization_2d,
               const Integer vector_element_0,
               const Integer vector_element_1,
               const Real step,
               const Tuple& coordinate_tuple_x,
               const Tuple& coordinate_tuple_y,
               const Real epsilon);

         void
         init (const Vector_Field_2D& vector_field_2d,
               const Scalarization_2d scalarization_2d,
               const Integer vector_element_0,
               const Integer vector_element_1,
               const Tuple& level_tuple,
               const Tuple& coordinate_tuple_x,
               const Tuple& coordinate_tuple_y,
               const Real epsilon);

      public:

         // 1
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Tuple& level_tuple,
                  const Real epsilon = GSL_NAN);

         // 2
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Real step,
                  const Real epsilon = GSL_NAN);

         // 3
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Tuple& level_tuple,
                  const Domain_2D& domain_2d,
                  const Real epsilon = GSL_NAN);

         // 4
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Real step,
                  const Domain_2D& domain_2d,
                  const Real epsilon = GSL_NAN);

         // 5
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Tuple& level_tuple,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon = GSL_NAN);

         // 6
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Real step,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon = GSL_NAN);

         // 7
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Tuple& level_tuple,
                  const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Real epsilon = GSL_NAN);

         // 8
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Real step,
                  const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Real epsilon = GSL_NAN);

         // 9
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Tuple& level_tuple,
                  const Real epsilon = GSL_NAN);

         // 10
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Real step,
                  const Real epsilon = GSL_NAN);

         // 11
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Tuple& level_tuple,
                  const Domain_2D& domain_2d,
                  const Real epsilon = GSL_NAN);

         // 12
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Real step,
                  const Domain_2D& domain_2d,
                  const Real epsilon = GSL_NAN);

         // 13
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Tuple& level_tuple,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon = GSL_NAN);

         // 14
         Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Real step,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon = GSL_NAN);

         // 15
         Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Tuple& level_tuple,
                  const Real epsilon = GSL_NAN);

         // 16
         Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Real step,
                  const Real epsilon = GSL_NAN);

         // 17
         Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Tuple& level_tuple,
                  const Domain_2D& domain_2d,
                  const Real epsilon = GSL_NAN);

         // 18
         Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Real step,
                  const Domain_2D& domain_2d,
                  const Real epsilon = GSL_NAN);

         // 19
         Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Tuple& level_tuple,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon = GSL_NAN);

         // 20
         Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Real step,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon = GSL_NAN);

         // 21
         Contour (const Vector_Field_2D& scalar_field_2d,
                  const Integer vector_element,
                  const Tuple& level_tuple,
                  const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Real epsilon = GSL_NAN);

         // 22
         Contour (const Vector_Field_2D& scalar_field_2d,
                  const Integer vector_element,
                  const Real step,
                  const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Real epsilon = GSL_NAN);

         ~Contour ();

         const Polygon&
         get_polygon (const Integer index) const;

         void
         work_out_polygons ();

         const Real&
         get_max_value () const;

         const Real&
         get_min_value () const;

         void
         render (const RefPtr<Context>& cr,
                 const Transform_2D& transform,
                 const Color_Chooser& color_chooser);

         void
         render_fill (const RefPtr<Context>& cr,
                      const Transform_2D& transform,
                      const Color_Chooser& color_chooser);

         void
         render_isoline (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Integer level_index,
                         const Dstring& format = "",
                         const Real label_multiplier = GSL_NAN,
                         const Real label_offset = 0,
                         const Real label_distance = 40,
                         const Integer label_stride = -1) const;

         void
         render_isolines (const RefPtr<Context>& cr,
                          const Transform_2D& transform,
                          const Dstring& format = "",
                          const Real label_multiplier = GSL_NAN,
                          const Real label_offset = 0,
                          const Real label_distance = 40,
                          const Integer label_stride = -1) const;

   };

   class Scalar_Raster : public Raster
   {

      private:

         const Transform_2D&
         transform;

      public:

         Scalar_Raster (const Transform_2D& transform,
                        const Scalar_Data_2D& scalar_data_2d,
                        const Box_2D& box_2d,
                        const Color_Chooser& color_chooser);

   };

   class Scalar_Renderer
   {

      public:

         static void
         render (const RefPtr<Context>& cr,
                 const Transform_2D& transform,
                 const Vector_Data_2D& vector_data_2d,
                 const Integer vector_element,
                 const Domain_2D& domain,
                 const Color_Chooser& color_chooser);

         static void
         render (const RefPtr<Context>& cr,
                 const Transform_2D& transform,
                 const Vector_Data_2D& vector_data_2d,
                 const Integer vector_element,
                 const Color_Chooser& color_chooser);

         static void
         render (const RefPtr<Context>& cr,
                 const Transform_2D& transform,
                 const Scalar_Data_2D& scalar_data_2d,
                 const Domain_2D& domain,
                 const Color_Chooser& color_chooser);

         static void
         render (const RefPtr<Context>& cr,
                 const Transform_2D& transform,
                 const Scalar_Data_2D& scalar_data_2d,
                 const Color_Chooser& color_chooser);

   };

}

#endif /* DENISE_VISUALIZE_H */

