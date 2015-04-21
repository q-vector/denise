//
// constream.h
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

#ifndef DENISE_CONSTREAM_H
#define DENISE_CONSTREAM_H

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_poly.h>
#include <denise/analysis.h>
#include <denise/cairo.h>
#include <denise/geometry.h>
#include <denise/util.h>

using namespace std;

namespace denise
{

   class Neighbor
   {

      public:

         Index_2D
         index_2d;

         Point_2D
         point_2d;

         Neighbor (const Index_2D& index_2d,
                   const Point_2D& point_2d);

   };

   class Isoline : public Simple_Polyline
   {

      public:

         bool
         looped;

         Point_1D
         contour_level;

         Index_2D
         front_si;

         Index_2D
         back_si;

         Index_2D
         front_si_2;

         Index_2D
         back_si_2;

         Isoline (const Point_1D contour_level,
                  const Point_2D point,
                  const Index_2D staggered_index);

         void
         render (const Cairo& cairo,
                 const Transform_2D& transform_2d) const;

   };

   class Isoline_Path : public Path
   {

      public:

         Duple
         range;

         Isoline_Path (const Duple& range);

         void
         render (const Cairo& cairo,
                 const Transform_2D& transform_2d,
                 const Color_Chooser& color_chooser) const;

   };

   class Contourer : public Grid_2D
   {

      private:

         Index_1D
         tuple_index;

         const Tuple_Field_2D&
         tuple_field_2d;

         Size_2D
         box_size_2d;

         Data_2D<Point_2D>*
         staggered_point_data_ptr;

         Scalar
         epsilon;

         list<Isoline*>
         isoline_ptr_list;

         list<Isoline_Path*>
         isoline_path_ptr_list;

         bool
         in_bounds (const Index_2D& staggered_index);

         bool
         on_boundary (const Index_2D& staggered_index);

         void
         grow_isoline (Isoline& isoline,
                       const vector<Neighbor>& neighbor_vector,
                       bool front);

         void
         elect_neighbors (vector<Neighbor>& nv,
                          const Direction direction,
                          const Index_2D& staggered_index);

         void
         grow_isoline (Isoline& isoline,
                       bool front);

         void
         grow_isoline (Isoline& isoline);

         void
         clear_staggered_array ();

         void
         set_staggered_array (Point_1D contour_level);

         Boundary
         get_boundary (const Index_2D& staggered_index);

         Index_2D
         get_box_index (const Index_2D& staggered_index);

         void
         append_to_path (bool forward,
                         Path& path,
                         Isoline& isoline,
                         Boundary& boundary,
                         Index_2D& box_index);

         Direction
         get_direction (const Index_2D& box_index,
                        Boundary boundary,
                        const Duple& range,
                        Index_2D& grid_index);

         bool
         grid_value_within_range (const Index_2D& grid_index,
                                  const Duple& range);

         void
         advance_along_boundary (Index_2D& grid_index,
                                 Direction& direction,
                                 Boundary& boundary);

         Isoline*
         search_isoline (list<Isoline*>& isoline_ptr_list,
                         const Index_2D& box_index,
                         bool& forward);

         void
         append_isoline (Path& path,
                         Isoline& isoline,
                         Boundary& boundary,
                         Index_2D& box_index,
                         bool forward);

         void
         fill_isoline_ptr_list (list<Isoline*>& isoline_ptr_list,
                                Point_1D contour_level);

         Isoline_Path*
         get_isoline_path_ptr (list<Isoline*> isoline_ptr_list,
                               const Duple& range);

         void
         append_isoline_ptrs (list<Isoline*>& destination,
                              list<Isoline*>& source);

         void
         free_isoline_ptr_list ();

         void
         free_isoline_path_ptr_list ();

         void
         analyze (const Tuple& contour_levels_tuple,
                  bool fill = true);

      public:

         Contourer (const Tuple_Data_2D& tuple_data_2d,
                    const Index_1D tuple_index,
                    const Tuple& contour_levels_tuple,
                    const bool fill = true,
                    const Scalar epsilon = 1e-9);

         Contourer (const Tuple_Data_2D& tuple_data_2d,
                    const Index_1D tuple_idnex,
                    const Tuple& contour_levels_tuple,
                    const Size_2D& size_2d,
                    const bool fill = true,
                    const Scalar epsilon = 1e-9);

         Contourer (const Tuple_Field_2D& tuple_field_2d,
                    const Index_1D tuple_index,
                    const Tuple& contour_levels_tuple,
                    const Size_2D& size_2d,
                    const Domain_1D& domain_x,
                    const Domain_1D& domain_y,
                    const bool fill = true,
                    const Scalar epsilon = 1e-9);

         Contourer (const Tuple_Field_2D& tuple_field_2d,
                    const Index_1D tuple_index,
                    const Tuple& contour_levels_tuple,
                    const Grid_2D& grid_2d,
                    const bool fill = true,
                    const Scalar epsilon = 1e-9);

         ~Contourer ();

         void
         render (const Cairo& cairo,
                 const Transform_2D& transform_2d,
                 bool render_lines = true,
                 const Color_Chooser* color_chooser_ptr = NULL) const;

   };

   class Scalar_Raster : public Raster_Cairo
   {

      private:

         const Transform_2D&
         transform_2d;

      public:

         Scalar_Raster (const Transform_2D& transform_2d,
                        const Scalar_Field_2D& scalar_field_2d,
                        const Box_2D& box_2d,
                        const Color_Chooser& color_chooser);

         Scalar_Raster (const Transform_2D& transform_2d,
                        const Tuple_Field_2D& tuple_field_2d,
                        const Index_1D tuple_index,
                        const Box_2D& box_2d,
                        const Color_Chooser& color_chooser);

   };

   class Scalar_Renderer
   {

      public:

         static void
         render (const Cairo& cairo,
                 const Transform_2D& transform_2d,
                 const Tuple_Data_2D& tuple_data_2d,
                 const Index_1D tuple_index,
                 const Color_Chooser& color_chooser);

         static void
         render (const Cairo& cairo,
                 const Transform_2D& transform_2d,
                 const Tuple_Data_2D& tuple_data_2d,
                 const Index_1D tuple_index,
                 const Color_Chooser& color_chooser,
                 const Size_2D& size_2d);

         static void
         render (const Cairo& cairo,
                 const Transform_2D& transform_2d,
                 const Tuple_Field_2D& tuple_field_2d,
                 const Index_1D tuple_index,
                 const Color_Chooser& color_chooser,
                 const Size_2D& size_2d,
                 const Domain_1D& domain_x,
                 const Domain_1D& domain_y);

         static void
         render (const Cairo& cairo,
                 const Transform_2D& transform_2d,
                 const Tuple_Field_2D& tuple_field_2d,
                 const Index_1D tuple_index,
                 const Color_Chooser& color_chooser,
                 const Grid_2D& grid_2d);


   };

/*
   typedef list<Point_2D>
   Streamline_Cell;

   class Streamline : public Polyline
   {

      public:

         bool
         is_separatrix;

         bool
         forward_ok;

         bool
         backward_ok;

         Streamline (const Point_2D& point);

         Streamline (const Point_2D& point_a,
                     const Point_2D& point_b,
                     bool forward);

         bool
         open () const;

   };

   class Streamline_Set : public vector<Streamline*>
   {

      private:

         Size_2D
         size_2d;

         Scalar
         delta_x;

         Scalar
         delta_y;

         Box
         bounding_box;

         Scalar
         separation;

         Scalar
         step_size;

         Scalar
         seed_separation;

         Size_1D
         max_count;

         Size_1D
         min_separatrix_count;
      
//         Streamline_Cell**
//         cell_array;
      
         Data_2D<Streamline_Cell>*
         cell_data_ptr;
      
         Vector_Data_2D*
         vector_data_ptr;
      
         set<Point_2D>*
         critical_point_set_ptr;
      
         vector<Point_2D>
         saddle_vector;
      
         vector<Point_2D>
         vortex_vector;
      
         vector<Point_2D>
         node_vector;

         void
         obtain_special_streamlines ();
      
         void
         seed_saddle (const Point_2D& saddle,
                      const Jacobian& jacobian,
                      Scalar B,
                      Scalar C);

         void
         obtain_ordinary_streamlines ();

         Index_2D
         get_index (const Point_2D& point) const;

         void
         add_to_cell (const Point_2D& point) const;

         void
         add_to_cell (Streamline& streamline) const;

         bool
         has_obstacle (const Point_2D& point,
                       const Vector& vector,
                       Scalar threshold) const;

         void
         grow (Streamline& streamline,
               bool backward,
               bool watch_ahead) const;

         void
         grow (Streamline& streamline,
               bool watch_ahead) const;

         Streamline*
         get_streamline_ptr (const Point_2D& point) const;

      public:

         Streamline_Set (Vector_Data_2D& vector_data,
                         const Box& bounding_box,
                         Scalar separation,
                         Size_1D step_size_granulity = 10,
                         Scalar seed_separation_ratio = 1.2,
                         Scalar max_streamline_length = 1000,
                         Scalar min_separatrix_length = 0.1);

         ~Streamline_Set ();

         const vector<Point_2D>&
         get_saddle_vector () const;

         const vector<Point_2D>&
         get_vortex_vector () const;


   };
*/

}

#endif /* DENISE_CONSTREAM_H */

