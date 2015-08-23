//
// constream.cc
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

#include "constream.h"

using namespace std;
using namespace denise;

Neighbor::Neighbor (const Index_2D& index_2d,
                    const Point_2D& point_2d)
        : index_2d (index_2d),
          point_2d (point_2d)
{
}

Isoline::Isoline (Point_1D contour_level,
                  Point_2D point,
                  Index_2D staggered_index)
{

   looped = true;
   add (point);

   this->contour_level = contour_level;
   this->front_si = staggered_index;
   this->back_si = staggered_index;
   this->front_si_2 = Index_2D (-1, -1);
   this->back_si_2 = Index_2D (-1, -1);

}

void
Isoline::render (const Cairo& cairo,
                 const Transform_2D& transform_2d) const
{
   cairo.simple_polyline (transform_2d, *this, looped);
   cairo.stroke ();
}

Isoline_Path::Isoline_Path (const Duple& range)
{
   this->range = range;
}

void
Isoline_Path::render (const Cairo& cairo,
                      const Transform_2D& transform_2d,
                      const Color_Chooser& color_chooser) const
{

   Point_1D level = (range.first + range.second) / 2;
   Color color = color_chooser.get_color (level);

   if (!color.is_nac ())
   {
      cairo.set_color (color);
      cairo.path (transform_2d, *this);
      cairo.fill ();
   }

}

bool
Contourer::in_bounds (const Index_2D& staggered_index)
{

   const Index_2D& si = staggered_index;
   const Size_2D& ss = staggered_point_data_ptr->get_size_2d ();

   return (si.i >= 0 && si.j >= 0 && si.j < ss.j &&
           (si.i < ss.i - 1 || (si.i == ss.i - 1 && si.j%2 == 1)));

}

bool
Contourer::on_boundary (const Index_2D& staggered_index)
{

   const Index_2D& si = staggered_index;
   const Size_2D& ss = staggered_point_data_ptr->get_size_2d ();

   return ((si.i == 0 && si.j%2 == 1) || (si.j == 0) ||
           (si.i == ss.i - 1) || (si.j == ss.j - 1));

}

void
Contourer::grow_isoline (Isoline& isoline,
                         const vector<Neighbor>& neighbor_vector,
                         bool front)
{

   const Neighbor& neighbor = neighbor_vector.front ();
   const Index_2D& si = neighbor.index_2d;
   const Point_2D& point = neighbor.point_2d;

   if (front)
   {
      isoline.front_si_2 = isoline.front_si;
      isoline.front_si = si;
      isoline.prepend (point);
   }
   else
   {
      isoline.back_si_2 = isoline.back_si;
      isoline.back_si = si;
      isoline.add (point);
   }

   grow_isoline (isoline, front);

}

void
Contourer::elect_neighbors (vector<Neighbor>& nv,
                            const Direction direction,
                            const Index_2D& staggered_index)
{

   Index_2D ni[3] = { staggered_index, staggered_index, staggered_index };

   switch (direction)
   {

      case UPWARD:
         ni[0].j += 1;
         ni[1].i += 1; ni[1].j += 1;
         ni[2].j += 2;
         break;

      case DOWNWARD:
         ni[0].i += 1; ni[0].j -= 1;
         ni[1].j -= 1;
         ni[2].j -= 2;
         break;

      case LEFTWARD:
         ni[0].i -= 1; ni[0].j -= 1;
         ni[1].i -= 1; ni[1].j += 1;
         ni[2].i -= 1;
         break;

      case RIGHTWARD:
         ni[0].j += 1;
         ni[1].j -= 1;
         ni[2].i += 1;
         break;

   }

   for (Index_1D i = 0; i < 3; i++)
   {
      const Index_2D& neighbor_index = ni[i];
      if (in_bounds (neighbor_index))
      {
         const Point_2D& p = staggered_point_data_ptr->get (neighbor_index);
         if (!p.is_nap ()) { nv.push_back (Neighbor (neighbor_index, p)); }
      }
   }

}

void
Contourer::grow_isoline (Isoline& isoline,
                         bool front)
{

   Index_2D& si = (front ? isoline.front_si : isoline.back_si);
   Index_2D& si_2 = (front ? isoline.front_si_2 : isoline.back_si_2);

   vector<Neighbor> nv;
   staggered_point_data_ptr->set (Scalar (GSL_NAN), si.i, si.j);

   if (si_2.i < 0 || si_2.j < 0)
   {

      if (si.j%2 == 0)
      {
         elect_neighbors (nv, UPWARD, si);
         if (nv.size () == 0) { elect_neighbors (nv, DOWNWARD, si); }
      }
      else
      {
         elect_neighbors (nv, LEFTWARD, si);
         if (nv.size () == 0) { elect_neighbors (nv, RIGHTWARD, si); }
      }

   }
   else
   {

      if (si.j%2 == 0)
      {
         if (si_2.j < si.j) { elect_neighbors (nv, UPWARD, si); }
         else { elect_neighbors (nv, DOWNWARD, si); }
      }
   
      else
      {
         if (si_2.i < si.i) { elect_neighbors (nv, RIGHTWARD, si); }
         else { elect_neighbors (nv, LEFTWARD, si); }
      }

   }
   
   if (nv.size () > 0 ) { grow_isoline (isoline, nv, front); }

}

void
Contourer::grow_isoline (Isoline& isoline)
{
   grow_isoline (isoline, true);
   grow_isoline (isoline, false);
}

void
Contourer::clear_staggered_array ()
{

   const Size_2D& ss = staggered_point_data_ptr->get_size_2d ();

   for (Index_1D i = 0; i < ss.i; i++)
   {
      for (Index_1D j = 0; j < ss.j; j++)
      {
         Point_2D& point = staggered_point_data_ptr->get (i, j);
         point.x = GSL_NAN;
         point.y = GSL_NAN;
      }
   }

}

void
Contourer::set_staggered_array (Point_1D contour_level)
{

   clear_staggered_array ();

   for (Index_1D i = 0; i < size_2d.i - 1; i++)
   {

      Point_1D a = get_x (i);
      Point_1D b = get_x (i+1);
      Scalar delta_x = b - a;

      for (Index_1D j = 0; j < size_2d.j - 1; j++)
      {

         Point_1D c = get_y (j);
         Point_1D d = get_y (j+1);
         Scalar delta_y = d - c;

         Scalar m = tuple_field_2d.evaluate_at (a, c, tuple_index);
         if (gsl_isnan (m)) { continue; }
         if (m == contour_level) { m += epsilon; }

         Scalar n = tuple_field_2d.evaluate_at (b, c, tuple_index);
         if (gsl_isnan (n)) { continue; }
         if (n == contour_level) { n += epsilon; }

         Scalar o = tuple_field_2d.evaluate_at (a, d, tuple_index);
         if (gsl_isnan (o)) { continue; }
         if (o == contour_level) { o += epsilon; }

         Scalar p = tuple_field_2d.evaluate_at (b, d, tuple_index);
         if (gsl_isnan (p)) { continue; }
         if (p == contour_level) { p += epsilon; }

         if ((contour_level - m) * (contour_level - o) < 0)
         {
            Point_2D& point = staggered_point_data_ptr->get (i, j*2+1);
            point.x = a;
            point.y = (contour_level - m) / (o - m) * delta_y + c;
         }

         if ((contour_level - n) * (contour_level - p) < 0)
         {
            Point_2D& point = staggered_point_data_ptr->get (i+1, j*2+1);
            point.x = b;
            point.y = (contour_level - n) / (p - n) * delta_y + c;
         }

         if ((contour_level - m) * (contour_level - n) < 0)
         {
            Point_2D& point = staggered_point_data_ptr->get (i, j*2);
            point.x = (contour_level - m) / (n - m) * delta_x + a;
            point.y = c;
         }

         if ((contour_level - o) * (contour_level - p) < 0)
         {
            Point_2D& point = staggered_point_data_ptr->get (i, j*2+2);
            point.x = (contour_level - o) / (p - o) * delta_x + a;
            point.y = d;
         }

      }
   }

}

Boundary
Contourer::get_boundary (const Index_2D& staggered_index)
{

   Boundary boundary;
   const Index_2D& si = staggered_index;
   const Index_2D& ss = staggered_point_data_ptr->get_size_2d ();

   // do now swap: j is to given preference over i
        if (si.j == 0)        { boundary = BOUNDARY_BOTTOM; }
   else if (si.j == ss.j - 1) { boundary = BOUNDARY_TOP;    }
   else if (si.i == ss.i - 1) { boundary = BOUNDARY_RIGHT;  }
   else if (si.i == 0)        { boundary = BOUNDARY_LEFT;   }

   return boundary;

}

Index_2D
Contourer::get_box_index (const Index_2D& staggered_index)
{

   Index_2D box_index;
   Boundary boundary = get_boundary (staggered_index);

   switch (boundary)
   {

      case BOUNDARY_LEFT: 
         box_index.i = 0;
         box_index.j = (staggered_index.j - 1)/2;
         break;

      case BOUNDARY_RIGHT: 
         box_index.i = box_size_2d.i - 1;
         box_index.j = (staggered_index.j - 1)/2;
         break;

      case BOUNDARY_TOP: 
         box_index.i = staggered_index.i;
         box_index.j = box_size_2d.j - 1;
         break;

      case BOUNDARY_BOTTOM: 
         box_index.i = staggered_index.i;
         box_index.j = 0;
         break;

   }

   return box_index;

}

void
Contourer::append_to_path (bool forward,
                           Path& path,
                           Isoline& isoline,
                           Boundary& boundary,
                           Index_2D& box_index)
{

   if (forward)
   {

      for (Isoline::iterator iterator = isoline.begin ();
           iterator != isoline.end (); iterator++)
      {
         path.line_to (*(iterator));
      }

      boundary = get_boundary (isoline.back_si);
      box_index = get_box_index (isoline.back_si);

   }

   else
   {

      for (Isoline::reverse_iterator iterator = isoline.rbegin ();
           iterator != isoline.rend (); iterator++)
      {
         path.line_to (*(iterator));
      }

      boundary = get_boundary (isoline.front_si);
      box_index = get_box_index (isoline.front_si);

   }

}

Direction
Contourer::get_direction (const Index_2D& box_index,
                          Boundary boundary,
                          const Duple& range,
                          Index_2D& grid_index)
{

   Direction direction_a, direction_b;
   Index_2D grid_index_a, grid_index_b;

   switch (boundary)
   {

      case BOUNDARY_LEFT:
         direction_a = UPWARD;
         direction_b = DOWNWARD;
         grid_index_a.i = 0;
         grid_index_a.j = box_index.j + 1;
         grid_index_b.i = 0;
         grid_index_b.j = box_index.j;
         break;

      case BOUNDARY_RIGHT:
         direction_a = UPWARD;
         direction_b = DOWNWARD;
         grid_index_a.i = box_index.i + 1;
         grid_index_a.j = box_index.j + 1;
         grid_index_b.i = box_index.i + 1;
         grid_index_b.j = box_index.j;
         break;

      case BOUNDARY_TOP:
         direction_a = LEFTWARD;
         direction_b = RIGHTWARD;
         grid_index_a.i = box_index.i;
         grid_index_a.j = box_index.j + 1;
         grid_index_b.i = box_index.i + 1;
         grid_index_b.j = box_index.j + 1;
         break;

      case BOUNDARY_BOTTOM:
         direction_a = LEFTWARD;
         direction_b = RIGHTWARD;
         grid_index_a.i = box_index.i;
         grid_index_a.j = 0;
         grid_index_b.i = box_index.i + 1;
         grid_index_b.j = 0;
         break;

   }

   Direction direction;

   if (grid_value_within_range (grid_index_a, range))
   {
      direction = direction_a;
      grid_index = grid_index_a;
   }

   else// if (grid_value_within_range (grid_index_b, range))
   {
      direction = direction_b;
      grid_index = grid_index_b;
   }

   return direction;

}

bool
Contourer::grid_value_within_range (const Index_2D& grid_index,
                                    const Duple& range)
{

   Point_1D x = get_x (grid_index.i);
   Point_1D y = get_y (grid_index.j);
   Scalar value = tuple_field_2d.evaluate_at (x, y, tuple_index);

   if (value == range.first || value == range.second) { value += epsilon; }
   return ((value - range.first) * (value - range.second) < 0);

}

void
Contourer::advance_along_boundary (Index_2D& grid_index, 
                                   Direction& direction,
                                   Boundary& boundary)
{

   if (grid_index.i == 0 && grid_index.j == 0)
   {

      if (direction == LEFTWARD)
      {
         direction = UPWARD;
         boundary = BOUNDARY_LEFT;
      }

      else
      if (direction == DOWNWARD)
      {
         direction = RIGHTWARD;
         boundary = BOUNDARY_BOTTOM;
      }

   }

   else
   if (grid_index.i == 0 && grid_index.j == size_2d.j - 1)
   {

      if (direction == UPWARD)
      {
         direction = RIGHTWARD;
         boundary = BOUNDARY_TOP;
      }

      else
      if (direction == LEFTWARD)
      {
         direction = DOWNWARD;
         boundary = BOUNDARY_LEFT;
      }

   }

   else
   if (grid_index.i == size_2d.i - 1 && grid_index.j == 0)
   {

      if (direction == RIGHTWARD)
      {
         direction = UPWARD;
         boundary = BOUNDARY_RIGHT;
      }

      else
      if (direction == DOWNWARD)
      {
         direction = LEFTWARD;
         boundary = BOUNDARY_BOTTOM;
      }

   }

   else
   if (grid_index.i == size_2d.i - 1 && grid_index.j == size_2d.j - 1)
   {

      if (direction == RIGHTWARD)
      {
         direction = DOWNWARD;
         boundary = BOUNDARY_RIGHT;
      }

      else
      if (direction == UPWARD)
      {
         direction = LEFTWARD;
         boundary = BOUNDARY_TOP;
      }

   }

   switch (direction)
   {
      case LEFTWARD:  grid_index.i--; break;
      case RIGHTWARD: grid_index.i++; break;
      case UPWARD:    grid_index.j++; break;
      case DOWNWARD:  grid_index.j--; break;
   }

}

Isoline*
Contourer::search_isoline (list<Isoline*>& isoline_ptr_list,
                           const Index_2D& box_index,
                           bool& forward)
{

   Isoline* isoline_ptr = NULL;

   for (list<Isoline*>::iterator iterator = isoline_ptr_list.begin ();
        iterator != isoline_ptr_list.end (); iterator++)
   {
      Isoline& isoline = **(iterator);
   }

   for (list<Isoline*>::iterator iterator = isoline_ptr_list.begin ();
        iterator != isoline_ptr_list.end (); iterator++)
   {

      Isoline& isoline = **(iterator);

      if (box_index == get_box_index (isoline.front_si))
      {
         forward = true;
         isoline_ptr = &isoline;
         iterator = isoline_ptr_list.erase (iterator);
         break;
      }

      if (box_index == get_box_index (isoline.back_si))
      {
         forward = false;
         isoline_ptr = &isoline;
         iterator = isoline_ptr_list.erase (iterator);
         break;
      }

   }

   return isoline_ptr;

}

void
Contourer::append_isoline (Path& path,
                           Isoline& isoline,
                           Boundary& boundary,
                           Index_2D& box_index,
                           bool forward)
{

   Index_2D& si = (forward ? isoline.back_si :
                             isoline.front_si);

   append_to_path (forward, path, isoline, boundary, box_index);
   box_index = get_box_index (si);

}

void
Contourer::fill_isoline_ptr_list (list<Isoline*>& isoline_ptr_list,
                                  Point_1D contour_level)
{

   Point_2D p;
   Index_2D si;
   set_staggered_array (contour_level);
   const Size_2D& ss = staggered_point_data_ptr->get_size_2d ();

   for (si.i = 0; si.i < ss.i; si.i++)
   {
      for (si.j = 0; si.j < ss.j; si.j++)
      {

         Point_2D& point = staggered_point_data_ptr->get (si.i, si.j);

         if (!point.is_nap ())
         {
            Isoline* isoline_ptr = new Isoline (contour_level, point, si);
            grow_isoline (*isoline_ptr);
            isoline_ptr_list.push_back (isoline_ptr);
            point.x = GSL_NAN; point.y = GSL_NAN;
         }

      }
   }

   for (list<Isoline*>::iterator iterator = isoline_ptr_list.begin ();
        iterator != isoline_ptr_list.end (); iterator++)
   {

      Isoline& isoline = **(iterator);

      if (on_boundary (isoline.front_si) && on_boundary (isoline.back_si))
      {
         isoline.looped = false;
      }


   }

}

Isoline_Path*
Contourer::get_isoline_path_ptr (list<Isoline*> isoline_ptr_list,
                                 const Duple& range)
{

   Isoline_Path* path_ptr = new Isoline_Path (range);
   Isoline_Path& path = *path_ptr;

   // Inner loops
   for (list<Isoline*>::iterator i = isoline_ptr_list.begin ();
        i != isoline_ptr_list.end (); i++)
   {

      Isoline* isoline_ptr = *(i);
      Isoline& isoline = *isoline_ptr;

      if (isoline_ptr->looped)
      {

         for (Isoline::iterator j = isoline.begin ();
              j != isoline.end (); j++)
         {
            if (j == isoline.begin ()) { path.move_to (*(j)); }
            else                       { path.line_to (*(j)); }
         }

         path.line_to ((*isoline_ptr).front ());
         path.close ();

         i = isoline_ptr_list.erase (i);
         i--;
      }

   }

   while (isoline_ptr_list.size () > 0)
   {

      bool forward = true;
      Direction direction;
      Index_2D grid_index;
      Point_2D latest_point;

      Isoline& first_isoline = *(isoline_ptr_list.front ());
      Index_2D box_index = get_box_index (first_isoline.front_si);
      Index_2D first_box_index = box_index;

      Isoline* isoline_ptr = isoline_ptr_list.front ();
      isoline_ptr_list.pop_front ();

      Boundary boundary = get_boundary (first_isoline.front_si);
      path.move_to (first_isoline.front ());

      while (isoline_ptr != NULL)
      {

         append_isoline (path, *isoline_ptr, boundary, box_index, forward);
         direction = get_direction (box_index, boundary, range, grid_index);

         while (grid_value_within_range (grid_index, range))
         {

            latest_point = get_point_2d (grid_index);
            path.line_to (latest_point);
            advance_along_boundary (grid_index, direction, boundary);

            box_index = grid_index;
            if (boundary == BOUNDARY_RIGHT) { box_index.i--; }
            if (boundary == BOUNDARY_TOP) { box_index.j--; }
            if (direction == RIGHTWARD) { box_index.i--; }
            if (direction == UPWARD) { box_index.j--; }

         }

         isoline_ptr = search_isoline (isoline_ptr_list, box_index, forward);
         if (isoline_ptr == NULL)
         {
            break;
         }

         if (get_box_index (isoline_ptr->front_si) ==
             get_box_index (isoline_ptr->back_si))
         {

            Point_2D& f = isoline_ptr->front ();
            Point_2D& b = isoline_ptr->back ();

            Scalar dx_f = f.x - latest_point.x;
            Scalar dy_f = f.y - latest_point.y;
            Scalar dx_b = b.x - latest_point.x;
            Scalar dy_b = b.y - latest_point.y;

            Scalar distance_f = sqrt ((dx_f*dx_f) + (dy_f*dy_f));
            Scalar distance_b = sqrt ((dx_b*dx_b) + (dy_b*dy_b));

            forward = (distance_f < distance_b);

         }

      }

      path.line_to (first_isoline.front ());
      path.close ();

   }

   return path_ptr;

}

void
Contourer::append_isoline_ptrs (list<Isoline*>& destination,
                                list<Isoline*>& source)
{

   for (list<Isoline*>::iterator iterator = source.begin ();
        iterator != source.end (); iterator++)
   {
      destination.push_back (*iterator);
   }

}

void
Contourer::free_isoline_ptr_list ()
{

   for (list<Isoline*>::iterator iterator = isoline_ptr_list.begin ();
        iterator != isoline_ptr_list.end (); iterator++)
   {
      Isoline* isoline_ptr = *(iterator);
      delete isoline_ptr;
   }

   isoline_ptr_list.clear ();

}

void
Contourer::free_isoline_path_ptr_list ()
{

   for (list<Isoline_Path*>::iterator iterator = isoline_path_ptr_list.begin ();
        iterator != isoline_path_ptr_list.end (); iterator++)
   {
      Isoline_Path* isoline_path_ptr = *(iterator);
      delete isoline_path_ptr;
   }

   isoline_path_ptr_list.clear ();

}

void
Contourer::analyze (const Tuple& contour_levels_tuple,
                    bool fill)
{

   free_isoline_ptr_list ();
   free_isoline_path_ptr_list ();

   if (!fill)
   {

      for (Tuple::const_iterator iterator = contour_levels_tuple.begin ();
           iterator != contour_levels_tuple.end (); iterator++)
      {
         const Point_1D& contour_level = *(iterator);
         fill_isoline_ptr_list (isoline_ptr_list, contour_level);
      }

   }

   else
   {

      Size_1D n = contour_levels_tuple.size ();
      list<Isoline*>* isoline_ptr_lists = new list<Isoline*>[n];

      for (Index_1D i = 0; i < n; i++)
      {

         const Point_1D contour_level = contour_levels_tuple[i];
         fill_isoline_ptr_list (isoline_ptr_lists[i], contour_level);

         if (i != 0)
         {

            list<Isoline*> ipl;
            append_isoline_ptrs (ipl, isoline_ptr_lists[i-1]);
            append_isoline_ptrs (ipl, isoline_ptr_lists[i]);

            Duple range (contour_levels_tuple[i-1], contour_levels_tuple[i]);
            Isoline_Path* isoline_path_ptr = get_isoline_path_ptr (ipl, range);
            isoline_path_ptr_list.push_back (isoline_path_ptr);

         }

      }

      for (Index_1D i = 0; i < n; i++)
      {
         append_isoline_ptrs (isoline_ptr_list, isoline_ptr_lists[i]);
      }

      delete[] isoline_ptr_lists;

   }

}

Contourer::Contourer (const Tuple_Data_2D& tuple_data_2d,
                      const Index_1D tuple_index,
                      const Tuple& contour_levels_tuple,
                      const bool fill,
                      const Scalar epsilon)
           : Grid_1D ((const Grid_1D&)tuple_data_2d),
             Grid_2D ((const Grid_2D&)tuple_data_2d),
      tuple_field_2d ((const Tuple_Field_2D&)tuple_data_2d),
         tuple_index (tuple_index),
             epsilon (epsilon),
         box_size_2d (size_2d.i - 1, size_2d.j - 1)
{
   Size_2D staggered_size_2d = Size_2D (size_2d.i, size_2d.j * 2 - 1);
   staggered_point_data_ptr = new Data_2D<Point_2D>(staggered_size_2d);
   analyze (contour_levels_tuple, fill);
}

Contourer::Contourer (const Tuple_Data_2D& tuple_data_2d,
                      const Index_1D tuple_index,
                      const Tuple& contour_levels_tuple,
                      const Size_2D& size_2d,
                      const bool fill,
                      const Scalar epsilon)
           : Grid_1D (size_2d.i, tuple_data_2d.get_domain_x ()),
             Grid_2D (size_2d, tuple_data_2d.get_domain_x (), tuple_data_2d.get_domain_y ()),
      tuple_field_2d ((const Tuple_Field_2D&)tuple_data_2d),
         tuple_index (tuple_index),
             epsilon (epsilon),
         box_size_2d (size_2d.i - 1, size_2d.j - 1)
{
   Size_2D staggered_size_2d = Size_2D (size_2d.i, size_2d.j * 2 - 1);
   staggered_point_data_ptr = new Data_2D<Point_2D>(staggered_size_2d);
   analyze (contour_levels_tuple, fill);
}

Contourer::Contourer (const Tuple_Field_2D& tuple_field_2d,
                      const Index_1D tuple_index,
                      const Tuple& contour_levels_tuple,
                      const Size_2D& size_2d,
                      const Domain_1D& domain_x,
                      const Domain_1D& domain_y,
                      const bool fill,
                      const Scalar epsilon)
           : Grid_1D (size_2d.i, domain_x),
             Grid_2D (size_2d, domain_x, domain_y),
      tuple_field_2d (tuple_field_2d),
         tuple_index (tuple_index),
             epsilon (epsilon),
         box_size_2d (size_2d.i - 1, size_2d.j - 1)
{
   Size_2D staggered_size_2d = Size_2D (size_2d.i, size_2d.j * 2 - 1);
   staggered_point_data_ptr = new Data_2D<Point_2D>(staggered_size_2d);
   analyze (contour_levels_tuple, fill);
}

Contourer::Contourer (const Tuple_Field_2D& tuple_field_2d,
                      const Index_1D tuple_index,
                      const Tuple& contour_levels_tuple,
                      const Grid_2D& grid_2d,
                      const bool fill,
                      const Scalar epsilon)
           : Grid_1D (grid_2d.get_x_tuple ()),
             Grid_2D (grid_2d),
      tuple_field_2d (tuple_field_2d),
         tuple_index (tuple_index),
             epsilon (epsilon),
         box_size_2d (size_2d.i - 1, size_2d.j - 1)
{
   Size_2D staggered_size_2d = Size_2D (size_2d.i, size_2d.j * 2 - 1);
   staggered_point_data_ptr = new Data_2D<Point_2D>(staggered_size_2d);
   analyze (contour_levels_tuple, fill);
}

Contourer::~Contourer ()
{
   free_isoline_ptr_list ();
   free_isoline_path_ptr_list ();
   delete staggered_point_data_ptr;
}

void
Contourer::render (const Cairo& cairo,
                   const Transform_2D& transform_2d,
                   const bool render_line,
                   const Color_Chooser* color_chooser_ptr) const
{

   if (!isoline_path_ptr_list.empty () && color_chooser_ptr != NULL)
   {

      for (list<Isoline_Path*>::const_iterator iterator =
              isoline_path_ptr_list.begin ();
           iterator != isoline_path_ptr_list.end (); iterator++)
      {
         const Isoline_Path& isoline_path = **(iterator);
         isoline_path.render (cairo, transform_2d, *color_chooser_ptr);
      }

   }

   if (render_line)
   {
      for (list<Isoline*>::const_iterator iterator = isoline_ptr_list.begin ();
           iterator != isoline_ptr_list.end (); iterator++)
      {
         const Isoline& isoline = **(iterator);
         isoline.render (cairo, transform_2d);
      }
   }

}

Scalar_Raster::Scalar_Raster (const Transform_2D& transform_2d,
                              const Scalar_Field_2D& scalar_field_2d,
                              const Box_2D& box_2d,
                              const Color_Chooser& color_chooser)
              : Raster_Cairo (box_2d),
                transform_2d (transform_2d)
{

   Point_2D p;
   const Size_2D& size_2d = box_2d.size_2d;
   const Index_2D& anchor = box_2d.index_2d;

   Index_2D end_index (anchor.i + size_2d.i, anchor.j  + size_2d.j);

   for (Index_1D i = anchor.i; i < end_index.i; i++)
   {

      for (Index_1D j = anchor.j; j < end_index.j; j++)
      {

         Point_2D p = transform_2d.transform (Scalar (i), Scalar (j), true);

         try
         {
            Scalar value = scalar_field_2d.evaluate_at (p.x, p.y);
            Color color = color_chooser.get_color (value);
            set_pixel (i - anchor.i, j - anchor.j, color);
         }
         catch (const Out_Of_Bounds_Exception& oobe) { cerr << oobe << endl; }
         catch (const Exception& e) { cerr << e << endl; }

      }
   }

}

Scalar_Raster::Scalar_Raster (const Transform_2D& transform_2d,
                              const Tuple_Field_2D& tuple_field_2d,
                              const Index_1D tuple_index,
                              const Box_2D& box_2d,
                              const Color_Chooser& color_chooser)
              : Raster_Cairo (box_2d),
                transform_2d (transform_2d)
{

   Point_2D p;
   const Size_2D& size_2d = box_2d.size_2d;
   const Index_2D& anchor = box_2d.index_2d;

   Index_2D end_index (anchor.i + size_2d.i, anchor.j  + size_2d.j);

   for (Index_1D i = anchor.i; i < end_index.i; i++)
   {

      for (Index_1D j = anchor.j; j < end_index.j; j++)
      {

         Point_2D p = transform_2d.transform (Scalar (i), Scalar (j), true);

         try
         {
            Scalar value = tuple_field_2d.evaluate_at (p.x, p.y, tuple_index);
            Color color = color_chooser.get_color (value);
            set_pixel (i - anchor.i, j - anchor.j, color);
         }
         catch (const Out_Of_Bounds_Exception& oobe) { }
         catch (const Exception& e) { }

      }
   }

}

void
Scalar_Renderer::render (const Cairo& cairo,
                         const Transform_2D& transform_2d,
                         const Tuple_Data_2D& tuple_data_2d,
                         const Index_1D tuple_index,
                         const Color_Chooser& color_chooser)
{
   render (cairo, transform_2d, tuple_data_2d,
      tuple_index, color_chooser, tuple_data_2d);
}

void
Scalar_Renderer::render (const Cairo& cairo,
                         const Transform_2D& transform_2d,
                         const Tuple_Data_2D& tuple_data_2d,
                         const Index_1D tuple_index,
                         const Color_Chooser& color_chooser,
                         const Size_2D& size_2d)
{

   Grid_2D grid_2d (size_2d, tuple_data_2d.get_domain_x (),
      tuple_data_2d.get_domain_y ());

   render (cairo, transform_2d, tuple_data_2d,
      tuple_index, color_chooser, grid_2d);

}

void
Scalar_Renderer::render (const Cairo& cairo,
                         const Transform_2D& transform_2d,
                         const Tuple_Field_2D& tuple_field_2d,
                         const Index_1D tuple_index,
                         const Color_Chooser& color_chooser,
                         const Size_2D& size_2d,
                         const Domain_1D& domain_x,
                         const Domain_1D& domain_y)
{
   Grid_2D grid_2d (size_2d, domain_x, domain_y);
   render (cairo, transform_2d, tuple_field_2d,
      tuple_index, color_chooser, grid_2d);
}

void
Scalar_Renderer::render (const Cairo& cairo,
                         const Transform_2D& transform_2d,
                         const Tuple_Field_2D& tuple_field_2d,
                         const Index_1D tuple_index,
                         const Color_Chooser& color_chooser,
                         const Grid_2D& grid_2d)
{

   const Size_2D& size_2d = grid_2d.get_size_2d ();

   for (Index_1D i = 0; i < size_2d.i; i++)
   {

      Point_1D x = grid_2d.get_x (i);
      Point_1D xb = x;
      Point_1D xe = x;

      if (i != 0) { xb += grid_2d.get_x (i-1); xb /= 2; }
      if (i != size_2d.i-1) { xe += grid_2d.get_x (i+1); xe /= 2; }

      for (Index_1D j = 0; j < size_2d.j; j++)
      {

         Point_1D y = grid_2d.get_y (j);
         Point_1D yb = y;
         Point_1D ye = y;

         if (j != 0) { yb += grid_2d.get_y (j-1); yb /= 2; }
         if (j != size_2d.j-1) { ye += grid_2d.get_y (j+1); ye /= 2; }

         Point_2D point = grid_2d.get_point_2d (Index_2D (i, j));

         Scalar value = tuple_field_2d.evaluate_at (x, y, tuple_index);
         if (gsl_isnan (value)) { continue; }

         Color color = color_chooser.get_color (value);
         if (color.is_nac ()) { continue; }

         cairo.set_color (color);

         Point_2D bb = transform_2d.transform (xb, yb);
         Point_2D be = transform_2d.transform (xb, ye);
         Point_2D eb = transform_2d.transform (xe, yb);
         Point_2D ee = transform_2d.transform (xe, ye);

         Polygon polygon;
         polygon.add (bb);
         polygon.add (be);
         polygon.add (ee);
         polygon.add (eb);

         cairo.polygon (polygon);
         cairo.fill ();

      }

   }

}

/*
Streamline::Streamline (const Point_2D& point)
{

   push_back (point);
   is_separatrix = false;

   forward_ok = true;
   backward_ok = true;

}

Streamline::Streamline (const Point_2D& point_a,
                        const Point_2D& point_b,
                        bool forward)
{
   push_back (point_a);
   is_separatrix = false;

   forward_ok = forward;
   backward_ok = !forward_ok;

   if (forward) { push_back (point_b); }
   else         { push_front (point_b); }
}

bool
Streamline::open () const
{
   return (forward_ok || backward_ok);
}

void
Streamline_Set::obtain_special_streamlines ()
{

   Size_1D norr;
   Vector nav (GSL_NAN, GSL_NAN);
   Scalar lambda_0, lambda_1;

   for (std::set<Point_2D>::const_iterator iterator =
           critical_point_set_ptr->begin ();
        iterator != critical_point_set_ptr->end (); iterator++)
   {

      const Point_2D point = *(iterator);
      if (bounding_box.is_out_of_bounds (point)) { continue; }

      Jacobian jacobian = vector_data_ptr->get_jacobian_at (point);
      Scalar b = -(jacobian.du_dx + jacobian.dv_dy);
      Scalar c = jacobian.get_determinant ();

      norr = gsl_poly_solve_quadratic (1, b, c, &lambda_0, &lambda_1);

      if (norr == 0) { vortex_vector.push_back (point); }
      else if (lambda_0 * lambda_1 < 0)
      {
         seed_saddle (point, jacobian, lambda_0, lambda_1);
      }

   }

   for (Index_1D count = 0; count < max_count; count++)
   {
      for (Streamline_Set::iterator iterator = begin ();
           iterator != end (); iterator++)
      {
         Streamline& streamline = **(iterator);
         bool watch_ahead = true;
//         bool watch_ahead = (count > min_separatrix_count ||
//                             !streamline.is_separatrix);
         grow (streamline, watch_ahead);
      }
   }

   for (vector<Point_2D>::iterator iterator = vortex_vector.begin ();
        iterator != vortex_vector.end (); iterator++)
   {

      Point_2D& point = *(iterator);

      for (Index_1D i = 0; i < 10; i++)
      {

         Point_2D p = point + Point_2D ((i + 0.5) * seed_separation, 0);
         if (bounding_box.is_out_of_bounds (p)) { continue; }
         if (has_obstacle (p, nav, separation)) { continue; }

         Streamline* streamline_ptr = get_streamline_ptr (p);

         if (streamline_ptr->size () < 5) { delete streamline_ptr; }
         else { push_back (streamline_ptr); }

      }

   }

}

void
Streamline_Set::seed_saddle (const Point_2D& saddle,
                             const Jacobian& jacobian,
                             Scalar lambda_0,
                             Scalar lambda_1)
{

   Scalar thetas[4];
   thetas[0] = atan (jacobian.dv_dx / (lambda_0 - jacobian.dv_dy));
   thetas[1] = atan (jacobian.dv_dx / (lambda_1 - jacobian.dv_dy));
   thetas[2] = thetas[0] + M_PI;
   thetas[3] = thetas[1] + M_PI;

   Point_2D seed;
   bool separatrix = true;
   Streamline* streamline_ptrs[4];

   for (Index_1D i = 0; i < 4; i++)
   {

      Scalar cos_theta = cos (thetas[i]);
      Scalar sin_theta = sin (thetas[i]);

      seed.x = saddle.x + 1.5 * step_size * cos_theta;
      seed.y = saddle.y + 1.5 * step_size * sin_theta;

      Vector vector = vector_data_ptr->get_aligned_vector_at (seed);
      vector.normalize ();

      Scalar dot_product = vector.u * cos_theta + vector.v * sin_theta;

      if (fabs (dot_product) < 0.0)
      {
         separatrix = false;
         streamline_ptrs[i] = NULL;
      }

      else
      {
         bool forward = (dot_product > 0);
         streamline_ptrs[i] = new Streamline (saddle, seed, forward);
      }

   }

   for (Index_1D i = 0; i < 4; i++)
   {

      Streamline* streamline_ptr = streamline_ptrs[i];

      if (separatrix)
      {
         streamline_ptr->is_separatrix = true;
         add_to_cell (*streamline_ptr);
         push_back (streamline_ptr);
         saddle_vector.push_back (saddle);
      }
      else { if (streamline_ptr != NULL) { delete streamline_ptr; } }

   }

}

void
Streamline_Set::obtain_ordinary_streamlines ()
{

   Vector nav (GSL_NAN, GSL_NAN);
   Point_2D left_point, right_point;
   vector<Streamline*> streamline_ptr_vector;

   if (size () == 0)
   {
      Point_2D point;
      point.x = (bounding_box.start_x + bounding_box.end_x) / 2;
      point.y = (bounding_box.start_y + bounding_box.end_y) / 2;
      Streamline* streamline_ptr = get_streamline_ptr (point);
      push_back (streamline_ptr);
   }

   for (Streamline_Set::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Streamline* streamline_ptr = *(iterator);
      streamline_ptr_vector.push_back (streamline_ptr);
   }

   for (Index_1D i = 0; i < streamline_ptr_vector.size (); i++)
   {

      Streamline* streamline_ptr = streamline_ptr_vector[i];

      Size_1D streamline_size = streamline_ptr->size ();
      Point_2D point_array[streamline_size];

      {

         Index_1D pi = 0;

         for (Streamline::iterator iterator = streamline_ptr->begin ();
              iterator != streamline_ptr->end (); iterator++)
         {
            Point_2D& point = *(iterator);
            point_array[pi] = point;

            pi++;

         }

      }

      shuffle_array (point_array, streamline_size);

      for (Index_1D pi = 0; pi < streamline_size; pi++)
      {

         Point_2D& point = point_array[pi];
         Vector vector = vector_data_ptr->get_aligned_vector_at (point);
         Scalar theta = atan2 (vector.v, vector.u);

         for (Index_1D j = -1; j <= 1; j++)
         {

            Point_2D p;
            p.x = point.x + j * seed_separation * cos (theta + M_PI_2);
            p.y = point.y + j * seed_separation * sin (theta + M_PI_2);

            if (bounding_box.is_out_of_bounds (p)) { continue; }
            if (has_obstacle (p, nav, separation)) { continue; }

            Streamline* streamline_ptr = get_streamline_ptr (p);
            if (streamline_ptr->size () < 5) { delete streamline_ptr; }
            else
            {
               streamline_ptr_vector.push_back (streamline_ptr);
               push_back (streamline_ptr);
            }

         }

      }

   }

}

Index_2D
Streamline_Set::get_index (const Point_2D& point) const
{
   Scalar dx = point.x - bounding_box.start_x;
   Scalar dy = point.y - bounding_box.start_y;
   Index_1D i = Index_1D (rint (dx / delta_x));
   Index_1D j = Index_1D (rint (dy / delta_y));
   return Index_2D (i, j);
}

void
Streamline_Set::add_to_cell (const Point_2D& point) const
{
   if (bounding_box.is_out_of_bounds (point)) { return; }
   Index_2D index_2d = get_index (point);
   cell_data_ptr->get (index_2d).push_back (point);
}

void
Streamline_Set::add_to_cell (Streamline& streamline) const
{

   for (Streamline::iterator iterator = streamline.begin ();
        iterator != streamline.end (); iterator++)
   {
      Point_2D& point = *(iterator);
      add_to_cell (point);
   }

}

bool
Streamline_Set::has_obstacle (const Point_2D& point,
                              const Vector& vector,
                              Scalar threshold) const
{

   bool obstacle = false;
   Index_2D index_2d = get_index (point);

   if (vector.u == 0 && vector.v == 0) { obstacle = true; }

   for (Index_1D i = -1; i <= 1 && !obstacle; i++)
   {

      Index_1D ii = index_2d.i + i;
      if (ii < 0 || ii >= size_2d.i) { continue; }

      for (Index_1D j = -1; j <= 1 && !obstacle; j++)
      {

         Index_1D jj = index_2d.j + j;
         if (jj < 0 || jj >= size_2d.j) { continue; }

         Streamline_Cell& cell = cell_data_ptr->get (ii, jj);

         for (Streamline_Cell::iterator iterator = cell.begin ();
              iterator != cell.end () && !obstacle; iterator++)
         {

            Point_2D& p = *(iterator);
            Scalar dx = p.x - point.x;
            Scalar dy = p.y - point.y;

            if (dx == 0 && dy == 0) { continue; }

            if (vector.is_nav () || vector.u * dx + vector.v * dy > 0)
            {
               Scalar d = sqrt (dx*dx + dy*dy);
               if (d < threshold) { obstacle = true; }
            }

         }

      }
   }

   return obstacle;

}

void
Streamline_Set::grow (Streamline& streamline,
                      bool backward,
                      bool watch_ahead) const
{

   Point_2D point = (backward ? streamline.front () : streamline.back ());

   Vector vector = vector_data_ptr->get_aligned_vector_at (point);
   if (backward) { vector *= -1; }
   bool zero_vector = (vector.u == 0 && vector.v == 0);
   if (!zero_vector) { vector.normalize (); }

   bool out_of_bounds = bounding_box.is_out_of_bounds (point);
   bool obstacle = has_obstacle (point, vector, separation);

   if (point.is_nap () || out_of_bounds ||
       (watch_ahead && obstacle) || zero_vector)
   {
      bool &ok = (backward ? streamline.backward_ok : streamline.forward_ok);
      ok = false;
   }

   else
   {

      Scalar h = step_size;
      Scalar theta = atan2 (vector.v, vector.u);
      point += Point_2D (h * cos (theta), h * sin (theta));

      add_to_cell (point);
      if (backward) { streamline.push_front (point); }
      else          { streamline.push_back (point); }

   }

}

void
Streamline_Set::grow (Streamline& streamline,
                      bool watch_ahead) const
{
   if (streamline.forward_ok) { grow (streamline, false, watch_ahead); }
   if (streamline.backward_ok) { grow (streamline, true, watch_ahead); }
}

Streamline*
Streamline_Set::get_streamline_ptr (const Point_2D& point) const
{

   Streamline* s_ptr = new Streamline (point);

   for (Index_1D count = 0; count < max_count && s_ptr->open (); count++)
   {
      grow (*s_ptr, true);
   }

   return s_ptr;

}

Streamline_Set::Streamline_Set (Vector_Data_2D& vector_data,
                                const Box& bounding_box,
                                Scalar separation,
                                Size_1D step_size_granulity,
                                Scalar seed_separation_ratio,
                                Scalar max_streamline_length,
                                Scalar min_separatrix_length)
{

   this->separation = separation;
   this->bounding_box = bounding_box;
   this->vector_data_ptr = &vector_data;

   size_2d.i = Size_1D (rint (bounding_box.get_width () / separation));
   size_2d.j = Size_1D (rint (bounding_box.get_height () / separation));

   delta_x = bounding_box.get_width () / (size_2d.i - 1);
   delta_y = bounding_box.get_height () / (size_2d.j - 1);

   critical_point_set_ptr =
      vector_data_ptr->get_critical_point_set_ptr (separation);

   cell_data_ptr = new Data_2D<Streamline_Cell> (size_2d);

   step_size = separation / 4;
   max_count = Size_1D (rint (100 / step_size));
   min_separatrix_count = Size_1D (rint (3 / step_size));
   seed_separation = seed_separation_ratio * separation;

   obtain_special_streamlines ();
   obtain_ordinary_streamlines ();

}

Streamline_Set::~Streamline_Set ()
{

   delete cell_data_ptr;
   delete critical_point_set_ptr;

   for (Streamline_Set::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Streamline* streamline_ptr = *(iterator);
      delete streamline_ptr;
   }

}

const vector<Point_2D>&
Streamline_Set::get_saddle_vector () const
{
   return saddle_vector;
}

const vector<Point_2D>&
Streamline_Set::get_vortex_vector () const
{
   return vortex_vector;
}
*/

