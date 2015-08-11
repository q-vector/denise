	//
// visualize.cc
// 
// Copyright (C) 2010 Simon E. Ching
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

#include "visualize.h"

using namespace std;
using namespace denise;

Contour::Label_Point_Set::Label_Point_Set (const Real tolerance)
   : tolerance (tolerance)
{
}

void
Contour::Label_Point_Set::add (const Point_2D& point)
{
   insert (point);
}

bool
Contour::Label_Point_Set::interfere_with (const Point_2D& point)
{

   for (Label_Point_Set::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Point_2D& lp = *(iterator);
      const Real dx = lp.x - point.x;
      const Real dy = lp.y - point.y;

      if (sqrt (dx*dx + dy*dy) < tolerance)
      {
         return true;
      }

   }

   return false;

}

Contour::Isoline::Isoline (const Integer level_index)
   : level_index (level_index),
     segment_ptr_a (NULL),
     segment_ptr_b (NULL)
{
}

Contour::Boundary::Segment*
Contour::Isoline::get_segment_ptr_a () const
{
   return segment_ptr_a;
}

Contour::Boundary::Segment*
Contour::Isoline::get_segment_ptr_b () const
{
   return segment_ptr_b;
}

void
Contour::Isoline::set_segment_ptr_a (Boundary::Segment* segment_ptr)
{
   this->segment_ptr_a = segment_ptr;
}

void
Contour::Isoline::set_segment_ptr_b (Boundary::Segment* segment_ptr)
{
   this->segment_ptr_b = segment_ptr;
}

bool
Contour::Isoline::is_closed () const
{
   return (segment_ptr_a == NULL || segment_ptr_b == NULL);
}

const bool
Contour::Isoline::is_soiled () const
{
   return soiled;
}

void
Contour::Isoline::mark_soiled ()
{
   soiled = true;
}

void
Contour::Isoline::reset_soiled ()
{
   soiled = false;
}

Contour::Boundary::Segment::Segment (const Integer index,
                                     const Real value_a,
                                     const Real value_b,
                                     const Point_2D& point_a,
                                     const Point_2D& point_b,
                                     const Integer number_of_levels)
   : index (index),
     value_a (value_a),
     value_b (value_b),
     point_a (point_a),
     point_b (point_b)
{

   typedef Isoline* Isoline_Ptr;

   isoline_ptrs = new Isoline_Ptr[number_of_levels];
   heads = new bool[number_of_levels];

   for (Integer i = 0; i < number_of_levels; i++)
   {
      isoline_ptrs[i] = NULL;
   }

}

Contour::Boundary::Segment::~Segment ()
{
   delete[] isoline_ptrs;
   delete[] heads;
}

Contour::Isoline*
Contour::Boundary::Segment::get_isoline_ptr (const Integer level_index) const
{
   return isoline_ptrs[level_index];
}

bool
Contour::Boundary::Segment::connected_to_head (const Integer level_index) const
{
   return heads[level_index];
}

bool
Contour::Boundary::Segment::connected_to_head (const Isoline& isoline) const
{
   return connected_to_head (isoline.level_index);
}

void
Contour::Boundary::Segment::register_isoline (const Isoline& isoline,
                                              const bool head) const
{
   const Integer& level_index = isoline.level_index;
   isoline_ptrs[level_index] = (Isoline*)(&isoline);
   heads[level_index] = head;
}

void
Contour::Boundary::Segment::detach_isoline (const Integer level_index)
{
   isoline_ptrs[level_index] = NULL;
}

const Integer&
Contour::Boundary::Segment::get_index () const
{
   return index;
}

const Real&
Contour::Boundary::Segment::get_value_a () const
{
   return value_a;
}

const Real&
Contour::Boundary::Segment::get_value_b () const
{
   return value_b;
}

const Point_2D&
Contour::Boundary::Segment::get_point_a () const
{
   return point_a;
}

const Point_2D&
Contour::Boundary::Segment::get_point_b () const
{
   return point_b;
}

Contour::Boundary::Boundary (const Scalar_Data_2D& scalar_data_2d,
                             const Tuple& level_tuple)
{

   const Integer level_tuple_size = level_tuple.size ();
   const Tuple& tuple_x = scalar_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = scalar_data_2d.get_coordinate_tuple (1);
   const Integer ni = tuple_x.size ();
   const Integer nj = tuple_y.size ();

   for (Integer i = 0; i < ni - 1; i++)
   {
      const Real value_a = scalar_data_2d.get_datum (i, 0);
      const Real value_b = scalar_data_2d.get_datum (i + 1, 0);
      const Point_2D point_a (tuple_x[i], tuple_y.front ());
      const Point_2D point_b (tuple_x[i + 1], tuple_y.front ());
      const Integer index = segment_ptr_vector.size ();
      Boundary::Segment* segment_ptr = new Boundary::Segment (index,
         value_a, value_b, point_a, point_b, level_tuple_size);
      segment_ptr_vector.push_back (segment_ptr);
   }

   for (Integer j = 0; j < nj - 1; j++)
   {
      const Real value_a = scalar_data_2d.get_datum (ni - 1, j);
      const Real value_b = scalar_data_2d.get_datum (ni - 1, j + 1);
      const Point_2D point_a (tuple_x.back (), tuple_y[j]);
      const Point_2D point_b (tuple_x.back (), tuple_y[j + 1]);
      const Integer index = segment_ptr_vector.size ();
      Boundary::Segment* segment_ptr = new Boundary::Segment (index,
         value_a, value_b, point_a, point_b, level_tuple_size);
      segment_ptr_vector.push_back (segment_ptr);
   }

   for (Integer i = ni - 2; i >= 0; i--)
   {
      const Real value_a = scalar_data_2d.get_datum (i + 1, nj - 1);
      const Real value_b = scalar_data_2d.get_datum (i, nj - 1);
      const Point_2D point_a (tuple_x[i + 1], tuple_y.back ());
      const Point_2D point_b (tuple_x[i], tuple_y.back ());
      const Integer index = segment_ptr_vector.size ();
      Boundary::Segment* segment_ptr = new Boundary::Segment (index,
         value_a, value_b, point_a, point_b, level_tuple_size);
      segment_ptr_vector.push_back (segment_ptr);
   }

   for (Integer j = nj - 2; j >= 0; j--)
   {
      const Real value_a = scalar_data_2d.get_datum (0, j + 1);
      const Real value_b = scalar_data_2d.get_datum (0, j);
      const Point_2D point_a (tuple_x.front (), tuple_y[j + 1]);
      const Point_2D point_b (tuple_x.front (), tuple_y[j]);
      const Integer index = segment_ptr_vector.size ();
      Boundary::Segment* segment_ptr = new Boundary::Segment (index,
         value_a, value_b, point_a, point_b, level_tuple_size);
      segment_ptr_vector.push_back (segment_ptr);
   }

}

Contour::Boundary::~Boundary ()
{
   for (vector<Segment*>::iterator iterator = segment_ptr_vector.begin ();
        iterator != segment_ptr_vector.end (); iterator++)
   {
      Segment* segment_ptr = *(iterator);
      delete segment_ptr;
   }
}

Integer
Contour::Boundary::size () const
{
   return segment_ptr_vector.size ();
}

Contour::Boundary::Segment*
Contour::Boundary::get_segment_ptr (const Integer index) const
{
   return segment_ptr_vector.at (index);
}

void
Contour::Boundary::advance_segment_ptr (Segment*& segment_ptr,
                                        const bool backward) const
{
   const Integer increment = (backward ? -1 : 1);
   const Integer index = segment_ptr->get_index () + increment;
   segment_ptr = get_segment_ptr (imodulo (index, segment_ptr_vector.size ()));
}

Contour::Edge::Edge (const Integer level_index,
                     const Point_2D& point_a,
                     const Point_2D& point_b,
                     const Side side_a,
                     const Side side_b)
   : denise::Edge (point_a, point_b),
     level_index (level_index),
     side_a (side_a),
     side_b (side_b)
{
}

bool
Contour::Edge::operator== (const Contour::Edge& edge) const
{
   return (level_index == edge.level_index &&
           side_a == edge.side_a && side_b == edge.side_b &&
           point_a == edge.point_a && point_b == edge.point_b);
}

Contour::Side_Edge::Side_Edge (const Integer lower_level_index,
                               const Real position_a,
                               const Real position_b,
                               const Side side,
                               const Side_End side_end_a,
                               const Side_End side_end_b)
   : lower_level_index (lower_level_index),
     position_a (position_a),
     position_b (position_b),
     side (side),
     side_end_a (side_end_a),
     side_end_b (side_end_b)
{
}

bool
Contour::Side_Edge::operator== (const Side_Edge& side_edge) const
{

   const Side_Edge& se = side_edge;

   return (lower_level_index == se.lower_level_index &&
           side == se.side &&
           side_end_a == se.side_end_a &&
           side_end_b == se.side_end_b &&
           position_a == se.position_a &&
           position_b == se.position_b);

}

Contour::Cell::Cell (const Tuple& level_tuple,
                     const Integer i,
                     const Integer j)
   : Index_2D (i, j),
     level_tuple (level_tuple)
{

   const Integer n = level_tuple.size ();
   const Integer n_1 = (n > 0 ? n - 1 : 0);

   typedef list<Edge> Edge_List;
   typedef list<Side_Edge> Side_Edge_List;
   typedef list<Isoline*> Isoline_Ptr_List;

   this->edge_lists = new Edge_List[n];
   this->side_edge_lists = new Side_Edge_List[n_1];
   this->start_isoline_ptr_lists = new Isoline_Ptr_List[n];
   this->end_isoline_ptr_lists = new Isoline_Ptr_List[n];

}

Contour::Cell::~Cell ()
{
   delete[] edge_lists;
   delete[] side_edge_lists;
   delete[] start_isoline_ptr_lists;
   delete[] end_isoline_ptr_lists;
}

bool
Contour::Cell::has_edge (const Integer level_index) const
{
   const Integer& li = level_index;
   list<Edge>& edge_list = edge_lists[li];
   return (edge_list.size () > 0);
}

bool
Contour::Cell::has_edge (const Integer level_index,
                         const Side side) const
{
   const Integer& li = level_index;
   list<Edge>& edge_list = edge_lists[li];

   switch (edge_list.size ())
   {

      case 0:
      {
         return false;
      }

      case 1:
      {
         const Edge& edge = edge_list.front ();
         return (edge.side_a == side || edge.side_b == side);
      }

      default: // 2 edges
      {
         return true;
      }

   }

}

bool
Contour::Cell::has_side_edge (const Integer lower_level_index) const
{
   const Integer& lli = lower_level_index;
   const list<Side_Edge>& sel = side_edge_lists[lli];
   return (sel.size () > 0);
}

Contour::Edge
Contour::Cell::get_edge (const Integer level_index) const
{
   const Integer& li = level_index;
   list<Edge>& edge_list = edge_lists[li];
   Edge& edge = edge_list.front ();
   edge_list.remove (edge);
   return edge;
}

Contour::Edge
Contour::Cell::get_edge (const Integer level_index,
                         const Side side) const
{
   const Integer& li = level_index;
   list<Edge>& edge_list = edge_lists[li];

   if (edge_list.size () == 1)
   {
      const Contour::Edge& edge = edge_list.front ();
      edge_list.remove (edge);
      return edge;
   }
   else // 2 edges
   {

      const Contour::Edge& e_a = edge_list.front ();
      const Contour::Edge& e_b = edge_list.back ();

      const bool is_a = (e_a.side_a == side || e_a.side_b == side);

      const Contour::Edge& edge = (is_a ? e_a : e_b);
      edge_list.remove (edge);
      return edge;

   }

}

Contour::Side_Edge
Contour::Cell::get_side_edge (const Integer lower_level_index) const
{
   const Integer& lli = lower_level_index;
   list<Side_Edge>& sel = side_edge_lists[lli];
   Contour::Side_Edge side_edge = sel.front ();
   sel.remove (side_edge);
   return side_edge;
}

void
Contour::Cell::insert (const Edge& edge)
{
   const Integer li = edge.level_index;
   list<Edge>& edge_list = edge_lists[li];
   edge_list.push_back (edge);
}

void
Contour::Cell::insert (const Side_Edge& side_edge)
{
   const Integer lli = side_edge.lower_level_index;
   list<Side_Edge>& sel = side_edge_lists[lli];
   sel.push_back (side_edge);
}

void
Contour::Cell::register_start_isoline (Isoline& isoline)
{
   const Integer li = isoline.level_index;
   list<Isoline*>& ipl = start_isoline_ptr_lists[li];
   ipl.push_back (&isoline);
}

void
Contour::Cell::register_end_isoline (Isoline& isoline)
{
   const Integer li = isoline.level_index;
   list<Isoline*>& ipl = end_isoline_ptr_lists[li];
   ipl.push_back (&isoline);
}

Contour::Valid_Edge::Valid_Edge (const Integer i,
                                 const Integer j,
                                 const bool along_x)
   : Index_2D (i, j),
     along_x (along_x)
{
}

bool
Contour::Valid_Edge::operator== (const Valid_Edge& valid_edge) const
{
   return ((Index_2D (*this) == Index_2D (valid_edge)) &&
           along_x == valid_edge.along_x);
}

bool
Contour::Valid_Edge::operator> (const Valid_Edge& valid_edge) const
{
   if (Index_2D (*this) > Index_2D (valid_edge)) { return true; }
   else if (Index_2D (*this) < Index_2D (valid_edge)) { return false; }
   else { return ((along_x ? 1 : 0) > (valid_edge.along_x ? 1 : 0)); }
}

bool
Contour::Valid_Edge::operator< (const Valid_Edge& valid_edge) const
{
   if (Index_2D (*this) < Index_2D (valid_edge)) { return true; }
   else if (Index_2D (*this) > Index_2D (valid_edge)) { return false; }
   else { return ((along_x ? 1 : 0) < (valid_edge.along_x ? 1 : 0)); }
}

bool
Contour::Valid_Edge_Set::contains (const Valid_Edge& valid_edge) const
{
   Contour::Valid_Edge_Set::const_iterator iterator = find (valid_edge);
   return (iterator != end ());
}

bool
Contour::Valid_Edge_Set::contains (const Integer i,
                                   const Integer j,
                                   const bool along_x) const
{
   return contains (Valid_Edge (i, j, along_x));
}

Contour::Valid_Data::Valid_Data (const Scalar_Data_2D& scalar_data_2d)
{

   const Tuple& tuple_x = scalar_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = scalar_data_2d.get_coordinate_tuple (1);
   size_2d.i = tuple_x.size () - 1;
   size_2d.j = tuple_y.size () - 1;

   buffer = new bool[size_2d.i * size_2d.j];

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real ne = scalar_data_2d.get_datum (i + 1, j + 1);
         const Real nw = scalar_data_2d.get_datum (i, j + 1);
         const Real se = scalar_data_2d.get_datum (i + 1, j);
         const Real sw = scalar_data_2d.get_datum (i, j);

         const bool ne_valid = (gsl_isnan (ne));
         const bool nw_valid = (gsl_isnan (nw));
         const bool se_valid = (gsl_isnan (se));
         const bool sw_valid = (gsl_isnan (sw));

         const bool valid = (ne_valid || nw_valid || se_valid || sw_valid);
//         const bool valid = (i > 5 && i < 20 && j > 15 && j < 30);
//         const bool valid = (i < 5 && j < 5);

         set (i, j, valid);

      }
   }

}

Contour::Valid_Data::~Valid_Data ()
{
   delete[] buffer;
}

const Size_2D&
Contour::Valid_Data::get_size_2d () const
{
   return size_2d;
}

bool
Contour::Valid_Data::get (const Integer i,
                          const Integer j) const
{
   if (i < 0 || i >= size_2d.i || j < 0 || j >= size_2d.j) { return false; }
   return buffer[i * size_2d.j + j];
}

void
Contour::Valid_Data::set (const Integer i,
                          const Integer j,
                          const bool b) const
{
   buffer[i * size_2d.j + j] = b;
}

bool
Contour::Valid_Data::is_valid_edge_x (const Integer i,
                                      const Integer j) const
{
   return get (i, j) != get (i - 1, j);
}

bool
Contour::Valid_Data::is_valid_edge_y (const Integer i,
                                      const Integer j) const
{
   return get (i, j) != get (i, j - 1);
}

bool
Contour::Valid_Data::is_valid (const Valid_Edge& valid_edge) const
{
   const Integer& i = valid_edge.i;
   const Integer& j = valid_edge.j;
   const bool along_x = valid_edge.along_x;
   return (along_x ? is_valid_edge_x (i, j) : is_valid_edge_y (i, j));
}

Contour::Valid_Path::Valid_Path (const Integer i,
                                 const Integer j,
                                 const bool along_x,
                                 const bool forward)
   : Valid_Edge (i, j, along_x),
     forward (forward)
{
}

Contour::Valid_Path::Valid_Path (const Valid_Edge& valid_edge,
                                 const bool forward)
   : Valid_Edge (valid_edge),
     forward (forward)
{
}

void
Contour::Valid_Path::set (const Integer i,
                          const Integer j,
                          const bool along_x,
                          const bool forward)
{
   this->i = i;
   this->j = j;
   this->along_x = along_x;
   this->forward = forward;
}

bool
Contour::Valid_Path::should_start_trace (const Valid_Data& valid_data,
                                         const Valid_Edge_Set& ves) const
{
   const bool not_processed_yet = !ves.contains (*this);
   return (valid_data.is_valid (*this) && not_processed_yet);
}

Contour::Valid_Path
Contour::Valid_Polygon::get_valid_path_a (const Valid_Path& valid_path) const
{

   const Integer& i = valid_path.i;
   const Integer& j = valid_path.j;
   const bool along_x = valid_path.along_x;
   const bool along_y = !along_x;
   const bool forward = valid_path.forward;
   const bool backward = !forward;

   if (along_x && forward) { return Valid_Path (i, j+1, false, true); }
   else if (along_x && backward) { return Valid_Path (i-1, j, false, false); }
   else if (along_y && forward) { return Valid_Path (i+1, j-1, true, false); }
   else if (along_y && backward) { return Valid_Path (i, j, true, true); }

}

Contour::Valid_Path
Contour::Valid_Polygon::get_valid_path_b (const Valid_Path& valid_path) const
{

   const Integer& i = valid_path.i;
   const Integer& j = valid_path.j;
   const bool along_x = valid_path.along_x;
   const bool along_y = !along_x;
   const bool forward = valid_path.forward;
   const bool backward = !forward;

   if (along_x && forward) { return Valid_Path (i, j+1, true, true); }
   else if (along_x && backward) { return Valid_Path (i, j-1, true, false); }
   else if (along_y && forward) { return Valid_Path (i+1, j, false, true); }
   else if (along_y && backward) { return Valid_Path (i-1, j, false, false); }

}

Contour::Valid_Path
Contour::Valid_Polygon::get_valid_path_c (const Valid_Path& valid_path) const
{

   const Integer& i = valid_path.i;
   const Integer& j = valid_path.j;
   const bool along_x = valid_path.along_x;
   const bool along_y = !along_x;
   const bool forward = valid_path.forward;
   const bool backward = !forward;

   if (along_x && forward) { return Valid_Path (i-1, j+1, false, false); }
   else if (along_x && backward) { return Valid_Path (i, j, false, true); }
   else if (along_y && forward) { return Valid_Path (i+1, j, true, true); }
   else if (along_y && backward) { return Valid_Path (i, j-1, true, false); }

}

Point_2D
Contour::Valid_Polygon::get_point_head (const Scalar_Data_2D& scalar_data_2d,
                                        const Valid_Path& valid_path) const
{

   const Integer& i = valid_path.i;
   const Integer& j = valid_path.j;
   const bool along_x = valid_path.along_x;
   const bool along_y = !along_x;
   const bool forward = valid_path.forward;
   const bool backward = !forward;

   Integer si, sj;

   if (along_x && forward) { si = i; sj = j; }
   else if (along_x && backward) { si = i; sj = j+1; }
   else if (along_y && forward) { si = i; sj = j; }
   else if (along_y && backward) { si = i+1; sj = j; }

   const Real x = scalar_data_2d.get_coordinate (0, si);
   const Real y = scalar_data_2d.get_coordinate (1, sj);
   return Point_2D (x, y);

}

Point_2D
Contour::Valid_Polygon::get_point_tail (const Scalar_Data_2D& scalar_data_2d,
                                        const Valid_Path& valid_path) const
{

   const Integer& i = valid_path.i;
   const Integer& j = valid_path.j;
   const bool along_x = valid_path.along_x;
   const bool along_y = !along_x;
   const bool forward = valid_path.forward;
   const bool backward = !forward;

   Integer si, sj;

   if (along_x && forward) { si = i; sj = j+1; }
   else if (along_x && backward) { si = i; sj = j; }
   else if (along_y && forward) { si = i+1; sj = j; }
   else if (along_y && backward) { si = i; sj = j; }

   const Real x = scalar_data_2d.get_coordinate (0, si);
   const Real y = scalar_data_2d.get_coordinate (1, sj);
   return Point_2D (x, y);

}

void
Contour::Valid_Polygon::trace (const Scalar_Data_2D& scalar_data_2d,
                               const Valid_Data& valid_data,
                               const Valid_Path& valid_path,
                               const bool new_handle)
{

   const Scalar_Data_2D& data_2d = scalar_data_2d;
   const Valid_Data& vd = valid_data;

   const Integer& i = valid_path.i;
   const Integer& j = valid_path.j;
   const bool along_x = valid_path.along_x;
   const bool forward = valid_path.forward;

   // should not start trace here, do nothing
   if (!valid_path.should_start_trace (vd, valid_edge_set)) { return; }

   valid_edge_set.insert (valid_path);
   
   if (new_handle) { add (get_point_head (data_2d, valid_path), true); }
   add (get_point_tail (data_2d, valid_path), false);

   const Valid_Path valid_path_a = get_valid_path_a (valid_path);
   const Valid_Path valid_path_b = get_valid_path_b (valid_path);
   const Valid_Path valid_path_c = get_valid_path_c (valid_path);

   if (vd.is_valid (valid_path_a)) { trace (data_2d, vd, valid_path_a); }
   else if (vd.is_valid (valid_path_b)) { trace (data_2d, vd, valid_path_b); }
   else if (vd.is_valid (valid_path_c)) { trace (data_2d, vd, valid_path_c); }

}

Contour::Valid_Polygon::Valid_Polygon (const Scalar_Data_2D& scalar_data_2d)
{

   const Valid_Data valid_data (scalar_data_2d);
   const Size_2D& size_2d = valid_data.get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Valid_Path valid_path (i, j, true, true);
         trace (scalar_data_2d, valid_data, valid_path, true);
      }
   }

}

Contour::Side
Contour::opposite (const Side side)
{

   Side opposite_side;

   switch (side)
   {
      case SIDE_N: opposite_side = SIDE_S; break;
      case SIDE_S: opposite_side = SIDE_N; break;
      case SIDE_E: opposite_side = SIDE_W; break;
      case SIDE_W: opposite_side = SIDE_E; break;
   }

   return opposite_side;

}

Contour::Cell&
Contour::get_cell (const Integer i,
                   const Integer j) const
{
   const Integer nj = scalar_data_2d_ptr->get_coordinate_tuple (1).size () - 1;
   return *(cell_ptrs[i * nj + j]);
}

Contour::Cell**
Contour::get_cell_ptrs () const
{

   const Tuple& tuple_x = scalar_data_2d_ptr->get_coordinate_tuple (0);
   const Tuple& tuple_y = scalar_data_2d_ptr->get_coordinate_tuple (1);

   const Integer ni = tuple_x.size () - 1;
   const Integer nj = tuple_y.size () - 1;

   typedef Cell* Cell_Ptr;
   Cell_Ptr* cell_ptrs = new Cell_Ptr[ni * nj];

   for (Integer i = 0; i < tuple_x.size (); i++)
   {
      for (Integer j = 0; j < tuple_y.size (); j++)
      {
         if (i != ni && j != nj)
         {
            const Tuple& clt = level_tuple;
            Cell* cell_ptr = new Cell (clt, i, j);
            cell_ptrs[i * nj + j] = cell_ptr;
         }
      }
   }

   return cell_ptrs;

}

void
Contour::substitude_nan (const Real substitude)
{

   const Size_2D& size_2d = scalar_data_2d_ptr->get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         Real& datum = scalar_data_2d_ptr->get_datum (i, j);
         if (gsl_isnan (datum)) { datum = substitude; }
      }
   }

}

void
Contour::work_out_edges ()
{
   for (Integer i = 0; i < level_tuple.size (); i++)
   {
      work_out_edges (i);
   }
}

void
Contour::work_out_side_edges ()
{
   if (level_tuple.size () == 0) { return; }
   for (Integer i = 0; i < level_tuple.size () - 1; i++)
   {
      work_out_side_edges (i, i + 1);
   }
}

void
Contour::trace_isolines ()
{

   const Integer n = level_tuple.size ();
   typedef list<Isoline*> Isoline_Ptr_List;
   this->isoline_ptr_lists = new Isoline_Ptr_List[n];

   for (Integer i = 0; i < n; i++)
   {
      trace_isolines (i);
   }

}

void
Contour::work_out_edges (const Integer level_index)
{

   const Integer li = level_index;
   const Real level = level_tuple[li];

   const bool outside = (((level - min_value) * (level - max_value)) > 0);
   if (outside) { return; }

   const Tuple& tuple_x = scalar_data_2d_ptr->get_coordinate_tuple (0);
   const Tuple& tuple_y = scalar_data_2d_ptr->get_coordinate_tuple (1);

   if (tuple_x.size () <= 1 || tuple_y.size () <= 1) { return; }

   for (Integer i = 0; i < tuple_x.size () - 1; i++)
   {
      for (Integer j = 0; j < tuple_y.size () - 1; j++)
      {

         Cell& cell = get_cell (i, j);

         const Real n = tuple_y[j + 1];
         const Real s = tuple_y[j];
         const Real e = tuple_x[i + 1];
         const Real w = tuple_x[i];

         const Real delta_x = e - w;
         const Real delta_y = n - s;

         const Real ne = scalar_data_2d_ptr->get_datum (i + 1, j + 1);
         const Real nw = scalar_data_2d_ptr->get_datum (i, j + 1);
         const Real se = scalar_data_2d_ptr->get_datum (i + 1, j);
         const Real sw = scalar_data_2d_ptr->get_datum (i, j);

         const Real residue_ne = level - ne;
         const Real residue_nw = level - nw;
         const Real residue_se = level - se;
         const Real residue_sw = level - sw;

         const bool cut_n = (residue_ne * residue_nw < 0);
         const bool cut_s = (residue_sw * residue_se < 0);
         const bool cut_e = (residue_se * residue_ne < 0);
         const bool cut_w = (residue_sw * residue_nw < 0);

         if (cut_n && cut_s && cut_e && cut_w)
         {

            const Real middle_value = (ne + nw + se + sw) / 4;
            const bool b = ((middle_value - level) * residue_ne) < 0;

            const Point_2D p_aa (residue_nw / (ne - nw) * delta_x + w, n);
            const Point_2D p_ab (residue_sw / (se - sw) * delta_x + w, s);
            const Point_2D p_ba (e, residue_se / (ne - se) * delta_y + s);
            const Point_2D p_bb (w, residue_sw / (nw - sw) * delta_y + s);

            if (((level - middle_value) * residue_ne) < 0)
            {
               cell.insert (Contour::Edge (li, p_aa, p_ba, SIDE_N, SIDE_E));
               cell.insert (Contour::Edge (li, p_ab, p_bb, SIDE_S, SIDE_W));
            }
            else
            {
               cell.insert (Contour::Edge (li, p_aa, p_bb, SIDE_N, SIDE_W));
               cell.insert (Contour::Edge (li, p_ab, p_ba, SIDE_S, SIDE_E));
            }

         }
         else
         if (cut_n && cut_s)
         {
            const Point_2D point_a (residue_nw / (ne - nw) * delta_x + w, n);
            const Point_2D point_b (residue_sw / (se - sw) * delta_x + w, s);
            cell.insert (Contour::Edge (li, point_a, point_b, SIDE_N, SIDE_S));
         }
         else
         if (cut_n && cut_e)
         {
            const Point_2D point_a (residue_nw / (ne - nw) * delta_x + w, n);
            const Point_2D point_b (e, residue_se / (ne - se) * delta_y + s);
            cell.insert (Contour::Edge (li, point_a, point_b, SIDE_N, SIDE_E));
         }
         else
         if (cut_n && cut_w)
         {
            const Point_2D point_a (residue_nw / (ne - nw) * delta_x + w, n);
            const Point_2D point_b (w, residue_sw / (nw - sw) * delta_y + s);
            cell.insert (Contour::Edge (li, point_a, point_b, SIDE_N, SIDE_W));
         }
         else
         if (cut_s && cut_e)
         {
            const Point_2D point_a (residue_sw / (se - sw) * delta_x + w, s);
            const Point_2D point_b (e, residue_se / (ne - se) * delta_y + s);
            cell.insert (Contour::Edge (li, point_a, point_b, SIDE_S, SIDE_E));
         }
         else
         if (cut_s && cut_w)
         {
            const Point_2D point_a (residue_sw / (se - sw) * delta_x + w, s);
            const Point_2D point_b (w, residue_sw / (nw - sw) * delta_y + s);
            cell.insert (Contour::Edge (li, point_a, point_b, SIDE_S, SIDE_W));
         }
         else
         if (cut_e && cut_w)
         {
            const Point_2D point_a (e, residue_se / (ne - se) * delta_y + s);
            const Point_2D point_b (w, residue_sw / (nw - sw) * delta_y + s);
            cell.insert (Contour::Edge (li, point_a, point_b, SIDE_E, SIDE_W));
         }

      }
   }

}

void
Contour::work_out_side_edges (const Integer lower_level_index,
                              const Integer upper_level_index)
{

   const Tuple& tuple_x = scalar_data_2d_ptr->get_coordinate_tuple (0);
   const Tuple& tuple_y = scalar_data_2d_ptr->get_coordinate_tuple (1);
   const Integer ni = tuple_x.size ();
   const Integer nj = tuple_y.size ();

   if (ni < 2 || nj < 2) { return; }

   // N Side
   for (Integer i = 0; i < ni - 1; i++)
   {

      Cell& cell = get_cell (i, 0);

      const Real left = tuple_x[i];
      const Real right = tuple_x[i + 1];

      const Real l_value = scalar_data_2d_ptr->get_datum (i, 0);
      const Real r_value = scalar_data_2d_ptr->get_datum (i + 1, 0);

      work_out_side_edge (cell, SIDE_N, left, right, l_value,
         r_value, lower_level_index, upper_level_index);

   }

   // S Side
   for (Integer i = 0; i < ni - 1; i++)
   {

      Cell& cell = get_cell (i, nj - 2);

      const Real left = tuple_x[i + 1];
      const Real right = tuple_x[i];

      const Real l_value = scalar_data_2d_ptr->get_datum (i + 1, nj - 1);
      const Real r_value = scalar_data_2d_ptr->get_datum (i, nj - 1);

      work_out_side_edge (cell, SIDE_S, left, right, l_value,
         r_value, lower_level_index, upper_level_index);

   }

   // W Side
   for (Integer j = 0; j < nj - 1; j++)
   {

      Cell& cell = get_cell (0, j);

      const Real left = tuple_y[j + 1];
      const Real right = tuple_y[j];

      const Real l_value = scalar_data_2d_ptr->get_datum (0, j + 1);
      const Real r_value = scalar_data_2d_ptr->get_datum (0, j);

      work_out_side_edge (cell, SIDE_W, left, right, l_value,
         r_value, lower_level_index, upper_level_index);

   }

   // E Side
   for (Integer j = 0; j < nj - 1; j++)
   {

      Cell& cell = get_cell (ni - 2, j);

      const Real left = tuple_y[j];
      const Real right = tuple_y[j + 1];

      const Real l_value = scalar_data_2d_ptr->get_datum (ni - 1, j);
      const Real r_value = scalar_data_2d_ptr->get_datum (ni - 1, j + 1);

      work_out_side_edge (cell, SIDE_E, left, right, l_value,
         r_value, lower_level_index, upper_level_index);

   }

}

void
Contour::work_out_side_edge (Cell& cell,
                             const Side side,
                             const Real left,
                             const Real right,
                             const Real l_value,
                             const Real r_value,
                             const Integer lower_level_index,
                             const Integer upper_level_index)
{

   const Integer lli = lower_level_index;
   const Integer uli = upper_level_index;
   const Real lower_level = level_tuple[lli];
   const Real upper_level = level_tuple[uli];
   const Real& lcl = lower_level;
   const Real& ucl = upper_level;

   const bool outside_l = (((lcl - min_value) * (lcl - max_value)) > 0);
   const bool outside_u = (((ucl - min_value) * (ucl - max_value)) > 0);
   if (outside_l && outside_u)
   {
      const Real centre = (max_value + min_value) / 2;
      if (((centre - lcl) * (centre - ucl)) > 0) { return; }
   }

   if (l_value > lcl && l_value < ucl && r_value > lcl && r_value < ucl)
   {
      cell.insert (Side_Edge (lli, left, right, side, SE_LEFT, SE_RIGHT));
   }

/*
   const Real delta = right - left;
   const Real delta_value = r_value - l_value;

   if (l_value < lcl && r_value > lcl && r_value < ucl)
   {
      const Real lower = (lcl - l_value) / delta_value * delta + left;
      //cell.insert (Side_Edge (lli, lower, right, side, SE_LOWER, SE_RIGHT));
   }
   else
   if (l_value < lcl && r_value > ucl)
   {
      const Real lower = (lcl - l_value) / delta_value * delta + left;
      const Real upper = (ucl - l_value) / delta_value * delta + left;
      //cell.insert (Side_Edge (lli, lower, upper, side, SE_LOWER, SE_UPPER));
   }
   else
   if (l_value > lcl && l_value < ucl && r_value < lcl)
   {
      const Real lower = (lcl - l_value) / delta_value * delta + left;
      //cell.insert (Side_Edge (lli, left, lower, side, SE_LEFT, SE_LOWER));
   }
   else
   if (l_value > lcl && l_value < ucl && r_value > lcl && r_value < ucl)
   {
      cell.insert (Side_Edge (lli, left, right, side, SE_LEFT, SE_RIGHT));
   }
   else
   if (l_value > lcl && l_value < ucl && r_value > ucl)
   {
      const Real upper = (ucl - l_value) / delta_value * delta + left;
      //cell.insert (Side_Edge (lli, left, upper, side, SE_LEFT, SE_UPPER));
   }
   else
   if (l_value > ucl && r_value < lcl)
   {
      const Real lower = (lcl - l_value) / delta_value * delta + left;
      const Real upper = (ucl - l_value) / delta_value * delta + left;
      //cell.insert (Side_Edge (lli, upper, lower, side, SE_UPPER, SE_LOWER));
   }
   else
   if (l_value > ucl && l_value > lcl && r_value < ucl)
   {
      const Real upper = (ucl - l_value) / delta_value * delta + left;
      //cell.insert (Side_Edge (lli, upper, right, side, SE_UPPER, SE_RIGHT));
   }
*/

}

void
Contour::work_out_polygon (const Integer lower_level_index)
{

   Polygon& polygon = *(polygon_ptrs[lower_level_index]);
   const Integer upper_level_index = lower_level_index + 1;

   const Integer lli = lower_level_index;
   const Integer uli = upper_level_index;
   const Real lcl = level_tuple[lli];
   const Real ucl = level_tuple[uli];

   list<Isoline*> isoline_ptr_list;
   list<Isoline*>& lower_isoline_ptr_list = isoline_ptr_lists[lli];
   list<Isoline*>& upper_isoline_ptr_list = isoline_ptr_lists[uli];

   if (isoline_ptr_lists == NULL) { return; }

   // Closed lower isolines
   for (list<Isoline*>::iterator iterator = lower_isoline_ptr_list.begin ();
        iterator != lower_isoline_ptr_list.end (); iterator++)
   {
      Isoline* isoline_ptr = *(iterator);
      isoline_ptr->reset_soiled ();
      isoline_ptr_list.push_back (isoline_ptr);
   }

   // Closed upper isolines
   for (list<Isoline*>::iterator iterator = upper_isoline_ptr_list.begin ();
        iterator != upper_isoline_ptr_list.end (); iterator++)
   {
      Isoline* isoline_ptr = *(iterator);
      isoline_ptr->reset_soiled ();
      isoline_ptr_list.push_back (isoline_ptr);
   }

   bool boundary_processed = false;

   for (list<Isoline*>::iterator iterator = isoline_ptr_list.begin ();
        iterator != isoline_ptr_list.end (); iterator++)
   {

      bool forward = true;
      Isoline* isoline_ptr = *(iterator);

      if (isoline_ptr->is_soiled ())
      {
         continue;
      }

      if (isoline_ptr->is_closed ())
      {
         isoline_ptr->mark_soiled ();
         follow_isoline (polygon, *isoline_ptr, forward, true);
         continue;
      }

      isoline_ptr->mark_soiled ();
      boundary_processed = true;
      bool new_handle = true;

      do
      {

         Boundary::Segment* s_ptr =
            follow_isoline (polygon, *isoline_ptr, forward, new_handle);

         isoline_ptr->mark_soiled ();
         new_handle = false;

         // preliminary search: another isoline in the same boundary segment?
         Integer li = isoline_ptr->level_index;

         li = (li == lli ? uli : lli);
         isoline_ptr = s_ptr->get_isoline_ptr (li);

         if (isoline_ptr != NULL)
         {
            // We found one!
            //    if soiled we've come back a loop. stop and proceed
            //    else follow that isoline by next iteration
            if (isoline_ptr->is_soiled ()) { break; }
            else
            {
               // But we still have to deterine this isoline on the same
               // boundary segemnt should go forward or backward
               forward = s_ptr->connected_to_head (*isoline_ptr);
               continue;
            }
         }

         // preliminary search == FAIL
         // now we have to follow the side-edges till we get to another one
         const Real value_a = s_ptr->get_value_a ();
         const Real value_b = s_ptr->get_value_b ();
         const bool backward = ((value_a - lcl) * (value_a - ucl) < 0);
         //const bool backward = ((value_b - lcl) * (value_b - ucl) < 0);

         while (isoline_ptr == NULL)
         {

            // polygon add point on bs depending on _backward_

            boundary_ptr->advance_segment_ptr (s_ptr, backward);
            const Point_2D& point_a = s_ptr->get_point_a ();
            const Point_2D& point_b = s_ptr->get_point_b ();
            const Point_2D& point = (backward ? point_a : point_b);

            polygon.add (point, false);

            Isoline* lower_isoline_ptr = s_ptr->get_isoline_ptr (lli);
            Isoline* upper_isoline_ptr = s_ptr->get_isoline_ptr (uli);
 
            if (lower_isoline_ptr != NULL)
            { 
               isoline_ptr = lower_isoline_ptr;
               forward = s_ptr->connected_to_head (*isoline_ptr);
            }
            else
            if (upper_isoline_ptr != NULL)
            {
               isoline_ptr = upper_isoline_ptr;
               forward = s_ptr->connected_to_head (*isoline_ptr);
            }

         }

      }
      while (!isoline_ptr->is_soiled ());

   }

   if (!boundary_processed)
   {

      const Real corner_datum = scalar_data_2d_ptr->get_datum (0, 0);

      if ((corner_datum - lcl) * (corner_datum - ucl) <= 0)
      {

         const Integer n = boundary_ptr->size ();

         for (Integer i = 0; i < n; i++)
         {
            Boundary::Segment* s_ptr = boundary_ptr->get_segment_ptr (i);
            if (i == 0) { polygon.add (s_ptr->get_point_a (), true); }
            polygon.add (s_ptr->get_point_b (), false);
         }

      }

   }

}

Contour::Boundary::Segment*
Contour::follow_isoline (Polygon& polygon,
                         const Isoline& isoline,
                         const bool forward,
                         const bool new_handle) const
{

   if (forward)
   {

      for (Isoline::const_iterator iterator = isoline.begin ();
           iterator != isoline.end (); iterator++)
      {
         const Point_2D& point = *(iterator);
         const bool nh = (new_handle && iterator == isoline.begin ());
         polygon.add (point, nh);
      }

      return isoline.get_segment_ptr_b ();

   }
   else
   {

      for (Isoline::const_reverse_iterator iterator = isoline.rbegin ();
           iterator != isoline.rend (); iterator++)
      {
         const Point_2D& point = *(iterator);
         const bool nh = (new_handle && iterator == isoline.rbegin ());
         polygon.add (point, nh);
      }

      return isoline.get_segment_ptr_a ();

   }

}

void
Contour::trace_isolines (const Integer level_index)
{

   const Integer li = level_index;

   const Scalar_Data_2D& scalar_data_2d = *scalar_data_2d_ptr;
   const Tuple& tuple_x = scalar_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = scalar_data_2d.get_coordinate_tuple (1);

   const Integer ni = tuple_x.size () - 1;
   const Integer nj = tuple_y.size () - 1;

   if (ni < 2 || nj < 2) { return; }

   for (Integer i = 0; i < ni; i++)
   {
      for (Integer j = 0; j < nj; j++)
      {

         Cell& cell = get_cell (i, j);

         while (cell.has_edge (li))
         {

            Edge edge = cell.get_edge (li);
            const Side& side_a = edge.side_a;
            const Side& side_b = edge.side_b;
            const Point_2D& point_a = edge.point_a;
            const Point_2D& point_b = edge.point_b;

            Isoline* isoline_ptr = new Isoline (li);
            list<Isoline*>& isoline_ptr_list = isoline_ptr_lists[li];
            isoline_ptr_list.push_back (isoline_ptr);

            Isoline& isoline = *isoline_ptr;
            isoline.push_back (point_a);
            isoline.push_front (point_b);

            if (!trace_isoline (isoline, edge, i, j, side_a, true))
            {
               trace_isoline (isoline, edge, i, j, side_b, false);
            }

         }

      }
   }

}

bool
Contour::trace_isoline (Isoline& isoline,
                        const Edge& edge,
                        const Integer i,
                        const Integer j,
                        const Side side,
                        const bool forward)
{

   Integer next_i = i;
   Integer next_j = j;

   const Integer li = edge.level_index;
   const Scalar_Data_2D& scalar_data_2d = *scalar_data_2d_ptr;
   const Tuple& tuple_x = scalar_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = scalar_data_2d.get_coordinate_tuple (1);

   const Integer ni = tuple_x.size ();
   const Integer nj = tuple_y.size ();

   bool out_of_bounds = false;

   switch (side)
   {
      case SIDE_N: if (j == nj - 2) { out_of_bounds = true; } next_j++; break;
      case SIDE_S: if (j == 0) { out_of_bounds = true; } next_j--; break;
      case SIDE_E: if (i == ni - 2) { out_of_bounds = true; } next_i++; break;
      case SIDE_W: if (i == 0) { out_of_bounds = true; } next_i--; break;
   }

   if (out_of_bounds)
   {

      // isoline is open

      Integer index;

      switch (side)
      {
         case SIDE_S: index = i; break;
         case SIDE_E: index = ni + j - 1; break;
         case SIDE_N: index = 2 * ni + nj - i - 4; break;
         case SIDE_W: index = 2 * (ni + nj) - j - 5; break;
      }

      Boundary::Segment* segment_ptr = boundary_ptr->get_segment_ptr (index);

      if (forward)
      {
         segment_ptr->register_isoline (isoline, false);
         isoline.set_segment_ptr_b (segment_ptr);
      }
      else
      {
         segment_ptr->register_isoline (isoline, true);
         isoline.set_segment_ptr_a (segment_ptr);
      }

      return false;

   }

   // Advance to next cell

   Cell& cell = get_cell (next_i, next_j);
   const Side opposite_side = opposite (side);

   // isoline is closed
   // do not have to go in another direction
   if (!cell.has_edge (li, opposite_side)) { return true; }

   Contour::Edge new_edge = cell.get_edge (li, opposite_side);
   
   const Side& side_a = new_edge.side_a;
   const Side& side_b = new_edge.side_b;
   const Point_2D& point_a = new_edge.point_a;
   const Point_2D& point_b = new_edge.point_b;

   const bool a_to_b = (opposite_side == side_a);
   const Side next_side = (a_to_b ? side_b : side_a);
   const Point_2D& next_point = (a_to_b ? point_b : point_a);

   if (forward) { isoline.push_back (next_point); }
   else { isoline.push_front (next_point); }

   return trace_isoline (isoline, new_edge, next_i, next_j, next_side, forward);

}

Polygon*
Contour::render_label (const RefPtr<Context>& cr,
                       const Transform_2D& transform,
                       const Dstring& format,
                       const Real label_multiplier,
                       const Real label_offset,
                       const Real label_distance,
                       const Integer label_stride) const
{

   const Tuple& clt = level_tuple;
   const Integer n = clt.size ();

//   Polygon* clip_polygon_ptr = new Polygon;
   Polygon* clip_polygon_ptr = NULL;
   Label_Point_Set label_point_set (label_distance);

   //cr->save ();
   //Color (1, 0, 0, 0.5).cairo (cr);

   for (Integer li = 0; li < n; li++)
   {

      const list<Isoline*>& ipl = isoline_ptr_lists[li];

      for (list<Isoline*>::const_iterator iterator = ipl.begin ();
           iterator != ipl.end (); iterator++)
      {

         Point_2D last_point;
         const Isoline& isoline = *(*(iterator));
         const Integer li = isoline.level_index;
         const Real cl = clt[li] * label_multiplier - label_offset;

         for (Isoline::const_iterator iterator = isoline.begin ();
              iterator != isoline.end (); iterator++)
         {

            const Integer d = distance (iterator, isoline.begin ());
            const Point_2D& this_point = *(iterator);

            if (!transform.is_out_of_domain (this_point))
            {

               if (d % label_stride == label_stride / 2)
               {

                  const Point_2D this_p = transform.transform (this_point);

                  if (this_p.is_nap ()) { continue; }
                  if (label_point_set.interfere_with (this_p)) { continue; }

                  const Point_2D last_p = transform.transform (last_point);
                  if (last_p.is_nap ()) { continue; }

                  const Real dx = this_p.x - last_p.x;
                  const Real dy = this_p.y - last_p.y;
                  const Real theta = atan (dy / dx);

                  if (gsl_isnan (theta)) { continue; }
                
                  if (format != "")
                  {

                     const Dstring str = Dstring::render (format, cl);
                     Label label (str, this_p, 'c', 'c');
                     label.set_text_angle (theta);
                     label.cairo (cr);

//                     const Rect& rect = (const Rect&)label;
//                     clip_polygon_ptr->add (rect);

                     //rect.cairo (cr);
                     //cr->stroke ();

                     label_point_set.add (this_p);

                  }

               }

            }

            last_point = this_point;

         }

      }

   }

   //cr->restore ();
   return clip_polygon_ptr;

}

// Basically patching up "naughty" grid values
void
Contour::init_a (const Vector_Field_2D& vector_field_2d,
                 const Integer vector_element,
                 const Tuple& coordinate_tuple_x,
                 const Tuple& coordinate_tuple_y)
{

   const Tuple& tuple_x = coordinate_tuple_x;
   const Tuple& tuple_y = coordinate_tuple_y;
   scalar_data_2d_ptr = new Scalar_Data_2D (tuple_x, tuple_y);

   Real datum;
   min_value = GSL_POSINF;
   max_value = GSL_NEGINF;

   for (Integer i = 0; i < tuple_x.size (); i++)
   {
      const Real& x = tuple_x[i];
      for (Integer j = 0; j < tuple_y.size (); j++)
      {
         const Real& y = tuple_y[j];
         datum = vector_field_2d.evaluate (vector_element, x, y, VALUE);
         if (gsl_finite (datum) && datum < min_value) { min_value = datum; }
         if (gsl_finite (datum) && datum > max_value) { max_value = datum; }
         scalar_data_2d_ptr->set_datum (i, j, datum);
      }
   }

}

// Basically patching up "naughty" grid values
void
Contour::init_a (const Vector_Field_2D& vector_field_2d,
                 const Scalarization_2d scalarization_2d,
                 const Integer vector_element_0,
                 const Integer vector_element_1,
                 const Tuple& coordinate_tuple_x,
                 const Tuple& coordinate_tuple_y)
{

   const Vector_Field_2D& vf_2d = vector_field_2d;
   const Tuple& tuple_x = coordinate_tuple_x;
   const Tuple& tuple_y = coordinate_tuple_y;
   scalar_data_2d_ptr = new Scalar_Data_2D (tuple_x, tuple_y);

   Real datum;
   max_value = GSL_NEGINF;
   min_value = GSL_POSINF;

   for (Integer i = 0; i < tuple_x.size (); i++)
   {

      const Real& x = tuple_x[i];

      for (Integer j = 0; j < tuple_y.size (); j++)
      {

         const Real& y = tuple_y[j];

         switch (scalarization_2d)
         {
            case MAGNITUDE:
            {
               Real u = vf_2d.evaluate (vector_element_0, x, y, VALUE);
               Real v = vf_2d.evaluate (vector_element_1, x, y, VALUE);
               datum = sqrt (u*u + v*v);

               break;
            }
            case DIVERGENCE:
            {
               Real u_x = vf_2d.evaluate (vector_element_0, x, y, DX);
               Real v_y = vf_2d.evaluate (vector_element_1, x, y, DY);
               datum = u_x + v_y;
               break;
            }
            case VORTICITY:
            {
               Real v_x = vf_2d.evaluate (vector_element_1, x, y, DX);
               Real u_y = vf_2d.evaluate (vector_element_0, x, y, DY);
               datum = v_x - u_y;
               break;
            }

         }

         if (gsl_finite (datum) && datum < min_value) { min_value = datum; }
         if (gsl_finite (datum) && datum > max_value) { max_value = datum; }
         scalar_data_2d_ptr->set_datum (i, j, datum);

      }

   }

}

void
Contour::init_b (const Tuple& level_tuple,
                 const Real epsilon)
{

   valid_polygon_ptr = NULL;

   const Tuple& tuple_x = scalar_data_2d_ptr->get_coordinate_tuple (0);
   const Tuple& tuple_y = scalar_data_2d_ptr->get_coordinate_tuple (1);

   if (tuple_x.size () < 2 || tuple_y.size () < 2) { return; }

   cell_ptrs = get_cell_ptrs ();
   max_value = GSL_NEGINF;
   min_value = GSL_POSINF;

   for (Integer i = 0; i < tuple_x.size (); i++)
   {
      for (Integer j = 0; j < tuple_y.size (); j++)
      {

         Real& datum = scalar_data_2d_ptr->get_datum (i, j);

         for (Tuple::const_iterator iterator = level_tuple.begin ();
              iterator != level_tuple.end (); iterator++)
         {
            const Real level = *(iterator);
            if (fabs (datum - level) < epsilon * 0.1)
            {
               datum += epsilon;
               break;
            }
         }

      }
   }

   //valid_polygon_ptr = new Valid_Polygon (*scalar_data_2d_ptr);

   const bool nil = (level_tuple.size () == 0);
   const Real substitude = (nil ? 0 :
      level_tuple.back () +
      M_PI * (level_tuple.back () -
       level_tuple.front ()));
   substitude_nan (substitude);

   for (Integer i = 0; i < tuple_x.size (); i++)
   {
      for (Integer j = 0; j < tuple_y.size (); j++)
      {

         Real& datum = scalar_data_2d_ptr->get_datum (i, j);

         if (!gsl_isnan (datum))
         {
            if (datum > max_value) { max_value = datum; }
            if (datum < min_value) { min_value = datum; }
         }

      }
   }

   this->boundary_ptr = new Boundary (*scalar_data_2d_ptr, level_tuple);

   work_out_edges ();
   work_out_side_edges ();
   trace_isolines ();


}

void
Contour::init (const Vector_Field_2D& vector_field_2d,
               const Integer vector_element,
               const Real step,
               const Tuple& coordinate_tuple_x,
               const Tuple& coordinate_tuple_y,
               const Real epsilon)
{

   const Vector_Field_2D& vf_2d = vector_field_2d;
   init_a (vf_2d, vector_element, coordinate_tuple_x, coordinate_tuple_y);

   const Real start = min_value - modulo (min_value, step);
   const Real end = max_value - modulo (max_value, step) + step;
   const Integer n = Integer (round ((end - start) / step)) + 1;
   const Real ep = (gsl_finite (epsilon) ? epsilon : step * 1e-3);

   if (n > 0 && n < 1000)
   {
      this->level_tuple = Tuple (n, start, end);
   }
   else
   {
      this->level_tuple = Tuple ();
   }

   init_b (level_tuple, ep);

}

void
Contour::init (const Vector_Field_2D& vector_field_2d,
               const Integer vector_element,
               const Tuple& level_tuple,
               const Tuple& coordinate_tuple_x,
               const Tuple& coordinate_tuple_y,
               const Real epsilon)
{

   const Vector_Field_2D& vf_2d = vector_field_2d;
   init_a (vf_2d, vector_element, coordinate_tuple_x, coordinate_tuple_y);

   const Real ep = (gsl_finite (epsilon) ? epsilon : 
      std::max (fabs (max_value), fabs (min_value)) * 1e-5);

   init_b (level_tuple, ep);

}

void
Contour::init (const Vector_Field_2D& vector_field_2d,
               const Scalarization_2d scalarization_2d,
               const Integer vector_element_0,
               const Integer vector_element_1,
               const Real step,
               const Tuple& coordinate_tuple_x,
               const Tuple& coordinate_tuple_y,
               const Real epsilon)
{

   init_a (vector_field_2d, scalarization_2d, vector_element_0,
      vector_element_1, coordinate_tuple_x, coordinate_tuple_y);

   const Real start = min_value - modulo (min_value, step);
   const Real end = max_value - modulo (max_value, step) + step;
   const Integer n = Integer (round ((end - start) / step)) + 1;
   const Real ep = (gsl_finite (epsilon) ? epsilon : step * 1e-3);

   if (n > 0 && n < 1000)
   {
      this->level_tuple = Tuple (n, start, end);
   }
   else
   {
      this->level_tuple = Tuple ();
   }

   init_b (level_tuple, epsilon);

}

void
Contour::init (const Vector_Field_2D& vector_field_2d,
               const Scalarization_2d scalarization_2d,
               const Integer vector_element_0,
               const Integer vector_element_1,
               const Tuple& level_tuple,
               const Tuple& coordinate_tuple_x,
               const Tuple& coordinate_tuple_y,
               const Real epsilon)
{

   init_a (vector_field_2d, scalarization_2d, vector_element_0,
      vector_element_1, coordinate_tuple_x, coordinate_tuple_y);

   const Real ep = (gsl_finite (epsilon) ? epsilon : 
      std::max (fabs (max_value), fabs (min_value)) * 1e-5);

   init_b (level_tuple, ep);

}

// 1
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Tuple& level_tuple,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Tuple& tuple_x = vector_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = vector_data_2d.get_coordinate_tuple (1);

   init (vector_data_2d, vector_element, level_tuple,
      tuple_x, tuple_y, epsilon);

}

// 2
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Real step,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Tuple& tuple_x = vector_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = vector_data_2d.get_coordinate_tuple (1);

   init (vector_data_2d, vector_element, step, tuple_x, tuple_y, epsilon);

}

// 3
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Tuple& level_tuple,
                  const Domain_2D& domain_2d,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x   = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y   = domain_2d.domain_y.end;

   if (start_x > end_x) { return; }
   if (start_y > end_y) { return; }

   const Tuple& orig_tuple_x = vector_data_2d.get_coordinate_tuple (0);
   const Tuple& orig_tuple_y = vector_data_2d.get_coordinate_tuple (1);

   Tuple tuple_x;
   Tuple tuple_y;

   for (Tuple::const_iterator iterator = orig_tuple_x.begin ();
        iterator != orig_tuple_x.end (); iterator++)
   {
      const Real& x = *(iterator);
      if ((x - start_x) * (x - end_x) <= 0)
      {
         tuple_x.push_back (x);
      }
   }

   for (Tuple::const_iterator iterator = orig_tuple_y.begin ();
        iterator != orig_tuple_y.end (); iterator++)
   {
      const Real& y = *(iterator);
      if ((y - start_y) * (y - end_y) <= 0)
      {
         tuple_y.push_back (y);
      }
   }

   init (vector_data_2d, vector_element, level_tuple,
      tuple_x, tuple_y, epsilon);

}

// 4
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Real step,
                  const Domain_2D& domain_2d,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x   = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y   = domain_2d.domain_y.end;

   if (start_x > end_x) { return; }
   if (start_y > end_y) { return; }

   const Tuple& orig_tuple_x = vector_data_2d.get_coordinate_tuple (0);
   const Tuple& orig_tuple_y = vector_data_2d.get_coordinate_tuple (1);

   Tuple tuple_x;
   Tuple tuple_y;

   for (Tuple::const_iterator iterator = orig_tuple_x.begin ();
        iterator != orig_tuple_x.end (); iterator++)
   {
      const Real& x = *(iterator);
      if ((x - start_x) * (x - end_x) <= 0)
      {
         tuple_x.push_back (x);
      }
   }

   for (Tuple::const_iterator iterator = orig_tuple_y.begin ();
        iterator != orig_tuple_y.end (); iterator++)
   {
      const Real& y = *(iterator);
      if ((y - start_y) * (y - end_y) <= 0)
      {
         tuple_y.push_back (y);
      }
   }

   init (vector_data_2d, vector_element, step, tuple_x, tuple_y, epsilon);

}

// 5
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Tuple& level_tuple,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{
   init (vector_data_2d, vector_element, level_tuple,
      coordinate_tuple_x, coordinate_tuple_y, epsilon);
}

// 6
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Real step,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{
   init (vector_data_2d, vector_element, step, coordinate_tuple_x,
      coordinate_tuple_y, epsilon);
}

// 7
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Tuple& level_tuple,
                  const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x   = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y   = domain_2d.domain_y.end;

   if (start_x > end_x) { return; }
   if (start_y > end_y) { return; }

   const Tuple tuple_x = Tuple (size_2d.i, start_x, end_x);
   const Tuple tuple_y = Tuple (size_2d.j, start_y, end_y);

   init (vector_data_2d, vector_element, level_tuple,
      tuple_x, tuple_y, epsilon);

}

// 8
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Integer vector_element,
                  const Real step,
                  const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x   = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y   = domain_2d.domain_y.end;

   if (start_x > end_x) { return; }
   if (start_y > end_y) { return; }

   const Tuple tuple_x = Tuple (size_2d.i, start_x, end_x);
   const Tuple tuple_y = Tuple (size_2d.j, start_y, end_y);

   init (vector_data_2d, vector_element, step,
      tuple_x, tuple_y, epsilon);

}

// 9
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Tuple& level_tuple,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Tuple& tuple_x = vector_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = vector_data_2d.get_coordinate_tuple (1);

   init (vector_data_2d, scalarization_2d, vector_element_0,
      vector_element_1, level_tuple, tuple_x, tuple_y,
      epsilon);

}

// 10
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Real step,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Tuple& tuple_x = vector_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = vector_data_2d.get_coordinate_tuple (1);

   init (vector_data_2d, scalarization_2d, vector_element_0,
      vector_element_1, step, tuple_x, tuple_y, epsilon);

}

// 11
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Tuple& level_tuple,
                  const Domain_2D& domain_2d,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x   = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y   = domain_2d.domain_y.end;

   if (start_x > end_x) { return; }
   if (start_y > end_y) { return; }

   const Tuple& orig_tuple_x = vector_data_2d.get_coordinate_tuple (0);
   const Tuple& orig_tuple_y = vector_data_2d.get_coordinate_tuple (1);

   Tuple tuple_x;
   Tuple tuple_y;

   for (Tuple::const_iterator iterator = orig_tuple_x.begin ();
        iterator != orig_tuple_x.end (); iterator++)
   {
      const Real& x = *(iterator);
      if ((x - start_x) * (x - end_x) <= 0)
      {
         tuple_x.push_back (x);
      }
   }

   for (Tuple::const_iterator iterator = orig_tuple_y.begin ();
        iterator != orig_tuple_y.end (); iterator++)
   {
      const Real& y = *(iterator);
      if ((y - start_y) * (y - end_y) <= 0)
      {
         tuple_y.push_back (y);
      }
   }

   init (vector_data_2d, scalarization_2d, vector_element_0,
      vector_element_1, level_tuple, tuple_x, tuple_y,
      epsilon);

}

// 12
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Real step,
                  const Domain_2D& domain_2d,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x   = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y   = domain_2d.domain_y.end;

   if (start_x > end_x) { return; }
   if (start_y > end_y) { return; }

   const Tuple& orig_tuple_x = vector_data_2d.get_coordinate_tuple (0);
   const Tuple& orig_tuple_y = vector_data_2d.get_coordinate_tuple (1);

   Tuple tuple_x;
   Tuple tuple_y;

   for (Tuple::const_iterator iterator = orig_tuple_x.begin ();
        iterator != orig_tuple_x.end (); iterator++)
   {
      const Real& x = *(iterator);
      if ((x - start_x) * (x - end_x) <= 0)
      {
         tuple_x.push_back (x);
      }
   }

   for (Tuple::const_iterator iterator = orig_tuple_y.begin ();
        iterator != orig_tuple_y.end (); iterator++)
   {
      const Real& y = *(iterator);
      if ((y - start_y) * (y - end_y) <= 0)
      {
         tuple_y.push_back (y);
      }
   }

   init (vector_data_2d, scalarization_2d, vector_element_0,
      vector_element_1, step, tuple_x, tuple_y,
      epsilon);

}

// 13
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Tuple& level_tuple,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{
   init (vector_data_2d, scalarization_2d, vector_element_0,
      vector_element_1, level_tuple, coordinate_tuple_x,
      coordinate_tuple_y, epsilon);
}

// 14
Contour::Contour (const Vector_Data_2D& vector_data_2d,
                  const Scalarization_2d scalarization_2d,
                  const Integer vector_element_0,
                  const Integer vector_element_1,
                  const Real step,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{
   init (vector_data_2d, scalarization_2d, vector_element_0,
      vector_element_1, step, coordinate_tuple_x,
      coordinate_tuple_y, epsilon);
}

// 15
Contour::Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Tuple& level_tuple,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Tuple& tuple_x = scalar_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = scalar_data_2d.get_coordinate_tuple (1);

   init ((const Vector_Field_2D&)scalar_data_2d, 0, level_tuple, tuple_x, tuple_y, epsilon);

}

// 16
Contour::Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Real step,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Tuple& tuple_x = scalar_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = scalar_data_2d.get_coordinate_tuple (1);

   init (scalar_data_2d, 0, step, tuple_x, tuple_y, epsilon);

}

// 17
Contour::Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Tuple& level_tuple,
                  const Domain_2D& domain_2d,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x   = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y   = domain_2d.domain_y.end;

   if (start_x > end_x) { return; }
   if (start_y > end_y) { return; }

   const Tuple& orig_tuple_x = scalar_data_2d.get_coordinate_tuple (0);
   const Tuple& orig_tuple_y = scalar_data_2d.get_coordinate_tuple (1);

   Tuple tuple_x;
   Tuple tuple_y;

   for (Tuple::const_iterator iterator = orig_tuple_x.begin ();
        iterator != orig_tuple_x.end (); iterator++)
   {
      const Real& x = *(iterator);
      if ((x - start_x) * (x - end_x) <= 0)
      {
         tuple_x.push_back (x);
      }
   }

   for (Tuple::const_iterator iterator = orig_tuple_y.begin ();
        iterator != orig_tuple_y.end (); iterator++)
   {
      const Real& y = *(iterator);
      if ((y - start_y) * (y - end_y) <= 0)
      {
         tuple_y.push_back (y);
      }
   }

   init (scalar_data_2d, 0, level_tuple, tuple_x, tuple_y, epsilon);

}

// 18
Contour::Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Real step,
                  const Domain_2D& domain_2d,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x   = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y   = domain_2d.domain_y.end;

   if (start_x > end_x) { return; }
   if (start_y > end_y) { return; }

   const Tuple& orig_tuple_x = scalar_data_2d.get_coordinate_tuple (0);
   const Tuple& orig_tuple_y = scalar_data_2d.get_coordinate_tuple (1);

   Tuple tuple_x;
   Tuple tuple_y;

   for (Tuple::const_iterator iterator = orig_tuple_x.begin ();
        iterator != orig_tuple_x.end (); iterator++)
   {
      const Real& x = *(iterator);
      if ((x - start_x) * (x - end_x) <= 0)
      {
         tuple_x.push_back (x);
      }
   }

   for (Tuple::const_iterator iterator = orig_tuple_y.begin ();
        iterator != orig_tuple_y.end (); iterator++)
   {
      const Real& y = *(iterator);
      if ((y - start_y) * (y - end_y) <= 0)
      {
         tuple_y.push_back (y);
      }
   }

   init (scalar_data_2d, 0, step, tuple_x, tuple_y, epsilon);

}

// 19
Contour::Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Tuple& level_tuple,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{
   init (scalar_data_2d, 0, level_tuple, coordinate_tuple_x,
      coordinate_tuple_y, epsilon);
}

// 20
Contour::Contour (const Scalar_Data_2D& scalar_data_2d,
                  const Real step,
                  const Tuple& coordinate_tuple_x,
                  const Tuple& coordinate_tuple_y,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{
   init (scalar_data_2d, 0, step, coordinate_tuple_x,
      coordinate_tuple_y, epsilon);
}

// 21
Contour::Contour (const Vector_Field_2D& vector_field_2d,
                  const Integer vector_element,
                  const Tuple& level_tuple,
                  const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Real epsilon)
   : level_tuple (level_tuple),
     scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x   = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y   = domain_2d.domain_y.end;

   if (start_x > end_x) { return; }
   if (start_y > end_y) { return; }

   const Tuple tuple_x = Tuple (size_2d.i, start_x, end_x);
   const Tuple tuple_y = Tuple (size_2d.j, start_y, end_y);

   init (vector_field_2d, vector_element, level_tuple,
      tuple_x, tuple_y, epsilon);

}

// 22
Contour::Contour (const Vector_Field_2D& vector_field_2d,
                  const Integer vector_element,
                  const Real step,
                  const Size_2D& size_2d,
                  const Domain_2D& domain_2d,
                  const Real epsilon)
   : scalar_data_2d_ptr (NULL),
     cell_ptrs (NULL),
     boundary_ptr (NULL),
     isoline_ptr_lists (NULL),
     polygon_ptrs (NULL),
     valid_polygon_ptr (NULL)
{

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x   = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y   = domain_2d.domain_y.end;

   if (start_x > end_x) { return; }
   if (start_y > end_y) { return; }

   const Tuple tuple_x = Tuple (size_2d.i, start_x, end_x);
   const Tuple tuple_y = Tuple (size_2d.j, start_y, end_y);

   init (vector_field_2d, vector_element, step, tuple_x, tuple_y, epsilon);

}

Contour::~Contour ()
{

   if (valid_polygon_ptr != NULL) { delete valid_polygon_ptr; }

   if (isoline_ptr_lists != NULL)
   {

      const Integer n = level_tuple.size ();

      for (Integer li = 0; li < n; li++)
      {

         list<Isoline*>& ipl = isoline_ptr_lists[li];

         for (list<Isoline*>::iterator iterator = ipl.begin ();
              iterator != ipl.end (); iterator++)
         {
            Isoline* isoline_ptr = (Isoline*)(*(iterator));
            delete isoline_ptr;
         }

      }

      delete[] isoline_ptr_lists;

   }


   if (cell_ptrs != NULL && scalar_data_2d_ptr != NULL)
   {
      const Tuple& tuple_x = scalar_data_2d_ptr->get_coordinate_tuple (0);
      const Tuple& tuple_y = scalar_data_2d_ptr->get_coordinate_tuple (1);
      if (tuple_x.size () > 1 && tuple_y.size () > 1)
      {
         const Integer n = (tuple_x.size () - 1) * (tuple_y.size () - 1);
         for (Integer i = 0; i < n; i++)
         {
            Cell* cell_ptr = cell_ptrs[i];
            delete cell_ptr;
         }
      }
      delete[] cell_ptrs;
   }

   if (boundary_ptr != NULL) { delete boundary_ptr; }

   if (scalar_data_2d_ptr != NULL) { delete scalar_data_2d_ptr; }

   if (polygon_ptrs != NULL)
   {
      const Integer n = level_tuple.size ();
      for (Integer i = 0; i < n - 1; i++)
      {
         delete polygon_ptrs[i];
      }
      delete[] polygon_ptrs;
   }

}

const Polygon&
Contour::get_polygon (const Integer index) const
{
   return *(polygon_ptrs[index]);
}

void
Contour::work_out_polygons ()
{

   if (level_tuple.size () == 0) { return; }

   typedef Polygon* Polygon_Ptr;
   const Integer n = level_tuple.size ();
   this->polygon_ptrs = new Polygon_Ptr[n - 1];

   for (Integer i = 0; i < n - 1; i++)
   {

      polygon_ptrs[i] = new Polygon ();
      work_out_polygon (i);

//      if (valid_polygon_ptr != NULL)
//      {
//         Polygon* candidate_polygon_ptr = polygon_ptrs[i];
//         Polygon* clipped_polygon_ptr = Polygon::boolean_op (INTERSECTION,
//            *candidate_polygon_ptr, *valid_polygon_ptr, 0.011);
//         delete clipped_polygon_ptr;
//      }

   }

}

const Real&
Contour::get_max_value () const
{
   return max_value;
}

const Real&
Contour::get_min_value () const
{
   return min_value;
}

void
Contour::render (const RefPtr<Context>& cr,
                 const Transform_2D& transform,
                 const Color_Chooser& color_chooser)
{
   render_fill (cr, transform, color_chooser);
   render_isolines (cr, transform);
}

void
Contour::render_fill (const RefPtr<Context>& cr,
                      const Transform_2D& transform,
                      const Color_Chooser& color_chooser)
{

   if (scalar_data_2d_ptr == NULL) { return; }

   if (polygon_ptrs == NULL)
   {
      work_out_polygons ();
   }

   cr->save ();

   cr->set_fill_rule (FILL_RULE_EVEN_ODD);
   const Integer n = level_tuple.size ();

   for (Integer li = 0; li < n - 1; li++)
   {

      const Real cl = level_tuple[li];
      const Real ncl = level_tuple[li + 1];
      const Real l = (cl + ncl) / 2;
      const Color& color = color_chooser.get_color (l);
      color.cairo (cr);

      const Polygon& polygon = *(polygon_ptrs[li]);
      transform.cairo (cr, polygon);
      cr->fill ();

   }

   cr->restore ();

}

void
Contour::render_isoline (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Integer level_index,
                         const Dstring& format,
                         const Real label_multiplier,
                         const Real label_offset,
                         const Real label_distance,
                         const Integer label_stride) const
{

   cr->save ();
   Point_2D p;
   Polygon* clip_polygon_ptr = NULL;

   if (gsl_finite (label_multiplier))
   {

      Integer ls = label_stride;
      if (ls < 1)
      {
         const Size_2D& size_2d = scalar_data_2d_ptr->get_size_2d ();
         ls = std::min (size_2d.i, size_2d.j) * 3 / 4;
      }

      clip_polygon_ptr = render_label (cr, transform, format,
         label_multiplier, label_offset, label_distance, ls);

   }

   const list<Isoline*>& ipl = isoline_ptr_lists[level_index];

   for (list<Isoline*>::const_iterator iterator = ipl.begin ();
        iterator != ipl.end (); iterator++)
   {

      const Isoline& isoline = *(*(iterator));
      if (isoline.size () < 2) { return; }

      if (clip_polygon_ptr != NULL)
      {
         transform.cairo (cr, isoline, *clip_polygon_ptr, true);
      }
      else
      {
         transform.cairo (cr, isoline);
      }

   }

   cr->stroke ();


   if (clip_polygon_ptr != NULL)
   {
      delete clip_polygon_ptr;
   }

   cr->restore ();

}

void
Contour::render_isolines (const RefPtr<Context>& cr,
                          const Transform_2D& transform,
                          const Dstring& format,
                          const Real label_multiplier,
                          const Real label_offset,
                          const Real label_distance,
                          const Integer label_stride) const
{

   if (scalar_data_2d_ptr == NULL) { return; }

   cr->save ();
   Point_2D p;
   Polygon* clip_polygon_ptr = NULL;

   if (gsl_finite (label_multiplier))
   {

      Integer ls = label_stride;
      if (ls < 1)
      {
         const Size_2D& size_2d = scalar_data_2d_ptr->get_size_2d ();
         ls = std::min (size_2d.i, size_2d.j) * 3 / 4;
      }

      clip_polygon_ptr = render_label (cr, transform, format,
         label_multiplier, label_offset, label_distance, ls);

   }

   const Integer n = level_tuple.size ();

   for (Integer li = 0; li < n; li++)
   {

      const list<Isoline*>& ipl = isoline_ptr_lists[li];

      for (list<Isoline*>::const_iterator iterator = ipl.begin ();
           iterator != ipl.end (); iterator++)
      {

         const Isoline& isoline = *(*(iterator));
         if (isoline.size () < 2) { return; }

         if (clip_polygon_ptr != NULL)
         {
            transform.cairo (cr, isoline, *clip_polygon_ptr, true);
         }
         else
         {
            transform.cairo (cr, isoline);
         }

      }

      cr->stroke ();

   }


   if (clip_polygon_ptr != NULL)
   {
      delete clip_polygon_ptr;
   }

   cr->restore ();

}

Scalar_Raster::Scalar_Raster (const Transform_2D& transform,
                              const Scalar_Data_2D& scalar_data_2d,
                              const Box_2D& box_2d,
                              const Color_Chooser& color_chooser)
   : Raster (box_2d),
     transform (transform)
{

   const Size_2D& size_2d = box_2d.size_2d;
   const Index_2D& anchor = box_2d.index_2d;

   const Index_2D end_index (anchor.i + size_2d.i, anchor.j  + size_2d.j);

   #pragma omp parallel for
   for (Integer i = anchor.i; i < end_index.i; i++)
   {

      for (Integer j = anchor.j; j < end_index.j; j++)
      {

         const Point_2D& p = transform.reverse (Real (i), Real (j));
         if (scalar_data_2d.out_of_bounds (p.x, p.y)) { continue; }

         try
         {
            const Real value = scalar_data_2d.evaluate (p.x, p.y, VALUE);
            const Color color = color_chooser.get_color (value);
            set_pixel (i - anchor.i, j - anchor.j, color);
         }
         catch (const Exception& e) { cerr << e << endl; }

      }
   }

}

void
Scalar_Renderer::render (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Vector_Data_2D& vector_data_2d,
                         const Integer vector_element,
                         const Domain_2D& domain,
                         const Color_Chooser& color_chooser)
{

   const Grid_nD& grid_nd = (const Grid_nD&)vector_data_2d;
   const Tuple& tuple_x = vector_data_2d.get_coordinate_tuple (0);
   const Tuple& tuple_y = vector_data_2d.get_coordinate_tuple (1);

   cr->save ();
   cr->set_antialias (ANTIALIAS_NONE);
   const Integer ve = vector_element;

   for (Integer i = 0; i < tuple_x.size (); i++)
   {

      const Real x = tuple_x[i];
      Real xb = x;
      Real xe = x;

      if (i != 0) { xb += tuple_x [i - 1]; xb /= 2; }
      if (i != tuple_x.size () - 1) { xe += tuple_x [i + 1]; xe /= 2; }

      for (Integer j = 0; j < tuple_y.size (); j++)
      {

         const Real y = tuple_y[j];
         Real yb = y;
         Real ye = y;

         if (domain.is_out_of_bounds (x, y)) { continue; }

         if (j != 0) { yb += tuple_y[j - 1]; yb /= 2; }
         if (j != tuple_y.size () - 1) { ye += tuple_y [j + 1]; ye /= 2; }

         const Real value = vector_data_2d.get_datum (ve, i, j);
         if (gsl_isnan (value)) { continue; }

         const Color color = color_chooser.get_color (value);
         if (color.is_nac ()) { continue; }

         color.cairo (cr);

         const Point_2D bb = transform.transform (xb, yb);
         const Point_2D be = transform.transform (xb, ye);
         const Point_2D eb = transform.transform (xe, yb);
         const Point_2D ee = transform.transform (xe, ye);

         Polygon polygon;
         polygon.add (bb);
         polygon.add (be);
         polygon.add (ee);
         polygon.add (eb);

         polygon.cairo (cr);
         cr->fill ();

      }

   }

   cr->restore ();

}

void
Scalar_Renderer::render (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Vector_Data_2D& vector_data_2d,
                         const Integer vector_element,
                         const Color_Chooser& color_chooser)
{
   const Domain_2D& domain = vector_data_2d.get_domain_2d ();
   render (cr, transform, vector_data_2d, vector_element,
      domain, color_chooser);
}

void
Scalar_Renderer::render (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Scalar_Data_2D& scalar_data_2d,
                         const Domain_2D& domain,
                         const Color_Chooser& color_chooser)
{
   const Vector_Data_2D& vector_data_2d = (const Vector_Data_2D&)scalar_data_2d;
   render (cr, transform, vector_data_2d, 0, domain, color_chooser);
}

void
Scalar_Renderer::render (const RefPtr<Context>& cr,
                         const Transform_2D& transform,
                         const Scalar_Data_2D& scalar_data_2d,
                         const Color_Chooser& color_chooser)
{
   const Vector_Data_2D& vector_data_2d = (const Vector_Data_2D&)scalar_data_2d;
   render (cr, transform, vector_data_2d, 0, color_chooser);
}

