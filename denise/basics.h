
// basics.h
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

#ifndef DENISE_BASICS_H
#define DENISE_BASICS_H

#if HAVE_IEEEFP_H
#include <ieeefp.h>
#endif

#include <complex>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <gsl/gsl_math.h>
#include <cairomm/context.h>

using namespace std;

namespace denise
{

   enum Extremity
   {
      MIN,
      MAX
   };

   enum Evaluate_Op
   {
      VALUE,
      DZ,
      DX,
      DY,
      DZ2,
      DX2,
      DY2,
      DZY,
      DZX,
      DXY,
      DZXY,
      LAPLACIAN
   };

   enum Coordinate_System
   {
      CARTESIAN,
      LAT_LONG,
      CYLINDRICAL,
      SPHERICAL
   };

   typedef int
   Integer;

   typedef double
   Real;

   const Real
   M_2_TIMES_PI = M_PI * 2;

   const Real
   M_4_TIMES_PI = M_PI * 4;

   const Real
   M_PI_3 = M_PI / 3;

   const Real
   M_SQRT1_3 = 1 / (M_SQRT3);

   const Real
   DEGREE_TO_RADIAN = M_PI/180;

   const Real
   RADIAN_TO_DEGREE = 180/M_PI;

   class Polygon;

   class Index_nD
   {

      public:

         const Integer
         n;

         Integer*
         buffer;

         Index_nD (const Integer n);

         Index_nD (const Index_nD& index_nd);

         ~Index_nD ();

         void
         operator= (const Index_nD& index_nd);

   };

   typedef Index_nD
   Size_nD;

   class Point_nD
   {

      public:

         const Integer
         n;

         Real*
         buffer;

         Point_nD (const Integer n);

         Point_nD (const Point_nD& point_nd);

         ~Point_nD ();

         void
         operator= (const Point_nD& point_nd);

   };

   class Box_nD
   {

      public:

         Index_nD
         index_nd;

         Size_nD
         size_nd;

         Box_nD (const Size_nD& size_nd);

         Box_nD (const Index_nD& index_nd,
                 const Size_nD& size_nd);

   };

   typedef Integer
   Index_1D;

   typedef Real
   Point_1D;

   typedef Real
   Scalar;

   typedef complex<Real>
   Cmplx;

   /// Represents a 2D index
   ///
   /// Index_2D consists of 2 Integer in i and j
   class Index_2D
   {

      public:

         Integer
         i;

         Integer
         j;

         Index_2D (const Integer i = 0,
                   const Integer j = 0);

         Index_2D (const string& str);

         string
         get_string () const;

         bool
         operator == (const Index_2D& index) const;

         bool
         operator != (const Index_2D& index) const;

         Index_2D
         operator = (const Index_2D& index);

         Index_2D
         operator + (const Index_2D& index) const;

         void
         operator += (const Index_2D& index);

         Index_2D
         operator - () const;

         Index_2D
         operator - (const Index_2D& index) const;

         void
         operator -= (const Index_2D& index);

         Index_2D
         operator * (Integer index) const;

         void
         operator *= (Integer index);

         bool
         operator < (const Index_2D& index_2d) const;

         bool
         operator > (const Index_2D& index_2d) const;

   };

   /// Represents a 3D index
   ///
   /// Index_3D consists of 3 Integer in k, i and j.
   /// It should be noted that the implicit order of the components is
   /// k-i-j is rather than the i-j-k.  This is, however, more intuitive
   /// because we usually consider a stack of vertical surfaces, and
   /// according to the C convention, the outmost index comes first.
   /// So anything associated with the vertical direction comes before
   /// the horizontal counterparts.
   class Index_3D
   {

      public:

         Integer
         k;

         Integer
         i;

         Integer
         j;

         Index_3D (const Integer k = 0,
                   const Integer i = 0,
                   const Integer j = 0);

         bool
         operator == (const Index_3D index);

         bool
         operator != (const Index_3D index);

         Index_3D
         operator = (const Index_3D index);

         Index_3D
         operator + (const Index_3D index);

         void
         operator += (const Index_3D index);

         Index_3D
         operator - ();

         Index_3D
         operator - (const Index_3D index);

         void
         operator -= (const Index_3D index);

         Index_3D
         operator * (const Integer index);

         void
         operator *= (const Integer index);

   };

   /// Represents a 1-D size
   ///
   /// Size_1D is typedef'ed from Index_1D
   typedef Index_1D
   Size_1D;

   /// Represents a 2-D size
   ///
   /// Size_2D is typedef'ed from Index_2D
   typedef Index_2D
   Size_2D;

   /// Represents a 3-D size
   ///
   /// Size_2D is typedef'ed from Index_3D
   typedef Index_3D
   Size_3D;

   class Box_2D
   {

      public:

         Index_2D
         index_2d;

         Size_2D
         size_2d;

         Box_2D ();

         Box_2D (const Box_2D& box_2d);

         Box_2D (const Size_2D& size_2d);

         Box_2D (const Box_2D& box_2d,
                 const Real margin);

         Box_2D (const Box_2D& box_2d,
                 const Real margin_h,
                 const Real margin_v);

         Box_2D (const Box_2D& box_2d,
                 const Real margin_l,
                 const Real margin_r,
                 const Real margin_t,
                 const Real margin_b);

         Box_2D (const Index_2D& index_2d,
                 const Size_2D& size_2d,
                 const Real anchor_x = 0,
                 const Real anchor_y = 0);

         Box_2D (const Index_2D& index_2d,
                 const Size_2D& size_2d,
                 const char justify_h,
                 const char justify_v);

         Real
         get_aspect () const;

         void
         shrink (const Real margin);

         void
         shrink (const Real margin_h,
                 const Real margin_v);

         void
         shrink (const Real margin_l,
                 const Real margin_r,
                 const Real margin_t,
                 const Real margin_b);

         Box_2D
         get_shrunk_box (const Real margin) const;

         Box_2D
         get_shrunk_box (const Real margin_h,
                         const Real margin_v) const;

         Box_2D
         get_shrunk_box (const Real margin_l,
                         const Real margin_r,
                         const Real margin_t,
                         const Real margin_b) const;

         bool
         contains (const Integer i,
                   const Integer j) const;

         bool
         contains (const Real x,
                   const Real y) const;

         Index_2D
         get_nw () const;

         Index_2D
         get_ne () const;

         Index_2D
         get_sw () const;

         Index_2D
         get_se () const;

   };

   class Box_3D
   {

      public:

         Index_3D
         index_3d;

         Size_3D
         size_3d;

         Box_3D (const Size_3D& size_3d);

         Box_3D (const Index_3D& index_3d,
                 const Size_3D& size_3d);

   };

   class Point_2D
   {

      public:

         Real
         x;

         Real
         y;

         Point_2D (const Real x = GSL_NAN,
                   const Real y = GSL_NAN);

         Point_2D (const complex<Real>& c);

         Point_2D (const Index_2D& index_2d);

         Point_2D (const string& str);

         complex<Real>
         get_complex () const;

         bool
         operator == (const Point_2D& point_2d) const;

         bool
         operator != (const Point_2D& point_2d) const;

         Point_2D
         operator = (const Point_2D& point_2d);

         Point_2D
         operator + (const Point_2D& point_2d) const;

         void
         operator += (const Point_2D& point_2d);

         Point_2D
         operator - () const;

         Point_2D
         operator - (const Point_2D& point_2d) const;

         void
         operator -= (const Point_2D& point_2d);

         Point_2D
         operator * (const Real scalar) const;

         void
         operator *= (const Real scalar);

         Point_2D
         operator / (const Real scalar) const;

         void
         operator /= (const Real scalar);

         bool
         operator < (const Point_2D& point_2d) const;

         bool
         operator > (const Point_2D& point_2d) const;

         void
         perturb (const Real magnitude = 1e-3);

         bool
         is_nap () const;

   };

   /// Represents a point on the \f$\mathbf{R}^3\f$ space.
   ///
   /// Point_3D consists of 3 Real in \f$(z, x, y)\f$.
   ///
   /// It should be noted that the implicit order of the components is
   /// z-x-y is rather than the x-y-z.  This is, however, more intuitive
   /// because we usually consider a stack of vertical surfaces, and
   /// according to the C convention, the outmost index comes first.
   /// So anything associated with the vertical direction comes before
   /// the horizontal counterparts.
   class Point_3D
   {

      public:

         Real
         z;

         Real
         x;

         Real
         y;

         Point_3D (const Real z = 0,
                   const Real x = 0,
                   const Real y = 0);

         bool
         operator == (const Point_3D& point_3d) const;

         bool
         operator != (const Point_3D& point_3d) const;

         Point_3D
         operator = (const Point_3D& point_3d);

         Point_3D
         operator + (const Point_3D& point_3d) const;

         void
         operator += (const Point_3D& point_3d);

         Point_3D
         operator - () const;

         Point_3D
         operator - (const Point_3D& point_3d) const;

         void
         operator -= (const Point_3D& point_3d);

         bool
         is_nap () const;

   };

   /// Extends an STL vector of Real
   ///
   class Tuple : public vector<Real>
   {

      public:

         Tuple () { }

         Tuple (const Tuple& tuple);

         Tuple (const set<Real>& set_real);

         Tuple (const set<Real>::const_iterator begin,
                const set<Real>::const_iterator end);

         Tuple (const string& tuple_string,
                const string& delimiter = ":");

         Tuple (const Integer number_of_values,
                const Real value);

         Tuple (const Integer number_of_values,
                const Real start_value,
                const Real end_value);

         void
         add_content (const set<Real>& set_real,
                      const bool clear_first = false);

         void
         add_content (const set<Real>::const_iterator begin,
                      const set<Real>::const_iterator end,
                      const bool clear_first = false);

         void
         add_content (const string& tuple_string,
                      const string& delimiter = ":",
                      const bool clear_first = false);

         void
         add_content (const Integer number_of_values,
                      const Real value,
                      const bool clear_first = false);

         void
         add_content (const Integer number_of_values,
                      const Real start_value,
                      const Real end_value,
                      const bool clear_first = false);

         void
         pop_front ();

         void
         push_front (const Real value);

         void
         cairo_dash (const Cairo::RefPtr<Cairo::Context> cr,
                     const Real offset) const;

         void
         operator *= (const Real value);

   };

   /// Extends an STL vector of Integer
   ///
   class Ituple : public vector<Integer>
   {

      public:

         Ituple () { }

         Ituple (const set<Integer>& set_integer);

         Ituple (const set<Integer>::const_iterator begin,
                 const set<Integer>::const_iterator end);

         Ituple (const string& ituple_string,
                 const string& delimiter = ":");

         Ituple (const Integer number_of_values,
                 const Integer value);

         Ituple (const Integer number_of_values,
                 const Integer start_value,
                 const Integer step_value);

         void
         add_content (const set<Integer>& set_integer,
                      const bool clear_first = false);

         void
         add_content (const set<Integer>::const_iterator begin,
                      const set<Integer>::const_iterator end,
                      const bool clear_first = false);

         void
         add_content (const string& ituple_string,
                      const string& delimiter = ":",
                      const bool clear_first = false);

         void
         add_content (const Integer number_of_values,
                      const Integer value,
                      const bool clear_first = false);

         void
         add_content (const Integer number_of_values,
                      const Integer start_value,
                      const Integer step_value,
                      const bool clear_first = false);

         void
         pop_front ();

         void
         push_front (const Integer value);

   };

   typedef pair<Real, Real>
   Duple;

   typedef pair<Integer, Integer>
   Iduple;

   typedef vector<Duple>
   Duple_Vector;

   typedef vector<Iduple>
   Iduple_Vector;

   typedef vector<Tuple>
   Tuple_Vector;

   typedef vector<Tuple>
   Ituple_Vector;

   class Real_Profile : public map<Real, Real>
   {

      public:

         Tuple
         get_abscissa_tuple () const;

         Tuple
         get_ordinate_tuple () const;

   };

   class Domain_1D
   {

      public:

         Real
         start;

         Real
         end;

         Domain_1D (const Real start = 0,
                    const Real end = 1);

         Domain_1D (const Domain_1D& domain_1d);

         bool
         is_reverse () const;

         void
         swap ();

         void
         swap_if_reverse ();

         Domain_1D
         get_swapped () const;

         Real
         get_span () const;

         bool
         snap (Real& x) const;

         bool
         is_out_of_bounds (const Real x) const;
         
         void
         translate (const Real delta);

         Real
         normalize (const Real x,
                    const Real gamma = 1) const;

   };

   class Domain_2D
   {  

      public:
         
         Domain_1D
         domain_x;
         
         Domain_1D         
         domain_y;         

         Domain_2D (const string& str);

         Domain_2D (const Real start_x = 0,
                    const Real end_x = 1,
                    const Real start_y = 0,
                    const Real end_y = 1);
                          
         Domain_2D (const Domain_1D& domain_x,
                    const Domain_1D& domain_y);
         
         Domain_2D (const Domain_2D& domain_2d);
         
         Real
         get_width () const;
         
         Real
         get_height () const;
         
         bool
         is_out_of_bounds (const Real x,
                           const Real y) const;
         
         bool
         is_out_of_bounds (const Point_2D& point) const;
         
         void
         translate (const Real delta_x,
                    const Real delta_y);

         void
         clip_line (Point_2D& point_0,
                    Point_2D& point_1) const;

         Point_2D
         get_random_point () const;

         Real
         get_span_x () const;

         Real
         get_span_y () const;

   };

   class Domain_3D
   {

      public:

         Domain_1D
         domain_z;

         Domain_1D
         domain_x;

         Domain_1D
         domain_y;

         Domain_3D (const Real start_z,
                    const Real end_z,
                    const Real start_x,
                    const Real end_x,
                    const Real start_y,
                    const Real end_y);

         Domain_3D (const Domain_1D& domain_z,
                    const Domain_1D& domain_x,
                    const Domain_1D& domain_y);

         Domain_3D (const Domain_3D& domain_3d);

         Real
         get_depth () const;

         Real
         get_width () const;

         Real
         get_height () const;

         bool
         is_out_of_bounds (const Real z,
                           const Real x,
                           const Real y) const;
         
         bool
         is_out_of_bounds (const Point_3D& point) const;
         
         void
         translate (const Real delta_z,
                    const Real delta_x,
                    const Real delta_y);

   };

   ostream&
   operator << (ostream& out_file,
                const Index_nD& index_nd);

   ostream&
   operator << (ostream& out_file,
                const Point_nD& point_nd);

   ostream&
   operator << (ostream &out_file,
                const Index_2D& index);

   ostream&
   operator << (ostream &out_file,
                const Index_3D& index);
 
   ostream&
   operator << (ostream &out_file,
                const Box_2D& box_2d);
 
   ostream&
   operator << (ostream &out_file,
                const Point_2D& point);

   ostream&
   operator << (ostream &out_file,
                const Point_3D& point);

   ostream&
   operator << (ostream &out_file,
                const Tuple& tuple);
   ostream&
   operator << (ostream &out_file,
                const Ituple& ituple);

   ostream&
   operator << (ostream &out_file,
                const Ituple& ituple);

   ostream&
   operator << (ostream &out_file,
                const Domain_1D& domain);

}

#endif /* DENISE_BASICS_H */

