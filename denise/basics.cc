//
// basics.cc
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

#include "analysis.h"
#include "basics.h"
#include "dstring.h"

using namespace denise;
using namespace Cairo;

Index_nD::Index_nD (const Integer n)
               : n (n)
{

   buffer = new Integer[n];

   for (Integer i = 0; i < n; i++)
   {
      buffer[i] = 0;
   }

}

Index_nD::Index_nD (const Index_nD& index_nd)
               : n (index_nd.n)
{
   buffer = new Integer[n];
   *this = index_nd;
}

Index_nD::~Index_nD ()
{
   delete[] buffer;
}

void
Index_nD::operator= (const Index_nD& index_nd)
{
   for (Integer i = 0; i < n; i++)
   {
      buffer[i] = index_nd.buffer[i];
   }
}

Point_nD::Point_nD (const Integer n)
               : n (n)
{

   buffer = new Real[n];

   for (Integer i = 0; i < n; i++)
   {
      buffer[i] = 0;
   }

}

Point_nD::Point_nD (const Point_nD& point_nd)
               : n (point_nd.n)
{
   buffer = new Real[n];
   *this = point_nd;
}

Point_nD::~Point_nD ()
{
   delete[] buffer;
}

void
Point_nD::operator= (const Point_nD& point_nd)
{
   for (Integer i = 0; i < n; i++)
   {
      buffer[i] = point_nd.buffer[i];
   }
}

Box_nD::Box_nD (const Size_nD& size_nd)
    : index_nd (size_nd.n),
       size_nd (size_nd)
{
}

Box_nD::Box_nD (const Index_nD& index_nd,
                const Size_nD& size_nd)
    : index_nd (index_nd),
       size_nd (size_nd)
{
}

Index_2D::Index_2D (const Index_1D i,
                    const Index_1D j)
               : i (i),
                 j (j)
{
}

Index_2D::Index_2D (const string& str)
{
   const Tokens tokens (str, "x");
   this->i = atoi (tokens[0].c_str ());
   this->j = (tokens.size () > 1 ? atoi (tokens[1].c_str ()) : 0); 
}

string
Index_2D::get_string () const
{
   return string_render ("(%d, %d)", i, j);
}

bool
Index_2D::operator == (const Index_2D& index_2d) const
{
   return (i == index_2d.i) && (j == index_2d.j);
}

bool
Index_2D::operator != (const Index_2D& index_2d) const
{
   return (i != index_2d.i) || (j != index_2d.j);
}

Index_2D
Index_2D::operator = (const Index_2D& index_2d)
{
   i = index_2d.i;
   j = index_2d.j;

   return *this;
}

Index_2D
Index_2D::operator + (const Index_2D& index_2d) const
{
   return Index_2D (i + index_2d.i, j + index_2d.j);
}

void
Index_2D::operator += (const Index_2D& index_2d)
{
   i += index_2d.i;
   j += index_2d.j;
}

Index_2D
Index_2D::operator - () const
{
   return Index_2D (-i, -j);
}

Index_2D
Index_2D::operator - (const Index_2D& index_2d) const
{
   return Index_2D (i - index_2d.i, j - index_2d.j);
}

void
Index_2D::operator -= (const Index_2D& index_2d)
{
   i -= index_2d.i;
   j -= index_2d.j;
}

Index_2D
Index_2D::operator * (const Index_1D index_1d) const
{
   return Index_2D (index_1d * i, index_1d * j);
}

void
Index_2D::operator *= (const Index_1D index_1d)
{
   i *= index_1d;
   j *= index_1d;
}

bool
Index_2D::operator < (const Index_2D& index_2d) const
{  
   if (i == index_2d.i) { return j < index_2d.j; }
   else { return i < index_2d.i; }
}

bool
Index_2D::operator > (const Index_2D& index_2d) const
{  
   if (i == index_2d.i) { return j > index_2d.j; }
   else { return i > index_2d.i; }
}

Index_3D::Index_3D (const Index_1D k,
                    const Index_1D i,
                    const Index_1D j)
               : k (k),
                 i (i),
                 j (j)
{
}

bool
Index_3D::operator == (const Index_3D index_3d)
{
   return (k == index_3d.k) && (i == index_3d.i) && (j == index_3d.j);
}

bool
Index_3D::operator != (const Index_3D index_3d)
{
   return (k != index_3d.k) || (i != index_3d.i) || (j != index_3d.j);
}

Index_3D
Index_3D::operator = (const Index_3D index_3d)
{
   k = index_3d.k;
   i = index_3d.i;
   j = index_3d.j;

   return *this;
}

Index_3D
Index_3D::operator + (const Index_3D index_3d)
{
   return Index_3D (k + index_3d.k, i + index_3d.i, j + index_3d.j);
}

void
Index_3D::operator += (const Index_3D index_3d)
{
   k += index_3d.k;
   i += index_3d.i;
   j += index_3d.j;
}

Index_3D
Index_3D::operator - ()
{
   return Index_3D (-k, -i, -j);
}

Index_3D
Index_3D::operator - (const Index_3D index_3d)
{
   return Index_3D (k - index_3d.k, i - index_3d.i, j - index_3d.j);
}

void
Index_3D::operator -= (const Index_3D index_3d)
{
   k -= index_3d.k;
   i -= index_3d.i;
   j -= index_3d.j;
}

Index_3D
Index_3D::operator * (const Index_1D index_1d)
{
   return Index_3D (index_1d * k, index_1d * i, index_1d * j);
}

void
Index_3D::operator *= (const Index_1D index_1d)
{
   k *= index_1d;
   i *= index_1d;
   j *= index_1d;
}

Box_2D::Box_2D ()
   : index_2d (0, 0),
     size_2d (0, 0)
{
}

Box_2D::Box_2D (const Box_2D& box_2d)
   : index_2d (box_2d.index_2d),
     size_2d (box_2d.size_2d)
{
}

Box_2D::Box_2D (const Size_2D& size_2d)
   : index_2d (0, 0),
     size_2d (size_2d)
{
}

Box_2D::Box_2D (const Box_2D& box_2d,
                const Real margin)
   : index_2d (box_2d.index_2d),
     size_2d (box_2d.size_2d)
{
   shrink (margin);
}

Box_2D::Box_2D (const Box_2D& box_2d,
                const Real margin_h,
                const Real margin_v)
   : index_2d (box_2d.index_2d),
     size_2d (box_2d.size_2d)
{
   shrink (margin_h, margin_v);
}

Box_2D::Box_2D (const Box_2D& box_2d,
                const Real margin_l,
                const Real margin_r,
                const Real margin_t,
                const Real margin_b)
   : index_2d (box_2d.index_2d),
     size_2d (box_2d.size_2d)
{
   shrink (margin_l, margin_r, margin_t, margin_b);
}

Box_2D::Box_2D (const Index_2D& index_2d,
                const Size_2D& size_2d,
                const Real anchor_x,
                const Real anchor_y)
   : index_2d (index_2d),
     size_2d (size_2d)
{
   this->index_2d.i -= Size_1D (rint (size_2d.i * anchor_x));
   this->index_2d.j -= Size_1D (rint (size_2d.j * anchor_y));
}

Box_2D::Box_2D (const Index_2D& index_2d,
                const Size_2D& size_2d,
                const char justify_h,
                const char justify_v)
    : index_2d (index_2d),
       size_2d (size_2d)
{

   if (justify_h == 'c') { this->index_2d.i -= size_2d.i/2; }
   else if (justify_h == 'r') { this->index_2d.i -= size_2d.i; }

   if (justify_v == 'c') { this->index_2d.j -= size_2d.j/2; }
   else if (justify_v == 't') { this->index_2d.j -= size_2d.j; }

}

Real
Box_2D::get_aspect () const
{
   return Real (size_2d.i) / Real (size_2d.j);
}

void
Box_2D::shrink (const Real margin)
{
   shrink (margin, margin, margin, margin);
}

void
Box_2D::shrink (const Real margin_h,
                const Real margin_v)
{
   shrink (margin_h, margin_h, margin_v, margin_v);
}

void
Box_2D::shrink (const Real margin_l,
                const Real margin_r,
                const Real margin_t,
                const Real margin_b)
{
   index_2d.i += Integer (margin_l);
   index_2d.j += Integer (margin_t);
   size_2d.i -= Integer (margin_l + margin_r);
   size_2d.j -= Integer (margin_t + margin_b);
}

Box_2D
Box_2D::get_shrunk_box (const Real margin) const
{
   Box_2D shrunk_box (*this);
   shrunk_box.shrink (margin);
   return shrunk_box;
}

Box_2D
Box_2D::get_shrunk_box (const Real margin_h,
                        const Real margin_v) const
{
   Box_2D shrunk_box (*this);
   shrunk_box.shrink (margin_h, margin_v);
   return shrunk_box;
}

Box_2D
Box_2D::get_shrunk_box (const Real margin_l,
                        const Real margin_r,
                        const Real margin_t,
                        const Real margin_b) const
{
   Box_2D shrunk_box (*this);
   shrunk_box.shrink (margin_l, margin_r, margin_t, margin_b);
   return shrunk_box;
}

bool
Box_2D::contains (const Integer i,
                  const Integer j) const
{
   return (i >= index_2d.i && i <= index_2d.i + size_2d.i &&
           j >= index_2d.j && j <= index_2d.j + size_2d.j);
}

bool
Box_2D::contains (const Real x,
                  const Real y) const
{
   return contains (Integer (rint (x)), Integer (rint (y)));
}

Box_3D::Box_3D (const Size_3D& size_3d)
    : index_3d (0, 0, 0),
       size_3d (size_3d)
{
}

Box_3D::Box_3D (const Index_3D& index_3d,
                const Size_3D& size_3d)
    : index_3d (index_3d),
       size_3d (size_3d)
{
}

Point_2D::Point_2D (const Real x,
                    const Real y)
   : x (x),
     y (y)
{
}

Point_2D::Point_2D (const complex<Real>& c)
   : x (std::real (c)),
     y (std::imag (c))
{
}

complex<Real>
Point_2D::get_complex () const
{
   return complex<Real> (x, y);
}

bool
Point_2D::operator == (const Point_2D& point_2d) const
{
   return (x == point_2d.x) && (y == point_2d.y);
}

bool
Point_2D::operator != (const Point_2D& point_2d) const
{
   return (x != point_2d.x) || (y != point_2d.y);
}

Point_2D
Point_2D::operator = (const Point_2D& point_2d)
{
   x = point_2d.x;
   y = point_2d.y;

   return *this;
}

Point_2D
Point_2D::operator + (const Point_2D& point_2d) const
{
   return Point_2D (x + point_2d.x, y + point_2d.y);
}

void
Point_2D::operator += (const Point_2D& point_2d)
{
   x += point_2d.x;
   y += point_2d.y;
}

Point_2D
Point_2D::operator - () const
{
   return Point_2D (-x, -y);
}

Point_2D
Point_2D::operator - (const Point_2D& point_2d) const
{
   return Point_2D (x - point_2d.x, y - point_2d.y);
}

void
Point_2D::operator -= (const Point_2D& point_2d)
{
   x -= point_2d.x;
   y -= point_2d.y;
}

Point_2D
Point_2D::operator * (const Real scalar) const
{
   return Point_2D (x * scalar, y * scalar);
}

void
Point_2D::operator *= (const Real scalar)
{
   x *= scalar;
   y *= scalar;
}

Point_2D
Point_2D::operator / (const Real scalar) const
{
   return Point_2D (x / scalar, y / scalar);
}

void
Point_2D::operator /= (const Real scalar)
{
   x /= scalar;
   y /= scalar;
}

bool
Point_2D::operator < (const Point_2D& point_2d) const
{  
   if (x == point_2d.x) { return y < point_2d.y; }
   else { return x < point_2d.x; }
}

bool
Point_2D::operator > (const Point_2D& point_2d) const
{  
   if (x == point_2d.x) { return y > point_2d.y; }
   else { return x > point_2d.x; }
}

void
Point_2D::perturb (const Real magnitude)
{
   x += magnitude * ((rand () % 3) - 1);
   y += magnitude * ((rand () % 3) - 1);
}

bool
Point_2D::is_nap () const
{
   return gsl_isnan (x) || gsl_isnan (y);
}

Point_3D::Point_3D (const Real z,
                    const Real x,
                    const Real y)
               : z (z),
                 x (x),
                 y (y)
{
}

bool
Point_3D::operator == (const Point_3D& point_3d) const
{
   return (z == point_3d.z) && (x == point_3d.x) && (y == point_3d.y);
}

bool
Point_3D::operator != (const Point_3D& point_3d) const
{
   return (z != point_3d.z) || (x != point_3d.x) || (y != point_3d.y);
}

Point_3D
Point_3D::operator = (const Point_3D& point_3d)
{
   z = point_3d.z;
   x = point_3d.x;
   y = point_3d.y;

   return *this;
}

Point_3D
Point_3D::operator + (const Point_3D& point_3d) const
{
   return Point_3D (z + point_3d.z, x + point_3d.x, y + point_3d.y);
}

void
Point_3D::operator += (const Point_3D& point_3d)
{
   z += point_3d.z;
   x += point_3d.x;
   y += point_3d.y;
}

Point_3D
Point_3D::operator - () const
{
   return Point_3D (-z, -x, -y);
}

Point_3D
Point_3D::operator - (const Point_3D& point_3d) const
{
   return Point_3D (z - point_3d.z, x - point_3d.x, y - point_3d.y);
}

void
Point_3D::operator -= (const Point_3D& point_3d)
{
   z -= point_3d.z;
   x -= point_3d.x;
   y -= point_3d.y;
}

bool
Point_3D::is_nap () const
{
   return gsl_isnan(z) || gsl_isnan (x) || gsl_isnan (y);
}

Tuple::Tuple (const Tuple& tuple)
   : vector<Real> (tuple)
{
}

Tuple::Tuple (const set<Real>& set_real)
{
   add_content (set_real.begin (), set_real.end ());
}

Tuple::Tuple (const set<Real>::const_iterator begin,
              const set<Real>::const_iterator end)
{
   add_content (begin, end);
}

Tuple::Tuple (const string& tuple_string,
              const string& delimiter)
{
   add_content (tuple_string, delimiter);
}

Tuple::Tuple (const Integer number_of_values,
              const Real value)
{
   add_content (number_of_values, value);
}

Tuple::Tuple (const Integer number_of_values,
              const Real start_value,
              const Real end_value)
{
   add_content (number_of_values, start_value, end_value);
}

void
Tuple::add_content (const set<Real>& set_real,
                    const bool clear_first)
{
   add_content (set_real.begin (), set_real.end (), clear_first);
}

void
Tuple::add_content (const set<Real>::const_iterator begin,
                    const set<Real>::const_iterator end,
                    const bool clear_first)
{

   if (clear_first) { clear (); }

   for (set<Real>::const_iterator iterator = begin;
        iterator != end; iterator++)
   {
      const Real element = *(iterator);
      push_back (element);
   }

}

void
Tuple::add_content (const string& tuple_string,
                    const string& delimiter,
                    const bool clear_first)
{

   if (clear_first) { clear (); }
   const vector<string>& token_vector = tokenize (tuple_string, delimiter);

   for (vector<string>::const_iterator iterator = token_vector.begin ();
        iterator != token_vector.end (); iterator++)
   {
      const string& token = *(iterator);
      push_back (atof (token.c_str ()));
   }

}

void
Tuple::add_content (const Integer number_of_values,
                    const Real value,
                    const bool clear_first)
{
   if (clear_first) { clear (); }
   for (Integer i = 0; i < number_of_values; i++) { push_back (value); }
}

void
Tuple::add_content (const Integer number_of_values,
                    const Real start_value,
                    const Real end_value,
                    const bool clear_first)
{

   if (clear_first) { clear (); }

   Real delta_value = (end_value - start_value) / (number_of_values-1);

   for (Integer i = 0; i < number_of_values - 1; i++)
   {
      push_back (start_value + i * delta_value);
   }

   push_back (end_value);

}

void
Tuple::pop_front ()
{
   erase (begin ());
}
                
void
Tuple::push_front (const Real value)
{
   insert (begin (), value);
}

void
Tuple::cairo_dash (const RefPtr<Context> cr,
                   const Real offset) const
{
   std::vector<double> vector;
   for (Tuple::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Real& datum = *(iterator);
      vector.push_back (datum);
   }
   cr->set_dash (vector, offset);
}

void
Tuple::operator *= (const Real value)
{
   for (Tuple::iterator iterator = begin (); iterator != end (); iterator++)
   {
      Real& datum = *(iterator);
      datum *= value;
   }
}

Ituple::Ituple (const set<Integer>& set_integer)
{
   add_content (set_integer.begin (), set_integer.end ());
}

Ituple::Ituple (const set<Integer>::const_iterator begin,
                const set<Integer>::const_iterator end)
{
   add_content (begin, end);
}

Ituple::Ituple (const string& ituple_string,
                const string& delimiter)
{
   add_content (ituple_string, delimiter);
}

Ituple::Ituple (const Integer number_of_values,
                const Integer value)
{
   add_content (number_of_values, value);
}

Ituple::Ituple (const Integer number_of_values,
                const Integer start_value,
                const Integer step_value)
{
   add_content (number_of_values, start_value, step_value);
}

void
Ituple::add_content (const set<Integer>& set_integer,
                     const bool clear_first)
{
   add_content (set_integer.begin (), set_integer.end (), clear_first);
}

void
Ituple::add_content (const set<Integer>::const_iterator begin,
                     const set<Integer>::const_iterator end,
                     const bool clear_first)
{

   if (clear_first) { clear (); }

   for (set<Integer>::const_iterator iterator = begin;
        iterator != end; iterator++)
   {
      const Integer element = *(iterator);
      push_back (element);
   }

}

void
Ituple::add_content (const string& ituple_string,
                     const string& delimiter,
                     const bool clear_first)
{

   if (clear_first) { clear (); }
   vector<string> token_vector = tokenize (ituple_string, delimiter);

   for (vector<string>::const_iterator iterator = token_vector.begin ();
        iterator != token_vector.end (); iterator++)
   {
      const string& token = *(iterator);
      push_back (atoi (token.c_str ()));
   }

}

void
Ituple::add_content (const Integer number_of_values,
                     const Integer value,
                     const bool clear_first)
{
   if (clear_first) { clear (); }
   for (Integer i = 0; i < number_of_values; i++) { push_back (value); }
}

void
Ituple::add_content (const Integer number_of_values,
                     const Integer start_value,
                     const Integer step_value,
                     const bool clear_first)
{

   if (clear_first) { clear (); }

   for (Integer i = 0; i < number_of_values; i++)
   {
      push_back (start_value + i * step_value);
   }

}

void
Ituple::pop_front ()
{
   erase (begin ());
}
                
void
Ituple::push_front (const Integer value)
{
   insert (begin (), value);
}

Domain_1D::Domain_1D (const Real start,
                      const Real end)
   : start (start),
     end (end)
{
}

Domain_1D::Domain_1D (const Domain_1D& domain_1d)
   : start (domain_1d.start),
     end (domain_1d.end)
{
}

bool
Domain_1D::is_reverse () const
{
   return (start > end);
}

void
Domain_1D::swap ()
{
   std::swap (start, end);
}

void
Domain_1D::swap_if_reverse ()
{
   if (is_reverse ()) { swap (); }
}

Real
Domain_1D::get_span () const
{
   return end - start;
}

bool
Domain_1D::snap (Real& x) const
{
   bool snapped = false;
        if (x < start) { x = start; snapped = true; }
   else if (x > end) { x = end; snapped = true; }
   return snapped;
}

bool
Domain_1D::is_out_of_bounds (const Real x) const
{
   return ((x - start) * (x - end)) > 0;
}
         
void
Domain_1D::translate (const Real delta)
{
   start += delta;
   end += delta;
}

Real
Domain_1D::normalize (const Real x,
                      const Real gamma) const
{
   const Real f = std::min (std::max ((x - start) / (end - start), 0.0), 1.0);
   return pow (f, gamma);
}

Domain_2D::Domain_2D (const Real start_x,
                      const Real end_x,
                      const Real start_y,
                      const Real end_y)
 : domain_x (start_x, end_x),
   domain_y (start_y, end_y)
{
}

Domain_2D::Domain_2D (const Domain_1D& domain_x,
                      const Domain_1D& domain_y)
 : domain_x (domain_x),
   domain_y (domain_y)
{
}

Domain_2D::Domain_2D (const Domain_2D& domain_2d)
   : domain_x (domain_2d.domain_x),
     domain_y (domain_2d.domain_y)
{
}

Real
Domain_2D::get_width () const
{
   return domain_x.get_span ();
}

Real
Domain_2D::get_height () const
{
   return domain_y.get_span ();
}

bool
Domain_2D::is_out_of_bounds (const Real x,
                             const Real y) const
{
   return domain_x.is_out_of_bounds (x) || domain_y.is_out_of_bounds (y);
}

bool
Domain_2D::is_out_of_bounds (const Point_2D& point) const
{
   return is_out_of_bounds (point.x, point.y);
}

void
Domain_2D::translate (const Real delta_x,
                      const Real delta_y)
{
   domain_x.translate (delta_x);
   domain_y.translate (delta_y);
}

void
Domain_2D::clip_line (Point_2D& p0,
                      Point_2D& p1) const
{

   Real t;
   Real t0 = 0;
   Real t1 = 1;
   Real dx = p1.x - p0.x;
   Real dy = p1.y - p0.y;

   t = (domain_x.start - p0.x) / dx;
   if (dx > 0)      { if (t > t0) { t0 = t; }  }
   else if (dx < 0) { if (t < t1) { t1 = t; } }

   t = (domain_x.end - p0.x) / dx;
   if (dx < 0)      { if (t > t0) { t0 = t; }  }
   else if (dx > 0) { if (t < t1) { t1 = t; } }

   t = (domain_y.start - p0.y) / dy;
   if (dy > 0)      { if (t > t0) { t0 = t; }  }
   else if (dy < 0) { if (t < t1) { t1 = t; } }

   t = (domain_y.end - p0.y) / dy;
   if (dy < 0)      { if (t > t0) { t0 = t; }  }
   else if (dy > 0) { if (t < t1) { t1 = t; } }

   if (t0 > t1)
   {
      p1.x = GSL_NAN;
      p1.y = GSL_NAN;
      p0.x = GSL_NAN;
      p0.y = GSL_NAN;
   }

   else
   {
      p1.x = p0.x + t1 * dx;
      p1.y = p0.y + t1 * dy;
      p0.x = p0.x + t0 * dx;
      p0.y = p0.y + t0 * dy;
   }

}

Point_2D
Domain_2D::get_random_point () const
{

   const Real start_x = domain_x.start;
   const Real end_x = domain_x.end;
   const Real start_y = domain_y.start;
   const Real end_y = domain_y.end;

   return Point_2D (random (start_x, end_x), random (start_y, end_y));

}

Real
Domain_2D::get_span_x () const
{
   return domain_x.get_span ();
}

Real
Domain_2D::get_span_y () const
{
   return domain_y.get_span ();
}

Domain_3D::Domain_3D (const Real start_z,
                      const Real end_z,
                      const Real start_x,
                      const Real end_x,
                      const Real start_y,
                      const Real end_y)
   : domain_z (start_z, end_z),
     domain_x (start_x, end_x),
     domain_y (start_y, end_y)
{
}

Domain_3D::Domain_3D (const Domain_1D& domain_z,
                      const Domain_1D& domain_x,
                      const Domain_1D& domain_y)
   : domain_z (domain_z),
     domain_x (domain_x),
     domain_y (domain_y)
{
}

Domain_3D::Domain_3D (const Domain_3D& domain_3d)
   : domain_z (domain_3d.domain_z),
     domain_x (domain_3d.domain_x),
     domain_y (domain_3d.domain_y)
{
}

Real
Domain_3D::get_depth () const
{
   return domain_z.get_span ();
}

Real
Domain_3D::get_width () const
{
   return domain_x.get_span ();
}

Real
Domain_3D::get_height () const
{
   return domain_y.get_span ();
}

bool
Domain_3D::is_out_of_bounds (const Real z,
                             const Real x,
                             const Real y) const
{
   return domain_z.is_out_of_bounds (x) ||
          domain_y.is_out_of_bounds (y) ||
          domain_y.is_out_of_bounds (y);
}

bool
Domain_3D::is_out_of_bounds (const Point_3D& point) const
{
   return is_out_of_bounds (point.z, point.x, point.y);
}

void
Domain_3D::translate (const Real delta_z,
                      const Real delta_x,
                      const Real delta_y)
{
   domain_z.translate (delta_z);
   domain_x.translate (delta_x);
   domain_y.translate (delta_y);
}

namespace denise
{

   ostream&
   operator << (ostream& out_file,
                const Index_nD& index_nd)
   {

      out_file << "(";

      for (Integer i = 0; i < index_nd.n; i++)
      {

         out_file << index_nd.buffer[i];

         if (i != index_nd.n - 1)
         {
            out_file << ", ";
         }

      }

      out_file << ")";
      return out_file;

   }

   ostream&
   operator << (ostream& out_file,
                const Point_nD& point_nd)
   {

      out_file << "(";

      for (Integer i = 0; i < point_nd.n; i++)
      {

         out_file << point_nd.buffer[i];

         if (i != point_nd.n - 1)
         {
            out_file << ", ";
         }

      }

      out_file << ")";
      return out_file;

   }

   ostream&
   operator << (ostream &out_file,
                const Index_2D& index)
   {
      out_file << "(" << index.i << ", " << index.j << ")";
      return out_file;
   }

   ostream&
   operator << (ostream &out_file,
                const Index_3D& index)
   {
      out_file << "(" << index.k << ", " << index.i << ", " << index.j << ")";
      return out_file;
   }

   ostream&
   operator << (ostream &out_file,
                const Box_2D& box_2d)
   {
      out_file << "(" << box_2d.index_2d << ", " << box_2d.size_2d << ")";
      return out_file;
   }

   ostream&
   operator << (ostream &out_file,
                const Point_2D& point)
   {
      out_file << "(" << point.x << ", " << point.y << ")";
      return out_file;
   }

   ostream&
   operator << (ostream &out_file,
                const Point_3D& point)
   {
      out_file << "(" << point.z << ", " << point.x << ", " << point.y << ")";
      return out_file;
   }

   ostream&
   operator << (ostream &out_file,
                const Tuple& tuple)
   {

      for (Tuple::const_iterator iterator = tuple.begin ();
           iterator != tuple.end (); iterator++)
      {
         const Real& component = *(iterator);
         if (distance (iterator, tuple.end ()) == 1) { out_file << component; }
         else { out_file << component << ":"; }
      }

      return out_file;

   }

   ostream&
   operator << (ostream &out_file,
                const Ituple& ituple)
   {

      for (Ituple::const_iterator iterator = ituple.begin ();
           iterator != ituple.end (); iterator++)
      {
         const Real& component = *(iterator);
         if (distance (iterator, ituple.end ()) == 1) { out_file << component; }
         else { out_file << component << ":"; }
      }

      return out_file;

   }

   ostream&
   operator << (ostream &out_file,
                const Domain_1D& domain)
   {
      out_file << "(" << domain.start << " -> " << domain.end << ")";
      return out_file;
   }

}
