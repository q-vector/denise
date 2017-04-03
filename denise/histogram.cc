//
// histogram.cc
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

#include "histogram.h"

using namespace std;
using namespace denise;

Histogram::Axis::Axis (const Real interval,
                       const Real offset)
   : interval (interval),
     offset (offset)
{
}
   
Histogram::Axis::Axis (const set<Real>& real_set)
   : set<Real> (real_set),
     interval (GSL_NAN),
     offset (GSL_NAN)
{
}

bool
Histogram::Axis::is_flexi () const
{
   return !gsl_isnan (interval);
}

void
Histogram::Axis::extend (const Real x)
{

   const bool is_empty = (size () == 0);

   if (is_empty)
   {
      const Real lower = x - modulo ((x - offset), interval);
      const Real upper = lower + interval;
      insert (lower);
      insert (upper);
      return;
   }

   // Note the asymmetry
   //    less than or equal to    vs   greater than

   const bool x_too_small = (x <= *(begin ()));
   if (x_too_small)
   {
      const Real half_interval = interval / 2;
      const Real lowest = x - modulo ((x - offset), interval);
      const Real end = *(begin ()) - half_interval;
      for (Real b = lowest; b < end; b += interval)
      {
         insert (b);
      }
      return;
   }

   const bool x_too_large = (x > *(rbegin ()));
   if (x_too_large)
   {
      const Real half_interval = interval / 2;
      const Real start = *(rbegin ()) + interval;
      const Real uppest = x - modulo ((x - offset), interval) + interval;
      for (Real b = start; b < uppest + half_interval; b += interval)
      {
         insert (b);
      }
      return;
   }

}

Real
Histogram::Axis::get_bin_top (const Real x) const
{
   return (*(lower_bound (x)));
}

Real
Histogram::Axis::get_bin_bottom (const Real x) const
{
   return (*(--(lower_bound (x))));
}

Histogram_1D::Histogram_1D (const Real interval,
                            const Real offset)
   : axis (interval, offset),
     max_value (0),
     min_value (0)
{
}
   
Histogram_1D::Histogram_1D (const set<Real>& real_set)
   : axis (real_set),
     max_value (0),
     min_value (0)
{
}

const Histogram::Axis&
Histogram_1D::get_axis () const
{
   return axis;
}

Real
Histogram_1D::get_bin_top (const Real x) const
{
   return axis.get_bin_top (x);
}

Real
Histogram_1D::get_bin_bottom (const Real x) const
{
   return axis.get_bin_bottom (x);
}

Real
Histogram_1D::get_value (const Real x) const
{
   const Real bin_top = axis.get_bin_top (x);
   Histogram_1D::const_iterator i = find (bin_top);
   return (i == end () ? 0 : i->second);
}

const Real&
Histogram_1D::get_max_value () const
{
   return max_value;
}

const Real&
Histogram_1D::get_min_value () const
{
   return min_value;
}

void
Histogram_1D::clear ()
{
   axis.clear ();
   std::map<Real, Real>::clear ();
   max_value = 0;
   min_value = 0;
}

void
Histogram_1D::increment (const Real x,
                         const Real weight)
{

   if (gsl_isnan (x)) { return; }
   if (axis.is_flexi ()) { axis.extend (x); }

   const Real bin_top = axis.get_bin_top (x);
   Histogram_1D::iterator i = find (bin_top);

   if (i == end ())
   {
      insert (make_pair (bin_top, weight));
      if (max_value < weight) { max_value = weight; }
      if (min_value < weight) { min_value = weight; }
   }
   else
   {
      i->second += weight;
      if (max_value < i->second) { max_value = i->second; }
      if (min_value < i->second) { min_value = i->second; }
   }

}

void
Histogram_1D::render (const RefPtr<Context> cr,
                      const Transform_2D& transform,
                      const Domain_1D& domain_y,
                      const Dstring& bin_format,
                      const Dstring& value_format,
                      const Color& color,
                      const Color& value_color,
                      const Color& baseline_color) const
{

   Polygon polygon;

   for (Axis::const_iterator iterator = axis.begin ();
        iterator != axis.end (); iterator++)
   {

      Axis::const_iterator next = iterator; next++;
      const bool first = (iterator == axis.begin ());
      const bool last = (next == axis.end ());

      const Real bin_bottom = *(iterator);
      const Real bin_top = *(next);
      const Real bin = (bin_top + bin_bottom) / 2;

      const Real value = get_value (bin_top);

      if (first)
      {

         polygon.add (Point_2D (0, bin_bottom));

         const Point_2D bin_bottom_p (0, bin_bottom);
         const Point_2D& tick_p = transform.transform (bin_bottom_p);
         const Edge tick (tick_p, tick_p + Point_2D (5, 0));

         baseline_color.cairo (cr);
         tick.cairo (cr);
         cr->stroke ();

      }

      if (!last)
      {

         const Dstring& bin_str = Dstring::render (bin_format, bin);
         const Dstring& value_str = Dstring::render (value_format, value);

         const Point_2D bin_point (0, bin);
         const Point_2D value_point (value, bin);

         baseline_color.cairo (cr);
         Label bin_label (bin_str, bin_point, 'r', 'c', 5);
         bin_label.cairo (cr, transform);

         value_color.cairo (cr);
         baseline_color.cairo (cr);
         Label value_label (value_str, value_point, 'l', 'c', 5);
         value_label.cairo (cr, transform);

         polygon.add (Point_2D (value, bin_bottom));
         polygon.add (Point_2D (value, bin_top));

         const Point_2D bin_top_p (0, bin_top);
         const Point_2D& tick_p = transform.transform (bin_top_p);
         const Edge tick (tick_p, tick_p + Point_2D (5, 0));

         baseline_color.cairo (cr);
         tick.cairo (cr);
         cr->stroke ();

      }
      else
      {
         polygon.add (Point_2D (0, bin_bottom));
      }

   }

   polygon.cairo (cr, transform);
   color.cairo (cr);
   cr->fill ();

   const Real start_y = *(axis.begin ());
   const Real end_y = *(axis.rbegin ());
   const Edge baseline (Point_2D (0, start_y), Point_2D (0, end_y));
   baseline_color.cairo (cr);
   baseline.cairo (cr, transform);
   cr->stroke ();

}

Histogram_2D::Histogram_2D (const Real interval_x,
                            const Real interval_y,
                            const Real offset_x,
                            const Real offset_y)
   : axis_x (interval_x, offset_x),
     axis_y (interval_y, offset_y),
     max_value (0),
     min_value (0)
{
}

Histogram_2D::Histogram_2D (const set<Real>& real_set_x,
                            const set<Real>& real_set_y)
   : axis_x (real_set_x),
     axis_y (real_set_y),
     max_value (0),
     min_value (0)
{
}

const Histogram::Axis&
Histogram_2D::get_axis_x () const
{
   return axis_x;
}

const Histogram::Axis&
Histogram_2D::get_axis_y () const
{
   return axis_y;
}

Real
Histogram_2D::get_bin_top_x (const Real x) const
{
   return axis_x.get_bin_top (x);
}

Real
Histogram_2D::get_bin_top_y (const Real y) const
{
   return axis_x.get_bin_top (y);
}

Real
Histogram_2D::get_bin_bottom_x (const Real x) const
{
   return axis_y.get_bin_bottom (x);
}

Real
Histogram_2D::get_bin_bottom_y (const Real y) const
{
   return axis_y.get_bin_bottom (y);
}

Real
Histogram_2D::get_value (const Point_2D& point_2d) const
{
   const Real bin_top_x = axis_x.get_bin_top (point_2d.x);
   const Real bin_top_y = axis_y.get_bin_top (point_2d.y);
   const Point_2D p (bin_top_x, bin_top_y);
   Histogram_2D::const_iterator i = find (p);
   return (i == end () ? 0 : i->second);
}

const Real&
Histogram_2D::get_max_value () const
{
   return max_value;
}

const Real&
Histogram_2D::get_min_value () const
{
   return min_value;
}

void
Histogram_2D::clear ()
{
   axis_x.clear ();
   axis_y.clear ();
   std::map<Point_2D, Real>::clear ();
   max_value = 0;
   min_value = 0;
}

void
Histogram_2D::increment (const Point_2D& point_2d,
                         const Real weight)
{

   const Real x = point_2d.x;
   const Real y = point_2d.y;

   if (gsl_isnan (x) || gsl_isnan (y)) { return; }
   if (axis_x.is_flexi ()) { axis_x.extend (x); }
   if (axis_y.is_flexi ()) { axis_y.extend (y); }

   const Real bin_top_x = axis_x.get_bin_top (x);
   const Real bin_top_y = axis_y.get_bin_top (y);
   const Point_2D p (bin_top_x, bin_top_y);
   Histogram_2D::iterator i = find (p);

   if (i == end ())
   {
      insert (make_pair (p, weight));
      if (max_value < weight) { max_value = weight; }
      if (min_value < weight) { min_value = weight; }
   }
   else
   {
      i->second += weight;
      if (max_value < i->second) { max_value = i->second; }
      if (min_value < i->second) { min_value = i->second; }
   }

}

void
Histogram_2D::render (const RefPtr<Context> cr,
                      const Transform_2D& transform,
                      const Color_Chooser& color_chooser,
                      const Dstring& format,
                      const Color& fg_color,
                      const Color& box_border_color) const
{

   cr->save ();
   cr->set_antialias (ANTIALIAS_NONE);

   Point_2D p;

   for (Axis::const_iterator i = axis_x.begin ();
        i != axis_x.end (); i++)
   {

      Axis::const_iterator next_i = i; next_i++;
      const bool last_x = (next_i == axis_x.end ());
      if (last_x) { continue; }

      const Real xb = *(i);
      const Real xe = *(next_i);
      p.x = (xb + xe) / 2;

      for (Axis::const_iterator j = axis_y.begin ();
           j != axis_y.end (); j++)
      {

         Axis::const_iterator next_j = j; next_j++;
         const bool last_y = (next_j == axis_y.end ());
         if (last_y) { continue; }

         const Real yb = *(j);
         const Real ye = *(next_j);
         p.y = (yb + ye) / 2;

         const Real value = get_value (p);
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

         if (value == 0) { cr->fill (); }
         else
         {
            cr->fill_preserve ();
            const Dstring& str = Dstring::render (format, value);
            box_border_color.cairo (cr);
            cr->stroke ();
            fg_color.cairo (cr);
            Label (str, transform.transform (p), 'c', 'c').cairo (cr);
         }

      }
   }

   cr->restore ();

}

