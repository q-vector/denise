//
// histogram.h
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

#ifndef DENISE_HISTOGRAM_H
#define DENISE_HISTOGRAM_H

#include <cmath>
#include <map>
#include <set>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <denise/basics.h>
#include <denise/geometry.h>
#include <denise/graphics.h>
#include <denise/linalg.h>
#include <denise/dataset.h>

using namespace std;

namespace denise
{

   class Histogram
   {

      protected:

         Integer
         number_of_points;

      public:

         class Axis : public set<Real>
         {

            protected:

               const Real
               interval;

               const Real
               offset;

            public:

               Axis (const Real interval,
                     const Real offset = 0);

               Axis (const set<Real>& real_set);

               bool
               is_flexi () const;

               void
               extend (const Real x);

               Real
               get_bin_top (const Real x) const;

               Real
               get_bin_bottom (const Real x) const;

         };

         Histogram ();

         virtual void
         increment ();

         const Integer&
         get_number_of_points () const;

   };

   class Histogram_1D : public map<Real, Real>,
                        public Histogram
   {

      protected:

         Axis
         axis;

         Real
         max_value;

         Real
         min_value;

      public:

         Histogram_1D (const Real interval,
                       const Real offset = 0);

         Histogram_1D (const set<Real>& real_set);

         const Histogram::Axis&
         get_axis () const;

         Real
         get_bin_top (const Real x) const;

         Real
         get_bin_bottom (const Real x) const;

         Real
         get_value (const Real x) const;

         const Real&
         get_max_value () const;

         const Real&
         get_min_value () const;

         void
         clear ();

         void
         increment (const Real x,
                    const Real weight = 1);

         void
         render (const RefPtr<Context>& cr,
                 const Transform_2D& transform,
                 const Dstring& bin_fmt,
                 const Dstring& value_format,
                 const Color& color,
                 const Color& value_color,
                 const Color& baseline_color) const;

         void
         render_outline (const RefPtr<Context>& cr,
                         const Transform_2D& transform) const;

   };

   class Histogram_2D : public map<Point_2D, Real>,
                        public Histogram
   {

      protected:

         Axis
         axis_x;

         Axis
         axis_y;

         Real
         max_value;

         Real
         min_value;

      public:

         Histogram_2D (const Real interval_x,
                       const Real interval_y,
                       const Real offset_x = 0,
                       const Real offset_y = 0);

         Histogram_2D (const set<Real>& real_set_x,
                       const set<Real>& real_set_y);

         const Histogram::Axis&
         get_axis_x () const;

         const Histogram::Axis&
         get_axis_y () const;

         Real
         get_bin_top_x (const Real x) const;

         Real
         get_bin_top_y (const Real y) const;

         Real
         get_bin_bottom_x (const Real x) const;

         Real
         get_bin_bottom_y (const Real y) const;

         Real
         get_value (const Point_2D& p) const;

         const Real&
         get_max_value () const;

         const Real&
         get_min_value () const;

         void
         clear ();

         void
         increment (const Point_2D& point_2d,
                    const Real weight = 1);

         void
         render (const RefPtr<Context> cr,
                 const Transform_2D& transform,
                 const Color_Chooser& color_chooser,
                 const Dstring& format,
                 const Color& fg_color,
                 const Color& box_border_color) const;

   };

}

#endif /* DENISE_HISTOGRAM_H */

