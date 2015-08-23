//
// fractal.h
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

#ifndef DENISE_FRACTAL_H
#define DENISE_FRACTAL_H

#include <cairomm/context.h>
#include <denise/analysis.h>
#include <denise/exception.h>
#include <denise/geometry.h>
#include <denise/graphics.h>

using namespace std;

namespace denise
{

   class Newton
   {

      private:

         static Point_2D
         power (const Point_2D& z,
                const Real power);

         static Real
         f (Real& x,
            Real& y,
            const Real power);

      public:

         class Raster : public denise::Raster
         {

            public:

               Raster (const Transform_2D& transform,
                       const Box_2D& box_2d,
                       const Real n,
                       const Color_Chooser& color_chooser,
                       const Integer maxiter);

         };

         static Integer
         get_iterations (const Point_2D& c,
                         const Real n,
                         const Integer maxiter,
                         const Real epsilon = 1e-6);

         static pair<Real, Integer>
         get_theta_iterations (const Point_2D& c,
                               const Real n,
                               const Integer maxiter,
                               const Real epsilon = 1e-6);

   };

   class Julia
   {

      private:

         static void
         f (Real& x,
            Real& y,
            const Point_2D& c);

      public:

         class Raster : public denise::Raster
         {

            public:

               Raster (const Transform_2D& transform,
                       const Box_2D& box_2d,
                       const Point_2D& c,
                       const Color_Chooser& color_chooser,
                       const Integer maxiter);

         };

         static Integer
         get_iterations (const Point_2D& z,
                         const Point_2D& c,
                         const Integer maxiter,
                         const Real bound);

   };

   class Mandelbrot
   {

      private:

         static void
         f (Real& x,
            Real& y,
            const Point_2D& c);

      public:

         class Raster : public denise::Raster
         {

            public:

               Raster (const Transform_2D& transform,
                       const Box_2D& box_2d,
                       const Color_Chooser& color_chooser,
                       const Integer maxiter);

         };

         static Real
         get_value (const Point_2D& c,
                    const Integer maxiter,
                    const Real bound);

         static Integer
         get_iterations (const Point_2D& c,
                         const Integer maxiter,
                         const Real bound);

   };

}

#endif /* DENISE_FRACTAL_H */

