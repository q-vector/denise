//
// fractal.cc
// 
// Copyright (C) 2010-2013 Simon E. Ching
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

#include "fractal.h"

using namespace std;
using namespace denise;

Point_2D
Newton::power (const Point_2D& z,
               const Real n)
{

   const Real theta = atan2 (z.y, z.x);
   const Real r = sqrt (z.x*z.x + z.y*z.y);

   const Real rho = pow (double (r), double (n));
   const Real phi = theta * n;

   return Point_2D (rho * cos (phi), rho * sin (phi));

}

Real
Newton::f (Real& x,
           Real& y,
           const Real n)
{

   const Point_2D z (x, y);
   const Point_2D& delta = (power (z, 1 - n) - z) / n;

   x += delta.x;
   y += delta.y;

   return sqrt (delta.x*delta.x + delta.y*delta.y);

}

Newton::Raster::Raster (const Transform_2D& transform,
                        const Box_2D& box_2d,
                        const Real n,
                        const Color_Chooser& color_chooser,
                        const Integer maxiter)
   : denise::Raster (box_2d)
{

   const Size_2D& size_2d = box_2d.size_2d;
   const Index_2D& anchor = box_2d.index_2d;
   const Index_2D end_index (anchor.i + size_2d.i, anchor.j  + size_2d.j);

   #pragma omp parallel for
   for (Integer i = anchor.i; i < end_index.i; i++)
   {
      for (Integer j = anchor.j; j < end_index.j; j++)
      {
         const Point_2D& z = transform.reverse (Real (i), Real (j));
         pair<Real, Integer> theta_iter = get_theta_iterations (z, n, maxiter);
         const Real theta = theta_iter.first;
         const Integer iter = theta_iter.second;
         const Color& color = color_chooser.get_color (Real (iter));
         //const Real hue = modulo ((theta + M_PI) / (2*M_PI) + 0.500001, 1)* 0.9;
         const Real hue = modulo ((theta) / (2*M_PI) + 0.0, 1);
         const Real saturation = pow (color.get_saturation (), 0.7);
         const Real brightness = pow (color.get_brightness (), 0.9);
         const Color& final_color = Color::hsb (hue, saturation, brightness);
         set_pixel (i - anchor.i, j - anchor.j, final_color);
      }
   }

}

Integer
Newton::get_iterations (const Point_2D& z,
                        const Real n,
                        const Integer maxiter,
                        const Real epsilon)
{

   Integer i = 0;
   Real x = z.x, y = z.y;
   Real step = GSL_POSINF;

   while (step > epsilon) 
   {
      step = f (x, y, n);
      i++;
      if (i == maxiter) { return maxiter; }
   }

   return i;

}

pair<Real, Integer>
Newton::get_theta_iterations (const Point_2D& z,
                              const Real n,
                              const Integer maxiter,
                              const Real epsilon)
{

   Integer i = 0;
   Real x = z.x, y = z.y;
   Real step = GSL_POSINF;

   while (step > epsilon)
   {
      step = f (x, y, n);
      i++;
      if (i == maxiter) { return make_pair (0, Real (maxiter)); }
   }

   return make_pair (atan2 (y, x), i);

}

void
Julia::f (Real& x,
          Real& y,
          const Point_2D& c)
{
   const Real xx = x*x - y*y + c.x;
   y = 2*x*y + c.y;
   x = xx;
}

Julia::Raster::Raster (const Transform_2D& transform,
                       const Box_2D& box_2d,
                       const Point_2D& c,
                       const Color_Chooser& color_chooser,
                       const Integer maxiter)
   : denise::Raster (box_2d)
{

   const Size_2D& size_2d = box_2d.size_2d;
   const Index_2D& anchor = box_2d.index_2d;
   const Index_2D end_index (anchor.i + size_2d.i, anchor.j  + size_2d.j);

   const Real bound = ((1 + sqrt (1 + 4 * sqrt (c.x*c.x + c.y*c.y))) / 2);

   const Color black (0, 0, 0);
   const Color white (1, 1, 1);

   #pragma omp parallel for
   for (Integer i = anchor.i; i < end_index.i; i++)
   {
      for (Integer j = anchor.j; j < end_index.j; j++)
      {
         const Point_2D& z = transform.reverse (Real (i), Real (j));
         const Real iter = get_iterations (z, c, maxiter, bound);
         const bool b = (iter == maxiter);
         const Color& color = (b ? black : color_chooser.get_color (iter));
         set_pixel (i - anchor.i, j - anchor.j, color);
      }
   }

}

Integer
Julia::get_iterations (const Point_2D& z,
                       const Point_2D& c,
                       const Integer maxiter,
                       const Real bound)
{

   Integer i = 0;
   Real x = z.x, y = z.y;
   const Real b = bound*bound;

   for (i = 0; (x*x + y*y) < b; i++)
   {
      if (i == maxiter) { return maxiter; }
      f (x, y, c);
   }

   return i;

}

void
Mandelbrot::f (Real& x,
               Real& y,
               const Point_2D& c)
{
   const Real xx = x*x - y*y + c.x;
   y = 2*x*y + c.y;
   x = xx;
}

Mandelbrot::Raster::Raster (const Transform_2D& transform,
                            const Box_2D& box_2d,
                            const Color_Chooser& color_chooser,
                            const Integer maxiter)
   : denise::Raster (box_2d)
{

   const Size_2D& size_2d = box_2d.size_2d;
   const Index_2D& anchor = box_2d.index_2d;
   const Index_2D end_index (anchor.i + size_2d.i, anchor.j  + size_2d.j);
   const Real bound = 4;

   const Color black (0, 0, 0);

   #pragma omp parallel for
   for (Integer i = anchor.i; i < end_index.i; i++)
   {
      for (Integer j = anchor.j; j < end_index.j; j++)
      {
         const Point_2D& c = transform.reverse (Real (i), Real (j));
         const Real iter = get_value (c, maxiter, bound);
         const bool b = (iter > (maxiter - 0.5));
         const Color& color = (b ? black : color_chooser.get_color (iter));
         set_pixel (i - anchor.i, j - anchor.j, color);
      }
   }

}

Real
Mandelbrot::get_value (const Point_2D& c,
                       const Integer maxiter,
                       const Real bound)
{

   Integer i = 0;
   Real x = c.x, y = c.y;
   const Real b = bound*bound;

   for (i = 0; x*x + y*y < b; i++)
   {
      if (i == maxiter) { return maxiter; }
      f (x, y, c);
   }

   const Real zn = sqrt (x*x + y*y);
   const Real nu = log (log (zn) / M_LN2) / M_LN2;
   return i + Real (1 - nu);

}

Integer
Mandelbrot::get_iterations (const Point_2D& c,
                            const Integer maxiter,
                            const Real bound)
{

   Integer i = 0;
   Real x = c.x, y = c.y;
   const Real b = bound*bound;

   for (i = 0; x*x + y*y < b; i++)
   {
      if (i == maxiter) { return maxiter; }
      f (x, y, c);
   }

   return i;

}

