//
// transform.h
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

#ifndef DENISE_TRANSFORM_H
#define DENISE_TRANSFORM_H

#include <cairomm/context.h>
#include <denise/basics.h>
#include <denise/dstring.h>

using namespace std;

namespace denise
{

   class Edge;
   class Circle;
   class Ellipse;
   class Simple_Polyline;
   class Polyline;
   class Polygon;

   class Transform_nD
   {

      public:

         virtual void
         transform (Point_nD& abscissa,
                    const Point_nD& ordinate) const = 0;

         virtual void
         reverse (Point_nD& ordinate,
                  const Point_nD& abscissa) const = 0;

   };

   class Transform_1D
   {

      public:

         virtual bool
         out_of_domain (const Real x) const = 0;

         virtual void
         transform (Real& transformed,
                    const Real x) const = 0;

         virtual Real
         transform (const Real x) const;

         virtual void
         reverse (Real& reversed,
                  const Real x) const = 0;

         virtual Real
         reverse (const Real x) const;

   };

   class Transform_2D
   {

      protected:

         Polygon*
         domain_polygon_ptr;

         void
         cr (const Cairo::RefPtr<Cairo::Context>& cr,
             const Simple_Polyline& simple_polyline) const;

         void
         cr (const Cairo::RefPtr<Cairo::Context>& cr,
             const Simple_Polyline& simple_polyline,
             const Polygon& clip_polygon,
             const bool negative) const;

         void
         cr (const Cairo::RefPtr<Cairo::Context>& cr,
             const Polyline& polyline) const;

         void
         cr (const Cairo::RefPtr<Cairo::Context>& cr,
             const Polyline& polyline,
             const Polygon& clip_polygon,
             const bool negative) const;

         void
         cr (const Cairo::RefPtr<Cairo::Context>& cr,
             const Polygon& polygon) const;

      public:

         Transform_2D ();

         static Transform_2D*
         get_transform_ptr (const Dstring& str);

         virtual bool
         out_of_domain (const Real x, 
                        const Real y) const;

         virtual bool
         is_out_of_domain (const Point_2D& point_2d) const;

         virtual void
         transform (Real& transformed_x,
                    Real& transformed_y,
                    const Real x,
                    const Real y) const;

         virtual Point_2D
         transform (const Real x,
                    const Real y) const;

         virtual Point_2D
         transform (const Point_2D& point_2d) const;

         virtual void
         reverse (Real& reversed_x,
                  Real& reversed_y,
                  const Real x,
                  const Real y) const;

         virtual Point_2D
         reverse (const Real x,
                  const Real y) const;

         virtual Point_2D
         reverse (const Point_2D& point_2d) const;

         virtual void
         transform_uv (Real& u,
                       Real& v,
                       const Real x,
                       const Real y) const;

         virtual void
         transform_uv (Real& u,
                       Real& v,
                       const Point_2D& point_2d) const;

         virtual Real
         get_theta (const Real u,
                    const Real v,
                    const Real x,
                    const Real y) const;

         virtual Real
         get_theta (const Real u,
                    const Real v,
                    const Point_2D& point_2d) const;

         virtual void
         cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                const Edge& edge) const;

         virtual void
         cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                const Circle& circle) const;

         virtual void
         cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                const Ellipse& ellipse) const;

         virtual void
         cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                const Simple_Polyline& simple_polyline) const;

         virtual void
         cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                const Simple_Polyline& simple_polyline,
                const Polygon& clip_polygon,
                const bool negative = false) const;

         virtual void
         cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                const Polyline& polyline) const;

         virtual void
         cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                const Polygon& polygon) const;

   };

   class Affine_Transform_1D : public Transform_1D
   {

      protected:

         Real
         s;

         Real
         o;

      public:

         Affine_Transform_1D ();

         Affine_Transform_1D (const Real scale,
                              const Real offset);

         Affine_Transform_1D (const Affine_Transform_1D& transform);

         virtual bool
         out_of_domain (const Real x) const;

         const Real&
         get_scale () const;

         const Real&
         get_offset () const;

         void
         scale (const Real scale);

         void
         scale (const Real scale,
                const Real pivot);

         void
         translate (const Real translation);

         virtual void
         transform (Real& transformed,
                    const Real x) const;

         virtual Real
         transform (const Real x) const;

         virtual void
         reverse (Real& reversed,
                  const Real x) const;

         virtual Real
         reverse (const Real x) const;

   };

   class Moebius_Transform : public Transform_2D
   {

      protected:

         complex<Real>
         a;

         complex<Real>
         b;

         complex<Real>
         c;

         complex<Real>
         d;

      public:

         Moebius_Transform ();

         Moebius_Transform (const Moebius_Transform& mt);

         Moebius_Transform (const Moebius_Transform& mt_a,
                            const Moebius_Transform& mt_b);

         Moebius_Transform (const complex<Real>& a,
                            const complex<Real>& b,
                            const complex<Real>& c,
                            const complex<Real>& d);

         Moebius_Transform (const Point_2D& a,
                            const Point_2D& b,
                            const Point_2D& c,
                            const Point_2D& d);

         Moebius_Transform (const complex<Real>& z_a,
                            const complex<Real>& z_b,
                            const complex<Real>& z_c);

         Moebius_Transform (const complex<Real>& z_a,
                            const complex<Real>& z_b,
                            const complex<Real>& z_c,
                            const complex<Real>& w_a,
                            const complex<Real>& w_b,
                            const complex<Real>& w_c);

         void
         set_identity ();

         Moebius_Transform
         get_inverse () const;

         Moebius_Transform
         get_compose (const Moebius_Transform& mt) const;

         pair<Point_2D, Point_2D>
         get_fixed_points () const;

         virtual void
         transform (Real& transformed_x,
                    Real& transformed_y,
                    const Real x,
                    const Real y) const;

         virtual Point_2D
         transform (const Real x,
                    const Real y) const;

         virtual Point_2D
         transform (const Point_2D& point_2d) const;

         virtual void
         reverse (Real& reversed_x,
                  Real& reversed_y,
                  const Real x,
                  const Real y) const;

         virtual Point_2D
         reverse (const Real x,
                  const Real y) const;

         virtual Point_2D
         reverse (const Point_2D& point_2d) const;

         void
         transform_uv (Real& u,
                       Real& v,
                       const Real x,
                       const Real y) const;

         Real
         get_theta (const Real u,
                    const Real v,
                    const Real x,
                    const Real y) const;

   };


   /// Represents a 2D affine transform that performs a linear mapping
   /// from 2D coordinates to other 2D coordinates that preserves 
   /// collinearity and ratios of distances.
   ///
   /// Affine transformations can be constructed using sequences of
   /// translations, scales, flips, rotations, and shears. Such a
   /// coordinate transformation can be represented by a 3 row by
   /// 3 column matrix with an implied last row of \f$[ 0 0 1 ]\f$.
   /// This matrix transforms source coordinates \f$(x, y)\f$ into
   /// destination coordinates \f$(x', y')\f$ by considering them
   /// to be a column vector and multiplying the coordinate vector
   /// by the matrix according to the following process:
   ///
   /// \f[
   ///    \left( \begin{array}{c} x\prime \\ y\prime \\ 1 \end{array} \right) =
   ///    \left( \begin{array}{ccc}
   ///       m_{00} & m_{01} & m_{02} \\ m_{10} & m_{11} & m_{12} \\ 0 & 0 & 1
   ///    \end{array} \right)
   ///    \left( \begin{array}{c} x \\ y \\ 1 \end{array} \right) =
   ///    \left( \begin{array}{ccc}
   ///       m_{00} x + m_{01} y + m_{02} \\ m_{10} x + m_{11} y + m_12} \\ 1
   ///    \end{array} \right)
   /// \f]
   ///
   class Affine_Transform_2D : public Transform_2D
   {

      private:

         Real
         m[2][3];

         Real
         r[2];

         Real
         determinant;

         void
         init ();

      public:

         Affine_Transform_2D ();

         Affine_Transform_2D (const Affine_Transform_2D& transform);

         void
         set_identity ();

         void
         rotate (const Real theta);

         void
         rotate (const Real theta,
                 const Point_2D& pivot);

         void
         scale (const Real scale_x,
                const Real scale_y);

         void
         scale (const Real scale_x,
                const Real scale_y,
                const Point_2D& pivot);

         void
         shear_x (const Real shear_x);

         void
         translate (const Real translation_x,
                    const Real translation_y);

         virtual Real
         get_jacobian () const;

         const Real&
         get_scale_x () const;

         const Real&
         get_scale_y () const;

         const Real&
         get_shear_x () const;

         const Real&
         get_shear_y () const;

         const Real&
         get_translation_x () const;

         const Real&
         get_translation_y () const;

         virtual void
         transform (Real& transformed_x,
                    Real& transformed_y,
                    const Real x,
                    const Real y) const;

         virtual Point_2D
         transform (const Real x,
                    const Real y) const;

         virtual Point_2D
         transform (const Point_2D& point_2d) const;

         virtual void
         reverse (Real& reversed_x,
                  Real& reversed_y,
                  const Real x,
                  const Real y) const;

         virtual Point_2D
         reverse (const Real x,
                  const Real y) const;

         virtual Point_2D
         reverse (const Point_2D& point_2d) const;

         void
         transform_uv (Real& u,
                       Real& v,
                       const Real x,
                       const Real y) const;

         Real
         get_theta (const Real u,
                    const Real v,
                    const Real x,
                    const Real y) const;

   };

}

#endif /* DENISE_TRANSFORM_H */

