//
// basics.cc
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

#include <cstring>
#include "analysis.h"
#include "transform.h"
#include "geodesy.h"
#include "geometry.h"

using namespace denise;

void
Transform_2D::cr (const Cairo::RefPtr<Cairo::Context>& cr,
                  const Simple_Polyline& simple_polyline) const
{

   Real x, y;

   for (Simple_Polyline::const_iterator iterator = simple_polyline.begin ();
        iterator != simple_polyline.end (); iterator++)
   {

      const Point_2D& point = *(iterator);
      transform (x, y, point.x, point.y);

      if (iterator == simple_polyline.begin ())
      {
         cr->move_to (x, y);
      }
      else
      {
         cr->line_to (x, y);
      }

   }

}

void
Transform_2D::cr (const Cairo::RefPtr<Cairo::Context>& cr,
                  const Simple_Polyline& simple_polyline,
                  const Polygon& clip_polygon,
                  const bool negative) const
{

   Point_2D p;
   Simple_Polyline sp;

   for (Simple_Polyline::const_iterator iterator = simple_polyline.begin ();
        iterator != simple_polyline.end (); iterator++)
   {
      const Point_2D& point = *(iterator);
      transform (p.x, p.y, point.x, point.y);
      sp.add (p);
   }

   Polyline* p_ptr = sp.clip (clip_polygon, negative);
   p_ptr->cairo (cr);
   delete p_ptr;

}

void
Transform_2D::cr (const Cairo::RefPtr<Cairo::Context>& cr,
                  const Polyline& polyline) const
{
   for (Polyline::const_iterator iterator = polyline.begin ();
        iterator != polyline.end (); iterator++)
   {
      const Simple_Polyline& simple_polyline = **(iterator);
      Transform_2D::cr (cr, simple_polyline);
   }
}

void
Transform_2D::cr (const Cairo::RefPtr<Cairo::Context>& cr,
                  const Polyline& polyline,
                  const Polygon& clip_polygon,
                  const bool negative) const
{
   for (Polyline::const_iterator iterator = polyline.begin ();
        iterator != polyline.end (); iterator++)
   {
      const Simple_Polyline& simple_polyline = **(iterator);
      Transform_2D::cr (cr, simple_polyline, clip_polygon, negative);
   }
}

void
Transform_2D::cr (const Cairo::RefPtr<Cairo::Context>& cr,
                  const Polygon& polygon) const
{

   Real x, y;

   if (polygon.size () <= 1) { return; }

   const Polygon::Vertex* first_handle_ptr = polygon.get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {

      const Integer n = current_handle_ptr->n;
      Polygon::Vertex* current_ptr = (Polygon::Vertex*)current_handle_ptr;

      for (Integer i = 0; i < n; i++)
      {
         const Point_2D& point = (const Point_2D&)(*current_ptr);
         transform (x, y, point.x, point.y);
         if (i == 0) { cr->move_to (x, y); }
         else        { cr->line_to (x, y); }
         current_ptr = current_ptr->next_ptr;
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);

   cr->close_path ();

}

Real
Transform_1D::transform (const Real x) const
{
   Real t;
   transform (t, x);
   return t;
}

Real
Transform_1D::reverse (const Real x) const
{
   Real r;
   reverse (r, x);
   return r;
}

Transform_2D::Transform_2D ()
   : domain_polygon_ptr (NULL)
{
}

Transform_2D*
Transform_2D::get_transform_ptr (const Dstring& str)
{

   const Tokens tokens (str, ":");

   if (tokens[0] == "CARTESIAN")
   {
      const Real domain_x_start = stof (tokens[1]);
      const Real domain_x_end = stof (tokens[2]);
      const Real domain_y_start = stof (tokens[3]);
      const Real domain_y_end = stof (tokens[4]);
      const Real width = stof (tokens[5]);
      const Real height = stof (tokens[6]);
      const Real origin_x = stof (tokens[7]);
      const Real origin_y = stof (tokens[8]);
      const Domain_1D domain_x (domain_x_start, domain_x_end);
      const Domain_1D domain_y (domain_y_start, domain_y_end);
      const Point_2D origin (origin_x, origin_y);
      return new Cartesian_Transform_2D (domain_x,
         domain_y, width, height, origin);
   }
   else
   if (tokens[0] == "LOG")
   {
      const Real domain_x_start = stof (tokens[1]);
      const Real domain_x_end = stof (tokens[2]);
      const Real domain_y_start = stof (tokens[3]);
      const Real domain_y_end = stof (tokens[4]);
      const Real width = stof (tokens[5]);
      const Real height = stof (tokens[6]);
      const Real origin_x = stof (tokens[7]);
      const Real origin_y = stof (tokens[8]);
      const bool log_x = (tokens[9] == "y");
      const bool log_y = (tokens[10] == "y");
      const Domain_1D domain_x (domain_x_start, domain_x_end);
      const Domain_1D domain_y (domain_y_start, domain_y_end);
      const Point_2D origin (origin_x, origin_y);
      return new Log_Transform_2D (domain_x, domain_y,
         width, height, log_x, log_y, origin);
   }
   else
   if (tokens[0] == "POLAR")
   {
      const Real x = stof (tokens[1]);
      const Real y = stof (tokens[2]);
      const Real scale = stof (tokens[3]);
      return new Polar_Transform_2D (Point_2D (x, y), scale);
   }
   else
   if (tokens[0] == "PARABOLIC")
   {
      const Real x = stof (tokens[1]);
      const Real y = stof (tokens[2]);
      return new Parabolic_Transform_2D (Point_2D (x, y));
   }
   else
   if (tokens[0] == "ELLIPTIC")
   {
      const Real x = stof (tokens[1]);
      const Real y = stof (tokens[2]);
      const Real scale = stof (tokens[3]);
      return new Elliptic_Transform_2D (Point_2D (x, y), scale);
   }
   else
   if (tokens[0] == "BIPOLAR")
   {
      const Real x = stof (tokens[1]);
      const Real y = stof (tokens[2]);
      const Real scale = stof (tokens[3]);
      return new Bipolar_Transform_2D (Point_2D (x, y), scale);
   }

   return NULL;

}

bool
Transform_2D::out_of_domain (const Real x,
                             const Real y) const
{
   return false;
}

bool
Transform_2D::is_out_of_domain (const Point_2D& point_2d) const
{
   return out_of_domain (point_2d.x, point_2d.y);
}

void
Transform_2D::transform (Real& transformed_x,
                         Real& transformed_y,
                         const Real x,
                         const Real y) const
{
   transformed_x = x;
   transformed_y = y;
}

Point_2D
Transform_2D::transform (const Real x,
                         const Real y) const
{
   Point_2D p;
   transform (p.x, p.y, x, y);
   return p;
}

Point_2D
Transform_2D::transform (const Point_2D& point_2d) const
{
   Point_2D p;
   transform (p.x, p.y, point_2d.x, point_2d.y);
   return p;
}

void
Transform_2D::reverse (Real& reversed_x,
                       Real& reversed_y,
                       const Real x,
                       const Real y) const
{
   reversed_x = x;
   reversed_y = y;
}

Point_2D
Transform_2D::reverse (const Real x,
                       const Real y) const
{
   Point_2D p;
   reverse (p.x, p.y, x, y);
   return p;
}

Point_2D
Transform_2D::reverse (const Point_2D& point_2d) const
{
   Point_2D p;
   reverse (p.x, p.y, point_2d.x, point_2d.y);
   return p;
}

void
Transform_2D::transform_uv (Real& u,
                            Real& v,
                            const Real x,
                            const Real y) const
{
}

void
Transform_2D::transform_uv (Real& u,
                            Real& v,
                            const Point_2D& point_2d) const
{
   transform_uv (u, v, point_2d.x, point_2d.y);
}

Real
Transform_2D::get_theta (const Real u,
                         const Real v,
                         const Real x,
                         const Real y) const
{
   Real uu = u, vv = v;
   transform_uv (uu, vv, x, y);
   return atan2 (vv, uu);
}

Real
Transform_2D::get_theta (const Real u,
                         const Real v,
                         const Point_2D& point_2d) const
{
   Real uu = u, vv = v;
   transform_uv (uu, vv, point_2d);
   return atan2 (v, u);
}

void
Transform_2D::cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                     const Edge& edge) const
{

   Real x, y;

   transform (x, y, edge.point_a.x, edge.point_a.y);
   cr->move_to (x, y);
   transform (x, y, edge.point_b.x, edge.point_b.y);
   cr->move_to (x, y);

}

void
Transform_2D::cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                     const Circle& circle) const
{

   const Point_2D& center = circle.get_center ();
   const Real radius = circle.get_radius ();

   Real x, y;
   transform (x, y, center.x, center.y);

   cr->move_to (x, y);
   cr->arc (x, y, radius, 0, 2*M_PI);

}

void
Transform_2D::cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                     const Ellipse& ellipse) const
{
   const Polygon polygon (ellipse);
   cairo (cr, polygon);
}

void
Transform_2D::cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                     const Simple_Polyline& simple_polyline) const
{
   if (domain_polygon_ptr != NULL)
   {
      Polyline* p_ptr = simple_polyline.clip (*domain_polygon_ptr);
      Transform_2D::cr (cr, *p_ptr);
      delete p_ptr;
   }
   else
   { 
      Transform_2D::cr (cr, simple_polyline);
   }
}

void
Transform_2D::cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                     const Simple_Polyline& simple_polyline,
                     const Polygon& clip_polygon,
                     const bool negative) const
{
   if (domain_polygon_ptr != NULL)
   {
      Polyline* p_ptr = simple_polyline.clip (*domain_polygon_ptr);
      Transform_2D::cr (cr, *p_ptr, clip_polygon, negative);
      delete p_ptr;
   }
   else
   { 
      Transform_2D::cr (cr, simple_polyline, clip_polygon, negative);
   }
}

void
Transform_2D::cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                     const Polyline& polyline) const
{
   if (domain_polygon_ptr != NULL)
   {
      Polyline* p_ptr = polyline.clip (*domain_polygon_ptr);
      Transform_2D::cr (cr, *p_ptr);
      delete p_ptr;
   }
   else
   { 
      Transform_2D::cr (cr, polyline);
   }
}

void
Transform_2D::cairo (const Cairo::RefPtr<Cairo::Context>& cr,
                     const Polygon& polygon) const
{
   if (domain_polygon_ptr != NULL)
   {
      const Polygon& dp = *domain_polygon_ptr;
      Polygon* p_ptr = Polygon::boolean_op (INTERSECTION, dp, polygon);
      Transform_2D::cr (cr, *p_ptr);
      delete p_ptr;
   }
   else
   {
      Transform_2D::cr (cr, polygon);
   }
}

Affine_Transform_1D::Affine_Transform_1D ()
   : s (1),
     o (0)
{
}

Affine_Transform_1D::Affine_Transform_1D (const Real scale,
                                          const Real offset)
   : s (scale),
     o (offset)
{
}

Affine_Transform_1D::Affine_Transform_1D (const Affine_Transform_1D& transform)
   : s (transform.s),
     o (transform.o)
{
}

bool
Affine_Transform_1D::out_of_domain (const Real x) const
{
   return false;
}

const Real&
Affine_Transform_1D::get_scale () const
{
   return s;
}

const Real&
Affine_Transform_1D::get_offset () const
{
   return o;
}

void
Affine_Transform_1D::scale (const Real scale)
{
   this->s *= scale;
}

void
Affine_Transform_1D::scale (const Real scale,
                            const Real pivot)
{
   translate (-pivot);
   this->scale (scale);
   translate (pivot);
}

void
Affine_Transform_1D::translate (const Real translation)
{
   o += translation;
}

void
Affine_Transform_1D::transform (Real& transformed,
                                const Real x) const
{
   transformed = s * x + o;
}

Real
Affine_Transform_1D::transform (const Real x) const
{
   return Transform_1D::transform (x);
}

void
Affine_Transform_1D::reverse (Real& reversed,
                              const Real x) const
{
   reversed = (x - o) / s;
}

Real
Affine_Transform_1D::reverse (const Real x) const
{
   return Transform_1D::reverse (x);
}

Moebius_Transform::Moebius_Transform ()
   : a (1, 0),
     b (0, 0),
     c (1, 0),
     d (0, 0)
{
}

Moebius_Transform::Moebius_Transform (const Moebius_Transform& mt)
   : a (mt.a),
     b (mt.b),
     c (mt.c),
     d (mt.d)
{
}

Moebius_Transform::Moebius_Transform (const Moebius_Transform& mt_a,
                                      const Moebius_Transform& mt_b)
   : a (mt_a.a * mt_b.a + mt_a.b * mt_b.c),
     b (mt_a.a * mt_b.b + mt_a.b * mt_b.d),
     c (mt_a.c * mt_b.a + mt_a.d * mt_b.c),
     d (mt_a.c * mt_b.b + mt_a.d * mt_b.d)
{
}

Moebius_Transform::Moebius_Transform (const complex<Real>& a,
                                      const complex<Real>& b,
                                      const complex<Real>& c,
                                      const complex<Real>& d)
   : a (a),
     b (b),
     c (c),
     d (d)
{
}

Moebius_Transform::Moebius_Transform (const Point_2D& point_a,
                                      const Point_2D& point_b,
                                      const Point_2D& point_c,
                                      const Point_2D& point_d)
   : a (point_a.get_complex ()),
     b (point_b.get_complex ()),
     c (point_c.get_complex ()),
     d (point_d.get_complex ())
{
}

Moebius_Transform::Moebius_Transform (const complex<Real>& z_a,
                                      const complex<Real>& z_b,
                                      const complex<Real>& z_c)
   : a (z_b - z_c),
     b (z_a * (z_c - z_b)),
     c (z_b - z_a),
     d (z_c * (z_a - z_b))
{
}

Moebius_Transform::Moebius_Transform (const complex<Real>& z_a,
                                      const complex<Real>& z_b,
                                      const complex<Real>& z_c,
                                      const complex<Real>& w_a,
                                      const complex<Real>& w_b,
                                      const complex<Real>& w_c)
{

   const Moebius_Transform f (z_a, z_b, z_c);
   const Moebius_Transform g (w_a, w_b, w_c);
   const Moebius_Transform mt (f.get_inverse (), g);

   this->a = mt.a;
   this->b = mt.b;
   this->c = mt.c;
   this->d = mt.d;

}

void
Moebius_Transform::set_identity ()
{
   a = Real (1.0);
   b = Real (0.0);
   c = Real (1.0);
   d = Real (0.0);
}

Moebius_Transform
Moebius_Transform::get_inverse () const
{
   return Moebius_Transform (d, -b, -c, a);
}

Moebius_Transform
Moebius_Transform::get_compose (const Moebius_Transform& mt) const
{
   return Moebius_Transform (*this, mt);
}

pair<Point_2D, Point_2D>
Moebius_Transform::get_fixed_points () const
{

   complex<Real> two_c = Real (2.0) * c;
   complex<Real> a_minus_d = a - d;
   complex<Real> bc4 = Real (4.0) * b * d;
   complex<Real> discriminant = (a_minus_d * a_minus_d - bc4);

   if (std::abs (c) == 0)
   {
      complex<Real> fixed = -b / a_minus_d;
      const Point_2D infinity (GSL_POSINF, GSL_POSINF);
      return make_pair (Point_2D (fixed), infinity);
   }
   else
   {

      complex<Real> sqrt_discriminant = sqrt (discriminant);
      complex<Real> fixed_a = (a_minus_d + sqrt_discriminant) / two_c;
      complex<Real> fixed_b = (a_minus_d - sqrt_discriminant) / two_c;

      return make_pair (Point_2D (fixed_a), Point_2D (fixed_b));

   }

}

void
Moebius_Transform::transform (Real& transformed_x,
                              Real& transformed_y,
                              const Real x,
                              const Real y) const
{

   const complex<Real> z (x, y);
   const complex<Real> w = (a*z + b) / (c*z + d);

   transformed_x = std::real (w);
   transformed_y = std::imag (w);

}

Point_2D
Moebius_Transform::transform (const Real x,
                              const Real y) const
{
   return Transform_2D::transform (x, y);
}

Point_2D
Moebius_Transform::transform (const Point_2D& point_2d) const
{
   return Transform_2D::transform (point_2d);
}

void
Moebius_Transform::reverse (Real& reversed_x,
                            Real& reversed_y,
                            const Real x,
                            const Real y) const
{

   const complex<Real> w (x, y);
   const complex<Real> z = (d*w - b) / (a - c*w);

   reversed_x = std::real (z);
   reversed_y = std::imag (z);

}

Point_2D
Moebius_Transform::reverse (const Real x,
                            const Real y) const
{
   return Transform_2D::reverse (x, y);
}

Point_2D
Moebius_Transform::reverse (const Point_2D& point_2d) const
{
   return Transform_2D::reverse (point_2d);
}

void
Moebius_Transform::transform_uv (Real& u,
                                 Real& v,
                                 const Real x,
                                 const Real y) const
{

   const complex<Real> w (u, v);
   const complex<Real> z (x, y);
   const complex<Real> az_plus_b = a*z + b;
   const complex<Real> cz_plus_d = c*z + d;
   const complex<Real> denominator = cz_plus_d * cz_plus_d;
   const complex<Real> M_prime = (a*cz_plus_d - az_plus_b*c) / denominator;

   const complex<Real> answer = M_prime * w;
   u = std::real (answer);
   v = std::real (answer);

/*
   const Real p_x = -d.x * x + d.y * y + b.x;
   const Real p_y = -d.y * x - d.x * y + b.y;
   const Real q_x = c.x * x - c.y * y - a.x;
   const Real q_y = c.y * x + c.x * y - a.y;

   const Real p_x_t = a.x * u - a.y * v;
   const Real p_y_t = a.y * u + a.x * v;
   const Real q_x_t = c.x * u - c.y * v;
   const Real q_y_t = c.y * u + c.x * v;

   const Real Q = q_x * q_x + q_y * q_y;
   const Real Q_t = 2 * (q_x * q_x_t + q_y + q_y_t);
   const Real QQ = Q * Q;

   const Real A = (p_x * q_x + p_y * q_y);
   const Real A_t = (p_x * q_x_t + p_x_t * q_x + p_y * q_y_t + p_y_t * q_y);

   const Real B = (p_y * q_x - p_x * q_y);
   const Real B_t = (p_y * q_x_t + p_y_t * q_x - p_x * q_y_t - p_x_t * q_y);

   u = (Q * A_t - A * Q_t) / QQ;
   v = (Q * B_t - B * Q_t) / QQ;
*/

}

Real
Moebius_Transform::get_theta (const Real u,
                              const Real v,
                              const Real x,
                              const Real y) const
{
   Real uu = u, vv = v;
   transform_uv (uu, vv, x, y);
   return atan2 (vv, uu);
}

void
Affine_Transform_2D::init ()
{
   r[0] = m[0][1] * m[1][2] - m[1][1] * m[0][2];
   r[1] = m[0][2] * m[1][0] - m[1][2] * m[0][0];
   determinant = m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

Affine_Transform_2D::Affine_Transform_2D ()
{
   set_identity ();
   init ();
}

Affine_Transform_2D::Affine_Transform_2D (const Affine_Transform_2D& transform)
{
   memcpy (m, transform.m, sizeof (Real) * 3);
   memcpy (m+1, transform.m+1, sizeof (Real) * 3);
   init ();
}

void
Affine_Transform_2D::set_identity ()
{
   m[0][0] = 1; m[0][1] = 0; m[0][2] = 0;
   m[1][0] = 0; m[1][1] = 1; m[1][2] = 0;
   init ();
}

void
Affine_Transform_2D::rotate (const Real theta)
{

   const Real c = cos (theta);
   const Real s = sin (theta);
   const Real m_00 = m[0][0], m_01 = m[0][1], m_02 = m[0][2];
   const Real m_10 = m[1][0], m_11 = m[1][1], m_12 = m[1][2];

   m[0][0] = c * m_00 - s * m_10;
   m[0][1] = c * m_01 - s * m_11;
   m[0][2] = c * m_02 - s * m_12;

   m[1][0] = s * m_00 + c * m_10;
   m[1][1] = s * m_01 + c * m_11;
   m[1][2] = s * m_02 + c * m_12;

   init ();

}

void
Affine_Transform_2D::rotate (const Real theta,
                             const Point_2D& pivot)
{
   translate (pivot.x, pivot.y);
   rotate (theta);
   translate (-pivot.x, -pivot.y);
}

void
Affine_Transform_2D::scale (const Real scale_x,
                            const Real scale_y)
{
   m[0][0] *= scale_x;
   m[0][1] *= scale_x;
   m[0][2] *= scale_x;
   m[1][0] *= scale_y;
   m[1][1] *= scale_y;
   m[1][2] *= scale_y;
   init ();
}

void
Affine_Transform_2D::scale (const Real scale_x,
                            const Real scale_y,
                            const Point_2D& pivot)
{
   translate (-pivot.x, -pivot.y);
   scale (scale_x, scale_y);
   translate (pivot.x, pivot.y);
}

void
Affine_Transform_2D::shear_x (const Real shear_x)
{

//   Real m_00 = m[0][0];
//   Real m_10 = m[1][0];
//
//   m[0][1] += m_00 * shear_x;
//   m[1][1] += m_10 * shear_x;

   m[0][0] += m[1][0] * shear_x;
   m[0][1] += m[1][1] * shear_x;
   m[0][2] += m[1][2] * shear_x;

   init ();

}

void
Affine_Transform_2D::translate (const Real translation_x,
                                const Real translation_y)
{
   m[0][2] += translation_x;
   m[1][2] += translation_y;
   init ();
}

Real
Affine_Transform_2D::get_jacobian () const
{
   return m[0][0] * m[1][1] - m[0][1] * m[1][0];
}

const Real&
Affine_Transform_2D::get_scale_x () const
{
   return m[0][0];
}

const Real&
Affine_Transform_2D::get_scale_y () const
{
   return m[1][1];
}

const Real&
Affine_Transform_2D::get_shear_x () const
{
   return m[0][1];
}

const Real&
Affine_Transform_2D::get_shear_y () const
{
   return m[1][0];
}

const Real&
Affine_Transform_2D::get_translation_x () const
{
   return m[0][2];
}

const Real&
Affine_Transform_2D::get_translation_y () const
{
   return m[1][2];
}

void
Affine_Transform_2D::transform (Real& transformed_x,
                                Real& transformed_y,
                                const Real x,
                                const Real y) const
{
   transformed_x = m[0][0] * x + m[0][1] * y + m[0][2];
   transformed_y = m[1][0] * x + m[1][1] * y + m[1][2];
}

Point_2D
Affine_Transform_2D::transform (const Real x,
                                const Real y) const
{
   return Transform_2D::transform (x, y);
}

Point_2D
Affine_Transform_2D::transform (const Point_2D& point_2d) const
{
   return Transform_2D::transform (point_2d);
}

void
Affine_Transform_2D::reverse (Real& reversed_x,
                              Real& reversed_y,
                              const Real x,
                              const Real y) const
{
   reversed_x = (m[1][1] * x - m[0][1] * y + r[0]) / determinant;
   reversed_y = (m[0][0] * y - m[1][0] * x + r[1]) / determinant;
}

Point_2D
Affine_Transform_2D::reverse (const Real x,
                              const Real y) const
{
   return Transform_2D::reverse (x, y);
}

Point_2D
Affine_Transform_2D::reverse (const Point_2D& point_2d) const
{
   return Transform_2D::reverse (point_2d);
}

void
Affine_Transform_2D::transform_uv (Real& u,
                                   Real& v,
                                   const Real x,
                                   const Real y) const
{
   const Real uu = m[0][0] * u + m[0][1] * v;
   const Real vv = m[1][0] * u + m[1][1] * v;
   u = uu;
   v = vv;
}

Real
Affine_Transform_2D::get_theta (const Real u,
                                const Real v,
                                const Real x,
                                const Real y) const
{
   return atan2 (m[1][0] * u + m[1][1] * v, m[0][0] * u + m[0][1] * v);
}

