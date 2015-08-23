//
// geometry.cc
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

#include "cairoable.h"
#include "geometry.h"
#include "transform.h"

using namespace denise;
using namespace Cairo;

Geometry_Exception::Geometry_Exception (const Dstring& description)
   : Exception ("Geometry_Exception", description)
{
}

void
Attractor::attract (Real& x,
                    Real& y) const
{
}

void
Attractor::attract (Point_2D& point) const
{
   attract (point.x, point.y);
}

Point_2D
Attractor::get_attraction (const Point_2D& point) const
{
   Real x = point.x;
   Real y = point.y;
   attract (x, y);
   return Point_2D (x, y);
}

Simple_Attractor::Simple_Attractor (const Real grid_x)
   : grid_x (grid_x),
     grid_y (grid_x)
{
}

Simple_Attractor::Simple_Attractor (const Real grid_x,
                                    const Real grid_y)
   : grid_x (grid_x),
     grid_y (grid_y)
{
}

void
Simple_Attractor::attract (Real& x,
                           Real& y) const
{
   x = round (x / grid_x) * grid_x;
   y = round (y / grid_y) * grid_y;
}

Real
Edge::get_determinant (const Real a00,
                       const Real a01,
                       const Real a10,
                       const Real a11) const
{
   return (a00 * a11) - (a01 * a10);
}

Edge::Edge (const Edge& edge)
   : point_a (edge.point_a),
     point_b (edge.point_b)
{
}

Edge::Edge (const Point_2D& point_a,
            const Point_2D& point_b)
   : point_a (point_a),
     point_b (point_b)
{
}

void
Edge::swap ()
{
   std::swap (point_a, point_b);
}

void
Edge::translate (const Real delta_x,
                 const Real delta_y)
{
   this->point_a.x += delta_x;
   this->point_a.y += delta_y;
   this->point_b.x += delta_x;
   this->point_b.y += delta_y;
}

void
Edge::translate (const Point_2D& point)
{
   translate (point.x, point.y);
}

Real
Edge::get_root_x (const Real y) const
{
   const Real xa = point_a.x, xb = point_b.x;
   const Real ya = point_a.y, yb = point_b.y;
   return (y - ya) / (yb - ya) * (xb - xa) + xa;
}

Real
Edge::get_root_y (const Real x) const
{
   const Real xa = point_a.x, xb = point_b.x;
   const Real ya = point_a.y, yb = point_b.y;
   return (x - xa) / (xb - xa) * (yb - ya) + ya;
}

Real
Edge::get_wec (const Point_2D& point,
               const Point_2D& point_a,
               const Point_2D& point_b)
{

   const Real a = point.x - point_a.x;
   const Real b = point.y - point_a.y;

   const Real c = point_b.x - point_a.x;
   const Real d = point_b.y - point_a.y;

   return (a * d + b * a);

}

bool
Edge::perturb_if_degenerate (Point_2D& point_aa,
                             Point_2D& point_ab,
                             Point_2D& point_ba,
                             Point_2D& point_bb,
                             const Real epsilon)
{

   bool perturbed = false;

   const Real A = point_ab.x - point_aa.x;
   const Real B = point_ba.x - point_bb.x;
   const Real C = point_ab.y - point_aa.y;
   const Real D = point_ba.y - point_bb.y;

   const Real determinant = A * D - B * C;

   const Real E = point_ba.x - point_aa.x;
   const Real F = point_ba.y - point_aa.y;

   const Real numerator_a = D * E - B * F;
   const Real numerator_b = A * F - C * E;

   if (determinant == 0 && numerator_a == 0)
   {

      const Real theta = atan2 (D, B) + random (-M_PI_4, M_PI_4);
      point_aa.x += epsilon * cos (theta);
      point_aa.y += epsilon * sin (theta);
      point_ab.x += epsilon * cos (theta);
      point_ab.y += epsilon * sin (theta);
      perturbed = true;
   }
   else
   {

      const Real alpha_a = numerator_a / determinant;
      const Real alpha_b = numerator_b / determinant;

      if (alpha_b >= 0 && alpha_b <= 1)
      {
         if (alpha_a == 1)
         {
            const Real theta = atan2 (C, A) + random (-M_PI_4, M_PI_4);
            point_ab.x += epsilon * cos (theta);
            point_ab.y += epsilon * sin (theta);
            perturbed = true;
         }
      }

      if (alpha_a >= 0 && alpha_a <= 1)
      {
         if (alpha_b == 1)
         {
            const Real theta = atan2 (D, B) + random (-M_PI_4, M_PI_4);
            point_bb.x += epsilon * cos (theta);
            point_bb.y += epsilon * sin (theta);
            perturbed = true;
         }
      }

   }

   return perturbed;

}

bool
Edge::intersection (Point_2D& intersection,
                    Real& alpha_a,
                    Real& alpha_b,
                    const Point_2D& point_aa,
                    const Point_2D& point_ab,
                    const Point_2D& point_ba,
                    const Point_2D& point_bb)
{

   const Real A = point_ab.x - point_aa.x;
   const Real B = point_ba.x - point_bb.x;
   const Real C = point_ab.y - point_aa.y;
   const Real D = point_ba.y - point_bb.y;

   const Real determinant = A * D - B * C;

   if (determinant == 0)
   {
      return false;
   }

   const Real E = point_ba.x - point_aa.x;
   const Real F = point_ba.y - point_aa.y;

   const Real numerator_a = D * E - B * F;
   const Real numerator_b = A * F - C * E;

   alpha_a = numerator_a / determinant;
   alpha_b = numerator_b / determinant;

   if (alpha_a < 0 || alpha_a > 1 || alpha_b < 0 || alpha_b > 1)
   {
      return false;
   }

   intersection.x = alpha_a * A + point_aa.x;
   intersection.y = alpha_a * C + point_aa.y;
   return true;

}

bool
Edge::intersection (Point_2D& intersection,
                    const Point_2D& point_aa,
                    const Point_2D& point_ab,
                    const Point_2D& point_ba,
                    const Point_2D& point_bb)
{
   Real alpha_a, alpha_b;
   return Edge::intersection (intersection, alpha_a, alpha_b, 
      point_aa, point_ab, point_ba, point_bb);
}

Point_2D
Edge::intersection_between (const Point_2D& point_aa,
                            const Point_2D& point_ab,
                            const Point_2D& point_ba,
                            const Point_2D& point_bb)
{
   Point_2D intersection;
   Edge::intersection (intersection, point_aa, point_ab, point_ba, point_bb);
   return intersection;
}

Point_2D
Edge::intersection_between (const Edge& edge_a,
                            const Edge& edge_b)
{
   return intersection_between (edge_a.point_a,
      edge_a.point_b, edge_b.point_a, edge_b.point_b);
}

Point_2D
Edge::get_intersection_with (const Edge& edge) const
{
   return intersection_between (*this, edge);
}

bool
Edge::intersects_with (const Edge& edge) const
{
   try { get_intersection_with (edge); return true; }
   catch (const Geometry_Exception& ge) { return false; }
}

Tuple
Edge::get_intersection_tuple (const Polygon& polygon) const
{

   const Real dx = point_b.x - point_a.x;
   const Real dy = point_b.y - point_a.y;
   const bool steep = fabs (dy / dx) > 1;

   Tuple t_tuple;
   const Polygon::Vertex* first_handle_ptr = polygon.get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   if (current_handle_ptr != NULL)
   {

      Point_2D p;
      Polygon::Vertex* current_ptr;
      Polygon::Vertex* next_ptr;

      // determine the intersections on the edge
      do
      {

         const Integer n = current_handle_ptr->n;

         current_ptr = (Polygon::Vertex*)current_handle_ptr;
         next_ptr = current_ptr->next_ptr;

         for (Integer i = 0; i < n; i++)
         {

            const Point_2D& point = (const Point_2D&)(*current_ptr);
            const Point_2D& next_point = (const Point_2D&)(*next_ptr);

            if (Edge::intersection (p, point_a, point_b, point, next_point))
            {
               const Real t_y = (p.y - point_a.y) / dy;
               const Real t_x = (p.x - point_a.x) / dx;
               const Real t = (steep ? t_y : t_x);
               t_tuple.push_back (t);
            }

            current_ptr = next_ptr;
            next_ptr = current_ptr->next_ptr;

         }

         current_handle_ptr = current_handle_ptr->next_handle_ptr;

      }
      while (current_handle_ptr != first_handle_ptr);

      sort (t_tuple.begin (), t_tuple.end ());

   }

   return t_tuple;

}

Real
Edge::distance_between (const Point_2D& point_a,
                        const Point_2D& point_b)
{
   const Real dx = point_b.x - point_a.x;
   const Real dy = point_b.y - point_a.y;
   return sqrt (dx*dx + dy*dy);
}

Real
Edge::distance (const Transform_2D& transform_2d,
                const Point_2D& point) const
{
   return Edge::distance (transform_2d, point_a, point_b, point);
}

Real
Edge::distance (const Point_2D& point) const
{
   return Edge::distance (point_a, point_b, point);
}

Real
Edge::distance (const Transform_2D& transform_2d,
                const Point_2D& edge_point_a,
                const Point_2D& edge_point_b,
                const Point_2D& point)
{
   const Point_2D& pa = transform_2d.transform (edge_point_a);
   const Point_2D& pb = transform_2d.transform (edge_point_b);
   return Edge::distance (pa, pb, point);
}

Real
Edge::distance (const Point_2D& edge_point_a,
                const Point_2D& edge_point_b,
                const Point_2D& point)
{

   const Real dx_ab = edge_point_b.x - edge_point_a.x;
   const Real dy_ab = edge_point_b.y - edge_point_a.y;
   const Real dx_a = point.x - edge_point_a.x;
   const Real dy_a = point.y - edge_point_a.y;
   const Real dx_b = point.x - edge_point_b.x;
   const Real dy_b = point.y - edge_point_b.y;

   const Real l2 = dx_ab * dx_ab + dy_ab * dy_ab;

   Real t = (dx_a * dx_ab + dy_a * dy_ab) / l2;
   if (gsl_isnan (t)) { return Edge::distance_between (edge_point_a, point); }

   if (t < 0) { t = 0; }
   if (t > 1) { t = 1; }

   const Real dx = edge_point_a.x + t * dx_ab - point.x;
   const Real dy = edge_point_a.y + t * dy_ab - point.y;

   return sqrt (dx * dx + dy * dy);

}

Real
Edge::length () const
{
   return distance_between (point_a, point_b);
}

bool
Edge::contains (const Point_2D& point) const
{
   return ((point.x - point_a.x) * (point.x - point_b.x) <= 0) &&
          ((point.y - point_a.y) * (point.y - point_b.y) <= 0);
}

Polyline*
Edge::clip (const Polygon& clip_polygon,
            const bool negative) const
{

   Point_2D p;
   Polyline* polyline_ptr = new Polyline ();

   bool start_segment = true;
   const Tuple& t_tuple = get_intersection_tuple (clip_polygon);

   const Real dx = point_b.x - point_a.x;
   const Real dy = point_b.y - point_a.y;

   // Determine the segments on the edge for rendering

   if (clip_polygon.contains (point_a, false) != negative)
   {
      polyline_ptr->add (point_a, start_segment);
      start_segment = !start_segment;
   }

   for (Tuple::const_iterator iterator = t_tuple.begin ();
        iterator != t_tuple.end (); iterator++)
   {

      const Real& t = *(iterator);
      p.x = point_a.x + t * dx;
      p.y = point_a.y + t * dy;

      polyline_ptr->add (p, start_segment);
      start_segment = !start_segment;

   }

   if (!start_segment)
   {
      polyline_ptr->add (point_b, start_segment);
   }

   return polyline_ptr;

}

void
Edge::cairo (const RefPtr<Context>& cr) const
{
   cr->move_to (point_a.x, point_a.y);
   cr->line_to (point_b.x, point_b.y);
}

void
Edge::cairo (const RefPtr<Context>& cr,
             const Transform_2D& transform_2d) const
{
   const Point_2D& transformed_point_a = transform_2d.transform (point_a);
   const Point_2D& transformed_point_b = transform_2d.transform (point_b);
   cr->move_to (transformed_point_a.x, transformed_point_a.y);
   cr->line_to (transformed_point_b.x, transformed_point_b.y);
}

void
Edge::cairo (const RefPtr<Context>& cr,
             const Polygon& clip_polygon,
             const bool negative) const
{

   Point_2D p;

   bool start_segment = true;
   const Tuple& t_tuple = get_intersection_tuple (clip_polygon);

   const Real dx = point_b.x - point_a.x;
   const Real dy = point_b.y - point_a.y;

   // Determine the segments on the edge for rendering

   if (clip_polygon.contains (point_a, false) != negative)
   {
      cr->move_to (point_a.x, point_a.y);
      start_segment = !start_segment;
   }

   for (Tuple::const_iterator iterator = t_tuple.begin ();
        iterator != t_tuple.end (); iterator++)
   {

      const Real& t = *(iterator);
      p.x = point_a.x + t * dx;
      p.y = point_a.y + t * dy;

      if (start_segment)
      {
         cr->move_to (p.x, p.y);
      }
      else
      {
         cr->line_to (p.x, p.y);
      }
      start_segment = !start_segment;

   }

   if (!start_segment)
   {
      cr->line_to (point_b.x, point_b.y);
   }

}

void
Circle::init (const Real A,
              const Real D,
              const Real E,
              const Real F)
{
   this->center = Point_2D (-D / (2*A), -E / (2*A));
   this->radius = sqrt (((D*D + E*E) / (4 *A*A)) - (F / A));
}
         
Circle::Circle ()
{
}
         
Circle::Circle (const Point_2D& center,
                const Real radius)
      : center (center),
        radius (radius)
{
}
         
Circle::Circle (const Real A,
                const Real D,
                const Real E,
                const Real F)
{
   init (A, D, E, F);
}

Circle::Circle (const Point_2D& point_a,
                const Point_2D& point_b,
                const Point_2D& point_c)
{

   const Real r_a = point_a.x*point_a.x + point_a.y*point_a.y;
   const Real r_b = point_b.x*point_b.x + point_b.y*point_b.y;
   const Real r_c = point_c.x*point_c.x + point_c.y*point_c.y;

   const Real A = determinant (point_a.x, point_a.y, 1,
      point_b.x, point_b.y, 1, point_c.x, point_c.y, 1);

   const Real D = - determinant (r_a, point_a.y, 1,
      r_b, point_b.y, 1, r_c, point_c.y, 1);

   const Real E = determinant (r_a, point_a.x, 1,
      r_b, point_b.x, 1, r_c, point_c.x, 1);

   const Real F = - determinant (r_a, point_a.x, point_a.y,
      r_b, point_b.x, point_b.y, r_c, point_c.x, point_c.y);

   init (A, D, E, F);

}

void
Circle::acquire_coefficients (Real& A,
                              Real& D,
                              Real& E,
                              Real& F) const
{
   A = 1;
   D = -2 * center.x;
   E = -2 * center.y;
   F = center.x*center.x + center.y*center.y - radius*radius;
}
                  
void     
Circle::translate (const Real delta_x,
                   const Real delta_y)
{
   this->center.x += delta_x;
   this->center.y += delta_y;
}
                  
void     
Circle::translate (const Point_2D& point)
{
   translate (point.x, point.y);
}

void     
Circle::scale (const Real scale)
{
   this->radius *= scale;
}

bool
Circle::contains (const Point_2D& point,
                  const bool border_included) const
{

   const Real dx = point.x - center.x;
   const Real dy = point.y - center.y;
   const Real r_square = dx*dx + dy*dy;
   const Real radius_square = radius*radius;
   if (r_square < radius_square) { return true; }
   if (r_square == radius_square && border_included) { return true; }

   return false;

}
                  
const Point_2D&
Circle::get_center () const
{
   return center;
}

Real
Circle::get_radius () const
{
   return radius;
}

Real
Circle::get_area () const
{
   return M_PI * radius * radius;
}

Real
Circle::get_perimeter () const
{
   return M_2_TIMES_PI * radius;
}

Point_2D
Circle::get_point_at_theta (const Real theta) const
{
   return Point_2D (center.x + radius * cos (theta),
                    center.y + radius * sin (theta));
}

pair<Point_2D, Point_2D>
Circle::get_tangent_point_pair (const Point_2D& point) const
{

   const Real dx = point.x - center.x;
   const Real dy = point.y - center.y;

   const Real theta = atan2 (dy, dx);
   const Real distance = sqrt (dx*dx + dy*dy);

   const Real t = acos (radius / distance);
   const Real t_0 = theta + t;
   const Real t_1 = theta - t;

   return make_pair (get_point_at_theta (t_0), get_point_at_theta (t_1));

}

bool
Circle::is_nac () const
{
   return (center.is_nap () || gsl_isnan (radius));
}

void
Circle::cairo (const RefPtr<Context>& cr) const
{
   cr->move_to (center.x, center.y);
   cr->arc (center.x, center.y, radius, 0, 2*M_PI);
}

void
Circle::cairo (const RefPtr<Context>& cr,
               const Transform_2D& transform_2d) const
{
   const Point_2D& transformed_center = transform_2d.transform (center);
   cr->move_to (transformed_center.x, transformed_center.y);
   cr->arc (transformed_center.x, transformed_center.y,
      radius, 0, 2*M_PI);
}

void
Ellipse::init (const Real A,
               const Real B,
               const Real C,
               const Real D,
               const Real E,
               const Real F)
{

   const Real f = 2 / (4*A*C - B*B);
   const Real d = (C*D*D - B*D*E + A*E*E) / (4*A*C - B*B) - F;

   this->center = Point_2D ((B*E/2 - C*D) * f, (B*D/2 - A*E) * f);

   gsl_matrix* M = gsl_matrix_alloc (2, 2);
   gsl_matrix_set (M, 0, 0, A / d);
   gsl_matrix_set (M, 0, 1, B / (2*d));
   gsl_matrix_set (M, 1, 0, B / (2*d));
   gsl_matrix_set (M, 1, 1, C / d);

   gsl_vector* eigenvalues = gsl_vector_alloc (2);
   gsl_matrix* eigenvectors = gsl_matrix_alloc (2, 2);

   gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc (2);
   gsl_eigen_symmv (M, eigenvalues, eigenvectors, w);
   gsl_eigen_symmv_free (w);

   this->a = 1.0 / sqrt (gsl_vector_get (eigenvalues, 0));
   this->b = 1.0 / sqrt (gsl_vector_get (eigenvalues, 1));
   this->tilt = asin (gsl_matrix_get (eigenvectors, 0, 1));

   if (a < b)
   {
      swap (a, b);
      this->tilt -= M_PI_2;
   }

   gsl_matrix_free (M);
   gsl_vector_free (eigenvalues);
   gsl_matrix_free (eigenvectors);

}

Ellipse::Ellipse ()
{
   this->center = Point_2D (GSL_NAN, GSL_NAN);
   this->a = GSL_NAN;
   this->b = GSL_NAN;
   this->tilt = GSL_NAN;
}

Ellipse::Ellipse (const Point_2D& center,
                  const Real major_radius,
                  const Real minor_radius,
                  const Real tilt)
        : center (center),
               a (major_radius),
               b (minor_radius),
            tilt (tilt)
{
//   this->tilt = fmod ((tilt + M_PI), M_2_TIMES_PI) - M_PI;
   this->tilt = fmodf ((tilt + M_PI), M_2_TIMES_PI) - M_PI;
}

Ellipse::Ellipse (const Real A,
                  const Real B,
                  const Real C,
                  const Real F)
{
   init (A, B, C, 0, 0, F);
}

Ellipse::Ellipse (const Real A,
                  const Real B,
                  const Real C,
                  const Real D,
                  const Real E,
                  const Real F)
{
   init (A, B, C, D, E, F);
}

void
Ellipse::acquire_coefficients (Real& A,
                               Real& B,
                               Real& C,
                               Real& D,
                               Real& E,
                               Real& F) const
{

   const Real s = sin (tilt);
   const Real c = cos (tilt);

   const Real a2 = a*a;
   const Real b2 = b*b;
   const Real c2 = c*c;
   const Real s2 = s*s;
   const Real sc = s*c;
   const Real x0 = center.x;
   const Real y0 = center.y;

   A = a2*s2 + b2*c2;
   B = 2*sc* (a2-b2);
   C = a2*c2 + b2*s2;
   D = 2*(y0*b2*sc - y0*a2*sc - x0*b2*c2 - x0*a2*s2);
   E = 2*(x0*b2*sc - x0*a2*sc - y0*a2*c2 - y0*b2*s2);
   F = x0*x0*b2*c2 + x0*x0*a2*s2 + y0*y0*a2*c2 + y0*y0*b2*s2 +
       2*x0*y0*a2*sc - 2*x0*y0*b2*sc - a2*b2;

}

void
Ellipse::translate (const Real delta_x,
                    const Real delta_y)
{
   this->center.x += delta_x;
   this->center.y += delta_y;
}

void
Ellipse::translate (const Point_2D& point)
{
   translate (point.x, point.y);
}

void
Ellipse::scale (const Real scale)
{
   this->a *= scale;
   this->b *= scale;
}

void
Ellipse::rotate (const Real tilt)
{
   this->tilt += tilt;
}

Point_2D
Ellipse::get_center () const
{
   return center;
}
         
Real
Ellipse::get_a () const
{
   return a;
}
         
Real
Ellipse::get_b () const
{
   return b;
}
   
Real
Ellipse::get_tilt () const
{
   return tilt;
}

Real
Ellipse::get_area () const
{
   return M_PI * a * b;
}

Real
Ellipse::get_perimeter (const gsl_mode_t gsl_mode) const
{
   const Real k = sqrt (1 - (a*a / b*b));
   return 4 * b * gsl_sf_ellint_Ecomp (k, gsl_mode);
/*
   Real gauss_kummer[] = { 1, 4, 64, 256, 16384, 65536, 1048576, 4194304 };

   Real h = (a - b) / (a + b);
   h *= h;

   Real perimeter = 1;
   for (Integer i = 0; i < 8; i++)
   {
      perimeter += 1.0 / gauss_kummer[i] * pow (h, i);
   }

   perimeter *= M_PI * (a + b);
   return perimeter;
*/
}

Point_2D
Ellipse::get_point_at_t (const Real t) const
{

   const Real c = cos (tilt);
   const Real s = sin (tilt);

   const Real A = a * c;
   const Real B = a * s;
   const Real C = b * c;
   const Real D = b * s;

   const Real radius_x = A * cos (t) - D * sin (t);
   const Real radius_y = B * cos (t) + C * sin (t);

   return Point_2D (center.x + radius_x, center.y + radius_y);

}

Point_2D
Ellipse::get_point_at_theta (const Real theta) const
{
   const Real t = atan2 (a * sin (theta - tilt), b * cos (theta - tilt));
   return get_point_at_t (t);
}

pair<Point_2D, Point_2D>
Ellipse::get_tangent_point_pair (const Point_2D& point) const
{

   const Real c = cos (tilt);
   const Real s = sin (tilt);

   const Real dx = point.x - center.x;
   const Real dy = point.y - center.y;

   const Real A = a * c;
   const Real B = -b * s;
   const Real C = a * s;
   const Real D = b * c;

   const Real psi = dy * B - dx * D;
   const Real phi = dx * C - dy * A;

   const Real lambda = atan2 (phi, psi);
   const Real amplitude = sqrt (psi*psi + phi*phi);

   const Real BCmAD = B*C - A*D;

   const Real K = acos (BCmAD / amplitude);

   const Real t_0 = K + lambda;
   const Real t_1 = M_2_TIMES_PI - K + lambda;

   const Point_2D point_0 = get_point_at_t (t_0);
   const Point_2D point_1 = get_point_at_t (t_1);

   Real angle_span = atan2 (point_0.y - point.y, point_0.x - point.x) -
                     atan2 (point_1.y - point.y, point_1.x - point.x);
   angle_span = fmod ((angle_span + M_4_TIMES_PI), M_2_TIMES_PI);
   if (angle_span > M_PI) { angle_span -= M_2_TIMES_PI; }

   if (angle_span > 0) { return make_pair (point_1, point_0); }
   else                { return make_pair (point_0, point_1); }

}

bool
Ellipse::is_nae () const
{
   return center.is_nap () || gsl_isnan (a) ||
      gsl_isnan (b) || gsl_isnan (tilt);
}

void
Ellipse::cairo (const RefPtr<Context>& cr) const
{
   const Polygon polygon (*this);
   polygon.cairo (cr);
}

void
Ellipse::cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d) const
{
   const Polygon polygon (*this);
   polygon.cairo (cr, transform_2d);
}

bool
Simple_Polyline::is_last (Simple_Polyline::iterator iterator) const
{
   return (++iterator) == end ();
}

bool
Simple_Polyline::is_last (Simple_Polyline::const_iterator iterator) const
{
   return (++iterator) == end ();
}

Simple_Polyline::Simple_Polyline (const bool closed)
   : closed (closed)
{
}

Simple_Polyline::Simple_Polyline (const Simple_Polyline& simple_polyline)
{
   set (simple_polyline);
}

void
Simple_Polyline::set (const Simple_Polyline& simple_polyline)
{

   clear ();
   this->closed = simple_polyline.closed;

   for (Simple_Polyline::const_iterator iterator = simple_polyline.begin ();
        iterator != simple_polyline.end (); iterator++)
   {
      push_back (*(iterator));
   }

}

void
Simple_Polyline::attract_by (const Attractor& attractor)
{
   for (Simple_Polyline::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Point_2D& p = *(iterator);
      attractor.attract (p);
   }
}

bool
Simple_Polyline::is_closed () const
{
   return closed;
}

void
Simple_Polyline::set_closed (const bool closed)
{
   this->closed = closed;
}

Edge
Simple_Polyline::get_edge (Simple_Polyline::iterator iterator) const
{
   if (is_last (iterator))
   {
      if (closed) { return Edge (*(iterator), *(begin ())); }
      else { throw Exception ("Last Node in Simple_Polyline"); }
   }
   else
   {
      Simple_Polyline::iterator next = iterator;
      return Edge (*(iterator), *(++next));
   }
}

Edge
Simple_Polyline::get_edge (Simple_Polyline::const_iterator iterator) const
{
   if (is_last (iterator))
   {
      if (closed) { return Edge (*(iterator), *(begin ())); }
      else { throw Exception ("Last Node in Simple_Polyline"); }
   }
   else
   {
      Simple_Polyline::const_iterator next = iterator;
      return Edge (*(iterator), *(++next));
   }
}

void
Simple_Polyline::add (const Point_2D& point)
{
   push_back (point);
}

void
Simple_Polyline::prepend (const Point_2D& point)
{
   push_front (point);
}

void
Simple_Polyline::translate (const Real dx,
                            const Real dy)
{
   for (Simple_Polyline::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Point_2D& p = *(iterator);
      p.x += dx;
      p.y += dy;
   }
}

void
Simple_Polyline::transform (const Transform_2D& transform)
{
   for (Simple_Polyline::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Point_2D& point = *(iterator);
      transform.transform (point.x, point.y, point.x, point.y);
   }
}

void
Simple_Polyline::reverse (const Transform_2D& transform)
{
   for (Simple_Polyline::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Point_2D& point = *(iterator);
      transform.reverse (point.x, point.y, point.x, point.y);
   }
}

Polyline*
Simple_Polyline::clip (const Polygon& clip_polygon,
                       const bool negative) const
{
   Polyline* polyline_ptr = new Polyline ();
   clip (*polyline_ptr, clip_polygon, negative);
   return polyline_ptr;
}

void
Simple_Polyline::clip (Polyline& polyline,
                       const Polygon& clip_polygon,
                       const bool negative) const
{


   Point_2D p;
   Edge edge;

   bool start_segment = true;
   const Point_2D& first_point = *(begin ());

   if (clip_polygon.contains (first_point, false) != negative)
   {
      polyline.add (first_point, start_segment);
      start_segment = !start_segment;
   }

   edge.point_b.x = first_point.x;
   edge.point_b.y = first_point.y;

   for (Simple_Polyline::const_iterator this_iterator = begin ();
        this_iterator != end (); this_iterator++)
   {

      Simple_Polyline::const_iterator next_iterator = this_iterator;
      next_iterator++;

      if (next_iterator == end ())
      {
         break;
      }

      edge.point_a.x = edge.point_b.x;
      edge.point_a.y = edge.point_b.y;
      edge.point_b.x = next_iterator->x;
      edge.point_b.y = next_iterator->y;

      const Tuple& t_tuple = edge.get_intersection_tuple (clip_polygon);

      const Real dx = edge.point_b.x - edge.point_a.x;
      const Real dy = edge.point_b.y - edge.point_a.y;

      // Determine the segments on the edge for rendering

      for (Tuple::const_iterator iterator = t_tuple.begin ();
           iterator != t_tuple.end (); iterator++)
      {

         const Real& t = *(iterator);
         p.x = edge.point_a.x + t * dx;
         p.y = edge.point_a.y + t * dy;

         polyline.add (p, start_segment);
         start_segment = !start_segment;

      }

      if (clip_polygon.contains (edge.point_b, false) != negative)
      {
         polyline.add (edge.point_b, start_segment);
      }

   }

}

void
Simple_Polyline::cairo (const RefPtr<Context>& cr) const
{

   for (Simple_Polyline::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Point_2D& point = *(iterator);

      if (iterator == begin ())
      {
         cr->move_to (point.x, point.y);
      }
      else
      {
         cr->line_to (point.x, point.y);
      }

   }

}

void
Simple_Polyline::cairo (const RefPtr<Context>& cr,
                        const Transform_2D& transform_2d) const
{

   for (Simple_Polyline::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Point_2D& point = *(iterator);
      const Point_2D& t_point = transform_2d.transform (point.x, point.y);

      if (iterator == begin ())
      {
         cr->move_to (t_point.x, t_point.y);
      }
      else
      {
         cr->line_to (t_point.x, t_point.y);
      }

   }
 
}

void
Simple_Polyline::cairo (const RefPtr<Context>& cr,
                        const Polygon& clip_polygon,
                        const bool negative) const
{

   if (size () < 2) { return; }

   Edge edge;

   for (Simple_Polyline::const_iterator this_iterator = begin ();
        this_iterator != end (); this_iterator++)
   {

      Simple_Polyline::const_iterator next_iterator = this_iterator;
      if ((++next_iterator) == end ()) { break; }

      const Point_2D& this_point = *(this_iterator);
      const Point_2D& next_point = *(next_iterator);

      edge.point_a.x = this_iterator->x;
      edge.point_a.y = this_iterator->y;
      edge.point_b.x = next_iterator->x;
      edge.point_b.y = next_iterator->y;

//         Polyline* p_ptr = edge.clip (clip_polygon, true);
//         p_ptr->cairo (cr);
//         delete p_ptr;


      edge.cairo (cr, clip_polygon, negative);
//      edge.cairo (cr);

   }

}

Simple_Polyline::iterator
Simple_Polyline::get_iterator (const Transform_2D& transform_2d,
                               const Point_2D& point_2d,
                               const Real threshold)
{

   for (Simple_Polyline::iterator i = begin (); i != end (); i++)
   {
      if (!closed && is_last (i)) { continue; }
      Edge edge = get_edge (i);
      if (edge.distance (transform_2d, point_2d) < threshold) { return i; }
   }

   return end ();

}

Simple_Polyline::const_iterator
Simple_Polyline::get_iterator (const Transform_2D& transform_2d,
                               const Point_2D& point_2d,
                               const Real threshold) const
{

   for (Simple_Polyline::const_iterator i = begin (); i != end (); i++)
   {
      if (!closed && is_last (i)) { continue; }
      Edge edge = get_edge (i);
      if (edge.distance (transform_2d, point_2d) < threshold) { return i; }
   }

   return end ();

}

Simple_Polyline::iterator
Simple_Polyline::get_iterator (const Point_2D& point_2d,
                               const Real threshold)
{

   for (Simple_Polyline::iterator i = begin (); i != end (); i++)
   {
      if (!closed && is_last (i)) { continue; }
      Edge edge = get_edge (i);
      if (edge.distance (point_2d) < threshold) { return i; }
   }

   return end ();

}

Simple_Polyline::const_iterator
Simple_Polyline::get_iterator (const Point_2D& point_2d,
                               const Real threshold) const
{

   for (Simple_Polyline::const_iterator i = begin (); i != end (); i++)
   {
      if (!closed && is_last (i)) { continue; }
      Edge edge = get_edge (i);
      if (edge.distance (point_2d) < threshold) { return i; }
   }

   return end ();

}

Simple_Polyline::iterator
Simple_Polyline::implant (const Transform_2D& transform_2d,
                          const Point_2D& point_2d,
                          const Real threshold)
{

   typedef Simple_Polyline::iterator Iterator;
   Iterator iterator = get_iterator (transform_2d, point_2d, threshold);
   if (iterator == end ()) { return iterator; }

   return insert (++iterator, transform_2d.reverse (point_2d));

}

Simple_Polyline::iterator
Simple_Polyline::implant (const Point_2D& point_2d,
                          const Real threshold)
{

   Simple_Polyline::iterator iterator = get_iterator (point_2d, threshold);
   if (iterator == end ()) { return iterator; }

   return insert (++iterator, point_2d);

}

Polygon::Vertex::Vertex (const Point_2D& point)
   : Point_2D (point),
     action (NOT_DECIDED),
     n (0),
     alpha (GSL_NAN),
     intersections (0),
     next_ptr (NULL),
     prev_ptr (NULL),
     handle_ptr (NULL),
     next_handle_ptr (NULL),
     prev_handle_ptr (NULL),
     neighbor_ptr (NULL),
     intersection (INTERSECTION_NOT)
{
}

void
Polyline::init (const Point_2D& point)
{
   Simple_Polyline* simple_polyline_ptr = new Simple_Polyline ();
   simple_polyline_ptr->push_back (point);
   push_back (simple_polyline_ptr);
}

Polyline::Polyline ()
{
   for (Polyline::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Simple_Polyline* simple_polyline_ptr = *(iterator);
      delete simple_polyline_ptr;
   }
}

Polyline::Polyline (const Simple_Polyline* simple_polyline_ptr)
{
   add (simple_polyline_ptr);
}

void
Polyline::attract_by (const Attractor& attractor)
{
   for (Polyline::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Simple_Polyline& simple_polyline = **(iterator);
      simple_polyline.attract_by (attractor);
   }
}

void
Polyline::add (const Simple_Polyline* simple_polyline_ptr)
{
   push_back ((Simple_Polyline*)simple_polyline_ptr);
}

void
Polyline::add (const Point_2D& point,
               const bool new_handle)
{
   if (size () == 0 || new_handle) { init (point); }
   else { back ()->push_back (point); }
}

void
Polyline::prepend (const Point_2D& point)
{
   if (size () == 0) { init (point); }
   else { front ()->push_front (point); }
}

void
Polyline::transform (const Transform_2D& transform)
{
   for (Polyline::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Simple_Polyline& simple_polyline = **(iterator);
      simple_polyline.transform (transform);
   }
}

void
Polyline::reverse (const Transform_2D& transform)
{
   for (Polyline::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Simple_Polyline& simple_polyline = **(iterator);
      simple_polyline.reverse (transform);
   }
}

Polyline*
Polyline::clip (const Polygon& clip_polygon,
                const bool negative) const
{

   Polyline* polyline_ptr = new Polyline ();

   for (Polyline::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Simple_Polyline& simple_polyline = **(iterator);
      simple_polyline.clip (*polyline_ptr, clip_polygon, negative);
   }

   return polyline_ptr;

}

void
Polyline::cairo (const RefPtr<Context>& cr) const
{
   for (Polyline::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Simple_Polyline& simple_polyline = **(iterator);
      simple_polyline.cairo (cr);
   }
}

void
Polyline::cairo (const RefPtr<Context>& cr,
                 const Transform_2D& transform_2d) const
{
   for (Polyline::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Simple_Polyline& simple_polyline = **(iterator);
      simple_polyline.cairo (cr, transform_2d);
   }
}

inline
Real
Polygon::determinant (const Point_2D& point,
                      const Point_2D& point_i,
                      const Point_2D& point_ip1)
{
   return (point_i.x - point.x) * (point_ip1.y - point.y) -
          (point_ip1.x - point.x) * (point_i.y - point.y);
}

inline
bool
Polygon::cross (const Point_2D& point,
                const Point_2D& point_i,
                const Point_2D& point_ip1)
{
   return ((point_i.y < point.y) != (point_ip1.y < point.y));
}

inline
bool
Polygon::cross_r (const Point_2D& point,
                  const Point_2D& point_i,
                  const Point_2D& point_ip1)
{
   const Real det = determinant (point, point_i, point_ip1);
   return ((det > 0) == (point_ip1.y > point_i.y));
}

inline
void
Polygon::modify (Integer& winding_number,
                 const Point_2D& point_i,
                 const Point_2D& point_ip1)
{
   winding_number += (2 * (Integer (point_ip1.y > point_i.y)) - 1);
}

void
Polygon::iterate_winding_number (Integer& w,
                                 const Polygon::Vertex* handle_ptr,
                                 const Point_2D& point)
{

   const Geometry_Exception gev ("point is vertex");
   const Geometry_Exception gee ("point on edge");

   const Integer n = handle_ptr->n;
   Polygon::Vertex* current_ptr = (Polygon::Vertex*)handle_ptr;
   if (((const Point_2D&)(*current_ptr)) == point) { throw gev; }

   for (Integer i = 0; i < n; i++)
   {

      const Polygon::Vertex* next_ptr = current_ptr->next_ptr;
      const Point_2D& p = (const Point_2D&)(*current_ptr);
      const Point_2D& np = (const Point_2D&)(*next_ptr);

      if (np.y == point.y)
      {

         if (np.x == point.x)
         {
            throw gev;
         }

         if ((p.y == point.y) && ((np.x > point.x) == (p.x < point.x)))
         {
            throw gee;
         }

      }

      if (cross (point, p, np))
      {

         if (p.x >= point.x)
         {
            if (np.x > point.x) { modify (w, p, np); }
            else if (cross_r (point, p, np)) { modify (w, p, np); }
         }
         else
         if ((np.x > point.x) && cross_r (point, p, np))
         {
            modify (w, p, np);
         }

      }

      current_ptr = current_ptr->next_ptr;

   }

}

Integer
Polygon::get_winding_number (const Polygon::Vertex* handle_ptr,
                             const Point_2D& point)
{
   Integer w = 0;
   iterate_winding_number (w, handle_ptr, point);
   return w;
}

bool
Polygon::contains (const Polygon::Vertex* handle_ptr,
                   const Point_2D& point,
                   bool border_included)
{
   try { return abs (get_winding_number (handle_ptr, point) % 2); }
   catch (const Geometry_Exception& ge) { return border_included; }
}

inline
void
Polygon::toggle_intersection (Polygon::Intersection& i)
{
        if (i == INTERSECTION_ENTRY) { i = INTERSECTION_EXIT; }
   else if (i == INTERSECTION_EXIT) { i = INTERSECTION_ENTRY; }
}

Polygon::Vertex*
Polygon::get_intersection_ptr () const
{

   Polygon::Vertex* intersection_ptr = NULL;
   Polygon::Vertex* current_handle_ptr = first_handle_ptr;

   do
   {

      Polygon::Vertex* current_ptr = current_handle_ptr;

      for (Integer i = 0; i < current_handle_ptr->n; i++)
      {
         if (current_ptr->neighbor_ptr != NULL)
         {
            intersection_ptr = current_ptr;
            break;
         }
         current_ptr = current_ptr->next_ptr;
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr && intersection_ptr == NULL);

   return intersection_ptr;

}

Point_2D
Polygon::get_intersection (Point_2D& point_aa,
                           Point_2D& point_ab,
                           Point_2D& point_ba,
                           Point_2D& point_bb,
                           Real perturbation)
{
   Point_2D ip;

   try
   {
      ip = Edge::intersection_between (point_aa, point_ab, point_ba, point_bb);
   }
   catch (const Geometry_Exception& ge)
   {
      point_ba.perturb (perturbation);
      point_bb.perturb (perturbation);
      return Polygon::get_intersection (point_aa, point_ab,
         point_ba, point_bb, perturbation);
   }

   return ip;

}

void
Polygon::boolean_op_phase_0 (const Polygon& polygon_a,
                             const Polygon& polygon_b,
                             const Real epsilon)
{

   typedef Polygon::Vertex Pv;

   const Pv* first_handle_b_ptr = polygon_b.first_handle_ptr;
   Pv* current_handle_b_ptr = (Polygon::Vertex*)first_handle_b_ptr;

   do
   {

      Pv* current_b_ptr = current_handle_b_ptr;

      do
      {

         const Pv* first_handle_a_ptr = polygon_a.first_handle_ptr;
         Pv* current_handle_a_ptr = (Pv*)first_handle_a_ptr;

         do
         {

            Pv* current_a_ptr = current_handle_a_ptr;

            do
            {

               Polygon::Vertex* next_a_ptr = current_a_ptr->next_ptr;
               Polygon::Vertex* next_b_ptr = current_b_ptr->next_ptr;

               Point_2D& pa = ((Point_2D&)(*current_a_ptr));
               Point_2D& npa = ((Point_2D&)(*next_a_ptr));
               Point_2D& pb = ((Point_2D&)(*current_b_ptr));
               Point_2D& npb = ((Point_2D&)(*next_b_ptr));

               Edge::perturb_if_degenerate (pa, npa, pb, npb, epsilon);

               current_a_ptr = current_a_ptr->next_ptr;
            }
            while (current_a_ptr != current_handle_a_ptr);

            current_handle_a_ptr = current_handle_a_ptr->next_handle_ptr;
         }
         while (current_handle_a_ptr != first_handle_a_ptr);

         current_b_ptr = current_b_ptr->next_ptr;
      }
      while (current_b_ptr != current_handle_b_ptr);

      current_handle_b_ptr = current_handle_b_ptr->next_handle_ptr;
   }
   while (current_handle_b_ptr != first_handle_b_ptr);

}

void
Polygon::boolean_op_phase_a (const Polygon& polygon_a,
                             const Polygon& polygon_b)
{

   Point_2D ip;
   Real alpha_a, alpha_b;

   typedef Polygon::Vertex Pv;
   typedef Polygon::Vertex* Pv_Ptr;
   typedef list<Pv_Ptr> Vertex_Ptr_List;

   Vertex_Ptr_List vertex_ptr_list;

   const Pv* first_handle_a_ptr = polygon_a.first_handle_ptr;
   Pv* current_handle_a_ptr = (Polygon::Vertex*)first_handle_a_ptr;

   do
   {

      Pv* current_a_ptr = current_handle_a_ptr;

      do
      {

         const Pv* first_handle_b_ptr = polygon_b.first_handle_ptr;
         Pv* current_handle_b_ptr = (Pv*)first_handle_b_ptr;

         do
         {

            Pv* current_b_ptr = current_handle_b_ptr;

            do
            {

               Polygon::Vertex* next_a_ptr = current_a_ptr->next_ptr;
               Polygon::Vertex* next_b_ptr = current_b_ptr->next_ptr;

               Point_2D& pa = ((Point_2D&)(*current_a_ptr));
               Point_2D& pb = ((Point_2D&)(*current_b_ptr));
               Point_2D& npa = ((Point_2D&)(*next_a_ptr));
               Point_2D& npb = ((Point_2D&)(*next_b_ptr));

               if (Edge::intersection (ip, alpha_a, alpha_b, pa, npa, pb, npb))
               {

                  Polygon::Vertex* iva_ptr = new Polygon::Vertex (ip);
                  Polygon::Vertex* ivb_ptr = new Polygon::Vertex (ip);

                  iva_ptr->neighbor_ptr = ivb_ptr;
                  iva_ptr->handle_ptr = current_handle_a_ptr;
                  iva_ptr->prev_ptr = current_a_ptr;
                  iva_ptr->next_ptr = next_a_ptr;
                  iva_ptr->alpha = alpha_a;

                  ivb_ptr->neighbor_ptr = iva_ptr;
                  ivb_ptr->handle_ptr = current_handle_b_ptr;
                  ivb_ptr->prev_ptr = current_b_ptr;
                  ivb_ptr->next_ptr = next_b_ptr;
                  ivb_ptr->alpha = alpha_b;

                  vertex_ptr_list.push_back (iva_ptr);
                  vertex_ptr_list.push_back (ivb_ptr);

               }


               current_b_ptr = current_b_ptr->next_ptr;
            }
            while (current_b_ptr != current_handle_b_ptr);

            current_handle_b_ptr = current_handle_b_ptr->next_handle_ptr;
         }
         while (current_handle_b_ptr != first_handle_b_ptr);

         current_a_ptr = current_a_ptr->next_ptr;
      }
      while (current_a_ptr != current_handle_a_ptr);

      current_handle_a_ptr = current_handle_a_ptr->next_handle_ptr;
   }
   while (current_handle_a_ptr != first_handle_a_ptr);



   for (Vertex_Ptr_List::iterator iterator = vertex_ptr_list.begin ();
        iterator != vertex_ptr_list.end (); iterator++)
   {

      Pv_Ptr iv_ptr = *(iterator);
      Real& alpha = iv_ptr->alpha;
      Point_2D& point = (Point_2D&)(*iv_ptr);

      Pv_Ptr edge_head_ptr = iv_ptr->prev_ptr;
      Pv_Ptr edge_tail_ptr = iv_ptr->next_ptr;
      Point_2D& prev_point = (Point_2D&)(*edge_head_ptr);
      Point_2D& next_point = (Point_2D&)(*edge_tail_ptr);

      Pv_Ptr prev_ptr = edge_head_ptr;
      Pv_Ptr next_ptr = edge_tail_ptr;

      if (next_ptr != prev_ptr->next_ptr)
      {

         edge_head_ptr->alpha = 0;
         edge_tail_ptr->alpha = 1;

         while (prev_ptr->next_ptr->alpha <= alpha)
         {
            prev_ptr = prev_ptr->next_ptr;
         }

         while (next_ptr->prev_ptr->alpha >= alpha)
         {
            next_ptr = next_ptr->prev_ptr;
         }

         if (prev_ptr == next_ptr)
         {

            const Real alpha_p = prev_ptr->alpha;
            const Real alpha_0 = next_ptr->alpha;
            const Real alpha_1 = next_ptr->next_ptr->alpha;
            const Real dx = edge_tail_ptr->x - edge_head_ptr->x;

            const Real dy = edge_tail_ptr->y - edge_head_ptr->y;
            alpha = (alpha_1 - alpha_0) * 1e-3 + alpha_0;

            iv_ptr->x = edge_head_ptr->x + alpha * dx;
            iv_ptr->y = edge_head_ptr->y + alpha * dy;
            next_ptr = next_ptr->next_ptr;
         }

         edge_head_ptr->alpha = GSL_NAN;
         edge_tail_ptr->alpha = GSL_NAN;

      }

      iv_ptr->next_ptr = next_ptr;
      iv_ptr->prev_ptr = prev_ptr;
      prev_ptr->next_ptr = iv_ptr;
      next_ptr->prev_ptr = iv_ptr;
      (iv_ptr->handle_ptr->n)++;
      (iv_ptr->handle_ptr->intersections)++;
      
   }


}

void
Polygon::boolean_op_phase_b (const Polygon& polygon_a,
                             const Polygon& polygon_b,
                             const bool flip)
{

   Polygon::Intersection intersection;

   const Polygon::Vertex* first_handle_ptr = polygon_a.first_handle_ptr;
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {

      Polygon::Vertex* current_ptr = current_handle_ptr;

      for (Integer i = 0; i < current_handle_ptr->n; i++)
      {

         const Point_2D& p = ((const Point_2D&)(*current_ptr));

         if (i == 0)
         {
            bool contains = polygon_b.contains (p);
            if (flip) { contains = !contains; }
            intersection = (contains ? INTERSECTION_EXIT: INTERSECTION_ENTRY);
         }

         if (current_ptr->neighbor_ptr != NULL)
         {
            current_ptr->intersection = intersection;
            toggle_intersection (intersection);
         }

         current_ptr = current_ptr->next_ptr;

      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);

}

void
Polygon::boolean_op_phase_c (Polygon& polygon,
                             const Polygon& polygon_a,
                             const Polygon& polygon_b,
                             const Boolean_Op boolean_op)
{

   Polygon::Vertex* first_intersection_ptr;

   vector<Polygon::Vertex*> remove_vector;

   while ((first_intersection_ptr = polygon_a.get_intersection_ptr ()) != NULL)
   {

      Polygon::Vertex* current_ptr = first_intersection_ptr;

      const Point_2D& point = ((const Point_2D&)(*current_ptr));
      polygon.add (point, true);

      do
      {

         Polygon::Intersection intersection = current_ptr->intersection;

         do
         {

            switch (intersection)
            {
               case INTERSECTION_ENTRY:
                  current_ptr = current_ptr->next_ptr;
                  break;

               case INTERSECTION_EXIT:
                  current_ptr = current_ptr->prev_ptr;
                  break;
            }

            if (current_ptr->neighbor_ptr != first_intersection_ptr)
            {
               const Point_2D& point = ((const Point_2D&)(*current_ptr));
               polygon.add (point, false);
            }

         }
         while (current_ptr->neighbor_ptr == NULL);

         current_ptr = current_ptr->neighbor_ptr;

         remove_vector.push_back (current_ptr->neighbor_ptr);
         remove_vector.push_back (current_ptr);

         current_ptr->neighbor_ptr->neighbor_ptr = NULL;
         current_ptr->neighbor_ptr = NULL;

      }
      while (current_ptr != first_intersection_ptr);//not closing the polygon

   }

   boolean_op_phase_c_tilde (polygon, polygon_a, polygon_b, boolean_op);

   for (vector<Polygon::Vertex*>::iterator iterator = remove_vector.begin ();
        iterator != remove_vector.end (); iterator++)
   {
      Polygon::Vertex* vertex_ptr = *(iterator);
      (vertex_ptr->handle_ptr->intersections)--;
      Polygon::remove (vertex_ptr);
   }

}

void
Polygon::boolean_op_phase_c_tilde (Polygon& polygon,
                                   const Polygon& polygon_a,
                                   const Polygon& polygon_b,
                                   const Boolean_Op boolean_op)
{

   const Polygon::Vertex* first_handle_a_ptr = polygon_a.first_handle_ptr;
   const Polygon::Vertex* first_handle_b_ptr = polygon_b.first_handle_ptr;

   Polygon::Vertex* current_handle_a_ptr = (Polygon::Vertex*)first_handle_a_ptr;

   do
   {

      if (current_handle_a_ptr->intersections == 0)
      {

         Polygon::Vertex* current_handle_b_ptr = (Polygon::Vertex*)first_handle_b_ptr;

         do
         {

            Polygon::Relation relation = get_relation (
               current_handle_a_ptr, current_handle_b_ptr);

            if ((boolean_op == INTERSECTION && relation == Polygon::A_ENTIRELY_IN_B) ||
                (boolean_op == DIFFERENCE && relation == Polygon::B_ENTIRELY_IN_A) ||
                (boolean_op == DIFFERENCE && relation == Polygon::DISJOINT) ||
                (boolean_op == UNION && relation == Polygon::B_ENTIRELY_IN_A) ||
                (boolean_op == UNION && relation == Polygon::DISJOINT))
            {
               if (boolean_op == INTERSECTION &&
                   current_handle_a_ptr->action == APPEND)
               {
                  current_handle_a_ptr->action = NOT_DECIDED;
               }
               else
               if (current_handle_a_ptr->action != DO_NOT_APPEND)
               {
                  current_handle_a_ptr->action = APPEND;
               }
            }
            else
            if ((boolean_op == DIFFERENCE && relation == Polygon::A_ENTIRELY_IN_B) ||
                (boolean_op == UNION && relation == Polygon::A_ENTIRELY_IN_B))
            {
               current_handle_a_ptr->action = DO_NOT_APPEND;
            }

            current_handle_b_ptr = current_handle_b_ptr->next_handle_ptr;

         }
         while (current_handle_b_ptr != first_handle_b_ptr);

      }

      current_handle_a_ptr = current_handle_a_ptr->next_handle_ptr;

   }
   while (current_handle_a_ptr != first_handle_a_ptr);

   Polygon::Vertex* current_handle_b_ptr = (Polygon::Vertex*)first_handle_b_ptr;

   do
   {

      if (current_handle_b_ptr->intersections == 0)
      {

         Polygon::Vertex* current_handle_a_ptr = (Polygon::Vertex*)first_handle_a_ptr;

         do
         {

            Polygon::Relation relation = get_relation (
               current_handle_a_ptr, current_handle_b_ptr);

            if ((boolean_op == INTERSECTION && relation == Polygon::B_ENTIRELY_IN_A) ||
                (boolean_op == DIFFERENCE && relation == Polygon::B_ENTIRELY_IN_A) ||
                (boolean_op == UNION && relation == Polygon::A_ENTIRELY_IN_B) ||
                (boolean_op == UNION && relation == Polygon::DISJOINT))
            {
               if (boolean_op == INTERSECTION &&
                   current_handle_b_ptr->action == APPEND)
               {
                  current_handle_b_ptr->action = NOT_DECIDED;
               }
               else
               if (current_handle_b_ptr->action != DO_NOT_APPEND)
               {
                  current_handle_b_ptr->action = APPEND;
               }
            }
            else
            if ((boolean_op == UNION && relation == Polygon::B_ENTIRELY_IN_A))
            {
               current_handle_a_ptr->action = DO_NOT_APPEND;
            }

            current_handle_a_ptr = current_handle_a_ptr->next_handle_ptr;

         }
         while (current_handle_a_ptr != first_handle_a_ptr);

      }

      current_handle_b_ptr = current_handle_b_ptr->next_handle_ptr;

   }
   while (current_handle_b_ptr != first_handle_b_ptr);

   {

      Polygon::Vertex* current_handle_a_ptr = (Polygon::Vertex*)first_handle_a_ptr;
      do
      {

         if (current_handle_a_ptr->action == APPEND)
         {
            polygon.append (current_handle_a_ptr);
         }

         current_handle_a_ptr->action = NOT_DECIDED;
         current_handle_a_ptr = current_handle_a_ptr->next_handle_ptr;

      }
      while (current_handle_a_ptr != first_handle_a_ptr);

      Polygon::Vertex* current_handle_b_ptr = (Polygon::Vertex*)first_handle_b_ptr;
      do
      {

         if (current_handle_b_ptr->action == APPEND)
         {
            polygon.append (current_handle_b_ptr);
         }

         current_handle_b_ptr->action = NOT_DECIDED;
         current_handle_b_ptr = current_handle_b_ptr->next_handle_ptr;

      }
      while (current_handle_b_ptr != first_handle_b_ptr);

   }

}

void
Polygon::simplify_phase_a (Polygon& polygon)
{

   Point_2D ip;
   set<Polygon::Vertex*> pv_ptr_set;
   Polygon::Vertex* first_handle_ptr = polygon.first_handle_ptr;

   Polygon::Vertex* a_current_handle_ptr = first_handle_ptr;
   do
   {

      Polygon::Vertex* a_current_ptr = a_current_handle_ptr;
      do
      {

         Polygon::Vertex* b_current_handle_ptr = a_current_handle_ptr;
         do
         {

            Polygon::Vertex* b_current_ptr = a_current_ptr;
            if (b_current_handle_ptr != a_current_handle_ptr)
            {
               b_current_ptr = b_current_handle_ptr;
            }

            do
            {

               Polygon::Vertex* a_next_ptr = a_current_ptr->next_ptr;
               Polygon::Vertex* b_next_ptr = b_current_ptr->next_ptr;

               const Point_2D& pa = ((Point_2D&)(*a_current_ptr));
               const Point_2D& npa = ((Point_2D&)(*a_next_ptr));
               const Point_2D& pb = ((Point_2D&)(*b_current_ptr));
               const Point_2D& npb = ((Point_2D&)(*b_next_ptr));

               if (pa != npa && pb != npb &&
                   pa != pb && pa != npb && npa != pb && npa != npb &&
                   a_current_ptr != b_current_ptr)
               {

                  if (Edge::intersection (ip, pa, npa, pb, npb))
                  {

                     Polygon::Vertex* iva_ptr = new Polygon::Vertex (ip);
                     Polygon::Vertex* ivb_ptr = new Polygon::Vertex (ip);

                     iva_ptr->intersection = INTERSECTION_SELF;
                     iva_ptr->prev_ptr = a_current_ptr;
                     iva_ptr->next_ptr = a_next_ptr;
                     iva_ptr->handle_ptr = a_current_handle_ptr;
                     iva_ptr->neighbor_ptr = ivb_ptr;

                     ivb_ptr->intersection = INTERSECTION_SELF;
                     ivb_ptr->prev_ptr = b_current_ptr;
                     ivb_ptr->next_ptr = b_next_ptr;
                     ivb_ptr->handle_ptr = b_current_handle_ptr;
                     ivb_ptr->neighbor_ptr = iva_ptr;

                     a_current_ptr->next_ptr = iva_ptr;
                     a_next_ptr->prev_ptr = iva_ptr;

                     b_current_ptr->next_ptr = ivb_ptr;
                     b_next_ptr->prev_ptr = ivb_ptr;

                     pv_ptr_set.insert (iva_ptr);

                  }

               }

               b_current_ptr = b_current_ptr->next_ptr;
            }
            while (b_current_ptr != b_current_handle_ptr);

            b_current_handle_ptr = b_current_handle_ptr->next_handle_ptr;
         }
         while (b_current_handle_ptr != first_handle_ptr);

         a_current_ptr = a_current_ptr->next_ptr;
      }
      while (a_current_ptr != a_current_handle_ptr);

      a_current_handle_ptr = a_current_handle_ptr->next_handle_ptr;
   }
   while (a_current_handle_ptr != first_handle_ptr);

   for (set<Polygon::Vertex*>::iterator iterator = pv_ptr_set.begin ();
        iterator != pv_ptr_set.end (); iterator++)
   {

      Polygon::Vertex* iva_ptr = *(iterator);
      Polygon::Vertex* ivb_ptr = iva_ptr->neighbor_ptr;

      Polygon::Vertex* next_a_ptr = iva_ptr->next_ptr;
      Polygon::Vertex* next_b_ptr = ivb_ptr->next_ptr;

      iva_ptr->next_ptr = next_b_ptr;
      ivb_ptr->next_ptr = next_a_ptr;

      iva_ptr->prev_ptr->next_ptr = iva_ptr;
      iva_ptr->next_ptr->prev_ptr = iva_ptr;

      ivb_ptr->prev_ptr->next_ptr = ivb_ptr;
      ivb_ptr->next_ptr->prev_ptr = ivb_ptr;

   }

}

void
Polygon::simplify_phase_b (Polygon& polygon)
{

   Polygon::Vertex* first_handle_ptr = polygon.first_handle_ptr;
   Polygon::Vertex* current_handle_ptr = first_handle_ptr;
   do
   {
      Polygon::Vertex* next_handle_ptr = current_handle_ptr->next_handle_ptr;

      current_handle_ptr->prev_handle_ptr = NULL;
      current_handle_ptr->next_handle_ptr = NULL;

      current_handle_ptr = next_handle_ptr;
   }
   while (current_handle_ptr != first_handle_ptr);

   first_handle_ptr->next_handle_ptr = first_handle_ptr;
   first_handle_ptr->prev_handle_ptr = first_handle_ptr;

   set<Polygon::Vertex*> pv_ptr_set;
   pv_ptr_set.insert (first_handle_ptr);

   while (pv_ptr_set.size () > 0)
   {

      Polygon::Vertex* handle_ptr = *(pv_ptr_set.begin ());
      pv_ptr_set.erase (handle_ptr);

      if (handle_ptr != first_handle_ptr)
      {

         Polygon::Vertex* last_handle_ptr = first_handle_ptr->prev_handle_ptr;

         handle_ptr->next_handle_ptr = first_handle_ptr;
         handle_ptr->prev_handle_ptr = last_handle_ptr;

         first_handle_ptr->prev_handle_ptr = handle_ptr;
         last_handle_ptr->next_handle_ptr = handle_ptr;

      }

      Polygon::Vertex* current_ptr = handle_ptr;

      do
      {

         if (current_ptr->intersection == INTERSECTION_SELF)
         {

            Polygon::Vertex* neighbor_ptr = current_ptr->neighbor_ptr;

            current_ptr->neighbor_ptr = NULL;
            neighbor_ptr->neighbor_ptr = NULL;
            current_ptr->intersection = INTERSECTION_NOT;
            neighbor_ptr->intersection = INTERSECTION_NOT;

            Polygon::Vertex* candidate_handle_ptr = neighbor_ptr;
            pv_ptr_set.insert (candidate_handle_ptr);

         }

         pv_ptr_set.erase (current_ptr);

         current_ptr->handle_ptr = handle_ptr;
         current_ptr = current_ptr->next_ptr;
      }
      while (current_ptr != handle_ptr);

   }

}

Polygon::Relation
Polygon::get_relation (const Polygon::Vertex* handle_a_ptr,
                       const Polygon::Vertex* handle_b_ptr)
{

   Polygon::Relation relation = Polygon::DISJOINT;

   if (entirely_within (handle_a_ptr, handle_b_ptr))
   {
      relation = Polygon::A_ENTIRELY_IN_B;
   }
   else
   if (entirely_within (handle_b_ptr, handle_a_ptr))
   {
      relation = Polygon::B_ENTIRELY_IN_A;
   }

   return relation;

}

bool
Polygon::entirely_within (const Polygon::Vertex* handle_a_ptr,
                          const Polygon::Vertex* handle_b_ptr)
{

   bool b = true;
   Polygon::Vertex* current_a_ptr = (Polygon::Vertex*)handle_a_ptr;

   do
   {
      const Point_2D& p = (const Point_2D&)(*current_a_ptr);
      if (!Polygon::contains (handle_b_ptr, p))
      {
         b = false;
         break;
      }
      current_a_ptr = current_a_ptr->next_ptr;
   }
   while (current_a_ptr != handle_a_ptr && b);

   return b;

}

bool
Polygon::append (Polygon::Vertex* handle_ptr)
{

   bool first_point = true;
   Polygon::Vertex* current_ptr = (Polygon::Vertex*)handle_ptr;

   do
   {
      const Point_2D& p = ((const Point_2D&)(*current_ptr));
      add (p, first_point);
      if (first_point) { first_point = false; }
      current_ptr = current_ptr->next_ptr;
   }
   while (current_ptr != handle_ptr);

}

Polygon::Polygon ()
   : first_handle_ptr (NULL),
     domain_x (GSL_POSINF, GSL_NEGINF),
     domain_y (GSL_POSINF, GSL_NEGINF),
     min_edge_length (GSL_POSINF)
{
}

Polygon::Polygon (const Polygon& polygon)
   : first_handle_ptr (NULL),
     domain_x (GSL_POSINF, GSL_NEGINF),
     domain_y (GSL_POSINF, GSL_NEGINF),
     min_edge_length (GSL_POSINF)
{
   add (polygon);
}

Polygon::Polygon (const Circle& circle,
                  const Integer n)
   : first_handle_ptr (NULL),
     domain_x (GSL_POSINF, GSL_NEGINF),
     domain_y (GSL_POSINF, GSL_NEGINF),
     min_edge_length (GSL_POSINF)
{

   const Real dt = M_2_TIMES_PI / Real (n);

   const Point_2D& center = circle.center;
   const Real& radius = circle.radius;

   for (Integer i = 0; i < n; i++)
   {

      const Real t = i * dt;

      const Real dx = radius * cos (t);
      const Real dy = radius * sin (t);

      add (Point_2D (center.x + radius * cos (t),
                     center.y + radius * sin (t)));

   }

}

Polygon::Polygon (const Ellipse& ellipse,
                  const Integer n)
   : first_handle_ptr (NULL),
     domain_x (GSL_POSINF, GSL_NEGINF),
     domain_y (GSL_POSINF, GSL_NEGINF),
     min_edge_length (GSL_POSINF)
{

   const Real dt = M_2_TIMES_PI / Real (n);

   const Point_2D& center = ellipse.center;
   const Real& a = ellipse.a;
   const Real& b = ellipse.b;
   const Real& tilt = ellipse.tilt;

   for (Integer i = 0; i < n; i++)
   {

      const Real t = i * dt;

      const Real dx = a * cos (t);
      const Real dy = b * sin (t);

      const Point_2D point (center.x + cos (tilt) * dx - sin (tilt) * dy,
                            center.y + sin (tilt) * dx + cos (tilt) * dy);

      add (point);

   }

}

Polygon::Polygon (const Domain_2D& domain_2d,
                  const Real delta_x,
                  Real delta_y)
   : first_handle_ptr (NULL),
     domain_x (GSL_POSINF, GSL_NEGINF),
     domain_y (GSL_POSINF, GSL_NEGINF),
     min_edge_length (GSL_POSINF)
{

   if (!gsl_finite (delta_y)) { delta_y = delta_x; }

   const Real start_x = domain_2d.domain_x.start;
   const Real end_x = domain_2d.domain_x.end;
   const Real start_y = domain_2d.domain_y.start;
   const Real end_y = domain_2d.domain_y.end;
   const Real span_x = end_x - start_x;
   const Real span_y = end_y - start_y;
   const Integer n_x = Integer (ceil (span_x / delta_x)) + 1;
   const Integer n_y = Integer (ceil (span_y / delta_y)) + 1;
   const Real d_x = span_x / (n_x - 1);
   const Real d_y = span_y / (n_y - 1);

   for (Integer i = 0; i < n_x; i++)
   {
      const Real x = start_x + i * d_x;
      add (Point_2D (x, start_y));
   }

   for (Integer j = 1; j < n_y; j++)
   {
      const Real y = start_y + j * d_y;
      add (Point_2D (end_x, y));
   }

   for (Integer i = n_x - 1; i > 0; i--)
   {
      const Real x = start_x  + i * d_x;
      add (Point_2D (x, end_y));
   }

   for (Integer j = n_y - 1; j > 0; j--)
   {
      const Real y = start_y  + j * d_y;
      add (Point_2D (start_x, y));
   }


}

Polygon::~Polygon ()
{
   clear ();
}

void
Polygon::operator= (const Polygon& polygon)
{

   first_handle_ptr = NULL;
   domain_x = Domain_1D (GSL_POSINF, GSL_NEGINF);
   domain_y = Domain_1D (GSL_POSINF, GSL_NEGINF);

   add (polygon);

}

void
Polygon::clear ()
{

   if (size () == 0) { return; }
   Polygon::Vertex* current_handle_ptr = first_handle_ptr;

   while (true)
   {

      Polygon::Vertex* next_handle_ptr = current_handle_ptr->next_handle_ptr;
      Polygon::Vertex* current_ptr = current_handle_ptr;

      const Integer n = current_handle_ptr->n;

      Polygon::Vertex* next_ptr = current_ptr->next_ptr;

      while (true)
      {

         delete current_ptr;
         current_ptr = next_ptr;

         if (current_ptr == current_handle_ptr) { break; }
         next_ptr = current_ptr->next_ptr;

      }

      if (current_handle_ptr == first_handle_ptr) { break; }
      current_handle_ptr = next_handle_ptr;

   }

   first_handle_ptr = NULL;

   domain_x.start = GSL_POSINF;
   domain_x.end = GSL_NEGINF;
   domain_y.start = GSL_POSINF;
   domain_y.end = GSL_NEGINF;
   min_edge_length = GSL_POSINF;

}

void
Polygon::attract_by (const Attractor& attractor)
{

   const Polygon::Vertex* first_handle_ptr = get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {

      const Integer n = current_handle_ptr->n;
      Polygon::Vertex* current_ptr = (Polygon::Vertex*)current_handle_ptr;

      for (Integer i = 0; i < n; i++)
      {
         Point_2D& point = (Point_2D&)(*current_ptr);
         attractor.attract (point);
         current_ptr = current_ptr->next_ptr;
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);

}

void
Polygon::add (const Point_2D& point,
              const bool new_handle)
{

   Polygon::Vertex* vertex_ptr = new Polygon::Vertex (point);

   if (first_handle_ptr == NULL)
   {

      vertex_ptr->prev_ptr = vertex_ptr;
      vertex_ptr->next_ptr = vertex_ptr;
      vertex_ptr->handle_ptr = vertex_ptr;

      vertex_ptr->prev_handle_ptr = vertex_ptr;
      vertex_ptr->next_handle_ptr = vertex_ptr;

      first_handle_ptr = vertex_ptr;
      vertex_ptr->n = 1;

   }
   else
   {

      Polygon::Vertex* last_handle_ptr = first_handle_ptr->prev_handle_ptr;

      if (new_handle) // main handle_ptr != NULL
      {

         vertex_ptr->next_ptr = vertex_ptr;
         vertex_ptr->prev_ptr = vertex_ptr;
         vertex_ptr->handle_ptr = vertex_ptr;

         vertex_ptr->prev_handle_ptr = last_handle_ptr;
         vertex_ptr->next_handle_ptr = first_handle_ptr;

         first_handle_ptr->prev_handle_ptr = vertex_ptr;
         last_handle_ptr->next_handle_ptr = vertex_ptr;
         vertex_ptr->n = 1;

      }
      else
      {

         Polygon::Vertex* last_ptr = last_handle_ptr->prev_ptr;

         if (point.x == last_ptr->x && point.y == last_ptr->y)
         {
            delete vertex_ptr;
            return;
         }

         vertex_ptr->next_ptr = last_handle_ptr;
         vertex_ptr->prev_ptr = last_ptr;

         last_handle_ptr->prev_ptr = vertex_ptr;
         last_ptr->next_ptr = vertex_ptr;

         vertex_ptr->handle_ptr = last_handle_ptr;
         (last_handle_ptr->n)++;

      }

   }

   if (point.x < domain_x.start) { domain_x.start = point.x; }
   else if (point.x < domain_x.end) { domain_x.end = point.x; }

   if (point.y < domain_y.start) { domain_y.start = point.y; }
   else if (point.y < domain_y.end) { domain_y.end = point.y; }

   if (vertex_ptr->prev_ptr != vertex_ptr)
   {

      const Point_2D& prev_point = (const Point_2D&)(*(vertex_ptr->prev_ptr));
      const Real dx = point.x - prev_point.x;
      const Real dy = point.y - prev_point.y;
      const Real edge_length = sqrt (dx*dx + dy*dy);

      if (edge_length < min_edge_length)
      {
         min_edge_length = edge_length;
      }

   }

}

void
Polygon::add (const Polygon& polygon)
{

   const Polygon::Vertex* first_handle_ptr = polygon.get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {

      Integer n = current_handle_ptr->n;
      Polygon::Vertex* current_ptr = (Polygon::Vertex*)current_handle_ptr;

      for (Integer i = 0; i < n; i++)
      {
         const Point_2D& point = (const Point_2D&)(*current_ptr);
         add (point, (i == 0));
         current_ptr = current_ptr->next_ptr;
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);

}

void
Polygon::add (const Polygon& polygon,
              const Point_2D& offset)
{

   const Polygon::Vertex* first_handle_ptr = polygon.get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {

      Integer n = current_handle_ptr->n;
      Polygon::Vertex* current_ptr = (Polygon::Vertex*)current_handle_ptr;

      for (Integer i = 0; i < n; i++)
      {
         const Point_2D& point = (const Point_2D&)(*current_ptr);
         add (point + offset, (i == 0));
         current_ptr = current_ptr->next_ptr;
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);
   
}

void
Polygon::add (const Simple_Polyline& simple_polyline,
              const bool new_handle)
{
   for (Simple_Polyline::const_iterator iterator = simple_polyline.begin ();
        iterator != simple_polyline.end (); iterator++)
   {
      const Point_2D& point = *(iterator);
      add (point, (new_handle && iterator == simple_polyline.begin ()));
   }
}

void
Polygon::transform (const Transform_2D& transform)
{

   const Polygon::Vertex* first_handle_ptr = get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {
      Polygon::Vertex* current_ptr = current_handle_ptr;
      do
      {
         Real& x = current_ptr->x;
         Real& y = current_ptr->y;
         transform.transform (x, y, x, y);
         current_ptr = current_ptr->next_ptr;
      }
      while (current_ptr != current_handle_ptr);
      current_handle_ptr = current_handle_ptr->next_handle_ptr;
   }
   while (current_handle_ptr != first_handle_ptr);

}

void
Polygon::reverse (const Transform_2D& transform)
{

   const Polygon::Vertex* first_handle_ptr = get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {
      Polygon::Vertex* current_ptr = current_handle_ptr;
      do
      {
         Real& x = current_ptr->x;
         Real& y = current_ptr->y;
         transform.reverse (x, y, x, y);
         current_ptr = current_ptr->next_ptr;
      }
      while (current_ptr != current_handle_ptr);
      current_handle_ptr = current_handle_ptr->next_handle_ptr;
   }
   while (current_handle_ptr != first_handle_ptr);

}

void
Polygon::simplify ()
{

   if (size () <= 1) { return; }

   simplify_phase_a (*this);
   simplify_phase_b (*this);
}

void
Polygon::translate (const Point_2D& point)
{

   if (size () <= 1) { return; }

   const Polygon::Vertex* first_handle_ptr = get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {

      const Integer n = current_handle_ptr->n;
      Polygon::Vertex* current_ptr = (Polygon::Vertex*)current_handle_ptr;

      for (Integer i = 0; i < n; i++)
      {
         Point_2D& p = (Point_2D&)(*current_ptr);
         p.x += point.x;
         p.y += point.y;
         current_ptr = current_ptr->next_ptr;
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);

}

void
Polygon::remove (Polygon::Vertex* vertex_ptr)
{

   (vertex_ptr->handle_ptr->n)--;

   Polygon::Vertex* prev_ptr = vertex_ptr->prev_ptr;
   Polygon::Vertex* next_ptr = vertex_ptr->next_ptr;

   if (vertex_ptr->handle_ptr == vertex_ptr) // is a handle
   {

      // next_ptr is new handle
      for (Polygon::Vertex* current_ptr = next_ptr;
           current_ptr != vertex_ptr;
           current_ptr = current_ptr->next_ptr)
      {
         current_ptr->handle_ptr = next_ptr;
      }

      next_ptr->n = vertex_ptr->n;
      next_ptr->prev_handle_ptr = vertex_ptr->prev_handle_ptr;
      next_ptr->next_handle_ptr = vertex_ptr->next_handle_ptr;

   }

   prev_ptr->next_ptr = next_ptr;
   next_ptr->prev_ptr = prev_ptr;

   delete vertex_ptr;

}

Integer
Polygon::size () const
{

   if (first_handle_ptr == NULL) { return 0; }
   else
   {

      Integer n = 0;

      Polygon::Vertex* current_handle_ptr = first_handle_ptr;
      do
      {
         n += current_handle_ptr->n;
         current_handle_ptr = current_handle_ptr->next_handle_ptr;
      }
      while (current_handle_ptr != first_handle_ptr);

      return n;

   }

}

Integer
Polygon::get_number_of_single_polygons () const
{

   if (first_handle_ptr == NULL) { return 0; }
   else
   {

      Integer n = 0;

      Polygon::Vertex* current_handle_ptr = first_handle_ptr;
      do
      {
         n++;
         current_handle_ptr = current_handle_ptr->next_handle_ptr;
      }
      while (current_handle_ptr != first_handle_ptr);

      return n;

   }

}

const Polygon::Vertex*
Polygon::get_first_handle_ptr () const
{
   return first_handle_ptr;
}

const Polygon::Vertex*
Polygon::get_last_handle_ptr () const
{
   return first_handle_ptr->prev_handle_ptr;
}

const Domain_1D&
Polygon::get_domain_x () const
{
   return domain_x;
}

const Domain_1D&
Polygon::get_domain_y () const
{
   return domain_y;
}

Real
Polygon::get_simple_polygon_area (const Polygon::Vertex* pv_ptr)
{

   Real area = 0;
   Polygon::Vertex* current_ptr = ((Polygon::Vertex*)pv_ptr);

   do
   {

      Polygon::Vertex* next_ptr = current_ptr->next_ptr;
      const Point_2D& current_point = (const Point_2D&)(*current_ptr);
      const Point_2D& next_point = (const Point_2D&)(*next_ptr);

      area += next_point.x * current_point.y;
      area -= next_point.y * current_point.x;

      current_ptr = next_ptr;

   }
   while (current_ptr != pv_ptr);

   return area / 2;

}

Real
Polygon::get_area () const
{

   Real area = 0;
   Polygon::Vertex* current_handle_ptr = first_handle_ptr;

   do
   {
      area += get_simple_polygon_area (current_handle_ptr);
      current_handle_ptr = current_handle_ptr->next_handle_ptr;
   }
   while (current_handle_ptr != first_handle_ptr);

   return area;

}

Real
Polygon::get_positive_area () const
{

   Real area = 0;
   Polygon::Vertex* current_handle_ptr = first_handle_ptr;

   do
   {
      Real sub_area = get_simple_polygon_area (current_handle_ptr);
      if (sub_area > 0) { area += sub_area; }
      current_handle_ptr = current_handle_ptr->next_handle_ptr;
   }
   while (current_handle_ptr != first_handle_ptr);

   return area;

}

Real
Polygon::get_negative_area () const
{

   Real area = 0;
   Polygon::Vertex* current_handle_ptr = first_handle_ptr;

   do
   {
      Real sub_area = get_simple_polygon_area (current_handle_ptr);
      if (sub_area < 0) { area += sub_area; }
      current_handle_ptr = current_handle_ptr->next_handle_ptr;
   }
   while (current_handle_ptr != first_handle_ptr);

   return area;

}

Point_2D
Polygon::get_centroid () const
{

   Integer n = 0;
   Point_2D centroid (0, 0);
   Polygon::Vertex* current_handle_ptr = first_handle_ptr;

   do
   {

      n += current_handle_ptr->n;
      Polygon::Vertex* current_ptr = current_handle_ptr;

      for (Integer i = 0; i < current_handle_ptr->n; i++)
      {
         const Point_2D& point = (const Point_2D&)(*current_ptr);
         centroid.x += point.x;
         centroid.y += point.y;
         current_ptr = current_ptr->next_ptr;
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);

   centroid.x /= Real (n);
   centroid.y /= Real (n);

   return centroid;

}

Integer
Polygon::get_winding_number (const Point_2D& point) const
{

   if (size () == 0) { return 0; }

   Integer w = 0;
   Polygon::Vertex* current_handle_ptr = first_handle_ptr;

   do
   {
      iterate_winding_number (w, current_handle_ptr, point);
      current_handle_ptr = current_handle_ptr->next_handle_ptr;
   }
   while (current_handle_ptr != first_handle_ptr);

   return w;

}

bool
Polygon::contains (const Point_2D& point,
                   bool border_included) const
{
   try { return abs (get_winding_number (point) % 2); }
   catch (const Geometry_Exception& ge) { return border_included; }
}

void
Polygon::boolean_op (Polygon& polygon,
                     const Boolean_Op boolean_op,
                     const Polygon& polygon_a,
                     const Polygon& polygon_b,
                     Real perturbation)
{

   bool exit_entry_a;
   bool exit_entry_b;

   switch (boolean_op)
   {
      case DIFFERENCE:   exit_entry_a = true;  exit_entry_b = false; break;
      case INTERSECTION: exit_entry_a = false; exit_entry_b = false; break;
      case UNION:        exit_entry_a = true;  exit_entry_b = true;  break;
   }

   if (polygon_a.size () > 2 && gsl_isnan (perturbation))
   {
      perturbation = min (polygon_a.min_edge_length, polygon_b.min_edge_length);
      perturbation *= 1e-1;
   }

   if (polygon_a.size () > 2 && polygon_b.size () > 2)
   {

      boolean_op_phase_0 (polygon_a, polygon_b, perturbation);

      boolean_op_phase_a (polygon_a, polygon_b);

//polygon_b.debug_print ("debug b ");
//      if (polygon_a.first_handle_ptr->intersections % 2 == 1 ||
//          polygon_b.first_handle_ptr->intersections % 2 == 1)
//      {
//         //throw Geometry_Exception ("Boolean Op: odd intersections");
//         Integer a = polygon_a.first_handle_ptr->intersections;
//         Integer b = polygon_b.first_handle_ptr->intersections;
//cout << "failure " << a << " " << b << endl;
//         return;
//      }

      boolean_op_phase_b (polygon_a, polygon_b, exit_entry_a);
      boolean_op_phase_b (polygon_b, polygon_a, exit_entry_b);

      boolean_op_phase_c (polygon, polygon_a, polygon_b, boolean_op);

   }

}

Polygon*
Polygon::boolean_op (const Boolean_Op boolean_op,
                     const Polygon& polygon_a,
                     const Polygon& polygon_b,
                     Real perturbation)
{
   Polygon* polygon_ptr = new Polygon ();
   Polygon::boolean_op (*polygon_ptr, boolean_op,
      polygon_a, polygon_b, perturbation);
   return polygon_ptr;
}

void
Polygon::cairo (const RefPtr<Context>& cr) const
{
   cairo (cr, Point_2D (0, 0));
}

void
Polygon::cairo (const RefPtr<Context>& cr,
                const Point_2D& offset) const
{

   if (size () <= 1) { return; }

   const Polygon::Vertex* first_handle_ptr = get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {

      const Integer n = current_handle_ptr->n;
      Polygon::Vertex* current_ptr = (Polygon::Vertex*)current_handle_ptr;

      for (Integer i = 0; i < n; i++)
      {
         const Point_2D& point = (const Point_2D&)(*current_ptr) + offset;
         if (i == 0) { cr->move_to (point.x, point.y); }
         else        { cr->line_to (point.x, point.y); }
         current_ptr = current_ptr->next_ptr;
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

      cr->close_path ();

   }
   while (current_handle_ptr != first_handle_ptr);

}

void
Polygon::cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d) const
{

   if (size () <= 1) { return; }

   const Polygon::Vertex* first_handle_ptr = get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {

      const Integer n = current_handle_ptr->n;
      Polygon::Vertex* current_ptr = (Polygon::Vertex*)current_handle_ptr;

      for (Integer i = 0; i < n; i++)
      {
         const Point_2D& point = (const Point_2D&)(*current_ptr);
         const Point_2D& t_point = transform_2d.transform (point);
         if (i == 0) { cr->move_to (t_point.x, t_point.y); }
         else        { cr->line_to (t_point.x, t_point.y); }
         current_ptr = current_ptr->next_ptr;
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);

   cr->close_path ();

}

void
Polygon::cairo (const RefPtr<Context>& cr,
                const Polygon::Vertex* polygon_vertex_ptr)
{

   Polygon::Vertex* current_ptr = (Polygon::Vertex*)polygon_vertex_ptr;

   do
   {

      const Point_2D& point = (const Point_2D&)(*current_ptr);

      if (current_ptr == polygon_vertex_ptr)
      {
         cr->move_to (point.x, point.y);
      }
      else
      {
         cr->line_to (point.x, point.y);
      }

      current_ptr = current_ptr->next_ptr;

   }
   while (current_ptr != polygon_vertex_ptr);
   cr->close_path ();

}

void
Polygon::debug_print (const Dstring& prefix,
                      ostream& out_stream) const
{

   const Polygon::Vertex* handle_ptr = first_handle_ptr;
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {

      out_stream << prefix << "new handle" << endl;
      Integer n = current_handle_ptr->n;
      Polygon::Vertex* current_ptr = current_handle_ptr;

      do
      {
         const Point_2D& point = (const Point_2D&)(*current_ptr);
         Dstring str ("");
         switch (current_ptr->intersection)
         {
            case INTERSECTION_NOT: str = "NOT_INTERSECTION"; break;
            case INTERSECTION_ENTRY: str = "ENTRY"; break;
            case INTERSECTION_EXIT: str = "EXIT"; break;
         }
         if (current_ptr->neighbor_ptr != NULL) { str += " HAS_NEIGHBOR"; }

         out_stream << prefix << current_ptr << " " << point << " " << str;
         out_stream << " | prev_ptr = " << current_ptr->prev_ptr;
         out_stream << " | next_ptr = " << current_ptr->next_ptr;
         out_stream << " | neighbor = " << current_ptr->neighbor_ptr << endl;

         current_ptr = current_ptr->next_ptr;

      }
      while (current_ptr != current_handle_ptr);

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);

   out_stream << prefix << "end polygon" << endl;

}

Real
Parametric::get_x (const Real t) const
{
   return t;
}

Point_2D
Parametric::get_point (const Real t) const
{
   const Real x = get_x (t);
   const Real y = get_y (t);
   return Point_2D (x, y);
}

Real
Parametric::get_weight (const Real t) const
{
   return 1;
}

Tuple
Parametric::get_t_tuple (const Integer n,
                         const Real start,
                         const Real end,
                         const Real w0) const
{

   const Real L = end - start;
   const Real delta = L / Real (n - 1);

   Real total_weight = 0;
   Tuple weight_tuple;
   for (Integer i = 0; i < n - 1; i++)
   {
      const Real t = start + i * delta;
      const Real weight = get_weight (t) + w0;
      weight_tuple.push_back (weight);
      total_weight += weight;
   }
   const Real k = L / total_weight;

   Real t = start;
   Tuple t_tuple (1, t);
   for (Tuple::const_iterator iterator = weight_tuple.begin ();
        iterator != weight_tuple.end (); iterator++)
   {
      const Real weight = *(iterator);
      const Real delta_t = k * weight;
      t += delta_t;
      t_tuple.push_back (t);
   }

   return t_tuple;

}

Simple_Polyline*
Parametric::get_simple_polyline_ptr (const Integer n,
                                     const Real start,
                                     const Real end,
                                     const Real w0) const
{

   Simple_Polyline* sp_ptr = new Simple_Polyline ();
   const Tuple& t_tuple = get_t_tuple (n, start, end, w0);

   for (Tuple::const_iterator iterator = t_tuple.begin ();
        iterator != t_tuple.end (); iterator++)
   {
      const Real& t = *(iterator);
      const Point_2D& point = get_point (t);
      sp_ptr->add (point);
   }

   return sp_ptr;

}

void
Parametric::cairo (const RefPtr<Context>& cr,
                   const Transform_2D& transform_2d,
                   const Integer n,
                   const Real start,
                   const Real end,
                   const Real w0) const
{
   Simple_Polyline* sp_ptr = get_simple_polyline_ptr (n, start, end, w0);
   sp_ptr->cairo (cr, transform_2d);
   delete sp_ptr;
}

Polar_Parametric::Polar_Parametric (const Point_2D& centre)
   : centre (centre)
{
}

Real
Polar_Parametric::get_x (const Real theta) const
{
   return centre.x + get_r (theta) * cos (theta);
}

Real
Polar_Parametric::get_y (const Real theta) const
{
   return centre.x + get_r (theta) * sin (theta);
}

Point_2D
Polar_Parametric::get_point (const Real theta) const
{
   const Real r = get_r (theta);
   const Real x = centre.x + r * cos (theta);
   const Real y = centre.y + r * sin (theta);
   return Point_2D (x, y);
}

Log_Spiral::Log_Spiral (const Real a,
                        const Real b,
                        const Point_2D& centre)
   : Polar_Parametric (centre),
     a (a),
     b (b)
{
}
 
Log_Spiral::Log_Spiral (const Point_2D& centre,
                        const Point_2D& point,
                        const Real pitch)
   : Polar_Parametric (centre),
     b (tan (pitch))
{

   const Real dx = point.x - centre.x;
   const Real dy = point.y - centre.y;
   const Real r = sqrt (dx * dx + dy * dy);
   const Real theta = atan2 (dy, dx);

   this->a = exp (log (r) - b * theta);

}

Real
Log_Spiral::get_r (const Real theta) const
{
   return a * exp (b * theta);
}

Descartes_Folium::Descartes_Folium (const Real a,
                                    const Point_2D& centre)
   : Polar_Parametric (centre),
     a (a)
{
}

Real
Descartes_Folium::get_r (const Real theta) const
{
   const Real s = sin (theta);
   const Real c = cos (theta);
   return 3 * a * s * c / (s * s * s + c * c * c);
}

Lemniscate::Lemniscate (const Real a,
                        const Point_2D& centre)
   : Polar_Parametric (centre),
     a (a)
{
}

Real
Lemniscate::get_r (const Real theta) const
{
   const Real rhs = a * a * cos (2 * theta);
   return sqrt (rhs < 0 ? 0 : rhs);
}

Pascal_Limacon::Pascal_Limacon (const Real a,
                                const Real b,
                                const Point_2D& centre)
   : Polar_Parametric (centre),
     a (a),
     b (b)
{
}

Real
Pascal_Limacon::get_r (const Real theta) const
{
   return b + a * cos (theta);
}

Cardoid::Cardoid (const Real a,
                  const Point_2D& centre)
   : Pascal_Limacon (a, a, centre)
{
}

Rose::Rose (const Real a,
            const Real n,
            const Real phi,
            const Point_2D& centre)
   : Polar_Parametric (centre),
     a (a),
     n (n),
     phi (phi)
{
}

Real
Rose::get_r (const Real theta) const
{
   return a * cos (theta * n + phi);
}

Trochoid::Trochoid (const Real a,
                    const Real b,
                    const Point_2D& centre)
   : a (a),
     b (b),
     centre (centre)
{
}

Real
Trochoid::get_x (const Real t) const
{
   return centre.x + a * t - b * sin (t);
}

Real
Trochoid::get_y (const Real t) const
{
   return centre.y + a - b * cos (t);
}

Cycloid::Cycloid (const Real a,
                  const Point_2D& centre)
   : Trochoid (a, a, centre)
{
}

void
Rect::translate_anchor (const Real translate_x,
                        const Real translate_y)
{

   Real dx = translate_x;
   Real dy = translate_y;

   if (tilt != 0)
   {

      const Real temp_dx = dx;
      //const Real temp_dy = -dy;
      const Real temp_dy = dy;

      const Real c = cos (tilt);
      const Real s = sin (tilt);

      dx = temp_dx * c - temp_dy * s;
      dy = temp_dx * s + temp_dy * c;

   }

   point_2d.x += dx;
   point_2d.y += dy;

}

void
Rect::construct_polygon ()
{

   const Real c = (tilt == 0 ? 1 : cos (tilt));
   const Real s = (tilt == 0 ? 0 : sin (tilt));

   const Real dx_sw = -height * s;
   const Real dy_sw = height * c;

   const Real dx_ne = width * c;
   const Real dy_ne = width * s;

   const Real dx_se = dx_sw + dx_ne;
   const Real dy_se = dy_sw + dy_ne;

   clear ();
   add (point_2d, true);
   add (Point_2D (point_2d.x + dx_sw, point_2d.y + dy_sw));
   add (Point_2D (point_2d.x + dx_se, point_2d.y + dy_se));
   add (Point_2D (point_2d.x + dx_ne, point_2d.y + dy_ne));

}

void
Rect::init (const Point_2D& point_2d,
            const Real width,
            const Real height,
            const Real translate_x,
            const Real translate_y,
            const Real tilt)
{

   this->tilt = tilt;
   this->width = width;
   this->height = height;
   this->point_2d = point_2d;

   translate_anchor (translate_x, translate_y);
   construct_polygon ();

}

Rect::Rect()
{
}

Rect::Rect(const Box_2D& box_2d)
{

   const Real width = box_2d.size_2d.i;
   const Real height = box_2d.size_2d.j;
   const Point_2D point_2d (box_2d.index_2d.i, box_2d.index_2d.j);

   init (point_2d, width, height, 0, 0, 0);

}

Rect::Rect(const Point_2D& point_a,
           const Point_2D& point_b)
{
   const Real width = point_b.x - point_a.x;
   const Real height = point_b.y - point_a.y;
   init (point_a, width, height, 0, 0, 0);
}

Rect::Rect(const Point_2D& point_2d,
           const Real width,
           const Real height,
           const Real tilt)
{
   init (point_2d, width, height, 0, 0, tilt);
}

Rect::Rect (const Point_2D& point_2d,
            const Real width,
            const Real height,
            const Real translate_x,
            const Real translate_y,
            const Real tilt)
{
   init (point_2d, width, height, translate_x, translate_y, tilt);
}

Rect::Rect (const Point_2D& point_2d,
            const Real width,
            const Real height,
            const char justify_h,
            const char justify_v,
            const Real tilt)
{

   Real translate_x = 0;
   Real translate_y = 0;

   switch (justify_h)
   {
      case 'c': translate_x = -0.5 * width; break;
      case 'r': translate_x = -width; break;
   }

   switch (justify_v)
   {
      case 'c': translate_y = -0.5 * height; break;
      case 'b': translate_y = -height; break;
   }

   init (point_2d, width, height, translate_x, translate_y, tilt);

}

Rect::Rect (const Rect& rect)
{
   init (rect.point_2d, rect.width, rect.height, 0, 0, rect.tilt);
}

void
Rect::set (const Rect& rect)
{
   init (rect.point_2d, rect.width, rect.height, 0, 0, rect.tilt);
}

void
Rect::operator= (const Rect& rect)
{
   init (rect.point_2d, rect.width, rect.height, 0, 0, rect.tilt);
}

void
Rect::grow (const Real margin_x,
            const Real margin_y)
{

   const Real c = (tilt == 0 ? 1 : cos (tilt));
   const Real s = (tilt == 0 ? 0 : sin (tilt));

   width += (margin_x + margin_x);
   height += (margin_y + margin_y);
   point_2d.x -= (margin_x * c - margin_y * s);
   point_2d.y -= (margin_x * s + margin_y * c);

   construct_polygon ();

}

/*
Path_Ob::Path_Ob (const Path_Ob_Genre path_ob_genre)
 : path_ob_genre (path_ob_genre)
{
}

Path_Ob::Path_Ob (const Path_Ob_Genre path_ob_genre,
                  const Point_2D& point)
 : path_ob_genre (path_ob_genre)
{
   push_back (point);
}

Path_Ob::Path_Ob (const Path_Ob_Genre path_ob_genre,
                  const Point_2D& point_0,
                  const Point_2D& point_1,
                  const Point_2D& point_2)
  : path_ob_genre (path_ob_genre)
{
   push_back (point_0);
   push_back (point_1);
   push_back (point_2);
}

Path::Path (const Simple_Polyline& simple_polyline)
{
   Path::simple_polyline (simple_polyline);
   if (simple_polyline.closed) { close (); }
}

Path::Path (const Polyline& polyline)
{

   for (Polyline::const_iterator iterator = polyline.begin ();
        iterator != polyline.end (); iterator++)
   {
      const Simple_Polyline& simple_polyline = **(iterator);
      Path::simple_polyline (simple_polyline);
   }

}

Path::Path ()
{
}

Path::Path (const Polygon& polygon)
{

   const Polygon::Vertex* first_handle_ptr = polygon.get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

   do
   {

      const Integer n = current_handle_ptr->n;
      Polygon::Vertex* current_ptr = (Polygon::Vertex*)current_handle_ptr;

      for (Integer i = 0; i < n; i++)
      {
         const Point_2D& point = (const Point_2D&)(*current_ptr);
         if (i == 0) { move_to (point); }
         else        { line_to (point); }
         current_ptr = current_ptr->next_ptr;
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);

   close ();

}

void
Path::simple_polyline (const Simple_Polyline& simple_polyline)
{

   for (Simple_Polyline::const_iterator iterator = simple_polyline.begin ();
        iterator != simple_polyline.end (); iterator++)
   {

      const Point_2D& point = *(iterator);

      if (iterator == simple_polyline.begin ()) { move_to (point); }
      else { line_to (point); }

   }

}

void
Path::move_to (const Point_2D& point)
{
   push_back (Path_Ob (MOVE_TO, point));
}
         
void
Path::line_to (const Point_2D& point)
{
   push_back (Path_Ob (LINE_TO, point));
}

void
Path::arc_to (const Point_2D& centre,
              const Real radius,
              const Real theta_0,
              const Real theta_1)
{
   push_back (Path_Ob (ARC_TO, centre, Point_2D (radius, radius),
      Point_2D (theta_0, theta_1)));
}

void
Path::curve_to (const Point_2D& point_0,
                const Point_2D& point_1,
                const Point_2D& point_2)
{
   push_back (Path_Ob (CURVE_TO, point_0, point_1, point_2));
}

void
Path::close ()
{
   push_back (Path_Ob (CLOSE_PATH));
}
*/

Arc::Arc (const Point_2D& point,
          const Real radius,
          const Real start_theta,
          const Real end_theta)
 : Simple_Polyline ()
{

   Point_2D p;
   const Real theta_span = end_theta - start_theta;
   const Real proportion = fabs (theta_span) / M_2_TIMES_PI;
   const Integer n = Integer ((log (fabs (radius) + 1) * 18) * proportion) + 5;
   const Real d_theta = theta_span / n;

   for (Integer i = 0; i < n + 1; i++)
   {
      const Real theta = start_theta + i * d_theta;
      p.x = point.x + radius * cos (theta);
      p.y = point.y + radius * sin (theta);
      add (p);
   }

}

Chord::Chord (const Point_2D& point,
              const Real radius,
              const Real start_theta,
              const Real end_theta)
      : Edge (Point_2D (point.x + radius * cos (start_theta), point.y + radius * sin (start_theta)),
              Point_2D (point.x + radius * cos (end_theta), point.y + radius * sin (end_theta)))
{
}

Sector::Sector (const Point_2D& point,
                const Real radius,
                const Real start_theta,
                const Real end_theta)
     : Polygon ()
{

   add (point);

   Point_2D p;
   const Real theta_span = end_theta - start_theta;
   const Real proportion = fabs (theta_span) / M_2_TIMES_PI;
   const Integer n = Integer ((log (fabs (radius) + 1) * 18) * proportion) + 5;
   const Real d_theta = theta_span / n;

   for (Integer i = 0; i < n + 1; i++)
   {
      const Real theta = start_theta + i * d_theta;
      p.x = point.x + radius * cos (theta);
      p.y = point.y + radius * sin (theta);
      add (p);
   }

}

Sector::Sector (const Point_2D& point,
                const Real outer_radius,
                const Real inner_radius,
                const Real start_theta,
                const Real end_theta)
     : Polygon ()
{

   Point_2D p;

   const Real theta_span = end_theta - start_theta;
   const Real proportion = fabs (theta_span) / M_2_TIMES_PI;
   const Integer n = Integer ((log (fabs (outer_radius) + 1) * 18) * proportion) + 5;
   const Real d_theta = theta_span / n;

   for (Integer i = 0; i < n + 1; i++)
   {
      const Real theta = start_theta + i * d_theta;
      p.x = point.x + outer_radius * cos (theta);
      p.y = point.y + outer_radius * sin (theta);
      add (p);
   }

   for (Integer i = n; i >= 0; i--)
   {
      const Real theta = start_theta + i * d_theta;
      p.x = point.x + inner_radius * cos (theta);
      p.y = point.y + inner_radius * sin (theta);
      add (p);
   }

}

Segment::Segment (const Point_2D& point,
                  const Real radius,
                  const Real start_theta,
                  const Real end_theta)
{

   Point_2D p;

   const Real theta_span = end_theta - start_theta;
   const Real proportion = fabs (theta_span) / M_2_TIMES_PI;
   const Integer n = Integer ((log (fabs (radius) + 1) * 18) * proportion) + 5;
   const Real d_theta = theta_span / n;

   for (Integer i = 0; i < n + 1; i++)
   {
      const Real theta = start_theta + i * d_theta;
      p.x = point.x + radius * cos (theta);
      p.y = point.y + radius * sin (theta);
      add (p);
   }

}

Symbol::Symbol ()
        : size (0)
{
}

Symbol::Symbol (const Real size)
        : size (size)
{
}

void
Symbol::add_circle (const bool new_handle,
                    const Point_2D& point,
                    const Real radius,
                    const bool positive_direction)
{

   Point_2D p;

   const Integer n = Integer ((log (fabs (radius) + 1) * 18)) + 5;
   const Real theta_span = M_2_TIMES_PI;
   const Real d_theta = (positive_direction ? 1 : -1) * theta_span / n;

   for (Integer i = 0; i < n; i++)
   {
      const Real theta = i * d_theta;
      p.x = point.x + radius * cos (theta);
      p.y = point.y + radius * sin (theta);
      add (p, (i == 0) && (new_handle));
   }

}

void
Symbol::add_arc (const bool new_handle,
                 const Real start_theta,
                 const Real end_theta)
{
   return add_arc (new_handle, Point_2D (0, 0), size, start_theta, end_theta);
}

void
Symbol::add_arc (const bool new_handle,
                 const Point_2D& point,
                 const Real radius,
                 const Real start_theta,
                 const Real end_theta)
{

   Point_2D p;

   const Real theta_span = end_theta - start_theta;
   const Real proportion = fabs (theta_span) / M_2_TIMES_PI;
   const Integer n = Integer ((log (fabs (radius) + 1) * 18) * proportion) + 5;
   const Real d_theta = theta_span / n;

   for (Integer i = 0; i < n + 1; i++)
   {
      const Real theta = start_theta + i * d_theta;
      p.x = point.x + radius * cos (theta);
      p.y = point.y + radius * sin (theta);
      add (p, (i == 0) && (new_handle));
   }

}

void
Symbol::add_sector (const bool new_handle,
                    const Real start_theta,
                    const Real end_theta)
{
   add (Point_2D (0, 0), new_handle);
   add_arc (false, start_theta, end_theta);
}

void
Symbol::add_sector (const bool new_handle,
                    const Point_2D& point,
                    const Real radius,
                    const Real start_theta,
                    const Real end_theta)
{
   add (point, new_handle);
   add_arc (false, point, radius, start_theta, end_theta);
}

Ring::Ring (const Real size)
  : Symbol (size)
{
   add_circle (true, Point_2D (0, 0), size);
}

void
Arrow::init (const Real theta,
             const Real size,
             const Real width_ratio)
{

   Real rs[7];
   Real thetas[7];

   const Real s = size / 2;
   const Real w = (!gsl_finite (width_ratio) ? size / 6 : size * width_ratio);

   thetas[0] = 0;
   thetas[1] = atan2 (2 * w, s - 3 * w);
   thetas[2] = atan2 (w, s - 3 * w);

   thetas[3] = atan2 (w, -s);
   thetas[4] = atan2 (-w, -s);
   thetas[5] = atan2 (-w, s - 3 * w);
   thetas[6] = atan2 (-2 * w, s - 3 * w);

   rs[0] = s;
   rs[1] = sqrt ((s - 3 * w) * (s - 3 * w) + 4 * w * w);
   rs[2] = sqrt ((s - 3 * w) * (s - 3 * w) + w * w);
   rs[3] = sqrt (s * s + w * w);
   rs[4] = sqrt (s * s + w * w);
   rs[5] = sqrt ((s - 3 * w) * (s - 3 * w) + w * w);
   rs[6] = sqrt ((s - 3 * w) * (s - 3 * w) + 4 * w * w);


   Point_2D vertex;

   for (Integer i = 0; i < 7; i++)
   {
      vertex.x = rs[i] * cos (thetas[i] + theta),
      vertex.y = rs[i] * sin (thetas[i] + theta);
      add (vertex);
   }

}

Arrow::Arrow (const Real theta,
              const Real size,
              const Real width_ratio)
    : Symbol (size)
{
   init (theta, size, width_ratio);
}

Arrow::Arrow (const Point_2D& point_from,
              const Point_2D& point_to,
              const Real width)
{

   const Real dx = (point_to.x - point_from.x);
   const Real dy = (point_to.y - point_from.y);
   const Real size = sqrt (dx * dx + dy * dy);
   const Real theta = atan2 (dy, dx);

   const Point_2D anchor (point_from.x + dx / 2, point_from.y + dy / 2);
   const Real width_ratio = width / size;

   Symbol::size = size;
   init (theta, size, width_ratio);
   translate (anchor);

}

Arrow
Arrow::left_arrow (const Real size,
                   const Real width_ratio)
{
   return Arrow (M_PI, size, width_ratio);
}

Arrow
Arrow::right_arrow (const Real size,
                    const Real width_ratio)
{
   return Arrow (0, size, width_ratio);
}

Arrow
Arrow::up_arrow (const Real size,
                 const Real width_ratio)
{
   return Arrow (M_PI_2, size, width_ratio);
}

Arrow
Arrow::down_arrow (const Real size,
                   const Real width_ratio)
{
   return Arrow (-M_PI_2, size, width_ratio);
}

Cross::Cross (const Real size,
              const Real width)
    : Symbol (size)
{

   const Real a = width / M_SQRT1_2;
   const Real b = size / M_SQRT1_2;
   const Real c = a + b;

   add (Point_2D (0, a));
   add (Point_2D (b, c));
   add (Point_2D (c, b));

   add (Point_2D (a, 0));
   add (Point_2D (c, -b));
   add (Point_2D (b, -c));

   add (Point_2D (0, -a));
   add (Point_2D (-b, -c));
   add (Point_2D (-c, -b));

   add (Point_2D (-a, 0));
   add (Point_2D (-c, b));
   add (Point_2D (-b, c));

}

Plus::Plus (const Real size,
            const Real width)
  : Symbol (size)
{

   const Real a = width / 2;
   const Real b = size;

   add (Point_2D (a, b));
   add (Point_2D (a, a));
   add (Point_2D (b, a));

   add (Point_2D (b, -a));
   add (Point_2D (a, -a));
   add (Point_2D (a, -b));

   add (Point_2D (-a, -b));
   add (Point_2D (-a, -a));
   add (Point_2D (-b, -a));

   add (Point_2D (-b, a));
   add (Point_2D (-a, a));
   add (Point_2D (-a, b));

}

Tee::Tee (const Real size,
          const Real width)
  : Symbol (size)
{

   const Real a = width / 2;
   const Real b = size;

   add (Point_2D (a, b));
   add (Point_2D (a, a));
   add (Point_2D (b / 2, a));

   add (Point_2D (b / 2, -a));
   add (Point_2D (-b / 2, -a));

   add (Point_2D (-b / 2, a));
   add (Point_2D (-a, a));
   add (Point_2D (-a, b));

}

Square::Square (const Real size)
      : Symbol (size)
{
   add (Point_2D (-size, -size));
   add (Point_2D (-size, size));
   add (Point_2D (size, size));
   add (Point_2D (size, -size));
}

Triangle::Triangle (const Real size,
                    const Real angle)
   : Symbol (size)
{

   const Real d_theta = M_2_TIMES_PI / 3;

   for (Integer i = 0; i < 3; i++)
   {
      const Real theta = angle + i * d_theta;
      const Real x = size * cos (theta);
      const Real y = size * sin (theta);
      add (Point_2D (x, y));
   }

}

Star::Star (const Real size)
   : Symbol (size)
{

   const Real radian_036 = 36 * DEGREE_TO_RADIAN;
   const Real radian_072 = 72 * DEGREE_TO_RADIAN;
   const Real radian_128 = 128 * DEGREE_TO_RADIAN;

   const Real O_309 = cos (radian_072);
   const Real O_809 = cos (radian_036);
   const Real O_587 = sin (radian_036);
   const Real O_951 = sin (radian_072);
   const Real O_788 = sin (radian_128);

   const Real size_l = size;
   const Real size_s = size * O_788 / O_309;

   add (Point_2D (0, size_l));
   add (Point_2D (size_s * O_587, size_s * O_809));
   add (Point_2D (size_l * O_951, size_l * O_309));
   add (Point_2D (size_s * O_951, -size_s * O_309));
   add (Point_2D (size_l * O_587, -size_l * O_809));
   add (Point_2D (0, -size_s));
   add (Point_2D (-size_l * O_587, -size_l * O_809));
   add (Point_2D (-size_s * O_951, -size_s * O_309));
   add (Point_2D (-size_l * O_951, size_l * O_309));
   add (Point_2D (-size_s * O_587, size_s * O_809));

}

Bowtie::Bowtie (const Real size)
   : Symbol (size)
{

   const Real s = size / M_SQRT2;

   add (Point_2D (s, -s));
   add (Point_2D (s, s));
   add (Point_2D (-s, -s));
   add (Point_2D (-s, s));

}

Cat_Head::Cat_Head (const Real size)
   : Symbol (size)
{

   const Real n = 360;
   const Real a = 37;
   const Real b = 33;
   const Real c = 10;
   const Real d = 28;

   const Real scale = size / a;

   const Real aa = a*a;
   const Real bb = b*b;
   const Real c_hat = sqrt (bb * (1 - c*c / aa));
   const Real d_hat = sqrt (bb * (1 - d*d / aa));
   const Real alpha = atan (c_hat / c);
   const Real beta = atan (d_hat / d);
   const Real gamma = alpha * 0.55 + beta * 0.45;
   const Real d_hat_hat = d * tan (gamma);

   const Real angle_rl = alpha;
   const Real angle_lr = M_PI - alpha;
   const Real angle_la = M_PI - gamma;
   const Real angle_ll = M_PI - beta;
   const Real angle_rr = M_2_TIMES_PI + beta;
   const Real angle_ra = M_2_TIMES_PI + gamma;

   const Real scalp_span = (M_PI - 2*alpha);
   const Real chin_span = (M_PI + 2*beta);

   const Integer n_scalp = Integer (Real ((scalp_span/M_2_TIMES_PI)*n) + 1);
   const Integer n_chin = Integer (Real ((chin_span/M_2_TIMES_PI)*n) + 1);
   const Real dt_scalp = scalp_span / (n_scalp - 1);
   const Real dt_chin = chin_span / (n_chin - 1);

   const Real major_radius = a * scale;
   const Real minor_radius = b * scale;
   const Real ear_x = d * scale;
   const Real ear_y = d_hat_hat * scale;

   for (Integer i = 0; i < n_scalp; i++)
   {
      const Real t = -(angle_rl + i * dt_scalp);
      add (Point_2D (major_radius * cos (t), minor_radius * sin (t)));
   }

   add (Point_2D (-ear_x, -ear_y));

   for (Integer i = 0; i < n_chin; i++)
   {
      const Real t = -(angle_ll + i * dt_chin);
      add (Point_2D (major_radius * cos (t), minor_radius * sin (t)));
   }

   add (Point_2D (ear_x, -ear_y));

}

Voronoi_Event::Voronoi_Event (const Voronoi_Event_Type type,
                              const Point_2D& point_2d)
   : type (type),
     Point_2D (point_2d)
{
}

bool
Voronoi_Event::operator < (const Point_2D& point_2d) const
{
   if (y == point_2d.y) { return x < point_2d.x; }
   else { return y < point_2d.y; }
}

bool
Voronoi_Event::operator == (const Point_2D& point_2d) const
{
   return (x == point_2d.x && y == point_2d.y);
}

bool
Voronoi_Event::operator > (const Point_2D& point_2d) const
{
   if (y == point_2d.y) { return x > point_2d.x; }
   else { return y > point_2d.y; }
}

Voronoi_Event_Queue::Voronoi_Event_Queue (const vector<Point_2D>& point_vector)
{

   for (vector<Point_2D>::const_iterator iterator = point_vector.begin ();
        iterator != point_vector.end (); iterator++)
   {
      const Point_2D& point_2d = *(iterator);
      Voronoi_Event event (VORONOI_SITE, point_2d);
      push (event);
   }

}

Real
Voronoi_Breakpoint::get_x (const Point_2D& site_0,
                           const Point_2D& site_1,
                           const Real y_b) const
{

   double x_0, x_1;
   Real d_0 = site_0.y - y_b;
   Real d_1 = site_1.y - y_b;
   Real y_b_2 = y_b * y_b;

   Real a_0 = 1 / (2 * d_0);
   Real a_1 = 1 / (2 * d_1);

   Real b_0 = - site_0.x / d_0;
   Real b_1 = - site_1.x / d_1;

   Real c_0 = (site_0.x*site_0.x + site_0.y*site_0.y) - y_b_2;
   Real c_1 = (site_1.x*site_1.x + site_1.y*site_1.y) - y_b_2;

   gsl_poly_solve_quadratic (a_0 - a_1, b_0 - b_1, c_0 - c_1, &x_0, &x_1);
   if (x_0 > x_1) { std::swap (x_0, x_1); }

   switch (breakpoint_genre)
   {
      case VORONOI_LEFT:  return x_0;
      case VORONOI_RIGHT: return x_1;
   }

}

Voronoi_Breakpoint::Voronoi_Breakpoint (const Voronoi_Breakpoint_Genre breakpoint_genre)
                    : breakpoint_genre (breakpoint_genre)
{
}

Voronoi_Arc::Voronoi_Arc (const Point_2D& site)
              : Point_2D (site),
               event_ptr (NULL),
        breakpoint_l_ptr (NULL),
        breakpoint_r_ptr (NULL)
{
}

Voronoi_Arc::Voronoi_Arc (const Voronoi_Arc& arc)
              : Point_2D (dynamic_cast<const Point_2D&> (arc)),
               event_ptr (arc.event_ptr),
        breakpoint_l_ptr (arc.breakpoint_l_ptr),
        breakpoint_r_ptr (arc.breakpoint_r_ptr)
{
}

Integer
Voronoi_Beachline::arc_match (const Real x,
                              const Integer index) const
{

   if (size () == 1) { return true; }
   else
   {

      const Voronoi_Arc& arc = at (index);
      const Point_2D& site = dynamic_cast<const Point_2D&> (arc);

      if (index == 0)
      {
         const Voronoi_Arc& arc_r = at (index + 1);
         const Voronoi_Breakpoint& breakpoint_r = *(arc.breakpoint_r_ptr);
         const Point_2D& site_r = dynamic_cast<const Point_2D&> (arc_r);
         Real x_r = breakpoint_r.get_x (site, site_r, y_b);
         if (x < x_r) { return true; }
      }
      else
      if (index == size () - 1)
      {
         const Voronoi_Arc& arc_l = at (index - 1);
         const Voronoi_Breakpoint& breakpoint_l = *(arc.breakpoint_l_ptr);
         const Point_2D& site_l = dynamic_cast<const Point_2D&> (arc_l);
         Real x_l = breakpoint_l.get_x (site_l, site, y_b);
         if (x >= x_l) { return true; }
      }
      else
      {
         const Voronoi_Arc& arc_l = at (index - 1);
         const Voronoi_Arc& arc_r = at (index + 1);
         const Voronoi_Breakpoint& breakpoint_l = *(arc.breakpoint_l_ptr);
         const Voronoi_Breakpoint& breakpoint_r = *(arc.breakpoint_r_ptr);
         const Point_2D& site_l = dynamic_cast<const Point_2D&> (arc_l);
         const Point_2D& site_r = dynamic_cast<const Point_2D&> (arc_r);
         Real x_l = breakpoint_l.get_x (site_l, site, y_b);
         Real x_r = breakpoint_r.get_x (site, site_r, y_b);
         if (x >= x_l && x < x_r) { return true; }
      }

   }

   return false;

}

Integer
Voronoi_Beachline::search_index (const Real x,
                                 const Integer start_index,
                                 const Integer end_index) const
{

   const Integer index = (start_index + end_index) / 2;
   const Integer match = arc_match (x, index);

   switch (match)
   {
      case 0:  return index;
      case -1: return search_index (x, start_index, index - 1);
      case +1: return search_index (x, index + 1, end_index);
   }

}

Voronoi_Beachline::iterator
Voronoi_Beachline::search_iterator (const Real x)
{
   const Integer index = search_index (x, 0, size () - 1);
   Voronoi_Beachline::iterator iterator = begin ();
   return (iterator + index);
}

Voronoi_Beachline::Voronoi_Beachline (const Real& y_b)
                               : y_b (y_b)
{
}

bool
Voronoi_Beachline::insert_site (const Point_2D& site)
{

   Voronoi_Arc arc (site);

   if (size () == 0) { push_back (arc); }
   else
   {

      Voronoi_Breakpoint* breakpoint_l_ptr;
      Voronoi_Breakpoint* breakpoint_r_ptr;

      breakpoint_l_ptr = new Voronoi_Breakpoint (VORONOI_LEFT);
      breakpoint_r_ptr = new Voronoi_Breakpoint (VORONOI_RIGHT);

      Voronoi_Beachline::iterator iterator = search_iterator (site.x);
      Voronoi_Arc& arc_r = *(iterator);
      Voronoi_Arc arc_l (arc_r);

      arc_l.breakpoint_r_ptr = breakpoint_l_ptr;
      arc_r.breakpoint_l_ptr = breakpoint_r_ptr;

      arc.breakpoint_l_ptr = breakpoint_l_ptr;
      arc.breakpoint_r_ptr = breakpoint_r_ptr;

      insert (iterator, arc);
      insert (iterator, arc_l);

   }

}

/*
void
Voronoi_Diagram::handle_event (const Voronoi_Event& event)
{
   switch (event.type)
   {
      case VORONOI_SITE:   handle_site_event (event); break;
      case VORONOI_CIRCLE: handle_circle_event (event); break;
   }
}

void
Voronoi_Diagram::handle_site_event (const Voronoi_Event& event)
{

   // Insert Site into Beachline

   // Create new half-edge records for edge separating p_i and p_j

   // Check the triple of consecutive arcs where the new arc for p_i
   // is the left arc to see if the breakpoints converge. If so,
   // insert the circle event into Q and add pointers between the node
   // in beachline and the node in event queue. Do the same for
   // the triple where the new arc is the right arc

}

void
Voronoi_Diagram::handle_circle_event (const Voronoi_Event& event)
{
}
*/

namespace denise
{

/*
   ostream&
   operator << (ostream& out_file,
                const Path_Ob& path_ob)
   {

      switch (path_ob.path_ob_genre)
      {

         case MOVE_TO:
            out_file << "[MOVE_TO " << path_ob[0] << "]";
            break;

         case LINE_TO:
            out_file << "[LINE_TO " << path_ob[0] << "]";
            break;

         case ARC_TO:
            out_file << "[ARC_TO " << path_ob[0] << " / " << path_ob[1].x
               << " FROM " << path_ob[2].x << " TO " << path_ob[2].y << "]";
            break;

         case CURVE_TO:
            out_file << "[CURVE_TO " << path_ob[0] << " TO " <<
               path_ob[1] << " TO " << path_ob[2] << "]";
            break;

         case CLOSE_PATH:
            out_file << "[CLOSE_PATH]";
            break;

      }

      return out_file;

   }
*/

   ostream&
   operator << (ostream& out_file,
                const Ellipse& ellipse)
   {
      Point_2D center = ellipse.get_center ();
      out_file << "(" << center << ", " << ellipse.get_a () << ", "
               << ellipse.get_b () << ", " << ellipse.get_tilt () << ")";
      return out_file;
   }

   ostream&
   operator << (ostream& out_file,
                const Simple_Polyline& simple_polyline)
   {
      for (Simple_Polyline::const_iterator iterator = simple_polyline.begin ();
           iterator != simple_polyline.end (); iterator++)
      {
         const Point_2D& point = *(iterator);
         if (iterator != simple_polyline.begin ()) { out_file << " "; }
         out_file << point;
      }
      return out_file;
   }

   ostream&
   operator << (ostream& out_file,
                const Polygon& polygon)
   {

      if (polygon.size () == 0) { return out_file; }

      const Polygon::Vertex* first_handle_ptr = polygon.get_first_handle_ptr ();
      Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)first_handle_ptr;

      do
      {

         Integer n = current_handle_ptr->n;
         Polygon::Vertex* current_ptr = (Polygon::Vertex*)current_handle_ptr;

         for (Integer i = 0; i < n; i++)
         {
            const Point_2D& point = (const Point_2D&)(*current_ptr);
            if (i == 0) { out_file << "N "; }
            out_file << point << " ";
            current_ptr = current_ptr->next_ptr;
         }

         current_handle_ptr = current_handle_ptr->next_handle_ptr;

      }
      while (current_handle_ptr != first_handle_ptr);

      return out_file;

   }

}

