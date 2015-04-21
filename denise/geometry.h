//
// geometry.h
// 
// Copyright (C) 2005-2013 Simon E. Ching
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

#ifndef DENISE_GEOMETRY_H
#define DENISE_GEOMETRY_H

#include <cmath>
#include <map>
#include <queue>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_ellint.h>
#include <cairomm/context.h>
#include <denise/basics.h>
#include <denise/cairoable.h>
#include <denise/exception.h>
#include <denise/linalg.h>
#include <denise/util.h>

using namespace std;
using namespace Cairo;

namespace denise
{

   class Polyline;
   class Polygon;
   class Transform_2D;

//   enum Path_Ob_Genre
//   {
//      MOVE_TO,
//      LINE_TO,
//      CURVE_TO,
//      ARC_TO,
//      CLOSE_PATH
//   };

   enum Boolean_Op
   {
      DIFFERENCE,
      INTERSECTION,
      UNION
   };

   enum Polygon_Intersection
   {
      INTERSECTION_NOT,
      INTERSECTION_ENTRY,
      INTERSECTION_EXIT,
      INTERSECTION_SELF
   };

   enum Polygon_Relation
   {
      A_ENTIRELY_IN_B,
      B_ENTIRELY_IN_A,
      DISJOINT
   };

   enum Polygon_Action
   {
      NOT_DECIDED,
      APPEND,
      DO_NOT_APPEND
   };

   class Geometry_Exception : public Exception
   {

      public:

         Geometry_Exception (const string& description = string (""));

   };

   class Attractor
   {

      private:

         const Real
         grid_x;

         const Real
         grid_y;

      public:

         Attractor ();

         Attractor (const Real grid_x);

         Attractor (const Real grid_x,
                    const Real grid_y);

         virtual void
         attract (Real& x,
                  Real& y) const;

         virtual void
         attract (Point_2D& point) const;

         virtual Point_2D
         get_attraction (const Point_2D& point) const;

   };

   /// A Line Segment connecting two Point_2D's
   class Edge : public Cairoable
   {

      private:

         Real
         get_determinant (const Real a00,
                          const Real a01,
                          const Real a10,
                          const Real a11) const;

      public:

         Point_2D
         point_a;

         Point_2D
         point_b;

         Edge () { }

         Edge (const Edge& edge);

         Edge (const Point_2D& point_a,
               const Point_2D& point_b);

         virtual void
         swap ();

         virtual void
         translate (const Real delta_x,
                    const Real delta_y);

         virtual void
         translate (const Point_2D& point);

         Real
         get_root_x (const Real y) const;

         Real
         get_root_y (const Real x) const;

         static Real
         get_wec (const Point_2D& point,
                  const Point_2D& point_a,
                  const Point_2D& point_b);

         static void
         list_polygon (const Polygon& polygon);

         static bool
         perturb_if_degenerate (Point_2D& point_aa,
                                Point_2D& point_ab,
                                Point_2D& point_ba,
                                Point_2D& point_bb,
                                const Real epsilon);

         static bool
         intersection (Point_2D& intersection,
                       Real& alpha_a,
                       Real& alpha_b,
                       const Point_2D& point_aa,
                       const Point_2D& point_ab,
                       const Point_2D& point_ba,
                       const Point_2D& point_bb);

         static bool
         intersection (Point_2D& intersection,
                       const Point_2D& point_aa,
                       const Point_2D& point_ab,
                       const Point_2D& point_ba,
                       const Point_2D& point_bb);

         static Point_2D
         intersection_between (const Point_2D& point_aa,
                               const Point_2D& point_ab,
                               const Point_2D& point_ba,
                               const Point_2D& point_bb);

         static Point_2D
         intersection_between (const Edge& edge_a,
                               const Edge& edge_b);

         Point_2D
         get_intersection_with (const Edge& edge) const;

         bool
         intersects_with (const Edge& edge) const;

         Tuple
         get_intersection_tuple (const Polygon& polygon) const;

         static Real
         distance_between (const Point_2D& point_a,
                           const Point_2D& point_b);

         static Real
         distance (const Transform_2D& transform_2d,
                   const Point_2D& edge_point_a,
                   const Point_2D& edge_point_b,
                   const Point_2D& point);

         static Real
         distance (const Point_2D& edge_point_a,
                   const Point_2D& edge_point_b,
                   const Point_2D& point);

         Real
         distance (const Transform_2D& transform_2d,
                   const Point_2D& point) const;

         Real
         distance (const Point_2D& point) const;

         Real
         length () const;

         bool
         contains (const Point_2D& point_2d) const;

         Polyline*
         clip (const Polygon& clip_polygon,
               const bool negative = false) const;

         void
         cairo (const RefPtr<Context>& cr) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Polygon& clip_polygon,
                const bool negative = false) const;

   };

   class Area
   {

      public:

         virtual bool
         contains (const Point_2D& point,
                   const bool border_included = false) const = 0;

   };

   /// Circle
   class Circle : public Area,
                  public Cairoable
   {

      private:

         void
         init (const Real A,
               const Real D,
               const Real E,
               const Real F);

      public:

         Point_2D
         center;

         Real
         radius;

         Circle ();

         Circle (const Point_2D& center,
                 const Real radius);

         Circle (const Real A,
                 const Real D,
                 const Real E,
                 const Real F);

         Circle (const Point_2D& point_a,
                 const Point_2D& point_b,
                 const Point_2D& point_c);

         void
         acquire_coefficients (Real& A,
                               Real& D,
                               Real& E,
                               Real& F) const;

         virtual void
         translate (const Real delta_x,
                    const Real delta_y);

         virtual void
         translate (const Point_2D& point);

         void
         scale (const Real scale);

         bool
         contains (const Point_2D& point,
                   const bool border_included = false) const;

         const Point_2D&
         get_center () const;

         Real
         get_radius () const;

         Real
         get_area () const;

         Real
         get_perimeter () const;

         Point_2D
         get_point_at_theta (const Real theta) const;

         pair<Point_2D, Point_2D>
         get_tangent_point_pair (const Point_2D& point) const;

         bool
         is_nac () const;

         void
         cairo (const RefPtr<Context>& cr) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d) const;

   };

   /// An Ellipse
   class Ellipse : public Cairoable
   {

      private:

         void
         init (const Real A,
               const Real B,
               const Real C,
               const Real D,
               const Real E,
               const Real F);

      public:

         Point_2D
         center;

         Real
         a;

         Real
         b;

         Real
         tilt;

         Ellipse ();

         Ellipse (const Point_2D& center,
                  const Real major_radius,
                  const Real minor_radius,
                  const Real slant);

         Ellipse (const Real A,
                  const Real B,
                  const Real C,
                  const Real F);

         Ellipse (const Real A,
                  const Real B,
                  const Real C,
                  const Real D,
                  const Real E,
                  const Real F);

         void
         acquire_coefficients (Real& A,
                               Real& B,
                               Real& C,
                               Real& D,
                               Real& E,
                               Real& F) const;

         virtual void
         translate (const Real delta_x,
                    const Real delta_y);

         virtual void
         translate (const Point_2D& point);

         void
         scale (const Real scale);

         void
         rotate (const Real tilt);

         Point_2D
         get_center () const;

         Real
         get_a () const;

         Real
         get_b () const;

         Real
         get_tilt () const;

         Real
         get_area () const;

         Real
         get_perimeter (const gsl_mode_t gsl_mode = GSL_PREC_APPROX) const;

         Point_2D
         get_point_at_t (const Real t) const;

         Point_2D
         get_point_at_theta (const Real theta) const;

         pair<Point_2D, Point_2D>
         get_tangent_point_pair (const Point_2D& point) const;

         bool
         is_nae () const;

         void
         cairo (const RefPtr<Context>& cr) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d) const;

   };

   class Simple_Polyline : public list<Point_2D>,
                           public Cairoable
   {

      protected:

         bool
         closed;

         bool
         is_last (Simple_Polyline::iterator iterator) const;

         bool
         is_last (Simple_Polyline::const_iterator iterator) const;

      public:

         Simple_Polyline (const bool closed = false);

         Simple_Polyline (const Simple_Polyline& simple_polyline);

         void
         set (const Simple_Polyline& simple_polyline);

         virtual void
         attract_by (const Attractor& attractor);

         bool
         is_closed () const;

         void
         set_closed (const bool closed);

         Edge
         get_edge (Simple_Polyline::iterator iterator) const;

         Edge
         get_edge (Simple_Polyline::const_iterator iterator) const;

         void
         add (const Point_2D& point);

         void
         prepend (const Point_2D& point);

         virtual void
         translate (const Real dx,
                    const Real dy);

         void
         transform (const Transform_2D& transform);

         void
         reverse (const Transform_2D& transform);

         Polyline*
         clip (const Polygon& clip_polygon,
               const bool negative = false) const;

         void
         clip (Polyline& polyline,
               const Polygon& clip_polygon,
               const bool negative = false) const;

         void
         cairo (const RefPtr<Context>& cr) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Polygon& clip_polygon,
                const bool negative = false) const;

         Simple_Polyline::iterator
         get_iterator (const Transform_2D& transform_2d,
                       const Point_2D& point_2d,
                       const Real threshold);

         Simple_Polyline::const_iterator
         get_iterator (const Transform_2D& transform_2d,
                       const Point_2D& point_2d,
                       const Real threshold) const;

         Simple_Polyline::iterator
         get_iterator (const Point_2D& point_2d,
                       const Real threshold);

         Simple_Polyline::const_iterator
         get_iterator (const Point_2D& point_2d,
                       const Real threshold) const;

         Simple_Polyline::iterator
         implant (const Transform_2D& transform_2d,
                  const Point_2D& point_2d,
                  const Real threshold);

         Simple_Polyline::iterator
         implant (const Point_2D& point_2d,
                  const Real threshold);

   };

   class Polyline : public list<Simple_Polyline*>,
                    public Cairoable
   {

      private:

         void
         init (const Point_2D& point);

      public:

         Polyline ();

         Polyline (const Simple_Polyline* simple_polyline_ptr);

         virtual void
         attract_by (const Attractor& attractor);

         void
         add (const Simple_Polyline* simple_polyline_ptr);

         void
         add (const Point_2D& point,
              const bool new_handle = false);

         void
         prepend (const Point_2D& point);

         void
         transform (const Transform_2D& transform);

         void
         reverse (const Transform_2D& transform);

         Polyline*
         clip (const Polygon& clip_polygon,
               const bool negative = false) const;

         void
         cairo (const RefPtr<Context>& cr) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d) const;

   };

   class Polygon_Vertex : public Point_2D
   {

      public:

         Polygon_Action
         action;

         Integer
         n;

         Real
         alpha;

         Integer
         intersections;

         Polygon_Vertex*
         next_ptr;

         Polygon_Vertex*
         prev_ptr;

         Polygon_Vertex*
         neighbor_ptr;

         Polygon_Vertex*
         handle_ptr;

         Polygon_Vertex*
         next_handle_ptr;

         Polygon_Vertex*
         prev_handle_ptr;

         Polygon_Intersection
         intersection;

         Polygon_Vertex (const Point_2D& point);

   };

   /// a list of point_2d's that closes
   class Polygon : public Area,
                   public Cairoable
   {

      private:

         Polygon_Vertex*
         first_handle_ptr;

         Real
         min_edge_length;

         Domain_1D
         domain_x;

         Domain_1D
         domain_y;

         static Real
         determinant (const Point_2D& point,
                      const Point_2D& point_i,
                      const Point_2D& point_ip1);

         static bool
         cross (const Point_2D& point,
                const Point_2D& point_i,
                const Point_2D& point_ip1);

         static bool
         cross_r (const Point_2D& point,
                  const Point_2D& point_i,
                  const Point_2D& point_ip1);

         static void
         modify (Integer& winding_number,
                 const Point_2D& point_i,
                 const Point_2D& point_ip1);

         static void
         iterate_winding_number (Integer& w,
                                 const Polygon_Vertex* handle_ptr,
                                 const Point_2D& point);

         static Integer
         get_winding_number (const Polygon_Vertex* handle_ptr,
                             const Point_2D& point);

         static bool
         contains (const Polygon_Vertex* handle_ptr,
                   const Point_2D& point,
                   bool border_included = false);

         static void
         toggle_intersection (Polygon_Intersection& intersection);

         Polygon_Vertex*
         get_intersection_ptr () const;

         static Point_2D
         get_intersection (Point_2D& point_aa,
                           Point_2D& point_ab,
                           Point_2D& point_ba,
                           Point_2D& point_bb,
                           Real perturbation);

         static void
         boolean_op_phase_0 (const Polygon& polygon_a,
                             const Polygon& polygon_b,
                             const Real perturbation);

         static void
         boolean_op_phase_a (const Polygon& polygon_a,
                             const Polygon& polygon_b);

         static void
         boolean_op_phase_b (const Polygon& polygon_a,
                             const Polygon& polygon_b,
                             const bool flip);

         static void
         boolean_op_phase_c (Polygon& polygon,
                             const Polygon& polygon_a,
                             const Polygon& polygon_b,
                             const Boolean_Op boolean_op);

         static void
         boolean_op_phase_c_tilde (Polygon& polygon,
                                   const Polygon& polygon_a,
                                   const Polygon& polygon_b,
                                   const Boolean_Op boolean_op);

         static void
         simplify_phase_a (Polygon& polygon);

         static void
         simplify_phase_b (Polygon& polygon);

         static Polygon_Relation
         get_polygon_relation (const Polygon_Vertex* handle_a_ptr,
                               const Polygon_Vertex* handle_b_ptr);

         static bool
         entirely_within (const Polygon_Vertex* handle_a_ptr,
                          const Polygon_Vertex* handle_b_ptr);

         bool
         append (Polygon_Vertex* handle_ptr);

      public:

         Polygon ();

         Polygon (const Polygon& polygon);

         Polygon (const Circle& circle,
                  const Integer n = 360);

         Polygon (const Ellipse& ellipse,
                  const Integer n = 360);

         ~Polygon ();

         void
         operator= (const Polygon& polygon);

         void
         clear ();

         virtual void
         attract_by (const Attractor& attractor);

         void
         add (const Point_2D& point,
              const bool new_handle = false);

         void
         add (const Polygon& polygon);

         void
         add (const Polygon& polygon,
              const Point_2D& offset);

         void
         add (const Simple_Polyline& simple_polyline,
              const bool new_handle = false);

         void
         transform (const Transform_2D& transform);

         void
         reverse (const Transform_2D& transform);

         void
         simplify ();

         virtual void
         translate (const Point_2D& point);

         static void
         remove (Polygon_Vertex* vertex_ptr);

         Integer
         size () const;

         Integer
         get_number_of_single_polygons () const;

         const Polygon_Vertex*
         get_first_handle_ptr () const;

         const Polygon_Vertex*
         get_last_handle_ptr () const;

         const Domain_1D&
         get_domain_x () const;

         const Domain_1D&
         get_domain_y () const;

         static Real
         get_simple_polygon_area (const Polygon_Vertex* pv_ptr);

         Real
         get_area () const;

         Real
         get_positive_area () const;

         Real
         get_negative_area () const;

         Point_2D
         get_centroid () const;

         Integer
         get_winding_number (const Real x,
                             const Real y) const;

         Integer
         get_winding_number (const Point_2D& point) const;

         bool
         contains (const Point_2D& point,
                   const bool border_included = false) const;

         static void
         boolean_op (Polygon& polygon,
                     const Boolean_Op op,
                     const Polygon& polygon_a,
                     const Polygon& polygon_b,
                     Real perturbation = GSL_NAN);

         static Polygon*
         boolean_op (const Boolean_Op op,
                     const Polygon& polygon_a,
                     const Polygon& polygon_b,
                     Real perturbation = GSL_NAN);

         void
         cairo (const RefPtr<Context>& cr) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Point_2D& offset) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d) const;

         static void
         cairo (const RefPtr<Context>& cr,
                const Polygon_Vertex* polygon_vertex_ptr);

         void
         debug_print (const string& prefix = "",
                      ostream& out_stream = cout) const;

   };

   class Parametric : public vector<Tuple>
   {

      public:

         virtual Real
         get_x (const Real t) const;

         virtual Real
         get_y (const Real t) const = 0;

         virtual Point_2D
         get_point (const Real t) const;

         virtual Real
         get_weight (const Real t) const;

         Tuple
         get_t_tuple (const Integer n,
                      const Real start,
                      const Real end,
                      const Real w0 = 0.1) const;

         Simple_Polyline*
         get_simple_polyline_ptr (const Integer n,
                                  const Real start,
                                  const Real end,
                                  const Real w0 = 0.1) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Transform_2D& transform_2d,
                const Integer n,
                const Real start,
                const Real end,
                const Real w0 = 0.1) const;

   };

   class Polar_Parametric : public Parametric
   {

      protected:

         Point_2D
         centre;

      public:

         Polar_Parametric (const Point_2D& centre = Point_2D (0, 0));

         virtual Real
         get_r (const Real theta) const = 0;

         Real
         get_x (const Real theta) const;

         Real
         get_y (const Real theta) const;

         Point_2D
         get_point (const Real theta) const;

   };

   class Log_Spiral : public Polar_Parametric
   {

      protected:

         Real
         a;

         Real
         b;

      public:

         Log_Spiral (const Real a,
                     const Real b,
                     const Point_2D& centre = Point_2D (0, 0));

         Log_Spiral (const Point_2D& centre,
                     const Point_2D& point,
                     const Real pitch);

         Real
         get_r (const Real theta) const;

   };

   class Descartes_Folium : public Polar_Parametric
   {

      protected:

         Real
         a;

      public:

         Descartes_Folium (const Real a,
                           const Point_2D& centre = Point_2D (0, 0));

         Real
         get_r (const Real theta) const;

   };

   class Lemniscate : public Polar_Parametric
   {

      protected:

         Real
         a;

      public:

         Lemniscate (const Real a,
                     const Point_2D& centre = Point_2D (0, 0));

         Real
         get_r (const Real theta) const;

   };

   class Pascal_Limacon : public Polar_Parametric
   {

      protected:

         Real
         a;

         Real
         b;

      public:

         Pascal_Limacon (const Real a,
                         const Real b,
                         const Point_2D& centre = Point_2D (0, 0));

         Real
         get_r (const Real theta) const;

   };

   class Cardoid : public Pascal_Limacon
   {

      public:

         Cardoid (const Real a,
                  const Point_2D& centre = Point_2D (0, 0));

   };

   class Rose : public Polar_Parametric
   {

      protected:

         Real
         a;

         Real
         n;

         Real
         phi;

      public:

         Rose (const Real a,
               const Real n,
               const Real phi = 0,
               const Point_2D& centre = Point_2D (0, 0));

         Real
         get_r (const Real theta) const;

   };

   class Trochoid : public Parametric
   {

      protected:

         Real
         a;

         Real
         b;

         Point_2D
         centre;

      public:

         Trochoid (const Real a,
                   const Real b,
                   const Point_2D& centre = Point_2D (0, 0));

         Real
         get_x (const Real t) const;

         Real
         get_y (const Real t) const;

   };

   class Cycloid : public Trochoid
   {

      public:

         Cycloid (const Real a,
                  const Point_2D& centre = Point_2D (0, 0));

   };

   class Rect : public Polygon
   {

      private:

         void
         translate_anchor (const Real translate_x,
                           const Real translate_y);

         void
         construct_polygon ();

         void
         init ();

         void
         init (const Point_2D& ponint_2d,
               const Real width,
               const Real height,
               const Real translate_x = 0, 
               const Real translate_y = 0, 
               const Real tilt = 0);

      public:

         Point_2D
         point_2d;

         Real
         width;

         Real
         height;

         Real
         tilt;

         Rect ();

         Rect (const Box_2D& box_2d);

         Rect (const Point_2D& point_a,
               const Point_2D& point_b);

         Rect (const Point_2D& point_2d,
               const Real width,
               const Real height,
               const Real tilt = 0);

         Rect (const Point_2D& point_2d,
               const Real width,
               const Real height,
               const Real translate_x,
               const Real translate_y,
               const Real tilt = 0);

         Rect (const Point_2D& point_2d,
               const Real width,
               const Real height,
               const char justify_h,
               const char justify_v,
               const Real tilt = 0);

         Rect (const Rect& rect);

         void
         set (const Rect& rect);

         void
         operator= (const Rect& rect);

         void
         grow (const Real margin_x,
               const Real margin_y);

   };

/*
   class Path_Ob : public vector<Point_2D>
   {

      public:

         Path_Ob_Genre
         path_ob_genre;

         Path_Ob () { }

         Path_Ob (const Path_Ob_Genre path_ob_genre);

         Path_Ob (const Path_Ob_Genre path_ob_genre,
                  const Point_2D& point);

         Path_Ob (const Path_Ob_Genre path_ob_genre,
                  const Point_2D& point_0,
                  const Point_2D& point_1,
                  const Point_2D& point_2);

   };

   /// Path
   class Path : public list<Path_Ob>
   {

      public:

         Path ();

         Path (const Simple_Polyline& simple_polyline);

         Path (const Polyline& polyline);

         Path (const Polygon& polygon);

         void
         simple_polyline (const Simple_Polyline& simple_polyline);

         void
         move_to (const Point_2D& point);

         void
         line_to (const Point_2D& point);

         void
         arc_to (const Point_2D& center,
                 const Real radius,
                 const Real theta_0,
                 const Real theta_1);

         void
         curve_to (const Point_2D& point_0,
                   const Point_2D& point_1,
                   const Point_2D& point_2);

         void
         close ();

   };
*/

   class Arc : public Simple_Polyline
   {

      public:

         Arc (const Point_2D& point,
              const Real radius,
              const Real start_theta,
              const Real end_theta);

   };

   class Chord : public Edge
   {

      public:

         Chord (const Point_2D& point,
                const Real radius,
                const Real start_theta,
                const Real end_theta);

   };

   class Sector : public Polygon
   {

      public:

         Sector (const Point_2D& point,
                 const Real radius,
                 const Real start_theta,
                 const Real end_theta);

         Sector (const Point_2D& point,
                 const Real outer_radius,
                 const Real inner_radius,
                 const Real start_theta,
                 const Real end_theta);

   };

   class Segment : public Polygon
   {

      public:

         Segment (const Point_2D& point,
                  const Real radius,
                  const Real start_theta,
                  const Real end_theta);

   };

   class Symbol : public Polygon
   {

      protected:

         Real
         size;

      public:

         Symbol ();

         Symbol (const Real size);

         void
         add_circle (const bool new_handle,
                     const Point_2D& point,
                     const Real radius,
                     const bool positive_direction = true);

         void
         add_arc (const bool new_handle,
                  const Real start_theta,
                  const Real end_theta);

         void
         add_arc (const bool new_handle,
                  const Point_2D& point,
                  const Real radius,
                  const Real start_theta,
                  const Real end_theta);

         void
         add_sector (const bool new_handle,
                     const Real start_theta,
                     const Real end_theta);

         void
         add_sector (const bool new_handle,
                     const Point_2D& point,
                     const Real radius,
                     const Real start_theta,
                     const Real end_theta);

   };

   class Ring : public Symbol
   {

      public:

         Ring (const Real size);

   };

   class Arrow : public Symbol
   {

      private:

         void
         init (const Real theta,
               const Real size,
               const Real width_ratio);

      public:

         Arrow (const Real theta,
                const Real size,
                const Real width_ratio = GSL_NAN);

         Arrow (const Point_2D& point_from,
                const Point_2D& point_to,
                const Real width = GSL_NAN);

         static Arrow
         left_arrow (const Real size,
                     const Real width_ratio = GSL_NAN);

         static Arrow
         right_arrow (const Real size,
                      const Real width_ratio = GSL_NAN);

         static Arrow
         up_arrow (const Real size,
                   const Real width_ratio = GSL_NAN);

         static Arrow
         down_arrow (const Real size,
                     const Real width_ratio = GSL_NAN);

   };

   class Cross : public Symbol
   {

      public:

         Cross (const Real size,
                const Real width = 1);

   };

   class Plus : public Symbol
   {

      public:

         Plus (const Real size,
               const Real width = 1);

   };

   class Tee : public Symbol
   {

      public:

         Tee (const Real size,
              const Real width = 1);

   };

   class Square : public Symbol
   {

      public:

         Square (const Real size);

   };

   class Triangle : public Symbol
   {

      public:

         Triangle (const Real size,
                   const Real angle = 0);

   };

   class Star : public Symbol
   {

      public:

         Star (const Real size);

   };

   class Bowtie : public Symbol
   {

      public:

         Bowtie (const Real size);

   };

   class Cat_Head : public Symbol
   {

      public:

         Cat_Head (const Real size);

   };

   enum Voronoi_Event_Type
   {
      VORONOI_SITE,
      VORONOI_CIRCLE
   };

   class Voronoi_Event : public Point_2D
   {

      public:

         Voronoi_Event_Type
         type;

         Voronoi_Event (const Voronoi_Event_Type type,
                        const Point_2D& point_2d);

         bool
         operator < (const Point_2D& point_2d) const;

         bool
         operator == (const Point_2D& point_2d) const;

         bool
         operator > (const Point_2D& point_2d) const;

   };

   class Voronoi_Event_Queue : public priority_queue<Voronoi_Event>
   {

      public:

         Voronoi_Event_Queue (const vector<Point_2D>& point_vector);

   };

   enum Voronoi_Breakpoint_Genre
   {
      VORONOI_LEFT,
      VORONOI_RIGHT
   };

   class Voronoi_Breakpoint
   {

      private:

         const Voronoi_Breakpoint_Genre
         breakpoint_genre;

      public:

         Voronoi_Breakpoint (const Voronoi_Breakpoint_Genre genre);

         Real
         get_x (const Point_2D& site_0,
                const Point_2D& site_1,
                const Real y_b) const;

   };

   class Voronoi_Arc : public Point_2D
   {

      public:

         Voronoi_Event*
         event_ptr;

         Voronoi_Breakpoint*
         breakpoint_l_ptr;

         Voronoi_Breakpoint*
         breakpoint_r_ptr;

         Voronoi_Arc (const Point_2D& site);

         Voronoi_Arc (const Voronoi_Arc& arc);

   };

   class Voronoi_Beachline : public vector<Voronoi_Arc>
   {

      private:

         const Real&
         y_b;

         Integer
         arc_match (const Real x,
                    const Integer index) const;

         Integer
         search_index (const Real x,
                       const Integer start_index,
                       const Integer end_index) const;

         Voronoi_Beachline::iterator
         search_iterator (const Real x);

      public:

         Voronoi_Beachline (const Real& y_b);

         bool
         insert_site (const Point_2D& site);

   };

/*
   class Voronoi_Diagram
   {

      private:

         Voronoi_Beachline
         beachline;

         Voronoi_Event_Queue
         event_queue;

         void
         handle_event (const Voronoi_Event& event);

         void
         handle_site_event (const Voronoi_Event& event);

         void
         handle_circle_event (const Voronoi_Event& event);

      public:

   };
*/

//   ostream&
//   operator << (ostream& out_file,
//                const Path_Ob& path_ob);

   ostream&
   operator << (ostream& out_file,
                const Ellipse& ellipse);

   ostream&
   operator << (ostream& out_file,
                const Simple_Polyline& simple_polylin);

   ostream&
   operator << (ostream& out_file,
                const Polygon& polygon);

}

#endif /* DENISE_GEOMETRY_H */

