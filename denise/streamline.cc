//
// streamline.cc
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

#include "streamline.h"

using namespace std;
using namespace denise;

/*
Duple_Field_2D::Duple_Field_2D (const Tuple_Field_2D& tuple_field_2d,
                                const Integer tuple_u_index,
                                const Integer tuple_v_index,
                                const Real u_multiplier,
                                const Real v_multiplier)
              : tuple_field_2d (tuple_field_2d),
                 tuple_u_index (tuple_u_index),
                 tuple_v_index (tuple_v_index),
                  u_multiplier (u_multiplier),
                  v_multiplier (v_multiplier)
{
}

Real
Duple_Field_2D::evaluate_u_at (const Real x,
                               const Real y,
                               const Evaluate_Op evaluate_op) const
{
   return tuple_field_2d.evaluate_at (x, y,
      tuple_u_index, evaluate_op) * u_multiplier;
}

Real
Duple_Field_2D::evaluate_u_at (const Point_2D& point,
                               const Evaluate_Op evaluate_op) const
{
   return evaluate_u_at (point.x, point.y, evaluate_op);
}

Real
Duple_Field_2D::evaluate_v_at (const Real x,
                               const Real y,
                               const Evaluate_Op evaluate_op) const
{
   return tuple_field_2d.evaluate_at (x, y,
      tuple_v_index, evaluate_op) * v_multiplier;
}

Real
Duple_Field_2D::evaluate_v_at (const Point_2D& point,
                               const Evaluate_Op evaluate_op) const
{
   return evaluate_v_at (point.x, point.y, evaluate_op);
}

*/

void
Streamliner::step (Point_2D& rk,
                   const Real x,
                   const Real y,
                   const Real h) const
{

   const Point_2D& p = transform.reverse (x, y);
   const Integer& ui = u_index;
   const Integer& vi = v_index;

   const Real u = vector_data_2d.evaluate_nocheck (ui, p.x, p.y, VALUE);
   const Real v = vector_data_2d.evaluate_nocheck (vi, p.x, p.y, VALUE);
   const Real theta = transform.get_theta (u, v * aspect, p.x, p.y);

   rk.x = h * cos (theta);
   rk.y = h * sin (theta);

}

void
Streamliner::integrate (Point_2D& rk,
                        const Point_2D& point_2d,
                        const Real h) const
{

   switch (integration_scheme)
   {

      case RUNGA_KUTTA:
      {

         Point_2D k1, k2, k3, k4;

         step (k1, point_2d.x, point_2d.y, h);
         step (k2, point_2d.x + k1.x*0.5, point_2d.y + k1.y*0.5, h);
         step (k3, point_2d.x + k2.x*0.5, point_2d.y + k2.y*0.5, h);
         step (k4, point_2d.x + k3.x, point_2d.y+k3.y, h);

         rk.x = (k1.x + k4.x) / 6 + (k2.x + k3.x) / 3 + point_2d.x;
         rk.y = (k1.y + k4.y) / 6 + (k2.y + k3.y) / 3 + point_2d.y;

         //Point_2D k5;
         //step (k5, rk, 1);
         //Real error_x = (k4.x - k5.x) / 6;
         //Real error_y = (k4.y - k5.y) / 6;
         //Real error_estimate = sqrt (error_x*error_x + error_y*error_y);

         break;

      }

      default:
      case EULER:
      {
         step (rk, point_2d.x, point_2d.y, h);
         rk.x += point_2d.x;
         rk.y += point_2d.y;
         break;
      }

   }

}

Streamliner::Streamliner (const Box_2D& box_2d,
                          const Transform_2D& transform,
                          const Vector_Data_2D& vector_data_2d,
                          const Integer u_index,
                          const Integer v_index,
                          const Integration_Scheme integration_scheme,
                          const Real aspect)
   : Box_2D (box_2d),
     transform (transform),
     vector_data_2d (vector_data_2d),
     u_index (u_index),
     v_index (v_index),
     integration_scheme (integration_scheme),
     aspect (aspect)
{
}

Intensity_Data::Intensity_Data (const Box_2D& box_2d)
   : anchor (box_2d.index_2d),
     size_2d (box_2d.size_2d)
{

   data = new Real[size_2d.i * size_2d.j];

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         data[i * size_2d.j + j] = 0;
      }
   }

}

Intensity_Data::~Intensity_Data ()
{
   delete[] data;
}

void
Intensity_Data::set_noise ()
{

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         set_raw (denise::random (), i, j);
      }
   }

}

void
Intensity_Data::set (const Real intensity,
                     const Integer i,
                     const Integer j)
{

   const Integer ii = i - anchor.i;
   const Integer jj = j - anchor.j;

   if (ii < 0 || ii >= size_2d.i || jj < 0 || jj >= size_2d.j)
   {
      throw Out_Of_Bounds_Exception ("Out Of Bounds");
   }

   set_raw (intensity, ii,  jj);

}

void
Intensity_Data::set_raw (const Real intensity,
                         const Integer i,
                         const Integer j)
{
   data[i * size_2d.j + j] = intensity;
}

Real&
Intensity_Data::get (const Integer i,
                     const Integer j) const
{

   const Integer ii = i - anchor.i;
   const Integer jj = j - anchor.j;

   if (ii < 0 || ii >= size_2d.i || jj < 0 || jj >= size_2d.j)
   {
      throw Out_Of_Bounds_Exception ();
   }

   return get_raw (ii, jj);
}   

Real&
Intensity_Data::get_raw (const Integer i,
                         const Integer j) const
{
   return data[i * size_2d.j + j];
}   

void
Intensity_Data::enhance (const Real min_intensity,
                         const Real max_intensity)
{

   Real span = max_intensity - min_intensity;

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         Real& intensity = get (i, j);
         intensity = ((intensity - min_intensity) / span);
      }
   }

}

Hits::Hits (const Box_2D& box_2d)
   : anchor (box_2d.index_2d),
     size_2d (box_2d.size_2d)
{

   data = new uint8_t[size_2d.i * size_2d.j];

   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         data[i * size_2d.j + j] = 0;
      }
   }

}

Hits::~Hits ()
{
   delete[] data;
}

uint8_t&
Hits::get (const Index_2D& index_2d) const
{

   const Integer ii = index_2d.i - anchor.i;
   const Integer jj = index_2d.j - anchor.j;

   if (ii < 0 || ii >= size_2d.i || jj < 0 || jj >= size_2d.j)
   {
      throw Out_Of_Bounds_Exception ();
   }

   return data[ii * size_2d.j + jj];

}   

void
Lic_Raster::grow (const Point_2D& point_2d,
                  Point_2D* point_array,
                  Real* noise_array,
                  const Intensity_Data& white_noise) const
{

   Index_2D index_2d;
   Integer L_plus_M_2 = 2 * L_plus_M;
   Integer array_size = 2 * L_plus_M + 1;

   for (Integer i = 0; i < array_size; i++)
   {
      point_array[i].x = GSL_NAN;
      point_array[i].y = GSL_NAN;
      noise_array[i] = GSL_NAN;
   }

   index_2d.i = Integer (rint (point_2d.x));
   index_2d.j = Integer (rint (point_2d.y));

   point_array[L_plus_M] = point_2d;
   noise_array[L_plus_M] = white_noise.get (index_2d.i, index_2d.j);    

   for (Integer i = L_plus_M; i < L_plus_M_2; i++)
   {

      try
      {
         Point_2D& next_point = point_array[i + 1];
         integrate (next_point, point_array[i], h);
         index_2d.i = Integer (rint (next_point.x));
         index_2d.j = Integer (rint (next_point.y));
         noise_array[i + 1] = white_noise.get (index_2d.i, index_2d.j);
      }
      catch (const Out_Of_Bounds_Exception& oobe) { break; }

   }

   for (Integer i = L_plus_M; i > 0; i--)
   {

      try
      {
         Point_2D& next_point = point_array[i - 1];
         integrate (next_point, point_array[i], -h);
         index_2d.i = Integer (rint (next_point.x));
         index_2d.j = Integer (rint (next_point.y));
         noise_array[i - 1] = white_noise.get (index_2d.i, index_2d.j);
      }
      catch (const Out_Of_Bounds_Exception& oobe) { break; }

   }

}

void
Lic_Raster::get_sequence (Index_2D* sequence_array,
                          const Box_2D& box_2d) const
{

   Integer sequence_size = box_2d.size_2d.i * box_2d.size_2d.j;

   for (Integer s = 0; s < sequence_size; s++)
   {
      Index_2D& i2d = sequence_array[s];
      i2d.i = box_2d.index_2d.i + s / box_2d.size_2d.j;
      i2d.j = box_2d.index_2d.j + s % box_2d.size_2d.j;
   }

   Integer n = sequence_size;

   // randomize sequence
   for (Integer i = 0; i < sequence_size; i++)
   {
      n--;
      Index_2D& index_2d = sequence_array[irandom (n)];
      const Index_2D temp = index_2d;
      index_2d = sequence_array[n];
      sequence_array[n] = temp;
   }

}

Intensity_Data*
Lic_Raster::get_slow_lic_intensity_data_ptr (Real& min_intensity,
                                             Real& max_intensity,
                                             const Box_2D& box_2d) const
{

/*
   Index_2D i2d;
   Size_2D size_2d = box_2d.size_2d;
   Index_2D anchor = box_2d.index_2d;

Intensity_Data white_noise (box_2d);
   white_noise.set_noise ();

   Intensity_Data* intensity_data_ptr = new Intensity_Data (box_2d);
   Intensity_Data& intensity_data = *intensity_data_ptr;

   min_intensity = GSL_POSINF;
   max_intensity = GSL_NEGINF;

   Index_2D index_2d;

   for (i2d.i = anchor.i; i2d.i < anchor.i + size_2d.i; i2d.i++)
   {
      for (i2d.j = anchor.j; i2d.j < anchor.j + size_2d.j; i2d.j++)
      {

         // for each pixel

         Integer n = 0;
         Real intensity = 0;

         Point_2D p, np;
         Point_2D point (i2d.i, i2d.j);

         const Point_2D& op = transform.reverse (point.x, point.y);
         const Real u = vector_data_2d.evaluate (u_index, op.x, op.y, VALUE);
         const Real v = vector_data_2d.evaluate (v_index, op.x, op.y, VALUE);
         const Real magnitude = sqrt (u*u + v*v);

//         Integer N = Integer (rint (magnitude / magnitude_0 * NN)) + 1;
         const Integer N = 5 * Integer (rint (magnitude)) + 5;

         p.x = point.x;
         p.y = point.y;

         // integrate a streamline, backward and forward
         for (Integer k = 0; k < N; k++)
         {

            try
            {
               integrate (np, p, h, EULER);
               index_2d.i = Integer (rint (np.x));
               index_2d.j = Integer (rint (np.y));
               intensity += white_noise.get (index_2d.i, index_2d.j);
               n++;
               p = np;
            }
            catch (const Out_Of_Bounds_Exception& oobe) { break; }

         }

         p.x = point.x;
         p.y = point.y;

         for (Integer k = N-1; k >= 0; k--)
         {

            try
            {
               integrate (np, p, -h, EULER);
               index_2d.i = Integer (rint (np.x));
               index_2d.j = Integer (rint (np.y));
               intensity += white_noise.get (index_2d.i, index_2d.j);
               n++;
               p = np;
            }
            catch (const Out_Of_Bounds_Exception& oobe) { break; }

         }

         // calculate the LIC
         intensity /= Real (n);

         if (intensity < min_intensity) { min_intensity = intensity; }
         if (intensity > max_intensity) { max_intensity = intensity; }

         intensity_data.set (intensity, i2d.i, i2d.j);

      }
   }

   return intensity_data_ptr;
*/

}

Intensity_Data*
Lic_Raster::get_lic_intensity_data_ptr (Real& min_intensity,
                                        Real& max_intensity,
                                        const Box_2D& box_2d) const
{

   Size_2D size_2d = box_2d.size_2d;
   Index_2D anchor = box_2d.index_2d;

   Integer min_hits = 1;
   Integer array_size = 2 * L_plus_M + 1;
   Integer standard_streamline_length = 2*L + 1;
   Integer sequence_size = size_2d.i * size_2d.j;

   Real* noise_array = new Real[array_size];
   Point_2D* point_array = new Point_2D[array_size];
   Index_2D* sequence_array = new Index_2D[sequence_size];

   Intensity_Data white_noise (box_2d);
   white_noise.set_noise ();

   Hits hits (box_2d);

   Intensity_Data* intensity_data_ptr = new Intensity_Data (box_2d);
   Intensity_Data& intensity_data = *intensity_data_ptr;

   get_sequence (sequence_array, box_2d);

   min_intensity = GSL_POSINF;
   max_intensity = GSL_NEGINF;

   for (Integer s = 0; s < sequence_size; s++)
   {

      Index_2D& ij = sequence_array[s];

      uint8_t& h = hits.get (ij);
      if (h >= min_hits) { continue; }

      const Point_2D point_2d (Real (ij.i), Real (ij.j));
      if (!box_2d.contains (ij.i, ij.j)) { continue; }

      grow (point_2d, point_array, noise_array, white_noise);

      // along the streamline
      //#pragma omp parallel for
      for (Integer k = -M; k <= M; k++)
      {

         try
         {

            Real intensity;
            const Point_2D& p = point_array[L_plus_M + k];
            if (p.is_nap ()) { continue; }

            Index_2D i2d (Integer (rint (p.x)), Integer (rint (p.y)));
            uint8_t& h = hits.get (i2d);
            if (h >= min_hits) { continue; }

            {
               Real v = 0;
               Integer n = standard_streamline_length;
               for (Integer l = -L; l <= L; l++)
               {
                  const Real& noise = noise_array[L_plus_M + k + l];
                  if (!gsl_isnan (noise)) { v += noise; } else { n--; }
               }
               intensity = v / n;
            }

            if (intensity < min_intensity) { min_intensity = intensity; }
            if (intensity > max_intensity) { max_intensity = intensity; }

            intensity_data.set (intensity, i2d.i, i2d.j);
            h++;

         }
         catch (const Out_Of_Bounds_Exception& oobe) { }

      }

   }

   delete[] noise_array;
   delete[] point_array;
   delete[] sequence_array;

   return intensity_data_ptr;

}

Lic_Raster::Lic_Raster (const Transform_2D& transform,
                        const Vector_Data_2D& vector_data_2d,
                        const Box_2D& box_2d,
                        const Integer u_index,
                        const Integer v_index,
                        const bool enhance,
                        const Real aspect,
                        const Integration_Scheme integration_scheme,
                        const Integer L,
                        const Integer M,
                        const Real h)
   : Raster (box_2d),
     Streamliner (box_2d, transform, vector_data_2d, u_index,
                  v_index, integration_scheme, aspect),
     L (L),
     M (M),
     h (h),
     L_plus_M (L + M)
{

   Real hue = 0.0;
   Real min_intensity, max_intensity;
   
   const Size_2D& size_2d = box_2d.size_2d;
   const Index_2D& anchor = box_2d.index_2d;
 
   Intensity_Data* intensity_data_ptr =
      get_lic_intensity_data_ptr (min_intensity, max_intensity, box_2d);
   Intensity_Data& intensity_data = *intensity_data_ptr;

   if (enhance)
   {
      intensity_data.enhance (min_intensity, max_intensity);
   }

   for (Integer i = anchor.i; i < anchor.i + size_2d.i; i++)
   {
      for (Integer j = anchor.j; j < anchor.j + size_2d.j; j++)
      {

         const Point_2D& p = transform.reverse (Real (i), Real (j));
         if (vector_data_2d.out_of_bounds (p.x, p.y)) { continue; }

         const Real u = vector_data_2d.evaluate (u_index, p.x, p.y, VALUE);
         const Real v = vector_data_2d.evaluate (v_index, p.x, p.y, VALUE);
         const Real magnitude = sqrt (u*u + v*v);

         const Real saturation = rint (magnitude / 10) / 8;

         const Real intensity = intensity_data.get (i, j);
         const Color& color = Color::hsb (hue, 0.0, intensity);
         //const Color& color = Color::hsb (intensity, 0.2, 1, alpha);
         //const Color& color = Color::hsb (saturation,
         //   saturation, intensity, 0.6);
         set_pixel (i - anchor.i, j - anchor.j, color);

      }
   }

   delete intensity_data_ptr;

}

Streamlines::Uv_Data::Uv_Data (const Vector_Data_2D& data,
                               const Integer u_index,
                               const Integer v_index,
                               const Real u_multiplier,
                               const Real v_multiplier)
   : data (data),
     u_index (u_index),
     v_index (v_index),
     u_multiplier (u_multiplier),
     v_multiplier (v_multiplier)
{
}

Streamlines::Singular_Point::Singular_Point (const Point_2D& point,
                                             const Jacobian_2D& jacobian)
   : Point_2D (point),
     jacobian (jacobian)
{
}

Streamlines::Singular_Point::Singular_Point (const Singular_Point& singular_point)
   : Point_2D (singular_point),
     jacobian (singular_point.jacobian)
{
}

bool
Streamlines::Singular_Point::is_saddle () const
{

   const Real u_x = jacobian.get_u_x ();
   const Real u_y = jacobian.get_u_y ();
   const Real v_x = jacobian.get_v_x ();
   const Real v_y = jacobian.get_v_y ();

   const Real a = (u_x * v_y - v_x * u_y);
   const Real b = u_x * u_x + v_x * v_x + u_y * u_y + v_y * v_y;
   return (a / b < 0);

}

pair<Real, Real>
Streamlines::Singular_Point::get_eigenvalue_pair () const
{

   const Real divergence = jacobian.get_divergence ();
   const Real jacobian_determinant = jacobian.get_determinant ();
   const Real determinant = divergence*divergence - 4*jacobian_determinant;

   const Real lambda_0 = (divergence + sqrt (determinant)) / 2;
   const Real lambda_1 = (divergence - sqrt (determinant)) / 2;
   return make_pair (lambda_0, lambda_1);

}

pair<Real, Real>
Streamlines::Singular_Point::get_eigentheta_pair () const
{

   const Real divergence = jacobian.get_divergence ();
   const Real jacobian_determinant = jacobian.get_determinant ();
   const Real determinant = divergence*divergence - 4*jacobian_determinant;

   const Real lambda_0 = (divergence + sqrt (determinant)) / 2;
   const Real lambda_1 = (divergence - sqrt (determinant)) / 2;
   const Real v_x = jacobian.get_v_x ();
   const Real v_y = jacobian.get_v_y ();

   const Real theta_0 = atan2 (v_x, (lambda_0 - v_y));
   const Real theta_1 = atan2 (v_x, (lambda_1 - v_y));

   return make_pair (theta_0, theta_1);

}

Real
Streamlines::Singular_Point::get_smaller_eigentheta () const
{

   const Real divergence = jacobian.get_divergence ();
   const Real jacobian_determinant = jacobian.get_determinant ();
   const Real determinant = divergence*divergence - 4*jacobian_determinant;

   const Real lambda_0 = (divergence + sqrt (determinant)) / 2;
   const Real lambda_1 = (divergence - sqrt (determinant)) / 2;
   const Real v_x = jacobian.get_v_x ();
   const Real v_y = jacobian.get_v_y ();

   const bool smaller_0 = fabs (lambda_0) < fabs (lambda_1);
   const Real lambda = (smaller_0 ? lambda_0 : lambda_1);
   return atan2 (v_x, (lambda - v_y));

}

int
Streamlines::Singular_Points::gsl_get_value (const gsl_vector* x,
                                             void* p,
                                             gsl_vector* f)
{

   const Vector_Data_2D* data_ptr = (const Vector_Data_2D*)p;
   const Vector_Data_2D& data = *data_ptr;

   const Real xx = gsl_vector_get (x, 0);
   const Real yy = gsl_vector_get (x, 1);

   gsl_vector_set (f, 0, data.evaluate_nocheck (0, xx, yy));
   gsl_vector_set (f, 1, data.evaluate_nocheck (1, xx, yy));

   return GSL_SUCCESS;

}

void
Streamlines::Singular_Points::potentially_add (const Uv_Data& uv_data,
                                               const Real a,
                                               const Real c,
                                               const Real delta_x,
                                               const Real delta_y,
                                               const Real s,
                                               const Real t)
{

   if (gsl_isnan (s) || gsl_isnan (t)) { return; }
   if (s < 0 || s >= 1 || t < 0 || t >= 1) { return; }

   const Vector_Data_2D& data = uv_data.data;
   const Integer u_index = uv_data.u_index;
   const Integer v_index = uv_data.v_index;

   const Point_2D p (s * delta_x + a, t * delta_y + c);
   const Jacobian_2D& j = data.get_jacobian_nocheck (u_index, v_index, p);
   const Singular_Point singular_point (p, j);
   push_back (singular_point);

}

Streamlines::Singular_Points::Singular_Points (const Uv_Data& uv_data,
                                               const Domain_2D& domain_2d,
                                               const Real residual,
                                               const Real epsilon)
   : uv_data (uv_data),
     residual (residual),
     epsilon (epsilon)
{

   const Vector_Data_2D& data = uv_data.data;
   const Integer u_index = uv_data.u_index;
   const Integer v_index = uv_data.v_index;

   const Tuple& tuple_x = data.get_coordinate_tuple (0);
   const Tuple& tuple_y = data.get_coordinate_tuple (1);
   const Integer ni = tuple_x.size ();
   const Integer nj = tuple_y.size ();

   for (Integer i = 0; i < ni - 1; i++)
   {

      const Real& a = tuple_x[i];
      const Real& b = tuple_x[i + 1];
      if (domain_2d.domain_x.is_out_of_bounds (a) ||
          domain_2d.domain_x.is_out_of_bounds (b)) { continue; }

      const Real delta_x = b - a;

      for (Integer j = 0; j < nj - 1; j++)
      {

         const Real& c = tuple_y[j];
         const Real& d = tuple_y[j + 1];
         if (domain_2d.domain_y.is_out_of_bounds (c) ||
             domain_2d.domain_y.is_out_of_bounds (d)) { continue; }

         const Real delta_y = d - c;

         const Real u_ac = data.get_datum (u_index, i, j);
         const Real u_ad = data.get_datum (u_index, i, j + 1);
         const Real u_bc = data.get_datum (u_index, i + 1, j);
         const Real u_bd = data.get_datum (u_index, i + 1, j + 1);

         const Real v_ac = data.get_datum (v_index, i, j);
         const Real v_ad = data.get_datum (v_index, i, j + 1);
         const Real v_bc = data.get_datum (v_index, i + 1, j);
         const Real v_bd = data.get_datum (v_index, i + 1, j + 1);

         const bool u_all_pos = (u_ac > 0 && u_ad > 0 && u_bc > 0 && u_bd > 0);
         const bool u_all_neg = (u_ac < 0 && u_ad < 0 && u_bc < 0 && u_bd < 0);
         const bool v_all_pos = (v_ac > 0 && v_ad > 0 && v_bc > 0 && v_bd > 0);
         const bool v_all_neg = (v_ac < 0 && v_ad < 0 && v_bc < 0 && v_bd < 0);
         const bool u_same_sign = (u_all_pos || u_all_neg);
         const bool v_same_sign = (v_all_pos || v_all_neg);
         if (u_same_sign && v_same_sign) { continue; }

         const Real a_u = (u_bc - u_ac);
         const Real b_u = (u_ad - u_ac);
         const Real c_u = (u_ac + u_bd - u_ad - u_bc);
         const Real d_u = u_ac;

         const Real a_v = (v_bc - v_ac);
         const Real b_v = (v_ad - v_ac);
         const Real c_v = (v_ac + v_bd - v_ad - v_bc);
         const Real d_v = v_ac;

         const Real A = (b_u * c_v - b_v * c_u);
         const Real B = (a_v * b_u - a_u * b_v + c_v * d_u - c_u * d_v);
         const Real C = (a_v * d_u - a_u * d_v);

         const Real determinant = B*B - 4*A*C;

         if (fabs (A) < 1e-9)
         {
            const Real t = -C / B;
            const Real s = -(d_v + b_v * t) / (a_v + c_v * t);
            potentially_add (uv_data, a, c, delta_x, delta_y, s, t);
         }
         else
         if (fabs (determinant) < 1e-9)
         {
            const Real t = -B / (A+A);
            const Real s = -(d_v + b_v * t) / (a_v + c_v * t);
            potentially_add (uv_data, a, c, delta_x, delta_y, s, t);
         }
         else
         if (determinant > 0)
         {

            const Real determinant_squared_2 = sqrt (determinant) / 2;

            const Real t_0 = (sqrt (determinant) - B) / (A+A);
            const Real s_0 = -(d_v + b_v * t_0) / (a_v + c_v * t_0);
            const Real t_1 = (-sqrt (determinant) - B) / (A+A);
            const Real s_1 = -(d_v + b_v * t_1) / (a_v + c_v * t_1);

            potentially_add (uv_data, a, c, delta_x, delta_y, s_0, t_0);
            potentially_add (uv_data, a, c, delta_x, delta_y, s_1, t_1);

         }

      }

   }

}

list<Streamlines::Singular_Point>
Streamlines::Singular_Points::get_saddle_point_list () const
{
   list<Singular_Point> saddle_point_list;
   for (Singular_Points::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Singular_Point& sp = *(iterator);
      if (sp.is_saddle ()) { saddle_point_list.push_back (sp); }
   }
   return saddle_point_list;
}

void
Streamlines::Streamline::set_nan (Real& dx,
                                  Real& dy) const
{
   dx = GSL_NAN;
   dy = GSL_NAN;
}

void
Streamlines::Streamline::euler (Real& dx,
                                Real& dy,
                                const Uv_Data& uv_data,
                                const Real x,
                                const Real y,
                                const Real h) const
{

   if (gsl_isnan (x) || gsl_isnan (y)) { set_nan (dx, dy); return; }

   const Vector_Data_2D& data = uv_data.data;
   const Integer& ui = uv_data.u_index;
   const Integer& vi = uv_data.v_index;

   const Real u = data.evaluate_nocheck (ui, x, y);
   const Real v = data.evaluate_nocheck (vi, x, y);
   const Real theta = atan2 (v, u);

   dx = h * cos (theta);
   dy = h * sin (theta);

}

void
Streamlines::Streamline::runga_kutta (Real& dx,
                                      Real& dy,
                                      const Uv_Data& uv_data,
                                      const Real x,
                                      const Real y,
                                      const Real h) const
{

   if (gsl_isnan (x) || gsl_isnan (y)) { set_nan (dx, dy); return; }

   Point_2D k1, k2, k3, k4;

   euler (k1.x, k1.y, uv_data, x, y, h);
   if (gsl_isnan (k1.x) || gsl_isnan (k1.y)) { set_nan (dx, dy); return; }

   euler (k2.x, k2.y, uv_data, x + k1.x*0.5, y + k1.y*0.5, h);
   if (gsl_isnan (k2.x) || gsl_isnan (k2.y)) { set_nan (dx, dy); return; }

   euler (k3.x, k3.y, uv_data, x + k2.x*0.5, y + k2.y*0.5, h);
   if (gsl_isnan (k3.x) || gsl_isnan (k3.y)) { set_nan (dx, dy); return; }

   euler (k4.x, k4.y, uv_data, x + k3.x, y + k3.y, h);
   if (gsl_isnan (k4.x) || gsl_isnan (k4.y)) { set_nan (dx, dy); return; }

   dx = (k1.x + k4.x) / 6 + (k2.x + k3.x) / 3;
   dy = (k1.y + k4.y) / 6 + (k2.y + k3.y) / 3;

}

Point_2D
Streamlines::Streamline::integrate (const Uv_Data& uv_data,
                                    const Real x,
                                    const Real y,
                                    const Real h,
                                    const Integration_Scheme scheme) const
{

   Real dx = 0, dy = 0;

   switch (scheme)
   {
      default:
      case EULER: euler (dx, dy, uv_data, x, y, h); break;
      case RUNGA_KUTTA: runga_kutta (dx, dy, uv_data, x, y, h); break;
   }

   return Point_2D (x + dx, y + dy);

}

Streamlines::Streamline::Streamline (const bool forward,
                                     const bool backward,
                                     const bool forward_stopped,
                                     const bool backward_stopped)
   : forward (forward),
     backward (backward),
     forward_stopped (forward_stopped),
     backward_stopped (backward_stopped)
{
}

Streamlines::Streamline::Streamline (const Point_2D& point,
                                     const bool forward,
                                     const bool backward)
   : forward (forward),
     backward (backward),
     forward_stopped (false),
     backward_stopped (false)
{
   add (point);
}

Streamlines::Streamline::Streamline (const Point_2D& point_0,
                                     const Point_2D& point_1,
                                     const bool forward,
                                     const bool backward)
   : forward (forward),
     backward (backward),
     forward_stopped (false),
     backward_stopped (false)
{

   if (!forward) { forward_stopped = true; }
   if (!backward) { backward_stopped = true; }

   add (point_0); add (point_1);
}

void
Streamlines::Streamline::stop ()
{
   forward_stop ();
   backward_stop ();
}

void
Streamlines::Streamline::forward_stop ()
{
   forward_stopped = true;
}

void
Streamlines::Streamline::backward_stop ()
{
   backward_stopped = true;
}

bool
Streamlines::Streamline::is_stopped () const
{
   return (is_forward_stopped () && is_backward_stopped ());
}

bool
Streamlines::Streamline::is_forward_stopped () const
{
   return forward_stopped;
}

bool
Streamlines::Streamline::is_backward_stopped () const
{
   return backward_stopped;
}

const Point_2D&
Streamlines::Streamline::tip () const
{
   return (forward ? back () : front ());
}

void
Streamlines::Streamline::add_to (Sand_Box& sand_box) const
{
   for (Streamline::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Point_2D& point = *(iterator);
      sand_box.add (point);
   }
}

void
Streamlines::Streamline::grow_step (const Uv_Data& uv_data,
                                    const Domain_2D& domain_2d,
                                    Sand_Box& sand_box,
                                    const Real h,
                                    const Integration_Scheme integration_scheme)
{

   // forward
   const bool f = (h > 0);

   const Point_2D& p = (f ? back () : front ());
   const Point_2D& np = integrate (uv_data, p.x, p.y, h, integration_scheme);
   if (np.is_nap ())
   {
      if (f) { forward_stop (); } else { backward_stop (); }
      return;
   }

   const Real theta = atan2 (np.y - p.y, np.x - p.x);
   const bool too_close = sand_box.too_close (np, theta);
   if (too_close)
   {
      if (f) { forward_stop (); } else { backward_stop (); }
      return;
   } 

   const bool out_of_bounds = domain_2d.is_out_of_bounds (np.x, np.y);
   if (out_of_bounds)
   {
      if (f) { forward_stop (); } else { backward_stop (); }
      return;
   }

   if (f) { add (np); } else { prepend (np); }
   sand_box.add (p);

}

void
Streamlines::Streamline::render (const RefPtr<Context>& cr,
                                 const Transform_2D& transform) const
{

   Simple_Polyline::cairo (cr, transform);
   cr->stroke ();

   const Integer d = 120;

   for (Streamline::const_iterator i = begin (); i != end (); i++)
   {
      if (std::distance (begin (), i) % d != d/10) { continue; }
      Streamline::const_iterator j = i; j--;
      Streamline::const_iterator k = j; k--;

      const Point_2D& p_a = transform.transform (*(k));
      const Point_2D& p_b = transform.transform (*(i));
      const Point_2D& p = transform.transform (*(j));

      const Real theta = atan2 ((p_b.y - p_a.y), (p_b.x - p_a.x));
      const Triangle triangle (5, theta);

      triangle.cairo (cr, p);
      cr->fill ();

   }

}

void
Streamlines::Streamline::grow (const Uv_Data& uv_data,
                               const Domain_2D& domain_2d,
                               Sand_Box& sand_box,
                               const Real step_size,
                               const Integration_Scheme scheme)
{

   const bool go_forward = (forward && !forward_stopped);
   const bool go_backward = (backward && !backward_stopped);

   if (go_forward)
   {
      grow_step (uv_data, domain_2d, sand_box, step_size, scheme);
   }

   if (go_backward)
   {
      grow_step (uv_data, domain_2d, sand_box, -step_size, scheme);
   }

}

void
Streamlines::Separatrix::grow_step (const Uv_Data& uv_data,
                                    const Domain_2D& domain_2d,
                                    Sand_Box& sand_box,
                                    const Real h,
                                    const Integration_Scheme integration_scheme)
{

   const bool hh = (forward ? fabs (h) : -fabs (h));

   const Point_2D& p = tip ();
   const Point_2D& np = integrate (uv_data, p.x, p.y, hh, integration_scheme);
   if (np.is_nap ()) { stop (); return; }

   const Real theta = atan2 (np.y - p.y, np.x - p.x);
   const bool too_close (sand_box.too_close (np));
   if (too_close) { stop (); return; }

   const bool out_of_bounds = domain_2d.is_out_of_bounds (np.x, np.y);
   if (out_of_bounds) { stop (); return; }

   if (forward) { add (np); } else { prepend (np); }
   sand_box.add (p);

}

Streamlines::Separatrix::Separatrix (const bool forward,
                                     const bool backward,
                                     const bool forward_stopped,
                                     const bool backward_stopped)
   : Streamline (forward, backward, forward_stopped, backward_stopped)
{
}

Streamlines::Separatrix::Separatrix (const Point_2D& point,
                                     const bool outgoing)
   : Streamline (outgoing, !outgoing, !outgoing, outgoing)
{
   add (point);
}

const Point_2D&
Streamlines::Separatrix::get_origin () const
{
   return (forward ? front () : back ());
}

const Point_2D&
Streamlines::Separatrix::get_tip () const
{
   return (forward ? back () : front ());
}

void
Streamlines::Separatrix::grow (const Uv_Data& uv_data,
                               const Domain_2D& domain_2d,
                               Sand_Box& sand_box,
                               const Real step_size,
                               const Integration_Scheme integration_scheme)
{

   const Real h = (forward ? fabs (step_size) : -fabs (step_size));

   const Point_2D& p = tip ();
   const Point_2D& np = integrate (uv_data, p.x, p.y, h, integration_scheme);
   if (np.is_nap ()) { stop (); return; }

   const Real theta = atan2 (np.y - p.y, np.x - p.x);
   const bool too_close (sand_box.too_close (np, theta));
   if (too_close) { stop (); return; }

   const bool out_of_bounds = domain_2d.is_out_of_bounds (np.x, np.y);
   if (out_of_bounds) { stop (); return; }

   if (forward) { add (np); } else { prepend (np); }
   sand_box.add (p);

}

void
Streamlines::Separatrix::render (const RefPtr<Context>& cr,
                                 const Transform_2D& transform) const
{
   Streamline::render (cr, transform);
}

void
Streamlines::Web_Stick::grow_step (const Uv_Data& uv_data,
                                   const Domain_2D& domain_2d,
                                   const Real step_size,
                                   const Integration_Scheme integration_scheme)
{

   const Real h = (forward ? fabs (step_size) : -fabs (step_size));

   const Point_2D& p = tip ();
   const Point_2D& np = integrate (uv_data, p.x, p.y, h, integration_scheme);
   if (np.is_nap ()) { stop (); return; }

   const bool out_of_bounds = domain_2d.is_out_of_bounds (np.x, np.y);
   if (out_of_bounds) { stop (); return; }

   if (forward) { add (np); } else { prepend (np); }

}

Streamlines::Web_Stick::Web_Stick (const Singular_Point& saddle_point,
                                   const Point_2D& point,
                                   const bool outgoing)
   : Separatrix (outgoing, !outgoing, !outgoing, outgoing)
{
   if (outgoing) { add (saddle_point); add (point); }
   else { prepend (saddle_point); prepend (point); }
}

void
Streamlines::Web_Stick::grow (const Uv_Data& uv_data,
                              const Domain_2D& domain_2d,
                              const Integer steps,
                              const Real h,
                              const Integration_Scheme integration_scheme)
{
   for (Integer count = 0; count < steps; count++)
   {
      grow_step (uv_data, domain_2d, h, integration_scheme);
   }
}

Streamlines::Saddle::Saddle (const Singular_Point& saddle_point,
                             const Uv_Data& uv_data,
                             const Domain_2D& domain_2d,
                             const Integer size,
                             const Real step_size,
                             const Integration_Scheme integration_scheme)
   : Singular_Point (saddle_point)
{

   const Singular_Point& sp = saddle_point;
   const Vector_Data_2D& data = uv_data.data;
   const Integer& ui = uv_data.u_index;
   const Integer& vi = uv_data.v_index;

   pair<Real, Real> thetas = sp.get_eigentheta_pair ();
   const Real theta_0 = thetas.first;
   const Real theta_1 = thetas.second;

   const Real dx_0 = step_size * cos (theta_0);
   const Real dy_0 = step_size * sin (theta_0);
   const Real dx_1 = step_size * cos (theta_1);
   const Real dy_1 = step_size * sin (theta_1);

   const Point_2D p_0 (sp.x + dx_0, sp.y + dy_0);
   const Point_2D p_1 (sp.x + dx_1, sp.y + dy_1);
   const Point_2D p_2 (sp.x - dx_0, sp.y - dy_0);
   const Point_2D p_3 (sp.x - dx_1, sp.y - dy_1);

   const Real u_0 = data.evaluate_nocheck (ui, p_0.x, p_0.y);
   const Real v_0 = data.evaluate_nocheck (vi, p_0.x, p_0.y);
   const Real outgoing_0 = (dx_0*u_0 + dy_0*v_0 > 0);

   push_back (Web_Stick (sp, p_0, outgoing_0));
   push_back (Web_Stick (sp, p_1, !outgoing_0));
   push_back (Web_Stick (sp, p_2, outgoing_0));
   push_back (Web_Stick (sp, p_3, !outgoing_0));

   for (Saddle::iterator iterator = begin (); iterator != end (); iterator++)
   {
      Web_Stick& web_stick = *(iterator);
      web_stick.grow (uv_data, domain_2d, size, step_size, integration_scheme);
   }

}

void
Streamlines::Saddle::add_to (Sand_Box& sand_box) const
{
   for (Saddle::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Web_Stick& web_stick = *(iterator);
      web_stick.add_to (sand_box);
   }
}

void
Streamlines::Saddle::render (const RefPtr<Context>& cr,
                             const Transform_2D& transform) const
{

   vector<Point_2D> tip_vector;

   for (Saddle::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Web_Stick& web_stick = *(iterator);
      const Point_2D& tip = web_stick.get_tip ();
      //const Edge edge (*this, tip);
      //edge.cairo (cr, transform);
      tip_vector.push_back (tip);
      web_stick.render (cr, transform);
   }

   const Point_2D& p_0 = tip_vector[0];
   const Point_2D& p_1 = tip_vector[1];
   const Point_2D& p_2 = tip_vector[2];
   const Point_2D& p_3 = tip_vector[3];

   Tuple tuple (2, this->x, this->y);
   Tuple tuple_0 (2, p_0.x, p_0.y);
   Tuple tuple_1 (2, p_1.x, p_1.y);
   Tuple tuple_2 (2, p_2.x, p_2.y);
   Tuple tuple_3 (2, p_3.x, p_3.y);

   vector<Tuple> tuple_vector_01;
   tuple_vector_01.push_back (tuple_0);
   tuple_vector_01.push_back (tuple);
   tuple_vector_01.push_back (tuple_1);

   vector<Tuple> tuple_vector_12;
   tuple_vector_12.push_back (tuple_1);
   tuple_vector_12.push_back (tuple);
   tuple_vector_12.push_back (tuple_2);

   vector<Tuple> tuple_vector_23;
   tuple_vector_23.push_back (tuple_2);
   tuple_vector_23.push_back (tuple);
   tuple_vector_23.push_back (tuple_3);

   vector<Tuple> tuple_vector_30;
   tuple_vector_30.push_back (tuple_3);
   tuple_vector_30.push_back (tuple);
   tuple_vector_30.push_back (tuple_0);

   list<Bezier_Curve> bezier_curve_list;
   bezier_curve_list.push_back (Bezier_Curve (tuple_vector_01));
   bezier_curve_list.push_back (Bezier_Curve (tuple_vector_12));
   bezier_curve_list.push_back (Bezier_Curve (tuple_vector_23));
   bezier_curve_list.push_back (Bezier_Curve (tuple_vector_30));

   for (list<Bezier_Curve>::const_iterator i = bezier_curve_list.begin ();
        i != bezier_curve_list.end (); i++)
   {
      const Tuple tuple_t (25, 0, 1);
      const Bezier_Curve& bezier_curve = *(i);
      Simple_Polyline simple_polyline;
      for (Tuple::const_iterator j = tuple_t.begin ();
           j != tuple_t.end (); j++)
      {
         const Real t = *(j);
         const Real x = bezier_curve.get_value_at (t, 0);
         const Real y = bezier_curve.get_value_at (t, 1);
         simple_polyline.push_back (Point_2D (x, y));
      }
      simple_polyline.cairo (cr, transform);
      cr->stroke ();
   }

}

void
Streamlines::Saddle_List::grow_separatrices (const Uv_Data& uv_data,
                                             const Domain_2D& domain_2d,
                                             Sand_Box& sand_box,
                                             const Real step_size,
                                             const Integration_Scheme scheme)
{
   
   const Integer n = separatrix_list.size ();

   for (Integer stopped_separatrices = 0; stopped_separatrices != n; )
   {

      stopped_separatrices = 0;

      for (list<Separatrix>::iterator iterator = separatrix_list.begin ();
           iterator != separatrix_list.end (); iterator++)
      {

         Separatrix& separatrix = *(iterator);
         if (separatrix.is_stopped ()) { stopped_separatrices++; continue; }

         separatrix.grow (uv_data, domain_2d, sand_box, step_size, scheme);
         if (separatrix.size () > 800) { separatrix.stop (); }

      }

   }

}

void
Streamlines::Saddle_List::add_to (Sand_Box& sand_box) const
{
   for (Saddle_List::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Saddle& saddle = *(iterator);
      saddle.add_to (sand_box);
   }

}

Streamlines::Saddle_List::Saddle_List (const list<Singular_Point>& sp_list,
                                       const Uv_Data& uv_data,
                                       const Domain_2D& domain_2d,
                                       Sand_Box& sand_box,
                                       const Integer size,
                                       const Real step_size,
                                       const Integration_Scheme scheme)
{

   for (list<Singular_Point>::const_iterator iterator = sp_list.begin ();
        iterator != sp_list.end (); iterator++)
   {
      const Singular_Point& sp = *(iterator);
      if (domain_2d.is_out_of_bounds (sp)) { continue; }
      push_back (Saddle (sp, uv_data, domain_2d, size, step_size, scheme));
   }

   for (Saddle_List::const_iterator i = begin (); i != end (); i++)
   {
      const Saddle& saddle = *(i);
      for (Saddle::const_iterator j = saddle.begin (); j != saddle.end (); j++)
      {
         const Web_Stick& web_stick = *(j);
         const Point_2D& tip = web_stick.tip ();
         separatrix_list.push_back (Separatrix (tip, web_stick.forward));
      }
   }

   grow_separatrices (uv_data, domain_2d, sand_box, step_size, scheme);
   add_to (sand_box);

}

void
Streamlines::Saddle_List::render (const RefPtr<Context>& cr,
                                  const Transform_2D& transform) const
{

   for (Saddle_List::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Saddle& saddle = *(iterator);
      saddle.render (cr, transform);
   }

   for (list<Separatrix>::const_iterator iterator = separatrix_list.begin ();
        iterator != separatrix_list.end (); iterator++)
   {
      const Separatrix& separatrix = *(iterator);
      separatrix.render (cr, transform);
   }

}

Integer
Streamlines::Sand_Box::get_offset (const Integer i,
                                   const Integer j) const
{
   return (i * size_nd.buffer[1]) + j;
}

Streamlines::Sand_Box::Sand_Box (const Domain_2D& domain_2d,
                                 const Real separation)
   : Grid_nD (2),
     separation (separation)
{

   const Domain_1D& domain_x = domain_2d.domain_x;
   const Domain_1D& domain_y = domain_2d.domain_y;
   const Real span_x = domain_x.get_span ();
   const Real span_y = domain_y.get_span ();
   const Integer ni = Integer (ceil (span_x / separation));
   const Integer nj = Integer (ceil (span_y / separation));

   spacings[0] = domain_x.get_span () / (ni - 1);
   spacings[1] = domain_y.get_span () / (nj - 1);
   periodics[0] = false;
   periodics[1] = false;
   coordinate_tuples[0] = Tuple (ni, domain_x.start, domain_x.end);
   coordinate_tuples[1] = Tuple (nj, domain_y.start, domain_y.end);

   size_nd.buffer[0] = ni;
   size_nd.buffer[1] = nj;

   const Integer n = ni * nj;
   typedef list<Point_2D> Point_List;
   point_lists = new Point_List[n];

}

Streamlines::Sand_Box::~Sand_Box ()
{
   delete[] point_lists;
}

bool
Streamlines::Sand_Box::same_box (const Point_2D& point_0,
                                 const Point_2D& point_1) const
{
   const Integer i_0 = get_node (0, point_0.x);
   const Integer j_0 = get_node (1, point_0.y);
   const Integer i_1 = get_node (0, point_1.x);
   const Integer j_1 = get_node (1, point_1.y);
   return (i_0 == i_1) && (j_0 == j_1);
}

void
Streamlines::Sand_Box::add (const Point_2D& point)
{
   const Integer i = get_node (0, point.x);
   const Integer j = get_node (1, point.y);
   const Integer offset = get_offset (i, j);
   list<Point_2D>& point_list = point_lists[offset];
   point_list.push_back (point);
}

bool
Streamlines::Sand_Box::too_close (const Point_2D& point,
                                  const Real theta) const
{

   const Real c = cos (theta);
   const Real s = sin (theta);
   const Integer i0 = get_node (0, point.x);
   const Integer j0 = get_node (1, point.y);
   const Real threshold = separation * separation;

   for (Integer i = i0 - 1; i <= i0 + 1; i++)
   {
      if (i < 0 || i > size_nd.buffer[0]) { continue; }
      for (Integer j = j0 - 1; j <= j0 + 1; j++)
      {
         if (j < 0 || j > size_nd.buffer[1]) { continue; }

         const Integer offset = get_offset (i, j);
         list<Point_2D>& point_list = point_lists[offset];

         for (list<Point_2D>::const_iterator iterator = point_list.begin ();
              iterator != point_list.end (); iterator++)
         {

            const Point_2D& p = *(iterator);

            const Real dx = p.x - point.x;
            const Real dy = p.y - point.y;
            if (!gsl_isnan (theta) && (dx*c + dy*s < 0)) { continue; }

            if (dx*dx + dy*dy < separation) { return true; }

         }

      }
   }

   return false;

}

void
Streamlines::construct_sides (Sand_Box& sand_box,
                              const Real separation,
                              const Real step_size)
{

   const Vector_Data_2D& data = uv_data.data;
   const Integer& ui = uv_data.u_index;
   const Integer& vi = uv_data.v_index;

   const Domain_1D& domain_x = domain_2d.domain_x;
   const Domain_1D& domain_y = domain_2d.domain_y;

   const Real sep = 5 * separation;

   for (Real y = domain_y.start; y < domain_y.end; y += sep)
   {

      const Point_2D point (domain_x.start, y);
      if (sand_box.too_close (point)) { continue; }
      Streamline streamline (point, true, true);

      while (!streamline.is_stopped ())
      {
         streamline.grow (uv_data, domain_2d, sand_box,
            step_size, integration_scheme);
         if (streamline.size () > 3000) { streamline.stop (); break; }
      }

      side_w_list.push_back (streamline);

      const Real u = data.evaluate_nocheck (ui, point.x, point.y);
      const Real v = data.evaluate_nocheck (vi, point.x, point.y);
      const Real theta = atan2 (v, u);
      const Real extra = std::min (2 * sep, sep / fabs (cos (theta)) - sep);
      y += extra;

   }

   for (Real x = domain_x.start; x < domain_x.end; x += sep)
   {

      const Point_2D point (x, domain_y.end);
      if (sand_box.too_close (point)) { continue; }
      Streamline streamline (point, true, true);

      while (!streamline.is_stopped ())
      {
         streamline.grow (uv_data, domain_2d, sand_box,
            step_size, integration_scheme);
         if (streamline.size () > 3000) { streamline.stop (); break; }
      }

      side_n_list.push_back (streamline);

      const Real u = data.evaluate_nocheck (ui, point.x, point.y);
      const Real v = data.evaluate_nocheck (vi, point.x, point.y);
      const Real theta = atan2 (u, v);
      const Real extra = std::min (2 * sep, sep / fabs (cos (theta)) - sep);
      x += extra;

   }

   for (Real y = domain_y.start; y < domain_y.end; y += sep)
   {

      const Point_2D point (domain_x.end, y);
      if (sand_box.too_close (point)) { continue; }
      Streamline streamline (point, true, true);

      while (!streamline.is_stopped ())
      {
         streamline.grow (uv_data, domain_2d, sand_box,
            step_size, integration_scheme);
         if (streamline.size () > 3000) { streamline.stop (); break; }
      }

      side_e_list.push_back (streamline);

      const Real u = data.evaluate_nocheck (ui, point.x, point.y);
      const Real v = data.evaluate_nocheck (vi, point.x, point.y);
      const Real theta = atan2 (v, u);
      const Real extra = std::min (2 * sep, sep / fabs (cos (theta)) - sep);
      y += extra;

   }

   for (Real x = domain_x.start; x < domain_x.end; x += sep)
   {

      const Point_2D point (x, domain_y.start);
      if (sand_box.too_close (point)) { continue; }
      Streamline streamline (point, true, true);

      while (!streamline.is_stopped ())
      {
         streamline.grow (uv_data, domain_2d, sand_box,
            step_size, integration_scheme);
         if (streamline.size () > 3000) { streamline.stop (); break; }
      }

      side_s_list.push_back (streamline);

      const Real u = data.evaluate_nocheck (ui, point.x, point.y);
      const Real v = data.evaluate_nocheck (vi, point.x, point.y);
      const Real theta = atan2 (u, v);
      const Real extra = std::min (2 * sep, sep / fabs (cos (theta)) - sep);
      x += extra;

   }

}

void
Streamlines::construct_streamlines (Sand_Box& sand_box,
                                    const Real separation,
                                    const Real step_size)
{

   const Vector_Data_2D& data = uv_data.data;
   const Poisson_Disk poisson_disk (domain_2d, separation * 2);

   for (Poisson_Disk::const_iterator iterator = poisson_disk.begin ();
        iterator != poisson_disk.end (); iterator++)
   {
      const Point_2D& point = *(iterator);
      if (sand_box.too_close (point)) { continue; }
      Streamline streamline (point, true, true);

      while (!streamline.is_stopped ())
      {
         streamline.grow (uv_data, domain_2d, sand_box,
            step_size, integration_scheme);
         if (streamline.size () > 3000) { streamline.stop (); break; }
      }

      streamline_list.push_back (streamline);

   }

}

Streamlines::Streamlines (const Uv_Data& uv_data,
                          const Real separation,
                          const Real step_size,
                          const Integer saddle_size,
                          const Integration_Scheme integration_scheme)
   : uv_data (uv_data),
     domain_2d (uv_data.data.get_domain_2d ()),
     integration_scheme (integration_scheme)
{

   Sand_Box sand_box (domain_2d, separation);
   Singular_Points singular_points (uv_data, domain_2d);

   list<Singular_Point> saddle_point_list =
      singular_points.get_saddle_point_list ();

   saddle_list_ptr = new Saddle_List (saddle_point_list, uv_data,
      domain_2d, sand_box, saddle_size, step_size, integration_scheme);

   construct_sides (sand_box, separation, step_size);
   construct_streamlines (sand_box, separation, step_size);

}

Streamlines::Streamlines (const Uv_Data& uv_data,
                          const Domain_2D& domain_2d,
                          const Real separation,
                          const Real step_size,
                          const Integer saddle_size,
                          const Integration_Scheme integration_scheme)
   : uv_data (uv_data),
     domain_2d (domain_2d),
     integration_scheme (integration_scheme)
{

   Sand_Box sand_box (domain_2d, separation);
   Singular_Points singular_points (uv_data, domain_2d);

   list<Singular_Point> saddle_point_list =
      singular_points.get_saddle_point_list ();

   saddle_list_ptr = new Saddle_List (saddle_point_list, uv_data,
      domain_2d, sand_box, saddle_size, step_size, integration_scheme);

   construct_sides (sand_box, separation, step_size);
   construct_streamlines (sand_box, separation, step_size);

}

Streamlines::~Streamlines ()
{
   delete saddle_list_ptr;
}

void  
Streamlines::render (const RefPtr<Context>& cr,
                     const Transform_2D& transform) const
{

   const Integer min_length = 25;

   cr->set_line_width (2);
   saddle_list_ptr->render (cr, transform);


   for (list<Streamline>::const_iterator iterator = vortex_line_list.begin ();
        iterator != vortex_line_list.end (); iterator++)
   {
      const Streamline& streamline = *(iterator);
      streamline.render (cr, transform);
      cr->stroke ();
   }

   for (list<Streamline>::const_iterator iterator = side_w_list.begin ();
        iterator != side_w_list.end (); iterator++)
   {
      const Streamline& streamline = *(iterator);
      if (streamline.size () < min_length) { continue; }
      streamline.render (cr, transform);
      cr->stroke ();
   }

   for (list<Streamline>::const_iterator iterator = side_n_list.begin ();
        iterator != side_n_list.end (); iterator++)
   {
      const Streamline& streamline = *(iterator);
      if (streamline.size () < min_length) { continue; }
      streamline.render (cr, transform);
      cr->stroke ();
   }

   for (list<Streamline>::const_iterator iterator = side_e_list.begin ();
        iterator != side_e_list.end (); iterator++)
   {
      const Streamline& streamline = *(iterator);
      if (streamline.size () < min_length) { continue; }
      streamline.render (cr, transform);
      cr->stroke ();
   }

   for (list<Streamline>::const_iterator iterator = side_s_list.begin ();
        iterator != side_s_list.end (); iterator++)
   {
      const Streamline& streamline = *(iterator);
      if (streamline.size () < min_length) { continue; }
      streamline.render (cr, transform);
      cr->stroke ();
   }

   for (list<Streamline>::const_iterator iterator = streamline_list.begin ();
        iterator != streamline_list.end (); iterator++)
   {
      const Streamline& streamline = *(iterator);
      if (streamline.size () < min_length) { continue; }
      streamline.render (cr, transform);
      cr->stroke ();
   }

}

