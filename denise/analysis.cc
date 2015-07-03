//
// analysis.cc
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

#include <typeinfo>
#include "analysis.h"
#include "util.h"

using namespace denise;

Real
Differentiation::diff (const Real y_0,
                       const Real y_1,
                       const Real h)
{
   return (y_1 - y_0) / h;
}

Real
Differentiation::diff_0 (const Real y_0,
                         const Real y_1,
                         const Real y_2,
                         const Real h)
{
   return (4 * y_1 - 3 * y_0 - y_2) / (h + h);
}

Real
Differentiation::diff_1 (const Real y_0,
                         const Real y_2,
                         const Real h)
{
   return (y_2 - y_0) / (h + h);
}

Real
Differentiation::diff_2 (const Real y_0,
                         const Real y_1,
                         const Real y_2,
                         const Real h)
{
   return (y_0 - 4 * y_1 + 3 * y_2) / (h + h);
}

Real
Differentiation::diff_2nd (const Real y_0,
                           const Real y_1,
                           const Real y_2,
                           const Real h)
{
   return (y_0 - 2 * y_1 + y_2) / (h * h);
}

Jacobian_2D::Jacobian_2D (const Real u_x,
                          const Real u_y,
                          const Real v_x,
                          const Real v_y)
   : u_x (u_x),
     u_y (u_y),
     v_x (v_x),
     v_y (v_y) 
{
}

Jacobian_2D::Jacobian_2D (const Jacobian_2D& jacobian_2d)
   : u_x (jacobian_2d.u_x),
     u_y (jacobian_2d.u_y),
     v_x (jacobian_2d.v_x),
     v_y (jacobian_2d.v_y) 
{
}

Real
Jacobian_2D::get_u_x () const
{
   return u_x;
}

Real
Jacobian_2D::get_u_y () const
{
   return u_y;
}

Real
Jacobian_2D::get_v_x () const
{
   return v_x;
}

Real
Jacobian_2D::get_v_y () const
{
   return v_y;
}

const Real
Jacobian_2D::get_divergence () const
{
   return u_x + v_y;
}

const Real
Jacobian_2D::get_determinant () const
{
   return u_x * v_y - u_y * v_x;
}

void
Grid_nD::init () throw (std::bad_alloc)
{
   spacings = new Real[n];
   periodics = new bool[n];
   coordinate_tuples = new Tuple[n];
}

void
Grid_nD::init (const Tuple* coordinate_tuples,
               const Real* spacings,
               const bool* periodics)
{

   init ();

   for (Integer i = 0 ; i < n; i++)
   {

      this->spacings[i] = spacings[i];
      this->periodics[i] = periodics[i];
      this->coordinate_tuples[i] = coordinate_tuples[i];

      size_nd.buffer[i] = coordinate_tuples[i].size ();

   }

}

Grid_nD::Grid_nD (const Integer n)
             : n (n),
         size_nd (n)
{
   init ();
}

Grid_nD::Grid_nD (const Grid_nD& grid_nd)
             : n (grid_nd.n),
         size_nd (grid_nd.n)
{
   init (grid_nd.coordinate_tuples, grid_nd.spacings, grid_nd.periodics);
}

Grid_nD::~Grid_nD ()
{
   delete[] spacings;
   delete[] periodics;
   delete[] coordinate_tuples;
}

Integer
Grid_nD::size () const
{

   Integer size = 1;

   for (Integer i = 0; i < n; i++)
   {
      size *= coordinate_tuples[i].size ();
   }

   return size;

}

Tuple&
Grid_nD::get_coordinate_tuple (const Integer dimension) const
{
   return coordinate_tuples[dimension];
}

const Real&
Grid_nD::get_coordinate (const Integer dimension,
                         const Integer node) const
{
   return coordinate_tuples[dimension][node];
}

Integer
Grid_nD::get_node (const Tuple& coordinate_tuple,
                   const Real spacing,
                   const Real coordinate)
{

   Integer node;
   const Real& x = coordinate;
   const Tuple& tuple = coordinate_tuple;

   Integer size = tuple.size ();
   Real span = tuple.back () - tuple.front ();
   Real start_diff = x - tuple.front ();
   Real end_diff   = x - tuple.back ();

   if (start_diff * end_diff >= 0)
   {
      if (fabs (start_diff) < fabs (end_diff))
      {
         return 0;
      }
      else
      {
         return size - 2;
      }
   }

   if (gsl_isnan (spacing))
   {
      for (node = 0; node < size - 2; node++)
      {
         if ((x - tuple[node]) * (x - tuple[node+1]) <= 0)
         {
            break;
         }
      }
   }
   else
   {
      Real u = (x - tuple.front ()) / span;
      node = Integer (floor (u * (size - 1)));
   }

   return node;

}

Integer
Grid_nD::get_node (const Integer dimension,
                   const Real coordinate) const
{
   const Real& spacing = get_spacing (dimension);
   const Tuple& coordinate_tuple = get_coordinate_tuple (dimension);
   return get_node (coordinate_tuple, spacing, coordinate);
}

Integer
Grid_nD::get_nearest_node (const Tuple& coordinate_tuple,
                           const Real spacing,
                           const Real coordinate)
{

   Integer node;
   const Real& x = coordinate;
   const Tuple& tuple = coordinate_tuple;
   
   Integer size = tuple.size ();
   const Real start_diff = x - tuple.front ();
   Real end_diff   = x - tuple.back ();
   
   if (start_diff * end_diff >= 0)
   {
      if (fabs (start_diff) < fabs (end_diff))
      {
         return 0;
      }
      else
      {
         return size - 1;
      }
   }
   
   if (gsl_isnan (spacing))
   {
      Real min_distance = GSL_POSINF;
      const Real distance = fabs (x - tuple[node]);
      for (node = 0; node < size - 1; node++)
      {

         const Real distance = fabs (x - tuple[node]);
         if (distance < min_distance)
         {
            min_distance = distance;
            continue;
         }
         else
         {
            return node - 1;
         }
      }
   }
   else
   {
      const Real span = tuple.back () - tuple.front ();
      const Real u = (x - tuple.front ()) / span;
      node = Integer (round (u * (size - 1)));
   }

   return node;

}

Integer
Grid_nD::get_nearest_node (const Integer dimension,
                           const Real coordinate) const
{
   const Real& spacing = get_spacing (dimension);
   const Tuple& coordinate_tuple = get_coordinate_tuple (dimension);
   return get_nearest_node (coordinate_tuple, spacing, coordinate);
}

const Real&
Grid_nD::get_spacing (const Integer dimension) const
{
   return spacings[dimension];
}

Real
Grid_nD::get_spacing (const Integer dimension,
                      const Integer node) const
{
   const Tuple& tuple = get_coordinate_tuple (dimension);
   return (tuple[node + 1] - tuple[node]);
}

bool
Grid_nD::node_out_of_bounds (const Integer dimension,
                             const Integer i) const
{
   return (i < 0 || i >= get_coordinate_tuple (dimension).size ());
}

bool
Grid_nD::out_of_bounds (const Integer dimension,
                        const Real coordinate) const
{
   if (periodics[dimension]) { return false; }
   const Tuple& tuple = get_coordinate_tuple (dimension);
   return ((coordinate - tuple.front ()) * (coordinate - tuple.back ())) > 0;
}

void
Grid_nD::translate (const Integer dimension,
                    const Real delta)
{

   Tuple& tuple = get_coordinate_tuple (dimension);

   for (Tuple::iterator i = tuple.begin (); i != tuple.end (); i++)
   {
      *i += delta;
   }

}

void
Grid_nD::standardize_node (const Integer dimension,
                           Integer& node) const
{
   const Tuple& coordinate_tuple = coordinate_tuples[dimension];
   node = imodulo (node, coordinate_tuple.size ());
}

void
Grid_nD::standardize_coordinate (const Integer dimension,
                                 Real& coordinate) const
{
   const Tuple& coordinate_tuple = coordinate_tuples[dimension];
   const Real start_coordinate = coordinate_tuple.front ();
   const Real end_coordinate = coordinate_tuple.back ();
   coordinate = modulo (coordinate, start_coordinate, end_coordinate);
}

Chunk::Chunk ()
   : chunk_size (0),
     buffer (NULL)
{
}

Chunk::~Chunk ()
{
   delete[] buffer;
   buffer = NULL;
}

void
Chunk::read (ifstream& file,
             const bool float_length)
{

   if (float_length)
   {

      float* temp_buffer = new float[chunk_size];
      file.read ((char*)temp_buffer, sizeof (float) * chunk_size);

      for (size_t i = 0; i < chunk_size; i++)
      {
         float& datum = temp_buffer[i];
#ifndef WORDS_BIGENDIAN
         swap_endian (&datum, sizeof (float));
#endif
         buffer[i] = Real (datum);
      }

      delete[] temp_buffer;

   }
   else
   {

#ifndef WORDS_BIGENDIAN
      Real* temp_buffer = new Real[chunk_size];
      file.read ((char*)temp_buffer, sizeof (Real) * chunk_size);

      for (size_t i = 0; i < chunk_size; i++)
      {
         Real& datum = temp_buffer[i];
         swap_endian (&datum, sizeof (Real));
         buffer[i] = datum;
      }

      delete[] temp_buffer;
#else
      file.read ((char*)temp_buffer, sizeof (Real) * chunk_size);
#endif

   }

}

void
Chunk::write (ofstream& file,
              const bool float_length) const
{


   if (float_length)
   {

      float* temp_buffer = new float[chunk_size];

      for (size_t i = 0; i < chunk_size; i++)
      {
         float datum = float (buffer[i]);
#ifndef WORDS_BIGENDIAN
         swap_endian (&datum, sizeof (float));
#endif
         temp_buffer[i] = datum;
      }

      file.write ((char*)temp_buffer, sizeof (float) * chunk_size);
      delete[] temp_buffer;

   }
   else
   {

#ifndef WORDS_BIGENDIAN
      Real* temp_buffer = new Real[chunk_size];

      for (size_t i = 0; i < chunk_size; i++)
      {
         Real datum = buffer[i];
         swap_endian (&datum, sizeof (Real));
         temp_buffer[i] = datum;
      }

      file.write ((char*)temp_buffer, sizeof (Real) * chunk_size);
      delete[] temp_buffer;
#else
      file.write ((char*)temp_buffer, sizeof (Real) * chunk_size);
#endif

   }

}

void
Chunk::init (const Integer chunk_size)
{
   this->chunk_size = chunk_size;
   buffer = new Real[chunk_size];
}

void
Chunk::copy (const Chunk& chunk,
             const Integer address)
{
   const Integer buffer_size = sizeof (Real) * chunk.chunk_size;
   memcpy (buffer + address, chunk.buffer, buffer_size);
}

void
Chunk::set (const Integer i,
            const Real datum)
{
   buffer[i] = datum;
}

const Real&
Chunk::get (const Integer i) const
{
   return buffer[i];
}

Real& 
Chunk::get (const Integer i)
{
   return buffer[i];
}

void
Chunk::initialize (const Real datum,
                   const size_t address,
                   const size_t size)
{

   const size_t ss = std::min (address, chunk_size - 1);
   const size_t nn = std::min (address + size, chunk_size);

   for (size_t i = ss; i < nn; i++)
   {
      buffer[i] = datum;
   }

}

void
Chunk::scale_offset (const Real scale,
                     const Real offset,
                     const size_t address,
                     const size_t size)
{

   const size_t ss = std::min (address, chunk_size - 1);
   const size_t nn = std::min (address + size, chunk_size);

   for (size_t i = ss; i < nn; i++)
   {
      buffer[i] = buffer[i] * scale + offset;
   }

}

Domain_1D
Chunk::get_max_min (const size_t address,
                    const size_t size) const
{

   Domain_1D max_min (GSL_POSINF, GSL_NEGINF);
   Real& min = max_min.start;
   Real& max = max_min.end;

   const size_t ss = std::min (address, chunk_size - 1);
   const size_t nn = std::min (address + size, chunk_size);

   for (size_t i = ss; i < nn; i++)
   {
      const Real& datum = buffer[i];
      if (!gsl_finite (datum)) { continue; }
      if (datum > max) { max = datum; }
      if (datum < min) { min = datum; }
   }

   return max_min;

}

Real
Chunk::get_mean (const size_t address,
                 const size_t size) const
{

   Integer nnn = 0;
   Real sigma = 0;

   const size_t ss = std::min (address, chunk_size - 1);
   const size_t nn = std::min (address + size, chunk_size);

   for (size_t i = ss; i < nn; i++)
   {
      const Real& datum = buffer[i];
      if (!gsl_finite (datum)) { continue; }
      sigma += datum;
      nnn++;
   }

   return sigma / nnn;

}

Real
Chunk::subtract_mean (const size_t address,
                      const size_t size)
{
   const Real mean = get_mean (address, size);
   const size_t ss = std::min (address, chunk_size - 1);
   const size_t nn = std::min (address + size, chunk_size);
   for (size_t i = ss; i < nn; i++) { buffer[i] -= mean; }
}

void
Vector_Data_nD::init (const Tuple* coordinate_tuples,
                      const Real* spacings,
                      const bool* periodics)
{
   Grid_nD::init (coordinate_tuples, spacings, periodics);
}

Vector_Data_nD::Vector_Data_nD (const Integer vector_size,
                                const Integer n)
   : Grid_nD (n),
     vector_size (vector_size)
{
}

Vector_Data_nD::~Vector_Data_nD ()
{
}

Chunk*
Vector_Data_nD::get_chunk_ptr (const Integer vector_element) const
{

   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;

   Chunk* chunk_ptr = new Chunk ();
   chunk_ptr->init (size);
   chunk_ptr->copy (*this, address);
   return chunk_ptr;

}

const Integer&
Vector_Data_nD::get_vector_size () const
{
   return vector_size;
}

void
Vector_Data_nD::initialize_all (const Real datum)
{
   for (Integer i = 0; i < vector_size; i++)
   {
      initialize (i, datum);
   }
}

Real
Vector_Data_nD::get_epsilon (const Integer vector_element) const
{
   const Domain_1D& max_min = get_max_min (vector_element);
   const Real max = max_min.start;
   const Real min = max_min.end;
   const Real abs = std::max (fabs (max), fabs (min));
   return abs * 1e-5;
}

void
Vector_Data_nD::subtract_mean (const Integer vector_element)
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;
   Chunk::subtract_mean (address, size);
}

void
Vector_Data_1D::init (const Tuple& coordinate_tuple,
                      const Real spacing,
                      const bool periodic)
{

   Grid_nD::init ();

   spacings[0] = spacing;
   periodics[0] = periodic;
   coordinate_tuples[0] = coordinate_tuple;

   gia_ptrs = new gsl_interp_accel*[vector_size];
   gs_ptrs = new gsl_spline*[vector_size];

   Vector_Data_nD::init (coordinate_tuples, spacings, periodics);
   Chunk::init (vector_size * Grid_nD::size ());

   for (Integer i = 0; i < vector_size; i++)
   {
      gia_ptrs[i] = NULL;
      gs_ptrs[i] = NULL;
   }

}

Vector_Data_1D::Vector_Data_1D (const Integer vector_size,
                                const Integer size_1d,
                                const Domain_1D& domain_1d,
                                const bool periodic)
              : Vector_Data_nD (vector_size, 1)
{

   Real spacing = domain_1d.get_span () / (size_1d - 1);
   Tuple coordinate_tuple (size_1d, domain_1d.start, domain_1d.end);

   init (coordinate_tuple, spacing, periodic);

}

Vector_Data_1D::Vector_Data_1D (const Integer vector_size,
                                const Tuple coordinate_tuple,
                                const bool periodic)
              : Vector_Data_nD (vector_size, 1)
{
   init (coordinate_tuple, GSL_NAN, periodic);
}

Vector_Data_1D::~Vector_Data_1D ()
{

   for (Integer i = 0; i < vector_size; i++)
   {
      gsl_spline* gs_ptr = gs_ptrs[i];
      gsl_interp_accel* gia_ptr = gia_ptrs[i];
      if (gs_ptr != NULL) { gsl_spline_free (gs_ptr); }
      if (gia_ptr != NULL) { gsl_interp_accel_free (gia_ptr); }
   }

   //this is useless... delete?
   //if (size () < 2) { return; }

}

void
Vector_Data_1D::modify_coordinate_tuple (const Integer node,
                                         const Real coordinate)
{
   Tuple& tuple = get_coordinate_tuple (0);
   tuple[node] = coordinate;
}

void
Vector_Data_1D::set_interpolation (const gsl_interp_type* interp_type)
{
   for (Integer i = 0; i < vector_size; i++)
   {
      set_interpolation (i, interp_type);
   }
}

void
Vector_Data_1D::set_interpolation (const Integer vector_element,
                                   const gsl_interp_type* interp_type)
{


   const Integer& n = Grid_nD::size ();

   //if (n < gsl_interp_min_size (interp_type)) { return; }

   double* x = new double[n];
   double* y = new double[n];

   gsl_interp_accel*& gia_ptr = gia_ptrs[vector_element];
   gsl_spline*& gs_ptr = gs_ptrs[vector_element];

   gia_ptrs[vector_element] = gsl_interp_accel_alloc ();
   gs_ptrs[vector_element] = gsl_spline_alloc (interp_type, n);

   for (Integer j = 0; j < n; j++)
   {
      x[j] = get_coordinate (0, j);
      y[j] = get_datum (vector_element, j);
   }

   gsl_spline_init (gs_ptrs[vector_element], x, y, n);

   delete[] x;
   delete[] y;

}

void
Vector_Data_1D::initialize (const Integer vector_element,
                            const Real datum)
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;
   Chunk::initialize (datum, address, size);
}

void
Vector_Data_1D::scale_offset (const Integer vector_element,
                              const Real scale,
                              const Real offset)
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;
   Chunk::scale_offset (scale, offset, address, size);
}

Domain_1D
Vector_Data_1D::get_max_min (const Integer vector_element) const
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;
   return Chunk::get_max_min (address, size);
}

void
Vector_Data_1D::set_datum (const Integer vector_element,
                           const Integer node,
                           const Real datum)
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size + node;
   Chunk::set (address, datum);
}

const Real&
Vector_Data_1D::get_datum (const Integer vector_element,
                           const Integer node) const
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size + node;
   return Chunk::get (address);
}

Real&
Vector_Data_1D::get_datum (const Integer vector_element,
                           const Integer node)
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size + node;
   return Chunk::get (address);
}

Real
Vector_Data_1D::evaluate (const Integer vector_element,
                          const Real coordinate,
                          const Evaluate_Op evaluate_op) const
{

   Real x = coordinate;
   if (periodics[0]) { standardize_coordinate (0, x); }

   switch (evaluate_op)
   {

      case VALUE:
      {
         return gsl_spline_eval (gs_ptrs[vector_element],
            x, gia_ptrs[vector_element]);
         break;
      }

      case DX:
      {
         return gsl_spline_eval_deriv (gs_ptrs[vector_element],
            x, gia_ptrs[vector_element]);
         break;
      }

      case DX2:
      {
         return gsl_spline_eval_deriv2 (gs_ptrs[vector_element],
            x, gia_ptrs[vector_element]);
         break;
      }

      default:
      {
         const string function_str ("Vector_Data_1D::evaluate (): ");
         const string code_str = string_render ("%d %f %d",
            vector_element, coordinate, evaluate_op);
         throw Exception (function_str + code_str);
      }

   }

}

Bicubic_Coefficients::Bicubic_Coefficients (const Size_2D& size_2d)
                                 : size_2d (size_2d)
{
   data = new Real**[size_2d.i];
   for (Integer i = 0; i < size_2d.i; i++)
   {
      data[i] = new Real*[size_2d.j];
      for (Integer j = 0; j < size_2d.j; j++)
      {
         data[i][j] = new Real[16];
      }
   }
}

Bicubic_Coefficients::~Bicubic_Coefficients ()
{
   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         delete[] data[i][j];
      }
      delete[] data[i];
   }
   delete[] data;
}

Real*
Bicubic_Coefficients::get_array (const Integer node_x,
                                 const Integer node_y) const
{
   return data[node_x][node_y];
}

Tricubic_Coefficients::Tricubic_Coefficients (const Size_3D& size_3d)
                                 : size_3d (size_3d)
{
   data = new Real***[size_3d.k];
   for (Integer k = 0; k < size_3d.k; k++)
   {
      data[k] = new Real**[size_3d.i];
      for (Integer i = 0; i < size_3d.i; i++)
      {
         data[k][i] = new Real*[size_3d.j];
         for (Integer j = 0; j < size_3d.j; j++)
         {
            data[k][i][j] = new Real[64];
         }
      }
   }
}

Tricubic_Coefficients::~Tricubic_Coefficients ()
{

   for (Integer k = 0; k < size_3d.k; k++)
   {
      for (Integer i = 0; i < size_3d.i; i++)
      {
         for (Integer j = 0; j < size_3d.j; j++)
         {
            delete[] data[k][i][j];
         }
         delete[] data[k][i];
      }
      delete[] data[k];
   }
   delete[] data;
}

Real*
Tricubic_Coefficients::get_array (const Integer node_z,
                                  const Integer node_x,
                                  const Integer node_y) const
{
   return data[node_z][node_x][node_y];
}

Bicubic_Coefficients*
Vector_Data_2D::get_bicubic_coefficients_ptr (const Chunk& chunk,
                                              const Chunk& chunk_x,
                                              const Chunk& chunk_y,
                                              const Chunk& chunk_xy) const
{

   const Integer size_i = coordinate_tuples[0].size ();
   const Integer size_j = coordinate_tuples[1].size ();

   Bicubic_Coefficients* bicubic_coefficients_ptr =
      new Bicubic_Coefficients (Size_2D (size_i - 1, size_j - 1));

   Real ff[16];
   const Real B[16][16] =
    { { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      {-3, 3, 0, 0,-2,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 2,-2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0,-3, 3, 0, 0,-2,-1, 0, 0 },
      { 0, 0, 0, 0, 0, 0, 0, 0, 2,-2, 0, 0, 1, 1, 0, 0 },
      {-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0,-3, 0, 3, 0, 0, 0, 0, 0,-2, 0,-1, 0 },
      { 9,-9,-9, 9, 6, 3,-6,-3, 6,-6, 3,-3, 4, 2, 2, 1 },
      {-6, 6, 6,-6,-3,-3, 3, 3,-4, 4,-2, 2,-2,-2,-1,-1 },
      { 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0 },
      { 0, 0, 0, 0, 2, 0,-2, 0, 0, 0, 0, 0, 1, 0, 1, 0 },
      {-6, 6, 6,-6,-4,-2, 4, 2,-3, 3,-3, 3,-2,-1,-2,-1 },
      { 4,-4,-4, 4, 2, 2,-2,-2, 2,-2, 2,-2, 1, 1, 1, 1 } };

   for (Integer i = 0; i < size_i - 1; i++)
   {

      const Integer ip1 = i + 1;
      const Real delta_x = get_spacing (0, i);

      for (Integer j = 0; j < size_j - 1; j++)
      {

         const Integer jp1 = j + 1;
         const Real delta_y = get_spacing (1, j);
         const Real delta_xy = delta_x * delta_y;

         const Integer i0j0 = i * size_j + j;
         const Integer i0j1 = i * size_j + jp1;
         const Integer i1j0 = ip1 * size_j + j;
         const Integer i1j1 = ip1 * size_j + jp1;

         ff[0] = chunk.get (i0j0);
         ff[1] = chunk.get (i0j1);
         ff[2] = chunk.get (i1j0);
         ff[3] = chunk.get (i1j1);

         ff[4]  = chunk_y.get (i0j0) * delta_y;
         ff[5]  = chunk_y.get (i0j1) * delta_y;
         ff[6]  = chunk_y.get (i1j0) * delta_y;
         ff[7]  = chunk_y.get (i1j1) * delta_y;

         ff[8]  = chunk_x.get (i0j0) * delta_x;
         ff[9]  = chunk_x.get (i0j1) * delta_x;
         ff[10] = chunk_x.get (i1j0) * delta_x;
         ff[11] = chunk_x.get (i1j1) * delta_x;

         ff[12] = chunk_xy.get (i0j0) * delta_xy;
         ff[13] = chunk_xy.get (i0j1) * delta_xy;
         ff[14] = chunk_xy.get (i1j0) * delta_xy;
         ff[15] = chunk_xy.get (i1j1) * delta_xy;

         Real* a = bicubic_coefficients_ptr->get_array (i, j);

         for (Integer ii = 0; ii < 16; ii++)
         {
            a[ii] = 0;
            for (Integer jj = 0; jj < 16; jj++)
            {
               const Real b = B[ii][jj];
               if (b != 0) { a[ii] += b * ff[jj]; }
            }
         }

      }
   }

   return bicubic_coefficients_ptr;

}

int
Vector_Data_2D::gsl_multiroot_f (const gsl_vector* x,
                                 void* params,
                                 gsl_vector* f)
{

   const Vector_Data_2D& data_2d = *((const Vector_Data_2D*)params);

   const Real xx = gsl_vector_get (x, 0);
   const Real yy = gsl_vector_get (x, 1);

   for (Integer i = 0; i < data_2d.n; i++)
   {
      gsl_vector_set (f, i, data_2d.evaluate (i, xx, yy, VALUE));
   }

   return GSL_SUCCESS;

}

int
Vector_Data_2D::gsl_multiroot_df (const gsl_vector* x,
                                  void* params,
                                  gsl_matrix* J)
{

   const Vector_Data_2D& data_2d = *((const Vector_Data_2D*)params);

   const Real xx = gsl_vector_get (x, 0);
   const Real yy = gsl_vector_get (x, 1);

   for (Integer i = 0; i < data_2d.n; i++)
   {
      gsl_matrix_set (J, i, 0, data_2d.evaluate (i, xx, yy, DX));
      gsl_matrix_set (J, i, 1, data_2d.evaluate (i, xx, yy, DY));
   }

   return GSL_SUCCESS;

}

int
Vector_Data_2D::gsl_multiroot_fdf (const gsl_vector* x,
                                   void* params,
                                   gsl_vector* f,
                                   gsl_matrix* J)
{

   const Vector_Data_2D& data_2d = *((const Vector_Data_2D*)params);

   const Real xx = gsl_vector_get (x, 0);
   const Real yy = gsl_vector_get (x, 1);

   for (Integer i = 0; i < data_2d.n; i++)
   {
      gsl_vector_set (f, i, data_2d.evaluate (i, xx, yy, VALUE));
      gsl_matrix_set (J, i, 0, data_2d.evaluate (i, xx, yy, DX));
      gsl_matrix_set (J, i, 1, data_2d.evaluate (i, xx, yy, DY));
   }

   return GSL_SUCCESS;

}

void
Vector_Data_2D::init (const Tuple& coordinate_tuple_x,
                      const Tuple& coordinate_tuple_y,
                      const Real spacing_x,
                      const Real spacing_y,
                      const bool periodic_x,
                      const bool periodic_y)

{

   Grid_nD::init ();

   spacings[0] = spacing_x;
   spacings[1] = spacing_y;
   periodics[0] = periodic_x;
   periodics[1] = periodic_y;
   coordinate_tuples[0] = coordinate_tuple_x;
   coordinate_tuples[1] = coordinate_tuple_y;

   Vector_Data_nD::init (coordinate_tuples, spacings, periodics);
   Chunk::init (vector_size * Grid_nD::size ());

}

Chunk*
Vector_Data_2D::get_chunk_x_ptr (const Chunk& chunk)
{

   const Tuple& tuple_x = coordinate_tuples[0];
   const Tuple& tuple_y = coordinate_tuples[1];
   const bool& periodic_x = periodics[0];

   const Integer ni = tuple_x.size ();
   const Integer nj = tuple_y.size ();
   const Integer n = ni * nj;
   Chunk* chunk_x_ptr = new Chunk ();
   chunk_x_ptr->init (n);

   for (Integer j = 0; j < nj; j++)
   {

      Scalar_Data_1D scalar_data_1d (tuple_x, periodic_x);

      for (Integer i = 0; i < ni; i++)
      {
         const Real& datum = chunk.get (i * nj + j);
         scalar_data_1d.set_datum (i, datum);
      }

      scalar_data_1d.set_interpolation ();

      for (Integer i = 0; i < ni; i++)
      {
         Real derivative = scalar_data_1d.evaluate (tuple_x[i], DX);
         chunk_x_ptr->set (i * nj + j, derivative);
      }

   }

   return chunk_x_ptr;

}

Chunk*
Vector_Data_2D::get_chunk_y_ptr (const Chunk& chunk)
{

   const Tuple& tuple_x = coordinate_tuples[0];
   const Tuple& tuple_y = coordinate_tuples[1];
   const bool& periodic_y = periodics[1];

   const Integer ni = tuple_x.size ();
   const Integer nj = tuple_y.size ();
   const Integer n = ni * nj;
   Chunk* chunk_y_ptr = new Chunk ();
   chunk_y_ptr->init (n);

   for (Integer i = 0; i < ni; i++)
   {

      Scalar_Data_1D scalar_data_1d (tuple_y, periodic_y);
      const Integer ii = i * nj;

      for (Integer j = 0; j < nj; j++)
      {
         const Real& datum = chunk.get (ii + j);
         scalar_data_1d.set_datum (j, datum);
      }

      scalar_data_1d.set_interpolation ();

      for (Integer j = 0; j < nj; j++)
      {
         Real derivative = scalar_data_1d.evaluate (tuple_y[j], DX);
         chunk_y_ptr->set (ii + j, derivative);
      }

   }

   return chunk_y_ptr; 

}

Real
Vector_Data_2D::evaluate_bilinear (const Integer vector_element,
                                   const Real coordinate_x,
                                   const Real coordinate_y,
                                   const Evaluate_Op evaluate_op) const
{

   if (evaluate_op != VALUE && evaluate_op != DX &&
       evaluate_op != DY && evaluate_op != DXY)
   {
      return 0;
   }

   const Integer i = get_node (0, coordinate_x);
   const Integer j = get_node (1, coordinate_y);

   const Real dx = get_spacing (0, i);
   const Real dy = get_spacing (1, j);

   const Real value_00 = get_datum (vector_element, i, j);
   const Real value_10 = get_datum (vector_element, i + 1, j);
   const Real value_01 = get_datum (vector_element, i, j + 1);
   const Real value_11 = get_datum (vector_element, i + 1, j + 1);

   switch (evaluate_op)
   {

      case VALUE:
      {
         const Real u = (coordinate_x - get_coordinate (0, i)) / dx;
         const Real v = (coordinate_y - get_coordinate (1, j)) / dy;
         const Real b0 = (value_10 - value_00) * u + value_00;
         const Real b1 = (value_11 - value_01) * u + value_01;
         return (b1 - b0) * v + b0;
      }

      case DX:
      {
         const Real v = (coordinate_y - get_coordinate (1, j)) / dy;
         const Real b0 = value_10 - value_00;
         const Real b1 = value_11 - value_01;
         return ((b1 - b0) * v + b0) / dx;
      }

      case DY:
      {
         const Real u = (coordinate_x - get_coordinate (0, i)) / dx;
         const Real b0 = value_01 - value_00;
         const Real b1 = value_11 - value_10;
         return ((b1 - b0) * u + b0) / dy;
      }

      case DXY:
      {
         return (value_11 + value_00 - value_10 - value_01) / (dx * dy);
      }

   }

}

Real
Vector_Data_2D::evaluate_bicubic (const Integer vector_element,
                                  const Real coordinate_x,
                                  const Real coordinate_y,
                                  const Evaluate_Op evaluate_op) const
{ 

   const Integer i = get_node (0, coordinate_x);
   const Integer j = get_node (1, coordinate_y);

   const Real dx = get_spacing (0, i);
   const Real dy = get_spacing (1, j);

   const Real u = (coordinate_x - get_coordinate (0, i)) / dx;
   const Real v = (coordinate_y - get_coordinate (1, j)) / dy;

   const Bicubic_Coefficients& bicubic_coefficients =
      *(bicubic_coefficients_ptrs[vector_element]);
   const Real* a = bicubic_coefficients.get_array (i, j);

   switch (evaluate_op)
   {

      case LAPLACIAN:
      {
         const Real dx2 = evaluate_bicubic (a, u, v, dx, dy, DX2);
         const Real dy2 = evaluate_bicubic (a, u, v, dx, dy, DY2);
         return dx2 + dy2;
      }

      default:
      {
         return evaluate_bicubic (a, u, v, dx, dy, evaluate_op);
      }

   }

}

Real
Vector_Data_2D::evaluate_bicubic (const Real* a,
                                  const Real u,
                                  const Real v,
                                  const Real dx,
                                  const Real dy,
                                  const Evaluate_Op evaluate_op) const
{

   switch (evaluate_op)
   {

      case VALUE:
      {
         const Real b0 = a[0]  + v * (a[1]  + v * (a[2]  + v * a[3]));
         const Real b1 = a[4]  + v * (a[5]  + v * (a[6]  + v * a[7]));
         const Real b2 = a[8]  + v * (a[9]  + v * (a[10] + v * a[11]));
         const Real b3 = a[12] + v * (a[13] + v * (a[14] + v * a[15]));
         return b0 + u * (b1 + u * (b2 + u * b3));
      }

      case DY:
      {
         const Real three_v = 3 * v;
         const Real b0 = a[1]  + v * (2 * a[2]  + three_v * a[3]);
         const Real b1 = a[5]  + v * (2 * a[6]  + three_v * a[7]);
         const Real b2 = a[9]  + v * (2 * a[10] + three_v * a[11]);
         const Real b3 = a[13] + v * (2 * a[14] + three_v * a[15]);
         return (b0 + u * (b1 + u * (b2 + u * b3))) / dy;
      }

      case DX:
      {
         const Real three_u = 3 * u;
         const Real b0 = a[4] + u * (2 * a[8]  + three_u * a[12]);
         const Real b1 = a[5] + u * (2 * a[9]  + three_u * a[13]);
         const Real b2 = a[6] + u * (2 * a[10] + three_u * a[14]);
         const Real b3 = a[7] + u * (2 * a[11] + three_u * a[15]);
         return (b0 + v * (b1 + v * (b2 + v * b3))) / dx;
      }

      case DY2:
      {
         const Real three_v = 3 * v;
         const Real b0 = (a[2]  + three_v * a[3]);
         const Real b1 = (a[6]  + three_v * a[7]);
         const Real b2 = (a[10] + three_v * a[11]);
         const Real b3 = (a[14] + three_v * a[15]);
         return 2 * (b0 + u * (b1 + u * (b2 + u * b3))) / (dy * dy);
      }

      case DX2:
      {
         const Real three_u = 3 * u;
         const Real b0 = (a[8]  + three_u * a[12]);
         const Real b1 = (a[9]  + three_u * a[13]);
         const Real b2 = (a[10] + three_u * a[14]);
         const Real b3 = (a[11] + three_u * a[15]);
         return 2 * (b0 + v * (b1 + v * (b2 + v * b3))) / (dx * dx);
      }

      case DXY:
      {
         const Real three_v = 3 * v;
         const Real b1 = a[5]  + v * (2 * a[6]  + three_v * a[7]);
         const Real b2 = a[9]  + v * (2 * a[10] + three_v * a[11]);
         const Real b3 = a[13] + v * (2 * a[14] + three_v * a[15]);
         return (b1 + u * (2 * b2 + u * 3 * b3)) / (dx * dy);
      }

   }

}

Real
Vector_Data_2D::get_dmagnitude_dx (const Integer vector_element_u,
                                   const Integer vector_element_v,
                                   const Integer node_x,
                                   const Integer node_y,
                                   const Real magnitude) const
{

   const Real hx = spacings[0];
   const Integer nx = size_nd.buffer[0];

   if (node_x == 0)
   {

      const Integer ip1 = node_x + 1;
      const Real u_ip1 = get_datum (vector_element_u, ip1, node_y);
      const Real v_ip1 = get_datum (vector_element_v, ip1, node_y);
      const Real magnitude_ip1 = sqrt (u_ip1*u_ip1 + v_ip1*v_ip1);

      return (magnitude_ip1 - magnitude) / hx;

   }
   else
   if (node_x == nx - 1)
   {

      const Integer im1 = node_x - 1;

      const Real u_im1 = get_datum (vector_element_u, im1, node_y);
      const Real v_im1 = get_datum (vector_element_v, im1, node_y);
      const Real magnitude_im1 = sqrt (u_im1*u_im1 + v_im1*v_im1);

      return (magnitude - magnitude_im1) / hx;

   }
   else
   {

      const Integer im1 = node_x - 1;
      const Integer ip1 = node_x + 1;

      const Real u_im1 = get_datum (vector_element_u, im1, node_y);
      const Real v_im1 = get_datum (vector_element_v, im1, node_y);
      const Real magnitude_im1 = sqrt (u_im1*u_im1 + v_im1*v_im1);

      const Real u_ip1 = get_datum (vector_element_u, ip1, node_y);
      const Real v_ip1 = get_datum (vector_element_v, ip1, node_y);
      const Real magnitude_ip1 = sqrt (u_ip1*u_ip1 + v_ip1*v_ip1);

      return (magnitude_ip1 - magnitude_im1) / hx;

   }

}

Real
Vector_Data_2D::get_dmagnitude_dy (const Integer vector_element_u,
                                   const Integer vector_element_v,
                                   const Integer node_x,
                                   const Integer node_y,
                                   const Real magnitude) const
{

   const Real hy = spacings[1];
   const Integer ny = size_nd.buffer[1];

   if (node_y == 0)
   {

      const Integer jp1 = node_y + 1;
      const Real u_jp1 = get_datum (vector_element_u, node_x, jp1);
      const Real v_jp1 = get_datum (vector_element_v, node_x, jp1);
      const Real magnitude_jp1 = sqrt (u_jp1*u_jp1 + v_jp1*v_jp1);

      return (magnitude_jp1 - magnitude) / hy;

   }
   else
   if (node_y == ny - 1)
   {

      const Integer jm1 = node_y - 1;
      const Real u_jm1 = get_datum (vector_element_u, node_x, jm1);
      const Real v_jm1 = get_datum (vector_element_v, node_x, jm1);
      const Real magnitude_jm1 = sqrt (u_jm1*u_jm1 + v_jm1*v_jm1);

      return (magnitude - magnitude_jm1) / hy;

   }
   else
   {

      const Integer jm1 = node_y - 1;
      const Integer jp1 = node_y + 1;

      const Real u_jm1 = get_datum (vector_element_u, node_x, jm1);
      const Real v_jm1 = get_datum (vector_element_v, node_x, jm1);
      const Real magnitude_jm1 = sqrt (u_jm1*u_jm1 + v_jm1*v_jm1);

      const Real u_jp1 = get_datum (vector_element_u, node_x, jp1);
      const Real v_jp1 = get_datum (vector_element_v, node_x, jp1);
      const Real magnitude_jp1 = sqrt (u_jp1*u_jp1 + v_jp1*v_jp1);

      return (magnitude_jp1 - magnitude_jm1) / hy;

   }

}

Real
Vector_Data_2D::get_dmagnitude_dx (const Integer vector_element_u,
                                   const Integer vector_element_v,
                                   const Real x,
                                   const Real y,
                                   const Real magnitude) const
{
   const Real u = evaluate (vector_element_u, x, y);
   const Real v = evaluate (vector_element_v, x, y);
   const Real du_dx = evaluate (vector_element_u, x, y, DX);
   const Real dv_dx = evaluate (vector_element_v, x, y, DX);
   return 2 * (u + v) + (du_dx + dv_dx);
}

Real
Vector_Data_2D::get_dmagnitude_dy (const Integer vector_element_u,
                                   const Integer vector_element_v,
                                   const Real x,
                                   const Real y,
                                   const Real magnitude) const
{
   const Real u = evaluate (vector_element_u, x, y);
   const Real v = evaluate (vector_element_v, x, y);
   const Real du_dy = evaluate (vector_element_u, x, y, DY);
   const Real dv_dy = evaluate (vector_element_v, x, y, DY);
   return 2 * (u + v) + (du_dy + dv_dy);
}

Vector_Data_2D::Vector_Data_2D (const Integer vector_size,
                                const Size_2D& size_2d,
                                const Domain_2D& domain_2d,
                                const bool periodic_x,
                                const bool periodic_y)
              : Vector_Data_nD (vector_size, 2)
{

   const Real& start_x = domain_2d.domain_x.start;
   const Real& end_x   = domain_2d.domain_x.end;
   const Real& start_y = domain_2d.domain_y.start;
   const Real& end_y   = domain_2d.domain_y.end;

   Real spacing_x = domain_2d.get_width () / (size_2d.i - 1);
   Real spacing_y = domain_2d.get_height () / (size_2d.j - 1);
   Tuple coordinate_tuple_x (size_2d.i, start_x, end_x);
   Tuple coordinate_tuple_y (size_2d.j, start_y, end_y);

   bicubic_coefficients_ptrs = new Bicubic_Coefficients*[vector_size];
   for (Integer v = 0; v < vector_size; v++)
   {
      bicubic_coefficients_ptrs[v] = NULL;
   }

   init (coordinate_tuple_x, coordinate_tuple_y,
      spacing_x, spacing_y, periodic_x, periodic_y);

}

Vector_Data_2D::Vector_Data_2D (const Integer vector_size,
                                const Tuple coordinate_tuple_x,
                                const Tuple coordinate_tuple_y,
                                const bool periodic_x,
                                const bool periodic_y)
              : Vector_Data_nD (vector_size, 2)
{

   bicubic_coefficients_ptrs = new Bicubic_Coefficients*[vector_size];
   for (Integer v = 0; v < vector_size; v++)
   {
      bicubic_coefficients_ptrs[v] = NULL;
   }

   init (coordinate_tuple_x, coordinate_tuple_y, gsl_nan (),
      gsl_nan (), periodic_x, periodic_y);

}

Vector_Data_2D::~Vector_Data_2D ()
{

   for (Integer v = 0; v < vector_size; v++)
   {
      delete bicubic_coefficients_ptrs[v];
   }

   delete[] bicubic_coefficients_ptrs;

}

void
Vector_Data_2D::set_bicubic_interpolation ()
{
   for (Integer i = 0; i < vector_size; i++)
   {
      set_bicubic_interpolation (i);
   }
}

void
Vector_Data_2D::set_bicubic_interpolation (const Integer vector_element)
{

   if (bicubic_coefficients_ptrs[vector_element] == NULL)
   {

      const Chunk* chunk_ptr = get_chunk_ptr (vector_element);
      const Chunk* chunk_x_ptr = get_chunk_x_ptr (*chunk_ptr);
      const Chunk* chunk_y_ptr = get_chunk_y_ptr (*chunk_ptr);
      const Chunk* chunk_xy_ptr = get_chunk_x_ptr (*chunk_y_ptr);

      bicubic_coefficients_ptrs[vector_element] =
         get_bicubic_coefficients_ptr (*chunk_ptr, *chunk_x_ptr,
            *chunk_x_ptr, *chunk_x_ptr);

      delete chunk_ptr;
      delete chunk_x_ptr;
      delete chunk_y_ptr;
      delete chunk_xy_ptr;

   }

}

bool
Vector_Data_2D::out_of_bounds (const Integer node_x,
                               const Integer node_y) const
{
   return Grid_nD::node_out_of_bounds (0, node_x) ||
          Grid_nD::node_out_of_bounds (1, node_y);
}

bool
Vector_Data_2D::out_of_bounds (const Real coordinate_x,
                               const Real coordinate_y) const
{
   return Grid_nD::out_of_bounds (0, coordinate_x) ||
          Grid_nD::out_of_bounds (1, coordinate_y);
}

Size_2D
Vector_Data_2D::get_size_2d () const
{
   return Size_2D (size_nd.buffer[0], size_nd.buffer[1]);
}

Domain_2D
Vector_Data_2D::get_domain_2d () const
{

   const Tuple& coordinate_tuple_x = get_coordinate_tuple (0);
   const Tuple& coordinate_tuple_y = get_coordinate_tuple (1);

   const Real start_x = coordinate_tuple_x.front ();
   const Real end_x = coordinate_tuple_x.back ();
   const Real start_y = coordinate_tuple_y.front ();
   const Real end_y = coordinate_tuple_y.back ();

   return Domain_2D (start_x, end_x, start_y, end_y);

}

void
Vector_Data_2D::initialize (const Integer vector_element,
                            const Real datum)
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;
   Chunk::initialize (datum, address, size);
}

void
Vector_Data_2D::scale_offset (const Integer vector_element,
                              const Real scale,
                              const Real offset)
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;
   Chunk::scale_offset (scale, offset, address, size);
}

Domain_1D
Vector_Data_2D::get_max_min (const Integer vector_element) const
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;
   return Chunk::get_max_min (address, size);
}

void
Vector_Data_2D::set_datum (const Integer vector_element,
                           const Integer node_x,
                           const Integer node_y,
                           const Real datum)
{
   const size_t index = node_x * size_nd.buffer[1] + node_y;
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size + index;
   Chunk::set (address, datum);
}

const Real&
Vector_Data_2D::get_datum (const Integer vector_element,
                           const Integer node_x,
                           const Integer node_y) const
{
   const size_t index = node_x * size_nd.buffer[1] + node_y;
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size + index;
   return Chunk::get (address);
}

Real&
Vector_Data_2D::get_datum (const Integer vector_element,
                           const Integer node_x,
                           const Integer node_y)
{
   const size_t index = node_x * size_nd.buffer[1] + node_y;
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size + index;
   return Chunk::get (address);
}

Real
Vector_Data_2D::evaluate (const Integer vector_element,
                          const Integer node_x,
                          const Integer node_y,
                          const Evaluate_Op evaluate_op) const
{

   if (out_of_bounds (node_x, node_y))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_2D: out_of_bounds Integer");
   }

   return evaluate_nocheck (vector_element, node_x, node_y, evaluate_op);

}

Real
Vector_Data_2D::evaluate (const Integer vector_element,
                          const Real coordinate_x,
                          const Real coordinate_y,
                          const Evaluate_Op evaluate_op) const
{

   if (out_of_bounds (coordinate_x, coordinate_y))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_2D: out_of_bounds Real");
   }

   return evaluate_nocheck (vector_element,
      coordinate_x, coordinate_y, evaluate_op);

}

Real
Vector_Data_2D::evaluate_2d (const Integer vector_element_u,
                             const Integer vector_element_v,
                             const Integer node_x,
                             const Integer node_y,
                             const Evaluate_Op_2D evaluate_op_2d) const
{

   if (out_of_bounds (node_x, node_y))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_2D: out_of_bounds Integer");
   }

   switch (evaluate_op_2d)
   {

      case MAGNITUDE_OP:
      {
         const Real u = evaluate (vector_element_u, node_x, node_y, VALUE);
         const Real v = evaluate (vector_element_v, node_x, node_y, VALUE);
         return sqrt (u*u + v*v);
         break;
      }

      case VORTICITY_OP:
      {
         const Real v_x = evaluate (vector_element_v, node_x, node_y, DX);
         const Real u_y = evaluate (vector_element_u, node_x, node_y, DY);
         return v_x - u_y;
         break;
      }

      case DIVERGENCE_OP:
      {
         const Real u_x = evaluate (vector_element_u, node_x, node_y, DX);
         const Real v_y = evaluate (vector_element_v, node_x, node_y, DY);
         return u_x + v_y;
         break;
      }

      case SHEAR_OP:
      {

         const Integer i = node_x;
         const Integer j = node_y;
         const Integer ve_u = vector_element_u;
         const Integer ve_v = vector_element_v;

         const Real u = get_datum (ve_u, i, j);
         const Real v = get_datum (ve_v, i, j);
         const Real speed = sqrt (u*u + v*v);

         const Real dspeed_dx = get_dmagnitude_dx (ve_u, ve_v, i, j, speed);
         const Real dspeed_dy = get_dmagnitude_dy (ve_u, ve_v, i, j, speed);
         return (v * dspeed_dx - u * dspeed_dy) / speed;

      }

      case CURVATURE_OP:
      {

         const Integer i = node_x;
         const Integer j = node_y;
         const Integer ve_u = vector_element_u;
         const Integer ve_v = vector_element_v;

         const Real u = get_datum (ve_u, i, j);
         const Real v = get_datum (ve_v, i, j);
         const Real speed = sqrt (u*u + v*v);

         const Real dspeed_dx = get_dmagnitude_dx (ve_u, ve_v, i, j, speed);
         const Real dspeed_dy = get_dmagnitude_dy (ve_u, ve_v, i, j, speed);
         const Real shear = (v * dspeed_dx - u * dspeed_dy) / speed;

         const Real v_x = evaluate (vector_element_v, node_x, node_y, DX);
         const Real u_y = evaluate (vector_element_u, node_x, node_y, DY);
         return v_x - u_y - shear;

      }

   }

   return GSL_NAN;

}

Real
Vector_Data_2D::evaluate_2d (const Integer vector_element_u,
                             const Integer vector_element_v,
                             const Real coordinate_x,
                             const Real coordinate_y,
                             const Evaluate_Op_2D evaluate_op_2d) const
{

   if (out_of_bounds (coordinate_x, coordinate_y))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_2D: out_of_bounds Real");
   }

   const Real& x = coordinate_x;
   const Real& y = coordinate_y;

   switch (evaluate_op_2d)
   {

      case MAGNITUDE_OP:
      {
         const Real u = evaluate (vector_element_u, x, y, VALUE);
         const Real v = evaluate (vector_element_v, x, y, VALUE);
         return sqrt (u*u + v*v);
         break;
      }

      case VORTICITY_OP:
      {
         const Real v_x = evaluate (vector_element_v, x, y, DX);
         const Real u_y = evaluate (vector_element_u, x, y, DY);
         return v_x - u_y;
         break;
      }

      case DIVERGENCE_OP:
      {
         const Real u_x = evaluate (vector_element_u, x, y, DX);
         const Real v_y = evaluate (vector_element_v, x, y, DY);
         return u_x + v_y;
         break;
      }

      case SHEAR_OP:
      {

         const Integer ve_u = vector_element_u;
         const Integer ve_v = vector_element_v;

         const Real u = evaluate (ve_u, x, y);
         const Real v = evaluate (ve_v, x, y);
         const Real speed = sqrt (u*u + v*v);

         const Real dspeed_dx = get_dmagnitude_dx (ve_u, ve_v, x, y, speed);
         const Real dspeed_dy = get_dmagnitude_dy (ve_u, ve_v, x, y, speed);
         return (v * dspeed_dx - u * dspeed_dy) / speed;

      }

      case CURVATURE_OP:
      {

         const Integer ve_u = vector_element_u;
         const Integer ve_v = vector_element_v;

         const Real u = evaluate (ve_u, x, y);
         const Real v = evaluate (ve_v, x, y);
         const Real speed = sqrt (u*u + v*v);

         const Real dspeed_dx = get_dmagnitude_dx (ve_u, ve_v, x, y, speed);
         const Real dspeed_dy = get_dmagnitude_dy (ve_u, ve_v, x, y, speed);
         const Real shear = (v * dspeed_dx - u * dspeed_dy) / speed;

         const Real v_x = evaluate (vector_element_v, x, y, DX);
         const Real u_y = evaluate (vector_element_u, x, y, DY);
         return v_x - u_y - shear;

      }

   }

   return GSL_NAN;

}

Real
Vector_Data_2D::evaluate_nocheck (const Integer vector_element,
                                  const Integer node_x,
                                  const Integer node_y,
                                  const Evaluate_Op evaluate_op) const
{

   Integer i = node_x;
   Integer j = node_y;

   if (periodics[0]) { standardize_node (0, i); }
   if (periodics[1]) { standardize_node (1, j); }

   if (evaluate_op != VALUE && evaluate_op != DX && evaluate_op != DY &&
       evaluate_op != DX2 && evaluate_op != DY2 && evaluate_op != DXY)
   {
      return 0;
   }

   const Real h = spacings[0];
   const Real k = spacings[1];
   const bool uniform_x = gsl_finite (h);
   const bool uniform_y = gsl_finite (k);

   switch (evaluate_op)
   {

      case VALUE:
      {
         if (periodics[0]) { standardize_node (0, i); }
         if (periodics[1]) { standardize_node (1, j); }
         return get_datum (vector_element, i, j);
      }

      case DX:
      {

         if (!uniform_x) { return GSL_NAN; }

         Integer ai = i - 1;
         Integer bi = i + 1;

         if (periodics[1]) { standardize_node (1, j); }

         if (!periodics[0])
         {
            const Integer n = get_coordinate_tuple (0).size ();
            if (i == 0)
            {
               const Real& a = get_datum (vector_element, 0, j);
               const Real& b = get_datum (vector_element, 1, j);
               const Real& c = get_datum (vector_element, 2, j);
               return Differentiation::diff_0 (a, b, c, h);
            }
            else
            if (i == n - 1)
            {
               const Real& a = get_datum (vector_element, n - 3, j);
               const Real& b = get_datum (vector_element, n - 2, j);
               const Real& c = get_datum (vector_element, n - 1, j);
               return Differentiation::diff_2 (a, b, c, h);
            }
         }
         else
         {
            standardize_node (0, ai);
            standardize_node (0, bi);
         }

         const Real a = get_datum (vector_element, ai, j);
         const Real b = get_datum (vector_element, bi, j);
         return Differentiation::diff_1 (a, b, h);

      }

      case DY:
      {

         if (!uniform_y) { return GSL_NAN; }

         Integer aj = j - 1;
         Integer bj = j + 1;

         if (periodics[0]) { standardize_node (0, i); }

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if (j == 0)
            {
               const Real& a = get_datum (vector_element, i, 0);
               const Real& b = get_datum (vector_element, i, 1);
               const Real& c = get_datum (vector_element, i, 2);
               return Differentiation::diff_0 (a, b, c, k);
            }
            else
            if (j == n - 1)
            {
               const Real& a = get_datum (vector_element, i, n - 3);
               const Real& b = get_datum (vector_element, i, n - 2);
               const Real& c = get_datum (vector_element, i, n - 1);
               return Differentiation::diff_2 (a, b, c, k);
            }
         }
         else
         {
            standardize_node (1, aj);
            standardize_node (1, bj);
         }

         const Real a = get_datum (vector_element, i, aj);
         const Real b = get_datum (vector_element, i, bj);
         return Differentiation::diff_1 (a, b, k);

      }

      case DX2:
      {

         if (!uniform_x) { return GSL_NAN; }

         Integer ai = i - 1;
         Integer bi = i;
         Integer ci = i + 1;

         if (periodics[1]) { standardize_node (1, j); }

         if (!periodics[0])
         {
            const Integer n = get_coordinate_tuple (0).size ();
            if (i == 0)
            {
               const Real& a = get_datum (vector_element, 0, j);
               const Real& b = get_datum (vector_element, 1, j);
               const Real& c = get_datum (vector_element, 2, j);
               return Differentiation::diff_2nd (a, b, c, h);
            }
            else
            if (i == n - 1)
            {
               const Real& a = get_datum (vector_element, n - 3, j);
               const Real& b = get_datum (vector_element, n - 2, j);
               const Real& c = get_datum (vector_element, n - 1, j);
               return Differentiation::diff_2nd (a, b, c, h);
            }
         }
         else
         {
            standardize_node (0, ai);
            standardize_node (0, bi);
            standardize_node (0, ci);
         }

         const Real a = get_datum (vector_element, ai, j);
         const Real b = get_datum (vector_element, bi, j);
         const Real c = get_datum (vector_element, ci, j);
         return Differentiation::diff_2nd (a, b, c, h);

      }

      case DY2:
      {

         if (!uniform_y) { return GSL_NAN; }

         Integer aj = j - 1;
         Integer bj = j;
         Integer cj = j + 1;

         if (periodics[0]) { standardize_node (0, i); }

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if (j == 0)
            {
               const Real& a = get_datum (vector_element, i, 0);
               const Real& b = get_datum (vector_element, i, 1);
               const Real& c = get_datum (vector_element, i, 2);
               return Differentiation::diff_2nd (a, b, c, k);
            }
            else
            if (j == n - 1)
            {
               const Real& a = get_datum (vector_element, i, n - 3);
               const Real& b = get_datum (vector_element, i, n - 2);
               const Real& c = get_datum (vector_element, i, n - 1);
               return Differentiation::diff_2nd (a, b, c, k);
            }
         }
         else
         {
            standardize_node (1, aj);
            standardize_node (1, bj);
            standardize_node (1, cj);
         }

         const Real a = get_datum (vector_element, i, aj);
         const Real b = get_datum (vector_element, i, bj);
         const Real c = get_datum (vector_element, i, cj);
         return Differentiation::diff_2nd (a, b, c, k);

      }

      case DXY:
      {

         if (!(uniform_x && uniform_y)) { return GSL_NAN; }

         Integer ai = i - 1;
         Integer bi = i + 1;
         Integer aj = j - 1;
         Integer bj = j + 1;

         if (!periodics[0])
         {
            const Integer n = get_coordinate_tuple (0).size ();
            if ((i == 0) || i == n - 1) { return 0; }
         }
         else
         {
            standardize_node (0, ai);
            standardize_node (0, bi);
         }

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if ((j == 0) || j == n - 1) { return 0; }
         }
         else
         {
            standardize_node (1, aj);
            standardize_node (1, bj);
         }

         const Real a = get_datum (vector_element, ai, aj);
         const Real b = get_datum (vector_element, ai, bj);
         const Real c = get_datum (vector_element, bi, aj);
         const Real d = get_datum (vector_element, bi, bj);
         return (a - b - c + d) / (4 * h * k);

      }

   }

}

Real
Vector_Data_2D::evaluate_nocheck (const Integer vector_element,
                                  const Real coordinate_x,
                                  const Real coordinate_y,
                                  const Evaluate_Op evaluate_op) const
{

   Real x = coordinate_x;
   Real y = coordinate_y;

   if (periodics[0]) { standardize_coordinate (0, x); }
   if (periodics[1]) { standardize_coordinate (1, y); }

   if (bicubic_coefficients_ptrs[vector_element] != NULL)
   {
      return evaluate_bicubic (vector_element, x, y, evaluate_op);
   }
   else
   {
      return evaluate_bilinear (vector_element, x, y, evaluate_op);
   }

}

Real
Vector_Data_2D::evaluate_nocheck (const Integer vector_element,
                                  const Point_2D& point_2d,
                                  const Evaluate_Op evaluate_op) const
{
   return evaluate_nocheck (vector_element,
      point_2d.x, point_2d.y, evaluate_op);
}

void
Vector_Data_2D::acquire_root (Real& root_coordinate_x,
                              Real& root_coordinate_y,
                              const Real residual) const
{

   gsl_multiroot_function f =
      { &Vector_Data_2D::gsl_multiroot_f, 2, (void*)this };

   gsl_vector* x = gsl_vector_alloc (2);
   gsl_vector_set (x, 0, root_coordinate_x);
   gsl_vector_set (x, 1, root_coordinate_y);

   const gsl_multiroot_fsolver_type* T = gsl_multiroot_fsolver_hybrid;
   gsl_multiroot_fsolver* s = gsl_multiroot_fsolver_alloc (T, 2);
   gsl_multiroot_fsolver_set (s, &f, x);

   for (int status = GSL_CONTINUE; status = GSL_CONTINUE; )
   {

      status = gsl_multiroot_fsolver_iterate (s);

      if (status != GSL_SUCCESS)
      {
         gsl_vector_set (s->x, 0, GSL_NAN);
         gsl_vector_set (s->x, 1, GSL_NAN);
      }

      status = gsl_multiroot_test_residual (s->f, residual);

   }

   gsl_multiroot_fsolver_free (s);
   gsl_vector_free (x);

}

Jacobian_2D
Vector_Data_2D::get_jacobian (const Integer u_index,
                              const Integer v_index,
                              const Point_2D& point_2d) const
{
   const Real u_x = evaluate (u_index, point_2d.x, point_2d.y, DX);
   const Real u_y = evaluate (u_index, point_2d.x, point_2d.y, DY);
   const Real v_x = evaluate (v_index, point_2d.x, point_2d.y, DX);
   const Real v_y = evaluate (v_index, point_2d.x, point_2d.y, DY);
   return Jacobian_2D (u_x, u_y, v_x, v_y);
}

Jacobian_2D
Vector_Data_2D::get_jacobian_nocheck (const Integer u_index,
                                      const Integer v_index,
                                      const Point_2D& point_2d) const
{
   const Real u_x = evaluate_nocheck (u_index, point_2d.x, point_2d.y, DX);
   const Real u_y = evaluate_nocheck (u_index, point_2d.x, point_2d.y, DY);
   const Real v_x = evaluate_nocheck (v_index, point_2d.x, point_2d.y, DX);
   const Real v_y = evaluate_nocheck (v_index, point_2d.x, point_2d.y, DY);
   return Jacobian_2D (u_x, u_y, v_x, v_y);
}

void
Vector_Data_2D::subtract_x_mean (const Integer vector_element)
{

   const Size_2D& size_2d = get_size_2d ();

   for (Integer j = 0; j < size_2d.j; j++)
   {
      
      Integer n = 0;
      Real sigma = 0;

      for (Integer i = 0; i < size_2d.i; i++)
      {
         const Real& datum = get_datum (vector_element, i, j);
         if (gsl_isnan (datum)) { continue; }
         sigma += datum;
         n++;
      }

      const Real mean = sigma / n;

      for (Integer i = 0; i < size_2d.i; i++)
      {
         Real& datum = get_datum (vector_element, i, j);
         datum -= mean;
      }

   }

}

void
Vector_Data_2D::subtract_y_mean (const Integer vector_element)
{

   const Size_2D& size_2d = get_size_2d ();

   for (Integer i = 0; i < size_2d.i; i++)
   {
      
      Integer n = 0;
      Real sigma = 0;

      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real& datum = get_datum (vector_element, i, j);
         if (gsl_isnan (datum)) { continue; }
         sigma += datum;
         n++;
      }

      const Real mean = sigma / n;

      for (Integer j = 0; j < size_2d.j; j++)
      {
         Real& datum = get_datum (vector_element, i, j);
         datum -= mean;
      }

   }

}

Tricubic_Coefficients*
Vector_Data_3D::get_tricubic_coefficients_ptr (const Chunk& chunk,
                                               const Chunk& chunk_y,
                                               const Chunk& chunk_x,
                                               const Chunk& chunk_xy,
                                               const Chunk& chunk_z,
                                               const Chunk& chunk_zy,
                                               const Chunk& chunk_zx,
                                               const Chunk& chunk_zxy) const
{

   Real f, fy, fx, fxy, fz, fzy, fzx, fzxy;
   const Integer size_k = coordinate_tuples[0].size ();
   const Integer size_i = coordinate_tuples[1].size ();
   const Integer size_j = coordinate_tuples[2].size ();
   const Integer size_ij = size_i * size_j;

   Tricubic_Coefficients* tricubic_coefficients_ptr =
      new Tricubic_Coefficients (Size_3D (size_k - 1, size_i - 1, size_j - 1));

   Real ff[64];
   const Real B[64][64] = {
      { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { -3, 3, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 2, -2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 9, -9, -9, 9, 0, 0, 0, 0, 6, 3, -6, -3, 0, 0, 0, 0, 6, -6, 3, -3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { -6, 6, 6, -6, 0, 0, 0, 0, -3, -3, 3, 3, 0, 0, 0, 0, -4, 4, -2, 2, 0, 0, 0, 0, -2, -2, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { -6, 6, 6, -6, 0, 0, 0, 0, -4, -2, 4, 2, 0, 0, 0, 0, -3, 3, -3, 3, 0, 0, 0, 0, -2, -1, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 4, -4, -4, 4, 0, 0, 0, 0, 2, 2, -2, -2, 0, 0, 0, 0, 2, -2, 2, -2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, -9, -9, 9, 0, 0, 0, 0, 6, 3, -6, -3, 0, 0, 0, 0, 6, -6, 3, -3, 0, 0, 0, 0, 4, 2, 2, 1, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 6, -6, 0, 0, 0, 0, -3, -3, 3, 3, 0, 0, 0, 0, -4, 4, -2, 2, 0, 0, 0, 0, -2, -2, -1, -1, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 6, -6, 0, 0, 0, 0, -4, -2, 4, 2, 0, 0, 0, 0, -3, 3, -3, 3, 0, 0, 0, 0, -2, -1, -2, -1, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -4, -4, 4, 0, 0, 0, 0, 2, 2, -2, -2, 0, 0, 0, 0, 2, -2, 2, -2, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0},
      { -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 9, -9, 0, 0, -9, 9, 0, 0, 6, 3, 0, 0, -6, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 3, -3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { -6, 6, 0, 0, 6, -6, 0, 0, -3, -3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 4, 0, 0, -2, 2, 0, 0, -2, -2, 0, 0, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, 0, 0, -1, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 9, -9, 0, 0, -9, 9, 0, 0, 6, 3, 0, 0, -6, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, -6, 0, 0, 3, -3, 0, 0, 4, 2, 0, 0, 2, 1, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 0, 0, 6, -6, 0, 0, -3, -3, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 4, 0, 0, -2, 2, 0, 0, -2, -2, 0, 0, -1, -1, 0, 0},
      { 9, 0, -9, 0, -9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0, -6, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, -6, 0, 3, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 9, 0, -9, 0, -9, 0, 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, 3, 0, -6, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 0, -6, 0, 3, 0, -3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 2, 0, 2, 0, 1, 0},
      { -27, 27, 27, -27, 27, -27, -27, 27, -18, -9, 18, 9, 18, 9, -18, -9, -18, 18, -9, 9, 18, -18, 9, -9, -12, -6, -6, -3, 12, 6, 6, 3, -18, 18, 18, -18, -9, 9, 9, -9, -12, -6, 12, 6, -6, -3, 6, 3, -12, 12, -6, 6, -6, 6, -3, 3, -8, -4, -4, -2, -4, -2, -2, -1},
      { 18, -18, -18, 18, -18, 18, 18, -18, 9, 9, -9, -9, -9, -9, 9, 9, 12, -12, 6, -6, -12, 12, -6, 6, 6, 6, 3, 3, -6, -6, -3, -3, 12, -12, -12, 12, 6, -6, -6, 6, 6, 6, -6, -6, 3, 3, -3, -3, 8, -8, 4, -4, 4, -4, 2, -2, 4, 4, 2, 2, 2, 2, 1, 1},
      { -6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, -3, 0, 3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, 4, 0, -2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, -1, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, -3, 0, 3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, 4, 0, -2, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -2, 0, -1, 0, -1, 0},
      { 18, -18, -18, 18, -18, 18, 18, -18, 12, 6, -12, -6, -12, -6, 12, 6, 9, -9, 9, -9, -9, 9, -9, 9, 6, 3, 6, 3, -6, -3, -6, -3, 12, -12, -12, 12, 6, -6, -6, 6, 8, 4, -8, -4, 4, 2, -4, -2, 6, -6, 6, -6, 3, -3, 3, -3, 4, 2, 4, 2, 2, 1, 2, 1},
      { -12, 12, 12, -12, 12, -12, -12, 12, -6, -6, 6, 6, 6, 6, -6, -6, -6, 6, -6, 6, 6, -6, 6, -6, -3, -3, -3, -3, 3, 3, 3, 3, -8, 8, 8, -8, -4, 4, 4, -4, -4, -4, 4, 4, -2, -2, 2, 2, -4, 4, -4, 4, -2, 2, -2, 2, -2, -2, -2, -2, -1, -1, -1, -1},
      { 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { -6, 6, 0, 0, 6, -6, 0, 0, -4, -2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 4, -4, 0, 0, -4, 4, 0, 0, 2, 2, 0, 0, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -6, 6, 0, 0, 6, -6, 0, 0, -4, -2, 0, 0, 4, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0, -2, -1, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, -4, 0, 0, -4, 4, 0, 0, 2, 2, 0, 0, -2, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0},
      { -6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, -2, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, -2, 0, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, -6, 0, 6, 0, 6, 0, -6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -4, 0, -2, 0, 4, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -3, 0, 3, 0, -3, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -2, 0, -1, 0, -2, 0, -1, 0},
      { 18, -18, -18, 18, -18, 18, 18, -18, 12, 6, -12, -6, -12, -6, 12, 6, 12, -12, 6, -6, -12, 12, -6, 6, 8, 4, 4, 2, -8, -4, -4, -2, 9, -9, -9, 9, 9, -9, -9, 9, 6, 3, -6, -3, 6, 3, -6, -3, 6, -6, 3, -3, 6, -6, 3, -3, 4, 2, 2, 1, 4, 2, 2, 1},
      { -12, 12, 12, -12, 12, -12, -12, 12, -6, -6, 6, 6, 6, 6, -6, -6, -8, 8, -4, 4, 8, -8, 4, -4, -4, -4, -2, -2, 4, 4, 2, 2, -6, 6, 6, -6, -6, 6, 6, -6, -3, -3, 3, 3, -3, -3, 3, 3, -4, 4, -2, 2, -4, 4, -2, 2, -2, -2, -1, -1, -2, -2, -1, -1},
      { 4, 0, -4, 0, -4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, -2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
      { 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, -4, 0, -4, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, -2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, -2, 0, 2, 0, -2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0},
      { -12, 12, 12, -12, 12, -12, -12, 12, -8, -4, 8, 4, 8, 4, -8, -4, -6, 6, -6, 6, 6, -6, 6, -6, -4, -2, -4, -2, 4, 2, 4, 2, -6, 6, 6, -6, -6, 6, 6, -6, -4, -2, 4, 2, -4, -2, 4, 2, -3, 3, -3, 3, -3, 3, -3, 3, -2, -1, -2, -1, -2, -1, -2, -1},
      { 8, -8, -8, 8, -8, 8, 8, -8, 4, 4, -4, -4, -4, -4, 4, 4, 4, -4, 4, -4, -4, 4, -4, 4, 2, 2, 2, 2, -2, -2, -2, -2, 4, -4, -4, 4, 4, -4, -4, 4, 2, 2, -2, -2, 2, 2, -2, -2, 2, -2, 2, -2, 2, -2, 2, -2, 1, 1, 1, 1, 1, 1, 1, 1} };

   for (Integer k = 0; k < size_k - 1; k++)
   {

      const Integer kp1 = k + 1;
      const Real delta_z = get_spacing (0, k);

      for (Integer i = 0; i < size_i - 1; i++)
      {

         const Integer ip1 = i + 1;
         const Real delta_x = get_spacing (1, i);
         const Real delta_zx = delta_z * delta_x;

         for (Integer j = 0; j < size_j - 1; j++)
         {

            const Integer jp1 = j + 1;
            const Real delta_y = get_spacing (2, j);
            const Real delta_zy = delta_z * delta_y;
            const Real delta_xy = delta_x * delta_y;
            const Real delta_zxy = delta_zx * delta_y;

            const Integer z0i0j0 = k * size_ij + i * size_j + j;
            const Integer z0i0j1 = k * size_ij + i * size_j + jp1;
            const Integer z0i1j0 = k * size_ij + ip1 * size_j + j;
            const Integer z0i1j1 = k * size_ij + ip1 * size_j + jp1;
            const Integer z1i0j0 = kp1 * size_ij + i * size_j + j;
            const Integer z1i0j1 = kp1 * size_ij + i * size_j + jp1;
            const Integer z1i1j0 = kp1 * size_ij + ip1 * size_j + j;
            const Integer z1i1j1 = kp1 * size_ij + ip1 * size_j + jp1;

            ff[0] = chunk.get (z0i0j0);
            ff[1] = chunk.get (z0i0j1);
            ff[2] = chunk.get (z0i1j0);
            ff[3] = chunk.get (z0i1j1);
            ff[4] = chunk.get (z1i0j0);
            ff[5] = chunk.get (z1i0j1);
            ff[6] = chunk.get (z1i1j0);
            ff[7] = chunk.get (z1i1j1);

            ff[8]  = chunk_y.get (z0i0j0) * delta_y;
            ff[9]  = chunk_y.get (z0i0j1) * delta_y;
            ff[10] = chunk_y.get (z0i1j0) * delta_y;
            ff[11] = chunk_y.get (z0i1j1) * delta_y;
            ff[12] = chunk_y.get (z1i0j0) * delta_y;
            ff[13] = chunk_y.get (z1i0j1) * delta_y;
            ff[14] = chunk_y.get (z1i1j0) * delta_y;
            ff[15] = chunk_y.get (z1i1j1) * delta_y;

            ff[16] = chunk_x.get (z0i0j0) * delta_x;
            ff[17] = chunk_x.get (z0i0j1) * delta_x;
            ff[18] = chunk_x.get (z0i1j0) * delta_x;
            ff[19] = chunk_x.get (z0i1j1) * delta_x;
            ff[20] = chunk_x.get (z1i0j0) * delta_x;
            ff[21] = chunk_x.get (z1i0j1) * delta_x;
            ff[22] = chunk_x.get (z1i1j0) * delta_x;
            ff[23] = chunk_x.get (z1i1j1) * delta_x;

            ff[24] = chunk_x.get (z0i0j0) * delta_xy;
            ff[25] = chunk_x.get (z0i0j1) * delta_xy;
            ff[26] = chunk_x.get (z0i1j0) * delta_xy;
            ff[27] = chunk_x.get (z0i1j1) * delta_xy;
            ff[28] = chunk_x.get (z1i0j0) * delta_xy;
            ff[29] = chunk_x.get (z1i0j1) * delta_xy;
            ff[30] = chunk_x.get (z1i1j0) * delta_xy;
            ff[31] = chunk_x.get (z1i1j1) * delta_xy;

            ff[32] = chunk_z.get (z0i0j0) * delta_z;
            ff[33] = chunk_z.get (z0i0j1) * delta_z;
            ff[34] = chunk_z.get (z0i1j0) * delta_z;
            ff[35] = chunk_z.get (z0i1j1) * delta_z;
            ff[36] = chunk_z.get (z1i0j0) * delta_z;
            ff[37] = chunk_z.get (z1i0j1) * delta_z;
            ff[38] = chunk_z.get (z1i1j0) * delta_z;
            ff[39] = chunk_z.get (z1i1j1) * delta_z;

            ff[40] = chunk_zy.get (z0i0j0) * delta_zy;
            ff[41] = chunk_zy.get (z0i0j1) * delta_zy;
            ff[42] = chunk_zy.get (z0i1j0) * delta_zy;
            ff[43] = chunk_zy.get (z0i1j1) * delta_zy;
            ff[44] = chunk_zy.get (z1i0j0) * delta_zy;
            ff[45] = chunk_zy.get (z1i0j1) * delta_zy;
            ff[46] = chunk_zy.get (z1i1j0) * delta_zy;
            ff[47] = chunk_zy.get (z1i1j1) * delta_zy;

            ff[48] = chunk_zx.get (z0i0j0) * delta_zx;
            ff[49] = chunk_zx.get (z0i0j1) * delta_zx;
            ff[50] = chunk_zx.get (z0i1j0) * delta_zx;
            ff[51] = chunk_zx.get (z0i1j1) * delta_zx;
            ff[52] = chunk_zx.get (z1i0j0) * delta_zx;
            ff[53] = chunk_zx.get (z1i0j1) * delta_zx;
            ff[54] = chunk_zx.get (z1i1j0) * delta_zx;
            ff[55] = chunk_zx.get (z1i1j1) * delta_zx;

            ff[56] = chunk_zxy.get (z0i0j0) * delta_zxy;
            ff[57] = chunk_zxy.get (z0i0j1) * delta_zxy;
            ff[58] = chunk_zxy.get (z0i1j0) * delta_zxy;
            ff[59] = chunk_zxy.get (z0i1j1) * delta_zxy;
            ff[60] = chunk_zxy.get (z1i0j0) * delta_zxy;
            ff[61] = chunk_zxy.get (z1i0j1) * delta_zxy;
            ff[62] = chunk_zxy.get (z1i1j0) * delta_zxy;
            ff[63] = chunk_zxy.get (z1i1j1) * delta_zxy;

            Real* a = tricubic_coefficients_ptr->get_array (k, i, j);

            for (Integer ii = 0; ii < 64; ii++)
            {
               a[ii] = 0;
               for (Integer jj = 0; jj < 64; jj++)
               {
                  const Real b = B[ii][jj];
                  if (b != 0) { a[ii] += b * ff[jj]; }
               }
            }

         }
      }
   }

   return tricubic_coefficients_ptr;

}

Integer
Vector_Data_3D::get_offset (const Integer vector_element,
                            const Integer node_z,
                            const Integer node_x,
                            const Integer node_y) const
{
   return ((vector_element * size_nd.buffer[0] + node_z) *
      size_nd.buffer[1] + node_x) * size_nd.buffer[2] + node_y;
}

void
Vector_Data_3D::init (const Tuple& coordinate_tuple_z,
                      const Tuple& coordinate_tuple_x,
                      const Tuple& coordinate_tuple_y,
                      const Real spacing_z,
                      const Real spacing_x,
                      const Real spacing_y,
                      const bool periodic_z,
                      const bool periodic_x,
                      const bool periodic_y)
{

   Grid_nD::init ();

   spacings[0] = spacing_z;
   spacings[1] = spacing_x;
   spacings[2] = spacing_y;
   periodics[0] = periodic_z;
   periodics[1] = periodic_x;
   periodics[2] = periodic_y;
   coordinate_tuples[0] = coordinate_tuple_z;
   coordinate_tuples[1] = coordinate_tuple_x;
   coordinate_tuples[2] = coordinate_tuple_y;

   Vector_Data_nD::init (coordinate_tuples, spacings, periodics);
   Chunk::init (vector_size * Grid_nD::size ());

}

Chunk*
Vector_Data_3D::get_chunk_z_ptr (const Chunk& chunk)
{

   const Tuple& tuple_z = coordinate_tuples[0];
   const Tuple& tuple_x = coordinate_tuples[1];
   const Tuple& tuple_y = coordinate_tuples[2];
   const bool& periodic_z = periodics[0];

   const Integer nk = tuple_z.size ();
   const Integer ni = tuple_x.size ();
   const Integer nj = tuple_y.size ();
   const Integer nkij = nk * ni * nj;
   const Integer nij = ni * nj;

   Chunk* chunk_z_ptr = new Chunk ();
   chunk_z_ptr->init (nkij);

   for (Integer i = 0; i < ni; i++)
   {
      for (Integer j = 0; j < nj; j++)
      {

         Scalar_Data_1D scalar_data_1d (tuple_z, periodic_z);

         for (Integer k = 0; k < nk; k++)
         {
            const Real& datum = chunk.get (k * nij + i * nj + j);
            scalar_data_1d.set_datum (k, datum);
         }

         scalar_data_1d.set_interpolation ();

         for (Integer k = 0; k < nk; k++)
         {
            const Real derivative = scalar_data_1d.evaluate (tuple_z[k], DX);
            chunk_z_ptr->set (k * nij + i * nj + j, derivative);
         }

      }
   }

   return chunk_z_ptr;

}

Chunk*
Vector_Data_3D::get_chunk_x_ptr (const Chunk& chunk)
{

   const Tuple& tuple_z = coordinate_tuples[0];
   const Tuple& tuple_x = coordinate_tuples[1];
   const Tuple& tuple_y = coordinate_tuples[2];
   const bool& periodic_x = periodics[1];

   const Integer nk = tuple_z.size ();
   const Integer ni = tuple_x.size ();
   const Integer nj = tuple_y.size ();
   const Integer nkij = nk * ni * nj;
   const Integer nij = ni * nj;

   Chunk* chunk_x_ptr = new Chunk ();
   chunk_x_ptr->init (nkij);

   for (Integer k = 0; k < nk; k++)
   {
      for (Integer j = 0; j < nj; j++)
      {

         Scalar_Data_1D scalar_data_1d (tuple_x, periodic_x);

         for (Integer i = 0; i < ni; i++)
         {
            const Real& datum = chunk.get (k * nij + i * nj + j);
            scalar_data_1d.set_datum (i, datum);
         }

         scalar_data_1d.set_interpolation ();

         for (Integer i = 0; i < ni; i++)
         {
            const Real derivative = scalar_data_1d.evaluate (tuple_x[i], DX);
            chunk_x_ptr->set (k * nij + i * nj + j, derivative);
         }

      }
   }

   return chunk_x_ptr;

}

Chunk*
Vector_Data_3D::get_chunk_y_ptr (const Chunk& chunk)
{

   const Tuple& tuple_z = coordinate_tuples[0];
   const Tuple& tuple_x = coordinate_tuples[1];
   const Tuple& tuple_y = coordinate_tuples[2];
   const bool& periodic_y = periodics[2];

   const Integer nk = tuple_z.size ();
   const Integer ni = tuple_x.size ();
   const Integer nj = tuple_y.size ();
   const Integer nkij = nk * ni * nj;
   const Integer nij = ni * nj;

   Chunk* chunk_y_ptr = new Chunk ();
   chunk_y_ptr->init (nkij);

   for (Integer k = 0; k < nk; k++)
   {
      for (Integer i = 0; i < ni; i++)
      {

         Scalar_Data_1D scalar_data_1d (tuple_y, periodic_y);

         for (Integer j = 0; j < nj; j++)
         {
            const Real& datum = chunk.get (k * nij + i * nj + j);
            scalar_data_1d.set_datum (j, datum);
         }

         scalar_data_1d.set_interpolation ();

         for (Integer j = 0; j < nj; j++)
         {
            const Real derivative = scalar_data_1d.evaluate (tuple_y[j], DX);
            chunk_y_ptr->set (k * nij + i * nj + j, derivative);
         }

      }
   }

   return chunk_y_ptr;

}

Real
Vector_Data_3D::evaluate_linear_z (const Integer vector_element,
                                   const Real coordinate_z,
                                   const Integer node_x,
                                   const Integer node_y) const
{

   const Integer ak = get_node (0, coordinate_z);
   const Integer bk = ak + 1;

   const Real dz = get_spacing (0, ak);
   const Real w = (coordinate_z - get_coordinate (0, ak)) / dz;

   const Real af = get_datum (vector_element, ak, node_x, node_y);
   const Real bf = get_datum (vector_element, bk, node_x, node_y);

   return af + (bf - af) * w;

}

Real
Vector_Data_3D::evaluate_bilinear (const Integer vector_element,
                                   const Integer node_z,
                                   const Real coordinate_x,
                                   const Real coordinate_y,
                                   const Evaluate_Op evaluate_op) const
{

   if (evaluate_op != VALUE && evaluate_op != DZ && evaluate_op != DX &&
       evaluate_op != DY && evaluate_op != DXY)
   {
      return 0;
   }

   const Integer k = node_z;
   const Integer i = get_node (1, coordinate_x);
   const Integer j = get_node (2, coordinate_y);

   const Real dx = get_spacing (1, i);
   const Real dy = get_spacing (2, j);

   switch (evaluate_op)
   {

      case VALUE:
      {
         const Real value_00 = get_datum (vector_element, k, i, j);
         const Real value_10 = get_datum (vector_element, k, i + 1, j);
         const Real value_01 = get_datum (vector_element, k, i, j + 1);
         const Real value_11 = get_datum (vector_element, k, i + 1, j + 1);
         const Real u = (coordinate_x - get_coordinate (1, i)) / dx;
         const Real v = (coordinate_y - get_coordinate (2, j)) / dy;
         const Real b0 = (value_10 - value_00) * u + value_00;
         const Real b1 = (value_11 - value_01) * u + value_01;
         return (b1 - b0) * v + b0;
      }

      case DZ:
      {

         const Real& x = coordinate_x;
         const Real& y = coordinate_y;
         const Integer ve = vector_element;
         const bool last_k = (k == size_nd.buffer[0] - 1);

         if (last_k)
         {
            const Real spacing = get_spacing (0, k - 1);
            const Real value = evaluate_bilinear (ve, k, x, y, VALUE);
            const Real alt_value = evaluate_bilinear (ve, k - 1, x, y, VALUE);
            return (value - alt_value) / spacing;
         }
         else
         {
            const Real spacing = get_spacing (0, k);
            const Real value = evaluate_bilinear (ve, k, x, y, VALUE);
            const Real alt_value = evaluate_bilinear (ve, k + 1, x, y, VALUE);
            return (alt_value - value) / spacing;
         }

      }

      case DX:
      {
         const Real value_00 = get_datum (vector_element, k, i, j);
         const Real value_10 = get_datum (vector_element, k, i + 1, j);
         const Real value_01 = get_datum (vector_element, k, i, j + 1);
         const Real value_11 = get_datum (vector_element, k, i + 1, j + 1);
         const Real v = (coordinate_y - get_coordinate (2, j)) / dy;
         const Real b0 = value_10 - value_00;
         const Real b1 = value_11 - value_01;
         return ((b1 - b0) * v + b0) / dx;
      }

      case DY:
      {
         const Real value_00 = get_datum (vector_element, k, i, j);
         const Real value_10 = get_datum (vector_element, k, i + 1, j);
         const Real value_01 = get_datum (vector_element, k, i, j + 1);
         const Real value_11 = get_datum (vector_element, k, i + 1, j + 1);
         const Real u = (coordinate_x - get_coordinate (1, i)) / dx;
         const Real b0 = value_01 - value_00;
         const Real b1 = value_11 - value_10;
         return ((b1 - b0) * u + b0) / dy;
      }

      case DXY:
      {
         const Real value_00 = get_datum (vector_element, k, i, j);
         const Real value_10 = get_datum (vector_element, k, i + 1, j);
         const Real value_01 = get_datum (vector_element, k, i, j + 1);
         const Real value_11 = get_datum (vector_element, k, i + 1, j + 1);
         return (value_11 + value_00 - value_10 - value_01) / (dx * dy);
      }

   }
}

Real
Vector_Data_3D::evaluate_trilinear (const Integer vector_element,
                                    const Real coordinate_z,
                                    const Real coordinate_x,
                                    const Real coordinate_y,
                                    const Evaluate_Op evaluate_op) const
{

   if (evaluate_op != VALUE && evaluate_op != DZ && evaluate_op != DX &&
       evaluate_op != DY && evaluate_op != DZY && evaluate_op != DZX &&
       evaluate_op != DXY && evaluate_op != DZXY)
   {
      return 0;
   }

   const Integer k = get_node (0, coordinate_z);
   const Integer i = get_node (1, coordinate_x);
   const Integer j = get_node (2, coordinate_y);

   const Real dz = get_spacing (0, k);
   const Real dx = get_spacing (1, i);
   const Real dy = get_spacing (2, j);

   const Integer kp1 = k + 1;
   const Integer ip1 = i + 1;
   const Integer jp1 = j + 1;

   const Real f_000 = get_datum (vector_element, k, i, j);
   const Real f_010 = get_datum (vector_element, k, ip1, j);
   const Real f_001 = get_datum (vector_element, k, i, jp1);
   const Real f_011 = get_datum (vector_element, k, ip1, jp1);
   const Real f_100 = get_datum (vector_element, kp1, i, j);
   const Real f_110 = get_datum (vector_element, kp1, ip1, j);
   const Real f_101 = get_datum (vector_element, kp1, i, jp1);
   const Real f_111 = get_datum (vector_element, kp1, ip1, jp1);

   switch (evaluate_op)
   {

      case VALUE:
      {

         const Real w = (coordinate_z - get_coordinate (0, k)) / dz;
         const Real u = (coordinate_x - get_coordinate (1, i)) / dx;
         const Real v = (coordinate_y - get_coordinate (2, j)) / dy;

         const Real a = f_100 + f_111 + f_001 + f_010 -
            f_101 - f_110 - f_000 - f_011;
         const Real b = f_000 + f_011 - f_001 - f_010;
         const Real c = f_101 + f_000 - f_100 - f_001;
         const Real d = f_110 + f_000 - f_100 - f_010;
         const Real e = f_001 - f_000;
         const Real f = f_010 - f_000;
         const Real g = f_100 - f_000;

         return a * w*u*v + b * u*v + c * w*v + d * w*u
              + e * v + f * u + g * w + f_000;

      }

      case DZ:
      {
         const Real u = (coordinate_x - get_coordinate (1, i)) / dx;
         const Real v = (coordinate_y - get_coordinate (2, j)) / dy;
         const Real a = f_100 + f_111 + f_001 + f_010 -
            f_101 - f_110 - f_000 - f_011;
         const Real c = f_101 + f_000 - f_100 - f_001;
         const Real d = f_110 + f_000 - f_100 - f_010;
         const Real g = f_100 - f_000;

         return (a * u*v + d * u + c * v + g) / dz;
      }

      case DX:
      {

         const Real w = (coordinate_z - get_coordinate (0, k)) / dz;
         const Real v = (coordinate_y - get_coordinate (2, j)) / dy;

         const Real a = f_100 + f_111 + f_001 + f_010 -
            f_101 - f_110 - f_000 - f_011;
         const Real b = f_000 + f_011 - f_001 - f_010;
         const Real d = f_110 + f_000 - f_100 - f_010;
         const Real f = f_010 - f_000;

         return (a * v*w + b * v + d * w + f) / dx;
      }

      case DY:
      {
         const Real w = (coordinate_z - get_coordinate (0, k)) / dz;
         const Real u = (coordinate_x - get_coordinate (1, i)) / dx;
         const Real a = f_100 + f_111 + f_001 + f_010 -
            f_101 - f_110 - f_000 - f_011;
         const Real b = f_000 + f_011 - f_001 - f_010;
         const Real c = f_101 + f_000 - f_100 - f_001;
         const Real e = f_001 - f_000;

         return (a * u*w + b * u + c * w + e) / dx;
      }

      case DZY:
      {
         const Real u = (coordinate_x - get_coordinate (1, i)) / dx;
         const Real a = f_100 + f_111 + f_001 + f_010 -
            f_101 - f_110 - f_000 - f_011;
         const Real c = f_101 + f_000 - f_100 - f_001;

         return (a * u + c) / (dz * dy);
      }

      case DZX:
      {
         const Real v = (coordinate_y - get_coordinate (2, j)) / dy;
         const Real a = f_100 + f_111 + f_001 + f_010 -
            f_101 - f_110 - f_000 - f_011;
         const Real d = f_110 + f_000 - f_100 - f_010;

         return (a * v + d) / (dz * dx);
      }

      case DXY:
      {
         const Real w = (coordinate_z - get_coordinate (0, k)) / dz;
         const Real a = f_100 + f_111 + f_001 + f_010 -
            f_101 - f_110 - f_000 - f_011;
         const Real b = f_000 + f_011 - f_001 - f_010;

         return (a * w + b) / (dx * dy);
      }

      case DZXY:
      {
         const Real a = f_100 + f_111 + f_001 + f_010 -
            f_101 - f_110 - f_000 - f_011;

         return a / (dz * dx * dy);
      }

   }

}

Real
Vector_Data_3D::evaluate_tricubic (const Integer vector_element,
                                   const Real coordinate_z,
                                   const Real coordinate_x,
                                   const Real coordinate_y,
                                   const Evaluate_Op evaluate_op) const
{ 

   const Integer k = get_node (0, coordinate_z);
   const Integer i = get_node (1, coordinate_x);
   const Integer j = get_node (2, coordinate_y);

   const Real dz = get_spacing (0, k);
   const Real dx = get_spacing (1, i);
   const Real dy = get_spacing (2, j);

   const Real w = (coordinate_z - get_coordinate (0, k)) / dz;
   const Real u = (coordinate_x - get_coordinate (1, i)) / dx;
   const Real v = (coordinate_y - get_coordinate (2, j)) / dy;

   const Tricubic_Coefficients& tricubic_coefficients =
      *(tricubic_coefficients_ptrs[vector_element]);
   const Real* a = tricubic_coefficients.get_array (k, i, j);

   switch (evaluate_op)
   {

      case LAPLACIAN:
      {
         const Real dz2 = evaluate_tricubic (a, w, u, v, dz, dx, dy, DZ2);
         const Real dx2 = evaluate_tricubic (a, w, u, v, dz, dx, dy, DX2);
         const Real dy2 = evaluate_tricubic (a, w, u, v, dz, dx, dy, DY2);
         return dz2 + dx2 + dy2;
      }

      default:
      {
         return evaluate_tricubic (a, w, u, v, dz, dx, dy, evaluate_op);
      }

   }

}

Real
Vector_Data_3D::evaluate_tricubic (const Real* a,
                                   const Real w,
                                   const Real u,
                                   const Real v,
                                   const Real dz,
                                   const Real dx,
                                   const Real dy,
                                   const Evaluate_Op evaluate_op) const
{

// all wrong
   switch (evaluate_op)
   {

      case VALUE:
      {
         const Real b0 = a[0]  + v * (a[1]  + v * (a[2]  + v * a[3]));
         const Real b1 = a[4]  + v * (a[5]  + v * (a[6]  + v * a[7]));
         const Real b2 = a[8]  + v * (a[9]  + v * (a[10] + v * a[11]));
         const Real b3 = a[12] + v * (a[13] + v * (a[14] + v * a[15]));
         return b0 + u * (b1 + u * (b2 + u * b3));
      }

      case DZ:
      {
         const Real three_u = 3 * u;
         const Real b0 = a[4] + u * (2 * a[8]  + three_u * a[12]);
         const Real b1 = a[5] + u * (2 * a[9]  + three_u * a[13]);
         const Real b2 = a[6] + u * (2 * a[10] + three_u * a[14]);
         const Real b3 = a[7] + u * (2 * a[11] + three_u * a[15]);
         return (b0 + v * (b1 + v * (b2 + v * b3))) / dz;
      }
      case DX:
      {
         const Real three_u = 3 * u;
         const Real b0 = a[4] + u * (2 * a[8]  + three_u * a[12]);
         const Real b1 = a[5] + u * (2 * a[9]  + three_u * a[13]);
         const Real b2 = a[6] + u * (2 * a[10] + three_u * a[14]);
         const Real b3 = a[7] + u * (2 * a[11] + three_u * a[15]);
         return (b0 + v * (b1 + v * (b2 + v * b3))) / dy;
      }

      case DY:
      {
         const Real three_v = 3 * v;
         const Real b0 = a[1]  + v * (2 * a[2]  + three_v * a[3]);
         const Real b1 = a[5]  + v * (2 * a[6]  + three_v * a[7]);
         const Real b2 = a[9]  + v * (2 * a[10] + three_v * a[11]);
         const Real b3 = a[13] + v * (2 * a[14] + three_v * a[15]);
         return (b0 + u * (b1 + u * (b2 + u * b3))) / dx;
      }

      case DZ2:
      {
         const Real three_u = 3 * u;
         const Real b0 = (a[8]  + three_u * a[12]);
         const Real b1 = (a[9]  + three_u * a[13]);
         const Real b2 = (a[10] + three_u * a[14]);
         const Real b3 = (a[11] + three_u * a[15]);
         return 2 * (b0 + v * (b1 + v * (b2 + v * b3))) / (dz * dz);
      }

      case DX2:
      {
         const Real three_u = 3 * u;
         const Real b0 = (a[8]  + three_u * a[12]);
         const Real b1 = (a[9]  + three_u * a[13]);
         const Real b2 = (a[10] + three_u * a[14]);
         const Real b3 = (a[11] + three_u * a[15]);
         return 2 * (b0 + v * (b1 + v * (b2 + v * b3))) / (dy * dy);
      }

      case DY2:
      {
         const Real three_v = 3 * v;
         const Real b0 = (a[2]  + three_v * a[3]);
         const Real b1 = (a[6]  + three_v * a[7]);
         const Real b2 = (a[10] + three_v * a[11]);
         const Real b3 = (a[14] + three_v * a[15]);
         return 2 * (b0 + u * (b1 + u * (b2 + u * b3))) / (dx * dx);
      }

      case DXY:
      {
         const Real b1 = a[5]  + v * (2 * a[6]  + 3 * v * a[7]);
         const Real b2 = a[9]  + v * (2 * a[10] + 3 * v * a[11]);
         const Real b3 = a[13] + v * (2 * a[14] + 3 * v * a[15]);
         return (b1 + u * (2 * b2 + u * 3 * b3)) / (dx * dy);
      }

      case DZY:
      {
         const Real b1 = a[5]  + v * (2 * a[6]  + 3 * v * a[7]);
         const Real b2 = a[9]  + v * (2 * a[10] + 3 * v * a[11]);
         const Real b3 = a[13] + v * (2 * a[14] + 3 * v * a[15]);
         return (b1 + u * (2 * b2 + u * 3 * b3)) / (dz * dy);
      }

      case DZX:
      {
         const Real b1 = a[5]  + v * (2 * a[6]  + 3 * v * a[7]);
         const Real b2 = a[9]  + v * (2 * a[10] + 3 * v * a[11]);
         const Real b3 = a[13] + v * (2 * a[14] + 3 * v * a[15]);
         return (b1 + u * (2 * b2 + u * 3 * b3)) / (dz * dx);
      }

      case DZXY:
      {
         const Real b1 = a[5]  + v * (2 * a[6]  + 3 * v * a[7]);
         const Real b2 = a[9]  + v * (2 * a[10] + 3 * v * a[11]);
         const Real b3 = a[13] + v * (2 * a[14] + 3 * v * a[15]);
         return (b1 + u * (2 * b2 + u * 3 * b3)) / (dz * dx * dy);
      }

   }

}



Vector_Data_3D::Vector_Data_3D (const Integer vector_size,
                                const Size_3D& size_3d,
                                const Domain_3D& domain_3d,
                                const bool periodic_z,
                                const bool periodic_x,
                                const bool periodic_y)
              : Vector_Data_nD (vector_size, 3)
{

   const Real& start_z = domain_3d.domain_z.start;
   const Real& end_z   = domain_3d.domain_z.end;
   const Real& start_x = domain_3d.domain_x.start;
   const Real& end_x   = domain_3d.domain_x.end;
   const Real& start_y = domain_3d.domain_y.start;
   const Real& end_y   = domain_3d.domain_y.end;

   Real spacing_z = domain_3d.get_depth () / (size_3d.k - 1);
   Real spacing_x = domain_3d.get_width () / (size_3d.i - 1);
   Real spacing_y = domain_3d.get_height () / (size_3d.j - 1);
   Tuple coordinate_tuple_z (size_3d.k, start_z, end_z);
   Tuple coordinate_tuple_x (size_3d.i, start_x, end_x);
   Tuple coordinate_tuple_y (size_3d.j, start_y, end_y);

   tricubic_coefficients_ptrs = new Tricubic_Coefficients*[vector_size];
   for (Integer v = 0; v < vector_size; v++)
   {
      tricubic_coefficients_ptrs[v] = NULL;
   }

   init (coordinate_tuple_z, coordinate_tuple_x, coordinate_tuple_y,
      spacing_z, spacing_x, spacing_y, periodic_z, periodic_x, periodic_y);

}

Vector_Data_3D::Vector_Data_3D (const Integer vector_size,
                                const Tuple coordinate_tuple_z,
                                const Tuple coordinate_tuple_x,
                                const Tuple coordinate_tuple_y,
                                const bool periodic_z,
                                const bool periodic_x,
                                const bool periodic_y)
              : Vector_Data_nD (vector_size, 3)
{

   tricubic_coefficients_ptrs = new Tricubic_Coefficients*[vector_size];
   for (Integer v = 0; v < vector_size; v++)
   {
      tricubic_coefficients_ptrs[v] = NULL;
   }

   const Real nan = gsl_nan ();

   init (coordinate_tuple_z, coordinate_tuple_x, coordinate_tuple_y,
      nan, nan, nan, periodic_z, periodic_x, periodic_y);

}

Vector_Data_3D::Vector_Data_3D (const Integer vector_size,
                                const Tuple coordinate_tuple_z,
                                const Size_2D& size_2d,
                                const Domain_2D& domain_2d,
                                const bool periodic_z,
                                const bool periodic_x,
                                const bool periodic_y)
              : Vector_Data_nD (vector_size, 3)
{

   const Real& start_x = domain_2d.domain_x.start;
   const Real& end_x   = domain_2d.domain_x.end;
   const Real& start_y = domain_2d.domain_y.start;
   const Real& end_y   = domain_2d.domain_y.end;

   Real spacing_x = domain_2d.get_width () / (size_2d.i - 1);
   Real spacing_y = domain_2d.get_height () / (size_2d.j - 1);
   Tuple coordinate_tuple_x (size_2d.i, start_x, end_x);
   Tuple coordinate_tuple_y (size_2d.j, start_y, end_y);

   tricubic_coefficients_ptrs = new Tricubic_Coefficients*[vector_size];
   for (Integer v = 0; v < vector_size; v++)
   {
      tricubic_coefficients_ptrs[v] = NULL;
   }

   const Real nan = gsl_nan ();

   init (coordinate_tuple_z, coordinate_tuple_x, coordinate_tuple_y,
      nan, spacing_x, spacing_y, periodic_z, periodic_x, periodic_y);

}

Vector_Data_3D::~Vector_Data_3D ()
{

   for (Integer v = 0; v < vector_size; v++)
   {
      delete tricubic_coefficients_ptrs[v];
   }

   delete[] tricubic_coefficients_ptrs;

}

void
Vector_Data_3D::set_tricubic_interpolation ()
{
   for (Integer i = 0; i < vector_size; i++)
   {
      set_tricubic_interpolation (i);
   }
}

void
Vector_Data_3D::set_tricubic_interpolation (const Integer vector_element)
{

   if (tricubic_coefficients_ptrs[vector_element] == NULL)
   {

      Chunk* chunk_ptr = get_chunk_ptr (vector_element);
      Chunk* chunk_z_ptr = get_chunk_z_ptr (*chunk_ptr);
      Chunk* chunk_y_ptr = get_chunk_y_ptr (*chunk_ptr);
      Chunk* chunk_x_ptr = get_chunk_x_ptr (*chunk_ptr);
      Chunk* chunk_xy_ptr = get_chunk_x_ptr (*chunk_y_ptr);
      Chunk* chunk_zx_ptr = get_chunk_z_ptr (*chunk_x_ptr);
      Chunk* chunk_zy_ptr = get_chunk_z_ptr (*chunk_y_ptr);
      Chunk* chunk_zxy_ptr = get_chunk_z_ptr (*chunk_xy_ptr);

      tricubic_coefficients_ptrs[vector_element] =
         get_tricubic_coefficients_ptr (*chunk_ptr, *chunk_y_ptr, *chunk_x_ptr,
            *chunk_xy_ptr, *chunk_z_ptr, *chunk_zy_ptr, *chunk_zx_ptr, *chunk_zxy_ptr);

      delete chunk_ptr;
      delete chunk_z_ptr;
      delete chunk_x_ptr;
      delete chunk_y_ptr;
      delete chunk_xy_ptr;
      delete chunk_zx_ptr;
      delete chunk_zy_ptr;
      delete chunk_zxy_ptr;

   }

}

bool
Vector_Data_3D::out_of_bounds (const Integer node_z,
                               const Integer node_x,
                               const Integer node_y) const
{
   return Grid_nD::node_out_of_bounds (0, node_z) ||
          Grid_nD::node_out_of_bounds (1, node_x) ||
          Grid_nD::node_out_of_bounds (2, node_y);
}

bool
Vector_Data_3D::out_of_bounds (const Integer node_z,
                               const Real coordinate_x,
                               const Real coordinate_y) const
{
   return Grid_nD::node_out_of_bounds (0, node_z) ||
          Grid_nD::out_of_bounds (1, coordinate_x) ||
          Grid_nD::out_of_bounds (2, coordinate_y);
}

bool
Vector_Data_3D::out_of_bounds (const Real coordinate_z,
                               const Integer node_x,
                               const Integer node_y) const
{
   return Grid_nD::out_of_bounds (0, coordinate_z) ||
          Grid_nD::node_out_of_bounds (1, node_x) ||
          Grid_nD::node_out_of_bounds (2, node_y);
}

bool
Vector_Data_3D::out_of_bounds (const Real coordinate_z,
                               const Real coordinate_x,
                               const Real coordinate_y) const
{
   return Grid_nD::out_of_bounds (0, coordinate_z) ||
          Grid_nD::out_of_bounds (1, coordinate_x) ||
          Grid_nD::out_of_bounds (2, coordinate_y);
}

Size_2D
Vector_Data_3D::get_size_2d () const
{
   return Size_2D (size_nd.buffer[1], size_nd.buffer[2]);
}

Domain_2D
Vector_Data_3D::get_domain_2d () const
{

   const Tuple& coordinate_tuple_x = get_coordinate_tuple (1);
   const Tuple& coordinate_tuple_y = get_coordinate_tuple (2);

   const Real start_x = coordinate_tuple_x.front ();
   const Real end_x = coordinate_tuple_x.back ();
   const Real start_y = coordinate_tuple_y.front ();
   const Real end_y = coordinate_tuple_y.back ();

   return Domain_2D (start_x, end_x, start_y, end_y);

}

Size_3D
Vector_Data_3D::get_size_3d () const
{
   return Size_3D (size_nd.buffer[0], size_nd.buffer[1], size_nd.buffer[2]);
}

void
Vector_Data_3D::initialize (const Integer vector_element,
                            const Real datum)
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;
   Chunk::initialize (datum, address, size);
}

void
Vector_Data_3D::scale_offset (const Integer vector_element,
                              const Real scale,
                              const Real offset)
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;
   Chunk::scale_offset (scale, offset, address, size);
}


Domain_1D
Vector_Data_3D::get_max_min (const Integer vector_element) const
{
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size;
   return Chunk::get_max_min (address, size);
}

void
Vector_Data_3D::set_datum (const Integer vector_element,
                           const Integer k,
                           const Integer i,
                           const Integer j,
                           const Real datum)
{
   const size_t ni = size_nd.buffer[1];
   const size_t nj = size_nd.buffer[2];
   const size_t index = (k * ni + i) * nj + j;
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size + index;
   Chunk::set (address, datum);
}

const Real&
Vector_Data_3D::get_datum (const Integer vector_element,
                           const Integer k,
                           const Integer i,
                           const Integer j) const
{
   const size_t ni = size_nd.buffer[1];
   const size_t nj = size_nd.buffer[2];
   const size_t index = (k * ni + i) * nj + j;
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size + index;
   return Chunk::get (address);
}

Real&
Vector_Data_3D::get_datum (const Integer vector_element,
                           const Integer k,
                           const Integer i,
                           const Integer j)
{
   const size_t ni = size_nd.buffer[1];
   const size_t nj = size_nd.buffer[2];
   const size_t index = (k * ni + i) * nj + j;
   const size_t size = Grid_nD::size ();
   const size_t address = vector_element * size + index;
   return Chunk::get (address);
}

Real
Vector_Data_3D::evaluate (const Integer vector_element,
                          const Integer k,
                          const Integer i,
                          const Integer j,
                          const Evaluate_Op evaluate_op) const
{

   if (out_of_bounds (k, i, j))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_3D: out_of_bounds");
   }

   return evaluate_nocheck (vector_element, k, i, j, evaluate_op);

}

Real
Vector_Data_3D::evaluate (const Integer vector_element,
                          const Integer k,
                          const Real x,
                          const Real y,
                          const Evaluate_Op evaluate_op) const
{

   if (out_of_bounds (k, x, y))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_3D: out_of_bounds");
   }

   return evaluate_nocheck (vector_element, k, x, y, evaluate_op);

}

Real
Vector_Data_3D::evaluate (const Integer vector_element,
                          const Real z,
                          const Integer i,
                          const Integer j,
                          const Evaluate_Op evaluate_op) const
{

   if (out_of_bounds (z, i, j))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_3D: out_of_bounds");
   }

   return evaluate_nocheck (vector_element, z, i, j, evaluate_op);

}

Real
Vector_Data_3D::evaluate (const Integer vector_element,
                          const Real z,
                          const Real x,
                          const Real y,
                          const Evaluate_Op evaluate_op) const
{

   if (out_of_bounds (z, x, y))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_3D: out_of_bounds");
   }

   return evaluate_nocheck (vector_element, z, x, y, evaluate_op);

}

Real
Vector_Data_3D::evaluate_uv (const Integer vector_element_u,
                             const Integer vector_element_v,
                             const Integer k,
                             const Integer i,
                             const Integer j,
                             const Evaluate_Op_2D evaluate_op_2d) const
{

   if (out_of_bounds (k, i, j))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_3D: out_of_bounds");
   }

   switch (evaluate_op_2d)
   {

      case MAGNITUDE_OP:
      {
         const Real u = evaluate (vector_element_u, k, i, j, VALUE);
         const Real v = evaluate (vector_element_v, k, i, j, VALUE);
         return sqrt (u*u + v*v);
         break;
      }

      case VORTICITY_OP:
      {
         const Real v_x = evaluate (vector_element_v, k, i, j, DX);
         const Real u_y = evaluate (vector_element_u, k, i, j, DY);
         return v_x - u_y;
         break;
      }

      case DIVERGENCE_OP:
      {
         const Real u_x = evaluate (vector_element_u, k, i, j, DX);
         const Real v_y = evaluate (vector_element_v, k, i, j, DY);
         return u_x + v_y;
         break;
      }

   }

   return GSL_NAN;

}

Real
Vector_Data_3D::evaluate_uv (const Integer vector_element_u,
                             const Integer vector_element_v,
                             const Integer k,
                             const Real x,
                             const Real y,
                             const Evaluate_Op_2D evaluate_op_2d) const
{

   if (out_of_bounds (k, x, y))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_3D: out_of_bounds");
   }

   switch (evaluate_op_2d)
   {

      case MAGNITUDE_OP:
      {
         const Real u = evaluate (vector_element_u, k, x, y, VALUE);
         const Real v = evaluate (vector_element_v, k, x, y, VALUE);
         return sqrt (u*u + v*v);
         break;
      }

      case VORTICITY_OP:
      {
         const Real v_x = evaluate (vector_element_v, k, x, y, DX);
         const Real u_y = evaluate (vector_element_u, k, x, y, DY);
         return v_x - u_y;
         break;
      }

      case DIVERGENCE_OP:
      {
         const Real u_x = evaluate (vector_element_u, k, x, y, DX);
         const Real v_y = evaluate (vector_element_v, k, x, y, DY);
         return u_x + v_y;
         break;
      }

   }

   return GSL_NAN;

}

Real
Vector_Data_3D::evaluate_uv (const Integer vector_element_u,
                             const Integer vector_element_v,
                             const Real z,
                             const Integer i,
                             const Integer j,
                             const Evaluate_Op_2D evaluate_op_2d) const
{

   if (out_of_bounds (z, i, j))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_3D: out_of_bounds");
   }

   switch (evaluate_op_2d)
   {

      case MAGNITUDE_OP:
      {
         const Real u = evaluate (vector_element_u, z, i, j, VALUE);
         const Real v = evaluate (vector_element_v, z, i, j, VALUE);
         return sqrt (u*u + v*v);
         break;
      }

      case VORTICITY_OP:
      {
         const Real v_x = evaluate (vector_element_v, z, i, j, DX);
         const Real u_y = evaluate (vector_element_u, z, i, j, DY);
         return v_x - u_y;
         break;
      }

      case DIVERGENCE_OP:
      {
         const Real u_x = evaluate (vector_element_u, z, i, j, DX);
         const Real v_y = evaluate (vector_element_v, z, i, j, DY);
         return u_x + v_y;
         break;
      }

   }

   return GSL_NAN;

}

Real
Vector_Data_3D::evaluate_uv (const Integer vector_element_u,
                             const Integer vector_element_v,
                             const Real z, 
                             const Real x, 
                             const Real y,
                             const Evaluate_Op_2D evaluate_op_2d) const
{

   if (out_of_bounds (z, x, y))
   {
      return GSL_NAN;
      //throw Out_Of_Bounds_Exception ("Vector_Data_3D: out_of_bounds");
   }

   switch (evaluate_op_2d)
   {

      case MAGNITUDE_OP:
      {
         const Real u = evaluate (vector_element_u, z, x, y, VALUE);
         const Real v = evaluate (vector_element_v, z, x, y, VALUE);
         return sqrt (u*u + v*v);
         break;
      }

      case VORTICITY_OP:
      {
         const Real v_x = evaluate (vector_element_v, z, x, y, DX);
         const Real u_y = evaluate (vector_element_u, z, x, y, DY);
         return v_x - u_y;
         break;
      }

      case DIVERGENCE_OP:
      {
         const Real u_x = evaluate (vector_element_u, z, x, y, DX);
         const Real v_y = evaluate (vector_element_v, z, x, y, DY);
         return u_x + v_y;
         break;
      }

   }

   return GSL_NAN;

}

Real
Vector_Data_3D::evaluate_nocheck (const Integer vector_element,
                                  const Integer kk,
                                  const Integer ii,
                                  const Integer jj,
                                  const Evaluate_Op evaluate_op) const
{

   Integer k = kk;
   Integer i = ii;
   Integer j = jj;

   if (periodics[0]) { standardize_node (0, k); }
   if (periodics[1]) { standardize_node (1, i); }
   if (periodics[2]) { standardize_node (2, j); }

   const Real hz = spacings[0];
   const Real hx = spacings[1];
   const Real hy = spacings[2];
   const bool uniform_z = gsl_finite (hz);
   const bool uniform_x = gsl_finite (hx);
   const bool uniform_y = gsl_finite (hy);

   switch (evaluate_op)
   {

      case VALUE:
      {
         return get_datum (vector_element, k, i, j);
      }

      case DZ:
      {

         if (!uniform_z) { return GSL_NAN; }

         Integer ak = k - 1;
         Integer bk = k + 1;

         if (!periodics[0])
         {
            const Integer n = get_coordinate_tuple (0).size ();
            if ((k == 0) || k == n - 1) { return 0; }
         }
         else
         {
            standardize_node (0, ak);
            standardize_node (0, bk);
         }

         const Real a = get_datum (vector_element, ak, i, j);
         const Real b = get_datum (vector_element, bk, i, j);
         return (b - a) / (hz + hz);

      }

      case DX:
      {

         if (!uniform_x) { return GSL_NAN; }

         Integer ai = i - 1;
         Integer bi = i + 1;

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if ((i == 0) || i == n - 1) { return 0; }
         }
         else
         {
            standardize_node (1, ai);
            standardize_node (1, bi);
         }

         const Real a = get_datum (vector_element, k, ai, j);
         const Real b = get_datum (vector_element, k, bi, j);
         return (b - a) / (hx + hx);

      }

      case DY:
      {

         if (!uniform_y) { return GSL_NAN; }

         Integer aj = j - 1;
         Integer bj = j + 1;

         if (!periodics[2])
         {
            const Integer n = get_coordinate_tuple (2).size ();
            if ((j == 0) || j == n - 1) { return 0; }
         }
         else
         {
            standardize_node (2, aj);
            standardize_node (2, bj);
         }

         const Real a = get_datum (vector_element, k, i, aj);
         const Real b = get_datum (vector_element, k, i, bj);
         return (b - a) / (hy + hy);

      }

      case DZ2:
      {

         if (!uniform_z) { return GSL_NAN; }

         Integer ak = k - 1;
         Integer bk = k;
         Integer ck = k + 1;

         if (!periodics[0])
         {
            const Integer n = get_coordinate_tuple (0).size ();
            if ((k == 0) || k == n - 1) { return 0; }
         }
         else
         {
            standardize_node (0, ak);
            standardize_node (0, bk);
            standardize_node (0, ck);
         }

         const Real a = get_datum (vector_element, ak, i, j);
         const Real b = get_datum (vector_element, bk, i, j);
         const Real c = get_datum (vector_element, ck, i, j);
         return (a - b - b + c) / (hz * hz);

      }

      case DX2:
      {

         if (!uniform_x) { return GSL_NAN; }

         Integer ai = i - 1;
         Integer bi = i;
         Integer ci = i + 1;

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if ((i == 0) || i == n - 1) { return 0; }
         }
         else
         {
            standardize_node (1, ai);
            standardize_node (1, bi);
            standardize_node (1, ci);
         }

         const Real a = get_datum (vector_element, k, ai, j);
         const Real b = get_datum (vector_element, k, bi, j);
         const Real c = get_datum (vector_element, k, ci, j);
         return (a - b - b + c) / (hx * hx);

      }

      case DY2:
      {

         if (!uniform_y) { return GSL_NAN; }

         Integer aj = j - 1;
         Integer bj = j;
         Integer cj = j + 1;

         if (!periodics[2])
         {
            const Integer n = get_coordinate_tuple (2).size ();
            if ((j == 0) || j == n - 1) { return 0; }
         }
         else
         {
            standardize_node (2, aj);
            standardize_node (2, bj);
            standardize_node (2, cj);
         }

         const Real a = get_datum (vector_element, k, i, aj);
         const Real b = get_datum (vector_element, k, i, bj);
         const Real c = get_datum (vector_element, k, i, cj);
         return (a - b - b + c) / (hy * hy);

      }

      case DZX:
      {

         if (!(uniform_z && uniform_x)) { return GSL_NAN; }

         Integer ak = k - 1;
         Integer bk = k + 1;
         Integer ai = i - 1;
         Integer bi = i + 1;

         if (!periodics[0])
         {
            const Integer n = get_coordinate_tuple (0).size ();
            if ((k == 0) || k == n - 1) { return 0; }
         }
         else
         {
            standardize_node (0, ak);
            standardize_node (0, bk);
         }

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if ((k == 0) || k == n - 1) { return 0; }
         }
         else
         {
            standardize_node (1, ai);
            standardize_node (1, bi);
         }

         const Real a = get_datum (vector_element, ak, ai, j);
         const Real b = get_datum (vector_element, ak, bi, j);
         const Real c = get_datum (vector_element, bk, ai, j);
         const Real d = get_datum (vector_element, bk, bi, j);
         return (a - b - c + d) / (4 * hz * hx);

      }

      case DZY:
      {

         if (!(uniform_z && uniform_y)) { return GSL_NAN; }

         Integer ak = k - 1;
         Integer bk = k + 1;
         Integer aj = j - 1;
         Integer bj = j + 1;

         if (!periodics[0])
         {
            const Integer n = get_coordinate_tuple (0).size ();
            if ((k == 0) || k == n - 1) { return 0; }
         }
         else
         {
            standardize_node (0, ak);
            standardize_node (0, bk);
         }

         if (!periodics[2])
         {
            const Integer n = get_coordinate_tuple (2).size ();
            if ((j == 0) || j == n - 1) { return 0; }
         }
         else
         {
            standardize_node (2, aj);
            standardize_node (2, bj);
         }

         const Real a = get_datum (vector_element, ak, i, aj);
         const Real b = get_datum (vector_element, ak, i, bj);
         const Real c = get_datum (vector_element, bk, i, aj);
         const Real d = get_datum (vector_element, bk, i, bj);
         return (a - b - c + d) / (4 * hz * hy);

      }

      case DXY:
      {

         if (!(uniform_x && uniform_y)) { return GSL_NAN; }

         Integer ai = i - 1;
         Integer bi = i + 1;
         Integer aj = j - 1;
         Integer bj = j + 1;

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if ((i == 0) || i == n - 1) { return 0; }
         }
         else
         {
            standardize_node (1, ai);
            standardize_node (1, bi);
         }

         if (!periodics[2])
         {
            const Integer n = get_coordinate_tuple (2).size ();
            if ((j == 0) || j == n - 1) { return 0; }
         }
         else
         {
            standardize_node (2, aj);
            standardize_node (2, bj);
         }

         const Real a = get_datum (vector_element, k, ai, aj);
         const Real b = get_datum (vector_element, k, ai, bj);
         const Real c = get_datum (vector_element, k, bi, aj);
         const Real d = get_datum (vector_element, k, bi, bj);
         return (a - b - c + d) / (4 * hx * hy);

      }

      case DZXY:
      {

         if (!(uniform_z && uniform_x && uniform_y)) { return GSL_NAN; }

         Integer ak = k - 1;
         Integer bk = k + 1;
         Integer ai = i - 1;
         Integer bi = i + 1;
         Integer aj = j - 1;
         Integer bj = j + 1;

         if (!periodics[0])
         {
            const Integer n = get_coordinate_tuple (0).size ();
            if ((k == 0) || k == n - 1) { return 0; }
         }
         else
         {
            standardize_node (0, ak);
            standardize_node (0, bk);
         }

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if ((i == 0) || i == n - 1) { return 0; }
         }
         else
         {
            standardize_node (1, ai);
            standardize_node (1, bi);
         }

         if (!periodics[2])
         {
            const Integer n = get_coordinate_tuple (2).size ();
            if ((j == 0) || j == n - 1) { return 0; }
         }
         else
         {
            standardize_node (2, aj);
            standardize_node (2, bj);
         }

         const Real aaa = get_datum (vector_element, ak, ai, aj);
         const Real aab = get_datum (vector_element, ak, ai, bj);
         const Real aba = get_datum (vector_element, ak, bi, aj);
         const Real abb = get_datum (vector_element, ak, bi, bj);
         const Real baa = get_datum (vector_element, bk, ai, aj);
         const Real bab = get_datum (vector_element, bk, ai, bj);
         const Real bba = get_datum (vector_element, bk, bi, aj);
         const Real bbb = get_datum (vector_element, bk, bi, bj);
         return 0;
         //return (a - b - c + d) / (8 * hz * hx * hy);

      }

   }

}

Real
Vector_Data_3D::evaluate_nocheck (const Integer vector_element,
                                  const Integer k,
                                  const Real xx,
                                  const Real yy,
                                  const Evaluate_Op evaluate_op) const
{

   Real x = xx;
   Real y = yy;

   if (periodics[1]) { standardize_coordinate (1, x); }
   if (periodics[2]) { standardize_coordinate (2, y); }

   return evaluate_bilinear (vector_element, k, x, y, evaluate_op);

}

Real
Vector_Data_3D::evaluate_nocheck (const Integer vector_element,
                                  const Real zz,
                                  const Integer ii,
                                  const Integer jj,
                                  const Evaluate_Op evaluate_op) const
{

   Real z = zz;
   Integer i = ii;
   Integer j = jj;

   if (periodics[0]) { standardize_coordinate (0, z); }

   if (evaluate_op == VALUE)
   {
      return evaluate_linear_z (vector_element, z, i, j);
   }

   const Integer ve = vector_element;

   if (periodics[1]) { standardize_node (1, i); }
   if (periodics[2]) { standardize_node (2, j); }

   const Real h = spacings[1];
   const Real k = spacings[2];
   const bool uniform_x = gsl_finite (h);
   const bool uniform_y = gsl_finite (k);

   switch (evaluate_op)
   {

      default:
      {
         return 0;
      }

      case DZ:
      {

         const Integer k = get_node (0, z);
         const Real dz = get_spacing (0, k);

         const Real value_0 = evaluate_nocheck (ve, k, i, j);
         const Real value_1 = evaluate_nocheck (ve, k + 1, i, j);

         return (value_1 - value_0) / dz;

      };

      case DX:
      {

         if (!uniform_x) { return GSL_NAN; }

         Integer ai = i - 1;
         Integer bi = i + 1;

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if ((i == 0) || i == n - 1) { return 0; }
         }
         else
         {
            standardize_node (1, ai);
            standardize_node (1, bi);
         }

         const Real a = evaluate_linear_z (vector_element, z, ai, j);
         const Real b = evaluate_linear_z (vector_element, z, bi, j);
         return (b - a) / (h + h);

      }

      case DY:
      {

         if (!uniform_y) { return GSL_NAN; }

         Integer aj = j - 1;
         Integer bj = j + 1;

         if (!periodics[2])
         {
            const Integer n = get_coordinate_tuple (2).size ();
            if ((j == 0) || j == n - 1) { return 0; }
         }
         else
         {
            standardize_node (2, aj);
            standardize_node (2, bj);
         }

         const Real a = evaluate_linear_z (vector_element, z, i, aj);
         const Real b = evaluate_linear_z (vector_element, z, i, bj);
         return (b - a) / (k + k);

      }

      case DX2:
      {

         if (!uniform_x) { return GSL_NAN; }

         Integer ai = i - 1;
         Integer bi = i;
         Integer ci = i + 1;

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if ((i == 0) || i == n - 1) { return 0; }
         }
         else
         {
            standardize_node (1, ai);
            standardize_node (1, bi);
            standardize_node (1, ci);
         }

         const Real a = evaluate_linear_z (vector_element, z, ai, j);
         const Real b = evaluate_linear_z (vector_element, z, bi, j);
         const Real c = evaluate_linear_z (vector_element, z, ci, j);
         return (a - b - b + c) / (h * h);

      }

      case DY2:
      {

         if (!uniform_y) { return GSL_NAN; }

         Integer aj = j - 1;
         Integer bj = j;
         Integer cj = j + 1;

         if (!periodics[2])
         {
            const Integer n = get_coordinate_tuple (2).size ();
            if ((j == 0) || j == n - 1) { return 0; }
         }
         else
         {
            standardize_node (2, aj);
            standardize_node (2, bj);
            standardize_node (2, cj);
         }

         const Real a = evaluate_linear_z (vector_element, z, i, aj);
         const Real b = evaluate_linear_z (vector_element, z, i, bj);
         const Real c = evaluate_linear_z (vector_element, z, i, cj);
         return (a - b - b + c) / (h * h);

      }

      case DXY:
      {

         if (!(uniform_x && uniform_y)) { return GSL_NAN; }

         Integer ai = i - 1;
         Integer bi = i + 1;
         Integer aj = j - 1;
         Integer bj = j + 1;

         if (!periodics[1])
         {
            const Integer n = get_coordinate_tuple (1).size ();
            if ((i == 0) || i == n - 1) { return 0; }
         }
         else
         {
            standardize_node (1, ai);
            standardize_node (1, bi);
         }

         if (!periodics[2])
         {
            const Integer n = get_coordinate_tuple (2).size ();
            if ((j == 0) || j == n - 1) { return 0; }
         }
         else
         {
            standardize_node (2, aj);
            standardize_node (2, bj);
         }

         const Real a = evaluate_linear_z (vector_element, z, ai, aj);
         const Real b = evaluate_linear_z (vector_element, z, ai, bj);
         const Real c = evaluate_linear_z (vector_element, z, bi, aj);
         const Real d = evaluate_linear_z (vector_element, z, bi, bj);
         return (a - b - c + d) / (4 * h * k);

      }

   }
}

Real
Vector_Data_3D::evaluate_nocheck (const Integer vector_element,
                                  const Real zz,
                                  const Real xx,
                                  const Real yy,
                                  const Evaluate_Op evaluate_op) const
{

   Real z = zz;
   Real x = xx;
   Real y = yy;

   if (periodics[0]) { standardize_coordinate (0, z); }
   if (periodics[1]) { standardize_coordinate (1, x); }
   if (periodics[2]) { standardize_coordinate (2, y); }

   if (tricubic_coefficients_ptrs[vector_element] != NULL)
   {
      return evaluate_tricubic (vector_element, z, x, y, evaluate_op);
   }
   else
   {
      return evaluate_trilinear (vector_element, z, x, y, evaluate_op);
   }

}

void
Vector_Data_3D::subtract_z_mean (const Integer vector_element)
{

   const Size_3D& size_3d = get_size_3d ();

   for (Integer i = 0; i < size_3d.i; i++)
   {
      
      for (Integer j = 0; j < size_3d.j; j++)
      {

         Integer n = 0;
         Real sigma = 0;

         for (Integer k = 0; k < size_3d.k; k++)
         {
            const Real& datum = get_datum (vector_element, k, i, j);
            if (gsl_isnan (datum)) { continue; }
            sigma += datum;
            n++;
         }

         const Real mean = sigma / n;

         for (Integer k = 0; k < size_3d.k; k++)
         {
            Real& datum = get_datum (vector_element, k, i, j);
            datum -= mean;
         }

      }
   }
      
}

void
Vector_Data_3D::subtract_x_mean (const Integer vector_element)
{

   const Size_3D& size_3d = get_size_3d ();

   for (Integer k = 0; k < size_3d.k; k++)
   {
      
      for (Integer j = 0; j < size_3d.j; j++)
      {

         Integer n = 0;
         Real sigma = 0;

         for (Integer i = 0; i < size_3d.i; i++)
         {
            const Real& datum = get_datum (vector_element, k, i, j);
            if (gsl_isnan (datum)) { continue; }
            sigma += datum;
            n++;
         }

         const Real mean = sigma / n;

         for (Integer i = 0; i < size_3d.i; i++)
         {
            Real& datum = get_datum (vector_element, k, i, j);
            datum -= mean;
         }

      }
   }

}

void
Vector_Data_3D::subtract_y_mean (const Integer vector_element)
{

   const Size_3D& size_3d = get_size_3d ();

   for (Integer k = 0; k < size_3d.k; k++)
   {
      
      for (Integer i = 0; i < size_3d.i; i++)
      {

         Integer n = 0;
         Real sigma = 0;

         for (Integer j = 0; j < size_3d.j; j++)
         {
            const Real& datum = get_datum (vector_element, k, i, j);
            if (gsl_isnan (datum)) { continue; }
            sigma += datum;
            n++;
         }

         const Real mean = sigma / n;

         for (Integer j = 0; j < size_3d.j; j++)
         {
            Real& datum = get_datum (vector_element, k, i, j);
            datum -= mean;
         }

      }
   }

}

void
Vector_Data_3D::subtract_zx_mean (const Integer vector_element)
{

   const Size_3D& size_3d = get_size_3d ();

   for (Integer j = 0; j < size_3d.j; j++)
   {
      
      Integer n = 0;
      Real sigma = 0;

      for (Integer k = 0; k < size_3d.k; k++)
      {
         for (Integer i = 0; i < size_3d.i; i++)
         {
            const Real& datum = get_datum (vector_element, k, i, j);
            if (gsl_isnan (datum)) { continue; }
            sigma += datum;
            n++;
         }
      }

      const Real mean = sigma / n;

      for (Integer k = 0; k < size_3d.k; k++)
      {
         for (Integer i = 0; i < size_3d.i; i++)
         {
            Real& datum = get_datum (vector_element, k, i, j);
            datum -= mean;
         }
      }

   }

}

void
Vector_Data_3D::subtract_zy_mean (const Integer vector_element)
{

   const Size_3D& size_3d = get_size_3d ();

   for (Integer i = 0; i < size_3d.i; i++)
   {
      
      Integer n = 0;
      Real sigma = 0;

      for (Integer k = 0; k < size_3d.k; k++)
      {
         for (Integer j = 0; j < size_3d.j; j++)
         {
            const Real& datum = get_datum (vector_element, k, i, j);
            if (gsl_isnan (datum)) { continue; }
            sigma += datum;
            n++;
         }
      }

      const Real mean = sigma / n;

      for (Integer k = 0; k < size_3d.k; k++)
      {
         for (Integer j = 0; j < size_3d.j; j++)
         {
            Real& datum = get_datum (vector_element, k, i, j);
            datum -= mean;
         }
      }

   }

}

void
Vector_Data_3D::subtract_xy_mean (const Integer vector_element)
{

   const Size_3D& size_3d = get_size_3d ();

   for (Integer k = 0; k < size_3d.k; k++)
   {
      
      Integer n = 0;
      Real sigma = 0;

      for (Integer i = 0; i < size_3d.i; i++)
      {
         for (Integer j = 0; j < size_3d.j; j++)
         {
            const Real& datum = get_datum (vector_element, k, i, j);
            if (gsl_isnan (datum)) { continue; }
            sigma += datum;
            n++;
         }
      }

      const Real mean = sigma / n;

      for (Integer i = 0; i < size_3d.i; i++)
      {
         for (Integer j = 0; j < size_3d.j; j++)
         {
            Real& datum = get_datum (vector_element, k, i, j);
            datum -= mean;
         }
      }

   }

}

Scalar_Data_1D::Scalar_Data_1D (const Integer size_1d,
                                const Domain_1D& domain_1d,
                                const bool periodic)
   : Vector_Data_1D (1, size_1d, domain_1d)
{
}

Scalar_Data_1D::Scalar_Data_1D (const Tuple coordinate_tuple,
                                const bool periodic)
   : Vector_Data_1D (1, coordinate_tuple, periodic)
{
}

void
Scalar_Data_1D::set_datum (const Integer node,
                           const Real datum)
{
   Vector_Data_1D::set_datum (0, node, datum);
}

const Real&
Scalar_Data_1D::get_datum (const Integer node) const
{
   return Vector_Data_1D::get_datum (0, node);
}

Real&
Scalar_Data_1D::get_datum (const Integer node)
{
   return Vector_Data_1D::get_datum (0, node);
}

Real
Scalar_Data_1D::evaluate (const Real coordinate,
                          const Evaluate_Op evaluate_op) const
{
   return Vector_Data_1D::evaluate (0, coordinate, evaluate_op);
}

Scalar_Data_2D::Scalar_Data_2D (const Size_2D& size_2d,
                                const Domain_2D& domain_2d,
                                const bool periodic_x,
                                const bool periodic_y)
              : Vector_Data_2D (1, size_2d, domain_2d, periodic_x, periodic_y)
{
}

Scalar_Data_2D::Scalar_Data_2D (const Tuple coordinate_tuple_x,
                                const Tuple coordinate_tuple_y,
                                const bool periodic_x,
                                const bool periodic_y)
              : Vector_Data_2D (1, coordinate_tuple_x, coordinate_tuple_y,
                                periodic_x, periodic_y)
{
}

Domain_1D
Scalar_Data_2D::get_max_min () const
{
   return Vector_Data_2D::get_max_min (0);
}

void
Scalar_Data_2D::set_datum (const Integer node_x,
                           const Integer node_y,
                           const Real datum)
{
   Vector_Data_2D::set_datum (0, node_x, node_y, datum);
}

const Real&
Scalar_Data_2D::get_datum (const Integer node_x,
                           const Integer node_y) const
{
   return Vector_Data_2D::get_datum (0, node_x, node_y);
}

Real&
Scalar_Data_2D::get_datum (const Integer node_x,
                           const Integer node_y)
{
   return Vector_Data_2D::get_datum (0, node_x, node_y);
}

Real
Scalar_Data_2D::evaluate (const Real coordinate_x,
                          const Real coordinate_y,
                          const Evaluate_Op evaluate_op) const
{
   return Vector_Data_2D::evaluate (0, coordinate_x, coordinate_y, evaluate_op);
}


Scalar_Data_3D::Scalar_Data_3D (const Size_3D& size_3d,
                                const Domain_3D& domain_3d,
                                const bool periodic_z,
                                const bool periodic_x,
                                const bool periodic_y)
   : Vector_Data_3D (1, size_3d, domain_3d, periodic_z, periodic_x, periodic_y)
{
}

Scalar_Data_3D::Scalar_Data_3D (const Tuple coordinate_tuple_z,
                                const Tuple coordinate_tuple_x,
                                const Tuple coordinate_tuple_y,
                                const bool periodic_z,
                                const bool periodic_x,
                                const bool periodic_y)
   : Vector_Data_3D (1, coordinate_tuple_z, coordinate_tuple_x,
                     coordinate_tuple_y, periodic_z, periodic_x, periodic_y)
{
}

Scalar_Data_3D::Scalar_Data_3D (const Tuple coordinate_tuple_z,
                                const Size_2D& size_2d,
                                const Domain_2D& domain_2d,
                                const bool periodic_z,
                                const bool periodic_x,
                                const bool periodic_y)
   : Vector_Data_3D (1, coordinate_tuple_z, size_2d, domain_2d,
                     periodic_z, periodic_x, periodic_y)
{
}

void
Scalar_Data_3D::set_datum (const Integer node_z,
                           const Integer node_x,
                           const Integer node_y,
                           const Real datum)
{
   Vector_Data_3D::set_datum (0, node_z, node_x, node_y, datum);
}

const Real&
Scalar_Data_3D::get_datum (const Integer node_z,
                           const Integer node_x,
                           const Integer node_y) const
{
   return Vector_Data_3D::get_datum (0, node_z, node_x, node_y);
}

Real&
Scalar_Data_3D::get_datum (const Integer node_z,
                           const Integer node_x,
                           const Integer node_y)
{
   return Vector_Data_3D::get_datum (0, node_z, node_x, node_y);
}

Real
Scalar_Data_3D::evaluate (const Real coordinate_z,
                          const Real coordinate_x,
                          const Real coordinate_y,
                          const Evaluate_Op evaluate_op) const
{
   return Vector_Data_3D::evaluate (0, coordinate_z,
      coordinate_x, coordinate_y, evaluate_op);
}

Real
Scalar_Data_3D::evaluate (const Integer node_z,
                          const Real coordinate_x,
                          const Real coordinate_y,
                          const Evaluate_Op evaluate_op) const
{
   return Vector_Data_3D::evaluate (0, node_z,
      coordinate_x, coordinate_y, evaluate_op);
}

Real
Bezier_Curve::bernstein (const Real t,
                         const Integer k,
                         const Integer n) const
{                  
   const Integer nk = gsl_sf_choose (n, k);
   return gsl_sf_choose (n, k) * gsl_pow_int (t, k) * gsl_pow_int (1-t, n-k);
}

Bezier_Curve::Bezier_Curve (const vector<Tuple>& tuple_vector)
   : vector<Tuple> (tuple_vector),
     tuple_size (tuple_vector.front ().size ())
{
}

Real
Bezier_Curve::get_value_at (const Real t,
                            const Integer tuple_index) const
{

   Real value = 0;
   const Integer n = size ();

   for (Integer k = 0; k < n; k++)
   {
      const Tuple tuple = at (k);
      value += tuple[tuple_index] * bernstein (t, k, n - 1);
   }

   return value;

}

Tuple
Bezier_Curve::get_value_at (const Real t) const
{

   Tuple tuple;

   for (Integer tuple_index = 0; tuple_index < tuple_size; tuple_index++)
   {
      tuple.push_back (get_value_at (t, tuple_index));
   }

   return tuple;

}

Scalar
B_Spline::cox_deboor (const Real t,
                      const Integer k,
                      const Integer d) const
{

   if (d == 1)
   {
      if (t >= knot_tuple[k] && t < knot_tuple[k+1]) { return 1; }
      else { return 0; }
   }
   else
   {
      Real a = (t - knot_tuple[k]) / (knot_tuple[k+d-1] - knot_tuple[k]);
      Real b = (knot_tuple[k+d] - t) / (knot_tuple[k+d] - knot_tuple[k+1]);
      if (!gsl_finite (a)) { a = 0; }
      if (!gsl_finite (b)) { b = 0; }
      return a * cox_deboor (t, k, d-1) + b * cox_deboor (t, k+1, d-1);
   }

}

B_Spline::B_Spline (const vector<Tuple>& tuple_vector,
                    const Integer degree,
                    const bool open,
                    const bool repeat)
   : vector<Tuple> (tuple_vector),
                 n (tuple_vector.front ().size ()),
                 d (degree + 1)
{

   if (repeat)
   {
      for (Integer i = 0; i < degree; i++) { push_back (at (i)); }
   }

   knot_tuple.resize (size () + d);
   Real delta_t = 1 / Real (size () - d + 1);

   for (Integer i = 0; i < knot_tuple.size (); i++)
   {
      knot_tuple[i] = (i - d + 1) * delta_t;
      if (open)
      {
         if (i < d) { knot_tuple[i] = 0; }
         if (i > size ()) { knot_tuple[i] = 1; }
      }
   }


}

B_Spline::B_Spline (const vector<Tuple>& tuple_vector,
                    const Tuple& knot_tuple)
   : vector<Tuple> (tuple_vector),
        knot_tuple (knot_tuple),
                 n (tuple_vector.front ().size ()),
                 d (knot_tuple.size () - size ())
{
}

Real
B_Spline::get_value_at (const Real t,
                        const Integer tuple_index) const
{

   Real value = 0;

   for (Integer j = 0; j < size (); j++)
   {
      Real cb = cox_deboor (t, j, d);
      value += at (j).at (tuple_index) * cb;
   }

   return value;

}

Tuple
B_Spline::get_value_at (const Real t) const
{

   Tuple tuple;

   for (Integer tuple_index = 0; tuple_index < n; tuple_index++)
   {
      tuple.push_back (get_value_at (t, tuple_index));
   }

   return tuple;

}

void
Cartesian_Transform_2D::init (const Domain_1D& domain_x,
                              const Domain_1D& domain_y,
                              const Real width,
                              const Real height,
                              const Point_2D& origin,
                              const bool transpose)
{

   if (transpose)
   {
      translate (-domain_x.start, -domain_y.start);
      rotate (M_PI/2);
      scale (height / domain_x.get_span (), -width / domain_y.get_span ());
      translate (origin.x, origin.y);
   }
   else
   {
      translate (-domain_x.start, -domain_y.start);
      scale (width / domain_x.get_span (), height / domain_y.get_span ());
      translate (origin.x, origin.y);
   }

}


Cartesian_Transform_2D::Cartesian_Transform_2D (const Domain_1D& domain_x,
                                                const Domain_1D& domain_y,
                                                const Box_2D& box_2d,
                                                const bool transpose)
{
   const Real width = Real (box_2d.size_2d.i);
   const Real height = Real (box_2d.size_2d.j);
   const Point_2D origin (box_2d.index_2d.i, box_2d.index_2d.j);
   init (domain_x, domain_y, width, height, origin, transpose);
}

Cartesian_Transform_2D::Cartesian_Transform_2D (const Domain_1D& domain_x,
                                                const Domain_1D& domain_y,
                                                const Real width,
                                                const Real height,
                                                const Point_2D& origin,
                                                const bool transpose)
{
   init (domain_x, domain_y, width, height, origin, transpose);
}

Log_Transform_2D::Log_Transform_2D (const Domain_1D& domain_x,
                                    const Domain_1D& domain_y,
                                    const Real width,
                                    const Real height,
                                    const bool log_x,
                                    const bool log_y,
                                    const Point_2D& origin)
   : Cartesian_Transform_2D (
        log_x ?
           Domain_1D (log (domain_x.start), log (domain_x.end)) : domain_x,
        log_y ?
           Domain_1D (log (domain_y.start), log (domain_y.end)) : domain_y,
        width, height, origin),
     log_x (log_x),
     log_y (log_y)
{
}

bool     
Log_Transform_2D::out_of_domain (const Real x,
                                 const Real y) const
{
   return false;
}

void
Log_Transform_2D::transform (Real& transformed_x,
                             Real& transformed_y,
                             const Real x,
                             const Real y) const
{
   Cartesian_Transform_2D::transform (transformed_x, transformed_y,
      (log_x ? log (x) : x), (log_y ? log (y) : y));
}

void
Log_Transform_2D::reverse (Real& reversed_x,
                           Real& reversed_y,
                           const Real x,
                           const Real y) const
{
   Cartesian_Transform_2D::reverse (reversed_x, reversed_y, x, y);
   if (log_x) { reversed_x = exp (reversed_x); }
   if (log_y) { reversed_y = exp (reversed_y); }
}

void
Log_Transform_2D::transform_uv (Real& u,
                                Real& v,
                                const Real x,
                                const Real y) const
{
   if (log_x) { u = log ((x + u) / x); }
   if (log_y) { v = log ((y + v) / y); }
   Cartesian_Transform_2D::transform_uv (u, v, x, y);
}

Real
Log_Transform_2D::get_theta (const Real u,
                             const Real v,
                             const Real x,
                             const Real y) const
{
   Real uu = log_x ? log ((x + u) / x) : u;
   Real vv = log_y ? log ((y + v) / y) : v;
   return Cartesian_Transform_2D::get_theta (uu, vv, x, y);
}

Polar_Transform_2D::Polar_Transform_2D ()
{
   set (Point_2D (0, 0), 1);
}

Polar_Transform_2D::Polar_Transform_2D (const Point_2D& origin,
                                        const Real scale)
                              : origin (origin),
                                 scale (scale)
{
   set (origin, scale);
}

void
Polar_Transform_2D::set (const Point_2D& origin,
                         const Real scale)
{
   this->origin = origin;
   this->scale = scale;
}

bool 
Polar_Transform_2D::out_of_domain (const Real r,
                                   const Real theta) const
{
   return false;
}

void
Polar_Transform_2D::transform (Real& x,
                               Real& y,
                               const Real r,
                               const Real theta) const
{
   const Real scaled_r = scale * r;
   x = origin.x + scaled_r * cos (theta);
   y = origin.y + scaled_r * sin (theta);
}

void
Polar_Transform_2D::reverse (Real& r,
                             Real& theta,
                             const Real x,
                             const Real y) const
{
   const Real dx = (x - origin.x);
   const Real dy = (y - origin.y);
   r = sqrt (dx*dx + dy*dy) / scale;
   theta = atan2 (dy, dx);
}

Parabolic_Transform_2D::Parabolic_Transform_2D (const Point_2D& origin)
   : origin (origin)
{
}

bool
Parabolic_Transform_2D::out_of_domain (const Real u,
                                       const Real v) const
{
   return false;
}

void
Parabolic_Transform_2D::transform (Real& x,
                                   Real& y,
                                   const Real u,
                                   const Real v) const
{
   x = origin.x + (u*u - v*v) / 2;
   y = origin.y + u*v;
}

void
Parabolic_Transform_2D::reverse (Real& u,
                                 Real& v,
                                 const Real x,
                                 const Real y) const
{
   u = GSL_NAN;
   v = GSL_NAN;
}

Elliptic_Transform_2D::Elliptic_Transform_2D (const Point_2D& origin,
                                              const Real scale)
   : origin (origin),
     scale (scale)
{
}

bool
Elliptic_Transform_2D::out_of_domain (const Real u,
                                      const Real v) const
{
   return (u < 0);
}

void
Elliptic_Transform_2D::transform (Real& x,
                                  Real& y,
                                  const Real u,
                                  const Real v) const
{
   x = origin.x + scale * cosh (u) * cos (v);
   y = origin.y + scale * sinh (u) * sin (v);
}

void
Elliptic_Transform_2D::reverse (Real& u,
                                Real& v,
                                const Real x,
                                const Real y) const
{
   u = GSL_NAN;
   v = GSL_NAN;
}

Bipolar_Transform_2D::Bipolar_Transform_2D (const Point_2D& origin,
                                            const Real scale)
   : origin (origin),
     scale (scale)
{
}

bool
Bipolar_Transform_2D::out_of_domain (const Real u,
                                     const Real v) const
{
   return false;
}

void
Bipolar_Transform_2D::transform (Real& x,
                                 Real& y,
                                 const Real u,
                                 const Real v) const
{
   const Real denominator = cosh (v) - cos (u);
   x = origin.x + scale * sinh (v) / denominator;
   y = origin.y - scale * sin (u) / denominator;
}

void
Bipolar_Transform_2D::reverse (Real& u,
                               Real& v,
                               const Real x,
                               const Real y) const
{
   u = GSL_NAN;
   v = GSL_NAN;
}

void
Poisson_Disk::Process_Vector::push (const Point_2D& point)
{
   push_back (point);
}

Point_2D
Poisson_Disk::Process_Vector::pop ()
{
   const Integer i = Integer (round (random () * size () - 0.5));
   Process_Vector::iterator iterator = begin ();
   std::advance (iterator, i);
   const Point_2D p = at (i);
   erase (iterator);
   return p;
}

Poisson_Disk::Grid::Grid (const Size_2D& size_2d,
                          const Domain_2D& domain_2d,
                          const Real r,
                          const Real h)
   : Vector_Data_2D (2, size_2d, domain_2d),
     r (r),
     h (h)
{
   initialize (0, GSL_NAN);
   initialize (1, GSL_NAN);
}

void
Poisson_Disk::Grid::set_point (const Point_2D& point)
{
   const Integer i = get_nearest_node (0, point.x);
   const Integer j = get_nearest_node (1, point.y);
   set_datum (0, i, j, point.x);
   set_datum (1, i, j, point.y);
}

bool
Poisson_Disk::Grid::has_neighbour (const Point_2D& point)
{

   const Real r2 = r*r;
   const Integer i0 = get_nearest_node (0, point.x);
   const Integer j0 = get_nearest_node (1, point.y);

   const Integer max_i = size_nd.buffer[0] - 1;
   const Integer max_j = size_nd.buffer[1] - 1;

   for (Integer i = i0 - 2; i <= i0 + 2; i++)
   {
      if (i < 0 || i > max_i) { continue; }
      for (Integer j = j0 - 2; j <= j0 + 2; j++)
      {
         if (j < 0 || j > max_j) { continue; }
         const Real x = get_datum (0, i, j);
         const Real y = get_datum (1, i, j);
         if (gsl_isnan (x) || gsl_isnan (y)) { continue; }
         const Real dx = x - point.x;
         const Real dy = y - point.y;
         if (dx*dx + dy*dy < r2) { return true; }
      }
   }

   return false;

}

Point_2D
Poisson_Disk::get_random_nearby_point (const Point_2D& p,
                                       const Real r) const
{
  const Real radius = random (r, r+r);
  const Real theta = random (0, 2 * M_PI);
  const Real x = p.x + radius * cos (theta);
  const Real y = p.y + radius * sin (theta);
  return Point_2D (x, y);
}

Poisson_Disk::Poisson_Disk (const Domain_2D& domain_2d,
                            const Real r)
{

   const Real h = r / sqrt (2);
   Process_Vector process_vector;

   const Integer ni = Integer (ceil (domain_2d.get_span_x () / h));
   const Integer nj = Integer (ceil (domain_2d.get_span_y () / h));
   const Size_2D size_2d (ni, nj);
   Grid grid (size_2d, domain_2d, r, h);

   const Point_2D& p = domain_2d.get_random_point ();
   process_vector.push (p);
   push_back (p);
   grid.set_point (p);

   const Integer n = 30;

   while (process_vector.size () != 0)
   {

      const Point_2D& p = process_vector.pop ();

      for (Integer index = 0; index < n; index++)
      {

         const Point_2D& np = get_random_nearby_point (p, r);
         if (grid.out_of_bounds (np.x, np.y)) { continue; }
         if (grid.has_neighbour (np)) { continue; }

         process_vector.push (np);
         push_back (np);
         grid.set_point (np);

      }

   }

}

Real
Graph_Raster::f (const Real x,
                 const Real y) const
{
   const Real xx = x*x;
   const Real yy = y*y;
   const Real xx_p_yy = xx + yy;
   const Real xx_m_yy = xx - yy;
   const Real a = 100;
   return (x*x + y*y) * (x*x + y*y) - a*a * (x*x - y*y);
   //return x*x + y*y - a*a;
   //return y - 100 * sin (x / 10);
}

Graph_Raster::Graph_Raster (const Transform_2D& transform,
                            const Box_2D& box_2d,
                            const Color& color,
                            const Real width)
   : Raster (box_2d)
{

   const Size_2D& size_2d = box_2d.size_2d;
   const Index_2D& anchor = box_2d.index_2d;
   const Index_2D end_index (anchor.i + size_2d.i, anchor.j  + size_2d.j);
   Real x, y;

   #pragma omp parallel for
   for (Integer i = anchor.i; i < end_index.i; i++)
   {
      for (Integer j = anchor.j; j < end_index.j; j++)
      {
         transform.reverse (x, y, Real (i), Real (j));
         const Real value = f (x, y);
         const Color& c = (fabs (value) < width ? color : transparent);
         set_pixel (i - anchor.i, j - anchor.j, c);
      }
   }

}

