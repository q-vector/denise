//
// linalg.h
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

#ifndef DENISE_LINALG_H
#define DENISE_LINALG_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <denise/basics.h>

using namespace std;

namespace denise
{

   enum Vector_Op
   {
      VO_ADD,
      VO_SUBTRACT,
      VO_MULTIPLY,
      VO_DIVIDE,
      VO_DOT
   };

   class Vector
   {

      protected:

         gsl_vector*
         gv;

      public:

         Vector (const Integer size,
                 const bool set_zero = false);

         Vector (const Vector& vector);

         Vector (const Vector& vector_a,
                 const Vector& vector_b,
                 const Vector_Op vector_op);

         Vector (gsl_vector* gv);

         ~Vector ();

         Integer
         get_size () const;

         gsl_vector*
         get_gv () const;

         const Real
         get_datum (const Integer i) const;

         double&
         get_datum (const Integer i);

         void
         set_datum (const Integer i,
                    const Real datum);

         void
         op (const Vector& vector,
             const Vector_Op vector_op);

         Real
         dot (const Vector& vector) const;

   };

   class Matrix
   {

      protected:

         gsl_matrix*
         gm;

         Real
         determinant;

         Matrix*
         inverse_ptr;

         void
         calculate_inverse ();

      public:

         Matrix (const Size_2D& size_2d,
                 const bool set_zero = false);

         Matrix (const Integer rows,
                 const Integer columns,
                 const bool set_zero = false);

         Matrix (const Matrix& matrix);

         Matrix (gsl_matrix* gm);

         ~Matrix ();

         bool
         is_square () const;

         Integer
         get_rows () const;

         Integer
         get_columns () const;

         const Matrix&
         get_inverse ();

         const Real&
         get_determinant ();

         gsl_matrix*
         get_gm () const;

         const Real
         get_datum (const Integer i,
                    const Integer j) const;

         double&
         get_datum (const Integer i,
                    const Integer j);

         void
         set_datum (const Integer i,
                    const Integer j,
                    const Real datum);

         void
         fill_row (const Integer row,
                   const Real* row_data);

         void
         fill_column (const Integer column,
                      const Real* column_data);

         void
         scale (const Real s);

   };

   class Real_Symmetric_Eigen
   {

      private:

         Vector*
         eigenvalue_vector_ptr;

         Matrix*
         eigenvector_matrix_ptr;

      public:

         Real_Symmetric_Eigen (const Matrix& matrix,
                               const bool sort = true);

         ~Real_Symmetric_Eigen ();

         const Vector&
         get_eigenvalue_vector () const;

         const Matrix&
         get_eigenvector_matrix () const;

   };

};

#endif /* DENISE_LINALG_H */
