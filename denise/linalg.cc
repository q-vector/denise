//
// linalg.cc
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

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <denise/linalg.h>
#include <denise/exception.h>

using namespace std;
using namespace denise;

Vector::Vector (const Integer size,
                const bool set_zero)
{
   if (set_zero)
   {
      gv = gsl_vector_alloc (size);
   }
   else
   {
      gv = gsl_vector_calloc (size);
   }
}

Vector::Vector (const Vector& vector)
{
   this->gv = gsl_vector_alloc (vector.get_size ());
   gsl_vector_memcpy (this->gv, vector.gv);
}

Vector::Vector (const Vector& vector_a,
                const Vector& vector_b,
                const Vector_Op vector_op)
{

   const Integer a = vector_a.get_size ();
   const Integer b = vector_b.get_size ();
   if (a != b) { throw Exception ("size mismatch"); }
   const Integer n = a;

   this->gv = gsl_vector_alloc (n);
   gsl_vector_memcpy (this->gv, vector_a.gv);

   op (vector_b, vector_op);

}

Vector::Vector (gsl_vector* gv)
{
   this->gv = gv;
}

Vector::~Vector ()
{
   gsl_vector_free (gv);
}

Integer
Vector::get_size () const
{
   return gv->size;
}
 
gsl_vector*
Vector::get_gv () const
{
   return gv;
}

const Real
Vector::get_datum (const Integer i) const
{
   return gsl_vector_get (gv, i);
}

double&
Vector::get_datum (const Integer i)
{
   return *gsl_vector_ptr (gv, i);
}

void
Vector::set_datum (const Integer i,
                   const Real datum)
{
   gsl_vector_set (gv, i, datum);
}

void
Vector::op (const Vector& vector,
            const Vector_Op vector_op)
{

   const Integer n = get_size ();
   if (n != vector.get_size ()) { throw Exception ("size mismatch"); }

   switch (vector_op)
   {

      default:
         throw Exception ("vector_op not supported");
         break;

      case VO_ADD:
         gsl_vector_add (this->gv, vector.gv);
         break;

      case VO_SUBTRACT:
         gsl_vector_sub (this->gv, vector.gv);
         break;

      case VO_MULTIPLY:
         gsl_vector_mul (this->gv, vector.gv);
         break;

      case VO_DIVIDE:
         gsl_vector_div (this->gv, vector.gv);
         break;

   }

}

Real
Vector::dot (const Vector& vector) const
{
   double result;
   gsl_blas_ddot (gv, vector.gv, &result);
   return Real (result);
}

void
Matrix::calculate_inverse ()
{

   if (!is_square ()) { throw Exception ("not square matrix"); }
   const Integer n = get_rows ();

   this->inverse_ptr = new Matrix (n, n);

   int sign;
   Matrix lu (*this);

   gsl_permutation* p = gsl_permutation_alloc (n);
   gsl_linalg_LU_decomp (lu.get_gm (), p, &sign);
   gsl_linalg_LU_invert (lu.get_gm (), p, inverse_ptr->get_gm ());
   this->determinant = gsl_linalg_LU_det (lu.get_gm (), sign);

   gsl_permutation_free (p);

}

Matrix::Matrix (const Size_2D& size_2d,
                const bool set_zero)
   : inverse_ptr (NULL),
     determinant (GSL_NAN)
{
   if (set_zero)
   {
      gm = gsl_matrix_alloc (size_2d.i, size_2d.j);
   }
   else
   {
      gm = gsl_matrix_calloc (size_2d.i, size_2d.j);
   }
}

Matrix::Matrix (const Integer rows,
                const Integer columns,
                const bool set_zero)
   : inverse_ptr (NULL),
     determinant (GSL_NAN)
{
   if (set_zero)
   {
      gm = gsl_matrix_alloc (rows, columns);
   }
   else
   {
      gm = gsl_matrix_calloc (rows, columns);
   }
}

Matrix::Matrix (const Matrix& matrix)
   : determinant (matrix.determinant)
{

   const Integer rows = matrix.get_rows ();
   const Integer columns = matrix.get_columns ();
   this->gm = gsl_matrix_alloc (rows, columns);
   gsl_matrix_memcpy (this->gm, matrix.gm);

   if (matrix.inverse_ptr == NULL) { this->inverse_ptr = NULL; }
   else
   {
      this->inverse_ptr = new Matrix (rows, columns);
      gsl_matrix* source = this->inverse_ptr->get_gm ();
      gsl_matrix* destination = this->inverse_ptr->get_gm ();
      gsl_matrix_memcpy (destination, source);
   }

}

Matrix::Matrix (gsl_matrix* gm)
   : inverse_ptr (NULL),
     determinant (GSL_NAN)
{
   this->gm = gm;
}

Matrix::~Matrix ()
{
   gsl_matrix_free (gm);
   if (inverse_ptr != NULL) { delete inverse_ptr; }
}

bool
Matrix::is_square () const
{
   return (get_rows () == get_columns ());
}

Integer
Matrix::get_rows () const
{
   return gm->size1;
}

Integer
Matrix::get_columns () const
{
   return gm->size2;
}

const Matrix&
Matrix::get_inverse ()
{
   if (inverse_ptr == NULL) { calculate_inverse (); }
   return *inverse_ptr;
}

const Real&
Matrix::get_determinant ()
{
   if (gsl_isnan (determinant)) { calculate_inverse (); }
   return determinant;
}

gsl_matrix*
Matrix::get_gm () const
{
   return gm;
}

const Real
Matrix::get_datum (const Integer i,
                   const Integer j) const
{
   return gsl_matrix_get (gm, i, j);
}

double&
Matrix::get_datum (const Integer i,
                   const Integer j)
{
   return *gsl_matrix_ptr (gm, i, j);
}

void
Matrix::set_datum (const Integer i,
                   const Integer j,
                   const Real datum)
{
   gsl_matrix_set (gm, i, j, datum);
}

void
Matrix::fill_row (const Integer row,
                  const Real* row_data)
{
   for (Integer j = 0; j < gm->size2; j++)
   {
      gsl_matrix_set (gm, row, j, row_data[j]);
   }
}

void
Matrix::fill_column (const Integer column,
                     const Real* column_data)
{
   for (Integer i = 0; i < gm->size1; i++)
   {
      gsl_matrix_set (gm, i, column, column_data[i]);
   }
}

void
Matrix::scale (const Real s)
{
   gsl_matrix_scale (gm, s);
}

Real_Symmetric_Eigen::Real_Symmetric_Eigen (const Matrix& matrix,
                                            const bool sort)
{

   if (!matrix.is_square ()) { throw Exception ("not square matrix"); }
   const Integer n = matrix.get_rows ();

   gsl_vector* eval = gsl_vector_alloc (n);
   gsl_matrix* evec = gsl_matrix_alloc (n, n);
   gsl_matrix* gm = gsl_matrix_alloc (n, n);
   gsl_matrix_memcpy (gm, matrix.get_gm ());

   gsl_eigen_symmv_workspace* w = gsl_eigen_symmv_alloc (n);
   gsl_eigen_symmv (gm, eval, evec, w);
   gsl_eigen_symmv_free (w);
   gsl_matrix_free (gm);

   if (sort) { gsl_eigen_symmv_sort (eval, evec, GSL_EIGEN_SORT_ABS_DESC); }

   eigenvalue_vector_ptr = new Vector (eval);
   eigenvector_matrix_ptr = new Matrix (evec);

}

Real_Symmetric_Eigen::~Real_Symmetric_Eigen ()
{
   delete eigenvalue_vector_ptr;
   delete eigenvector_matrix_ptr;
}

const Vector&
Real_Symmetric_Eigen::get_eigenvalue_vector () const
{
   return *eigenvalue_vector_ptr;
}

const Matrix&
Real_Symmetric_Eigen::get_eigenvector_matrix () const
{
   return *eigenvector_matrix_ptr;
}

