//
// dataset.cc
// 
// Copyright (C) 2011 Simon E. Ching
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
#include <gsl/gsl_statistics_double.h>
#include "dataset.h"

using namespace std;
using namespace denise;

Dataset::Dataset (const Integer dimension,
                  const Integer number_of_points)
   : data_matrix_ptr (new Matrix (number_of_points, dimension)),
     covariance_matrix_ptr (NULL),
     mean_vector_ptr (NULL),
     residue_matrix_ptr (NULL),
     principal_component_ptr (NULL)
{
}

Dataset::Dataset (Matrix* data_matrix_ptr)
   : data_matrix_ptr (data_matrix_ptr),
     covariance_matrix_ptr (NULL),
     mean_vector_ptr (NULL),
     residue_matrix_ptr (NULL),
     principal_component_ptr (NULL)
{
}

Dataset::~Dataset ()
{
   delete data_matrix_ptr;
   delete mean_vector_ptr;
   delete residue_matrix_ptr;
   delete covariance_matrix_ptr;
}

void
Dataset::set_data (const Integer i,
                   const Real* data)
{
   data_matrix_ptr->fill_row (i, data);
}

void
Dataset::set_datum (const Integer i,
                    const Integer j,
                    const Real datum)
{
   data_matrix_ptr->set_datum (i, j, datum);
}

Integer
Dataset::get_number_of_points () const
{
   return data_matrix_ptr->get_rows ();
}

Integer
Dataset::get_dimension () const
{
   return data_matrix_ptr->get_columns ();
}

Matrix&
Dataset::get_data_matrix ()
{
   return *data_matrix_ptr;
}

const Matrix&
Dataset::get_data_matrix () const
{
   return *data_matrix_ptr;
}

Matrix&
Dataset::get_covariance_matrix ()
{
   compute_covariance (true);
   return *covariance_matrix_ptr;
}

const Vector&
Dataset::get_mean_vector ()
{
   compute_mean (true);
   return *mean_vector_ptr;
}

const Matrix&
Dataset::get_residue_matrix ()
{
   compute_residue (true);
   return *residue_matrix_ptr;
}

const Vector&
Dataset::get_pca_coefficient_vector ()
{
   compute_principal_component (true);
   return principal_component_ptr->get_eigenvalue_vector ();
}

const Matrix&
Dataset::get_pca_vector_matrix ()
{
   compute_principal_component (true);
   return principal_component_ptr->get_eigenvector_matrix ();
}

const Real
Dataset::get_mahalonobis (const Vector& vector)
{

   const Integer size = get_dimension ();

   const Vector& mju = get_mean_vector ();
   Matrix& sigma = get_covariance_matrix ();
   const Matrix& inverse_sigma = sigma.get_inverse ();
   const Real& det_sigma = sigma.get_determinant ();

   Vector inverse_sigma_deviation (size);
   Vector deviation (vector, mju, VO_SUBTRACT);

   gsl_blas_dgemv (CblasNoTrans, 1, inverse_sigma.get_gm (),
      deviation.get_gv (), 0, inverse_sigma_deviation.get_gv ());

   return deviation.dot (inverse_sigma_deviation);

}

void
Dataset::compute_residue (const bool idle_if_present)
{

   if (residue_matrix_ptr != NULL && idle_if_present) { return; }
   if (residue_matrix_ptr != NULL) { delete residue_matrix_ptr; }

   compute_mean (idle_if_present);

   const Integer rows = get_number_of_points ();
   const Integer k = get_dimension ();
   const Vector& mean_vector = get_mean_vector ();

   residue_matrix_ptr = new Matrix (rows, k);
   for (Integer j = 0; j < k; j++)
   {

      const Real mean = mean_vector.get_datum (j);

      for (Integer i = 0; i < rows; i++)
      {
         const Real datum = data_matrix_ptr->get_datum (i, j);
         const Real residue = datum - mean;
         residue_matrix_ptr->set_datum (i, j, residue); 
      }

   }

}

void
Dataset::compute_covariance (const bool idle_if_present)
{

   if (idle_if_present && covariance_matrix_ptr != NULL) { return; }
   if (covariance_matrix_ptr != NULL) { delete covariance_matrix_ptr; }
   compute_residue (idle_if_present);

   const Integer rows = get_number_of_points ();
   const Integer columns = get_dimension ();
   const Real one_over_n_minus_1 = 1.0 / Real (rows - 1);

   gsl_matrix* residue_gm = residue_matrix_ptr->get_gm ();
   gsl_matrix* covariance_gm = gsl_matrix_alloc (columns, columns);
   gsl_blas_dgemm (CblasTrans, CblasNoTrans, one_over_n_minus_1,
      residue_gm, residue_gm, 0, covariance_gm);

   covariance_matrix_ptr = new Matrix (covariance_gm);

}

void
Dataset::compute_mean (const bool idle_if_present)
{

   if (idle_if_present && mean_vector_ptr != NULL) { return; }
   if (mean_vector_ptr != NULL) { delete mean_vector_ptr; }

   const Integer rows = get_number_of_points ();
   const Integer k = get_dimension ();

   mean_vector_ptr = new Vector (k);
   for (Integer j = 0; j < k; j++)
   {

      Real sigma = 0;

      for (Integer i = 0; i < rows; i++)
      {
         sigma += data_matrix_ptr->get_datum (i, j);
      }

      mean_vector_ptr->set_datum (j, sigma / rows);

   }

}

void
Dataset::compute_principal_component (const bool sort,
                                      const bool idle_if_present)
{

   if (idle_if_present && principal_component_ptr != NULL) { return; }
   if (principal_component_ptr != NULL) { delete principal_component_ptr; }

   compute_covariance (idle_if_present);
   principal_component_ptr = new Real_Symmetric_Eigen (
      *covariance_matrix_ptr, sort);

}

int
Dataset_1D::compare (const void* a,
                     const void* b)
{
   return *(double*)a > *(double*)b;
}

Dataset_1D::Dataset_1D (const Integer number_of_points)
   : Dataset (1, number_of_points)
{
}

void
Dataset_1D::set_datum (const Integer i,
                       const Real datum)
{
   Dataset::set_datum (i, 0, datum);
}

void
Dataset_1D::sort ()
{
   double* data = get_data_matrix ().get_gm ()->data;
   const Integer n = get_number_of_points ();
   qsort (data, n, sizeof (double), (Dataset_1D::compare));
}

Real
Dataset_1D::get_mean ()
{
   return get_mean_vector ().get_datum (0);
}

Real
Dataset_1D::get_variance ()
{
   return get_covariance_matrix ().get_datum (0, 0);
}

Real
Dataset_1D::get_mean_square () const
{
   const Integer n = get_number_of_points ();
   double* squares = new double[n];
   double* data = get_data_matrix ().get_gm ()->data;
   for (Integer i = 0; i < n; i++) { squares[i] = gsl_pow_2 (data[i]); }
   const Real mean_square = gsl_stats_mean (squares, 1, n);
   delete[] squares;
   return mean_square;
}
   
Real
Dataset_1D::get_root_mean_square () const
{
   return get_rms ();
}
   
Real
Dataset_1D::get_rms () const
{
   return sqrt (get_mean_square ());
}

Real      
Dataset_1D::get_sd ()
{
   return get_standard_deviation ();
}
         
Real
Dataset_1D::get_standard_deviation ()
{
   return sqrt (get_variance ());
}
         
Real
Dataset_1D::get_max () const
{
   const Integer n = get_number_of_points ();
   const double* data = get_data_matrix ().get_gm ()->data;
   return gsl_stats_max (data, 1, n); 
}
         
Real
Dataset_1D::get_min () const
{
   const Integer n = get_number_of_points ();
   const double* data = get_data_matrix ().get_gm ()->data;
   return gsl_stats_min (data, 1, n); 
}

Real
Dataset_1D::get_median () const
{
   const Integer n = get_number_of_points ();
   const double* data = get_data_matrix ().get_gm ()->data;
   return gsl_stats_median_from_sorted_data (data, 1, n); 
}

Real
Dataset_1D::get_first_quartile () const
{
   return get_percentile (25);
}

Real
Dataset_1D::get_second_quartile () const
{
   return get_median ();
}

Real
Dataset_1D::get_third_quartile () const
{
   return get_percentile (75);
}

Real
Dataset_1D::get_percentile (const Real percentage) const
{
   const Integer n = get_number_of_points ();
   const double* data = get_data_matrix ().get_gm ()->data;
   return gsl_stats_quantile_from_sorted_data (data, 1, n, percentage*1e-2);
}

Dataset_2D::Dataset_2D (const Integer number_of_points)
   : Dataset (2, number_of_points)
{
}

void
Dataset_2D::set_data (const Integer i,
                      const Real datum_x,
                      const Real datum_y)
{
   Dataset::set_datum (i, 0, datum_x);
   Dataset::set_datum (i, 1, datum_y);
}

Real
Dataset_2D::get_mean_x ()
{
   return get_mean_vector ().get_datum (0);
}

Real
Dataset_2D::get_mean_y ()
{
   return get_mean_vector ().get_datum (1);
}

Real
Dataset_2D::get_variance_x ()
{
   return get_covariance_matrix ().get_datum (0, 0);
}

Real
Dataset_2D::get_variance_y ()
{
   return get_covariance_matrix ().get_datum (1, 1);
}

Real
Dataset_2D::get_covariance ()
{
   return get_covariance_matrix ().get_datum (0, 1);
}

