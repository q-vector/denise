//
// stat.cc
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

#include "stat.h"

using namespace std;
using namespace denise;

int
Treatment::compare (const void* a,
                    const void* b)
{
   return *(double*)a > *(double*)b;
}

Treatment::Treatment (const double* data,
                      const Integer size)
   : size (size)
{

   this->data = new double[size];

   for (Integer i = 0; i < size; i++)
   {
      this->data[i] = data[i];
   }

}

Treatment::Treatment (const Tuple& tuple)
   : size (tuple.size ())
{

   data = new double[size];

   for (Integer i = 0; i < size; i++)
   {
      data[i] = tuple[i];
   }

}

Treatment::~Treatment ()
{
   delete[] data;
}

void
Treatment::sort ()
{
   qsort (data, size, sizeof (double), (Treatment::compare));
}

const double*
Treatment::get_data () const
{
   return data;
}

Integer
Treatment::get_size () const
{
   return size;
}

Real
Treatment::get_mean () const
{
   return gsl_stats_mean (data, 1, size);
}

Real
Treatment::get_mean_square () const
{
   double* squares = new double[size];
   for (Integer i = 0; i < size; i++) { squares[i] = gsl_pow_2 (data[i]); }
   const Real mean_square = gsl_stats_mean (squares, 1, size);
   delete[] squares;
   return mean_square;
}

Real
Treatment::get_root_mean_square () const
{
   return get_rms ();
}

Real
Treatment::get_rms () const
{
   return sqrt (get_mean_square ());
}


Real
Treatment::get_variance () const
{
   return gsl_stats_variance (data, 1, size);
}

Real
Treatment::get_sd () const
{
   return get_standard_deviation ();
}

Real
Treatment::get_standard_deviation () const
{
   return sqrt (get_variance ());
}

Real 
Treatment::get_max () const
{
   return gsl_stats_max (data, 1, size);
}

Real
Treatment::get_min () const
{
   return gsl_stats_min (data, 1, size);
}

Real 
Treatment::get_median () const
{
   return gsl_stats_median_from_sorted_data (data, 1, size);
}

Real
Treatment::get_first_quartile () const
{
   return get_percentile (25);
}

Real
Treatment::get_second_quartile () const
{
   return get_median ();
}

Real
Treatment::get_third_quartile () const
{
   return get_percentile (75);
}

Real
Treatment::get_percentile (const Real percentage) const
{
   return gsl_stats_quantile_from_sorted_data (data, 1, size, percentage*1e-2);
}

void
Lognormal_Distribution::init (const Real mean,
                              const Real variance)
{
   this->mean = mean;
   this->variance = variance;
}

void
Lognormal_Distribution::init (const double* array,
                              const Integer n)
{

   const Real mean = gsl_stats_mean (array, 1, n);
   const Real variance = gsl_stats_variance_m (array, 1, n, mean);

   init (mean, variance);

}

Lognormal_Distribution::Lognormal_Distribution ()
{
   init (0, 1);
}

Lognormal_Distribution::Lognormal_Distribution (const Real mean,
                                                const Real variance)
{
   init (mean, variance);
}

Lognormal_Distribution::Lognormal_Distribution (const double* array,
                                                const Integer n)
{
   init (array, n);
}

Lognormal_Distribution::Lognormal_Distribution (const Tuple& tuple)
{

   Integer n = tuple.size ();
   double* array = new double[n];

   for (Integer i = 0; i < n; i++)
   {
      array[i] = tuple[i];
   }

   init (array, n);
   delete[] array;

}

Real
Lognormal_Distribution::get_pdf (const Real point) const
{
   return get_pdf (point, mean, variance);
}

Real
Lognormal_Distribution::get_P (const Real point) const
{
   return get_P (point, mean, variance);
}

Real
Lognormal_Distribution::get_Q (const Real point) const
{
   return get_Q (point, mean, variance);
}

Real
Lognormal_Distribution::get_P_inv (const Real P) const
{
   return get_P_inv (P, mean, variance);
}

Real
Lognormal_Distribution::get_Q_inv (const Real Q) const
{
   return get_Q_inv (Q, mean, variance);
}

Real
Lognormal_Distribution::get_mean () const
{
   return mean;
}

Real
Lognormal_Distribution::get_variance () const
{
   return variance;
}

Real
Lognormal_Distribution::get_sigma () const
{
   return sqrt (variance);
}

Real
Lognormal_Distribution::get_pdf (const Real point,
                                 const Real mean,
                                 const Real variance)
{
   return gsl_ran_lognormal_pdf (point, mean, sqrt (variance));
}

Real
Lognormal_Distribution::get_P (const Real point,
                               const Real mean,
                               const Real variance)
{
   return gsl_cdf_lognormal_P (point, mean, sqrt (variance));
}

Real
Lognormal_Distribution::get_Q (const Real point,
                               const Real mean,
                               const Real variance)
{
   return gsl_cdf_lognormal_Q (point, mean, sqrt (variance));
}

Real
Lognormal_Distribution::get_P_inv (const Real P,
                                   const Real mean,
                                   const Real variance)
{
  return gsl_cdf_lognormal_Pinv (P, mean, sqrt (variance));
}

Real
Lognormal_Distribution::get_Q_inv (const Real Q,
                                   const Real mean,
                                   const Real variance)
{
  return gsl_cdf_lognormal_Qinv (Q, mean, sqrt (variance));
}

Rayleigh_Distribution::Rayleigh_Distribution (const Real sigma)
                                     : sigma (sigma)
{
}

Real
Rayleigh_Distribution::get_pdf (const Real point) const
{
   return get_pdf (point, sigma);
}

Real
Rayleigh_Distribution::get_P (const Real point) const
{
   return get_P (point, sigma);
}

Real
Rayleigh_Distribution::get_Q (const Real point) const
{
   return get_Q (point, sigma);
}

Real
Rayleigh_Distribution::get_P_inv (const Real P) const
{
   return get_P_inv (P, sigma);
}

Real
Rayleigh_Distribution::get_Q_inv (const Real Q) const
{
   return get_Q_inv (Q, sigma);
}

Real
Rayleigh_Distribution::get_sigma () const
{
   return this->sigma;
}

Real
Rayleigh_Distribution::get_pdf (const Real point,
                                const Real sigma)
{
   return gsl_ran_rayleigh_pdf (point, sigma);
}

Real
Rayleigh_Distribution::get_P (const Real point,
                              const Real sigma)
{
   return gsl_cdf_rayleigh_P (point, sigma);
}

Real
Rayleigh_Distribution::get_Q (const Real point,
                              const Real sigma)
{
   return gsl_cdf_rayleigh_Q (point, sigma);
}

Real
Rayleigh_Distribution::get_P_inv (const Real P,
                                  const Real sigma)
{
   return gsl_cdf_rayleigh_Pinv (P, sigma);
}

Real
Rayleigh_Distribution::get_Q_inv (const Real Q,
                                  const Real sigma)
{
   return gsl_cdf_rayleigh_Qinv (Q, sigma);
}

void
Gaussian_Distribution::init (const Real mean,
                             const Real variance)
{
   this->mean = mean;
   this->variance = variance;
}

void
Gaussian_Distribution::init (const double* array,
                             const Integer n)
{

   const Real mean = gsl_stats_mean (array, 1, n);
   const Real variance = gsl_stats_variance_m (array, 1, n, mean);

   init (mean, variance);

}

Gaussian_Distribution::Gaussian_Distribution ()
{
   init (0, 1);
}

Gaussian_Distribution::Gaussian_Distribution (const Real mean,
                                              const Real variance)
{
   init (mean, variance);
}

Gaussian_Distribution::Gaussian_Distribution (Dataset_1D& dataset_1d)
{
   const Real mean = dataset_1d.get_mean ();
   const Real variance = dataset_1d.get_variance ();
   init (mean, variance);
}

Gaussian_Distribution::Gaussian_Distribution (const double* array,
                                              const Integer n)
{
   init (array, n);
}

Gaussian_Distribution::Gaussian_Distribution (const Tuple& tuple)
{

   const Integer n = tuple.size ();
   double* array = new double[n];

   for (Integer i = 0; i < n; i++)
   {
      array[i] = tuple[i];
   }

   delete[] array;
   init (array, n);

}

Real
Gaussian_Distribution::get_pdf (const Real point) const
{
   return get_pdf (point, mean, variance);
}

Real
Gaussian_Distribution::get_P (const Real point) const
{
   return get_P (point, mean, variance);
}

Real
Gaussian_Distribution::get_Q (const Real point) const
{
   return get_Q (point, mean, variance);
}

Real
Gaussian_Distribution::get_P_inv (const Real P) const
{
   return get_P_inv (P, mean, variance);
}

Real
Gaussian_Distribution::get_Q_inv (const Real Q) const
{
   return get_Q_inv (Q, mean, variance);
}

Real
Gaussian_Distribution::get_mean () const
{
   return mean;
}
  
Real
Gaussian_Distribution::get_variance () const
{
   return variance;
}

Real
Gaussian_Distribution::get_sigma () const
{
   return sqrt (variance);
}

Real
Gaussian_Distribution::get_pdf (const Real point,
                                const Real mean,
                                const Real variance)
{
   return gsl_ran_gaussian_pdf (point - mean, sqrt (variance));
}

Real
Gaussian_Distribution::get_P (const Real point,
                              const Real mean,
                              const Real variance)
{
   return gsl_cdf_gaussian_P (point - mean, sqrt (variance));
}

Real
Gaussian_Distribution::get_Q (const Real point,
                              const Real mean,
                              const Real variance)
{
   return gsl_cdf_gaussian_Q (point - mean, sqrt (variance));
}

Real
Gaussian_Distribution::get_P_inv (const Real P,
                                  const Real mean,
                                  const Real variance)
{
  return gsl_cdf_gaussian_Pinv (P, sqrt (variance)) + mean;
}

Real
Gaussian_Distribution::get_Q_inv (const Real Q,
                                  const Real mean,
                                  const Real variance)
{
  return gsl_cdf_gaussian_Qinv (Q, sqrt (variance)) + mean;
}

void
Bivariate_Gaussian_Distribution::init (const Real mean_x,
                                       const Real mean_y,
                                       const Real variance_x,
                                       const Real variance_y,
                                       const Real covariance)
{
   this->mean_x = mean_x;
   this->mean_y = mean_y;
   this->variance_x = variance_x;
   this->variance_y = variance_y;
   this->covariance = covariance;
}

void
Bivariate_Gaussian_Distribution::init (const double* x_array,
                                       const double* y_array,           
                                       const Integer n)
{

   Real mean_x = gsl_stats_mean (x_array, 1, n);
   Real mean_y = gsl_stats_mean (y_array, 1, n);
   Real variance_x = gsl_stats_variance_m (x_array, 1, n, mean_x);
   Real variance_y = gsl_stats_variance_m (y_array, 1, n, mean_y);
   Real covariance = gsl_stats_covariance_m (x_array, 1, y_array, 1, n,
                                               mean_x, mean_y);

   init (mean_x, mean_y, variance_x, variance_y, covariance);

}

Bivariate_Gaussian_Distribution::Bivariate_Gaussian_Distribution ()
{
   init (0, 0, 1, 1, 0);
}

Bivariate_Gaussian_Distribution::Bivariate_Gaussian_Distribution (const Real mean_x,
                                                                  const Real mean_y,
                                                                  const Real variance_x,
                                                                  const Real variance_y,
                                                                  const Real covariance)
{
   init (mean_x, mean_y, variance_x, variance_y, covariance);
}

Bivariate_Gaussian_Distribution::Bivariate_Gaussian_Distribution (Dataset_2D& dataset_2d)
{
   const Real mean_x = dataset_2d.get_mean_x ();
   const Real mean_y = dataset_2d.get_mean_y ();
   const Real variance_x = dataset_2d.get_variance_x ();
   const Real variance_y = dataset_2d.get_variance_y ();
   const Real covariance = dataset_2d.get_covariance ();
   init (mean_x, mean_y, variance_x, variance_y, covariance);
}

Bivariate_Gaussian_Distribution::Bivariate_Gaussian_Distribution (const double* x_array,
                                                                  const double* y_array,
                                                                  const Integer n)
{
   init (x_array, y_array, n);
}

Bivariate_Gaussian_Distribution::Bivariate_Gaussian_Distribution (const vector<Point_2D>& point_vector)
{

   const Integer n = point_vector.size ();

   double* x_array = new double[n];
   double* y_array = new double[n];

   for (Integer i = 0; i < n; i++)
   {
      x_array[i] = point_vector[i].x;
      y_array[i] = point_vector[i].y;
   }

   init (x_array, y_array, n);
   delete[] x_array;
   delete[] y_array;

}

Bivariate_Gaussian_Distribution::Bivariate_Gaussian_Distribution (const Duple_Vector& duple_vector)
{

   const Integer n = duple_vector.size ();

   double* x_array = new double[n];
   double* y_array = new double[n];

   for (Integer i = 0; i < n; i++)
   {
      x_array[i] = duple_vector[i].first;
      y_array[i] = duple_vector[i].second;
   }

   init (x_array, y_array, n);
   delete[] x_array;
   delete[] y_array;

}

Real
Bivariate_Gaussian_Distribution::get_pdf (const Point_2D& point_2d) const
{
   return get_pdf (point_2d.x, point_2d.y);
}

Real
Bivariate_Gaussian_Distribution::get_pdf (const Real x,
                                          const Real y) const
{
   return get_pdf (x, y, mean_x, mean_y, variance_x, variance_y, covariance);
}

Point_2D
Bivariate_Gaussian_Distribution::get_mean () const
{
   return Point_2D (mean_x, mean_y);
}

Real
Bivariate_Gaussian_Distribution::get_mean_x () const
{
   return mean_x;
}
   
Real
Bivariate_Gaussian_Distribution::get_mean_y () const
{
   return mean_y;
}
         
Real
Bivariate_Gaussian_Distribution::get_variance_x () const
{
   return variance_x;
}

Real
Bivariate_Gaussian_Distribution::get_variance_y () const
{
   return variance_y;
}

Real
Bivariate_Gaussian_Distribution::get_covariance () const
{
   return covariance;
}

Real
Bivariate_Gaussian_Distribution::get_sigma_x () const
{
   return sqrt (variance_x);
}

Real
Bivariate_Gaussian_Distribution::get_sigma_y () const
{
   return sqrt (variance_y);
}

Real
Bivariate_Gaussian_Distribution::get_rho () const
{
   return covariance / sqrt (variance_x * variance_y);
}

denise::Matrix*
Bivariate_Gaussian_Distribution::get_covariance_matrix_ptr () const
{

   Matrix* matrix_ptr = new Matrix (2, 2);

   matrix_ptr->set_datum (0, 0, variance_x);
   matrix_ptr->set_datum (1, 1, variance_y);
   matrix_ptr->set_datum (0, 1, covariance);
   matrix_ptr->set_datum (1, 0, covariance);

   return matrix_ptr;

}

denise::Matrix*
Bivariate_Gaussian_Distribution::get_inverse_covariance_matrix_ptr () const
{

   const Real determinant = variance_x*variance_y - covariance*covariance;
   Matrix* matrix_ptr = new Matrix (2, 2);

   const Real a_00 = variance_y / determinant;
   const Real a_11 = variance_x / determinant;
   const Real a_01 = -covariance / determinant;

   matrix_ptr->set_datum (0, 0, a_00);
   matrix_ptr->set_datum (1, 1, a_11);
   matrix_ptr->set_datum (0, 1, a_01);
   matrix_ptr->set_datum (1, 0, a_01);

   return matrix_ptr;

}

Ellipse
Bivariate_Gaussian_Distribution::get_ellipse (const Real probability) const
{

   const Real threshold = gsl_cdf_exponential_Pinv (probability, 2);

   const Matrix* inverse_sigma_ptr = get_inverse_covariance_matrix_ptr ();
   const Matrix& inverse_sigma = *inverse_sigma_ptr;

   const Real A = inverse_sigma.get_datum (0, 0);
   const Real B = inverse_sigma.get_datum (0, 1) * -2;
   const Real C = inverse_sigma.get_datum (1, 1);
   const Real F = -threshold;

   Ellipse ellipse (A, B, C, F);
   ellipse.translate (mean_x, mean_y);

   delete inverse_sigma_ptr;
   return ellipse;

}

Real
Bivariate_Gaussian_Distribution::get_pdf (const Point_2D& point_2d,
                                          const Real mean_x,
                                          const Real mean_y,        
                                          const Real variance_x,      
                                          const Real variance_y,      
                                          const Real covariance)
{
   return get_pdf (point_2d.x, point_2d.y, mean_x, mean_y,
                   variance_x, variance_y, covariance);
}

Real
Bivariate_Gaussian_Distribution::get_pdf (const Real x, 
                                          const Real y,
                                          const Real mean_x,
                                          const Real mean_y,
                                          const Real variance_x,
                                          const Real variance_y,
                                          const Real covariance)
{

   const Real sigma_x = sqrt (variance_x);
   const Real sigma_y = sqrt (variance_y);
   Real rho = covariance / (sigma_x * sigma_y);

   if (rho > 1) { rho = 0.9999; }
   if (rho < -1) { rho = -0.9999; }

   const Real pdf = gsl_ran_bivariate_gaussian_pdf (
      x - mean_x, y - mean_y, sigma_x, sigma_y, rho);
   return pdf;

}

Nvariate_Gaussian_Distribution::Nvariate_Gaussian_Distribution (Dataset& dataset)
   : dataset (dataset)
{
}

/*
Nvariate_Gaussian_Distribution::Nvariate_Gaussian_Distribution (const vector<Tuple>& tuple_vector)
{

   Integer vector_length = tuple_vector.size ();
   Integer vector_width = tuple_vector[0].size ();

   gsl_matrix* data = gsl_matrix_alloc (vector_length, vector_width);

   for (Integer i = 0; i < vector_length; i++)
   {

      Tuple tuple = tuple_vector[i];

      for (Integer j = 0; j < vector_width; j++)
      {
         gsl_matrix_set (data, i, j, tuple[j]);
      }

   }

   init (data);

   gsl_matrix_free (data);

}
*/

Integer
Nvariate_Gaussian_Distribution::get_size () const
{
   return dataset.get_dimension ();
}

Real
Nvariate_Gaussian_Distribution::get_pdf (const Vector& vector)
{
   const Integer size = get_size ();
   Matrix& sigma = dataset.get_covariance_matrix ();
   const Real& det_sigma = sigma.get_determinant ();
   const Real exponent = dataset.get_mahalonobis (vector) / -2;
   return exp (exponent) / sqrt (pow (M_2_TIMES_PI, size) * det_sigma);
}

Regression::Regression (const Integer n,
                        const Integer degrees_of_freedom)
   : n (n),
     degrees_of_freedom (degrees_of_freedom)
{
}
         
Integer
Regression::get_n () const
{
   return n;
}

Integer
Regression::get_degrees_of_freedom () const
{
   return degrees_of_freedom;
}

Real
Regression::get_chi_square () const
{
   return chi_square;
}

Real
Regression::get_error_variance () const
{
   return chi_square / (n - degrees_of_freedom - 1);
}

Real
Regression::get_regression_sum_of_squares () const
{
   return regression_sum_of_squares;
}

Real
Regression::get_SSr () const
{
   return get_regression_sum_of_squares ();
}

Real
Regression::get_SSe () const
{
   return get_chi_square ();
}

Real
Regression::get_Syy () const
{
   return get_SSr () + get_SSe ();
}

Real
Regression::get_MSr () const
{
   return get_SSr () / degrees_of_freedom;
}

Real
Regression::get_MSe () const
{
   return get_error_variance ();
}

Real
Regression::get_F0 () const
{
   return get_MSr () / get_MSe ();
}

Real
Regression::get_R_square () const
{
   return get_SSr () / get_Syy ();
}

void
Regression_nD::init (const Matrix& X,
                     const Vector& y)
{

   double chi_square_d;
   const Integer size = X.get_columns ();

   gsl_vector* gv = coefficient_vector.get_gv ();
   gsl_matrix* gm = covariance_matrix.get_gm ();

   gsl_multifit_linear_workspace* w = gsl_multifit_linear_alloc (n, size);
   gsl_multifit_linear (X.get_gm (), y.get_gv (), gv, gm, &chi_square_d, w);
   gsl_multifit_linear_free (w);

   this->chi_square = chi_square_d;

   this->regression_sum_of_squares = 0;
   const Real mean_y = gsl_stats_mean (y.get_gv ()->data, 1, n);
   for (Integer i = 0; i < n; i++)
   {
      gsl_vector_view x = gsl_matrix_row ((gsl_matrix*)X.get_gm (), i);
      const Real deviation = get_estimate (&(x.vector)) - mean_y;
      regression_sum_of_squares += deviation * deviation;
   }

}

Regression_nD::Regression_nD (const Integer n,
                              const Integer degrees_of_freedom)
   : Regression (n, degrees_of_freedom),
     coefficient_vector (degrees_of_freedom + 1),
     covariance_matrix (degrees_of_freedom + 1, degrees_of_freedom + 1)
{
}
         
Regression_nD::Regression_nD (const Vector& coefficient_vector,
                              const Matrix& covariance_matrix,
                              const Real chi_square,
                              const Real regression_sum_of_squares,
                              const Integer n)
   : Regression (n, coefficient_vector.get_size ()),
     coefficient_vector (coefficient_vector),
     covariance_matrix (covariance_matrix)
{
   this->chi_square = chi_square;
   this->regression_sum_of_squares = regression_sum_of_squares;
}

Regression_nD::Regression_nD (const Matrix& X,
                              const Vector& y)
   : Regression (X.get_rows (), X.get_columns () - 1),
     coefficient_vector (X.get_columns ()),
     covariance_matrix (X.get_columns (), X.get_columns ())
{
   init (X, y);
}

Regression_nD::Regression_nD (const vector<Tuple>& tuple_vector)
   : Regression (tuple_vector.size (), tuple_vector[0].size () - 1),
     coefficient_vector (tuple_vector[0].size ()),
     covariance_matrix (tuple_vector[0].size (), tuple_vector[0].size ())
{

   const Integer size = coefficient_vector.get_size ();
   Matrix X (n, size);
   Vector y (n);

   for (Integer i = 0; i < n; i++)
   {

      y.set_datum (i, tuple_vector[i][size-1]);
      X.set_datum (i, 0, 1);

      for (Integer j = 1; j < size; j++)
      {
         X.set_datum (i, j, tuple_vector[i][j-1]);
      }

   }

   init (X, y);

}

Regression_nD::~Regression_nD ()
{
}

const Vector&
Regression_nD::get_coefficient_vector () const
{
   return coefficient_vector;
}

const denise::Matrix&
Regression_nD::get_covariance_matrix () const
{
   return covariance_matrix;
}

Real
Regression_nD::get_ci_range (const Real alpha,
                             const Integer i) const
{
   const Real covariance = covariance_matrix.get_datum (i, i);
   return gsl_cdf_tdist_Qinv (alpha/2, n - 2) * sqrt (covariance);
}

Real
Regression_nD::get_t0 (const Integer i) const
{
   const Vector& cv = coefficient_vector;
   const Matrix& cm = covariance_matrix;
   return cv.get_datum (i) * sqrt (cm.get_datum (i, i));
}

Real
Regression_nD::get_estimate (const Vector& x) const
{
   return coefficient_vector.dot (x);
}

Real
Regression_nD::get_estimate_variance (const Vector& x) const
{

   const Integer n = covariance_matrix.get_rows ();
   Vector covariance_matrix_times_x (n, true);

   gsl_blas_dsymv (CblasUpper, 1, covariance_matrix.get_gm (),
      x.get_gv (), 0, covariance_matrix_times_x.get_gv ());

   const Real standard_error = x.dot (covariance_matrix_times_x);
   return standard_error;

}

Real
Regression_nD::get_estimate_standard_error (const Vector& x) const
{
   return sqrt (get_estimate_variance (x));
}

Gaussian_Distribution
Regression_nD::get_estimate_distribution (const Vector& x) const
{
   return Gaussian_Distribution (get_estimate (x), get_estimate_variance (x));
}

void
Regression_1D::init (const double* x,
                     const double* y,
                     const Integer n)
{

   double cs;
   double b_0, b_1;
   double cov_00, cov_01, cov_11;
   gsl_fit_linear (x, 1, y, 1, n, &b_0, &b_1, &cov_00, &cov_01, &cov_11, &cs);

   const Real mean_y = gsl_stats_mean (y, 1, n);

   this->beta_0 = b_0;
   this->beta_1 = b_1;
   this->chi_square = cs;
   this->covariance_00 = cov_00;
   this->covariance_01 = cov_01;
   this->covariance_11 = cov_11;
   this->regression_sum_of_squares = 0;

   for (Integer i = 0; i < n; i++)
   {
      const Real deviation = get_estimate (x[i]) - mean_y;
      regression_sum_of_squares += deviation * deviation;
   }

}

Regression_1D::Regression_1D (const double* x,
                              const double* y,
                              const Integer n)
   : Regression (n, 1)
{
   init (x, y, n);
}

Regression_1D::Regression_1D (const vector<Point_2D>& point_vector,
                              bool y_at_front)
   : Regression (n, 1)
{

   Integer n = point_vector.size ();
   double* x = new double[n];
   double* y = new double[n];

   for (Integer i = 0; i < n; i++)
   {
      x[i] = point_vector[i].x;
      y[i] = point_vector[i].y;
      if (y_at_front) { swap (x[i], y[i]); }
   }

   init (x, y, n);
   delete[] x;
   delete[] y;

}

Real
Regression_1D::get_beta_0 () const
{
   return beta_0;
}

Real
Regression_1D::get_beta_1 () const
{
   return beta_1;
}

Real
Regression_1D::get_covariance_00 () const
{
   return covariance_00;
}

Real
Regression_1D::get_covariance_01 () const
{
   return covariance_01;
}

Real
Regression_1D::get_covariance_11 () const
{
   return covariance_11;
}

Real
Regression_1D::get_ci_range_0 (const Real alpha) const
{
   return gsl_cdf_tdist_Qinv (alpha/2, n - 2) * sqrt (covariance_00);
}

Real
Regression_1D::get_ci_range_1 (const Real alpha) const
{
   return gsl_cdf_tdist_Qinv (alpha/2, n - 2) * sqrt (covariance_11);
}

Real
Regression_1D::get_estimate (const Real x) const
{
   return beta_0 + beta_1 * x;
}

Real
Regression_1D::get_estimate_variance (const Real x) const
{
   return covariance_00 + x * (2 * covariance_01 + covariance_11 * x);
}

Real
Regression_1D::get_estimate_standard_error (const Real x) const
{
   return sqrt (get_estimate_variance (x));
}

Gaussian_Distribution
Regression_1D::get_estimate_distribution (Real x) const
{
   return Gaussian_Distribution (get_estimate (x), get_estimate_variance (x));
}

Regression_2D::Regression_2D (const double* x0,
                              const double* x1,
                              const double* y,
                              const Integer n)
   : Regression_nD (n, 2)
{

   const Integer size = 3;

   Matrix X_matrix (n, size);
   Vector y_vector (n);

   for (Integer i = 0; i < n; i++)
   {
      y_vector.set_datum (i, y[i]);
      X_matrix.set_datum (i, 0, 1);
      X_matrix.set_datum (i, 1, x0[i]);
      X_matrix.set_datum (i, 2, x1[i]);
   }

   Regression_nD::init (X_matrix, y_vector);

}

Regression_2D::Regression_2D (const vector<Point_2D>& predictor_vector,
                              const vector<Real>& predictant_vector)
   : Regression_nD (n, 2)
{

   const Integer size = 3;
   const Integer n = predictor_vector.size ();

   Matrix X_matrix (n, size);
   Vector y_vector (n);

   for (Integer i = 0; i < n; i++)
   {
      y_vector.set_datum (i, predictant_vector[i]);
      X_matrix.set_datum (i, 0, 1);
      X_matrix.set_datum (i, 1, predictor_vector[i].x);
      X_matrix.set_datum (i, 2, predictor_vector[i].y);
   }

   Regression_nD::init (X_matrix, y_vector);

}

Regression_2D::Regression_2D (const vector<Point_3D>& point_vector)
   : Regression_nD (n, 2)
{

   const Integer size = 3;
   const Integer n = point_vector.size ();

   Matrix X_matrix (n, size);
   Vector y_vector (n);

   for (Integer i = 0; i < n; i++)
   {
      y_vector.set_datum (i, point_vector[i].y);
      X_matrix.set_datum (i, 0, 1);
      X_matrix.set_datum (i, 1, point_vector[i].z);
      X_matrix.set_datum (i, 2, point_vector[i].x);
   }

   Regression_nD::init (X_matrix, y_vector);

}

Real
Regression_2D::get_estimate (const Real x0,
                             const Real x1) const
{

   gsl_vector* x = gsl_vector_alloc (3);
   gsl_vector_set (x, 0, 1);
   gsl_vector_set (x, 1, x0);
   gsl_vector_set (x, 2, x1);

   Real estimate = Regression_nD::get_estimate (x);
   gsl_vector_free (x);

   return estimate;

}
   
Real
Regression_2D::get_estimate_variance (const Real x0,
                                      const Real x1) const
{

   gsl_vector* x = gsl_vector_alloc (3);
   gsl_vector_set (x, 0, 1);
   gsl_vector_set (x, 1, x0);
   gsl_vector_set (x, 2, x1);

   Real estimate = Regression_nD::get_estimate_variance (x);
   gsl_vector_free (x);

   return estimate;

}
   
Real
Regression_2D::get_estimate_standard_error (const Real x0,
                                            const Real x1) const
{
   return sqrt (get_estimate_variance (x0, x1));
}

Gaussian_Distribution
Regression_2D::get_estimate_distribution (const Real x0,
                                          const Real x1) const
{
   Real mean = get_estimate (x0, x1);
   Real variance = get_estimate_variance (x0, x1);
   return Gaussian_Distribution (mean, variance);
}

namespace denise
{

   ostream&
   operator << (ostream &out_file,
                Bgd& bgd)
   {
      out_file << "(" << bgd.get_mean () << ", " << bgd.get_variance_x ()
               << ", " << bgd.get_variance_y () << ", "
               << bgd.get_covariance () << ")";
      return out_file;
   }

}
