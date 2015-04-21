//
// stat.h
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

#ifndef DENISE_STAT_H
#define DENISE_STAT_H

#include <cmath>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_statistics_double.h>
#include <denise/basics.h>
#include <denise/geometry.h>
#include <denise/linalg.h>
#include <denise/dataset.h>

using namespace std;

namespace denise
{

   class Ellipse;

   class Treatment : public Tuple
   {

      private:

         double*
         data;

         Integer
         size;

         static int
         compare (const void* a,
                  const void* b);

      public:

         Treatment (const double* data,
                    const Integer size);

         Treatment (const Tuple& tuple);

         ~Treatment ();

         void
         sort ();

         const double*
         get_data () const;

         Integer
         get_size () const;

         Real
         get_mean () const;

         Real
         get_mean_square () const;

         Real
         get_root_mean_square () const;

         Real
         get_rms () const;

         Real
         get_variance () const;

         Real
         get_sd () const;

         Real
         get_standard_deviation () const;

         Real
         get_max () const;

         Real
         get_min () const;

         Real
         get_median () const;

         Real
         get_first_quartile () const;

         Real
         get_second_quartile () const;

         Real
         get_third_quartile () const;

         Real
         get_percentile (const Real percentage) const;

   };

   class Lognormal_Distribution
   {

      private:

         Real
         mean;

         Real
         variance;

         void
         init (Real mean,
               Real variace);

         void
         init (const double* array,
               Integer n);

      public:

         Lognormal_Distribution ();

         Lognormal_Distribution (const Real mean,
                                 const Real variance);

         Lognormal_Distribution (const double* array,
                                 const Integer n);

         Lognormal_Distribution (const Tuple& tuple);

         Real
         get_pdf (const Real point) const;

         Real
         get_P (const Real point) const;

         Real
         get_Q (const Real point) const;

         Real
         get_P_inv (const Real P) const;

         Real
         get_Q_inv (const Real Q) const;

         Real
         get_mean () const;

         Real
         get_variance () const;

         Real
         get_sigma () const;

         static Real
         get_pdf (const Real point,
                  const Real mean,
                  const Real variance);

         static Real
         get_P (const Real point,
                const Real mean,
                const Real variance);

         static Real
         get_Q (const Real point,
                const Real mean,
                const Real variance);

         static Real
         get_P_inv (const Real P,
                    const Real mean,
                    const Real variance);

         static Real
         get_Q_inv (const Real Q,
                    const Real mean,
                    const Real variance);

   };

   class Rayleigh_Distribution
   {

      private:

         Real
         sigma;

      public:

         Rayleigh_Distribution (const Real sigma);

         Real
         get_pdf (const Real point) const;

         Real
         get_P (const Real point) const;

         Real
         get_Q (const Real point) const;

         Real
         get_P_inv (const Real P) const;

         Real
         get_Q_inv (const Real Q) const;

         Real
         get_sigma () const;

         static Real
         get_pdf (const Real point,
                  const Real sigma);

         static Real
         get_P (const Real point,
                const Real sigma);

         static Real
         get_Q (const Real point,
                const Real sigma);

         static Real
         get_P_inv (const Real P,
                    const Real sigma);

         static Real
         get_Q_inv (const Real Q,
                    const Real sigma);

   };

   /// Represents the Gaussian Distribution
   ///
   /// A Gaussian Distribution in a variate \f$X\f$ with a mean \f$\mu\f$
   /// and variance \f$\sigma^2\f$ has the probability density function
   ///
   /// \f[
   ///   \frac{1}{\sqrt {2 \pi \sigma^2}} \exp (-(x - \mu)^2 / 2 \sigma^2)
   /// \f]
   class Gaussian_Distribution
   {

      private:

         Real
         mean;

         Real
         variance;

         void
         init (const Real mean,
               const Real variance);

         void
         init (const double* array,
               const Integer n);

      public:

         Gaussian_Distribution ();

         /// Constructor
         ///
         /// \param mean mean
         /// \param variance variance
         Gaussian_Distribution (const Real mean,
                                const Real variance);

         Gaussian_Distribution (Dataset_1D& dataset_1d);

         Gaussian_Distribution (const double* array,
                                const Integer n);

         /// Constructor
         ///
         /// \param tuple tuple holding the sample data
         Gaussian_Distribution (const Tuple& tuple);

         /// Evaulates the probability density function
         Real
         get_pdf (const Real point) const;

         /// Evaluates the lower tail cumulative probability function
         Real
         get_P (const Real point) const;

         /// Evaluates the upper tail cumulative probability function
         Real
         get_Q (const Real point) const;

         /// Evaluates the inverse of the lower tail cumulative
         /// probability function
         Real
         get_P_inv (const Real P) const;

         /// Evaluates the inverse of the upper tail cumulative
         /// probability function
         Real
         get_Q_inv (const Real Q) const;

         /// Returns the mean of the Gaussian Distribution
         Real
         get_mean () const;

         /// Returns the variance of the Gaussian Distribution
         Real
         get_variance () const;

         /// Returns the standard deviation of the Gaussian Distribution
         Real
         get_sigma () const;

         static Real
         get_pdf (const Real point,
                  const Real mean,
                  const Real variance);

         static Real
         get_P (const Real point,
                const Real mean,
                const Real variance);

         static Real
         get_Q (const Real point,
                const Real mean,
                const Real variance);

         static Real
         get_P_inv (const Real P,
                    const Real mean,
                    const Real variance);

         static Real
         get_Q_inv (const Real Q,
                    const Real mean,
                    const Real variance);

   };

   typedef Gaussian_Distribution Normal_Distribution;

   /// Represents the Gaussian Distribution
   ///
   /// A Bivariate Gaussian Distribution can be specified by 5 attributes:
   /// 2 means \f$\mu_x\f$, \f$\mu_y\f$, 2 variances \f$\sigma_x^2\f$,
   /// \f$\sigma_y_2\f$ and a covariance \f$\sigma_{xy}\f$.  It has
   /// the probability denisty function:
   ///
   /// \f[
   ///   \frac{1}{2 \pi \sigma_x \sigma_y \sqrt {1 - \rho^2}}
   ///    \exp (- \frac{z}{2 (1 - \rho^2)})
   /// \f]
   ///
   /// where
   ///
   /// \f[
   ///    z = \frac{(x - \mu_x)^2}{\sigma_x^2}
   ///      - \frac{2 \rho (x - \mu_x) (y - \mu_y)}{\sigma_x \sigma_y}
   ///      + \frac{(y - \mu_y)^2}{\sigma_y^2}
   /// \f]
   ///
   /// and
   ///
   /// \rho = \frac{\sigma_{xy}}{\sigma_x \sigma_y}
   class Bivariate_Gaussian_Distribution
   {

      private:

         Real
         mean_x;

         Real
         mean_y;

         Real
         variance_x;

         Real
         variance_y;

         Real
         covariance;

         void
         init (const Real mean_x,
               const Real mean_y,
               const Real variance_x,
               const Real variance_y,
               const Real covariance);

         void
         init (const double* x_array,
               const double* y_array,
               const Integer n);

      public:

         Bivariate_Gaussian_Distribution ();

         /// Constructor
         ///
         /// \param mean_x mean in the x direction
         /// \param mean_y mean in the y direction
         /// \param variance_x variance in the x direction
         /// \param variance_y variance in the y direction
         /// \param covariance covariance
         Bivariate_Gaussian_Distribution (const Real mean_x,
                                          const Real mean_y,
                                          const Real variance_x,
                                          const Real variance_y,
                                          const Real covariance);

         Bivariate_Gaussian_Distribution (Dataset_2D& dataset_2d);

         Bivariate_Gaussian_Distribution (const double* x_array,
                                          const double* y_array,
                                          const Integer n);

         /// Constructor
         ///
         /// \param point_vector vector holding the sample data
         Bivariate_Gaussian_Distribution (const vector<Point_2D>& point_vector);

         Bivariate_Gaussian_Distribution (const Duple_Vector& duple_vector);

         /// Evaulates the probability density function
         Real
         get_pdf (const Point_2D& point_2d) const;

         /// Evaulates the probability density function
         Real
         get_pdf (const Real x,
                  const Real y) const;

         /// Returns the mean of the Bivariate Gaussian Distribution
         Point_2D
         get_mean () const;

         /// Returns the x-mean of the Bivariate Gaussian Distribution
         Real
         get_mean_x () const;

         /// Returns the y-mean of the Bivariate Gaussian Distribution
         Real
         get_mean_y () const;

         /// Returns the x-variance of the Bivariate Gaussian Distribution
         Real
         get_variance_x () const;

         /// Returns the y-variance of the Bivariate Gaussian Distribution
         Real
         get_variance_y () const;

         /// Returns the covariance of the Bivariate Gaussian Distribution
         Real
         get_covariance () const;

         /// Returns the x standard deviation of the
         /// Bivariate Gaussian Distribution
         Real
         get_sigma_x () const;

         /// Returns the y standard deviation of the
         /// Bivariate Gaussian Distribution
         Real
         get_sigma_y () const;

         /// Returns the correlation of the Bivariate Gaussian Distribution
         Real
         get_rho () const;

         /// Returns the covariance matrix of the
         /// Bivariate Gaussian Distribution
         Matrix*
         get_covariance_matrix_ptr () const;

         /// Returns the inverse of the covariance matrix of the
         /// Bivariate Gaussian Distribution
         Matrix*
         get_inverse_covariance_matrix_ptr () const;

         /// Returns an ellipse that encloses the given probability
         /// of mass of the Bivariate Gaussian Distribution
         Ellipse
         get_ellipse (const Real probability) const;

         static Real
         get_pdf (const Point_2D& point_2d,
                  const Real mean_x,
                  const Real mean_y,
                  const Real variance_x,
                  const Real variance_y,
                  const Real covariance);

         static Real
         get_pdf (const Real x,
                  const Real y,
                  const Real mean_x,
                  const Real mean_y,        
                  const Real variance_x,      
                  const Real variance_y,      
                  const Real covariance);

   };

   typedef Bivariate_Gaussian_Distribution Bivariate_Normal_Distribution;
   typedef Bivariate_Gaussian_Distribution Bgd;
   typedef Bivariate_Normal_Distribution Bnd;

   /// Represents the N-variate Gaussian Distribution, or the
   /// Multivariate Gaussian Distribution.
   ///
   /// N-variate Gaussian Distribution is a generalization of the
   /// bivariate gaussian distribution.  It is specified by a mean
   /// vector \f$\mathbf{\mu}\f$ and covariance matrix \f$\mathbf{\Sigma}\f$.
   class Nvariate_Gaussian_Distribution
   {

      private:

         Dataset&
         dataset;

      public:

         /// Constructor
         Nvariate_Gaussian_Distribution (Dataset& dataset);

         /// Constructor
         ///
         /// \param tuple_vector vector holding the sample data
         //Nvariate_Gaussian_Distribution (const vector<Tuple>& tuple_vector);

         ~Nvariate_Gaussian_Distribution ();

         /// Returns the number of dimensions of the
         /// N-variate Gaussian Distribution
         Integer
         get_size () const;

         /// Evaulates the probability density function
         Real
         get_pdf (const Vector& vector);

   };

   typedef Nvariate_Gaussian_Distribution Nvariate_Normal_Distribution;
   typedef Nvariate_Gaussian_Distribution Ngd;
   typedef Nvariate_Normal_Distribution Nnd;

   class Regression
   {

      protected:

         Integer
         n;

         Real
         chi_square;

         Integer
         degrees_of_freedom;

         Real
         regression_sum_of_squares;

      public:

         //Regression ();

         Regression (const Integer n,
                     const Integer degrees_of_freedom);

         Integer
         get_n () const;

         Integer
         get_degrees_of_freedom () const;

         Real 
         get_chi_square () const;

         Real
         get_error_variance () const;

         Real
         get_regression_sum_of_squares () const;

         Real
         get_SSr () const;

         Real
         get_SSe () const;

         Real
         get_Syy () const;

         Real
         get_MSr () const;

         Real
         get_MSe () const;

         Real
         get_F0 () const;

         Real
         get_R_square () const;

   };

   class Regression_nD : public Regression
   {

      protected:

         Vector
         coefficient_vector;

         Matrix
         covariance_matrix;

         void
         init (const Matrix& X,
               const Vector& y);

      public:

         //Regression_nD () { }

         Regression_nD (const Integer n,
                        const Integer degrees_of_freedom);

         Regression_nD (const Vector& coefficient_vector,
                        const Matrix& covariance_matrix,
                        const Real chi_square,
                        const Real regression_sum_of_squares,
                        const Integer n);

         Regression_nD (const Matrix& X,
                        const Vector& y);

         Regression_nD (const vector<Tuple>& tuple_vector);

         ~Regression_nD ();

         const Vector&
         get_coefficient_vector () const;

         const Matrix&
         get_covariance_matrix () const;

         Real
         get_ci_range (const Real alpha,
                       const Integer i) const;

         Real
         get_t0 (const Integer i) const;

         Real
         get_estimate (const Vector& x) const;

         Real
         get_estimate_variance (const Vector& x) const;

         Real
         get_estimate_standard_error (const Vector& y) const;

         Gaussian_Distribution
         get_estimate_distribution (const Vector& x) const;

   };

   class Regression_1D : public Regression
   {

      private:

         Real
         beta_0;

         Real
         beta_1;

         Real
         covariance_00;

         Real
         covariance_01;

         Real
         covariance_11;

         void
         init (const double* x,
               const double* y,
               const Integer n);

      public:

         Regression_1D (const double* x,
                        const double* y,
                        const Integer n);

         Regression_1D (const vector<Point_2D>& point_vector,
                        bool y_at_front = false);

         Real
         get_beta_0 () const;

         Real
         get_beta_1 () const;

         Real
         get_covariance_00 () const;

         Real
         get_covariance_01 () const;

         Real
         get_covariance_11 () const;

         Real
         get_ci_range_0 (const Real alpha) const;

         Real
         get_ci_range_1 (const Real alpha) const;

         Real
         get_estimate (const Real x) const;

         Real
         get_estimate_variance (const Real x) const;

         Real
         get_estimate_standard_error (const Real x) const;

         Gaussian_Distribution
         get_estimate_distribution (const Real x) const;

   };

   class Regression_2D : public Regression_nD
   {

      public:

         Regression_2D (const double* x0,
                        const double* x1,
                        const double* y,
                        const Integer n);

         Regression_2D (const vector<Point_2D>& predictor_vector,
                        const vector<Real>& predictant_vector);

         Regression_2D (const vector<Point_3D>& point_vector);

         Real
         get_estimate (const Real x0,
                       const Real x1) const;

         Real
         get_estimate_variance (const Real x0,
                                const Real x1) const;

         Real
         get_estimate_standard_error (const Real x0,
                                      const Real x1) const;

         Gaussian_Distribution
         get_estimate_distribution (const Real x0,
                                    const Real x1) const;

   };

}

#endif /* DENISE_STAT_H */

