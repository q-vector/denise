//
// kalman.h
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

#ifndef DENISE_KALMAN_H
#define DENISE_KALMAN_H

#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

using namespace std;

namespace denise
{

   /// Implements the Kalman Filter. 
   ///
   /// Quoting from Welch and Bishop (2002),
   ///
   /// The Kalman Filter addresses the general problem of estimating
   /// the state x (element of \f$R^n\f$) of a discrete-time controlled
   /// process that is governed by the linear stochastic difference equation
   ///
   /// \f[
   /// x_k = A x_{k-1} + B u_k + w_{k-1}
   /// \f]
   ///
   ///  with a measurement z (elemennt of \f$R^m\f$), i.e.
   ///
   /// \f[
   /// z_k = H x_k + v_k
   /// \f]
   ///
   /// The random variables \f$w_k\f$ and \f$v_k\f$ represent the process
   /// and measurement noise respectively.  They are assumed to be
   /// independent of each other, white, and with normal probability
   /// distributions
   ///
   /// \f[
   ///     p(w) \sim N (0, Q)
   /// \f]
   /// \f[
   ///     p(v) \sim N (0, R)
   /// \f]
   ///
   /// The n-by-n matrix A in the difference equation (1) relates the
   ///  state at the previous timp step k-1 to the state at the current
   ///  step k, in the absence of either a driving function or process
   ///  noise.  The n-by-l matrix B relates matrix B relates the optional
   ///  control input u (element of \f$R^l\f$) to the state x.  The m-by-n
   ///  matrix H in the measurement equation (2) relates the state to
   ///  the measurement \f$z_k\f$.
   ///
   class Kalman_Filter
   {

      private:

         int
         n;

         int
         m;

         gsl_matrix*
         H;

         gsl_matrix*
         Q;

         gsl_matrix*
         R;

         void
         time_update (gsl_vector* x,
                      gsl_matrix* P);

         void
         time_update (gsl_vector* x,
                      gsl_matrix* P,
                      const gsl_matrix* A);

         void
         time_update (gsl_vector* x,
                      gsl_matrix* P,
                      const gsl_matrix* B,
                      const gsl_vector* u);

         void
         time_update (gsl_vector* x,
                      gsl_matrix* P,
                      const gsl_matrix* A,
                      const gsl_matrix* B,
                      const gsl_vector* u);

         void
         measurement_update (gsl_vector* x,
                             gsl_matrix* P,
                             const gsl_vector* z);

      public:

         /// Constructor
         ///
         /// \param H m-by-n matrix
         /// \param Q n-by-n matrix
         /// \param R m-by-m matrix
         Kalman_Filter (const gsl_matrix* H,
                        const gsl_matrix* Q,
                        const gsl_matrix* R);

         /// Destructor
         ~Kalman_Filter ();

         /// Iterates kalman_filter by one step.
         ///
         /// A is assumed to be I
         /// B is assumed to be 0
         ///
         /// \param x n vector
         /// \param P n-by-n matrix
         /// \param z m vector
         void
         filter (gsl_vector* x,
                 gsl_matrix* P,
                 const gsl_vector* z);

         /// Iterates kalman_filter by one step.
         ///
         /// B is assumed to be 0
         ///
         /// \param x n vector
         /// \param P n-by-n matrix
         /// \param A n-by-n matrix
         /// \param z m vector
         void
         filter (gsl_vector* x,
                 gsl_matrix* P,
                 const gsl_matrix* A,
                 const gsl_vector* z);

         /// Iterates kalman_filter by one step.
         ///
         /// A assumed to be I.
         ///
         /// \param x n vector
         /// \param P n-by-n matrix
         /// \param z m vector
         /// \param B n-by-l matrix
         /// \param u l vector
         void
         filter (gsl_vector* x,
                 gsl_matrix* P,
                 const gsl_matrix* B,
                 const gsl_vector* u,
                 const gsl_vector* z);

         /// Iterates kalman_filter by one step.
         ///
         /// \param x n vector
         /// \param P n-by-n matrix
         /// \param A n-by-n matrix
         /// \param z m vector
         /// \param B n-by-l matrix
         /// \param u l vector
         void
         filter (gsl_vector* x,
                 gsl_matrix* P,
                 const gsl_matrix* A,
                 const gsl_matrix* B,
                 const gsl_vector* u,
                 const gsl_vector* z);

   };

}

#endif /* DENISE_KALMAN_H */

