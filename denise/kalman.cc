//
// kalman.cc
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

#include <kalman.h>

using namespace std;
using namespace denise;

void
Kalman_Filter::time_update (gsl_vector* x,
                            gsl_matrix* P)
{

   // P = P + Q
   //
   // This step is composed by the following steps:
   // - P += Q
   gsl_matrix_add (P, Q);

}

void
Kalman_Filter::time_update (gsl_vector* x,
                            gsl_matrix* P,
                            const gsl_matrix* A)
{

   // Allocate
   gsl_vector* Ax = gsl_vector_alloc (n);
   gsl_matrix* AP = gsl_matrix_alloc (n, n);
   gsl_matrix* APAt = gsl_matrix_alloc (n, n);

   // x = Ax
   //
   // This step is composed by the following steps:
   // - Ax = A * x
   // - x = Ax
   gsl_blas_dgemv (CblasNoTrans, 1, A, x, 0, Ax);
   gsl_vector_memcpy (x, Ax);

   // P = APAt + Q
   //
   // This step is composed by the following steps:
   // - AP = A * P
   // - APAt = AP*At
   // - P = APAt
   // - P += Q
   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, A, P, 0, AP);
   gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1, AP, A, 0, APAt);
   gsl_matrix_memcpy (P, APAt);
   gsl_matrix_add (P, Q);

   // Clean-up
   gsl_vector_free (Ax);
   gsl_matrix_free (AP);
   gsl_matrix_free (APAt);

}

void
Kalman_Filter::time_update (gsl_vector* x,
                            gsl_matrix* P,
                            const gsl_matrix* B,
                            const gsl_vector* u)
{

   // Allocate
   gsl_vector* Bu = gsl_vector_alloc (n);
   gsl_matrix* AP = gsl_matrix_alloc (n, n);

   // x = Ax + Bu
   //
   // This step is composed by the following steps:
   // - Bu = B * u
   // - x += Bu
   gsl_blas_dgemv (CblasNoTrans, 1, B, u, 0, Bu);
   gsl_vector_add (x, Bu);

   // P = P + Q
   //
   // This step is composed by the following steps:
   // - P += Q
   gsl_matrix_add (P, Q);

   // Clean-up
   gsl_vector_free (Bu);
   gsl_matrix_free (AP);

}

void
Kalman_Filter::time_update (gsl_vector* x,
                            gsl_matrix* P,
                            const gsl_matrix* A,
                            const gsl_matrix* B,
                            const gsl_vector* u)
{

   // Allocate
   gsl_vector* Ax = gsl_vector_alloc (n);
   gsl_vector* Bu = gsl_vector_alloc (n);
   gsl_matrix* AP = gsl_matrix_alloc (n, n);
   gsl_matrix* APAt = gsl_matrix_alloc (n, n);

   // x = Ax + Bu
   //
   // This step is composed by the following steps:
   // - Ax = A * x
   // - Bu = B * u
   // - x = Ax
   // - x += Bu
   gsl_blas_dgemv (CblasNoTrans, 1, A, x, 0, Ax);
   gsl_blas_dgemv (CblasNoTrans, 1, B, u, 0, Bu);
   gsl_vector_memcpy (x, Ax);
   gsl_vector_add (x, Bu);

   // P = APAt + Q
   //
   // This step is composed by the following steps:
   // - AP = A * P
   // - APAt = AP*At
   // - P = APAt
   // - P += Q
   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, A, P, 0, AP);
   gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1, AP, A, 0, APAt);
   gsl_matrix_memcpy (P, APAt);
   gsl_matrix_add (P, Q);

   // Clean-up
   gsl_vector_free (Ax);
   gsl_vector_free (Bu);
   gsl_matrix_free (AP);
   gsl_matrix_free (APAt);

}

void
Kalman_Filter::measurement_update (gsl_vector* x,
                                   gsl_matrix* P,
                                   const gsl_vector* z)
{

   int signum;

   // Allocate
   gsl_vector* Hx = gsl_vector_alloc (m);
   gsl_vector* ozmHxc = gsl_vector_alloc (m);
   gsl_vector* KozmHxc = gsl_vector_alloc (n);
   gsl_matrix* K = gsl_matrix_alloc (n, m);
   gsl_matrix* LU = gsl_matrix_alloc (m, m);
   gsl_matrix* KH = gsl_matrix_alloc (n, n);
   gsl_matrix* PHt = gsl_matrix_alloc (n, m);
   gsl_matrix* HPHt = gsl_matrix_alloc (m, m);
   gsl_matrix* oImKHc = gsl_matrix_alloc (n, n);
   gsl_matrix* oImKHcP = gsl_matrix_alloc (n, n);
   gsl_matrix* oHPHtpRc = gsl_matrix_alloc (m, m);
   gsl_matrix* oHPHtpRci = gsl_matrix_alloc (m, m);
   gsl_permutation* permutation = gsl_permutation_alloc (m);

   // K = P * Ht * oHPHtpRci
   //
   // This step is composed by the following steps:
   // - PHt = P * Ht
   // - HPHt = H * PHt
   // - oHPHtpRc = HPHt
   // - oHPHtpRc += R
   // - LU = oHPHtpRc
   // - LU = LU_decomp (oHPHtpRc)
   //   calculate oHPHtpRci
   // - K = oHPHtpRci
   gsl_blas_dgemm (CblasNoTrans, CblasTrans, 1, P, H, 0, PHt);
   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, H, PHt, 0, HPHt);
   gsl_matrix_memcpy (oHPHtpRc, HPHt);
   gsl_matrix_add (oHPHtpRc, R);
   gsl_matrix_memcpy (LU, oHPHtpRc);
   gsl_linalg_LU_decomp (LU, permutation, &signum);
   gsl_linalg_LU_invert (LU, permutation, oHPHtpRci);
   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, PHt, oHPHtpRci, 0, K);

   // x = x + KozmHxc
   //
   // This step is composed by the following steps:
   // - Hx = H * x
   // - ozmHxc = z 
   // - ozmHxc -= Hx
   // - KozmHxc = K * omzHxc
   // - x += KozmHxc

   gsl_blas_dgemv (CblasNoTrans, 1, H, x, 0, Hx);
   gsl_vector_memcpy (ozmHxc, z);
   gsl_vector_sub (ozmHxc, Hx);
   gsl_blas_dgemv (CblasNoTrans, 1, K, ozmHxc, 0, KozmHxc);
   gsl_vector_add (x, KozmHxc);

   // P = oImKHc * P
   //
   // This step is composed by the following steps:
   // - KH = K * H
   // - oImKHc = I
   // - oImKHc -= KH
   // - oImKHcP = oImKHc * P
   // - P = oImKHcP
   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, K, H, 0, KH);
   gsl_matrix_set_identity (oImKHc);
   gsl_matrix_sub (oImKHc, KH);
   gsl_blas_dgemm (CblasNoTrans, CblasNoTrans, 1, oImKHc, P, 0, oImKHcP);
   gsl_matrix_memcpy (P, oImKHcP);

   // Clean-up
   gsl_vector_free (Hx);
   gsl_vector_free (ozmHxc);
   gsl_vector_free (KozmHxc);
   gsl_matrix_free (K);
   gsl_matrix_free (LU);
   gsl_matrix_free (KH);
   gsl_matrix_free (PHt);
   gsl_matrix_free (HPHt);
   gsl_matrix_free (oImKHc);
   gsl_matrix_free (oImKHcP);
   gsl_matrix_free (oHPHtpRc);
   gsl_matrix_free (oHPHtpRci);
   gsl_permutation_free (permutation);

}

Kalman_Filter::Kalman_Filter (const gsl_matrix* H,
                              const gsl_matrix* Q,
                              const gsl_matrix* R)
{

   this->n = Q->size1;
   this->m = R->size1;

   this->H = gsl_matrix_alloc (m, n);
   this->Q = gsl_matrix_alloc (n, n);
   this->R = gsl_matrix_alloc (m, m);

   gsl_matrix_memcpy (this->H, H);
   gsl_matrix_memcpy (this->Q, Q);
   gsl_matrix_memcpy (this->R, R);

}

Kalman_Filter::~Kalman_Filter ()
{
   gsl_matrix_free (this->H);
   gsl_matrix_free (this->Q);
   gsl_matrix_free (this->R);
}

void
Kalman_Filter::filter (gsl_vector* x,
                       gsl_matrix* P,
                       const gsl_vector* z)
{
   time_update (x, P);
   measurement_update (x, P, z);
}

void
Kalman_Filter::filter (gsl_vector* x,
                       gsl_matrix* P,
                       const gsl_matrix* A,
                       const gsl_vector* z)
{
   time_update (x, P, A);
   measurement_update (x, P, z);
}

void
Kalman_Filter::filter (gsl_vector* x,
                       gsl_matrix* P,
                       const gsl_matrix* B,
                       const gsl_vector* u,
                       const gsl_vector* z)
{
   time_update (x, P, B, u);
   measurement_update (x, P, z);
}

void
Kalman_Filter::filter (gsl_vector* x,
                       gsl_matrix* P,
                       const gsl_matrix* A,
                       const gsl_matrix* B,
                       const gsl_vector* u,
                       const gsl_vector* z)
{
   time_update (x, P, A, B, u);
   measurement_update (x, P, z);
}

