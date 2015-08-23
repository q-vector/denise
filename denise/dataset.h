//
// dataset.h
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

#ifndef DENISE_DATASET_H
#define DENISE_DATASET_H

#include <denise/basics.h>
#include <denise/linalg.h>

using namespace std;

namespace denise
{

   class Dataset
   {

      protected:

         Matrix*
         data_matrix_ptr;

         Matrix*
         covariance_matrix_ptr;

         Vector*
         mean_vector_ptr;

         Matrix*
         residue_matrix_ptr;

         Real_Symmetric_Eigen*
         principal_component_ptr;

      public:

         Dataset (const Integer dimension,
                  const Integer number_of_points);

         Dataset (Matrix* data_matrix_ptr);

         ~Dataset ();

         void
         set_data (const Integer i,
                   const Real* data);

         void
         set_datum (const Integer i,
                    const Integer j,
                    const Real datum);

         Integer
         get_number_of_points () const;

         Integer
         get_dimension () const;

         Matrix&
         get_data_matrix ();

         const Matrix&
         get_data_matrix () const;

         Matrix&
         get_covariance_matrix ();

         const Vector&
         get_mean_vector ();

         const Matrix&
         get_residue_matrix ();

         const Vector&
         get_pca_coefficient_vector ();

         const Matrix&
         get_pca_vector_matrix ();

         const Real
         get_mahalonobis (const Vector& vector);

         void
         compute_residue (const bool idle_if_present = true);

         void
         compute_covariance (const bool idle_if_present = true);

         void
         compute_mean (const bool idle_if_present = true);

         void
         compute_principal_component (const bool sort = true,
                                      const bool idle_if_present = true);

   };

   class Dataset_1D : public Dataset
   {

      private:

         static int
         compare (const void* a,
                  const void* b);

      public:

         Dataset_1D (const Integer number_of_points);

         void
         set_datum (const Integer i,
                    const Real datum);

         void
         sort ();

         Real
         get_mean ();

         Real
         get_variance ();

         Real
         get_mean_square () const;

         Real
         get_root_mean_square () const;

         Real
         get_rms () const;

         Real
         get_sd ();

         Real
         get_standard_deviation ();

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

   class Dataset_2D : public Dataset
   {

      public:

         Dataset_2D (const Integer number_of_points);

         void
         set_data (const Integer i,
                   const Real datum_x,
                   const Real datum_y);

         Real
         get_mean_x ();

         Real
         get_mean_y ();

         Real
         get_variance_x ();

         Real
         get_variance_y ();

         Real
         get_covariance ();

   };

};

#endif /* DENISE_DATASET_H */
