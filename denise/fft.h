//
// fft.h
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

#ifndef DENISE_FFT_H
#define DENISE_FFT_H

#include <fftw3.h>
#include <denise/basics.h>

using namespace std;
using namespace denise;

namespace denise
{

   enum Fft_Scaling
   {
      UNITARY,
      FORWARD,
      BACKWARD
   };

   class Fft_Real_1D
   {

      private:

         fftw_plan
         plan;

         Fft_Scaling
         fft_scaling;

         void
         scale ();

         void
         scale_array_t (const Real scale);

         void
         scale_array_f (const Real scale);

      public:

         const Integer
         n;

         const bool
         forward;

         double*
         array_t;

         fftw_complex*
         array_f;

         Fft_Real_1D (const Integer n,
                      const bool forward = true,
                      const Fft_Scaling fft_scaling = FORWARD,
                      const unsigned plan_vigor = FFTW_ESTIMATE);

         ~Fft_Real_1D ();

         void
         set_array_t (const Tuple& tuple_t);

         void
         set_array_t (const double* array_t);

         void
         set_coefficient (const Integer k,
                          const Cmplx& coefficient);

         void
         set_amplitude_phase (const Integer wavenumber,
                              const Real amplitude,
                              const Real phase);

         void
         set_trig_coefficients (const Integer wavenumber,
                                const Real cosine_coefficient,
                                const Real sine_coefficient);

         void
         transform ();

         Cmplx
         get_coefficient (const Integer k) const;

         void
         acquire_amplitude_phase (const Integer wavenumber,
                                  Real& amplitude,
                                  Real& phase) const;

         void
         acquire_trig_coefficients (const Integer wavenumber,
                                    Real& cosine_coefficient,
                                    Real& sine_coefficient) const;

   };

}

#endif /* DENISE_FFT_H */

