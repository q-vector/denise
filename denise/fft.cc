//
// fft.cc
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

#include "fft.h"
#include "util.h"

using namespace std;
using namespace denise;

void
Fft_Real_1D::scale ()
{

   if (forward)
   {
      switch (fft_scaling)
      {
         case UNITARY: scale_array_f (sqrt (Real (n))); break;
         case FORWARD: scale_array_f (Real (n));        break;
      }
   }
   else
   {
      switch (fft_scaling)
      {
         case UNITARY:  scale_array_t (sqrt (Real (n))); break;
         case BACKWARD: scale_array_t (Real (n));        break;
      }
   }

}

void
Fft_Real_1D::scale_array_t (const Real scale)
{
   for (Integer i = 0; i < n; i++)
   {
      array_t[i] /= scale;
   }
}

void
Fft_Real_1D::scale_array_f (const Real scale)
{
   for (Integer i = 0; i < n/2+1; i++)
   {
      array_f[i][0] /= scale;
      array_f[i][1] /= scale;
   }
}

Fft_Real_1D::Fft_Real_1D (const Integer n,
                          const bool forward,
                          const Fft_Scaling fft_scaling,
                          const unsigned plan_vigor)
                     : n (n),
                 forward (forward),
             fft_scaling (fft_scaling)
{

   array_t = (double*)fftw_malloc (sizeof (double) * n); 
   array_f = (fftw_complex*)fftw_malloc (sizeof (fftw_complex) * (n/2+1)); 

   const bool& f = forward;
   if (f) { plan = fftw_plan_dft_r2c_1d (n, array_t, array_f, plan_vigor); }
   else   { plan = fftw_plan_dft_c2r_1d (n, array_f, array_t, plan_vigor); }

}

Fft_Real_1D::~Fft_Real_1D ()
{
   fftw_free (array_f);
   fftw_free (array_t); 
   fftw_destroy_plan (plan); 
}

void
Fft_Real_1D::set_array_t (const Tuple& tuple_t)
{
   for (Integer i = 0; i < n; i++) { array_t[i] = tuple_t[i]; }
}

void
Fft_Real_1D::set_array_t (const double* array_t)
{
   memcpy (this->array_t, array_t, sizeof (double) * n);
}

void
Fft_Real_1D::set_coefficient (const Integer k,
                              const Cmplx& coefficient)
{
   fftw_complex& fftwc = array_f[abs (k)];
   fftwc[0] = coefficient.real ();
   fftwc[1] = (k < 0 ? -1 : 1) * coefficient.imag ();
}

void
Fft_Real_1D::set_amplitude_phase (const Integer wavenumber,
                                  const Real amplitude,
                                  const Real phase)
{
   const Real r = amplitude * cos (phase) / 2;
   const Real i = amplitude * sin (phase) / -2;
   set_coefficient (wavenumber, Cmplx (r, i));
}

void
Fft_Real_1D::set_trig_coefficients (const Integer wavenumber,
                                    const Real cosine_coefficient,
                                    const Real sine_coefficient)
{
   const Real r = cosine_coefficient / 2;
   const Real i = sine_coefficient / -2;
   set_coefficient (wavenumber, Cmplx (r, i));
}

void
Fft_Real_1D::transform ()
{
   fftw_execute (plan);
   scale ();
}

Cmplx
Fft_Real_1D::get_coefficient (const Integer k) const
{
   const fftw_complex& fftwc = array_f[abs (k)];
   return Cmplx (fftwc[0], (k < 0 ? -1 : 1) * fftwc[1]);
}

void
Fft_Real_1D::acquire_amplitude_phase (const Integer wavenumber,
                                      Real& amplitude,
                                      Real& phase) const
{
   const Cmplx z = get_coefficient (-imodulo (wavenumber, n/2+1));
   const Real real = z.real ();
   const Real imag = z.imag ();
   amplitude = 2 * sqrt (real*real + imag*imag);
   phase = atan2 (imag, real);
}

void
Fft_Real_1D::acquire_trig_coefficients (const Integer wavenumber,
                                        Real& cosine_coefficient,
                                        Real& sine_coefficient) const
{
   const Cmplx z = get_coefficient (imodulo (wavenumber, n/2+1));
   cosine_coefficient = 2 * z.real ();
   sine_coefficient = -2 * z.imag ();
}

