//
// gts.h
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

#ifndef DENISE_GTS_H
#define DENISE_GTS_H

#include <map>
#include <denise/geodesy.h>
#include <denise/thermo.h>

using namespace std;

namespace denise
{

   class Gts : public Tokens
   {

      protected:

         const Dtime
         time_hint;

         Dstring
         yygg;

         Integer
         wmo_id;

      public:

         Gts (const Dstring& message,
              const Dtime& time_hint = Dtime::now ());

         Dstring
         get_yygg () const;

         Dstring
         get_time_string () const;

         const Dstring& 
         get_wmo_id () const;

         Dstring
         get_key () const;

         static Tokens*
         parse_file (const Dstring& file_path);

   };

   class Gts_Ptr_Store : public map<Dstring, Gts*>
   {

      public:

         ~Gts_Ptr_Store ();

         void
         ingest (const Dstring& message,
                 const Dtime& time_hint);

         const Gts&
         get_gts (const Dstring& gts_key) const;

   };

   class Gts_Ptr_Grandstore : public map<Dstring, Gts_Ptr_Store>
   {

      private:

         void
         ingest (const Dstring& genre,
                 const Dstring& message,
                 const Dtime& time_hint);

      public:

         void
         ingest (const Dstring& message,
                 const Dtime& time_hint);

         void
         ingest_file (const Dstring& file_path);

         void
         ingest_dir (const Dstring& dir_path,
                     const Reg_Exp& file_reg_exp);

         const Gts_Ptr_Store&
         get_gts_ptr_store (const Dstring& genre) const;

         const Gts&
         get_gts (const Dstring& genre,
                  const Dstring& gts_key) const;

   };

   class Temp : public Gts
   {

      protected:

         pair<Real, Real>
         parse_tttdd (const Dstring& tttdd) const;

         Wind
         parse_ddfff (const Dstring& ddfff) const;

         bool
         use_knots () const;

      public:

         Temp (const Dstring& message,
               const Dtime& time_hint);

         virtual void
         parse_to (Sounding& sounding) const = 0;

   };

   class Ttac : public Temp
   {

      protected:

         const wchar_t
         get_highest_wind_code () const;

         Integer
         parse_special_level (Sounding& sounding,
                              Integer& index) const;

         Integer
         parse_standard_levels (Sounding& sounding,
                                Integer& index) const;

         pair<Real, Real>
         parse_pphhh (const Dstring& pphhh) const;

         virtual void
         interpret_p_z (const Integer pp,
                        Real& pressure,
                        Real& geopotential_height) const = 0;

         virtual Integer
         get_number_of_standard_levels () const = 0;

      public:

         Ttac (const Dstring& message,
               const Dtime& time_hint);

   };

   class Ttbd : public Temp
   {

      protected:

         virtual Real
         parse_nnppp (const Dstring& nnppp) const = 0;

         bool
         parse_significant_temperature_levels (Sounding& sounding,
                                               Integer& index) const;

         bool
         parse_significant_wind_levels (Sounding& sounding,
                                        Integer& index) const;

      public:

         Ttbd (const Dstring& message,
               const Dtime& time_hint);

   };

   class Ttaa : public Ttac
   {

      protected:

         void
         interpret_p_z (const Integer pp,
                        Real& pressure,
                        Real& geopotential_height) const;

         Integer
         get_number_of_standard_levels () const;

      public:

         Ttaa (const Dstring& message,
               const Dtime& time_hint);

         void
         parse_to (Sounding& sounding) const;

   };

   class Ttbb : public Ttbd
   {

      protected:

         Real
         parse_nnppp (const Dstring& nnppp) const;

      public:

         Ttbb (const Dstring& message,
               const Dtime& time_hint);

         void
         parse_to (Sounding& sounding) const;

   };

   class Ttcc : public Ttac
   {

      protected:

         void
         interpret_p_z (const Integer pp,
                        Real& pressure,
                        Real& geopotential_height) const;

         Integer
         get_number_of_standard_levels () const;

      public:

         Ttcc (const Dstring& message,
               const Dtime& time_hint);

         void
         parse_to (Sounding& sounding) const;

   };

   class Ttdd : public Ttbd
   {

      protected:

         Real
         parse_nnppp (const Dstring& nnppp) const;

      public:

         Ttdd (const Dstring& message,
               const Dtime& time_hint);

         void
         parse_to (Sounding& sounding) const;

   };

   class Temp_Sounding : public Sounding
   {
   };

   class Temp_Sounding_Store : public map<Dstring, Temp_Sounding>
   {
   };

   ostream&
   operator << (ostream &out_file,
                const Gts& gts);

};

#endif /* DENISE_GTS_H */ 
