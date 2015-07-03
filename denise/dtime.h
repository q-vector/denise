//
// dtime.h
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

#ifndef DENISE_DTIME_H
#define DENISE_DTIME_H

#include <set>
#include <ctime>
#include <denise/basics.h>
#include <denise/exception.h>
#include <denise/dstring.h>

using namespace std;

namespace denise
{

   /// Represent time.
   ///
   /// The internal representation of Dtime is double t which
   /// represents numbers of hours since epoch (January 1, 1970
   /// 00UTC).
   class Dtime
   {

      private:

         void
         init (const string& time_string,
               bool is_local);

         void
         init (const string& time_string,
               const string& format,
               bool is_local);

         void
         init (const struct tm* time_struct,
               bool is_local = false);

         void
         init (Integer year,
               Integer month,
               Integer day,
               Integer hour = 0,
               Integer minute = 0,
               Integer second = 0,
               bool is_local = false);

         struct tm
         get_time_struct (bool is_local = false) const;

      public:

         double
         t;

         Dtime ();

         /// Constructor that accepts t
         Dtime (const double t,
                const bool snap_to_minute = true);

         Dtime (const Dtime& time);

         /// Constructor that accepts a string representation of time
         /// of the format %Y%m%d%H or %Y%m%d%H%M or %Y%m%d%H%M%S.
         Dtime (const string& time_string,
                const bool is_local = false);

         /// Constructor that accepts a string representation of time
         /// of the given format %Y%m%d%H or %Y%m%d%H%M or %Y%m%d%H%M%S.
         Dtime (const string& time_string,
                const string& format,
                const bool is_local = false);

         /// Constructor that accepts integral values of year
         /// month, day, hour, minute and second.
         Dtime (const Integer year,
                const Integer month,
                const Integer day,
                const Integer hour = 0,
                const Integer minute = 0,
                const Integer second = 0,
                bool is_local = false);

         static Dtime
         now ();

         static Dtime
         last_sub_synoptic ();

         static Dtime
         last_synoptic ();

         static Dtime
         last_synoptic_12 ();

         static Dtime
         hours_later (const double hours);

         /// Returns a string representation of the time with
         /// the given format.
         string
         get_string (const string& format = string ("%Y%m%d%H"),
                     const bool is_local = false) const;

         /// Returns the integral value of year.
         Integer
         get_year () const;

         /// Returns the integral value of month.
         Integer
         get_month () const;

         /// Returns the integral value of day of month.
         Integer
         get_day () const;

         /// Returns the integral value of hour.
         Integer
         get_hour () const;

         /// Returns the integral value of minute.
         Integer
         get_minute () const;

         /// Returns the integral value of second.
         Integer
         get_second () const;

         /// Returns the integral value of day of year.
         Integer
         get_day_of_year () const;

         Dtime
         get_start_of_year () const;

         Dtime
         get_start_of_month () const;

         /// Returns true if internal value is NAN.
         bool
         is_nat () const;

         /// Returns current Dtime.
         static Dtime
         get_current_time ();

         /// Returns true if given year is a leap year.
         static bool
         is_leap (Integer year);

         /// Returns the Julian day of the given year, month and day.
         static Integer
         get_julian_day (Integer year,
                         Integer month,
                         Integer day);

         static Tuple
         get_time_tuple (const string& time_str,
                         const string& delimiter = string (":"));

         static Tuple
         get_yearly_time_tuple (const Dtime& start_time,
                                const Dtime& end_time);

         static Tuple
         get_monthly_time_tuple (const Dtime& start_time,
                                 const Dtime& end_time);

         static Tuple
         get_trimonthly_time_tuple (const Dtime& start_time,
                                    const Dtime& end_time);

         bool
         operator == (const Dtime& time) const;

         bool
         operator != (const Dtime& time) const;

         bool
         operator < (const Dtime& time) const;

         bool
         operator > (const Dtime& time) const;

         bool
         operator <= (const Dtime& time) const;

         bool
         operator >= (const Dtime& time) const;

         class Span
         {

            public:

               Real
               start_t;

               Real
               end_t;

            public:

               Span ();

               Span (const string& str);

               Span (const Dtime& start,
                     const Dtime& end);

               Dtime
               get_start (const bool snap_to_minute = false) const; 

               Dtime
               get_end (const bool snap_to_minute = false) const; 

               bool
               match (const Dtime& dtime) const;

               bool
               operator == (const Span& span) const;

               bool
               operator != (const Span& span) const;

               bool
               operator < (const Span& span) const;

               bool
               operator > (const Span& span) const;

         };

         class Set : public set<Span>
         {

            public: 

               Set (const string& str);

               Set (const Dtime& start,
                    const Dtime& end);

               virtual bool
               match (const Dtime& dtime) const;

         };

   };

   ostream&
   operator << (ostream& out,
                const Dtime& time);

   ostream&
   operator << (ostream& out,
                const Dtime::Span& span);

   ostream&
   operator << (ostream& out,
                const Dtime::Set& set);

}

#endif /* DENISE_DTIME_H */

