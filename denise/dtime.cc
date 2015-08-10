//
// dtime.cc
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

#include "dstring.h"
#include "dtime.h"
#include "util.h"

using namespace std;
using namespace denise;

void
Dtime::init (const Dstring& time_string,
             bool is_local)
{

   if (time_string == L"big_bang") { this->t = GSL_NEGINF; return; }
   if (time_string == L"big_crunch") { this->t = GSL_POSINF; return; }
   if (time_string == L"nat") { this->t = GSL_NAN; return; }

   Integer year  = stoi (time_string.substr (0, 4));
   Integer month = stoi (time_string.substr (4, 2));
   Integer day   = stoi (time_string.substr (6, 2));

   Integer n = time_string.size ();
   Integer hour = 0, minute = 0, second = 0;

   if (n >= 10)
   {
      hour = stoi (time_string.substr (8, 2));
      if (n >= 12)
      {
         minute = stoi (time_string.substr (10, 2));
         if (n >= 14)
         {
            second = stoi (time_string.substr (12, 2));
         }
      }
   }

   init (year, month, day, hour, minute, second, is_local);

}

void
Dtime::init (const Dstring& time_string,
             const Dstring& format,
             bool is_local)
{

   if (time_string == L"big_bang") { this->t = GSL_NEGINF; return; }
   if (time_string == L"big_crunch") { this->t = GSL_POSINF; return; }
   if (time_string == L"nat") { this->t = GSL_NAN; return; }

   struct tm tm_struct;
   tm_struct.tm_min = 0; tm_struct.tm_sec = 0; tm_struct.tm_sec = 0;
   tm_struct.tm_min = 0; tm_struct.tm_hour = 0; tm_struct.tm_mday = 0;
   tm_struct.tm_mon = 0; tm_struct.tm_year = 0; tm_struct.tm_wday = 0;
   tm_struct.tm_yday = 0; tm_struct.tm_isdst = 0;

   const string& ts = time_string.get_string ();
   const string& fmt = format.get_string ();
   strptime (ts.c_str (), fmt.c_str (), &tm_struct);

   if (tm_struct.tm_min < 0 || tm_struct.tm_min > 59) { tm_struct.tm_min = 0; }
   if (tm_struct.tm_sec < 0 || tm_struct.tm_sec > 61) { tm_struct.tm_sec = 0; }

   init (&tm_struct, is_local);

}

void
Dtime::init (const struct tm* time_struct,
             bool is_local)
{

   Integer year = time_struct->tm_year + 1900;
   Integer month = time_struct->tm_mon + 1;
   Integer day   = time_struct->tm_mday;
   Integer hour  = time_struct->tm_hour;
   Integer minute = time_struct->tm_min;
   Integer second = time_struct->tm_sec;

   return init (year, month, day, hour, minute, second, is_local);

}

void
Dtime::init (Integer year,
             Integer month,
             Integer day,
             Integer hour,
             Integer minute,
             Integer second,
             bool is_local)
{

   tzset ();
   const Integer epoch_day = 2440588;
   Integer jd = get_julian_day (year, month, day);

   t = double (jd - epoch_day) * 24 + double (hour);
   t += double (minute * 60 + second) / 3600;
//   t = (double (jd - epoch_day) * 24 + double (hour)) * 3600;
//   t += double (minute * 60 + second);

   if (is_local)
   {
#ifndef HAVE_TIMEZONE
      t += double (timezone) / 3600;
//      t += double (timezone);
#else
      time_t now = time (NULL);
      struct tm ts = *localtime (&now);
      t -= double (ts.tm_gmtoff) / 3600;
//      t -= double (ts.tm_gmtoff);
#endif
   }

}

struct tm
Dtime::get_time_struct (bool is_local) const
{

   tzset ();
   struct tm time_struct;
   const Integer epoch_day = 2440588;

   double t = this->t;

   if (is_local)
   {
#ifndef HAVE_TIMEZONE
      t -= double (timezone / 3600);
//      t -= double (timezone);
#else
      time_t now = time (NULL);
      struct tm ts = *localtime (&now);
      t += double (ts.tm_gmtoff) / 3600;
//      t += double (ts.tm_gmtoff);
#endif
   }

   double ut = fmod (t, 24);
   while (ut < 0) { ut += 24; }
   Integer jd = Integer (rint (t - ut)) / 24 + epoch_day;
//   double ut = fmod (t, 86400);
//   while (ut < 0) { ut += 86400; }
//   Integer jd = Integer (rint (t - ut)) / 86400 + epoch_day;

   Integer a = jd + 32044;
   Integer b = (4*a + 3) / 146097;
   Integer c = a - (146097*b) / 4;
   Integer d = (4*c + 3) / 1461;
   Integer e = c - (1461*d) / 4;
   Integer m = (5*e + 2) / 153;

   Integer day = e - (153*m + 2)/5 + 1;
   Integer month = m + 3 - 12 * (m/10);
   Integer year = 100 * b + d - 4800 + (m/10);

   Integer hour = Integer (floor (ut));
//   Integer hour = Integer (floor (ut / 3600));

   ut -= hour;
//   ut -= (hour * 3600);

   Integer second = Integer (rint (ut * 3600));
//   Integer second = Integer (rint (ut));
   Integer minute = second / 60;
   second -= minute * 60;

   time_struct.tm_year = year - 1900;
   time_struct.tm_mon  = month - 1;
   time_struct.tm_mday = day;
   time_struct.tm_hour = hour;
   time_struct.tm_min  = minute;
   time_struct.tm_sec  = second;

   {
      Integer a = (14 - month) / 12;
      Integer y = year - a;
      Integer m = month + 12*a - 2;
      time_struct.tm_wday = (day + y + y/4 - y/100 + y/400 + (31*m)/12) % 7;
   }

   {
      Integer jd_0 = get_julian_day (year, 1, 1);
      time_struct.tm_yday = jd - jd_0;
   }

   time_struct.tm_isdst = false;

   return time_struct;

}

Dtime::Dtime ()
         : t (double (time (NULL)) / 3600)
{
}

Dtime::Dtime (const double t,
              const bool snap_to_minute)
         : t (snap_to_minute ? round (t*60) / 60 : t)
{
}

Dtime::Dtime (const Dtime& time)
         : t (time.t)
{
}

Dtime::Dtime (const Dstring& time_string,
              const bool is_local)
{
   init (time_string, is_local);
}

Dtime::Dtime (const Dstring& time_string,
              const Dstring& format,
              const bool is_local)
{
   if (format.empty ()) { init (time_string, is_local); }
   else { init (time_string, format, is_local); }
}

Dtime::Dtime (const Integer year,
              const Integer month,
              const Integer day,
              const Integer hour,
              const Integer minute,
              const Integer second,
              const bool is_local)
{
   init (year, month, day, hour, minute, second, is_local);
}

Dtime
Dtime::now ()
{
   return Dtime ();
}

Dtime
Dtime::last_sub_synoptic ()
{
   const Dtime now;
   return Dtime (now.t - modulo (now.t, 3));
}

Dtime
Dtime::last_synoptic ()
{
   const Dtime now;
   return Dtime (now.t - modulo (now.t, 6));
}

Dtime
Dtime::last_synoptic_12 ()
{
   const Dtime now;
   return Dtime (now.t - modulo (now.t, 12));
}

Dtime
Dtime::hours_later (const double hours)
{
   return Dtime (Dtime::now ().t + hours);
}

Dstring
Dtime::get_string (const Dstring& format,
                   const bool is_local) const
{

   int inf;

   if (gsl_isnan (t))
   {
      return Dstring (L"nat");
   }
   else
   if ((inf = gsl_isinf (t)) != 0)
   {
      const bool future = (inf > 0);
      return Dstring (future ? L"big_crunch" : L"big_bang");
   }
   else
   {

      struct tm time_struct = get_time_struct (is_local);
      const string fmt (format.begin (), format.end ());

      Integer n = format.size () * 5;
      char* c_time_string = new char[n];
      strftime (c_time_string, n, fmt.c_str (), &time_struct);

      const Dstring str (c_time_string);
      delete[] c_time_string;
      return str;

   }

}

Integer
Dtime::get_year () const
{
   return stoi (get_string (L"%Y"));
}

Integer
Dtime::get_month () const
{
   return stoi (get_string (L"%m"));
}

Integer
Dtime::get_day () const
{
   return stoi (get_string (L"%d"));
}

Integer 
Dtime::get_hour () const
{
   return stoi (get_string (L"%H"));
}
         
Integer
Dtime::get_minute () const
{
   return stoi (get_string (L"%M"));
}
         
Integer
Dtime::get_second () const
{
   return stoi (get_string (L"%S"));
}
         
Integer
Dtime::get_day_of_year () const
{
   return stoi (get_string (L"%j"));
}

Dtime
Dtime::get_start_of_year () const
{
   return Dtime (get_string (L"%Y") + L"010100");
}

Dtime
Dtime::get_start_of_month () const
{
   return Dtime (get_string (L"%Y%m") + L"0100");
}

bool
Dtime::is_nat () const
{
   return gsl_isnan (t);
}

Dtime
Dtime::get_current_time ()
{
   return Dtime (double (time (NULL)) / 3600);
//   return Dtime (double (time (NULL)));
}
         
bool
Dtime::is_leap (Integer year)
{
   if (year % 400 == 0)      { return true; }
   else if (year % 100 == 0) { return false; }
   else if (year % 4 == 0)   { return true; }
   else                      { return false; }
}

Integer
Dtime::get_julian_day (Integer year,
                       Integer month,
                       Integer day)
{

   Integer a = (14 - month) / 12;
   Integer y = year + 4800 - a;
   Integer m = month + 12 * a - 3;

   return day + (153*m + 2) / 5 + 365 * y + (y/4) - (y/100) + (y/400) - 32045;

}

Tuple
Dtime::get_time_tuple (const Dstring& time_string,
                       const Dstring& delimiter)
{

   Tuple tuple;
   const Tokens tokens (time_string, delimiter);

   for (Tokens::const_iterator iterator = tokens.begin ();
        iterator != tokens.end (); iterator++)
   {
      const Dstring& token = *(iterator);
      const Dtime dtime (token);
      tuple.push_back (dtime.t);
   }

   return tuple;

}

Tuple
Dtime::get_yearly_time_tuple (const Dtime& start_time,
                              const Dtime& end_time)
{

   Tuple tuple;
   Integer year = start_time.get_year ();

   for (Dtime dtime = Dtime (year, 1, 1);
        dtime <= end_time;
        dtime = Dtime (++year, 1, 1))
   {
      if (dtime < start_time) { continue; }
      tuple.push_back (dtime.t);
   }

}

Tuple
Dtime::get_monthly_time_tuple (const Dtime& start_time,
                               const Dtime& end_time)
{

   Tuple tuple;
   Integer start_year = start_time.get_year ();
   Integer end_year = end_time.get_year ();

   for (Integer year = start_year; year <= end_year; year++)
   {
      for (Integer month = 1; month <= 12; month++)
      {
         const Dtime dtime (year, month, 1);
         if (dtime < start_time) { continue; }
         if (dtime > end_time) { continue; }
         tuple.push_back (dtime.t);
      }
   }

   return tuple;

}

Tuple
Dtime::get_trimonthly_time_tuple (const Dtime& start_time,
                                  const Dtime& end_time)
{

   Tuple tuple;
   Integer start_year = start_time.get_year ();
   Integer end_year = end_time.get_year ();

   for (Integer year = start_year; year <= end_year; year++)
   {
      for (Integer month = 1; month <= 12; month++)
      {
         for (Integer day = 1; day < 25; day += 10)
         {
            const Dtime dtime (year, month, day);
            if (dtime < start_time) { continue; }
            if (dtime > end_time) { continue; }
            tuple.push_back (dtime.t);
         }
      }
   }

   return tuple;

}

bool
Dtime::operator == (const Dtime& time) const
{
   return fabs (time.t - t) < 0.01;
}

bool
Dtime::operator != (const Dtime& time) const
{
   return fabs (time.t - t) >= 0.01;
}

bool
Dtime::operator < (const Dtime& time) const
{
   return t < time.t;
}

bool
Dtime::operator > (const Dtime& time) const
{
   return t > time.t;
}

bool
Dtime::operator <= (const Dtime& time) const
{
   return (t - time.t) < 0.01;
}

bool
Dtime::operator >= (const Dtime& time) const
{
   return (t - time.t - t) > -0.01;
}

Dtime::Span::Span ()
   : start_t (GSL_NEGINF),
     end_t (GSL_POSINF)
{
   this->start_t = GSL_NEGINF;
   this->end_t = GSL_POSINF;
}

Dtime::Span::Span (const Dstring& str)
   : start_t (GSL_NEGINF),
     end_t (GSL_POSINF)
{
   const Tokens tokens (str, L"-");
   if (tokens.size () == 1)
   {
      const Real t = Dtime (tokens[0]).t;
      const bool dash_at_front = (str[0] == '-');
      const bool dash_at_back = (str[str.length () - 1] == '-');
           if (dash_at_front) { end_t = t; }
      else if (dash_at_back)  { start_t = t; }
      else                    { start_t = t; end_t = t; }
   }
   else
   if (tokens.size () == 2)
   { 
      start_t = Dtime (tokens[0]).t;
      end_t = Dtime (tokens[1]).t;
   }
}

Dtime::Span::Span (const Dtime& start,
                   const Dtime& end)
   : start_t (start.t),
     end_t (end.t)
{
}

Dtime
Dtime::Span::get_start (const bool snap_to_minute) const
{
   return Dtime (start_t, snap_to_minute);
}

Dtime
Dtime::Span::get_end (const bool snap_to_minute) const
{
   return Dtime (end_t, snap_to_minute);
}

bool
Dtime::Span::match (const Dtime& dtime) const
{
   return (dtime.t >= start_t) && (dtime.t <= end_t);
}

bool
Dtime::Span::operator == (const Dtime::Span& span) const
{
   return (start_t == span.start_t) && (end_t == span.end_t);
}

bool
Dtime::Span::operator != (const Dtime::Span& span) const
{
   return (start_t != span.start_t) || (end_t != span.end_t);
}

bool
Dtime::Span::operator < (const Dtime::Span& span) const
{
   if (start_t == span.start_t) { return end_t < span.end_t; }
   else                         { return start_t < span.start_t; }
}

bool
Dtime::Span::operator > (const Dtime::Span& span) const
{
   if (start_t == span.start_t) { return end_t > span.end_t; }
   else                         { return start_t > span.start_t; }
}

bool
Dtime::Set::match (const Dtime& dtime) const
{
   if (size () == 0) { return true; }
   for (auto iterator = begin (); iterator != end (); iterator++)
   {
      const Dtime::Span& span = *(iterator);
      if (span.match (dtime)) { return true; }
   }
   return false;
}

Dtime::Set::Set (const Dstring& str)
{

   const Tokens tokens (str, L":");

   for (auto iterator = tokens.begin (); iterator != tokens.end (); iterator++)
   {
      const Dstring& token = *(iterator);
      insert (Dtime::Span (token));
   }

}

Dtime::Set::Set (const Dtime& start,
                 const Dtime& end)
{
   insert (Dtime::Span (start, end));
}

namespace denise
{

   wostream&
   operator << (wostream &out,
                const Dtime& time)
   {
      out << time.get_string ();
      return out;
   }

   wostream&
   operator << (wostream &out,
                const Dtime::Span& span)
   {
      const Dtime& start = span.get_start ();
      const Dtime& end = span.get_end ();
      if (start == end) { out << start; }
      else { out << start << "-" << end; }
      return out;
   }

   wostream&
   operator << (wostream &out,
                const Dtime::Set& set)
   {
      for (auto iterator = set.begin (); iterator != set.end (); iterator++)
      {
         if (iterator != set.begin ()) { out << ":"; }
         out << *(iterator);
      }
      return out;
   }

}

