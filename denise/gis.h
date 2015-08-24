//
// gis.h
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

#ifndef DENISE_GIS_H
#define DENISE_GIS_H

#include <stdint.h>
#include <denise/basics.h>
#include <denise/geodesy.h>
#include <denise/geometry.h>
#include <denise/graphics.h>

using namespace std;
using namespace Cairo;

namespace denise
{

   class Gshhs : public Polygon,
                 public Geodetic_Cairoable
   {

      protected:

         class Point
         {

            public:

               int32_t
               x;

               int32_t
               y;

               void
               swap_endian ();

               void
               read (FILE* file);

               void
               read (ifstream& file);

         };

         class Header
         {

            protected:

               Integer
               buffer_size;

               uint8_t*
               buffer;

               virtual void
               swap_endian ();

               uint8_t
               get_version () const;

               const int32_t&
               get_flag () const;

            public:

               Header ();

               Header (const Header& header);

               ~Header ();

               void
               read (FILE* file);

               void
               read (ifstream& file);

               void
               reset (Polygon::Vertex* handle_ptr);

               const int32_t&
               get_int32_t (const Integer position) const;

               int32_t&
               get_int32_t (const Integer position);

               const int16_t&
               get_int16_t (const Integer position) const;

               int16_t&
               get_int16_t (const Integer position);

               const int8_t&
               get_int8_t (const Integer position) const;

               int8_t&
               get_int8_t (const Integer position);

               void
               write (FILE* file,
                      const bool preserve_header);

               virtual int8_t
               get_level () const; 

               virtual bool
               is_cross_greenwich () const; 

               virtual bool
               is_western () const; 

         };

         virtual Header*
         get_header_ptr () const;

         void
         add_polygon (FILE* file,
                      const Header& header,
                      const Real max_longitude);

         void
         add_polygon (ifstream& file,
                      const Header& header,
                      const Real max_longitude);

         void
         write_simple_polygon (FILE* file,
                               const Polygon::Vertex* vertex_ptr) const;

      public:

         Gshhs ();

         Gshhs (const Dstring& file_path,
                const Integer level = 1);

         ~Gshhs ();

         void
         save (const Dstring& file_path) const;

         Gshhs*
         clip (const Polygon& polygon);

         Gshhs*
         clip (const Domain_1D& domain_latitude,
               const Domain_1D& domain_longitude);

         virtual void
         cairo (const RefPtr<Context> cr,
                const Geodetic_Transform& transform) const;

         virtual void
         cairo (const RefPtr<Context> cr,
                const Geodetic_Transform& transform,
                const Polygon& clip_polygon);

   };

   class Gshhs2 : public Gshhs
   {

      protected:

         class Header : public Gshhs::Header
         {

            public:

               int32_t
               id;

               int32_t
               n;

               int32_t
               flag;

               int32_t
               west;

               int32_t
               east;

               int32_t
               south;

               int32_t
               north;

               int32_t
               area;

               int32_t
               area_full;

               int16_t
               container;

               int16_t
               ancestor;

               void
               swap_endian ();

         };

   };

   class Blue_Marble : public Raster
   {

      private:

         static Color
         get_color (const Lat_Long& lat_long,
                    const Real start_latitude,
                    const Real start_longitude,
                    FILE* west_file,
                    FILE* east_file);

         void
         fill_raster (const Dstring& blue_marble_path,
                      const Transform_2D& transform_2d);

      public:

         class Exception : public denise::Exception
         {

            public:

               Exception (const Dstring& description = "");

         };

         Blue_Marble (const Dstring& blue_marble_path,
                      const Transform_2D& transform,
                      const Size_2D& size_2d);

         ~Blue_Marble ();

   };

   class Gtopo30 : public Raster
   {

      private:

         enum Tile
         {
            W180S60 = 0,
            W120S60 = 1,
            W060S60 = 2,
            W000S60 = 3,
            E060S60 = 4,
            E120S60 = 5,
            W180S10 = 6,
            W140S10 = 7,
            W100S10 = 8,
            W060S10 = 9,
            W020S10 = 10,
            E020S10 = 11,
            E060S10 = 12,
            E100S10 = 13,
            E140S10 = 14,
            W180N40 = 15,
            W140N40 = 16,
            W100N40 = 17,
            W060N40 = 18,
            W020N40 = 19,
            E020N40 = 20,
            E060N40 = 21,
            E100N40 = 22,
            E140N40 = 23,
            W180N90 = 24,
            W140N90 = 25,
            W100N90 = 26,
            W060N90 = 27,
            W020N90 = 28,
            E020N90 = 29,
            E060N90 = 30,
            E100N90 = 31,
            E140N90 = 32,
            NUMBER_OF_TILES = 33
         };

         const Real
         delta;

         const Real
         start_latitude;

         const Real
         start_longitude;

         const Dstring
         gtopo30_path;

         FILE**
         file_array;

         Size_2D*
         tile_size_array;

         Index_2D*
         tile_index_array;

         static Tile
         get_tile (const Integer gtopo30_i,
                   const Integer gtopo30_j);

         static Dstring
         get_tile_string (const Tile tile);

         static Size_2D
         get_tile_size (const Tile tile);

         static Index_2D
         get_tile_index (const Tile tile);

         static Domain_2D
         get_tile_domain (const Tile tile);

         static Domain_1D
         get_tile_domain_latitude (const Tile tile);

         static Domain_1D
         get_tile_domain_longitude (const Tile tile);

         Real*
         get_data (const Transform_2D& transform);

         void
         fill_raster (const Real* data,
                      const Blue_Marble* blue_marble_ptr,
                      const Real scale_factor);

      public:

         class Exception : public denise::Exception
         {

            public:

               Exception (const Dstring& description = "");

         };

         Gtopo30 (const Dstring& gtopo30_path,
                  const Transform_2D& transform,
                  const Size_2D& size_2d,
                  const Dstring& blue_marble_path = "",
                  const Real scale_factor = 0.0025);

         ~Gtopo30 ();

         Real
         get_datum (const Lat_Long& lat_long) const;

   };

   class Land_Mask
   {

      private:

         uint8_t*
         buffer;

         uint16_t
         n;

         int16_t
         start_i;

         int16_t
         start_j;

         uint16_t
         ni;

         uint16_t
         nj;

         uint16_t
         get_degree_box_size () const;

         uint8_t
         get_byte_mask (const uint8_t bj) const;

         void
         acquire_lat_long (Real& latitude,
                           Real& longitude,
                           const int16_t i,
                           const int16_t j,
                           const uint8_t ii,
                           const uint8_t jj) const;

         void
         acquire_indices (int16_t& i,
                          int16_t& j,
                          uint8_t& ii,
                          uint8_t& jj,
                          const Lat_Long& lat_long) const;

      public:

         Land_Mask (const Polygon& polygon,
                    const uint16_t n = 20,
                    const int16_t start_i = -90,
                    const int16_t start_j = -180,
                    const uint16_t ni = 180,
                    const uint16_t nj = 360);

         Land_Mask (const Gtopo30& gtopo30,
                    const uint16_t n = 20,
                    const int16_t start_i = -90,
                    const int16_t start_j = -180,
                    const uint16_t ni = 180,
                    const uint16_t nj = 360);

         Land_Mask (const Dstring& land_mask_file_path);

         ~Land_Mask ();

         void
         populate (const Polygon& polygon,
                   const int16_t start_i,
                   const int16_t start_j,
                   const uint16_t ni,
                   const uint16_t nj);

         void
         populate (const Gtopo30& gtopo30,
                   const int16_t start_i,
                   const int16_t start_j,
                   const uint16_t ni,
                   const uint16_t nj);

         void
         save (const Dstring& land_mask_file_path) const;

         bool
         is_land (const Lat_Long& lat_long) const;

         Real
         get_fetch (const Lat_Long& lat_long,
                    const Real bearing) const;

         Real
         get_effective_fetch (const Lat_Long& lat_long,
                              const Real bearing) const;

   };

   class Fetch
   {

      private:

         uint8_t*
         buffer;

         uint16_t
         nb;

         uint16_t
         n;

         int16_t
         start_i;

         int16_t
         start_j;

         uint16_t
         ni;

         uint16_t
         nj;

         uint8_t&
         get_byte (const Lat_Long& lat_long,
                   const Real bearing);

         void
         populate (const Land_Mask& land_mask);

      public:

         Fetch (const Gshhs& gshhs,
                const uint16_t nb = 36,
                const uint16_t n = 20,
                const int16_t start_i = -90,
                const int16_t start_j = -180,
                const uint16_t ni = 180,
                const uint16_t nj = 360);

         Fetch (const Gtopo30& gtopo30,
                const uint16_t nb = 36,
                const uint16_t n = 20,
                const int16_t start_i = -90,
                const int16_t start_j = -180,
                const uint16_t ni = 180,
                const uint16_t nj = 360);

         Fetch (const Land_Mask& land_mask,
                const uint16_t nb = 36,
                const uint16_t n = 20,
                const int16_t start_i = -90,
                const int16_t start_j = -180,
                const uint16_t ni = 180,
                const uint16_t nj = 360);

         Fetch (const Dstring& fetch_file_path);

         ~Fetch ();

         void
         save (const Dstring& fetch_file_path) const;

         Real
         get_fetch (const Real& bearing,
                    const Lat_Long& lat_long) const;

   };

   class Nsd
   {

      private:

         static Real
         get_lat_long_number (const Dstring& lat_long_string);

         void
         wmo_init (const Dstring& str);

         void
         icao_init (const Dstring& str);

      public:

         Integer
         wmo_id;

         Dstring
         icao_id;

         Dstring
         name;

         Dstring
         country;

         Lat_Long
         lat_long;

         Real
         elevation;

         Nsd (const Dstring& nsd_string,
              const bool icao);

   };

   class Nsd_Wmo : private map<Integer, Nsd>
   {

      public:

         Nsd_Wmo (const Dstring& file_path);

         const Nsd&
         get_nsd (const Integer wmo_id) const;

   };

   class Nsd_Icao : private map<Dstring, Nsd>
   {

      public:

         Nsd_Icao (const Dstring& file_path);

         const Nsd&
         get_nsd (const Dstring& icao_id) const;

   };

   class Nsd_Exception : public Exception
   {

      public:

         Nsd_Exception (const Dstring& description = "");

   };

   class Outl : public Polyline,
                public Geodetic_Cairoable
   {

      private:

         void
         read (const Dstring& file_path);

      public:

         Outl (const Dstring& file_path);

         virtual void
         cairo (const RefPtr<Context> cr,
                const Geodetic_Transform& transform) const;

   };

   class Kidney : public map<Real, Polygon*>,
                  public Geodetic_Cairoable
   {

      private:

         void
         read (const Dstring& file_path);

      public:

         Kidney (const Dstring& file_path);

         ~Kidney ();

         virtual void
         cairo (const RefPtr<Context> cr,
                const Geodetic_Transform& transform) const;

   };

}

#endif /* DENISE_GIS_H */ 
