//
// gis.cc
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

#include <fstream>
#include "gis.h"

using namespace std;
using namespace denise;

void
Gshhs::Point::swap_endian ()
{
   denise::swap_endian (&x, sizeof (int32_t));
   denise::swap_endian (&y, sizeof (int32_t));
}

void
Gshhs::Point::read (FILE* file)
{

   fread (this, sizeof (class Point), 1, file);

#ifndef WORDS_BIGENDIAN
   swap_endian ();
#endif

   if (x % 100000 == 0) { x += 1; }
   if (y % 100000 == 0) { y += 1; }
   
}

void
Gshhs::Header::swap_endian ()
{

   const Integer size_4 = sizeof (int32_t);
   const Integer size_2 = sizeof (int16_t);

   // swapped or unswapped, old or new version, 10th byte must = 0
   const bool old = (buffer[10] == 0);

   denise::swap_endian (buffer, size_4);  // id
   denise::swap_endian (buffer + 4, size_4);  // n
   denise::swap_endian (buffer + 8, size_4);  // level / flag
   denise::swap_endian (buffer + 12, size_4); // west
   denise::swap_endian (buffer + 16, size_4); // east
   denise::swap_endian (buffer + 20, size_4); // south
   denise::swap_endian (buffer + 24, size_4); // north
   denise::swap_endian (buffer + 28, size_4); // area
   denise::swap_endian (buffer + 32, size_4); // version / area_full

   if (old)
   {
      denise::swap_endian (buffer + 36, size_2); // greenwich
      denise::swap_endian (buffer + 38, size_2); // source
   }
   else
   {
      denise::swap_endian (buffer + 36, size_4); // container
      denise::swap_endian (buffer + 40, size_4); // ancestor
   }

}

uint8_t
Gshhs::Header::get_version () const
{
   return (get_flag () >> 8) & 0xff;
}

const int32_t&
Gshhs::Header::get_flag () const
{
   return get_int32_t (8);
}

Gshhs::Header::Header ()
   : buffer (NULL)
{
}

Gshhs::Header::Header (const Header& header)
{
   const bool is_empty = (header.buffer == NULL);
   if (is_empty) { this->buffer = NULL; }
   else
   {
      buffer_size = header.buffer_size;
      this->buffer = new uint8_t[buffer_size];
      memcpy (this->buffer, header.buffer, buffer_size);
   }
}

Gshhs::Header::~Header ()
{
   if (buffer != NULL)
   {
      delete[] buffer;
   }
}

void
Gshhs::Header::read (FILE* file)
{

   const size_t old_size = 40;
   const size_t new_size = 44;
   const size_t size_difference = new_size - old_size;

   uint8_t* temp_buffer = new uint8_t[old_size];

   if (fread (temp_buffer, old_size, 1, file) == 0)
   {
      throw IO_Exception ("EOF reached");
   }

   const bool old = (temp_buffer[10] == 0);
   buffer_size = (old ? old_size : new_size);
   buffer = new uint8_t[buffer_size];

   memcpy (buffer, temp_buffer, old_size);
   delete[] temp_buffer;

   if (!old)
   {
      if (fread (buffer + old_size, size_difference, 1, file) == 0)
      {
         throw IO_Exception ("EOF reached");
      }
   }

#ifndef WORDS_BIGENDIAN
   swap_endian ();
#endif
   
}

void
Gshhs::Header::reset (Polygon::Vertex* handle_ptr)
{

   Real start_latitude  = GSL_POSINF;
   Real end_latitude    = GSL_NEGINF;
   Real start_longitude = GSL_POSINF;
   Real end_longitude   = GSL_NEGINF;

   Integer n = handle_ptr->n;
   Polygon::Vertex* current_ptr = (Polygon::Vertex*)handle_ptr;

   for (Integer i = 0; i < n; i++)
   {

      const Lat_Long& lat_long = (const Lat_Long&)(*current_ptr);
      const Real& latitude = lat_long.latitude;
      const Real& longitude = lat_long.longitude;

      if (latitude < start_latitude) { start_latitude = latitude; }
      if (latitude > end_latitude) { end_latitude = latitude; }
      if (longitude < start_longitude) { start_longitude = longitude; }
      if (longitude > end_longitude) { end_longitude = longitude; }

      current_ptr = current_ptr->next_ptr;

   }

   get_int32_t (0) = 0;
   get_int32_t (4) = handle_ptr->n;
   get_int32_t (8) = 1;
   get_int32_t (12) = int32_t (start_longitude * 1e6);
   get_int32_t (16) = int32_t (end_longitude * 1e6);
   get_int32_t (20) = int32_t (start_latitude * 1e6);
   get_int32_t (24) = int32_t (end_latitude* 1e6);
   get_int32_t (28) = -1;
   get_int32_t (32) = -1;
   get_int16_t (36) = 0;
   get_int16_t (38) = 0;

}

const int32_t&
Gshhs::Header::get_int32_t (const Integer position) const
{
   return *((int32_t*)(buffer + position));
}

int32_t&
Gshhs::Header::get_int32_t (const Integer position)
{
   return *((int32_t*)(buffer + position));
}

const int16_t&
Gshhs::Header::get_int16_t (const Integer position) const
{
   return *((int16_t*)(buffer + position));
}

int16_t&
Gshhs::Header::get_int16_t (const Integer position)
{
   return *((int16_t*)(buffer + position));
}

void
Gshhs::Header::write (FILE* file,
                      const bool preserve_header)
{

   if (preserve_header)
   {
      Header temp_header (*this);
      temp_header.write (file, false);
   }
   else
   {
#ifndef WORDS_BIGENDIAN
      swap_endian ();
#endif
      fwrite (buffer, buffer_size, 1, file);
   }

}

int8_t
Gshhs::Header::get_level () const
{
   const bool old = (get_version () == 0);
   const int32_t& flag = get_flag ();
   return (old ? int8_t (flag) : (flag & 0xff));
}

bool
Gshhs::Header::is_cross_greenwich () const
{
   const bool old = (get_version () == 0);
   return (old ? (get_int16_t (36) >= 1) : ((get_flag () >> 16) & 1));
}

bool
Gshhs::Header::is_western () const
{
   const int32_t w = get_int32_t (12);
   const int32_t e = get_int32_t (16);
   return (w > 180000000 && e > 180000000);
}

Gshhs::Header*
Gshhs::get_header_ptr () const
{
   return new Header ();
}

void
Gshhs::add_polygon (FILE* file,
                    const Header& header,
                    const Real max_longitude)
{

   Point point;
   Point_2D p;
   Real last_y;

   const uint32_t n = header.get_int32_t (4);
   const bool antartica = (header.get_int32_t (20) * 1e-6 < -89);
   const bool cross_greenwich = header.is_cross_greenwich ();
   const bool western = header.is_western ();

   if (antartica)
   {
//      add (Point_2D (-68.9257, 0), true);
//      add (Point_2D (-89, 0), false);
//      add (Point_2D (-89, 360), false);
//      add (Point_2D (-68.9257, 360), false);
   }

   for (Integer i = 0; i < n; i++)
   {

      point.read (file);

      p.x = point.y * 1e-6;
      p.y = point.x * 1e-6;

      if (cross_greenwich && p.y > max_longitude) { p.y -= 360; }
      if ((antartica || western) && p.y > 180) { p.y -= 360; }

      if (antartica && (fabs (last_y - p.y) > 350))
      {
         add (Point_2D (-78.5, -180), true);
         add (Point_2D (-89.9, -180), false);
         add (Point_2D (-89.9, 180), false);
         add (Point_2D (-78.5, 180), false);
      }

      if (i < n - 1)
      {
         //const bool new_handle = (!antartica && (i == 0));
         const bool new_handle = (i == 0);
         add (p, new_handle);
      }

      last_y = p.y;

   }

}

void
Gshhs::write_simple_polygon (FILE* file,
                             const Polygon::Vertex* handle_ptr) const
{

   Point point;
   Integer n = handle_ptr->n;
   Polygon::Vertex* current_ptr = (Polygon::Vertex*)handle_ptr;

   for (Integer i = 0; i < n; i++)
   {

      const Point_2D& p = (const Point_2D&)(*current_ptr);

      point.x = int32_t (point.y * 1e6);
      point.y = int32_t (point.x * 1e6);

#ifndef WORDS_BIGENDIAN
      point.swap_endian ();
#endif

      fwrite (&point, sizeof (class Point), 1, file);

      current_ptr = current_ptr->next_ptr;

   }

}

Gshhs::Gshhs ()
{
}

Gshhs::Gshhs (const Dstring& file_path,
              const Integer level)
{

   Real max_longitude = 270;
   size_t point_size = sizeof (class Point);

   Header* header_ptr = get_header_ptr ();
   Header& header = *header_ptr;

   FILE* file = get_input_file (file_path);

   while (true)
   {

      try { header.read (file); } catch (const IO_Exception& ioe) { break; }

      const int32_t n = header.get_int32_t (4);
      const bool fail_level = (header.get_level () > level);
      if (fail_level) { fseek (file, point_size * n, SEEK_CUR); continue; }

      add_polygon (file, header, max_longitude);
      max_longitude = 180;

   }

   delete header_ptr;
   fclose (file);

}

Gshhs::~Gshhs ()
{
}

void
Gshhs::save (const Dstring& file_path) const
{

   Header* header_ptr = get_header_ptr ();
   Header& header = *header_ptr;

   FILE* file = get_output_file (file_path);

   const Polygon::Vertex* first_handle_ptr = this->get_first_handle_ptr ();
   Polygon::Vertex* current_handle_ptr = (Polygon::Vertex*)(first_handle_ptr);

   do
   {

      if (current_handle_ptr->n > 2)
      {
         header.reset (current_handle_ptr);
         header.write (file, false);
         write_simple_polygon (file, current_handle_ptr);
      }

      current_handle_ptr = current_handle_ptr->next_handle_ptr;

   }
   while (current_handle_ptr != first_handle_ptr);
   
   delete header_ptr;
   fclose (file);

}

Gshhs*
Gshhs::clip (const Polygon& clip_polygon)
{
   Gshhs* gshhs_ptr = new Gshhs ();
   Polygon::boolean_op (*gshhs_ptr, INTERSECTION, *this, clip_polygon);
   return gshhs_ptr;
}

Gshhs*
Gshhs::clip (const Domain_1D& domain_latitude,
             const Domain_1D& domain_longitude)
{

   Size_2D size_2d (0, 0);

   const Real start_latitude = domain_latitude.start;
   const Real end_latitude   = domain_latitude.end;
   const Real start_longitude = domain_longitude.start;
   const Real end_longitude   = domain_longitude.end;

   const Real latitude_span = end_latitude - start_latitude;
   const Real longitude_span = end_longitude - start_longitude;

   const Integer nominal_n = Integer (latitude_span) + 2;
   const Integer nominal_m = Integer (longitude_span) + 2;

   const Integer n = (size_2d.i >= 2 ? size_2d.i : nominal_n);
   const Integer m = (size_2d.j >= 2 ? size_2d.j : nominal_m);

   const Real delta_latitude = latitude_span / (n - 1);
   const Real delta_longitude = longitude_span / (m - 1);

   Polygon clip_polygon;

   for (Integer i = 0; i < n; i++)
   {
      Real latitude = start_latitude + i * delta_latitude;
      clip_polygon.add (Point_2D (latitude, start_longitude));
   }

   for (Integer i = 0; i < m; i++)
   {
      Real longitude = start_longitude + i * delta_longitude;
      clip_polygon.add (Point_2D (end_latitude, longitude));
   }

   for (Integer i = 0; i < n; i++)
   {
      Real latitude = end_latitude - i * delta_latitude;
      clip_polygon.add (Point_2D (latitude, end_longitude));
   }

   for (Integer i = 0; i < m; i++)
   {
      Real longitude = end_longitude - i * delta_longitude;
      clip_polygon.add (Point_2D (start_latitude, longitude));
   }

   return clip (clip_polygon);

}

void
Gshhs::cairo (const RefPtr<Context> cr,
              const Geodetic_Transform& transform) const
{
   transform.cairo (cr, *this);
}

void
Gshhs::cairo (const RefPtr<Context> cr,
              const Geodetic_Transform& transform,
              const Polygon& clip_polygon)
{
   Gshhs* clipped_gshhs_ptr = this->clip (clip_polygon);
   transform.cairo (cr, *clipped_gshhs_ptr);
   delete clipped_gshhs_ptr;
}

Dstring
Blue_Marble::get_tile_string (const Blue_Marble::Tile tile)
{

   Dstring tile_string;

   switch (tile)
   {
      case EAST: tile_string += "blue_marble_east"; break;
      case WEST: tile_string += "blue_marble_west"; break;
      default: throw Blue_Marble::Exception ("Invalid Tile");
   }

   return tile_string;

}

Color
Blue_Marble::get_color (const Lat_Long& lat_long,
                        const Real start_latitude,
                        const Real start_longitude,
                        FILE** blue_marble_files)
{

   uint8_t raw_datum[3];
   Integer pixel_size = sizeof (uint8_t) * 3;

   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;

   // 1 point = 1/120 degrees
   const Integer gi = Integer (rint ((longitude - start_longitude) * 120));
   const Integer gj = Integer (rint ((latitude - start_latitude) * 120));

   const Blue_Marble::Tile tile = (gi >= 21600 ? EAST : WEST);
   FILE* file = blue_marble_files[tile];

   const Integer tile_i = gi - (tile == EAST ? 21600 : 0);
   const Integer tile_j = 21600 - gj - 1;

   const long position = (tile_j * 21600 + tile_i) * pixel_size;
   fseek (file, position, SEEK_SET);
   fread (raw_datum, pixel_size, 1, file);

   return Color (Real (raw_datum[0]) / 255,
                 Real (raw_datum[1]) / 255,
                 Real (raw_datum[2]) / 255);

}

void
Blue_Marble::fill_raster (const Dstring& blue_marble_path,
                          const Transform_2D& transform_2d)

{

   const Dstring& east_file_path = blue_marble_path + "/blue_marble_east.bin";
   const Dstring& west_file_path = blue_marble_path + "/blue_marble_west.bin";

   FILE** blue_marble_files = new FILE*[2];
   blue_marble_files[EAST] = get_input_file (east_file_path);
   blue_marble_files[WEST] = get_input_file (west_file_path);

   const Real delta_2 = Real (1) / Real (240);
   const Real start_latitude = -90 + delta_2;
   const Real start_longitude = -180 + delta_2;

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real x = Real (i);

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real y = Real (j);
         const Lat_Long lat_long = transform_2d.reverse (x, y);

         const Color& color = get_color (lat_long, start_latitude,
            start_longitude, blue_marble_files);

         set_pixel (i, j, color);

      }

   }

   fclose (blue_marble_files[EAST]);
   fclose (blue_marble_files[WEST]);

}

Blue_Marble::Blue_Marble (const Dstring& blue_marble_path,
                          const Transform_2D& transform_2d,
                          const Size_2D& size_2d)
   : Raster (size_2d)
{
   fill_raster (blue_marble_path, transform_2d);
}

Blue_Marble::~Blue_Marble ()
{
}

Blue_Marble::Exception::Exception (const Dstring& description)
   : denise::Exception ("Blue_Marble::Exception", description)
{
}

Gtopo30::Tile
Gtopo30::get_tile (const Integer gtopo30_i,
                   const Integer gtopo30_j)
{

   Integer i = gtopo30_i;
   Integer j = denise::imodulo (gtopo30_j, 43200);

   Integer tile_index;

   if (j >= 3600)
   {
      tile_index = ((j - 3600) / 6000) * 9 + (i / 4800) + 6;
   }
   else
   {
      tile_index = i / 7200;
   }

   return Gtopo30::Tile (tile_index);

}

Dstring
Gtopo30::get_tile_string (const Gtopo30::Tile tile)
{

   Dstring tile_string;

   switch (tile)
   {
      case W180S60: tile_string += "W180S60"; break;
      case W120S60: tile_string += "W120S60"; break;
      case W060S60: tile_string += "W060S60"; break;
      case W000S60: tile_string += "W000S60"; break;
      case E060S60: tile_string += "E060S60"; break;
      case E120S60: tile_string += "E120S60"; break;
      case W180S10: tile_string += "W180S10"; break;
      case W140S10: tile_string += "W140S10"; break;
      case W100S10: tile_string += "W100S10"; break;
      case W060S10: tile_string += "W060S10"; break;
      case W020S10: tile_string += "W020S10"; break;
      case E020S10: tile_string += "E020S10"; break;
      case E060S10: tile_string += "E060S10"; break;
      case E100S10: tile_string += "E100S10"; break;
      case E140S10: tile_string += "E140S10"; break;
      case W180N40: tile_string += "W180N40"; break;
      case W140N40: tile_string += "W140N40"; break;
      case W100N40: tile_string += "W100N40"; break;
      case W060N40: tile_string += "W060N40"; break;
      case W020N40: tile_string += "W020N40"; break;
      case E020N40: tile_string += "E020N40"; break;
      case E060N40: tile_string += "E060N40"; break;
      case E100N40: tile_string += "E100N40"; break;
      case E140N40: tile_string += "E140N40"; break;
      case W180N90: tile_string += "W180N90"; break;
      case W140N90: tile_string += "W140N90"; break;
      case W100N90: tile_string += "W100N90"; break;
      case W060N90: tile_string += "W060N90"; break;
      case W020N90: tile_string += "W020N90"; break;
      case E020N90: tile_string += "E020N90"; break;
      case E060N90: tile_string += "E060N90"; break;
      case E100N90: tile_string += "E100N90"; break;
      case E140N90: tile_string += "E140N90"; break;
      default: throw Gtopo30::Exception ("Invalid Tile");
   }

   return tile_string;

}

Size_2D
Gtopo30::get_tile_size (const Gtopo30::Tile tile)
{

   Size_2D size_2d;
   Integer tile_index = Integer (tile);

   if (tile_index >= 6 && tile_index <= 32)
   {
      size_2d.i = 4800;
      size_2d.j = 6000;
   }
   else
   if (tile_index >= 0 && tile_index <= 5)
   {
      size_2d.i = 7200;
      size_2d.j = 3600;
   }
   else
   {
      throw Gtopo30::Exception ("Invalid Tile");
   }

   return size_2d;

}

Index_2D
Gtopo30::get_tile_index (const Gtopo30::Tile tile)
{

   Index_2D index_2d;
   Integer tile_index = Integer (tile);

   if (tile_index >= 6 && tile_index <= 32)
   {
      index_2d.i = ((tile_index - 6) % 9) * 4800;
      index_2d.j = ((tile_index - 6) / 9) * 6000 + 3600;
   }
   else
   if (tile_index >= 0 && tile_index <= 5)
   {
      index_2d.i = tile_index * 7200;
      index_2d.j = 0;
   }
   else
   {
      throw Gtopo30::Exception ("Invalid Tile");
   }

   return index_2d;

}

Domain_2D
Gtopo30::get_tile_domain (const Gtopo30::Tile tile)
{

   Domain_1D domain_latitude = get_tile_domain_latitude (tile);
   Domain_1D domain_longitude = get_tile_domain_longitude (tile);

   return Domain_2D (domain_latitude, domain_longitude);

}

Domain_1D
Gtopo30::get_tile_domain_latitude (const Gtopo30::Tile tile)
{

   Domain_1D domain_1d;
   Integer tile_index = Integer (tile);

   const Real delta = Real (1) / Real (120);

   if (tile_index >= 0 && tile_index <= 32)
   {

      Integer j = tile_index / 9;

      if (j == 3)
      {
         domain_1d.start = -90 + delta;
         domain_1d.end = -60 - delta;
      }
      else
      {
         domain_1d.end = 90 - j * 50 - delta;
         domain_1d.start = 90 - (j+1) * 50 + delta;
      }

   }
   else
   {
      throw Gtopo30::Exception ("Invalid Tile");
   }

   return domain_1d;

}

Domain_1D
Gtopo30::get_tile_domain_longitude (const Gtopo30::Tile tile)
{

   Domain_1D domain_1d;
   Integer tile_index = Integer (tile);

   const Real delta = Real (1) / Real (120);

   if (tile_index >= 0 && tile_index <= 26)
   {
      Integer j = tile_index % 9;
      domain_1d.start = -180 + j * 40 + delta;
      domain_1d.end = -180 + (j + 1) * 40 - delta;
   }
   else
   if (tile_index >= 27 && tile_index <= 32)
   {
      Integer j = tile_index - 27;
      domain_1d.start = -180 + j * 60 + delta;
      domain_1d.end = -180 + (j + 1) * 60 - delta;
   }
   else
   {
      throw Gtopo30::Exception ("Invalid Tile");
   }

   return domain_1d;

}

Real*
Gtopo30::get_data (const Transform_2D& transform)
{

   file_array = new FILE*[33];
   tile_size_array = new Size_2D[33];
   tile_index_array = new Index_2D[33];

   for (Integer tile = 0; tile < 33; tile++)
   {  
      file_array[tile] = NULL;
   }   

   Lat_Long lat_long;
   Real& latitude = lat_long.latitude;
   Real& longitude = lat_long.longitude;

   Real* data = new Real[size_2d.i * size_2d.j];

   for (Integer i = 0; i < size_2d.i; i++)
   {

      const Real x = Real (i);

      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real y = Real (j);

         transform.reverse (latitude, longitude, x, y);
         const Real& datum = get_datum (lat_long);
         data[i*size_2d.j + j] = datum;

      }

   }

   return data;

}

void
Gtopo30::fill_raster (const Real* data,
                      const Real max_gradient)
{

   //const Color& land = Color::hsb (0.5, 0.27, 0.3);
   //const Color& sea = Color::hsb (0.67, 0.67, 0.3);
   const Color& land = Color::land ();
   const Color& sea = Color::sea ();
   const Real mg = max_gradient;

   for (Integer i = 0; i < size_2d.i; i++)
   {

      for (Integer j = 0; j < size_2d.j; j++)
      {

         const Real datum = data[i * size_2d.j + j];

         const Real h = std::min (std::max (datum / 1200.0, 0.0), 1.0);
         const Real hue = 0.45 - h * 0.4;
         //const Real brightness = h * 0.7 + 0.28;
         const Real brightness = h * 0.3 + 0.2;
         Color color = (datum > 0 ? Color::hsb (hue, 0.54, brightness) : sea);
         //Color color = (datum > 0 ? land : sea);
         //Color color = color_chooser.get_color (datum);

         Real gradient = 0;

         if (i != 0 && i != size_2d.i - 1 && j != 0 && i != size_2d.j - 1)
         {

            const Real datum_e = data[(i+1) * size_2d.j + j];
            const Real datum_w = data[(i-1) * size_2d.j + j];
            const Real datum_n = data[i * size_2d.j + (j+1)];
            const Real datum_s = data[i * size_2d.j + (j-1)];

            gradient = datum_e - datum_n - datum_w + datum_s;
            if (!gsl_finite (gradient)) { gradient = 0; }

            Real scale = 1 + gradient * 0.0025;
            scale = std::max (0.7, std::min (1.3, scale));
            color.scale_brightness (scale);

         }

         set_pixel (i, j, color);

      }
   }

}

Gtopo30::Gtopo30 (const Dstring& gtopo30_path,
                  const Transform_2D& transform_2d,
                  const Size_2D& size_2d,
                  const Real max_gradient)
   : Raster (size_2d),
     gtopo30_path (gtopo30_path),
     delta (Real (1) / Real (120)),
     start_latitude (-90 + delta/2),
     start_longitude (-180 + delta/2)
{
   Real* data = get_data (transform_2d);
   fill_raster (data, max_gradient);
   delete[] data;
}

Gtopo30::~Gtopo30 ()
{

   for (Integer tile = 0; tile < 33; tile++)
   {
      FILE*& file = file_array[tile];
      if (file != NULL) { fclose (file); }
   }

   delete[] file_array;
   delete[] tile_size_array;
   delete[] tile_index_array;

}

Real
Gtopo30::get_datum (const Lat_Long& lat_long) const
{

   int16_t datum;

   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;

   Integer gi = Integer (rint ((longitude - start_longitude) / delta));
   Integer gj = Integer (rint ((latitude - start_latitude) / delta));

   Gtopo30::Tile tile = get_tile (gi, gj);
   FILE*& file = file_array[tile];
   Size_2D& tile_size = tile_size_array[tile];
   Index_2D& tile_index = tile_index_array[tile];

   if (file == NULL)
   { 

      Dstring tile_string = get_tile_string (tile);
      Dstring file_path = gtopo30_path + "/" + tile_string + ".DEM";
      file = get_input_file (file_path);

      tile_size = get_tile_size (tile);
      tile_index = get_tile_index (tile);

   }

   const Integer tile_i = gi - tile_index.i;
   const Integer tile_j = tile_size.j - (gj - tile_index.j) - 1;

   long position = (tile_j * tile_size.i + tile_i) * sizeof (int16_t);
   fseek (file, position, SEEK_SET);
   fread (&datum, sizeof (int16_t), 1, file);
#ifndef WORDS_BIGENDIAN
   swap_endian (&datum, sizeof (int16_t));
#endif

   return (datum < -1000 ? GSL_NAN : Real (datum));

}

Gtopo30::Exception::Exception (const Dstring& description)
   : denise::Exception ("Gtopo30::Exception", description)
{
}

uint16_t
Land_Mask::get_degree_box_size () const
{
   return (n * n) / 8;
}

uint8_t
Land_Mask::get_byte_mask (const uint8_t bj) const
{
//   uint8_t byte_mask = 0x80;
//   return (byte_mask >> bj);
   switch (bj)
   {
      default: return 0;
      case 0: return 1;
      case 1: return 2;
      case 2: return 4;
      case 3: return 8;
      case 4: return 16;
      case 5: return 32;
      case 6: return 64;
      case 7: return 128;
   }
}

void
Land_Mask::acquire_lat_long (Real& latitude,
                             Real& longitude,
                             const int16_t i,
                             const int16_t j,
                             const uint8_t ii,
                             const uint8_t jj) const
{
   const Real h = Real (1) / Real (n);
   const Real h_2 = h / 2;
   latitude = Real (i) + Real (ii) * h + h_2 + Real (start_i);
   longitude = Real (j) + Real (jj) * h + h_2 + Real (start_j);
}

void
Land_Mask::acquire_indices (int16_t& i,
                            int16_t& j,
                            uint8_t& ii,
                            uint8_t& jj,
                            const Lat_Long& lat_long) const
{

   const Real h = Real (1) / Real (n);
   const Real epsilon = h / 5000;
   const Real h_2 = h / 2 - epsilon;

   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;

   const Real floor_latitude = floor (latitude + epsilon);
   const Real floor_longitude = floor (longitude + epsilon);

   i = int16_t (Integer (floor_latitude) - start_i);
   j = int16_t (Integer (floor_longitude) - start_j);

   ii = Integer (round ((latitude - floor_latitude - h_2) / h));
   jj = Integer (round ((longitude - floor_longitude - h_2) / h));

}

Land_Mask::Land_Mask (const Polygon& polygon,
                      const uint16_t n,
                      const int16_t start_i,
                      const int16_t start_j,
                      const uint16_t ni,
                      const uint16_t nj)
   : n (n),
     start_i (start_i),
     start_j (start_j),
     ni (ni),
     nj (nj)
{

   const uint16_t degree_box_size = get_degree_box_size ();
   const uint32_t buffer_size = ni * nj * degree_box_size;
   buffer = new uint8_t[buffer_size];

   populate (polygon, start_i, start_j, ni, nj);

}
 
Land_Mask::Land_Mask (const Gtopo30& gtopo30,
                      const uint16_t n,
                      const int16_t start_i,
                      const int16_t start_j,
                      const uint16_t ni,
                      const uint16_t nj)
   : n (n),
     start_i (start_i),
     start_j (start_j),
     ni (ni),
     nj (nj)
{

   const uint16_t degree_box_size = get_degree_box_size ();
   const uint32_t buffer_size = ni * nj * degree_box_size;
   buffer = new uint8_t[buffer_size];

   populate (gtopo30, start_i, start_j, ni, nj);

}
 
Land_Mask::Land_Mask (const Dstring& land_mask_file_path)
{


   FILE* file = get_input_file (land_mask_file_path);

   fread (&n, 1, sizeof (uint16_t), file);
   fread (&start_i, 1, sizeof (int16_t), file);
   fread (&start_j, 1, sizeof (int16_t), file);
   fread (&ni, 1, sizeof (uint16_t), file);
   fread (&nj, 1, sizeof (uint16_t), file);

#ifndef WORDS_BIGENDIAN
   denise::swap_endian (&n, sizeof (uint16_t));
   denise::swap_endian (&start_i, sizeof (int16_t));
   denise::swap_endian (&start_j, sizeof (int16_t));
   denise::swap_endian (&ni, sizeof (uint16_t));
   denise::swap_endian (&nj, sizeof (uint16_t));
#endif

   const uint16_t degree_box_size = get_degree_box_size ();
   const uint32_t buffer_size = ni * nj * degree_box_size;
   buffer = new uint8_t[buffer_size];

   fread (buffer, buffer_size, sizeof (uint8_t), file);

   fclose (file);

}

Land_Mask::~Land_Mask ()
{
   delete[] buffer;
}

void
Land_Mask::populate (const Polygon& polygon,
                     const int16_t start_i,
                     const int16_t start_j,
                     const uint16_t ni,
                     const uint16_t nj)
{

   const uint16_t n_tilde = n / 4;
   const uint16_t n_tilde_2 = n_tilde / 2;
   const bool n_tilde_odd = (n_tilde % 2);
   const uint32_t box_size = get_degree_box_size ();

   uint16_t bi;
   uint8_t bj;
   Point_2D point;
   Real& latitude = point.x;
   Real& longitude = point.y;

   for (int16_t i = 0; i < ni; i++)
   {
      for (int16_t j = 0; j < nj; j++)
      {

         for (uint8_t ii = 0; ii < n; ii++)
         {

            bool ii_odd = (ii % 2);
            uint16_t iii = (n_tilde_odd ? (ii / 2) * n_tilde : ii * n_tilde_2);

            for (uint8_t jj = 0; jj < n; jj++)
            {

               acquire_lat_long (latitude, longitude, i, j, ii, jj);
               const bool is_land = polygon.contains (point);

               if (n_tilde_odd)
               {
                  bi = (ii_odd ? iii + (jj+4)/8 + n_tilde_2 : iii + jj/8);
                  bj = (ii_odd ? (jj+4)%8 : jj%8);
               }
               else
               {
                  bi = iii + jj/8;
                  bj = jj%8;
               }

               uint8_t& byte = buffer[(i * nj + j) * box_size + bi];
               const uint8_t byte_mask = get_byte_mask (bj);

               byte = (is_land ? (byte | byte_mask) : (byte & !byte_mask));

            }
         }
      }
   }
 
}

void
Land_Mask::populate (const Gtopo30& gtopo30,
                     const int16_t start_i,
                     const int16_t start_j,
                     const uint16_t ni,
                     const uint16_t nj)
{

   const uint16_t n_tilde = n / 4;
   const uint16_t n_tilde_2 = n_tilde / 2;
   const bool n_tilde_odd = (n_tilde % 2);
   const uint32_t box_size = get_degree_box_size ();

   uint16_t bi;
   uint8_t bj;
   Lat_Long lat_long;
   Real& latitude = lat_long.latitude;
   Real& longitude = lat_long.longitude;

   for (int16_t i = 0; i < ni; i++)
   {
      for (int16_t j = 0; j < nj; j++)
      {

         for (uint8_t ii = 0; ii < n; ii++)
         {

            bool ii_odd = (ii % 2);
            uint16_t iii = (n_tilde_odd ? (ii / 2) * n_tilde : ii * n_tilde_2);

            for (uint8_t jj = 0; jj < n; jj++)
            {

               acquire_lat_long (latitude, longitude, i, j, ii, jj);
               bool is_land = !gsl_isnan (gtopo30.get_datum (lat_long));

               if (n_tilde_odd)
               {
                  bi = (ii_odd ? iii + (jj+4)/8 + n_tilde_2 : iii + jj/8);
                  bj = (ii_odd ? (jj+4)%8 : jj%8);
               }
               else
               {
                  bi = iii + jj/8;
                  bj = jj%8;
               }

               uint8_t& byte = buffer[(i * nj + j) * box_size + bi];
               uint8_t byte_mask = get_byte_mask (bj);

               byte = (is_land ? (byte | byte_mask) : (byte & ~byte_mask));

            }
         }
      }
   }
 
}

void
Land_Mask::save (const Dstring& land_mask_file_path) const
{

   FILE* file = get_output_file (land_mask_file_path);

   const uint32_t box_size = get_degree_box_size ();
   const uint32_t buffer_size = ni * nj * box_size;

   uint16_t n_copy = n;
   int16_t start_i_copy = start_i;
   int16_t start_j_copy = start_j;
   uint16_t ni_copy = ni;
   uint16_t nj_copy = nj;

#ifndef WORDS_BIGENDIAN
   denise::swap_endian (&n_copy, sizeof (uint16_t));
   denise::swap_endian (&start_i_copy, sizeof (int16_t));
   denise::swap_endian (&start_j_copy, sizeof (int16_t));
   denise::swap_endian (&ni_copy, sizeof (uint16_t));
   denise::swap_endian (&nj_copy, sizeof (uint16_t));
#endif

   fwrite (&n_copy, 1, sizeof (uint16_t), file);
   fwrite (&start_i_copy, 1, sizeof (uint16_t), file);
   fwrite (&start_j_copy, 1, sizeof (uint16_t), file);
   fwrite (&ni_copy, 1, sizeof (uint16_t), file);
   fwrite (&nj_copy, 1, sizeof (uint16_t), file);
   fwrite (buffer, buffer_size, sizeof (uint8_t), file);

   fclose (file);

}

bool
Land_Mask::is_land (const Lat_Long& lat_long) const
{

   const uint16_t n_tilde = n / 4;
   const uint16_t n_tilde_2 = n_tilde / 2;
   const bool n_tilde_odd = (n_tilde % 2);
   const uint32_t box_size = get_degree_box_size ();

   int16_t i, j;
   uint8_t ii, jj;
   acquire_indices (i, j, ii, jj, lat_long);
   if (i < 0 || i >= ni || j < 0 || j >= nj) { return false; }

   uint16_t bi;
   uint8_t bj;

   if (n_tilde_odd)
   {
      const bool ii_odd = (ii % 2);
      const uint8_t iii = (ii / 2) * n_tilde;
      bi = (ii_odd ? iii + (jj+4)/8 + n_tilde_2 : iii + jj/8);
      bj = (ii_odd ? (jj+4)%8 : jj%8);
   }
   else
   {
      const uint16_t iii = (ii * n_tilde_2);
      bi = iii + jj/8;
      bj = jj%8;
   }

   const uint8_t& byte = buffer[(i * nj + j) * box_size + bi];
   const uint8_t byte_mask = get_byte_mask (bj);

   return (byte & byte_mask);

}

Real
Land_Mask::get_fetch (const Lat_Long& lat_long,
                      const Real bearing) const
{

   const Geodesy geodesy;
   const Real max_d = 1000e3;
   const Real dd = LATITUDE_LENGTH / n;

   if (is_land (lat_long)) { return 0; }

   for (Real d = dd; d < max_d; d += dd)
   {
      const Journey::Simple simple_journey (lat_long, d, bearing, geodesy);
      const Lat_Long& ll = simple_journey.get_destination ();
      if (is_land (ll)) { return d; }
   }

   return max_d;

}

Real
Land_Mask::get_effective_fetch (const Lat_Long& lat_long,
                                const Real bearing) const
{

   const Integer n = 17;
   Tuple deviation_bearing_tuple (n, 5, 85);

   Real fetch = get_fetch (lat_long, bearing);
   Real total_weight = 1;

   for (Tuple::const_iterator iterator = deviation_bearing_tuple.begin ();
        iterator != deviation_bearing_tuple.end (); iterator++)
   {

      const Real deviation_bearing = *(iterator);
      const Real weight = cos (deviation_bearing * DEGREE_TO_RADIAN);

      const Real b_plus = bearing + deviation_bearing;
      const Real b_minus = bearing - deviation_bearing;
      const Real raw_fetch_plus = get_fetch (lat_long, b_plus);
      const Real raw_fetch_minus = get_fetch (lat_long, b_minus);

      fetch += weight * (raw_fetch_plus + raw_fetch_minus);
      total_weight += 2 * weight;

   }

   return (fetch / total_weight);

}

void
Fetch::populate (const Land_Mask& land_mask)
{

   Lat_Long lat_long;
   Real& latitude = lat_long.latitude;
   Real& longitude = lat_long.longitude;

   const Real h = Real (1) / Real (n);
   const Real h_2 = h / 2;

   const Real d_bearing = 360.0 / nb;

   const uint16_t nii = ni * n;
   const uint16_t njj = nj * n;

   const uint32_t block_size = nii * njj;
   buffer = new uint8_t[nb * block_size];

   for (int16_t ii = 0; ii < nii; ii++)
   {

      latitude = Real (start_i) + ii * h + h_2;

      for (int16_t jj = 0; jj < njj; jj++)
      {

         longitude = Real (start_j) + jj * h + h_2;

         for (Integer b = 0; b < nb; b++)
         {

            const Real bearing = b * d_bearing;
            const Real f = land_mask.get_effective_fetch (lat_long, bearing);
//            const Real f = land_mask.get_fetch (lat_long, bearing);
            const Real fetch = std::min (f, Real (510e3));

            uint8_t& byte = buffer[b * block_size + ii * njj + jj];
            byte = uint8_t (round (fetch / 2e3));

         }

      }
   }

}

Fetch::Fetch (const Gshhs& gshhs,
              const uint16_t nb,
              const uint16_t n,
              const int16_t start_i,
              const int16_t start_j,
              const uint16_t ni,
              const uint16_t nj)
   : nb (nb),
     n (n),
     start_i (start_i),
     start_j (start_j),
     ni (ni),
     nj (nj)
{
   const Land_Mask land_mask (gshhs, n, start_i, start_j, ni, nj);
   populate (land_mask);
}

Fetch::Fetch (const Gtopo30& gtopo30,
              const uint16_t nb,
              const uint16_t n,
              const int16_t start_i,
              const int16_t start_j,
              const uint16_t ni,
              const uint16_t nj)
   : nb (nb),
     n (n),
     start_i (start_i),
     start_j (start_j),
     ni (ni),
     nj (nj)
{
   const Land_Mask land_mask (gtopo30, n, start_i, start_j, ni, nj);
   populate (land_mask);
}

Fetch::Fetch (const Land_Mask& land_mask,
              const uint16_t nb,
              const uint16_t n,
              const int16_t start_i,
              const int16_t start_j,
              const uint16_t ni,
              const uint16_t nj)
   : nb (nb),
     n (n),
     start_i (start_i),
     start_j (start_j),
     ni (ni),
     nj (nj)
{
   populate (land_mask);
}

Fetch::Fetch (const Dstring& fetch_file_path)
{

   int i = 0;
   //FILE* file = get_input_file (fetch_file_path);
   gzFile file = get_gzfile (fetch_file_path);

   //fread (&nb, 1, sizeof (uint16_t), file);
   //fread (&n, 1, sizeof (uint16_t), file);
   //fread (&start_i, 1, sizeof (int16_t), file);
   //fread (&start_j, 1, sizeof (int16_t), file);
   //fread (&ni, 1, sizeof (uint16_t), file);
   //fread (&nj, 1, sizeof (uint16_t), file);
   gzread (file, &nb, sizeof (uint16_t));
   gzread (file, &n, sizeof (uint16_t));
   gzread (file, &start_i, sizeof (int16_t));
   gzread (file, &start_j, sizeof (int16_t));
   gzread (file, &ni, sizeof (uint16_t));
   gzread (file, &nj, sizeof (uint16_t));

#ifndef WORDS_BIGENDIAN
   denise::swap_endian (&nb, sizeof (uint16_t));
   denise::swap_endian (&n, sizeof (uint16_t));
   denise::swap_endian (&start_i, sizeof (int16_t));
   denise::swap_endian (&start_j, sizeof (int16_t));
   denise::swap_endian (&ni, sizeof (uint16_t));
   denise::swap_endian (&nj, sizeof (uint16_t));
#endif

   const uint16_t nii = ni * n;
   const uint16_t njj = nj * n;
   const uint32_t buffer_size = nb * nii * njj;

   int r;

   buffer = new uint8_t[buffer_size];
   //fread (buffer, buffer_size, sizeof (uint8_t), file);
   r = gzread (file, buffer, buffer_size);


   //fclose (file);
   gzclose (file);

}

Fetch::~Fetch ()
{
   delete[] buffer;
}

void
Fetch::save (const Dstring& fetch_file_path) const
{

   const uint16_t nii = ni * n;
   const uint16_t njj = nj * n;
   const uint32_t buffer_size = nb * nii * njj;

   FILE* file = get_output_file (fetch_file_path);

   uint16_t nb_copy = nb;
   uint16_t n_copy = n;
   int16_t start_i_copy = start_i;
   int16_t start_j_copy = start_j;
   uint16_t ni_copy = ni;
   uint16_t nj_copy = nj;

#ifndef WORDS_BIGENDIAN
   denise::swap_endian (&nb_copy, sizeof (uint16_t));
   denise::swap_endian (&n_copy, sizeof (uint16_t));
   denise::swap_endian (&start_i_copy, sizeof (int16_t));
   denise::swap_endian (&start_j_copy, sizeof (int16_t));
   denise::swap_endian (&ni_copy, sizeof (uint16_t));
   denise::swap_endian (&nj_copy, sizeof (uint16_t));
#endif

   fwrite (&nb_copy, 1, sizeof (uint16_t), file);
   fwrite (&n_copy, 1, sizeof (uint16_t), file);
   fwrite (&start_i_copy, 1, sizeof (int16_t), file);
   fwrite (&start_j_copy, 1, sizeof (int16_t), file);
   fwrite (&ni_copy, 1, sizeof (uint16_t), file);
   fwrite (&nj_copy, 1, sizeof (uint16_t), file);
   fwrite (buffer, buffer_size, sizeof (uint8_t), file);

   fclose (file);

}

Real
Fetch::get_fetch (const Real& bearing,
                  const Lat_Long& lat_long) const
{

   const Real d_bearing = 360.0 / nb;
   const uint16_t b = uint16_t (round (bearing / d_bearing)) % nb;

   const Real h = Real (1) / Real (n);
   const Real epsilon = h / 5000;
   const Real h_2 = h / 2 - epsilon;

   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;

   const uint16_t nii = ni * n;
   const uint16_t njj = nj * n;

   const int16_t ii = int16_t (round ((latitude - Real (start_i) - h_2) / h));
   if (ii < 0 || ii >= nii) { return 510e3; }

   const int16_t jj = int16_t (round ((longitude - Real (start_j) - h_2) / h));
   if (jj < 0 || jj >= njj) { return 510e3; }

   const uint32_t block_size = nii * njj;
   const uint8_t& byte = buffer[b * block_size + ii * njj + jj];

   const Real fetch = Real (byte) * 2e3;
   return fetch;

}

Real
Nsd::get_lat_long_number (const Dstring& lat_long_string)
{

   Real x = 0;

   Dstring str (lat_long_string);
   str.chop ();

   const Tokens tokens (str, "-");
   Integer n = tokens.size ();

   x += stof (tokens[0]);
   if (n > 1) { x += stof (tokens[1]) / 60; }
   if (n > 2) { x += stof (tokens[2]) / 3600; }

   char last_char = *(lat_long_string.rbegin ());
   if (last_char == 'S' || last_char == 'W') { x *= -1; }

   return x;

}

void
Nsd::wmo_init (const Dstring& str)
{

   const Tokens tokens (str, ";");

   wmo_id = stoi (tokens[0]) * 1000;
   wmo_id += stoi (tokens[1]);
   icao_id = tokens[2];

   name = tokens[3];
   country = tokens[5];

   lat_long.latitude = get_lat_long_number (tokens[7]);
   lat_long.longitude = get_lat_long_number (tokens[8]);

   elevation = stof (tokens[11]);

}

void
Nsd::icao_init (const Dstring& str)
{

   const Tokens tokens (str, ";");

   icao_id = tokens[0];
   wmo_id = stoi (tokens[1]) * 1000;
   wmo_id += stoi (tokens[2]);

   name = tokens[3];
   country = tokens[5];

   lat_long.latitude = get_lat_long_number (tokens[7]);
   lat_long.longitude = get_lat_long_number (tokens[8]);

   elevation = stof (tokens[11]);

}

Nsd::Nsd (const Dstring& nsd_string,
          const bool icao)
{
   if (icao) { icao_init (nsd_string); }
   else { wmo_init (nsd_string); }
}

Nsd_Wmo::Nsd_Wmo (const Dstring& file_path)
{

   string il;
   Dstring input_line;
   
   ifstream file (file_path.get_string ());

   while (file)
   {
      std::getline (file, il);
      if (il == "") { continue; }
      Dstring input_line (il);
      Nsd nsd (input_line, false);
      insert (make_pair (nsd.wmo_id, nsd));
   }

}

const Nsd&
Nsd_Wmo::get_nsd (const Integer wmo_id) const
{

   Nsd_Wmo::const_iterator iterator = find (wmo_id);

   if (iterator != end ())
   {
      return (iterator->second);
   }
   else
   {
      Dstring str = Dstring::render ("Nsd_Wmo cannot find: %05d", wmo_id);
      throw Nsd_Exception (str);
   }

}

Nsd_Icao::Nsd_Icao (const Dstring& file_path)
{

   ifstream file (file_path.get_string ().c_str ());

   string il;

   while (file)
   {
      std::getline (file, il);
      if (il == "") { continue; }
      const Dstring input_line (il);
      Nsd nsd (input_line, true);
      insert (make_pair (nsd.icao_id, nsd));
   }

}

const Nsd&
Nsd_Icao::get_nsd (const Dstring& icao_id) const
{

   Nsd_Icao::const_iterator iterator = find (icao_id);

   if (iterator != end ())
   {
      return (iterator->second);
   }
   else
   {
      Dstring str = "Nsd_Icao cannot find: " + icao_id;
      throw Nsd_Exception (str);
   }

}

Nsd_Exception::Nsd_Exception (const Dstring& description)
   : Exception ("Nsd_Exception", description)
{
}

void
Outl::read (const Dstring& file_path)
{

   long header_pos;
   Lat_Long lat_long;
   int32_t raw_latitude, raw_longitude;
   int32_t pos, number_of_lines, number_of_points;

   FILE* file = get_input_file (file_path);
   fread (&number_of_lines, sizeof (int32_t), 1, file);

#ifndef WORDS_BIGENDIAN
   denise::swap_endian (&number_of_lines, sizeof (int32_t));
#endif

   for (Integer i = 0; i < number_of_lines; i++)
   {

      header_pos = (i * 6 + 5) * sizeof (int32_t);
      fseek (file, header_pos, SEEK_SET);

      fread (&pos, sizeof (int32_t), 1, file);
      fread (&number_of_points, sizeof (int32_t), 1, file);

#ifndef WORDS_BIGENDIAN
      denise::swap_endian (&pos, sizeof (int32_t));
      denise::swap_endian (&number_of_points, sizeof (int32_t));
#endif

      pos *= sizeof (int32_t);
      number_of_points /= 2;

      fseek (file, pos, SEEK_SET);

      for (Integer j = 0; j < number_of_points; j++)
      {

         fread (&raw_latitude, sizeof (int32_t), 1, file);
         fread (&raw_longitude, sizeof (int32_t), 1, file);

#ifndef WORDS_BIGENDIAN
         denise::swap_endian (&raw_latitude, sizeof (int32_t));
         denise::swap_endian (&raw_longitude, sizeof (int32_t));
#endif

         lat_long.latitude = Real (raw_latitude) * 1e-4;
         lat_long.longitude = -Real (raw_longitude) * 1e-4;

         add (lat_long, (j == 0));

      }

   }

   fclose (file);

}

Outl::Outl (const Dstring& file_path)
{
   read (file_path);
}

void
Outl::cairo (const RefPtr<Context> cr,
             const Geodetic_Transform& transform) const
{
   transform.cairo (cr, *this);
}

void
Kidney::read (const Dstring& file_path)
{

   ifstream file (file_path.get_string ());
   string is;

   for (string is; std::getline (file, is); )
   {

      const Dstring input_string (is);
      const Tokens tokens (input_string, ":");
      const Real level = tokens.real (0);

      Polygon* polygon_ptr = new Polygon ();
      insert (make_pair (level, polygon_ptr));

      for (Integer i = 1; i < tokens.size (); i++)
      {
         const Tokens lat_long_tokens (tokens[i]);
         const Real latitude = lat_long_tokens.real (0);
         const Real longitude = lat_long_tokens.real (1);
         const Point_2D point (latitude, longitude);
         polygon_ptr->add (Point_2D (latitude, longitude));
      }

   }

   file.close ();

}

Kidney::Kidney (const Dstring& file_path)
{
   read (file_path);
}

Kidney::~Kidney ()
{
   for (Kidney::iterator iterator = begin (); iterator != end (); iterator++)
   {
      Polygon* polygon_ptr = iterator->second;
      delete polygon_ptr;
   }
}

void
Kidney::cairo (const RefPtr<Context> cr,
               const Geodetic_Transform& transform) const
{

   RefPtr<Pattern> pattern = cr->get_source ();
   const Color color (cr);
   const Real line_width = cr->get_line_width ();

   for (Kidney::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Real threshold = iterator->first;
      const Polygon* polygon_ptr = iterator->second;

      const bool highlight = fabs (threshold - 50) < 0.5;

      if (highlight)
      {
         const Real lw = line_width * 2;
         const Color c (color.r, color.g, color.b, color.a * 2);
         cr->save ();
         cr->set_line_width (line_width * 2);
         Color (color.r, color.g, color.b, color.a * 2).cairo (cr);
      }

      polygon_ptr->cairo (cr, transform);
      cr->stroke ();

      if (highlight)
      {
         cr->restore ();
      }

   }

}

