//
// grib.cc
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
// but WITHOUT ANY WARRANTY), without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with libdenise.  If not, see <http://www.gnu.org/licenses/>.

#include "grib.h"
#include <denise/thermo.h>
#include <jasper/jasper.h>

using namespace std;
using namespace denise;

Grib::Block::Block ()
{
}

Grib::Block::Block (const uint32_t n,
                    const bool set_zero)
   : n (n)
{
   buffer = new uint8_t[n];
   if (set_zero) { for (uint16_t i = 0; i < n; i++) { buffer[i] = 0; } }
}

Grib::Block::Block (FILE* file,
                    const uint32_t address)
{

   n = get_size (file, address);

   buffer = new uint8_t[n];
   fseek (file, address, SEEK_SET);
   fread (buffer, n, 1, file);

}

Grib::Block::Block (const Block& block)
   : n (block.n)
{
   buffer = new uint8_t[n];
   memcpy (buffer, block.buffer, n);
}

Grib::Block::~Block ()
{
   delete[] buffer;
}

uint32_t
Grib::Block::get_size (FILE* file,
                       const uint32_t address)
{
   return get_uint (file, address, 3);
}

bool
Grib::Block::get_bit (const uint32_t position,
                      const uint8_t bit_number) const
{

   bool boolean;
   const uint8_t c = buffer[position];

   switch (bit_number)
   {
      case 0: boolean = c & 0x80; break;
      case 1: boolean = c & 0x40; break;
      case 2: boolean = c & 0x20; break;
      case 3: boolean = c & 0x10; break;
      case 4: boolean = c & 0x08; break;
      case 5: boolean = c & 0x04; break;
      case 6: boolean = c & 0x02; break;
      case 7: boolean = c & 0x01; break;
   }

   return boolean;

}

void
Grib::Block::set_bit (const uint32_t position,
                      const uint8_t bit_number)
{

   uint8_t& c = buffer[position];

   switch (bit_number)
   {
      case 0: c |= 0x80; break;
      case 1: c |= 0x40; break;
      case 2: c |= 0x20; break;
      case 3: c |= 0x10; break;
      case 4: c |= 0x08; break;
      case 5: c |= 0x04; break;
      case 6: c |= 0x02; break;
      case 7: c |= 0x01; break;
   }

}

void
Grib::Block::unset_bit (const uint32_t position,
                        const uint8_t bit_number)
{

   uint8_t& c = buffer[position];

   switch (bit_number)
   {
      case 0: c &= 0x7f; break;
      case 1: c &= 0xbf; break;
      case 2: c &= 0xdf; break;
      case 3: c &= 0xef; break;
      case 4: c &= 0xf7; break;
      case 5: c &= 0xfb; break;
      case 6: c &= 0xfd; break;
      case 7: c &= 0xfe; break;
   }

}

uint32_t
Grib::Block::get_uint (FILE* file,
                       const uint32_t address,
                       const uint32_t n)
{
   Block block (n);
   fseek (file, address, SEEK_SET);
   fread (block.buffer, n, 1, file);
   return block.get_uint (0, n);
}

int32_t
Grib::Block::get_int (FILE* file,
                      const uint32_t address,
                      const uint32_t n)
{
   Block block (n);
   fseek (file, address, SEEK_SET);
   fread (block.buffer, n, 1, file);
   return block.get_int (0, n);
}

uint32_t
Grib::Block::get_uint (const uint32_t position,
                       const uint32_t n) const
{

   uint32_t integer = 0;
   uint8_t* pointer = buffer + position;

   for (uint32_t i = 0; i < n; i++)
   {
      if (i != 0) { integer <<= 8; }
      integer += pointer[i];
   }

   return integer;

}

int
Grib::Block::get_int (const uint32_t position,
                      const uint32_t n) const
{

   int integer = 0;
   uint8_t* pointer = buffer + position;

   for (uint32_t i = 0; i < n; i++)
   {
      if (i != 0) { integer <<= 8; }
      integer += (i == 0 ? (0x7f & pointer[i]) : pointer[i]);
   }

   if (*pointer & 0x80) { integer *= -1; }

   return integer;

}

uint32_t
Grib::Block::get_uint (const uint32_t position,
                       const uint8_t start_bit,
                       const uint8_t number_of_bits) const
{

   uint8_t* pointer = buffer + position;

   uint8_t number_of_bytes = number_of_bits / 8;
   if (number_of_bits & 0x7) { number_of_bytes++; }

   Block block (number_of_bytes, true);
   const uint8_t redundant_bits = number_of_bytes * 8 - number_of_bits;

   for (uint8_t i = 0; i < number_of_bits; i++)
   {

      uint8_t block_offset_byte  = (i + redundant_bits) / 8;
      uint8_t block_offset_bit   = (i + redundant_bits) & 0x7;
      uint8_t offset_byte = (i + start_bit) / 8;
      uint8_t offset_bit  = (i + start_bit) & 0x7;

      const bool b = get_bit (position + offset_byte, offset_bit);
      if (b) { block.set_bit (block_offset_byte, block_offset_bit); }

   }

   return block.get_uint (0, number_of_bytes);

}

int32_t
Grib::Block::get_int (const uint32_t position,
                      const uint8_t start_bit,
                      const uint8_t number_of_bits) const
{

   uint8_t* pointer = buffer + position;

   uint8_t number_of_bytes = number_of_bits / 8;
   if (number_of_bits & 0x7) { number_of_bytes++; }

   Block block (number_of_bytes, true);
   const uint8_t redundant_bits = number_of_bytes * 8 - number_of_bits;

   for (uint8_t i = 0; i < number_of_bits; i++)
   {

      uint8_t block_offset_byte  = (i + redundant_bits) / 8;
      uint8_t block_offset_bit   = (i + redundant_bits) & 0x7;
      uint8_t offset_byte = (i + start_bit) / 8;
      uint8_t offset_bit  = (i + start_bit) & 0x7;

      const bool b = get_bit (position + offset_byte, offset_bit);
      if (b) { block.set_bit (block_offset_byte, block_offset_bit); }

   }

   return block.get_int (0, number_of_bytes);

}

float
Grib::Block::get_float (const uint32_t position) const
{

   uint8_t* pointer = buffer + position;

   int A = (pointer[0] & 0x7f);
   int B = (pointer[1] << 16) + (pointer[2] << 8) + (pointer[3]);

   float number = pow (float (2), -24) * B * pow (float (16), A-64);
   if (pointer[0] & 0x80) { number *= -1; }

   return number;

}

bool
Grib::Block::operator== (const Block& block) const
{

   for (uint32_t i = 0; i < n; i++)
   {
      if (buffer[i] != block.buffer[i]) { return false; }
   }

   return true;

}

bool
Grib::Block::operator> (const Block& block) const
{

   for (uint32_t i = 0; i < n; i++)
   {
      if (buffer[i] > block.buffer[i]) { return true; }
      if (buffer[i] < block.buffer[i]) { return false; }
   }

   return false;

}

bool
Grib::Block::operator< (const Block& block) const
{

   for (uint32_t i = 0; i < n; i++)
   {
      if (buffer[i] < block.buffer[i]) { return true; }
      if (buffer[i] > block.buffer[i]) { return false; }
   }

   return false;

}

Grib::Section::Section (const Block& block,
                        const uint32_t offset,
                        const uint32_t length)
   : block (block),
     offset (offset),
     length (length)
{
}

uint32_t
Grib::Section::get_uint (const uint32_t position,
                         const uint32_t n) const
{
   return block.get_uint (offset + position, n);
}

int32_t
Grib::Section::get_int (const uint32_t position,
                        const uint32_t n) const
{
   return block.get_int (offset + position, n);
}

uint32_t
Grib::Section::get_uint (const uint32_t position,
                         const uint8_t start_bit,
                         const uint8_t number_of_bits) const
{
   return block.get_uint (offset + position, start_bit, number_of_bits);
}

int32_t
Grib::Section::get_int (const uint32_t position,
                        const uint8_t start_bit,
                        const uint8_t number_of_bits) const
{
   return block.get_int (offset + position, start_bit, number_of_bits);
}

float
Grib::Section::get_float (const uint32_t position) const
{
   return block.get_float (offset + position);
}

Grib::Pds::Level::Level (const Pds& pds)
   : Section (pds, 9, 3)
{
}

Grib::Pds::Forecast_Time::Forecast_Time (const Pds& pds)
   : Section (pds, 17, 6)
{
}

uint8_t
Grib::Pds::Forecast_Time::get_p1 () const
{
   return get_uint (1, 1);
}

uint8_t
Grib::Pds::Forecast_Time::get_p2 () const
{
   return get_uint (2, 1);
}

Grib::Pds::Pds (FILE* file,
                const uint32_t address)
   : Block (file, address)
{
}

uint8_t
Grib::Pds::get_center_id () const
{
   return buffer[4];
}

uint8_t
Grib::Pds::get_grid_id () const
{
   return buffer[6];
}

uint8_t
Grib::Pds::get_element_id () const
{
   return buffer[8];
}

bool
Grib::Pds::is_gds_present () const
{
   const uint8_t c = buffer[7];
   return (c & 0x80);
}

bool
Grib::Pds::is_bms_present () const
{
   const uint8_t c = buffer[7];
   return (c & 0x40);
}

Dtime
Grib::Pds::get_base_time () const
{

   uint8_t century = buffer[24];
   uint8_t year    = buffer[12];
   uint8_t month   = buffer[13];
   uint8_t day     = buffer[14];
   uint8_t hour    = buffer[15];
   uint8_t minute  = buffer[16];

   return Dtime ((century - 1) * 100 + year, month, day, hour, minute);

}

Grib::Pds::Level
Grib::Pds::get_level () const
{
   return Grib::Pds::Level (*this);
}

Grib::Pds::Forecast_Time
Grib::Pds::get_forecast_time () const
{
   return Grib::Pds::Forecast_Time (*this);
}

int16_t
Grib::Pds::get_D () const
{
   return get_int (26, 2);
}

Grib::Gds::Gds (FILE* file,
                const uint32_t address)
   : Block (file, address)
{
}

const Size_2D
Grib::Gds::get_size_2d () const
{

   const uint8_t grid_type = buffer[5];
   const bool regular_lat_long = (grid_type == 0);

   const uint16_t n_x = get_uint (6, 2);
   const uint16_t n_y = get_uint (8, 2);

   if (regular_lat_long) { return Size_2D (n_y, n_x); }
   else                  { return Size_2D (n_x, n_y); }

}

Grib::Bms::Bms (FILE* file,
                const uint32_t address)
   : Block (file, address)
{
}

Grib::Bds::Bds (FILE* file,
                const uint32_t address)
   : Block (file, address)
{
}

float
Grib::Bds::get_float (const uint32_t position,
                      const uint8_t start_bit,
                      const uint8_t number_of_bits,
                      const int16_t E,
                      const float R,
                      const int16_t D) const
{

   uint32_t X = get_uint (position, start_bit, number_of_bits);
   float datum = (R + (X * pow (float (2), E))) / pow (float (10), D);

   return datum;

}

Grib::Key::Key ()
   : Block (16)
{
}

Grib::Key::Key (const Pds& pds)
   : Block (16)
{
   memcpy (buffer, pds.buffer + 24, 1);     // Base Century
   memcpy (buffer + 1, pds.buffer + 12, 5); // Base YYMMDDHHMM
   memcpy (buffer + 6, pds.buffer + 17, 6); // Forecast_Hour
   memcpy (buffer + 12, pds.buffer + 8, 4); // Parameter + Layer
}

Grib::Key::Key (const Key& key)
   : Block (key)
{
}

ostream&
Grib::Key::operator << (ostream& out)
{
   for (Integer i = 0; i < 16; i++)
   {
      const uint8_t byte = buffer[i];
      out << byte;
      if (i != 15) { out << " "; }
   }
}

Grib::Header::Header (FILE* file,
                      uint32_t& address)
   : address (address)
{

   this->record_size = Block::get_uint (file, address + 4, 3);
   this->pds_ptr = new Pds (file, address + 8);
   this->gds_position = pds_ptr->n + 8;
   this->gds_ptr = new Gds (file, address + gds_position);

   const bool gds_present = pds_ptr->is_gds_present ();
   const uint32_t gds_address = address + gds_position;
   const uint32_t gds_size = (gds_present ? Block::get_size (file, gds_address) : 0);

   this->bms_position = gds_position + gds_size;
   const bool bms_present = pds_ptr->is_bms_present ();
   const uint32_t bms_address = address + bms_position;
   const uint32_t bms_size = (bms_present ? Block::get_size (file, bms_address) : 0);

   this->bds_position = bms_position + bms_size;
   //const uint32_t bds_address = address + bds_position;
   //const uint32_t bds_size = get_size (file, bbs_address);

   address += record_size;

}

Grib::Header::~Header ()
{
   delete pds_ptr;
   delete gds_ptr;
}

const Grib::Pds&
Grib::Header::get_pds () const
{
   return *pds_ptr;
}

const Grib::Gds&
Grib::Header::get_gds () const
{
   return *gds_ptr;
}

Grib::Grib (const Dstring& file_path)
   : file_path (file_path)
{

   FILE* file = get_input_file (file_path);
   uint32_t file_size = get_file_size (file_path);

   char buffer[5];
   buffer[4] = '\0';

   for (uint32_t offset = 0; offset < file_size; )
   {

      if ((offset + 4) >= file_size) { break; }

      fseek (file, offset, SEEK_SET);
      fread (buffer, 4, 1, file);

      if (strcmp ("GRIB", buffer) == 0)
      {
         Header* header_ptr = new Header (file, offset);
         const Grib::Key key (header_ptr->get_pds ());
         header_ptr_map.insert (make_pair (key, header_ptr));
      }
      else
      {
         offset++;
      }

   }

   fclose (file);

}

Grib::~Grib ()
{

   typedef map<Grib::Key, Grib::Header*>::iterator Iterator;

   for (Iterator iterator = header_ptr_map.begin ();
        iterator != header_ptr_map.end (); iterator++)
   {
      Grib::Header* header_ptr = (iterator->second);
      delete header_ptr;
   }

}

const map<Grib::Key, Grib::Header*>&
Grib::get_header_ptr_map () const
{
   return header_ptr_map;
}

Grib::Gds*
Grib::get_gds_ptr (FILE* file,
                   const Key& key) const
{
   const Header& header = *(header_ptr_map.find (key)->second);
   if (!header.get_pds ().is_gds_present ()) { return NULL; }
   return new Gds (file, header.address + header.gds_position);
}

Grib::Bms*
Grib::get_bms_ptr (FILE* file,
                   const Key& key) const
{
   const Header& header = *(header_ptr_map.find (key)->second);
   if (!header.get_pds ().is_bms_present ()) { return NULL; }
   return new Bms (file, header.address + header.bms_position);
}

Grib::Bds*
Grib::get_bds_ptr (FILE* file,
                   const Key& key) const
{
   const Header& header = *(header_ptr_map.find (key)->second);
   return new Bds (file, header.address + header.bds_position);
}

void
Grib::fill_data (Geodetic_Vector_Data_3D& gvd_3d,
                 const Integer element_index,
                 const Integer k,
                 const Key& key) const
{

   FILE* file = get_input_file (file_path);

   typedef map<Grib::Key, Grib::Header*>::const_iterator Iterator;
   Iterator iterator = header_ptr_map.find (key);

   //for (Iterator iterator = header_ptr_map.begin ();
   //     iterator != header_ptr_map.end (); iterator++)
   //{
   //   const Key& key = iterator->first;
   //   cout << key << endl;
   //   cout << k << endl;
   //   cout << endl;
   //}

   if (iterator == header_ptr_map.end ())
   {
      throw Nwp_Exception ("Grib::fill_data 3d Data Not Available");
   }

   const Header& header = *(iterator->second);
   const Pds& pds = header.get_pds ();
   const Gds gds (file, header.address + header.gds_position);
   const Bds bds (file, header.address + header.bds_position);

   const uint16_t n_x = gds.get_uint (6, 2);
   const uint16_t n_y = gds.get_uint (8, 2);

   const int16_t D = pds.get_D ();
   const int16_t E = bds.get_uint (4, 2);
   const float R = static_cast<Block> (bds).get_float (6);

   const uint8_t bits_per_datum = bds.buffer[10];
   const uint8_t grid_scan_flag = (gds.buffer[27] & 0xe0) >> 5;
   const bool c_notation = gds.buffer[27] & 0x20;

   uint8_t start_bit = 0;
   uint32_t position = 11;

   const bool swap = !c_notation;
   const uint32_t grid_size = n_x * n_y;

   for (int32_t serial_number = 0; serial_number < grid_size; serial_number++)
   {

      float datum = bds.get_float (position,
         start_bit, bits_per_datum, E, R, D);

      const uint8_t parameter = key.buffer[12];
      if (parameter == 52 || parameter == 60 ||
          (parameter >= 71 && parameter <=75))
      {
         datum /= 100;
      }

      Integer i, j;
      switch (grid_scan_flag)
      {

         case 0:
            i = serial_number % n_x;
            j = n_y - (serial_number / n_x) - 1;
            break;

         case 2:
            i = serial_number % n_x;
            j = serial_number / n_x;
            break;

         case 4:
            i = n_x - (serial_number % n_x) - 1;
            j = n_y - (serial_number / n_x) - 1;
            break;

         case 6:
            i = n_x - (serial_number % n_x) - 1;
            j = serial_number / n_x;
            break;

         case 1:
            i = (serial_number / n_y);
            j = n_y - (serial_number % n_y) - 1;
            break;

         case 3:
            i = serial_number / n_y;
            j = serial_number % n_y;
            break;

         case 5:
            i = n_x - (serial_number / n_y) - 1;
            j = n_y - (serial_number % n_y) - 1;
            break;

         case 7:
            i = n_x - (serial_number / n_y) - 1;
            j = serial_number % n_y;
            break;

      }

      if (swap) { gvd_3d.set_datum (element_index, k, j, i, datum); }
      else      { gvd_3d.set_datum (element_index, k, i, j, datum); }

      position += (start_bit + bits_per_datum) / 8;
      start_bit = (start_bit + bits_per_datum) & 0x7;

   }

   fclose (file);

}

void
Grib::fill_data (Geodetic_Vector_Data_2D& gvd_2d,
                 const Integer element_index,
                 const Key& key) const
{

   FILE* file = get_input_file (file_path);

   typedef map<Grib::Key, Grib::Header*>::const_iterator Iterator;
   Iterator iterator = header_ptr_map.find (key);

   if (iterator == header_ptr_map.end ())
   {
      throw Nwp_Exception ("Grib::fill_data 2d Data Not Available");
   }

   const Header& header = *(iterator->second);
   const Pds& pds = header.get_pds ();
   const Gds gds (file, header.address + header.gds_position);
   const Bds bds (file, header.address + header.bds_position);

   const uint16_t n_x = gds.get_uint (6, 2);
   const uint16_t n_y = gds.get_uint (8, 2);

   const int16_t D = pds.get_D ();
   const int16_t E = bds.get_int (4, 2);
   const float R = static_cast<Block> (bds).get_float (6);

   const uint8_t bits_per_datum = bds.buffer[10];
   const uint8_t grid_scan_flag = (gds.buffer[27] & 0xe0) >> 5;
   const bool c_notation = gds.buffer[27] & 0x20;

   uint8_t start_bit = 0;
   uint32_t position = 11;

   const bool swap = !c_notation;
   const uint32_t grid_size = n_x * n_y;

   for (int32_t serial_number = 0; serial_number < grid_size; serial_number++)
   {

      float datum = bds.get_float (position,
         start_bit, bits_per_datum, E, R, D);

      const uint8_t parameter = key.buffer[12];
      if (parameter == 52 || parameter == 60 ||
          (parameter >= 71 && parameter <= 75))
      {
         datum /= 100;
      }

      Integer i, j;
      switch (grid_scan_flag)
      {

         case 0:
            i = serial_number % n_x;
            j = n_y - (serial_number / n_x) - 1;
            break;

         case 2:
            i = serial_number % n_x;
            j = serial_number / n_x;
            break;

         case 4:
            i = n_x - (serial_number % n_x) - 1;
            j = n_y - (serial_number / n_x) - 1;
            break;

         case 6:
            i = n_x - (serial_number % n_x) - 1;
            j = serial_number / n_x;
            break;

         case 1:
            i = (serial_number / n_y);
            j = n_y - (serial_number % n_y) - 1;
            break;

         case 3:
            i = serial_number / n_y;
            j = serial_number % n_y;
            break;

         case 5:
            i = n_x - (serial_number / n_y) - 1;
            j = n_y - (serial_number % n_y) - 1;
            break;

         case 7:
            i = n_x - (serial_number / n_y) - 1;
            j = serial_number % n_y;
            break;

      }

      if (swap) { gvd_2d.set_datum (element_index, j, i, datum); }
      else      { gvd_2d.set_datum (element_index, i, j, datum); }

      position += (start_bit + bits_per_datum) / 8;
      start_bit = (start_bit + bits_per_datum) & 0x7;

   }

   fclose (file);

} 

Grib2::Block::Block ()
{
}

Grib2::Block::Block (const uint32_t n,
                     const bool set_zero)
   : n (n)
{
   buffer = new uint8_t[n];
   if (set_zero) { for (uint16_t i = 0; i < n; i++) { buffer[i] = 0; } }
}

Grib2::Block::Block (FILE* file,
                     uint32_t& address)
{

   n = get_size (file, address);

   buffer = new uint8_t[n];
   fseek (file, address, SEEK_SET);
   fread (buffer, n, 1, file);

   address += n;

}

Grib2::Block::Block (const Block& block)
   : n (block.n)
{
   buffer = new uint8_t[n];
   memcpy (buffer, block.buffer, n);
}

Grib2::Block::~Block ()
{
   delete[] buffer;
}

uint32_t
Grib2::Block::get_size (FILE* file,
                        const uint32_t address)
{
   return get_uint (file, address, 4);
}

bool
Grib2::Block::get_bit (const uint32_t position,
                       const uint8_t bit_number) const
{

   bool boolean;
   const uint8_t c = buffer[position];

   switch (bit_number)
   {
      case 0: boolean = c & 0x80; break;
      case 1: boolean = c & 0x40; break;
      case 2: boolean = c & 0x20; break;
      case 3: boolean = c & 0x10; break;
      case 4: boolean = c & 0x08; break;
      case 5: boolean = c & 0x04; break;
      case 6: boolean = c & 0x02; break;
      case 7: boolean = c & 0x01; break;
   }

   return boolean;

}

void
Grib2::Block::set_bit (const uint32_t position,
                       const uint8_t bit_number)
{

   uint8_t& c = buffer[position];

   switch (bit_number)
   {
      case 0: c |= 0x80; break;
      case 1: c |= 0x40; break;
      case 2: c |= 0x20; break;
      case 3: c |= 0x10; break;
      case 4: c |= 0x08; break;
      case 5: c |= 0x04; break;
      case 6: c |= 0x02; break;
      case 7: c |= 0x01; break;
   }

}

void
Grib2::Block::unset_bit (const uint32_t position,
                         const uint8_t bit_number)
{

   uint8_t& c = buffer[position];

   switch (bit_number)
   {
      case 0: c &= 0x7f; break;
      case 1: c &= 0xbf; break;
      case 2: c &= 0xdf; break;
      case 3: c &= 0xef; break;
      case 4: c &= 0xf7; break;
      case 5: c &= 0xfb; break;
      case 6: c &= 0xfd; break;
      case 7: c &= 0xfe; break;
   }

}

uint64_t
Grib2::Block::get_uint (FILE* file,
                        const uint32_t address,
                        const uint32_t n)
{
   Block block (n);
   fseek (file, address, SEEK_SET);
   fread (block.buffer, n, 1, file);
   return block.get_uint (0, n);
}

int64_t
Grib2::Block::get_int (FILE* file,
                       const uint32_t address,
                       const uint32_t n)
{
   Block block (n);
   fseek (file, address, SEEK_SET);
   fread (block.buffer, n, 1, file);
   return block.get_int (0, n);
}

uint64_t
Grib2::Block::get_uint (const uint32_t position,
                        const uint32_t n) const
{

   uint64_t integer = 0;
   uint8_t* pointer = buffer + position;

   for (uint32_t i = 0; i < n; i++)
   {
      if (i != 0) { integer <<= 8; }
      integer += pointer[i];
   }

   return integer;

}

int64_t
Grib2::Block::get_int (const uint32_t position,
                       const uint32_t n) const
{

   int integer = 0;
   uint8_t* pointer = buffer + position;

   for (uint32_t i = 0; i < n; i++)
   {
      if (i != 0) { integer <<= 8; }
      integer += (i == 0 ? (0x7f & pointer[i]) : pointer[i]);
   }

   if (*pointer & 0x80) { integer *= -1; }

   return integer;

}

uint32_t
Grib2::Block::get_uint (const uint32_t position,
                        const uint8_t start_bit,
                        const uint8_t number_of_bits) const
{

   uint8_t* pointer = buffer + position;

   uint8_t number_of_bytes = number_of_bits / 8;
   if (number_of_bits & 0x7) { number_of_bytes++; }

   Block block (number_of_bytes, true);
   const uint8_t redundant_bits = number_of_bytes * 8 - number_of_bits;

   for (uint8_t i = 0; i < number_of_bits; i++)
   {

      uint8_t block_offset_byte  = (i + redundant_bits) / 8;
      uint8_t block_offset_bit   = (i + redundant_bits) & 0x7;
      uint8_t offset_byte = (i + start_bit) / 8;
      uint8_t offset_bit  = (i + start_bit) & 0x7;

      const bool b = get_bit (position + offset_byte, offset_bit);
      if (b) { block.set_bit (block_offset_byte, block_offset_bit); }

   }

   return block.get_uint (0, number_of_bytes);

}

int32_t
Grib2::Block::get_int (const uint32_t position,
                       const uint8_t start_bit,
                       const uint8_t number_of_bits) const
{

   uint8_t* pointer = buffer + position;

   uint8_t number_of_bytes = number_of_bits / 8;
   if (number_of_bits & 0x7) { number_of_bytes++; }

   Block block (number_of_bytes, true);
   const uint8_t redundant_bits = number_of_bytes * 8 - number_of_bits;

   for (uint8_t i = 0; i < number_of_bits; i++)
   {

      uint8_t block_offset_byte  = (i + redundant_bits) / 8;
      uint8_t block_offset_bit   = (i + redundant_bits) & 0x7;
      uint8_t offset_byte = (i + start_bit) / 8;
      uint8_t offset_bit  = (i + start_bit) & 0x7;

      const bool b = get_bit (position + offset_byte, offset_bit);
      if (b) { block.set_bit (block_offset_byte, block_offset_bit); }

   }

   return block.get_int (0, number_of_bytes);

}

float
Grib2::Block::get_float (const uint32_t position) const
{

   uint8_t* pointer = buffer + position;

   int A = (pointer[0] & 0x7f);
   int B = (pointer[1] << 16) + (pointer[2] << 8) + (pointer[3]);

   float number = pow (float (2), -24) * B * pow (float (16), A-64);
   if (pointer[0] & 0x80) { number *= -1; }

   return number;

}

bool
Grib2::Block::operator== (const Block& block) const
{

   for (uint32_t i = 0; i < n; i++)
   {
      if (buffer[i] != block.buffer[i]) { return false; }
   }

   return true;

}

bool
Grib2::Block::operator> (const Block& block) const
{

   for (uint32_t i = 0; i < n; i++)
   {
      if (buffer[i] > block.buffer[i]) { return true; }
      if (buffer[i] < block.buffer[i]) { return false; }
   }

   return false;

}

bool
Grib2::Block::operator< (const Block& block) const
{

   for (uint32_t i = 0; i < n; i++)
   {
      if (buffer[i] < block.buffer[i]) { return true; }
      if (buffer[i] > block.buffer[i]) { return false; }
   }

   return false;

}

Grib2::Section::Section (const Block& block,
                        const uint32_t offset,
                        const uint32_t length)
   : block (block),
     offset (offset),
     length (length)
{
}

uint32_t
Grib2::Section::get_uint (const uint32_t position,
                         const uint32_t n) const
{
   return block.get_uint (offset + position, n);
}

int32_t
Grib2::Section::get_int (const uint32_t position,
                        const uint32_t n) const
{
   return block.get_int (offset + position, n);
}

uint32_t
Grib2::Section::get_uint (const uint32_t position,
                          const uint8_t start_bit,
                          const uint8_t number_of_bits) const
{
   return block.get_uint (offset + position, start_bit, number_of_bits);
}

int32_t
Grib2::Section::get_int (const uint32_t position,
                         const uint8_t start_bit,
                         const uint8_t number_of_bits) const
{
   return block.get_int (offset + position, start_bit, number_of_bits);
}

float
Grib2::Section::get_float (const uint32_t position) const
{
   return block.get_float (offset + position);
}

Grib2::Block_1::Base_Time::Base_Time (const Block_1& block_1)
   : Section (block_1, 12, 7)
{
}

Grib2::Block_1::Block_1 (FILE* file,
                         uint32_t& address)
   : Block (file, address)
{
}

Dtime
Grib2::Block_1::get_reference_time () const
{

   uint16_t century = get_uint (12, 2);
   uint8_t year    = buffer[13];
   uint8_t month   = buffer[14];
   uint8_t day     = buffer[15];
   uint8_t hour    = buffer[16];
   uint8_t minute  = buffer[17];
   uint8_t second  = buffer[18];

   return Dtime (century * 100 + year, month, day, hour, minute, second);

}

Grib2::Block_1::Base_Time
Grib2::Block_1::get_base_time () const
{
   return Grib2::Block_1::Base_Time (*this);
}

Grib2::Block_2::Block_2 (FILE* file,
                         uint32_t& address)
   : Block (file, address)
{
}

uint16_t
Grib2::Block_3::get_template_number () const
{
   return get_uint (12, 2);
}

uint8_t
Grib2::Block_3::get_scanning_mode () const
{
   return buffer[71];
}

Grib2::Block_3::Block_3 (FILE* file,
                         uint32_t& address)
   : Block (file, address)
{
}

uint32_t
Grib2::Block_3::get_number_of_points () const
{
   return get_uint (6, 4);
}

const Size_2D
Grib2::Block_3::get_size_2d () const
{

   const uint16_t template_number = get_template_number ();

   switch (template_number)
   {

      case 0:
      {
         const uint32_t n_x = get_uint (34, 4);
         const uint32_t n_y = get_uint (30, 4);
         return Size_2D (n_x, n_y);
      }

      default:
         throw Exception ("Not implemented yet");
   }

}

const uint8_t
Grib2::Block_3::get_scan_flag () const
{

   const uint16_t template_number = get_template_number ();

   switch (template_number)
   {
      default: throw Exception ("Not implemented yet");
      case 0: return ((buffer[71] & 0xf0) >> 4);
      case 10: return ((buffer[59] & 0xf0) >> 4);
      case 20: return ((buffer[64] & 0xf0) >> 4);
      case 30: return ((buffer[64] & 0xf0) >> 4);
      case 31: return ((buffer[64] & 0xf0) >> 4);
      case 40: return ((buffer[71] & 0xf0) >> 4);
      case 41: return ((buffer[71] & 0xf0) >> 4);
      case 42: return ((buffer[71] & 0xf0) >> 4);
      case 43: return ((buffer[71] & 0xf0) >> 4);
      case 44: return ((buffer[71] & 0xf0) >> 4);
   }

}

Grib2::Block_4::Parameter::Parameter (const Block_4& block_4)
   : Section (block_4, 9, 2)
{
}

Grib2::Block_4::Level::Level (const Block_4& block_4)
   : Section (block_4, 22, 12)
{
}

uint8_t
Grib2::Block_4::Level::get_first_level_type () const
{
   return block.buffer[offset];
}

uint32_t
Grib2::Block_4::Level::get_first_level () const
{
   int8_t scale_factor = block.buffer[offset + 1];
   uint32_t value = get_uint (2, 4);
   uint32_t scale = uint32_t (round (pow (10.0, double (scale_factor))));
   return value * scale;
}

uint8_t
Grib2::Block_4::Level::get_second_level_type () const
{
   return block.buffer[offset + 6];
}

uint32_t
Grib2::Block_4::Level::get_second_level () const
{
   int8_t scale_factor = block.buffer[offset + 7];
   uint32_t value = get_uint (8, 4);
   uint32_t scale = uint32_t (round (pow (10.0, double (scale_factor))));
   return value * scale;
}

Grib2::Block_4::Forecast_Time::Forecast_Time (const Block_4& block_4)
   : Section (block_4, 17, 5)
{
}

uint16_t
Grib2::Block_4::get_template_number () const
{
   return get_uint (7, 2);
}

Grib2::Block_4::Block_4 (FILE* file,
                         uint32_t& address)
   : Block (file, address)
{
}

Grib2::Block_4::Parameter
Grib2::Block_4::get_parameter () const
{
   return Grib2::Block_4::Parameter (*this);
}

Grib2::Block_4::Level
Grib2::Block_4::get_level () const
{
   return Grib2::Block_4::Level (*this);
}

Grib2::Block_4::Forecast_Time
Grib2::Block_4::get_forecast_time () const
{
   return Grib2::Block_4::Forecast_Time (*this);
}

bool
Grib2::Block_4::parameter_is_percentage () const
{

   const uint8_t octet_9 = buffer[9];
   const uint8_t octet_10 = buffer[10];

   if (octet_9 == 1)
   {
      if (octet_10 == 1 || octet_10 == 27 || octet_10 == 39 ||
          octet_10 == 42 || octet_10 == 198 || octet_10 == 201 ||
          octet_10 == 242)
      {
         return true;
      }
   }
   else
   if (octet_9 == 6)
   {
      if (octet_10 == 1 || octet_10 == 2 || octet_10 == 3 ||
          octet_10 == 4 || octet_10 == 5 || octet_10 == 7 ||
          octet_10 == 14 || octet_10 == 22 || octet_10 == 25 ||
          octet_10 == 192)
      {
         return true;
      }
   }

   return false;

}

Grib2::Block_5::Block_5 (FILE* file,
                         uint32_t& address)
   : Block (file, address)
{
}

uint16_t
Grib2::Block_5::get_template_number () const
{
   return get_uint (9, 2);
}

uint32_t
Grib2::Block_5::get_number_of_points () const
{
   return get_uint (5, 4);
}

float
Grib2::Block_5::get_reference_value () const
{

   const uint8_t template_number = get_template_number ();

   switch (template_number)
   {
      default:
         throw Exception ("Not implemented");
      case 0:
      case 1:
      case 2:
      case 3:
      case 40:
      case 41:
      case 50:
      case 61:
      {
         float value;
         memcpy (&value, buffer + 11, sizeof (float));
#ifndef WORDS_BIGENDIAN
         swap_endian (&value, sizeof (float));
#endif
         return value;
      }
   }

}

int16_t
Grib2::Block_5::get_binary_scale_factor () const
{

   const uint8_t template_number = get_template_number ();

   switch (template_number)
   {
      default:
         throw Exception ("Not implemented");
      case 0:
      case 1:
      case 2:
      case 3:
      case 40:
      case 41:
      case 50:
      case 61:
         return get_int (15, 2);
   }

}

int16_t
Grib2::Block_5::get_decimal_scale_factor () const
{

   const uint8_t template_number = get_template_number ();

   switch (template_number)
   {
      default:
         throw Exception ("Not implemented");
      case 0:
      case 1:
      case 2:
      case 3:
      case 40:
      case 41:
      case 50:
      case 61:
         return get_int (17, 2);
   }

}

uint8_t
Grib2::Block_5::get_bits_per_datum () const
{

   const uint8_t template_number = get_template_number ();

   switch (template_number)
   {
      default:
         throw Exception ("Not implemented");
      case 0:
      case 1:
      case 2:
      case 3:
      case 40:
      case 41:
      case 50:
      case 61:
         return buffer[19];
   }

}

Grib2::Block_6::Block_6 (FILE* file,
                         uint32_t& address)
   : Block (file, address)
{
}

uint8_t
Grib2::Block_6::get_bitmap_indicator () const
{
   return buffer[5];
}

bool
Grib2::Block_6::bitmap_present () const
{
   return (get_bitmap_indicator () == 0);
}

Grib2::Block_7::Block_7 (FILE* file,
                         uint32_t& address)
   : Block (file, address)
{
}

Grib2::Data::Data (const uint32_t n)
   : n (n),
     buffer (new float[n])
{
}

Grib2::Data::~Data ()
{
   delete[] buffer;
}

const uint32_t
Grib2::Data::get_n () const
{
   return n;
}

float*
Grib2::Data::get_buffer () const
{
   return buffer;
}

void
Grib2::Data::simple (const Block_5& block_5,
                     const Block_7& block_7)
{

   uint8_t start_bit = 0;
   uint32_t position = 0;

   const uint8_t bits_per_datum = block_5.get_bits_per_datum ();
   const int16_t D = block_5.get_decimal_scale_factor ();
   const int16_t E = block_5.get_binary_scale_factor ();
   const float R = block_5.get_reference_value ();
   const uint32_t n = block_5.get_number_of_points ();

   for (uint32_t index = 0; index < n; index++)
   {
      const uint32_t X = block_7.get_uint (position, start_bit, bits_per_datum);
      buffer[index] = (R + (X * pow (float (2), E))) / pow (float (10), D);
      position += (start_bit + bits_per_datum) / 8;
      start_bit = (start_bit + bits_per_datum) & 0x7;
   }

}

void
Grib2::Data::jpeg2000 (const Block_5& block_5,
                       const Block_7& block_7)
{

   const uint8_t bits_per_datum = block_5.get_bits_per_datum ();
   const int16_t D = block_5.get_decimal_scale_factor ();
   const int16_t E = block_5.get_binary_scale_factor ();
   const float R = block_5.get_reference_value ();

   //jas_init();
   char* encoded_buffer = (char*)block_7.buffer + 5;
   const uint32_t encoded_buffer_n = block_7.n - 5;
   jas_stream_t* stream = jas_stream_memopen (encoded_buffer, encoded_buffer_n);

   char* opts = 0;
   jas_image_t* image = jpc_decode (stream, opts);
   if (image == 0) { throw Exception ("JPEG2000 decode error"); }
   const int number_of_components = image->numcmpts_;
   const bool grayscale = (number_of_components == 1);
   if (!grayscale) { throw Exception ("JPEG2000 Color image"); }

   const Integer image_height = jas_image_height (image);
   const Integer image_width = jas_image_width (image);
   jas_matrix_t* matrix_data = jas_matrix_create (image_height, image_width);
   jas_image_readcmpt (image, 0, 0, 0, image_width, image_height, matrix_data);

   const uint32_t n = block_5.get_number_of_points ();
   const jas_image_cmpt_t* component_ptr = image->cmpts_[0];

   const Integer w = component_ptr->width_;
   const Integer h = component_ptr->height_;

   for (Integer row = 0; row < h; row++) 
   {
      for (Integer column = 0; column < w; column++)
      {
         const Integer index = row * w + column;
         const uint32_t X = matrix_data->rows_[row][column];
         float datum = (R + (X * powf (float (2), E))) / powf (float (10), D);
         buffer[index] = datum;
      }
   }

   jas_matrix_destroy (matrix_data);
   jas_stream_close (stream);
   jas_image_destroy (image);

}

Grib2::Key::Key ()
   : Block (7 + 5 + 2 + 12)
{
}

Grib2::Key::Key (const Header& header)
   : Block (7 + 5 + 2 + 12)
{

   const Grib2::Block_4::Level l = header.block_4.get_level ();
   const Grib2::Block_1::Base_Time bt = header.block_1.get_base_time ();
   const Grib2::Block_4::Parameter p = header.block_4.get_parameter ();
   const Grib2::Block_4::Forecast_Time ft = header.block_4.get_forecast_time ();

   memcpy (buffer + 0, bt.block.buffer + bt.offset, bt.length);
   memcpy (buffer + 7, ft.block.buffer + ft.offset, ft.length);
   memcpy (buffer + 12, p.block.buffer + p.offset, p.length);
   memcpy (buffer + 14, l.block.buffer + l.offset, l.length);

}

Grib2::Key::Key (const Key& key)
   : Block (key)
{
}

Grib2::Header::Header (const Block_1& block_1,
                       const Block_3& block_3,
                       const Block_4& block_4,
                       const Block_5& block_5,
                       const Block_6& block_6,
                       const uint32_t block_7_address)
   : block_1 (block_1),
     block_3 (block_3),
     block_4 (block_4),
     block_5 (block_5),
     block_6 (block_6),
     block_7_address (block_7_address)
{
}

const Grib2::Header& 
Grib2::get_header (const Key& key) const
{

   typedef map<Grib2::Key, Grib2::Header*>::const_iterator Iterator;
   Iterator iterator = header_ptr_map.find (key);

   //for (Iterator i = header_ptr_map.begin ();
   //     i != header_ptr_map.end (); i++)
   //{
   //   const Grib2::Key& k = i->first;
   //   cout << key << "| " << k << endl;
   //}

   if (iterator == header_ptr_map.end ())
   {
      throw Nwp_Exception ("Data Not Available");
   }

   const Header& header = *(iterator->second);
   return header;

}

Grib2::Grib2 (const Dstring& file_path)
   : file_path (file_path)
{

   FILE* file = get_input_file (file_path);
   uint32_t file_size = get_file_size (file_path);

   char buffer[5];
   buffer[4] = '\0';

Integer count = 0;
   for (uint32_t offset = 0; offset < file_size; )
   {

      const uint32_t address = offset;
      if ((offset + 4) >= file_size) { break; }

      fseek (file, offset, SEEK_SET);
      fread (buffer, 4, 1, file);

      const bool grib_head = (strcmp ("GRIB", buffer) == 0);

      if (!grib_head) { offset++; continue; }
count++;

      const uint32_t this_offset = offset;
      uint64_t grib_size = Block::get_uint (file, offset + 8, 8);

      offset += 16; // Skipping Block_0
      const Block_1 block_1 (file, offset);
      const uint8_t next_section = Block::get_uint (file, offset + 4, 1);
      const bool block_2_present = (next_section == 2);

      if (block_2_present)
      {
         const uint32_t block_2_n = Block::get_uint (file, offset, 4);
         offset += block_2_n;
      }

      Block_3* block_3_ptr = new Block_3 (file, offset);
      Block_4* block_4_ptr = new Block_4 (file, offset);
      Block_5* block_5_ptr = new Block_5 (file, offset);
      Block_6* block_6_ptr = new Block_6 (file, offset);
      const uint32_t block_7_address = offset;

      Header* header_ptr = new Header (block_1, *block_3_ptr,
         *block_4_ptr, *block_5_ptr, *block_6_ptr, block_7_address);
      const Grib2::Key key (*header_ptr);
      header_ptr_map.insert (make_pair (key, header_ptr));

      const uint32_t block_7_n = Block::get_uint (file, offset, 4);
      offset += block_7_n; // Block 7

   const uint8_t scan_flag = (block_3_ptr->get_scan_flag ());

const bool bit_0 = scan_flag & 0x80;
const bool bit_1 = scan_flag & 0x40;
const bool bit_2 = scan_flag & 0x20;
const bool bit_3 = scan_flag & 0x10;
const Dstring b0 (bit_0 ? "yes" : "no");
const Dstring b1 (bit_1 ? "yes" : "no");
const Dstring b2 (bit_2 ? "yes" : "no");
const Dstring b3 (bit_3 ? "yes" : "no");
      while (offset - this_offset + 4 < grib_size)
      {

         const uint32_t next_block_n = Block::get_uint (file, offset, 4);
         const uint8_t block_number = Block::get_uint (file, offset + 4, 1);

         Block_3* b_3_ptr = NULL;
         Block_4* b_4_ptr = NULL;

         // intended no breaks here in the switch construct
         //    presence of block {n}, implies presese of block_{n+1}
         switch (block_number)
         {
            default: throw Exception ("Grib2 Coding Error");
            case 2: offset += Block::get_uint (file, offset, 4);
            case 3: b_3_ptr = new Block_3 (file, offset);
            case 4: b_4_ptr = new Block_4 (file, offset);
         }

         Block_5* b_5_ptr = new Block_5 (file, offset);
         Block_6* b_6_ptr = new Block_6 (file, offset);

         const uint32_t block_7_address = offset;
         Header* header_ptr = new Header (block_1,
            (b_3_ptr == NULL ? *block_3_ptr : *b_3_ptr),
            (b_4_ptr == NULL ? *block_4_ptr : *b_4_ptr),
            (b_5_ptr == NULL ? *block_5_ptr : *b_5_ptr),
            (b_6_ptr == NULL ? *block_6_ptr : *b_6_ptr),
            block_7_address);
         const Grib2::Key key (*header_ptr);
         header_ptr_map.insert (make_pair (key, header_ptr));

         if (b_3_ptr != NULL) { delete b_3_ptr; }
         if (b_4_ptr != NULL) { delete b_4_ptr; }
         if (b_5_ptr != NULL) { delete b_5_ptr; }
         if (b_6_ptr != NULL) { delete b_6_ptr; }

         const uint32_t block_7_n = Block::get_uint (file, offset, 4);
         offset += block_7_n; // Block 7

      }

      delete block_3_ptr;
      delete block_4_ptr;
      delete block_5_ptr;
      delete block_6_ptr;

      offset += 4;

   }

   fclose (file);

   const Header& header = *(header_ptr_map.begin ()->second);
   const Integer n = header.block_5.get_number_of_points ();
   data_ptr = new Data (n);

}

Grib2::~Grib2 ()
{

   delete data_ptr;

   typedef map<Grib2::Key, Grib2::Header*>::iterator Iterator;

   for (Iterator iterator = header_ptr_map.begin ();
        iterator != header_ptr_map.end (); iterator++)
   {
      Grib2::Header* header_ptr = (iterator->second);
      delete header_ptr;
   }

}

const map<Grib2::Key, Grib2::Header*>&
Grib2::get_header_ptr_map () const
{
   return header_ptr_map;
}

void
Grib2::fill_data (Geodetic_Vector_Data_3D& gvd_3d,
                  const Integer element_index,
                  const Integer k,
                  const Key& key) const
{
/*
*/

   const Header& header = get_header (key);
   const Grib2::Key kkk (header);

   const Block_3& block_3 = header.block_3;
   const Block_5& block_5 = header.block_5;

   FILE* file = get_input_file (file_path);
   uint32_t block_7_address = header.block_7_address;
   const Block_7 block_7 (file, block_7_address);
   fclose (file);

   const uint32_t n = block_5.get_number_of_points ();
//   Data& data = *data_ptr;
   Data data (n);
   switch (block_5.get_template_number ())
   {
      default: throw Exception ("Packing Type not supported.");
      case 0:  data.simple (block_5, block_7); break;
      case 40: data.jpeg2000 (block_5, block_7); break;
   }

   float* buffer = data.get_buffer ();
   const Size_2D& size_2d = block_3.get_size_2d ();
   const uint8_t scan_flag = (block_3.get_scan_flag ());

   const bool c_notation = (scan_flag & 0x2);
   const bool swap = !c_notation;

   for (uint32_t index = 0; index < n; index++)
   {

      const Real datum = buffer[index];

      Integer i, j;
      switch (scan_flag >> 1)
      {

         case 0:
            i = size_2d.i - index / size_2d.j - 1;
            j = index % size_2d.j;
            break;

         case 2: // untested
            i = index % size_2d.i;
            j = index / size_2d.i;
            break;

         case 4: // untested
            i = size_2d.i - (index % size_2d.i) - 1;
            j = size_2d.j - (index / size_2d.i) - 1;
            break;

         case 6: // untested
            i = size_2d.i - (index % size_2d.i) - 1;
            j = index / size_2d.i;
            break;

         case 1: // untested
            i = (index / size_2d.j);
            j = size_2d.j - (index % size_2d.j) - 1;
            break;

         case 3: // untested
            i = index / size_2d.j;
            j = index % size_2d.j;
            break;

         case 5: // untested
            i = size_2d.i - (index / size_2d.j) - 1;
            j = size_2d.j - (index % size_2d.j) - 1;
            break;

         case 7: // untested
            i = size_2d.i - (index / size_2d.j) - 1;
            j = index % size_2d.j;
            break;

      }

//      if (swap) { std::swap (i, j); }
      gvd_3d.set_datum (element_index, k, i, j, datum);

   }

}

void
Grib2::fill_data (Geodetic_Vector_Data_2D& gvd_2d,
                  const Integer element_index,
                  const Key& key) const
{

   const Header& header = get_header (key);
   const Grib2::Key kkk (header);

   const Block_3& block_3 = header.block_3;
   const Block_4& block_4 = header.block_4;
   const Block_5& block_5 = header.block_5;

   FILE* file = get_input_file (file_path);
   uint32_t block_7_address = header.block_7_address;
   const Block_7 block_7 (file, block_7_address);
   fclose (file);

   const uint32_t n = block_5.get_number_of_points ();
//   Data& data = *data_ptr;
   Data data (n);
   switch (block_5.get_template_number ())
   {
      default: throw Exception ("Packing Type not supported.");
      case 0:  data.simple (block_5, block_7); break;
      case 40: data.jpeg2000 (block_5, block_7); break;
   }

   float* buffer = data.get_buffer ();
   const Size_2D& size_2d = block_3.get_size_2d ();
   const uint8_t scan_flag = (block_3.get_scan_flag ());

   const bool c_notation = (scan_flag & 0x2);
   const bool swap = !c_notation;

   const Real multiplier = (block_4.parameter_is_percentage () ? 0.01 : 1);
   for (uint32_t index = 0; index < n; index++)
   {

      const Real datum = buffer[index] * multiplier;

      Integer i, j;
      switch (scan_flag >> 1)
      {

         case 0:
            i = size_2d.i - index / size_2d.j - 1;
            j = index % size_2d.j;
            break;

         case 2: // untested
            i = index % size_2d.i;
            j = index / size_2d.i;
            break;

         case 4: // untested
            i = size_2d.i - (index % size_2d.i) - 1;
            j = size_2d.j - (index / size_2d.i) - 1;
            break;

         case 6: // untested
            i = size_2d.i - (index % size_2d.i) - 1;
            j = index / size_2d.i;
            break;

         case 1: // untested
            i = (index / size_2d.j);
            j = size_2d.j - (index % size_2d.j) - 1;
            break;

         case 3: // untested
            i = index / size_2d.j;
            j = index % size_2d.j;
            break;

         case 5: // untested
            i = size_2d.i - (index / size_2d.j) - 1;
            j = size_2d.j - (index % size_2d.j) - 1;
            break;

         case 7: // untested
            i = size_2d.i - (index / size_2d.j) - 1;
            j = index % size_2d.j;
            break;

      }

//      if (swap) { std::swap (i, j); }
      gvd_2d.set_datum (element_index, i, j, datum);

   }

} 

namespace denise
{

   ostream&
   operator << (ostream &out,
                const Grib::Key& key)
   {
      for (Integer i = 0; i < 16; i++)
      {
         out << Integer (key.buffer[i]) << " ";
      }
      return out;
   }

   ostream&
   operator << (ostream &out,
                const Grib2::Key& key)
   {
      for (Integer i = 0; i < 26; i++)
      {
         Dstring sep (" ");
         if (i == 6 || i == 11 || i == 13) { sep = " : "; }
         out << Integer (key.buffer[i]) << sep;
      }
      return out;
   }

};

Access::Data_3D::Data_3D (const vector<Met_Element>& met_element_vector,
                          const Key& key)
   : Nwp::Data_3D (met_element_vector, key)
{
}

Real
Access::Data_3D::evaluate (const Met_Element element,
                           const Real p,
                           const Real latitude,
                           const Real longitude,
                           const Evaluate_Op evaluate_op) const
{

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& RH = denise::RELATIVE_HUMIDITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real rh = Nwp::Data_3D::evaluate (RH, p, latitude, longitude);
         const Thermo_Medium thermo_medium = (t < 0 ? ICE : WATER);
         return Moisture::get_t_d (t, rh, WATER);
      }

      case denise::OMEGA:
      {
         const Met_Element& T = TEMPERATURE;
         const Met_Element& W = VERTICAL_VELOCITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real w = Nwp::Data_3D::evaluate (W, p, latitude, longitude);
         const Real rho = p / (R_d * t);
         return -rho * g * w;
      }

   }

   return Nwp::Data_3D::evaluate (element, p,
      latitude, longitude, evaluate_op);

}

Grib::Key
Access::get_grib_key (const Key& key,
                      const Met_Element met_element,
                      const Level& level) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   Grib::Key grib_key;
   set_grib_key (grib_key, met_element, base_time, forecast_hour);
   set_grib_key (grib_key, met_element, level);

   return grib_key;

}

void
Access::set_grib_key (Grib::Key& grib_key,
                      const Met_Element met_element,
                      const Dtime& base_time,
                      const Integer forecast_hour) const
{

   uint8_t* buffer = grib_key.buffer;

   const Integer yyyy = base_time.get_year ();
   const Integer mm = base_time.get_month ();
   const Integer dd = base_time.get_day ();
   const Integer HH = base_time.get_hour ();
   const Integer MM = base_time.get_minute ();

   buffer[0] = uint8_t (yyyy / 100) + 1;
   buffer[1] = uint8_t (yyyy % 100);
   buffer[2] = uint8_t (mm);
   buffer[3] = uint8_t (dd);
   buffer[4] = uint8_t (HH);
   buffer[5] = uint8_t (MM);
   buffer[6] = 1; // hour
   buffer[7] = uint8_t (forecast_hour);
   buffer[8] = 0;
   buffer[9] = 0;
   buffer[10] = 0;
   buffer[11] = 0;

   if (met_element == RAINFALL_CUMULATIVE)
   {
      buffer[7] = 0;
      buffer[8] = uint8_t (forecast_hour);
      buffer[9] = 4;
      buffer[10] = 0;
      buffer[11] = 1;
   }

}

void
Access::set_grib_key (Grib::Key& grib_key,
                      const Met_Element met_element,
                      const denise::Level& level) const
{

   uint8_t* buffer = grib_key.buffer;
   const uint8_t n = (omega_as_w ? 135 : 134);

   if (level.type == Level::PRESSURE)
   {
      switch (met_element)
      {
         case denise::ZONAL_WIND:          { buffer[12] = 131; break; }
         case denise::MERIDIONAL_WIND:     { buffer[12] = 132; break; }
         case denise::TEMPERATURE:         { buffer[12] = 130; break; }
         case denise::RELATIVE_HUMIDITY:   { buffer[12] = 157; break; }
         case denise::VERTICAL_VELOCITY:   { buffer[12] = n; break; }
         case denise::GEOPOTENTIAL_HEIGHT: { buffer[12] = 156; break; }
      }
      buffer[13] = 100;
      uint16_t hpa = uint16_t (round (level.get_value () * 1e-2)); 
#ifndef WORDS_BIGENDIAN
      swap_endian (&hpa, sizeof (uint16_t));
#endif
      memcpy (buffer + 14, &hpa, sizeof (uint16_t));
   }
   else
   {

      switch (met_element)
      {
         case MEAN_SEA_LEVEL_PRESSURE:
            buffer[12] = 151;
            buffer[13] = 102;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case PRESSURE:
            buffer[12] = 134;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case PPTN:
         case RAINFALL_CUMULATIVE:
            buffer[12] = 61;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case HIGH_CLOUD:
            buffer[12] = 186;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case MIDDLE_CLOUD:
            buffer[12] = 187;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case LOW_CLOUD:
            buffer[12] = 188;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case ZONAL_WIND:
            buffer[12] = 165;
            buffer[13] = 105;
            buffer[14] = 0;
            buffer[15] = 10;
            break;
         case MERIDIONAL_WIND:
            buffer[12] = 166;
            buffer[13] = 105;
            buffer[14] = 0;
            buffer[15] = 10;
            break;
         case TEMPERATURE:
            buffer[12] = 167;
            buffer[13] = 105;
            buffer[14] = 0;
            buffer[15] = 2;
            break;
         case DEW_POINT:
            buffer[12] = 168;
            buffer[13] = 105;
            buffer[14] = 0;
            buffer[15] = 2;
            break;
      }

   }

}

void
Access::initialize_3d_data (const Key& key)
{
   typedef Access::Data_3D Ad_3d;
   Ad_3d* ad_3d_ptr = new Ad_3d (met_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, ad_3d_ptr));
}

void          
Access::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Met_Element>::const_iterator Iterator;
   const Key& key = data_3d.key;

   for (Iterator iterator = met_element_vector.begin ();
        iterator != met_element_vector.end (); iterator++)
   {
      const Met_Element& met_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = get_gvd_3d_ptr (met_element, key);
      data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
   }

   data_3d.set_available ();

}

Geodetic_Vector_Data_2D*
Access::get_initialized_vd_2d (const Integer vector_size) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   Gvd_2d* data_ptr = new Gvd_2d (vector_size, size_2d, domain_2d, false);
   return data_ptr;

}

void 
Access::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Level& level,
                                const Met_Element met_element)
{

   const Integer vi = vector_index;

   switch (met_element)
   {

      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& surface = Level::surface_level ();
         fill_grib_data (gvd_2d, vi, met_element, key, surface);
         return;
      }

   }

   Nwp::fill_ts_diagnosis_data (gvd_2d, vi, key, level, met_element);

}

void
Access::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                        const Integer vector_index,
                        const Key& key,
                        const Met_Element met_element)
{

   switch (met_element)
   {

      case RAINFALL_CUMULATIVE:
      {

         if (key.forecast_hour == 0)
         {
            gvd_2d.initialize (vector_index, 0);
            return;
         }

         const Level& nil_level = Level::nil_level ();
         fill_grib_data (gvd_2d, vector_index,
            RAINFALL_CUMULATIVE, key, nil_level);
         return;

      }

   }

   Nwp::fill_rain_data (gvd_2d, vector_index, key, met_element);

}

void
Access::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                         const Integer vector_index,
                         const Key& key,
                         const Met_Element met_element)
{
   const Level& nil_level = Level::nil_level ();
   fill_grib_data (gvd_2d, vector_index, met_element, key, nil_level);
}

void
Access::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Met_Element met_element)
{

   const Level& screen = Level::screen_level ();

   switch (met_element)
   {

      case TEMPERATURE:
      case DEW_POINT:
      {
         fill_grib_data (gvd_2d, vector_index, met_element, key, screen);
         return;
      }

      case RELATIVE_HUMIDITY:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, DEW_POINT);

         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real t_d = data_ptr->get_datum (1, i, j);
               const Real rh = Moisture::get_rh (t - K, t_d - K);
               gvd_2d.set_datum (vector_index, i, j, rh);
            }
         }

         delete data_ptr;
         return;

      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, met_element);

}

void
Access::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                             const Integer vector_index,
                             const Key& key,
                             const Met_Element met_element)
{

   const Level& ten = Level::ten_metre_level ();

   switch (met_element)
   {

      case ZONAL_WIND:
      case MERIDIONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, met_element, key, ten);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, met_element);

}

void
Access::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Key& key,
                       const Met_Element met_element)
{

   const Level& msl = Level::mean_sea_level ();

   switch (met_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         const Met_Element mslp = MEAN_SEA_LEVEL_PRESSURE;
         fill_grib_data (gvd_2d, vector_index, mslp, key, msl);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, met_element);

}

void
Access::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key,
                                 const Met_Element met_element)
{

   const Level& surface = Level::surface_level ();

   switch (met_element)
   {

      case PRESSURE:
      {
         fill_grib_data (gvd_2d, vector_index, PRESSURE, key, surface);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, met_element);

}

Geodetic_Vector_Data_3D*
Access::get_gvd_3d_ptr (const Met_Element met_element,
                        const Key& key) const
{

   Geodetic_Vector_Data_3D* gvd_3d_ptr =
      new Geodetic_Vector_Data_3D (1, tuple_p, size_2d, domain_2d);
   Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   const bool is_rh = (met_element == RELATIVE_HUMIDITY);
   const bool is_zonal_wind = (met_element == ZONAL_WIND);
   const bool is_meridional_wind = (met_element == MERIDIONAL_WIND);
   const bool is_wind = (is_zonal_wind || is_meridional_wind);

   for (Integer k = 0; k < tuple_p.size (); k++)
   {

      const Real p = tuple_p[k];
      const denise::Level level (Level::PRESSURE, p);
      const Grib::Key& grib_key = get_grib_key (key, met_element, level);

      Access::const_iterator iterator = find (grib_key);
      if (iterator == end ())
      {
         cout << key.base_time << " " << key.forecast_hour << " " <<
                 met_element << " " << level.get_string () <<
                 " throw exception " << endl;
         throw Nwp_Exception ("Access::gvd_3d_ptr Not Available");
      }
      const Grib& grib = *(iterator->second);

      Geodetic_Vector_Data_2D* grib_data_ptr =
         get_grib_data_ptr (grib, grib_key);

      if (is_rh) { grib_data_ptr->scale_offset (0, 0.01, 0); }

      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Integer jj = (j == size_2d.j - 1 ? j - 1 : j);
            Real datum = grib_data_ptr->get_datum (0, i, jj);
            if (is_wind && fabs (datum) > 500) { datum = 0; }
            gvd_3d.set_datum (0, k, i, j, datum);
         }
      }

      delete grib_data_ptr;

   }

   return gvd_3d_ptr;

}

Geodetic_Vector_Data_2D*
Access::get_grib_data_ptr (const Grib& grib,
                           const Grib::Key& grib_key) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = new Gvd_2d (1, size_2d, domain_2d);
   grib.fill_data (*data_ptr, 0, grib_key);

   const Integer last_j = size_2d.j - 1;

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real datum = data_ptr->get_datum (0, i, 0);
      data_ptr->set_datum (0, i, last_j, datum);
   }

   return data_ptr;

}

void
Access::fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                        const Integer vector_index,
                        const Met_Element met_element,
                        const Key& key,
                        const Level& level) const
{

   const Grib::Key& grib_key = get_grib_key (key, met_element, level);
   const Nwp_Exception exception ("Access::fill_grib_data Not Available");

   Access::const_iterator iterator = find (grib_key);
   if (iterator == end ()) { throw exception; }
   const Grib& grib = *(iterator->second);

   Geodetic_Vector_Data_2D* grib_data_ptr = get_grib_data_ptr (grib, grib_key);

   const bool is_zonal_wind = (met_element == ZONAL_WIND);
   const bool is_meridional_wind = (met_element == MERIDIONAL_WIND);
   const bool is_wind = (is_zonal_wind || is_meridional_wind);

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Integer jj = (j == size_2d.j - 1 ? j - 1 : j);
         Real datum = grib_data_ptr->get_datum (0, i, jj);
         if (is_wind && fabs (datum) > 500) { datum = 0; }
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete grib_data_ptr;

}

Access::Access (const Dstring& description,
                const Dstring& data_path,
                const Dstring& search_string,
                const bool omega_as_w)
   : Nwp (description, data_path),
     data_path (data_path),
     search_string (search_string),
     omega_as_w (omega_as_w)
{

   this->status = "Unloaded";

   met_element_vector.push_back (denise::TEMPERATURE);
   met_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   met_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   met_element_vector.push_back (denise::ZONAL_WIND);
   met_element_vector.push_back (denise::MERIDIONAL_WIND);
   met_element_vector.push_back (denise::VERTICAL_VELOCITY);

}

Access::~Access ()
{
   clean_up ();
}

void
Access::survey ()
{

   typedef map<Grib::Key, Grib::Header*> Header_Ptr_Map;

   const vector<Dstring>& dir_listing = get_dir_listing (path, search_string);
   for (vector<Dstring>::const_iterator iterator = dir_listing.begin ();
        iterator != dir_listing.end (); iterator++)
   {

      const Dstring& file_name = *(iterator);
      const Dstring& file_path = path + "/" + file_name;

cout << "Grib file_path " << file_path << endl;
      Grib* grib_ptr = new Grib (file_path);
      const Header_Ptr_Map& header_ptr_map = grib_ptr->get_header_ptr_map ();

      for (Header_Ptr_Map::const_iterator i = header_ptr_map.begin ();
           i != header_ptr_map.end (); i++)
      {
         const Grib::Header& header = *(i->second);
         const Grib::Pds& pds = header.get_pds ();
         const Dtime base_time = pds.get_base_time ();
         const Integer forecast_hour = pds.get_forecast_time ().get_p1 ();

         const Grib::Key grib_key (pds);
         const Nwp::Key key (base_time, forecast_hour);

         key_multimap.add (key);
         insert (make_pair (grib_key, grib_ptr));

      }

      grib_ptr_set.insert (grib_ptr);

   }

   if (size () > 0)
   {

      set<uint16_t> set_p;
      const Grib& grib = *(begin ()->second);
      const Header_Ptr_Map& header_ptr_map = grib.get_header_ptr_map ();

      const Grib::Header& first_header = *(header_ptr_map.begin ()->second);
      const Grib::Gds& gds = first_header.get_gds ();
      size_2d = gds.get_size_2d ();
      const Real latitude_0 = gds.get_int (10, 3) * 1e-3;
      const Real longitude_0 = gds.get_int (13, 3) * 1e-3;
      const Real latitude_1 = gds.get_int (17, 3) * 1e-3;
      const Real longitude_1 = gds.get_int (20, 3) * 1e-3;

      Domain_1D& domain_latitude = domain_2d.domain_x;
      Domain_1D& domain_longitude = domain_2d.domain_y;
      domain_latitude.start = std::min (latitude_0, latitude_1);
      domain_latitude.end = std::max (latitude_0, latitude_1);
      domain_longitude.start = std::min (longitude_0, longitude_1);
      domain_longitude.end = std::max (longitude_0, longitude_1);

      for (Header_Ptr_Map::const_iterator iterator = header_ptr_map.begin ();
           iterator != header_ptr_map.end (); iterator++)
      {

         const Grib::Header& header = *(iterator->second);
         const Grib::Pds& pds = header.get_pds ();
         const Grib::Pds::Level& level = pds.get_level ();

         if (level.get_uint (0, 1) == 100)
         {
            const uint16_t p = level.get_uint (1, 2);
            if (p < 100) { continue; }
            set_p.insert (p);
         }

      }

      for (set<uint16_t>::const_iterator iterator = set_p.begin ();
           iterator != set_p.end (); iterator++)
      {
         const Real p = Real (*(iterator)) * 100;
         tuple_p.push_back (p);
      }

   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

   for (Key_Multimap::const_iterator iterator = key_multimap.begin ();
        iterator != key_multimap.end (); iterator++)
   {
      const Key key (iterator->first, iterator->second);
      initialize_3d_data (key);
   }

}

void
Access::clean_up ()
{

   for (set<Grib*>::iterator iterator = grib_ptr_set.begin ();
        iterator != grib_ptr_set.end (); iterator++)
   {
      Grib* grib_ptr = *(iterator);
      delete grib_ptr;
   }

   key_multimap.clear ();
   tuple_p.clear ();
   clear_data_3d_ptr_map ();

   clear ();
   grib_ptr_set.clear ();

}

Ecmwf::Data_3D::Data_3D (const vector<Met_Element>& met_element_vector,
                         const Key& key)
   : Nwp::Data_3D (met_element_vector, key)
{
}

Real
Ecmwf::Data_3D::evaluate (const Met_Element element,
                          const Real p,
                          const Real latitude,
                          const Real longitude,
                          const Evaluate_Op evaluate_op) const
{

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& R = denise::MIXING_RATIO;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real q = Nwp::Data_3D::evaluate (R, p, latitude, longitude);
         const Real rh = 0.00263 * p*q / (exp ((17.67 * (t-K) / (t-29.65))));
         return Moisture::get_t_d (t, rh);
      }

      case denise::RELATIVE_HUMIDITY:
      {
         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& R = denise::MIXING_RATIO;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real q = Nwp::Data_3D::evaluate (R, p, latitude, longitude);
         return 0.00263 * p * q / (exp ((17.67 * (t - K) / (t - 29.65))));
      }

      case denise::OMEGA:
      {
         const Met_Element& T = TEMPERATURE;
         const Met_Element& W = VERTICAL_VELOCITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real w = Nwp::Data_3D::evaluate (W, p, latitude, longitude);
         const Real rho = p / (R_d * t);
         return -rho * g * w;
      }

   }

   return Nwp::Data_3D::evaluate (element, p,
      latitude, longitude, evaluate_op);

}

Grib::Key
Ecmwf::get_grib_key (const Key& key,
                     const Met_Element met_element,
                     const Level& level) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   Grib::Key grib_key;
   set_grib_key (grib_key, met_element, base_time, forecast_hour);
   set_grib_key (grib_key, met_element, level);

   return grib_key;

}

void
Ecmwf::set_grib_key (Grib::Key& grib_key,
                     const Met_Element met_element,
                     const Dtime& base_time,
                     const Integer forecast_hour) const
{

   uint8_t* buffer = grib_key.buffer;

   const Integer yyyy = base_time.get_year ();
   const Integer mm = base_time.get_month ();
   const Integer dd = base_time.get_day ();
   const Integer HH = base_time.get_hour ();
   const Integer MM = base_time.get_minute ();

   buffer[0] = uint8_t (yyyy / 100) + 1;
   buffer[1] = uint8_t (yyyy % 100);
   buffer[2] = uint8_t (mm);
   buffer[3] = uint8_t (dd);
   buffer[4] = uint8_t (HH);
   buffer[5] = uint8_t (MM);
   buffer[6] = 1; // hour
   buffer[7] = uint8_t (forecast_hour);
   buffer[8] = 0;
   buffer[9] = 0;
   buffer[10] = 0;
   buffer[11] = 0;

   if (met_element == RAINFALL_CUMULATIVE)
   {
      buffer[7] = 0;
      buffer[8] = uint8_t (forecast_hour);
      buffer[9] = 4;
      buffer[10] = 0;
      buffer[11] = 1;
   }

//cout << "Ecmwf::set_grib_key " << grib_key << endl;

}

void
Ecmwf::set_grib_key (Grib::Key& grib_key,
                     const Met_Element met_element,
                     const denise::Level& level) const
{

   uint8_t* buffer = grib_key.buffer;
   //const uint8_t n = (omega_as_w ? 135 : 134);

   if (level.type == Level::PRESSURE)
   {
      switch (met_element)
      {
         case denise::ZONAL_WIND:          { buffer[12] = 131; break; }
         case denise::MERIDIONAL_WIND:     { buffer[12] = 132; break; }
         case denise::TEMPERATURE:         { buffer[12] = 130; break; }
         case denise::MIXING_RATIO:        { buffer[12] = 133; break; }
         case denise::VERTICAL_VELOCITY:   { buffer[12] = 135; break; }
         case denise::GEOPOTENTIAL_HEIGHT: { buffer[12] = 156; break; }
      }
      buffer[13] = 100;
      uint16_t hpa = uint16_t (round (level.get_value () * 1e-2)); 
#ifndef WORDS_BIGENDIAN
      swap_endian (&hpa, sizeof (uint16_t));
#endif
      memcpy (buffer + 14, &hpa, sizeof (uint16_t));
   }
   else
   {

      switch (met_element)
      {
         case MEAN_SEA_LEVEL_PRESSURE:
            buffer[12] = 151;
            buffer[13] = 102;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case PRESSURE:
            buffer[12] = 134;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case PPTN:
         case RAINFALL_CUMULATIVE:
            buffer[12] = 61;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case HIGH_CLOUD:
            buffer[12] = 186;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case MIDDLE_CLOUD:
            buffer[12] = 187;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case LOW_CLOUD:
            buffer[12] = 188;
            buffer[13] = 200;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case ZONAL_WIND:
            buffer[12] = 165;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case MERIDIONAL_WIND:
            buffer[12] = 166;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case TEMPERATURE:
            buffer[12] = 167;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
         case DEW_POINT:
            buffer[12] = 168;
            buffer[13] = 1;
            buffer[14] = 0;
            buffer[15] = 0;
            break;
      }

//cout << "!!Ecmwf::set_grib_key " << grib_key << endl;
   }

}

void
Ecmwf::initialize_3d_data (const Key& key)
{
   typedef Ecmwf::Data_3D Ed_3d;
   Ed_3d* ed_3d_ptr = new Ed_3d (met_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, ed_3d_ptr));
}

void          
Ecmwf::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Met_Element>::const_iterator Iterator;
   const Key& key = data_3d.key;

   for (Iterator iterator = met_element_vector.begin ();
        iterator != met_element_vector.end (); iterator++)
   {
      const Met_Element& met_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = get_gvd_3d_ptr (met_element, key);
      data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
   }

   data_3d.set_available ();

}

Geodetic_Vector_Data_2D*
Ecmwf::get_initialized_vd_2d (const Integer vector_size) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   const Size_2D& size_2d = size_2d_map.at (1000);
   Gvd_2d* data_ptr = new Gvd_2d (vector_size, size_2d, domain_2d, false);
   return data_ptr;

}

void 
Ecmwf::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Level& level,
                               const Met_Element met_element)
{

   const Integer vi = vector_index;

   switch (met_element)
   {

      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& surface = Level::surface_level ();
         fill_grib_data (gvd_2d, vi, met_element, key, surface);
         return;
      }

   }

   Nwp::fill_ts_diagnosis_data (gvd_2d, vi, key, level, met_element);

}

void
Ecmwf::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Key& key,
                       const Met_Element met_element)
{

   switch (met_element)
   {

      case RAINFALL_CUMULATIVE:
      {

         if (key.forecast_hour == 0)
         {
            gvd_2d.initialize (vector_index, 0);
            return;
         }

         const Level& nil_level = Level::nil_level ();
         fill_grib_data (gvd_2d, vector_index,
            RAINFALL_CUMULATIVE, key, nil_level);
         return;

      }

   }

   Nwp::fill_rain_data (gvd_2d, vector_index, key, met_element);

}

void
Ecmwf::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                        const Integer vector_index,
                        const Key& key,
                        const Met_Element met_element)
{
   const Level& nil_level = Level::nil_level ();
   fill_grib_data (gvd_2d, vector_index, met_element, key, nil_level);
}

void
Ecmwf::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Met_Element met_element)
{

   const Level& screen = Level::screen_level ();

   switch (met_element)
   {

      case TEMPERATURE:
      case DEW_POINT:
      {
         fill_grib_data (gvd_2d, vector_index, met_element, key, screen);
         return;
      }

      case RELATIVE_HUMIDITY:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, DEW_POINT);

         const Size_2D& size_2d = gvd_2d.get_size_2d ();
         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            const Real latitude = gvd_2d.get_coordinate (0, i);
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real longitude = gvd_2d.get_coordinate (1, j);
               const Real t = data_ptr->evaluate (0, latitude, longitude);
               const Real t_d = data_ptr->evaluate (1, latitude, longitude);
               const Real rh = Moisture::get_rh (t - K, t_d - K);
               gvd_2d.set_datum (vector_index, i, j, rh);
            }
         }

         delete data_ptr;
         return;

      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, met_element);

}

void
Ecmwf::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                            const Integer vector_index,
                            const Key& key,
                            const Met_Element met_element)
{

   const Level& ten = Level::ten_metre_level ();

   switch (met_element)
   {

      case ZONAL_WIND:
      case MERIDIONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, met_element, key, ten);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, met_element);

}

void
Ecmwf::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Met_Element met_element)
{

   const Level& msl = Level::mean_sea_level ();

   switch (met_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         const Met_Element mslp = MEAN_SEA_LEVEL_PRESSURE;
         fill_grib_data (gvd_2d, vector_index, mslp, key, msl);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, met_element);

}

void
Ecmwf::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                                const Integer vector_index,
                                const Key& key,
                                const Met_Element met_element)
{

   const Level& surface = Level::surface_level ();

   switch (met_element)
   {

      case PRESSURE:
      {
         fill_grib_data (gvd_2d, vector_index, PRESSURE, key, surface);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, met_element);

}

Geodetic_Vector_Data_3D*
Ecmwf::get_gvd_3d_ptr (const Met_Element met_element,
                       const Key& key) const
{

   const Size_2D& size_2d = size_2d_map.at (1000);
cout << "GVD 3D " << size_2d << " " << tuple_p << " " << met_element << endl;
   Geodetic_Vector_Data_3D* gvd_3d_ptr =
      new Geodetic_Vector_Data_3D (1, tuple_p, size_2d, domain_2d);
   Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   const bool is_rh = (met_element == RELATIVE_HUMIDITY);
   const bool is_zonal_wind = (met_element == ZONAL_WIND);
   const bool is_meridional_wind = (met_element == MERIDIONAL_WIND);
   const bool is_wind = (is_zonal_wind || is_meridional_wind);

   for (Integer k = 0; k < tuple_p.size (); k++)
   {

      const Real p = tuple_p[k];
      const denise::Level level (Level::PRESSURE, p);
      const Grib::Key& grib_key = get_grib_key (key, met_element, level);

      Ecmwf::const_iterator iterator = find (grib_key);
      if (iterator == end ())
      {
         //cout << key.base_time << " " << key.forecast_hour << " " <<
         //        met_element << " " << level.get_string () <<
         //        " throw exception " << endl;
         throw Nwp_Exception ("Ecmwf::gvd_3d_ptr Not Available");
      }
      const Grib& grib = *(iterator->second);

      Geodetic_Vector_Data_2D* grib_data_ptr =
         get_grib_data_ptr (grib, grib_key);

      if (is_rh) { grib_data_ptr->scale_offset (0, 0.01, 0); }

      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         const Real latitude = gvd_3d_ptr->get_coordinate (1, i);
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real longitude = gvd_3d_ptr->get_coordinate (2, j);
            const Integer jj = (j == size_2d.j - 1 ? j - 1 : j);
            Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
            if (is_wind && fabs (datum) > 500) { datum = 0; }
            gvd_3d.set_datum (0, k, i, j, datum);
         }
      }

      delete grib_data_ptr;

   }

   return gvd_3d_ptr;

}

Geodetic_Vector_Data_2D*
Ecmwf::get_grib_data_ptr (const Grib& grib,
                          const Grib::Key& grib_key) const
{

   //const Size_2D& size_2d = size_2d_map.at (1000);

   const Grib::Gds& gds = grib.get_header_ptr_map ().at (grib_key)->get_gds ();
   const Size_2D& size_2d = gds.get_size_2d ();

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = new Gvd_2d (1, size_2d, domain_2d);
   grib.fill_data (*data_ptr, 0, grib_key);

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   const Integer last_j = size_2d.j - 1;

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real datum = data_ptr->get_datum (0, i, 0);
      data_ptr->set_datum (0, i, last_j, datum);
   }

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   return data_ptr;

}

void
Ecmwf::fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Met_Element met_element,
                       const Key& key,
                       const Level& level) const
{

   const Grib::Key& grib_key = get_grib_key (key, met_element, level);
   const Nwp_Exception exception ("Ecmwf::fill_grib_data Not Available");

   Ecmwf::const_iterator iterator = find (grib_key);
   if (iterator == end ()) { throw exception; }
   const Grib& grib = *(iterator->second);

   Geodetic_Vector_Data_2D* grib_data_ptr = get_grib_data_ptr (grib, grib_key);

   const bool is_zonal_wind = (met_element == ZONAL_WIND);
   const bool is_meridional_wind = (met_element == MERIDIONAL_WIND);
   const bool is_wind = (is_zonal_wind || is_meridional_wind);

   const Size_2D& size_2d = gvd_2d.get_size_2d ();

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_coordinate (0, i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_coordinate (1, j);
         Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
         if (is_wind && fabs (datum) > 500) { datum = 0; }
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete grib_data_ptr;

}

Ecmwf::Ecmwf (const Dstring& description,
              const Dstring& data_path,
              const Dstring& search_string,
              const bool omega_as_w)
   : Nwp (description, data_path),
     data_path (data_path),
     search_string (search_string),
     omega_as_w (omega_as_w)
{

   this->status = "Unloaded";

   met_element_vector.push_back (denise::TEMPERATURE);
   met_element_vector.push_back (denise::MIXING_RATIO);
   met_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   met_element_vector.push_back (denise::ZONAL_WIND);
   met_element_vector.push_back (denise::MERIDIONAL_WIND);
   met_element_vector.push_back (denise::VERTICAL_VELOCITY);

}

Ecmwf::~Ecmwf ()
{
   clean_up ();
}

void
Ecmwf::survey ()
{

   typedef map<Grib::Key, Grib::Header*> Header_Ptr_Map;

   const vector<Dstring>& dir_listing = get_dir_listing (path, search_string);
   for (vector<Dstring>::const_iterator iterator = dir_listing.begin ();
        iterator != dir_listing.end (); iterator++)
   {

      const Dstring& file_name = *(iterator);
      const Dstring& file_path = path + "/" + file_name;

      Grib* grib_ptr = new Grib (file_path);
      const Header_Ptr_Map& header_ptr_map = grib_ptr->get_header_ptr_map ();

      for (Header_Ptr_Map::const_iterator i = header_ptr_map.begin ();
           i != header_ptr_map.end (); i++)
      {
         const Grib::Header& header = *(i->second);
         const Grib::Pds& pds = header.get_pds ();
         const Dtime base_time = pds.get_base_time ();
         const Integer forecast_hour = pds.get_forecast_time ().get_p1 ();

         const Grib::Key grib_key (pds);
         const Nwp::Key key (base_time, forecast_hour);

         key_multimap.add (key);
         insert (make_pair (grib_key, grib_ptr));

      }

      grib_ptr_set.insert (grib_ptr);

   }

   if (size () > 0)
   {

      set<uint16_t> set_p;
      const Grib& grib = *(begin ()->second);
      const Header_Ptr_Map& header_ptr_map = grib.get_header_ptr_map ();

      const Grib::Header& first_header = *(header_ptr_map.begin ()->second);
      const Grib::Gds& gds = first_header.get_gds ();
      const Real latitude_0 = gds.get_int (10, 3) * 1e-3;
      const Real longitude_0 = gds.get_int (13, 3) * 1e-3;
      const Real latitude_1 = gds.get_int (17, 3) * 1e-3;
      const Real longitude_1 = gds.get_int (20, 3) * 1e-3;

      Domain_1D& domain_latitude = domain_2d.domain_x;
      Domain_1D& domain_longitude = domain_2d.domain_y;
      domain_latitude.start = std::min (latitude_0, latitude_1);
      domain_latitude.end = std::max (latitude_0, latitude_1);
      domain_longitude.start = std::min (longitude_0, longitude_1);
      domain_longitude.end = std::max (longitude_0, longitude_1);

      for (Header_Ptr_Map::const_iterator iterator = header_ptr_map.begin ();
           iterator != header_ptr_map.end (); iterator++)
      {

         const Grib::Header& header = *(iterator->second);
         const Grib::Pds& pds = header.get_pds ();
         const Grib::Pds::Level& level = pds.get_level ();

         const Grib::Gds& gds = header.get_gds ();

         if (level.get_uint (0, 1) == 100)
         {
            const uint16_t p = level.get_uint (1, 2);
            if (p < 100) { continue; }
            size_2d_map[p] = gds.get_size_2d ();
            set_p.insert (p);
         }

      }

      for (set<uint16_t>::const_iterator iterator = set_p.begin ();
           iterator != set_p.end (); iterator++)
      {
         const Real p = Real (*(iterator)) * 100;
         tuple_p.push_back (p);
      }

   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

   for (Key_Multimap::const_iterator iterator = key_multimap.begin ();
        iterator != key_multimap.end (); iterator++)
   {
      const Key key (iterator->first, iterator->second);
      initialize_3d_data (key);
   }

}

void
Ecmwf::clean_up ()
{

   for (set<Grib*>::iterator iterator = grib_ptr_set.begin ();
        iterator != grib_ptr_set.end (); iterator++)
   {
      Grib* grib_ptr = *(iterator);
      delete grib_ptr;
   }

   key_multimap.clear ();
   tuple_p.clear ();
   clear_data_3d_ptr_map ();

   clear ();
   grib_ptr_set.clear ();

}

Gfs3::Data_3D::Data_3D (const vector<Met_Element>& met_element_vector,
                        const Key& key)
   : Nwp::Data_3D (met_element_vector, key)
{
}

Real
Gfs3::Data_3D::evaluate (const Met_Element element,
                         const Real p,
                         const Real latitude,
                         const Real longitude,
                         const Evaluate_Op evaluate_op) const
{

   typedef Nwp::Data_3D Nd_3d;
   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& RH = denise::RELATIVE_HUMIDITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real rh = Nwp::Data_3D::evaluate (RH, p, latitude, longitude);
         return Moisture::get_t_d (t, rh);
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& O = denise::OMEGA;

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Evaluate_Op& eo = evaluate_op;
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real ts = Nd_3d::evaluate (T, p, latitude, longitude, eo);
            const Real os = Nd_3d::evaluate (O, p, latitude, longitude, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nd_3d::evaluate (element, p, latitude, longitude, eo);

}

Grib::Key
Gfs3::get_grib_key (const Key& key,
                    const Met_Element met_element,
                    const Level& level) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   Grib::Key grib_key;
   set_grib_key (grib_key, base_time);
   set_grib_key (grib_key, met_element, forecast_hour);
   set_grib_key (grib_key, met_element);
   set_grib_key (grib_key, met_element, level);

   return grib_key;

}

void
Gfs3::set_grib_key (Grib::Key& grib_key,
                    const Dtime& base_time) const
{

   uint8_t* buffer = grib_key.buffer;

   const Integer yyyy = base_time.get_year ();
   const Integer mm = base_time.get_month ();
   const Integer dd = base_time.get_day ();
   const Integer HH = base_time.get_hour ();
   const Integer MM = base_time.get_minute ();

   buffer[0] = uint8_t (yyyy / 100) + 1;
   buffer[1] = uint8_t (yyyy % 100);
   buffer[2] = uint8_t (mm);
   buffer[3] = uint8_t (dd);
   buffer[4] = uint8_t (HH);
   buffer[5] = uint8_t (MM);

}

void
Gfs3::set_grib_key (Grib::Key& grib_key,
                    const Met_Element met_element,
                    const Integer forecast_hour) const
{

   uint8_t* buffer = grib_key.buffer;

   buffer[6] = 1;
   buffer[10] = 0;
   buffer[11] = 0;

   switch (met_element)
   {

      default:
         buffer[7] = uint8_t (forecast_hour >> 8);
         buffer[8] = uint8_t (forecast_hour % 256);
         buffer[9] = 10;
         buffer[10] = 0;
         buffer[11] = 0;
         break;

      case PPT3:
         buffer[7] = uint8_t (forecast_hour - 3);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;

      case PPT6:
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;

      case PPTN:
         buffer[7] = 0;
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;

      case TOTAL_CLOUD:
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;

      case HIGH_CLOUD:
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;

      case MIDDLE_CLOUD:
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;

      case LOW_CLOUD:
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;

   }

}

void
Gfs3::set_grib_key (Grib::Key& grib_key,
                    const Met_Element met_element) const
{

   uint8_t* buffer = grib_key.buffer;
   uint8_t& octet = buffer[12];

   switch (met_element)
   {
      case PRESSURE:
         octet = 1;
         break;
      case TEMPERATURE:
         octet = 11;
         break;
      case TEMPERATURE_CELCIUS:
         octet = 0;
         break;
      case MIX_DOWN_TEMPERATURE:
         octet = 0;
         break;
      case ZONAL_WIND:
         octet = 33;
         break;
      case MERIDIONAL_WIND:
         octet = 34;
         break;
      case WIND_SPEED:
         octet = 32;
         break;
      case VERTICAL_VELOCITY:
         octet = 40;
         break;
      case GEOPOTENTIAL_HEIGHT:
         octet = 7;
         break;
      case DEW_POINT:
         octet = 17;
         break;
      case MEAN_SEA_LEVEL_PRESSURE:
         octet = 2;
         break;
      case HIGH_CLOUD:
         octet = 71;
         break;
      case MIDDLE_CLOUD:
         octet = 71;
         break;
      case LOW_CLOUD:
         octet = 71;
         break;
      case TOTAL_CLOUD:
         octet = 71;
         break;
      case OMEGA:
         octet = 39;
         break;
      case PPT3:
         octet = 61;
         break;
      case PPT6:
         octet = 61;
         break;
      case PPTN:
         octet = 61;
         break;
      case RAINFALL_STEP:
         octet = 0;
         break;
      case RELATIVE_HUMIDITY:
         octet = 52;
         break;
      case DEW_POINT_DEPRESSION:
         octet = 0;
         break;
      case POTENTIAL_TEMPERATURE:
         octet = 13;
         break;
      case THETA_E:
         octet = 14;
         break;
      case MIXING_RATIO:
         octet = 53;
         break;
      case MONTGOMERY:
         octet = 37;
         break;
      case SLI:
         octet = 131;
         break;
      case SHOWALTER:
         octet = 0;
         break;
      case LI_700:
         octet = 0;
         break;
      case LI_THUNDER:
         octet = 0;
         break;
      case K_INDEX:
         octet = 133;
         break;
      case TOTAL_TOTALS:
         octet = 0;
         break;
      case CAPE:
         octet = 157;
         break;
      case PRECIPITABLE_WATER:
         octet = 54;
         break;
      case FOG_FRACTION:
         octet = 0;
         break;
      case THICKNESS:
         octet = 0;
         break;
      case POTENTIAL_VORTICITY:
         octet = 149;
         break;
      case ABSOLUTE_VORTICITY:
         octet = 41;
         break;
      case SHEAR_VORTICITY:
         octet = 0;
         break;
      case CURVATURE_VORTICITY:
         octet = 0;
         break;
   }

}

void
Gfs3::set_grib_key (Grib::Key& grib_key,
                    const Met_Element met_element,
                    const denise::Level& level) const
{

   uint8_t* buffer = grib_key.buffer;

   switch (met_element)
   {
      case TOTAL_CLOUD:
         buffer[13] = 200;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
      case HIGH_CLOUD:
         buffer[13] = 214;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
      case MIDDLE_CLOUD:
         buffer[13] = 224;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
      case LOW_CLOUD:
         buffer[13] = 234;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
      case PRECIPITABLE_WATER:
         buffer[13] = 200;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
   }

   switch (level.type)
   {
      case Level::PRESSURE:
      {
         const uint16_t p = uint16_t (round (level.value * 1e-2));
         buffer[13] = 100;
         buffer[14] = uint8_t (p >> 8);
         buffer[15] = uint8_t (p % 256);
         break;
      }
      case Level::THETA:
      {
         const uint16_t theta = uint16_t (round (level.value));
         buffer[13] = 113;
         buffer[14] = uint8_t (theta >> 8);
         buffer[15] = uint8_t (theta % 256);
         break;
      }
      case Level::SIGMA:
      {
         const uint16_t sigma = uint16_t (round (level.value * 10000));
         buffer[13] = 107;
         buffer[14] = uint8_t (sigma >> 8);
         buffer[15] = uint8_t (sigma % 256);
         break;
      }
      case Level::SCREEN:
      {
         const uint16_t z = uint16_t (2);
         buffer[13] = 105;
         buffer[14] = uint8_t (z >> 8);
         buffer[15] = uint8_t (z % 256);
         break;
      }
      case Level::FIFTY_METRE:
      {
         const uint16_t z = uint16_t (50);
         buffer[13] = 105;
         buffer[14] = uint8_t (z >> 8);
         buffer[15] = uint8_t (z % 256);
         break;
      }
      case Level::TEN_METRE:
      {
         const uint16_t z = uint16_t (10);
         buffer[13] = 105;
         buffer[14] = uint8_t (z >> 8);
         buffer[15] = uint8_t (z % 256);
         break;
      }
      case Level::MEAN_SEA:
      {
         buffer[13] = 102;
         buffer[14] = 0;
         buffer[15] = 0;
         break;
      }
      case Level::SURFACE:
      {
         buffer[13] = 1;
         buffer[14] = 0;
         buffer[15] = 0;
         break;
      }
      case Level::NIL:
      case Level::NAL:
      {
         buffer[13] = 0;
         buffer[14] = 0;
         buffer[15] = 0;
         break;
      }
   }

}

void
Gfs3::initialize_3d_data (const Key& key)
{
   typedef Gfs3::Data_3D G3d_3d;
   G3d_3d* g3d_3d_ptr = new G3d_3d (met_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, g3d_3d_ptr));
cout << "data_3d_ptr_map.size () = " << &data_3d_ptr_map << ": " << data_3d_ptr_map.size () << "  key = " << key.base_time.get_string () << " " << key.forecast_hour << endl;
}

void          
Gfs3::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Met_Element>::const_iterator Iterator;
   const Key& key = data_3d.key;

   for (Iterator iterator = met_element_vector.begin ();
        iterator != met_element_vector.end (); iterator++)
   {
      const Met_Element& met_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = get_gvd_3d_ptr (met_element, key);
      data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
   }

   data_3d.set_available ();

}

Geodetic_Vector_Data_2D*
Gfs3::get_initialized_vd_2d (const Integer vector_size) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   Gvd_2d* data_ptr = new Gvd_2d (vector_size, size_2d, domain_2d, false);
   return data_ptr;

}

void 
Gfs3::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Level& level,
                              const Met_Element met_element)
{

   const Integer vi = vector_index;

   switch (met_element)
   {

      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& surface = Level::surface_level ();
         fill_grib_data (gvd_2d, vi, met_element, key, surface);
         return;
      }

   }

   Nwp::fill_ts_diagnosis_data (gvd_2d, vi, key, level, met_element);

}

void
Gfs3::fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key)
{

   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour < 0)
   {
      throw Nwp_Exception ("Forecast Hour < 0");
      return;
   }

   gvd_2d.initialize (vector_index, 0);
   if (forecast_hour == 0) { return; }

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();

   Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
   Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

   for (Integer fh = 0; fh < forecast_hour; fh += 6)
   {
      fill_grib_data (precip_data, 0, PPT6, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (precip_data, 0, PPT3, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   delete precip_data_ptr;

}

void
Gfs3::fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();
   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour == 0)
   {
      gvd_2d.initialize (vector_index, 0);
      return;
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (gvd_2d, vector_index, PPT3, key, surface_level);
   }
   else
   {

      fill_grib_data (gvd_2d, vector_index, PPT6, key, surface_level);

      Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
      Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

      try
      {
         const Dtime& base_time = key.base_time;
         const Integer forecast_hour = key.forecast_hour;
         const Key prev_key (base_time, forecast_hour - 3);
         fill_grib_data (precip_data, 0, PPT3, prev_key, surface_level);
         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real ppt3 = precip_data.get_datum (0, i, j);
               gvd_2d.get_datum (vector_index, i, j) -= ppt3;
            }
         }
      }
      catch (const Nwp_Exception& ne)
      {
      }

      delete precip_data_ptr;

   }

}

void
Gfs3::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Met_Element met_element)
{

   switch (met_element)
   {

      case RAINFALL_CUMULATIVE:
      {
         fill_cumulative_rain_data (gvd_2d, vector_index, key);
         return;
      }

      case RAINFALL_STEP:
      {
         fill_step_rain_data (gvd_2d, vector_index, key);
         return;
      }

   }

   Nwp::fill_rain_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs3::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Key& key,
                       const Met_Element met_element)
{
   const Level& nil_level = Level::nil_level ();
   fill_grib_data (gvd_2d, vector_index, met_element, key, nil_level);
}

void
Gfs3::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Met_Element met_element)
{

   const Level& screen = Level::screen_level ();
   const Level& surface = Level::surface_level ();

   switch (met_element)
   {

      case TEMPERATURE:
      {
         fill_grib_data (gvd_2d, vector_index, TEMPERATURE, key, screen);
         return;
      }

      case DEW_POINT:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, RELATIVE_HUMIDITY);

         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real rh = data_ptr->get_datum (1, i, j);
               const Real t_d = Moisture::get_t_d (t, rh);
               gvd_2d.set_datum (vector_index, i, j, t_d);
            }
         }

         delete data_ptr;
         return;

      }

      case RELATIVE_HUMIDITY:
      {
         fill_grib_data (gvd_2d, vector_index, RELATIVE_HUMIDITY, key, screen);
         return;
      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs3::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key,
                           const Met_Element met_element)
{

   const Level& ten = Level::ten_metre_level ();

   switch (met_element)
   {

      case ZONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, ZONAL_WIND, key, ten);
         return;
      }

      case MERIDIONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, MERIDIONAL_WIND, key, ten);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs3::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                     const Integer vector_index,
                     const Key& key,
                     const Met_Element met_element)
{

   const Level& msl = Level::mean_sea_level ();

   switch (met_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         const Met_Element mslp = MEAN_SEA_LEVEL_PRESSURE;
         fill_grib_data (gvd_2d, vector_index, mslp, key, msl);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs3::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Met_Element met_element)
{

   const Level& surface = Level::surface_level ();

   switch (met_element)
   {

      case PRESSURE:
      {
         fill_grib_data (gvd_2d, vector_index, PRESSURE, key, surface);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, met_element);

}

/*
void
Gfs3::fill_grib_data (Geodetic_Vector_Data_3D& gvd_3d,
                      const Integer vector_index,
                      const Met_Element met_element,
                      const Key& key) const
{

   Gfs3::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib& grib = *(iterator->second);

   for (Integer k = 0; k < tuple_p.size (); k++)
   {
      const Real p = tuple_p[k];
      const Level level (Level::PRESSURE, p);
      const Grib::Key& key = get_grib_key (key, met_element, level);
      grib.fill_data (gvd_3d, vector_index, k, key);
   }

}
*/

Geodetic_Vector_Data_3D*
Gfs3::get_gvd_3d_ptr (const Met_Element met_element,
                      const Key& key) const
{

   Geodetic_Vector_Data_3D* gvd_3d_ptr =
      new Geodetic_Vector_Data_3D (1, tuple_p, size_2d, domain_2d);
   Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   for (Integer k = 0; k < tuple_p.size (); k++)
   {

      const Real p = tuple_p[k];
      const denise::Level level (Level::PRESSURE, p);
      const Grib::Key& grib_key = get_grib_key (key, met_element, level);

      Gfs3::const_iterator iterator = find (key);
      if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
      const Grib& grib = *(iterator->second);

      Geodetic_Vector_Data_2D* grib_data_ptr =
         get_grib_data_ptr (grib, grib_key);

      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         const Real latitude = grib_data_ptr->get_latitude (i);
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real longitude = grib_data_ptr->get_longitude (j);
            const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
            gvd_3d.set_datum (0, k, i, j, datum);
         }
      }

      delete grib_data_ptr;

   }

   return gvd_3d_ptr;

}

Geodetic_Vector_Data_2D*
Gfs3::get_grib_data_ptr (const Grib& grib,
                         const Grib::Key& grib_key) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = new Gvd_2d (1, size_2d, domain_2d, true);
   grib.fill_data (*data_ptr, 0, grib_key);

   const Integer last_j = size_2d.j - 1;

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real datum = data_ptr->get_datum (0, i, 0);
      data_ptr->set_datum (0, i, last_j, datum);
   }

   return data_ptr;

}

void
Gfs3::fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Met_Element met_element,
                      const Key& key,
                      const Level& level) const
{

   Gfs3::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib& grib = *(iterator->second);

   const Grib::Key& grib_key = get_grib_key (key, met_element, level);
   Geodetic_Vector_Data_2D* grib_data_ptr = get_grib_data_ptr (grib, grib_key);

   Geodetic_Vector_Data_2D grib_data (1, size_2d, domain_2d, true);

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_latitude (i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_longitude (j);
         const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete grib_data_ptr;

}

Gfs3::Gfs3 (const Dstring& data_path)
   : Nwp ("Gfs3", data_path),
     data_path (data_path),
     size_2d (181, 361),
     domain_2d (Domain_1D (-90, 90), Domain_1D (0, 360))
{

   this->status = "Unloaded";

   met_element_vector.push_back (denise::TEMPERATURE);
   met_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   met_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   met_element_vector.push_back (denise::ZONAL_WIND);
   met_element_vector.push_back (denise::MERIDIONAL_WIND);
   met_element_vector.push_back (denise::OMEGA);

}

Gfs3::~Gfs3 ()
{
   clean_up ();
}

void
Gfs3::survey ()
{

   const Dstring re_ym ("[0-9].....");
   const Dstring re_ymd ("[0-9].......");
   const Dstring fmt ("gfs_3_[0-9]......._[0-9]..._[0-9]..\\.grb");

   typedef map<Grib::Key, Grib::Header*> Header_Ptr_Map;

   const vector<Dstring>& dir_ym = get_dir_listing (path, re_ym);
   for (vector<Dstring>::const_iterator i = dir_ym.begin ();
        i != dir_ym.end (); i++)
   {

      const Dstring& ym = *(i);
      const Dstring& ym_path = path + "/" + ym;

      const vector<Dstring>& dir_ymd = get_dir_listing (ym_path, re_ymd);
      for (vector<Dstring>::const_iterator j = dir_ymd.begin ();
           j != dir_ymd.end (); j++)
      {

         const Dstring& ymd = *(j);
         const Dstring& ymd_path = ym_path + "/" + ymd;

         const vector<Dstring>& dir_listing = get_dir_listing (ymd_path, fmt);
         for (vector<Dstring>::const_iterator iterator = dir_listing.begin ();
              iterator != dir_listing.end (); iterator++)
         {

            // fn = filename
            const Dstring& fn = *(iterator);
            const Dstring& file_path = ymd_path + "/" + fn;
            const Dstring& bt_str = fn.substr (6, 8) + fn.substr (15, 2);
            const Dstring& fh_str = fn.substr (20, 3);

            const Dtime base_time (bt_str, Dstring ("%Y%m%d%H"));
            const Integer forecast_hour = stoi (fh_str);
            const Key key (base_time, forecast_hour);

            key_multimap.add (key);
            Grib* grib_ptr = new Grib (file_path);
            insert (make_pair (key, grib_ptr));

         }

      }

   }

   if (size () > 0)
   {

      set<uint16_t> set_p;
      const Grib& grib = *(begin ()->second);
      const Header_Ptr_Map& header_ptr_map = grib.get_header_ptr_map ();

      for (Header_Ptr_Map::const_iterator iterator = header_ptr_map.begin ();
           iterator != header_ptr_map.end (); iterator++)
      {

         const Grib::Header& header = *(iterator->second);
         const Grib::Pds& pds = header.get_pds ();
         const Grib::Pds::Level& level = pds.get_level ();

         if (level.get_uint (0, 1) == 100)
         {
            const uint16_t p = level.get_uint (1, 2);
            if (p < 200) { continue; }
            set_p.insert (p);
         }

      }

      for (set<uint16_t>::const_iterator iterator = set_p.begin ();
           iterator != set_p.end (); iterator++)
      {
         const Real p = Real (*(iterator)) * 100;
         tuple_p.push_back (p);
      }

   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

   for (Key_Multimap::const_iterator iterator = key_multimap.begin ();
        iterator != key_multimap.end (); iterator++)
   {
      const Key key (iterator->first, iterator->second);
      initialize_3d_data (key);
   }

}

void
Gfs3::clear_3d_data ()
{
}

void
Gfs3::clean_up ()
{
   for (Gfs3::iterator iterator = begin (); iterator != end (); iterator++)
   {
      Grib* grib_ptr = iterator->second;
      delete grib_ptr;
   }
   clear ();
}

vector<Dtime>
Gfs3::get_valid_time_vector () const
{

   vector<Dtime> valid_time_vector;

   for (Gfs3::const_iterator iterator = begin ();
       iterator != end (); iterator++)
   {
      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime dtime (bt.t + fh);
      valid_time_vector.push_back (dtime);
   }

   return valid_time_vector;

/*
   vector<Dtime> valid_time_vector;
   Dtime start_time ("2012011300");
   Dtime end_time ("2012011400");
   for (Real t = start_time.t; t < end_time.t; t += 3)
   {
      Dtime dtime (t);
      valid_time_vector.push_back (dtime);
   }
   return valid_time_vector;
*/
}

Nwp::Key
Gfs3::get_key (const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs3::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());

   const Dtime& base_time = base_time_vector.back ();
   const Integer forecast_hour = Integer (round (dtime.t - base_time.t));

   return Key (base_time, forecast_hour);

}

void
Gfs3::acquire_base_time_forecast_hour (Dtime& base_time,
                                       Integer& forecast_hour,
                                       const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs3::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());
   base_time.t = base_time_vector.back ().t;
   forecast_hour = Integer (round (dtime.t - base_time.t));

}

Gfs4::Data_3D::Data_3D (const vector<Met_Element>& met_element_vector,
                        const Key& key)
   : Nwp::Data_3D (met_element_vector, key)
{
}

Real
Gfs4::Data_3D::evaluate (const Met_Element element,
                         const Real p,
                         const Real latitude,
                         const Real longitude,
                         const Evaluate_Op evaluate_op) const
{

   typedef Nwp::Data_3D Nd_3d;
   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& RH = denise::RELATIVE_HUMIDITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real rh = Nwp::Data_3D::evaluate (RH, p, latitude, longitude);
         return Moisture::get_t_d (t, rh);
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& O = denise::OMEGA;

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Evaluate_Op& eo = evaluate_op;
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real ts = Nd_3d::evaluate (T, p, latitude, longitude, eo);
            const Real os = Nd_3d::evaluate (O, p, latitude, longitude, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nd_3d::evaluate (element, p, latitude, longitude, eo);

}

Grib2::Key
Gfs4::get_grib_key (const Key& key,
                    const Met_Element met_element,
                    const Level& level) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   Grib2::Key grib_key;
   set_grib_key (grib_key, base_time);
   set_grib_key (grib_key, met_element, forecast_hour);
   set_grib_key (grib_key, met_element);
   set_grib_key (grib_key, met_element, level);
   return grib_key;

}

void
Gfs4::set_grib_key (Grib2::Key& grib_key,
                    const Dtime& base_time) const
{

   uint8_t* buffer = grib_key.buffer;

   const Integer yyyy = base_time.get_year ();
   const Integer mm = base_time.get_month ();
   const Integer dd = base_time.get_day ();
   const Integer HH = base_time.get_hour ();
   const Integer MM = base_time.get_minute ();
   const Integer SS = base_time.get_second ();

   buffer[0] = uint8_t (yyyy / 256);
   buffer[1] = uint8_t (yyyy % 256);
   buffer[2] = uint8_t (mm);
   buffer[3] = uint8_t (dd);
   buffer[4] = uint8_t (HH);
   buffer[5] = uint8_t (MM);
   buffer[6] = uint8_t (SS);

}

void
Gfs4::set_grib_key (Grib2::Key& grib_key,
                    const Met_Element met_element,
                    const Integer forecast_hour) const
{

   uint8_t* buffer = grib_key.buffer;
   buffer[7] = 1;

   switch (met_element)
   {

      default:
      {
         uint32_t value = forecast_hour;
#ifndef WORDS_BIGENDIAN
         swap_endian (&value, sizeof (uint32_t));
#endif
         memcpy (buffer + 8, &value, sizeof (uint32_t));
         break;
      }

      case PPT3:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 3);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case PPT6:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case PPTN:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = 0;
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case TOTAL_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case HIGH_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case MIDDLE_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case LOW_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

   }

}

void
Gfs4::set_grib_key (Grib2::Key& grib_key,
                    const Met_Element met_element) const
{

   uint8_t* buffer = grib_key.buffer;

   switch (met_element)
   {
      case PRESSURE:
         buffer[12] = 3;
         buffer[13] = 0;
         break;
      case TEMPERATURE:
         buffer[12] = 0;
         buffer[13] = 0;
         break;
      case TEMPERATURE_CELCIUS:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case MIX_DOWN_TEMPERATURE:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case ZONAL_WIND:
         buffer[12] = 2;
         buffer[13] = 2;
         break;
      case MERIDIONAL_WIND:
         buffer[12] = 2;
         buffer[13] = 3;
         break;
      case WIND_SPEED:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case VERTICAL_VELOCITY:
         buffer[12] = 2;
         buffer[13] = 9;
         break;
      case GEOPOTENTIAL_HEIGHT:
         buffer[12] = 3;
         buffer[13] = 5;
         break;
      case DEW_POINT:
         buffer[12] = 0;
         buffer[13] = 6;
         break;
      case MEAN_SEA_LEVEL_PRESSURE:
         buffer[12] = 3;
         buffer[13] = 1;
         break;
      case HIGH_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case MIDDLE_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case LOW_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case TOTAL_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case OMEGA:
         buffer[12] = 2;
         buffer[13] = 8;
         break;
      case PPT3:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case PPT6:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case PPTN:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case RAINFALL_STEP:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case RELATIVE_HUMIDITY:
         buffer[12] = 1;
         buffer[13] = 1;
         break;
      case DEW_POINT_DEPRESSION:
         buffer[12] = 0;
         buffer[13] = 7;
         break;
      case POTENTIAL_TEMPERATURE:
         buffer[12] = 0;
         buffer[13] = 2;
         break;
      case THETA_E:
         buffer[12] = 0;
         buffer[13] = 3;
         break;
      case MIXING_RATIO:
         buffer[12] = 1;
         buffer[13] = 2;
         break;
      case MONTGOMERY:
         buffer[12] = 2;
         buffer[13] = 6;
         break;
      case SLI:
         buffer[12] = 7;
         buffer[13] = 10;
         break;
      case SHOWALTER:
         buffer[12] = 7;
         buffer[13] = 13;
         break;
      case LI_700:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case LI_THUNDER:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case K_INDEX:
         buffer[12] = 7;
         buffer[13] = 2;
         break;
      case TOTAL_TOTALS:
         buffer[12] = 7;
         buffer[13] = 4;
         break;
      case CAPE:
         buffer[12] = 7;
         buffer[13] = 6;
         break;
      case PRECIPITABLE_WATER:
         buffer[12] = 1;
         buffer[13] = 3;
         break;
      case FOG_FRACTION:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case THICKNESS:
         buffer[12] = 3;
         buffer[13] = 12;
         break;
      case POTENTIAL_VORTICITY:
         buffer[12] = 2;
         buffer[13] = 14;
         break;
      case ABSOLUTE_VORTICITY:
         buffer[12] = 2;
         buffer[13] = 10;
         break;
      case SHEAR_VORTICITY:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case CURVATURE_VORTICITY:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
   }

}

void
Gfs4::set_grib_key (Grib2::Key& grib_key,
                    const Met_Element met_element,
                    const denise::Level& level) const
{

   uint8_t* buffer = grib_key.buffer;

   switch (met_element)
   {
      case TOTAL_CLOUD:
         buffer[14] = 200;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case HIGH_CLOUD:
         buffer[14] = 214;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case MIDDLE_CLOUD:
         buffer[14] = 224;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case LOW_CLOUD:
         buffer[14] = 234;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case PRECIPITABLE_WATER:
         buffer[13] = 200;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
   }

   switch (level.type)
   {
      case Level::PRESSURE:
      {
         uint32_t p = uint32_t (round (level.value));
         buffer[14] = 100;
         buffer[15] = 0;
#ifndef WORDS_BIGENDIAN
         swap_endian (&p, sizeof (uint32_t));
#endif
         memcpy (buffer + 16, &p, sizeof (uint32_t));
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::THETA:
      {
         uint32_t theta = uint32_t (round (level.value));
         buffer[14] = 107;
         buffer[15] = 0;
#ifndef WORDS_BIGENDIAN
         swap_endian (&theta, sizeof (uint32_t));
#endif
         memcpy (buffer + 16, &theta, sizeof (uint32_t));
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::SIGMA:
      {
         uint32_t sigma = uint32_t (round (level.value * 10000));
         buffer[14] = 113;
         buffer[15] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::SCREEN:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 2;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::FIFTY_METRE:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 50;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::TEN_METRE:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 10;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::MEAN_SEA:
      {
         buffer[14] = 101;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::SURFACE:
      {
         buffer[14] = 1;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::NIL:
      case Level::NAL:
      {
         buffer[14] = 1;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
   }

}

void
Gfs4::initialize_3d_data (const Key& key)
{
   typedef Gfs4::Data_3D G3d_3d;
   G3d_3d* g3d_3d_ptr = new G3d_3d (met_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, g3d_3d_ptr));
}

void          
Gfs4::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Met_Element>::const_iterator Iterator;
   const Key& key = data_3d.key;

   for (Iterator iterator = met_element_vector.begin ();
        iterator != met_element_vector.end (); iterator++)
   {
      const Met_Element& met_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = get_gvd_3d_ptr (met_element, key);
      data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
   }

   data_3d.set_available ();

}

Geodetic_Vector_Data_2D*
Gfs4::get_initialized_vd_2d (const Integer vector_size) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   Gvd_2d* data_ptr = new Gvd_2d (vector_size, size_2d, domain_2d, false);
   return data_ptr;

}

void 
Gfs4::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Level& level,
                              const Met_Element met_element)
{

   const Integer vi = vector_index;

   switch (met_element)
   {

      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& surface = Level::surface_level ();
         fill_grib_data (gvd_2d, vi, met_element, key, surface);
         return;
      }

   }

   Nwp::fill_ts_diagnosis_data (gvd_2d, vi, key, level, met_element);

}

void
Gfs4::fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key)
{

   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour < 0)
   {
      throw Nwp_Exception ("Forecast Hour < 0");
      return;
   }

   gvd_2d.initialize (vector_index, 0);
   if (forecast_hour == 0) { return; }

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();

   Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
   Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

   for (Integer fh = 0; fh < forecast_hour; fh += 6)
   {
      fill_grib_data (precip_data, 0, PPT6, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (precip_data, 0, PPT3, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   delete precip_data_ptr;

}

void
Gfs4::fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();
   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour == 0)
   {
      gvd_2d.initialize (vector_index, 0);
      return;
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (gvd_2d, vector_index, PPT3, key, surface_level);
   }
   else
   {

      fill_grib_data (gvd_2d, vector_index, PPT6, key, surface_level);

      Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
      Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

      try
      {

         const Dtime& base_time = key.base_time;
         const Integer forecast_hour = key.forecast_hour;
         const Key prev_key (base_time, forecast_hour - 3);

         fill_grib_data (precip_data, 0, PPT3, prev_key, surface_level);
         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real ppt3 = precip_data.get_datum (0, i, j);
               gvd_2d.get_datum (vector_index, i, j) -= ppt3;
            }
         }
      }
      catch (const Nwp_Exception& ne)
      {
      }

      delete precip_data_ptr;

   }

}

void
Gfs4::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Met_Element met_element)
{

   switch (met_element)
   {

      case RAINFALL_CUMULATIVE:
      {
         fill_cumulative_rain_data (gvd_2d, vector_index, key);
         return;
      }

      case RAINFALL_STEP:
      {
         fill_step_rain_data (gvd_2d, vector_index, key);
         return;
      }

   }

   Nwp::fill_rain_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs4::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                       const Integer vector_index,
                       const Key& key,
                       const Met_Element met_element)
{
   const Level& nil_level = Level::nil_level ();
   fill_grib_data (gvd_2d, vector_index, met_element, key, nil_level);
}

void
Gfs4::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Met_Element met_element)
{

   const Level& screen = Level::screen_level ();
   const Level& surface = Level::surface_level ();

   switch (met_element)
   {

      case TEMPERATURE:
      {
         fill_grib_data (gvd_2d, vector_index, TEMPERATURE, key, screen);
         return;
      }

      case DEW_POINT:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, RELATIVE_HUMIDITY);

         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real rh = data_ptr->get_datum (1, i, j);
               const Real t_d = Moisture::get_t_d (t, rh);
               gvd_2d.set_datum (vector_index, i, j, t_d);
            }
         }

         delete data_ptr;
         return;

      }

      case RELATIVE_HUMIDITY:
      {
         fill_grib_data (gvd_2d, vector_index, RELATIVE_HUMIDITY, key, screen);
         return;
      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs4::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key,
                           const Met_Element met_element)
{

   const Level& ten = Level::ten_metre_level ();

   switch (met_element)
   {

      case ZONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, ZONAL_WIND, key, ten);
         return;
      }

      case MERIDIONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, MERIDIONAL_WIND, key, ten);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs4::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                     const Integer vector_index,
                     const Key& key,
                     const Met_Element met_element)
{

   const Level& msl = Level::mean_sea_level ();

   switch (met_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         const Met_Element mslp = MEAN_SEA_LEVEL_PRESSURE;
         fill_grib_data (gvd_2d, vector_index, mslp, key, msl);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs4::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                               const Integer vector_index,
                               const Key& key,
                               const Met_Element met_element)
{

   const Level& surface = Level::surface_level ();

   switch (met_element)
   {

      case PRESSURE:
      {
         fill_grib_data (gvd_2d, vector_index, PRESSURE, key, surface);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, met_element);

}

/*
void
Gfs4::fill_grib_data (Geodetic_Vector_Data_3D& gvd_3d,
                      const Integer vector_index,
                      const Met_Element met_element,
                      const Key& key) const
{

   const Key key (key);

   Gfs4::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib2& grib = *(iterator->second);

   for (Integer k = 0; k < tuple_p.size (); k++)
   {
      const Real p = tuple_p[k];
      const Level level (Level::PRESSURE, p);
      const Grib2::Key& key = get_grib_key (key, met_element, level);
      grib.fill_data (gvd_3d, vector_index, k, key);
   }

}
*/

Geodetic_Vector_Data_3D*
Gfs4::get_gvd_3d_ptr (const Met_Element met_element,
                      const Key& key) const
{

   Geodetic_Vector_Data_3D* gvd_3d_ptr =
      new Geodetic_Vector_Data_3D (1, tuple_p, size_2d, domain_2d);
   Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   for (Integer k = 0; k < tuple_p.size (); k++)
   {

      const Real p = tuple_p[k];
      const denise::Level level (Level::PRESSURE, p);
      const Grib2::Key& grib_key = get_grib_key (key, met_element, level);

      Gfs4::const_iterator iterator = find (key);
      if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
      const Grib2& grib = *(iterator->second);

      Geodetic_Vector_Data_2D* grib_data_ptr =
         get_grib_data_ptr (grib, grib_key);

      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         const Real latitude = grib_data_ptr->get_latitude (i);
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real longitude = grib_data_ptr->get_longitude (j);
            const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
            gvd_3d.set_datum (0, k, i, j, datum);
         }
      }

      delete grib_data_ptr;

   }

   return gvd_3d_ptr;

}

Geodetic_Vector_Data_2D*
Gfs4::get_grib_data_ptr (const Grib2& grib,
                         const Grib2::Key& grib_key) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = new Gvd_2d (1, grib_size_2d, grib_domain_2d, true);
   grib.fill_data (*data_ptr, 0, grib_key);

   const Integer last_j = grib_size_2d.j - 1;

   #pragma omp parallel for
   for (Integer i = 0; i < grib_size_2d.i; i++)
   {
      const Real datum = data_ptr->get_datum (0, i, 0);
      data_ptr->set_datum (0, i, last_j, datum);
   }

   return data_ptr;

}

void
Gfs4::fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Met_Element met_element,
                      const Key& key,
                      const Level& level) const
{

   Gfs4::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib2& grib = *(iterator->second);

   const Grib2::Key& grib_key = get_grib_key (key, met_element, level);
   Geodetic_Vector_Data_2D* grib_data_ptr = get_grib_data_ptr (grib, grib_key);

   Geodetic_Vector_Data_2D grib_data (1, grib_size_2d, grib_domain_2d, true);

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_latitude (i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_longitude (j);
         const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete grib_data_ptr;

}

Gfs4::Gfs4 (const Dstring& data_path)
   : Nwp ("Gfs4", data_path),
     data_path (data_path),
     grib_size_2d (361, 720),
     grib_domain_2d (Domain_1D (-90, 90), Domain_1D (0, 360)),
     size_2d (51, 71),
     domain_2d (Domain_1D (-60, 0), Domain_1D (80, 160))
{

   this->status = "Unloaded";

   met_element_vector.push_back (denise::TEMPERATURE);
   met_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   met_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   met_element_vector.push_back (denise::ZONAL_WIND);
   met_element_vector.push_back (denise::MERIDIONAL_WIND);
   met_element_vector.push_back (denise::OMEGA);

}

Gfs4::~Gfs4 ()
{
   clean_up ();
}

void
Gfs4::survey ()
{

   const Dstring re_ym ("[0-9].....");
   const Dstring re_ymd ("[0-9].......");
   const Dstring fmt ("gfs_4_[0-9]......._[0-9]..._[0-9]..\\.grb");

   typedef map<Grib2::Key, Grib2::Header*> Header_Ptr_Map;

   const vector<Dstring>& dir_ym = get_dir_listing (path, re_ym);
   for (vector<Dstring>::const_iterator i = dir_ym.begin ();
        i != dir_ym.end (); i++)
   {

      const Dstring& ym = *(i);
      const Dstring& ym_path = path + "/" + ym;

      const vector<Dstring>& dir_ymd = get_dir_listing (ym_path, re_ymd);
      for (vector<Dstring>::const_iterator j = dir_ymd.begin ();
           j != dir_ymd.end (); j++)
      {

         const Dstring& ymd = *(j);
         const Dstring& ymd_path = ym_path + "/" + ymd;

         const vector<Dstring>& dir_listing = get_dir_listing (ymd_path, fmt);
         for (vector<Dstring>::const_iterator iterator = dir_listing.begin ();
              iterator != dir_listing.end (); iterator++)
         {

            // fn = filename
            const Dstring& fn = *(iterator);
            const Dstring& file_path = ymd_path + "/" + fn;
            const Dstring& bt_str = fn.substr (6, 8) + fn.substr (15, 2);
            const Dstring& fh_str = fn.substr (20, 3);

            const Dtime base_time (bt_str, "%Y%m%d%H");
            const Integer forecast_hour = stoi (fh_str);

            const Key key (base_time, forecast_hour);
            key_multimap.add (key);
            Grib2* grib_ptr = new Grib2 (file_path);
            insert (make_pair (key, grib_ptr));
            //valid_time_set.insert (valid_time);
            //key_set.insert (Key (base_time, fh));

         }

      }

   }

   if (size () > 0)
   {

      set<uint32_t> set_p;
      const Grib2& grib = *(begin ()->second);
      const Header_Ptr_Map& header_ptr_map = grib.get_header_ptr_map ();

      for (Header_Ptr_Map::const_iterator iterator = header_ptr_map.begin ();
           iterator != header_ptr_map.end (); iterator++)
      {

         const Grib2::Header& header = *(iterator->second);
         const Grib2::Block_4& block_4 = header.block_4;
         const Grib2::Block_4::Level& level = block_4.get_level ();

         if (level.get_first_level_type () == 100)
         {
            const uint32_t p = level.get_first_level ();
            if (p < 200) { continue; }
            set_p.insert (p);
         }

      }

      for (set<uint32_t>::const_iterator iterator = set_p.begin ();
           iterator != set_p.end (); iterator++)
      {
         const Real p = Real (*(iterator));
         tuple_p.push_back (p);
      }

   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

}

void
Gfs4::clean_up ()
{
   for (Gfs4::iterator iterator = begin (); iterator != end (); iterator++)
   {
      Grib2* grib_ptr = iterator->second;
      delete grib_ptr;
   }
}

void
Gfs4::set_domain_2d (const Domain_2D& domain_2d)
{

   this->domain_2d.domain_x.start = ceil (domain_2d.domain_x.start);
   this->domain_2d.domain_x.end = floor (domain_2d.domain_x.end);
   this->domain_2d.domain_y.start = ceil (domain_2d.domain_y.start);
   this->domain_2d.domain_y.end = floor (domain_2d.domain_y.end);

   const Real& start_latitude = this->domain_2d.domain_x.start;
   const Real& end_latitude = this->domain_2d.domain_x.end;
   const Real& start_longitude = this->domain_2d.domain_y.start;
   const Real& end_longitude = this->domain_2d.domain_y.end;
 
   size_2d.i = Integer (round (end_latitude - start_latitude));
   size_2d.j = Integer (round (end_longitude - start_longitude));

}

vector<Dtime>
Gfs4::get_valid_time_vector () const
{

   vector<Dtime> valid_time_vector;

   for (Gfs4::const_iterator iterator = begin ();
       iterator != end (); iterator++)
   {
      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime dtime (bt.t + fh);
      valid_time_vector.push_back (dtime);
   }

   return valid_time_vector;

/*
   vector<Dtime> valid_time_vector;
   Dtime start_time ("2012011300");
   Dtime end_time ("2012011400");
   for (Real t = start_time.t; t < end_time.t; t += 3)
   {
      Dtime dtime (t);
      valid_time_vector.push_back (dtime);
   }
   return valid_time_vector;
*/
}

Nwp::Key
Gfs4::get_key (const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs4::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());

   const Dtime& base_time = base_time_vector.back ();
   const Integer forecast_hour = Integer (round (dtime.t - base_time.t));

   return Key (base_time, forecast_hour);

}

void
Gfs4::acquire_base_time_forecast_hour (Dtime& base_time,
                                       Integer& forecast_hour,
                                       const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs4::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());
   base_time.t = base_time_vector.back ().t;
   forecast_hour = Integer (round (dtime.t - base_time.t));

}

Gfs::Data_3D::Data_3D (const vector<Met_Element>& met_element_vector,
                       const Key& key)
   : Nwp::Data_3D (met_element_vector, key)
{
}

Real
Gfs::Data_3D::evaluate (const Met_Element element,
                        const Real p,
                        const Real latitude,
                        const Real longitude,
                        const Evaluate_Op evaluate_op) const
{

   typedef Nwp::Data_3D Nd_3d;
   const Evaluate_Op& eo = evaluate_op;

   switch (element)
   {

      case denise::DEW_POINT:
      {
         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& RH = denise::RELATIVE_HUMIDITY;
         const Real t = Nwp::Data_3D::evaluate (T, p, latitude, longitude);
         const Real rh = Nwp::Data_3D::evaluate (RH, p, latitude, longitude);
         return Moisture::get_t_d (t, rh);
      }

      case denise::VERTICAL_VELOCITY:
      {

         const Met_Element& T = denise::TEMPERATURE;
         const Met_Element& O = denise::OMEGA;

         if (evaluate_op == DX || evaluate_op == DY)
         {
            const Evaluate_Op& eo = evaluate_op;
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real ts = Nd_3d::evaluate (T, p, latitude, longitude, eo);
            const Real os = Nd_3d::evaluate (O, p, latitude, longitude, eo);
            const Real rho = p / (R_d * t);
            const Real oRdp = o * R_d / p;
            return (oRdp * ts + os / rho) / -g;
         }
         else
         {
            const Real t = Nd_3d::evaluate (T, p, latitude, longitude);
            const Real o = Nd_3d::evaluate (O, p, latitude, longitude);
            const Real rho = p / (R_d * t);
            return o / (-rho * g);
         }

      }

   }

   return Nd_3d::evaluate (element, p, latitude, longitude, eo);

}

Grib2::Key
Gfs::get_grib_key (const Key& key,
                   const Met_Element met_element,
                   const Level& level) const
{

   const Dtime& base_time = key.base_time;
   const Integer forecast_hour = key.forecast_hour;

   Grib2::Key grib_key;
   set_grib_key (grib_key, base_time);
   set_grib_key (grib_key, met_element, forecast_hour);
   set_grib_key (grib_key, met_element);
   set_grib_key (grib_key, met_element, level);
   return grib_key;

}

void
Gfs::set_grib_key (Grib2::Key& grib_key,
                   const Dtime& base_time) const
{

   uint8_t* buffer = grib_key.buffer;

   const Integer yyyy = base_time.get_year ();
   const Integer mm = base_time.get_month ();
   const Integer dd = base_time.get_day ();
   const Integer HH = base_time.get_hour ();
   const Integer MM = base_time.get_minute ();
   const Integer SS = base_time.get_second ();

   buffer[0] = uint8_t (yyyy / 256);
   buffer[1] = uint8_t (yyyy % 256);
   buffer[2] = uint8_t (mm);
   buffer[3] = uint8_t (dd);
   buffer[4] = uint8_t (HH);
   buffer[5] = uint8_t (MM);
   buffer[6] = uint8_t (SS);

}

void
Gfs::set_grib_key (Grib2::Key& grib_key,
                   const Met_Element met_element,
                   const Integer forecast_hour) const
{

   uint8_t* buffer = grib_key.buffer;
   buffer[7] = 1;

   switch (met_element)
   {

      default:
      {
         uint32_t value = forecast_hour;
#ifndef WORDS_BIGENDIAN
         swap_endian (&value, sizeof (uint32_t));
#endif
         memcpy (buffer + 8, &value, sizeof (uint32_t));
         break;
      }

      case PPT3:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 3);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case PPT6:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case PPTN:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = 0;
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 4;
         break;
      }

      case TOTAL_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case HIGH_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case MIDDLE_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

      case LOW_CLOUD:
      {
         throw Nwp_Exception ("Not Available");
         buffer[7] = uint8_t (forecast_hour - 6);
         buffer[8] = uint8_t (forecast_hour);
         buffer[9] = 3;
         break;
      }

   }

}

void
Gfs::set_grib_key (Grib2::Key& grib_key,
                   const Met_Element met_element) const
{

   uint8_t* buffer = grib_key.buffer;

   switch (met_element)
   {
      case PRESSURE:
         buffer[12] = 3;
         buffer[13] = 0;
         break;
      case TEMPERATURE:
         buffer[12] = 0;
         buffer[13] = 0;
         break;
      case TEMPERATURE_CELCIUS:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case MIX_DOWN_TEMPERATURE:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case ZONAL_WIND:
         buffer[12] = 2;
         buffer[13] = 2;
         break;
      case MERIDIONAL_WIND:
         buffer[12] = 2;
         buffer[13] = 3;
         break;
      case WIND_SPEED:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case VERTICAL_VELOCITY:
         buffer[12] = 2;
         buffer[13] = 9;
         break;
      case GEOPOTENTIAL_HEIGHT:
         buffer[12] = 3;
         buffer[13] = 5;
         break;
      case DEW_POINT:
         buffer[12] = 0;
         buffer[13] = 6;
         break;
      case MEAN_SEA_LEVEL_PRESSURE:
         buffer[12] = 3;
         buffer[13] = 1;
         break;
      case HIGH_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case MIDDLE_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case LOW_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case TOTAL_CLOUD:
         buffer[12] = 6;
         buffer[13] = 1;
         break;
      case OMEGA:
         buffer[12] = 2;
         buffer[13] = 8;
         break;
      case PPT3:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case PPT6:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case PPTN:
         buffer[12] = 1;
         buffer[13] = 8;
         break;
      case RAINFALL_STEP:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case RELATIVE_HUMIDITY:
         buffer[12] = 1;
         buffer[13] = 1;
         break;
      case DEW_POINT_DEPRESSION:
         buffer[12] = 0;
         buffer[13] = 7;
         break;
      case POTENTIAL_TEMPERATURE:
         buffer[12] = 0;
         buffer[13] = 2;
         break;
      case THETA_E:
         buffer[12] = 0;
         buffer[13] = 3;
         break;
      case MIXING_RATIO:
         buffer[12] = 1;
         buffer[13] = 2;
         break;
      case MONTGOMERY:
         buffer[12] = 2;
         buffer[13] = 6;
         break;
      case SLI:
         buffer[12] = 7;
         buffer[13] = 10;
         break;
      case SHOWALTER:
         buffer[12] = 7;
         buffer[13] = 13;
         break;
      case LI_700:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case LI_THUNDER:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case K_INDEX:
         buffer[12] = 7;
         buffer[13] = 2;
         break;
      case TOTAL_TOTALS:
         buffer[12] = 7;
         buffer[13] = 4;
         break;
      case CAPE:
         buffer[12] = 7;
         buffer[13] = 6;
         break;
      case PRECIPITABLE_WATER:
         buffer[12] = 1;
         buffer[13] = 3;
         break;
      case FOG_FRACTION:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case THICKNESS:
         buffer[12] = 3;
         buffer[13] = 12;
         break;
      case POTENTIAL_VORTICITY:
         buffer[12] = 2;
         buffer[13] = 14;
         break;
      case ABSOLUTE_VORTICITY:
         buffer[12] = 2;
         buffer[13] = 10;
         break;
      case SHEAR_VORTICITY:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
      case CURVATURE_VORTICITY:
         buffer[12] = 255;
         buffer[13] = 255;
         break;
   }

}

void
Gfs::set_grib_key (Grib2::Key& grib_key,
                   const Met_Element met_element,
                   const denise::Level& level) const
{

   uint8_t* buffer = grib_key.buffer;

   switch (met_element)
   {
      case TOTAL_CLOUD:
         buffer[14] = 200;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case HIGH_CLOUD:
         buffer[14] = 214;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case MIDDLE_CLOUD:
         buffer[14] = 224;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case LOW_CLOUD:
         buffer[14] = 234;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         return;
      case PRECIPITABLE_WATER:
         buffer[13] = 200;
         buffer[14] = 0;
         buffer[15] = 0;
         return;
   }

   switch (level.type)
   {
      case Level::PRESSURE:
      {
         uint32_t p = uint32_t (round (level.value));
         buffer[14] = 100;
         buffer[15] = 0;
#ifndef WORDS_BIGENDIAN
         swap_endian (&p, sizeof (uint32_t));
#endif
         memcpy (buffer + 16, &p, sizeof (uint32_t));
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::THETA:
      {
         uint32_t theta = uint32_t (round (level.value));
         buffer[14] = 107;
         buffer[15] = 0;
#ifndef WORDS_BIGENDIAN
         swap_endian (&theta, sizeof (uint32_t));
#endif
         memcpy (buffer + 16, &theta, sizeof (uint32_t));
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::SIGMA:
      {
         uint32_t sigma = uint32_t (round (level.value * 10000));
         buffer[14] = 113;
         buffer[15] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::SCREEN:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 2;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::FIFTY_METRE:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 50;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::TEN_METRE:
      {
         buffer[14] = 103;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 10;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::MEAN_SEA:
      {
         buffer[14] = 101;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::SURFACE:
      {
         buffer[14] = 1;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
      case Level::NIL:
      case Level::NAL:
      {
         buffer[14] = 1;
         buffer[15] = 0;
         buffer[16] = 0;
         buffer[17] = 0;
         buffer[18] = 0;
         buffer[19] = 0;
         buffer[20] = 255;
         buffer[21] = 0;
         buffer[22] = 0;
         buffer[23] = 0;
         buffer[24] = 0;
         buffer[25] = 0;
         break;
      }
   }

}

void
Gfs::initialize_3d_data (const Key& key)
{
   typedef Gfs::Data_3D G3d_3d;
   G3d_3d* g3d_3d_ptr = new G3d_3d (met_element_vector, key);
   data_3d_ptr_map.insert (make_pair (key, g3d_3d_ptr));
}

void          
Gfs::load_3d_data (Nwp::Data_3D& data_3d)
{

   typedef vector<Met_Element>::const_iterator Iterator;
   const Key& key = data_3d.key;

   for (Iterator iterator = met_element_vector.begin ();
        iterator != met_element_vector.end (); iterator++)
   {
      const Met_Element& met_element = *(iterator);
      Geodetic_Vector_Data_3D* gvd_3d_ptr = get_gvd_3d_ptr (met_element, key);
      data_3d.set_gvd_3d_ptr (met_element, gvd_3d_ptr);
   }

   data_3d.set_available ();

}

Geodetic_Vector_Data_2D*
Gfs::get_initialized_vd_2d (const Integer vector_size) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;

   Gvd_2d* data_ptr = new Gvd_2d (vector_size, size_2d, domain_2d, false);
   return data_ptr;

}

void 
Gfs::fill_ts_diagnosis_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Level& level,
                              const Met_Element met_element)
{

   const Integer vi = vector_index;

   switch (met_element)
   {

      case CAPE:
      case PRECIPITABLE_WATER:
      {
         const Level& surface = Level::surface_level ();
         fill_grib_data (gvd_2d, vi, met_element, key, surface);
         return;
      }

   }

   Nwp::fill_ts_diagnosis_data (gvd_2d, vi, key, level, met_element);

}

void
Gfs::fill_cumulative_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                                 const Integer vector_index,
                                 const Key& key)
{

   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour < 0)
   {
      throw Nwp_Exception ("Forecast Hour < 0");
      return;
   }

   gvd_2d.initialize (vector_index, 0);
   if (forecast_hour == 0) { return; }

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();

   Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
   Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

   for (Integer fh = 0; fh < forecast_hour; fh += 6)
   {
      fill_grib_data (precip_data, 0, PPT6, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (precip_data, 0, PPT3, key, surface_level);
      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real precip = precip_data.get_datum (0, i, j);
            gvd_2d.get_datum (vector_index, i, j) += precip;
         }
      }
   }

   delete precip_data_ptr;

}

void
Gfs::fill_step_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                           const Integer vector_index,
                           const Key& key)
{

   const Size_2D& size_2d = gvd_2d.get_size_2d ();
   const Level& surface_level = Level::surface_level ();
   const Integer forecast_hour = key.forecast_hour;

   if (forecast_hour == 0)
   {
      gvd_2d.initialize (vector_index, 0);
      return;
   }

   if (forecast_hour % 6 != 0)
   {
      fill_grib_data (gvd_2d, vector_index, PPT3, key, surface_level);
   }
   else
   {

      fill_grib_data (gvd_2d, vector_index, PPT6, key, surface_level);

      Geodetic_Vector_Data_2D* precip_data_ptr = get_initialized_vd_2d (1);
      Geodetic_Vector_Data_2D& precip_data = *precip_data_ptr;

      try
      {

         const Dtime& base_time = key.base_time;
         const Integer forecast_hour = key.forecast_hour;
         const Key prev_key (base_time, forecast_hour - 3);

         fill_grib_data (precip_data, 0, PPT3, prev_key, surface_level);
         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real ppt3 = precip_data.get_datum (0, i, j);
               gvd_2d.get_datum (vector_index, i, j) -= ppt3;
            }
         }
      }
      catch (const Nwp_Exception& ne)
      {
      }

      delete precip_data_ptr;

   }

}

void
Gfs::fill_rain_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Met_Element met_element)
{

   switch (met_element)
   {

      case RAINFALL_CUMULATIVE:
      {
         fill_cumulative_rain_data (gvd_2d, vector_index, key);
         return;
      }

      case RAINFALL_STEP:
      {
         fill_step_rain_data (gvd_2d, vector_index, key);
         return;
      }

   }

   Nwp::fill_rain_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs::fill_cloud_data (Geodetic_Vector_Data_2D& gvd_2d,
                      const Integer vector_index,
                      const Key& key,
                      const Met_Element met_element)
{
   const Level& nil_level = Level::nil_level ();
   fill_grib_data (gvd_2d, vector_index, met_element, key, nil_level);
}

void
Gfs::fill_screen_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                             const Integer vector_index,
                             const Key& key,
                             const Met_Element met_element)
{

   const Level& screen = Level::screen_level ();
   const Level& surface = Level::surface_level ();

   switch (met_element)
   {

      case TEMPERATURE:
      {
         fill_grib_data (gvd_2d, vector_index, TEMPERATURE, key, screen);
         return;
      }

      case DEW_POINT:
      {

         typedef Geodetic_Vector_Data_2D Gvd_2d;
         Gvd_2d* data_ptr = get_initialized_vd_2d (2);
         fill_data (*data_ptr, 0, key, screen, TEMPERATURE);
         fill_data (*data_ptr, 1, key, screen, RELATIVE_HUMIDITY);

         #pragma omp parallel for
         for (Integer i = 0; i < size_2d.i; i++)
         {
            for (Integer j = 0; j < size_2d.j; j++)
            {
               const Real t = data_ptr->get_datum (0, i, j);
               const Real rh = data_ptr->get_datum (1, i, j);
               const Real t_d = Moisture::get_t_d (t, rh);
               gvd_2d.set_datum (vector_index, i, j, t_d);
            }
         }

         delete data_ptr;
         return;

      }

      case RELATIVE_HUMIDITY:
      {
         fill_grib_data (gvd_2d, vector_index, RELATIVE_HUMIDITY, key, screen);
         return;
      }

   }

   Nwp::fill_screen_level_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs::fill_10m_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                          const Integer vector_index,
                          const Key& key,
                          const Met_Element met_element)
{

   const Level& ten = Level::ten_metre_level ();

   switch (met_element)
   {

      case ZONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, ZONAL_WIND, key, ten);
         return;
      }

      case MERIDIONAL_WIND:
      {
         fill_grib_data (gvd_2d, vector_index, MERIDIONAL_WIND, key, ten);
         return;
      }

   }

   Nwp::fill_10m_level_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs::fill_msl_data (Geodetic_Vector_Data_2D& gvd_2d,
                    const Integer vector_index,
                    const Key& key,
                    const Met_Element met_element)
{

   const Level& msl = Level::mean_sea_level ();

   switch (met_element)
   {

      case PRESSURE:
      case MEAN_SEA_LEVEL_PRESSURE:
      {
         const Met_Element mslp = MEAN_SEA_LEVEL_PRESSURE;
         fill_grib_data (gvd_2d, vector_index, mslp, key, msl);
         return;
      }

   }

   Nwp::fill_msl_data (gvd_2d, vector_index, key, met_element);

}

void
Gfs::fill_surface_level_data (Geodetic_Vector_Data_2D& gvd_2d,
                              const Integer vector_index,
                              const Key& key,
                              const Met_Element met_element)
{

   const Level& surface = Level::surface_level ();

   switch (met_element)
   {

      case PRESSURE:
      {
         fill_grib_data (gvd_2d, vector_index, PRESSURE, key, surface);
         return;
      }

   }

   Nwp::fill_surface_level_data (gvd_2d, vector_index, key, met_element);

}

/*
void
Gfs::fill_grib_data (Geodetic_Vector_Data_3D& gvd_3d,
                     const Integer vector_index,
                     const Met_Element met_element,
                     const Key& key) const
{

   const Key key (key);

   Gfs::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib2& grib = *(iterator->second);

   for (Integer k = 0; k < tuple_p.size (); k++)
   {
      const Real p = tuple_p[k];
      const Level level (Level::PRESSURE, p);
      const Grib2::Key& key = get_grib_key (key, met_element, level);
      grib.fill_data (gvd_3d, vector_index, k, key);
   }

}
*/

Geodetic_Vector_Data_3D*
Gfs::get_gvd_3d_ptr (const Met_Element met_element,
                     const Key& key) const
{

   Geodetic_Vector_Data_3D* gvd_3d_ptr =
      new Geodetic_Vector_Data_3D (1, tuple_p, size_2d, domain_2d);
   Geodetic_Vector_Data_3D& gvd_3d = *gvd_3d_ptr;

   for (Integer k = 0; k < tuple_p.size (); k++)
   {

      const Real p = tuple_p[k];
      const denise::Level level (Level::PRESSURE, p);
      const Grib2::Key& grib_key = get_grib_key (key, met_element, level);

      Gfs::const_iterator iterator = find (key);
      if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
      const Grib2& grib = *(iterator->second);

      Geodetic_Vector_Data_2D* grib_data_ptr =
         get_grib_data_ptr (grib, grib_key);

      #pragma omp parallel for
      for (Integer i = 0; i < size_2d.i; i++)
      {
         const Real latitude = grib_data_ptr->get_latitude (i);
         for (Integer j = 0; j < size_2d.j; j++)
         {
            const Real longitude = grib_data_ptr->get_longitude (j);
            const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
            gvd_3d.set_datum (0, k, i, j, datum);
         }
      }

      delete grib_data_ptr;

   }

   return gvd_3d_ptr;

}

Geodetic_Vector_Data_2D*
Gfs::get_grib_data_ptr (const Grib2& grib,
                        const Grib2::Key& grib_key) const
{

   typedef Geodetic_Vector_Data_2D Gvd_2d;
   Gvd_2d* data_ptr = new Gvd_2d (1, grib_size_2d, grib_domain_2d, true);
   grib.fill_data (*data_ptr, 0, grib_key);

   const Integer last_j = grib_size_2d.j - 1;

   #pragma omp parallel for
   for (Integer i = 0; i < grib_size_2d.i; i++)
   {
      const Real datum = data_ptr->get_datum (0, i, 0);
      data_ptr->set_datum (0, i, last_j, datum);
   }

   return data_ptr;

}

void
Gfs::fill_grib_data (Geodetic_Vector_Data_2D& gvd_2d,
                     const Integer vector_index,
                     const Met_Element met_element,
                     const Key& key,
                     const Level& level) const
{

   Gfs::const_iterator iterator = find (key);
   if (iterator == end ()) { throw Nwp_Exception ("Not Available"); }
   const Grib2& grib = *(iterator->second);

   const Grib2::Key& grib_key = get_grib_key (key, met_element, level);
   Geodetic_Vector_Data_2D* grib_data_ptr = get_grib_data_ptr (grib, grib_key);

   Geodetic_Vector_Data_2D grib_data (1, grib_size_2d, grib_domain_2d, true);

   #pragma omp parallel for
   for (Integer i = 0; i < size_2d.i; i++)
   {
      const Real latitude = gvd_2d.get_latitude (i);
      for (Integer j = 0; j < size_2d.j; j++)
      {
         const Real longitude = gvd_2d.get_longitude (j);
         const Real datum = grib_data_ptr->evaluate (0, latitude, longitude);
         gvd_2d.set_datum (vector_index, i, j, datum);
      }
   }

   delete grib_data_ptr;

}

Gfs::Gfs (const Dstring& data_path)
   : Nwp ("Gfs", data_path),
     data_path (data_path),
     grib_size_2d (361, 720),
     grib_domain_2d (Domain_1D (-90, 90), Domain_1D (0, 360)),
     size_2d (51, 71),
     domain_2d (Domain_1D (-60, 0), Domain_1D (80, 160))
{

   this->status = "Unloaded";

   met_element_vector.push_back (denise::TEMPERATURE);
   met_element_vector.push_back (denise::RELATIVE_HUMIDITY);
   met_element_vector.push_back (denise::GEOPOTENTIAL_HEIGHT);
   met_element_vector.push_back (denise::ZONAL_WIND);
   met_element_vector.push_back (denise::MERIDIONAL_WIND);
   met_element_vector.push_back (denise::OMEGA);

}

Gfs::~Gfs ()
{
   clean_up ();
}

void
Gfs::survey ()
{

   const Dstring re_ym ("[0-9].....");
   const Dstring re_ymd ("[0-9].......");
   const Dstring fmt ("gfs_4_[0-9]......._[0-9]..._[0-9]..\\.grb");

   typedef map<Grib2::Key, Grib2::Header*> Header_Ptr_Map;

   const vector<Dstring>& dir_ym = get_dir_listing (path, re_ym);
   for (vector<Dstring>::const_iterator i = dir_ym.begin ();
        i != dir_ym.end (); i++)
   {

      const Dstring& ym = *(i);
      const Dstring& ym_path = path + "/" + ym;

      const vector<Dstring>& dir_ymd = get_dir_listing (ym_path, re_ymd);
      for (vector<Dstring>::const_iterator j = dir_ymd.begin ();
           j != dir_ymd.end (); j++)
      {

         const Dstring& ymd = *(j);
         const Dstring& ymd_path = ym_path + "/" + ymd;

         const vector<Dstring>& dir_listing = get_dir_listing (ymd_path, fmt);
         for (vector<Dstring>::const_iterator iterator = dir_listing.begin ();
              iterator != dir_listing.end (); iterator++)
         {

            // fn = filename
            const Dstring& fn = *(iterator);
            const Dstring& file_path = ymd_path + "/" + fn;
            const Dstring& bt_str = fn.substr (6, 8) + fn.substr (15, 2);
            const Dstring& fh_str = fn.substr (20, 3);

            const Dtime base_time (bt_str, Dstring ("%Y%m%d%H"));
            const Integer forecast_hour = stoi (fh_str);

            const Key key (base_time, forecast_hour);
            key_multimap.add (key);
            Grib2* grib_ptr = new Grib2 (file_path);
            insert (make_pair (key, grib_ptr));
            //valid_time_set.insert (valid_time);
            //key_set.insert (Key (base_time, fh));

         }

      }

   }

   if (size () > 0)
   {

      set<uint32_t> set_p;
      const Grib2& grib = *(begin ()->second);
      const Header_Ptr_Map& header_ptr_map = grib.get_header_ptr_map ();

      for (Header_Ptr_Map::const_iterator iterator = header_ptr_map.begin ();
           iterator != header_ptr_map.end (); iterator++)
      {

         const Grib2::Header& header = *(iterator->second);
         const Grib2::Block_4& block_4 = header.block_4;
         const Grib2::Block_4::Level& level = block_4.get_level ();

         if (level.get_first_level_type () == 100)
         {
            const uint32_t p = level.get_first_level ();
            if (p < 200) { continue; }
            set_p.insert (p);
         }

      }

      for (set<uint32_t>::const_iterator iterator = set_p.begin ();
           iterator != set_p.end (); iterator++)
      {
         const Real p = Real (*(iterator));
         tuple_p.push_back (p);
      }

   }

   status = "";
   const set<Dtime>& base_time_set = key_multimap.get_base_time_set ();
   for (set<Dtime>::const_iterator iterator = base_time_set.begin ();
        iterator != base_time_set.end (); iterator++)
   {
      const Dtime& base_time = *(iterator);
      status += base_time.get_string () + " ";
   }

}

void
Gfs::clean_up ()
{
   for (Gfs::iterator iterator = begin (); iterator != end (); iterator++)
   {
      Grib2* grib_ptr = iterator->second;
      delete grib_ptr;
   }
}

void
Gfs::set_domain_2d (const Domain_2D& domain_2d)
{

   this->domain_2d.domain_x.start = ceil (domain_2d.domain_x.start);
   this->domain_2d.domain_x.end = floor (domain_2d.domain_x.end);
   this->domain_2d.domain_y.start = ceil (domain_2d.domain_y.start);
   this->domain_2d.domain_y.end = floor (domain_2d.domain_y.end);

   const Real& start_latitude = this->domain_2d.domain_x.start;
   const Real& end_latitude = this->domain_2d.domain_x.end;
   const Real& start_longitude = this->domain_2d.domain_y.start;
   const Real& end_longitude = this->domain_2d.domain_y.end;
 
   size_2d.i = Integer (round (end_latitude - start_latitude));
   size_2d.j = Integer (round (end_longitude - start_longitude));

}

vector<Dtime>
Gfs::get_valid_time_vector () const
{

   vector<Dtime> valid_time_vector;

   for (Gfs::const_iterator iterator = begin ();
       iterator != end (); iterator++)
   {
      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime dtime (bt.t + fh);
      valid_time_vector.push_back (dtime);
   }

   return valid_time_vector;

/*
   vector<Dtime> valid_time_vector;
   Dtime start_time ("2012011300");
   Dtime end_time ("2012011400");
   for (Real t = start_time.t; t < end_time.t; t += 3)
   {
      Dtime dtime (t);
      valid_time_vector.push_back (dtime);
   }
   return valid_time_vector;
*/
}

Nwp::Key
Gfs::get_key (const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());

   const Dtime& base_time = base_time_vector.back ();
   const Integer forecast_hour = Integer (round (dtime.t - base_time.t));

   return Key (base_time, forecast_hour);

}

void
Gfs::acquire_base_time_forecast_hour (Dtime& base_time,
                                      Integer& forecast_hour,
                                      const Dtime& dtime) const
{

   vector<Dtime> base_time_vector;

   for (Gfs::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Key& key = iterator->first;
      const Dtime& bt = key.base_time;
      const Integer fh = key.forecast_hour;
      const Dtime t (bt.t + fh);

      if (fabs (t.t - dtime.t) < 0.5)
      {
         base_time_vector.push_back (bt);
      }

   }

   if (base_time_vector.size () == 0)
   {
      throw Nwp_Exception ("timestep not available");
   }

   sort (base_time_vector.begin (), base_time_vector.end ());
   base_time.t = base_time_vector.back ().t;
   forecast_hour = Integer (round (dtime.t - base_time.t));

}

