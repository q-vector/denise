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

Grib::Grib (const string& file_path)
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

Grib2::Grib2 (const string& file_path)
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
const string b0 (bit_0 ? "yes" : "no");
const string b1 (bit_1 ? "yes" : "no");
const string b2 (bit_2 ? "yes" : "no");
const string b3 (bit_3 ? "yes" : "no");
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
   operator << (ostream &out_file,
                const Grib::Key& key)
   {
      for (Integer i = 0; i < 16; i++)
      {
         out_file << Integer (key.buffer[i]) << " ";
      }
      return out_file;
   }

};


namespace denise
{

   ostream&
   operator << (ostream &out_file,
                const Grib2::Key& key)
   {
      for (Integer i = 0; i < 26; i++)
      {
         string sep (" ");
         if (i == 6 || i == 11 || i == 13) { sep = " : "; }
         out_file << Integer (key.buffer[i]) << sep;
      }
      return out_file;
   }

};


