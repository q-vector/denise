//
// ds.h
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

#ifndef DENISE_DS_H
#define DENISE_DS_H

#include <denise/basics.h>
#include <denise/exception.h>

using namespace std;

namespace denise
{

   template <class T> class Heap
   {

      private:

         bool
         delete_on_destruction;

         T*
         array;

         Integer
         n;

         const Integer
         max_size;

         Heap (T* array,
               const Integer n)
          : n (n),
            max_size (n),
            array (array),
            delete_on_destruction (false)
         {
            for (Integer i = n/2 - 1; i >= 0; i--) { heapify (i, n); }
         }

         void
         heapify (const Integer i,
                  const Integer heap_size)
         {

            Integer l = 2*i + 1;
            Integer r = 2*i + 2;
            Integer largest;

            if (l < heap_size && array[l] > array[i])
            { largest = l; } else { largest = i; }

            if (r < heap_size && array[r] > array[largest])
            { largest = r; }

            if (largest != i)
            {
               std::swap (array[i], array[largest]);
               heapify (largest, heap_size);
            }

         }

         void
         sort ()
         {
            for (Integer i = n-1; i > 0; i--)
            {
               std::swap (array[0], array[i]);
               heapify (0, i);
            }
         }

      public:

         Heap (const Integer max_size) throw (std::bad_alloc)
          : n (0),
            max_size (max_size),
            delete_on_destruction (true)
         {
            array = new T[max_size];
         }

         ~Heap ()
         {
            if (delete_on_destruction) { delete[] array; }
         }

         const Integer&
         size () const
         {
            return n;
         }

         const T&
         get_max () const
         {
            return array[0];
         }

         void
         insert (const T& datum)
         {
            if (n == max_size)
            {
               throw Exception ("heap overflow");
            }
            else
            {
               n++;
               Integer i = n - 1;
               Integer p = (i - 1) / 2;
               while (i > 0 && array[p] < datum)
               {
                  array[i] = array[p];
                  i = p;
                  p = (i - 1) / 2;
               }
               array[i] = datum;
            }
         }

         T
         pop_max ()
         {
            if (n == 0)
            {
               throw Exception ("pop from empty heap");
            }
            else
            {
               const T datum = array[0];
               array[0] = array[n - 1];
               n--;
               heapify (0, n);
               return datum;
            }
         }

         static void
         sort (T* array,
               const Integer n)
         {
            Heap<T> heap (array, n);
            heap.sort ();
         }

         void
         print () const
         {
            for (Integer i = 0; i < n; i++) { cout << array[i] << " "; }
            cout << endl;
         }

   };

}

#endif /* DENISE_DS_H */

