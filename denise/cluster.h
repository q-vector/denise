//
// cluster.h
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

#ifndef DENISE_CLUSTER_H
#define DENISE_CLUSTER_H

#include <map>
#include <denise/basics.h>
#include <denise/linalg.h>
#include <denise/dataset.h>

using namespace std;

namespace denise
{

   class Cluster
   {

      private:

         enum Distance_Genre
         {
            EUCLIDEAN,
            SQUARED_EUCLIDEAN,
            ANGLE,
            MANHATTAN
         };

         enum Method
         {
            MIN_DISTANCE,
            MAX_DISTANCE,
            MEAN,
            WEIGHTED_MEAN,
            CENTROID,
            WEIGHTED_CENTROID
         };

         class Node
         {

            public:

               Integer
               index_a;

               Integer
               index_b;

               Real
               distance;

               Node ();

               Node (const Integer index_a,
                     const Integer index_b,
                     const Real distance);

         };

         const Dataset&
         dataset;

         Method
         method;

         Distance_Genre
         distance_genre;

         Matrix*
         distance_matrix_ptr;

         Node*
         node_data;

         void
         construct_distance_matrix ();

         void
         construct_nodes ();

         Real
         squared_euclidean_distance (const Matrix& data_matrix,
                                     const Integer index_a,
                                     const Integer index_b) const;

         Real
         angle_distance (const Matrix& data_matrix,
                         const Integer index_a,
                         const Integer index_b) const;

         Real
         manhattan_distance (const Matrix& data_matrix,
                             const Integer index_a,
                             const Integer index_b) const;

         Real
         distance (const Matrix& data_matrix,
                   const Integer index_a,
                   const Integer index_b) const;

         Real
         leaf_distance (const Integer leaf_a,
                        const Integer leaf_b) const;

         Real
         cluster_distance_min (const Integer cluster_a,
                               const Integer cluster_b) const;

         Real
         cluster_distance_max (const Integer cluster_a,
                               const Integer cluster_b) const;

         Real
         cluster_distance_mean (const Integer cluster_a,
                                const Integer cluster_b) const;

         Real
         cluster_distance_centroid (const Integer cluster_a,
                                    const Integer cluster_b) const;

         Real
         get_distance_cluster (const Integer cluster_a,
                               const Integer cluster_b) const;

         void
         fill_cluster_set (set<Integer>& leaf_set,
                           const Integer cluster_index,
                           const Real distance) const;

         void
         fill_leaf_set (set<Integer>& leaf_set,
                        const Integer cluster_index) const;

      public:

         Cluster (const Dataset& dataset,
                  const Method method = MEAN,
                  const Distance_Genre distance_genre = EUCLIDEAN);

         ~Cluster ();

         class Multimap : public multimap<Integer, Integer>
         {

            private:

               set<Integer>
               cluster_set;

            public:

               void
               insert (const Integer cluster_index,
                       const Integer leaf_index);

               Integer
               get_number_of_clusters () const;

               const set<Integer>&
               get_cluster_set () const;

               set<Integer>
               get_leaf_index_set (const Integer cluster_index) const;

         };

         Multimap
         get_multimap (const Real distance) const;

   };

/*
*/

}

#endif /* DENISE_CLUSTER_H */
