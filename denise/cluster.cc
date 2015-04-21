//
// cluster.cc
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

#include "cluster.h"

using namespace std;
using namespace denise;

Cluster::Node::Node ()
   : index_a (-1),
     index_b (-1),
     distance (GSL_NAN)
{
}

Cluster::Node::Node (const Integer index_a,
                     const Integer index_b,
                     const Real distance)
   : index_a (index_a),
     index_b (index_b),
     distance (distance)
{
}

void
Cluster::construct_distance_matrix ()
{

   const Integer n = dataset.get_number_of_points ();
   distance_matrix_ptr = new Matrix (n, n, true);

   gsl_matrix* gm = distance_matrix_ptr->get_gm ();

   for (Integer i = 0; i < n; i++)
   {
      for (Integer j = i; j < n; j++)
      {
         if (i == j) { continue; }
         const Real d = leaf_distance (i, j);
         gsl_matrix_set (gm, i, j, d);
         gsl_matrix_set (gm, j, i, d);
      }
   }

}

void
Cluster::construct_nodes ()
{

   const Integer n = dataset.get_number_of_points ();
   const Integer number_of_nodes = 2*n - 1;
   node_data = new Node[number_of_nodes];

   set<Integer> cluster_index_set;
   for (Integer i = 0; i < n; i++)
   {
      cluster_index_set.insert (i);
   }

   // For n leaves, there are n - 1 groupings to be done
   for (Integer c = n; c < number_of_nodes; c++)
   {

      set<Integer>::iterator i;
      set<Integer>::iterator a, b;
      set<Integer>::iterator closest_a, closest_b;
      Integer closest_aa, closest_bb;
      Real min_distance = GSL_POSINF;

      for (a = cluster_index_set.begin (); a != cluster_index_set.end (); a++)
      {

         const Integer aa = *(a);

         for (b = a; b != cluster_index_set.end (); b++)
         {

            if (a == b) { continue; }
            const Integer bb = *(b);
            const Real d = get_distance_cluster (aa, bb);

            if (d < min_distance)
            {
               min_distance = d;
               closest_a = a;
               closest_b = b;
               closest_aa = aa;
               closest_bb = bb;
            }

         }
      }


      Node& new_node = node_data[c];
      new_node.index_a = *(closest_a);
      new_node.index_b = *(closest_b);
      new_node.distance = min_distance;

      cluster_index_set.erase (closest_aa);
      cluster_index_set.erase (closest_bb);
      cluster_index_set.insert (c);

   }

}

Real
Cluster::squared_euclidean_distance (const Matrix& data_matrix,
                                     const Integer index_a,
                                     const Integer index_b) const
{

   const Integer k = data_matrix.get_columns ();
   gsl_matrix* gm = data_matrix.get_gm ();

   Real sigma_d2 = 0;

   for (Integer i = 0; i < k; i++)
   {

      const Real a = gsl_matrix_get (gm, index_a, i);
      const Real b = gsl_matrix_get (gm, index_b, i);
      const Real d = (a - b);

      sigma_d2 += d * d;

   }

   return sigma_d2;

}

Real
Cluster::angle_distance (const Matrix& data_matrix,
                         const Integer index_a,
                         const Integer index_b) const
{

   const Integer k = data_matrix.get_columns ();
   gsl_matrix* gm = data_matrix.get_gm ();

   Real a_dot_b = 0;
   Real magnitude_a = 0;
   Real magnitude_b = 0;

   for (Integer i = 0; i < k; i++)
   {

      const Real a = gsl_matrix_get (gm, index_a, i);
      const Real b = gsl_matrix_get (gm, index_b, i);
      const Real d = fabs (a - b);

      a_dot_b += a * b;
      magnitude_a += a * a;
      magnitude_b += b * b;

   }

   return acos (a_dot_b / sqrt (magnitude_a * magnitude_b));

}

Real
Cluster::manhattan_distance (const Matrix& data_matrix,
                             const Integer index_a,
                             const Integer index_b) const
{

   const Integer k = data_matrix.get_columns ();
   gsl_matrix* gm = data_matrix.get_gm ();

   Real distance = 0;

   for (Integer i = 0; i < k; i++)
   {

      const Real a = gsl_matrix_get (gm, index_a, i);
      const Real b = gsl_matrix_get (gm, index_b, i);
      const Real d = fabs (a - b);

      distance += d * d;

   }

   return distance;

}

Real
Cluster::distance (const Matrix& data_matrix,
                   const Integer index_a,
                   const Integer index_b) const
{

   const Matrix& dm = data_matrix;

   switch (distance_genre)
   {

      default:
      case EUCLIDEAN:
      {
         const Real d = squared_euclidean_distance (dm, index_a, index_b);
         return sqrt (d);
      }

      case SQUARED_EUCLIDEAN:
      {
         const Real d = squared_euclidean_distance (dm, index_a, index_b);
         return d;
      }

      case ANGLE:
      {
         const Real d = angle_distance (dm, index_a, index_b);
         return d;
      }

      case MANHATTAN:
      {
         const Real d = manhattan_distance (dm, index_a, index_b);
         return d;
      }

   }

}

Real
Cluster::leaf_distance (const Integer leaf_a, 
                        const Integer leaf_b) const
{
   return distance (dataset.get_data_matrix (), leaf_a, leaf_b);
}

Real
Cluster::cluster_distance_min (const Integer cluster_a,
                               const Integer cluster_b) const
{

   set<Integer> leaf_set_a;
   set<Integer> leaf_set_b;
   fill_leaf_set (leaf_set_a, cluster_a);
   fill_leaf_set (leaf_set_b, cluster_b);

   gsl_matrix* gm = distance_matrix_ptr->get_gm ();
   typedef set<Integer>::const_iterator Iterator;

   Real min_d = GSL_POSINF;

   for (Iterator a = leaf_set_a.begin (); a != leaf_set_a.end (); a++)
   {
      const Integer aa = *(a);
      for (Iterator b = leaf_set_b.begin (); b != leaf_set_b.end (); b++)
      {
         const Integer bb = *(b);
         const Real d = gsl_matrix_get (gm, aa, bb);
         if (d < min_d) { min_d = d; }
      }
   }

   return min_d;

}

Real
Cluster::cluster_distance_max (const Integer cluster_a,
                               const Integer cluster_b) const
{

   set<Integer> leaf_set_a;
   set<Integer> leaf_set_b;
   fill_leaf_set (leaf_set_a, cluster_a);
   fill_leaf_set (leaf_set_b, cluster_b);

   gsl_matrix* gm = distance_matrix_ptr->get_gm ();
   typedef set<Integer>::const_iterator Iterator;

   Real max_d = GSL_NEGINF;

   for (Iterator a = leaf_set_a.begin (); a != leaf_set_a.end (); a++)
   {
      const Integer aa = *(a);
      for (Iterator b = leaf_set_b.begin (); b != leaf_set_b.end (); b++)
      {
         const Integer bb = *(b);
         const Real d = gsl_matrix_get (gm, aa, bb);
         if (d > max_d) { max_d = d; }
      }
   }

   return max_d;

}

Real
Cluster::cluster_distance_mean (const Integer cluster_a,
                                const Integer cluster_b) const
{

   set<Integer> leaf_set_a;
   set<Integer> leaf_set_b;
   fill_leaf_set (leaf_set_a, cluster_a);
   fill_leaf_set (leaf_set_b, cluster_b);

   gsl_matrix* gm = distance_matrix_ptr->get_gm ();
   typedef set<Integer>::const_iterator Iterator;

   Real sigma_d = 0;

   for (Iterator a = leaf_set_a.begin (); a != leaf_set_a.end (); a++)
   {
      const Integer aa = *(a);
      for (Iterator b = leaf_set_b.begin (); b != leaf_set_b.end (); b++)
      {
         const Integer bb = *(b);
         const Real d = gsl_matrix_get (gm, aa, bb);
         sigma_d += d;
      }
   }

   return sigma_d / (leaf_set_a.size () * leaf_set_b.size ());

}

Real
Cluster::cluster_distance_centroid (const Integer cluster_a,
                                    const Integer cluster_b) const
{

   set<Integer> leaf_set_a;
   set<Integer> leaf_set_b;
   fill_leaf_set (leaf_set_a, cluster_a);
   fill_leaf_set (leaf_set_b, cluster_b);

   const Integer k = dataset.get_dimension ();
   Matrix centroid_matrix (2, k, true);

   const Matrix& data_matrix = dataset.get_data_matrix ();
   gsl_matrix* dm = data_matrix.get_gm ();
   gsl_matrix* cm = centroid_matrix.get_gm ();

   typedef set<Integer>::const_iterator Iterator;

   for (Integer i = 0; i < k; i++)
   {
      // calcuate_centroid

      double* a_ptr = gsl_matrix_ptr (cm, 0, i);
      double* b_ptr = gsl_matrix_ptr (cm, 1, i);

      for (Iterator a = leaf_set_a.begin (); a != leaf_set_a.end (); a++)
      {
         const Integer aa = *(a);
         *a_ptr += gsl_matrix_get (dm, aa, i);
      }

      for (Iterator b = leaf_set_b.begin (); b != leaf_set_b.end (); b++)
      {
         const Integer bb = *(b);
         *b_ptr += gsl_matrix_get (dm, bb, i);
      }

      *a_ptr /= leaf_set_a.size ();
      *b_ptr /= leaf_set_b.size ();

   }

   return distance (centroid_matrix, 0, 1);

}

Real
Cluster::get_distance_cluster (const Integer cluster_a,
                               const Integer cluster_b) const
{

   switch (method)
   {

      case MIN_DISTANCE:
      {
         const Real d = cluster_distance_min (cluster_a, cluster_b);
         return d;
      }

      case MAX_DISTANCE:
      {
         const Real d = cluster_distance_max (cluster_a, cluster_b);
         return d;
      }

      default:
      case MEAN:
      case WEIGHTED_MEAN:
      {
         const Real d = cluster_distance_mean (cluster_a, cluster_b);
         return d;
      }

      case CENTROID:
      case WEIGHTED_CENTROID:
      {
         const Real d = cluster_distance_centroid (cluster_a, cluster_b);
         return d;
      }

   }

}

void
Cluster::fill_cluster_set (set<Integer>& cluster_set,
                           const Integer cluster_index,
                           const Real distance) const
{

   const Integer n = dataset.get_number_of_points ();
   const bool is_leaf = (cluster_index < n);

   if (is_leaf)
   {
      cluster_set.insert (cluster_index);
      return;
   }

   const Node& node = node_data[cluster_index];
   const bool node_distance_small_enough = (node.distance < distance);
   if (node_distance_small_enough)
   {
      cluster_set.insert (cluster_index);
      return;
   }

   fill_cluster_set (cluster_set, node.index_a, distance);
   fill_cluster_set (cluster_set, node.index_b, distance);

   return;

}

void
Cluster::fill_leaf_set (set<Integer>& leaf_set,
                        const Integer cluster_index) const
{

   const Integer n = dataset.get_number_of_points ();

   if (cluster_index < n)
   {
      // is a leaf
      leaf_set.insert (cluster_index);
      return;
   }
   else
   {
      // is a node
      const Node& node = node_data[cluster_index];
      fill_leaf_set (leaf_set, node.index_a);
      fill_leaf_set (leaf_set, node.index_b);
   }

}

Cluster::Cluster (const Dataset& dataset,
                  const Method method,
                  const Distance_Genre distance_genre)
   : dataset (dataset),
     method (method),
     distance_genre (distance_genre)
{
   construct_distance_matrix ();
   construct_nodes ();
}

Cluster::~Cluster ()
{
   delete distance_matrix_ptr;
   delete[] node_data;
}

void
Cluster::Multimap::insert (const Integer cluster_index,
                           const Integer leaf_index)
{
   cluster_set.insert (cluster_index);
   multimap<Integer, Integer>::insert (make_pair (cluster_index, leaf_index));
}

Integer
Cluster::Multimap::get_number_of_clusters () const
{
   return cluster_set.size ();
}

const set<Integer>&
Cluster::Multimap::get_cluster_set () const
{
   return cluster_set;
}

set<Integer>
Cluster::Multimap::get_leaf_index_set (const Integer cluster_index) const
{

   typedef Cluster::Multimap::const_iterator Iterator;
   typedef pair<Iterator, Iterator> Range;
   set<Integer> leaf_index_set;

   const Range& range = equal_range (cluster_index);

   for (Iterator j = range.first; j != range.second; j++)
   {
      const Integer leaf_index = j->second;
      leaf_index_set.insert (leaf_index);
   }

   return leaf_index_set;

}

Cluster::Multimap
Cluster::get_multimap (const Real distance) const
{

   set<Integer> cluster_set;
   Multimap multimap;

   const Integer n = dataset.get_number_of_points ();
   const Integer top_cluster_index = n + n - 2;
   fill_cluster_set (cluster_set, top_cluster_index, distance);

   for (set<Integer>::const_iterator i = cluster_set.begin ();
        i != cluster_set.end (); i++)
   {

      const Integer cluster_index = *(i);

      set<Integer> leaf_set;
      fill_leaf_set (leaf_set, cluster_index);

      for (set<Integer>::const_iterator j = leaf_set.begin ();
           j != leaf_set.end (); j++)
      {
         const Integer leaf_index = *(j);
         multimap.insert (cluster_index, leaf_index);
      }

   }

   return multimap;

}

/*
*/
