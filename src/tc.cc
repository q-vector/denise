//
// tc.cc
// 
// Copyright (C) 2005-2015 Simon E. Ching
// 
// This file is part of the denise suite.
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

#include "andrea.h"
#include "tc.h"

using namespace std;
using namespace denise;
using namespace andrea;

Best_Tracks_Package::Best_Tracks_Package (Andrea& andrea)
   : Package (andrea)
{
}

Best_Tracks_Package::~Best_Tracks_Package ()
{
}

void
Best_Tracks_Package::best_tracks_ingest (const Dstring& identifier,
                                         const Dstring& file_path)
{

   map<Dstring, Best_Tracks>::iterator i;

   if ((i = best_tracks_map.find (identifier)) == best_tracks_map.end ())
   {
      i = best_tracks_map.insert (make_pair (identifier, Best_Tracks ())).first;
   }

   Best_Tracks& best_tracks = i->second;
   best_tracks.ingest_jma (file_path);

}

void
Best_Tracks_Package::best_tracks_parse (const Tokens& tokens)
{

   const Integer n = tokens.size ();

   if (tokens[0] == "ingest")
   {
      const Dstring& identifier = tokens[1];
      const Dstring& file_path = tokens[2];
      best_tracks_ingest (identifier, file_path);
   }

}

const map<Dstring, Best_Tracks>&
Best_Tracks_Package::get_best_tracks_map () const
{
   return best_tracks_map;
}

const Best_Tracks&
Best_Tracks_Package::get_best_tracks (const Dstring& identifier) const
{
   Exception e ("best_track not found: " + identifier);
   try { return best_tracks_map.at (identifier); }
   catch (const std::out_of_range& oor) { throw e; }
}

