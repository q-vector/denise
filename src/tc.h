//
// tc.h
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

#ifndef DENISE_SRC_TC_H
#define DENISE_SRC_TC_H

#include <map>
#include <denise/tc.h>
#include "package.h"

using namespace std;
using namespace denise;

namespace andrea
{

   class Best_Tracks_Package : public Package
   {

      protected:

         map<Dstring, Best_Tracks>
         best_tracks_map;

         Best_Tracks_Package (Andrea& andrea);

         ~Best_Tracks_Package ();

         void
         best_tracks_ingest (const Dstring& identifier,
                             const Dstring& file_path);

         void
         best_tracks_parse (const Tokens& tokens);

      public:

         const map<Dstring, Best_Tracks>&
         get_best_tracks_map () const;

         const Best_Tracks&
         get_best_tracks (const Dstring& identifier) const;

   };

};

#endif /* DENISE_TC_GSHHS_H */
