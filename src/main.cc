#include <cstdlib>
#include <iostream>
#include <readline/readline.h>
#include <readline/history.h>
#include "denise.h"

using namespace std;
using namespace denise;

int
main (int argc,
      char** argv)
{

   bool done;
   char* line;

   Denise denise;
   //char prompt[] = "";
   char prompt[] = "\ndenise> ";

   for ( ; !done; )
   {

      line = readline (prompt);
      if (!line) { break; }

      try
      {

         const Tokens tokens (line);
         if (tokens.size () == 0) { continue; }

         denise.parse (tokens);

      }
      catch (const Exception& e)
      {
         cerr << "denise: " << e << " " << line << endl;
      }
      catch (const std::out_of_range& oor)
      {
         cerr << "denise: out_of_range " << oor.what () << " " << line << endl;
      }

      add_history (line);
      free (line);

   }

}

