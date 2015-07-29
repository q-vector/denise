#include <cstdlib>
#include <iostream>
#include <readline/readline.h>
#include <readline/history.h>
#include "andrea.h"

using namespace std;
using namespace denise;

int
main (int argc,
      char** argv)
{

   bool done;
   char* line;

   Andrea andrea;
   //char prompt[] = "";
   char prompt[] = "\nandrea> ";

   for ( ; !done; )
   {

      line = readline (prompt);
      if (!line) { break; }

      const string& input_line = get_trimmed (string (line));
      if (input_line[0] == '#') { continue; }

      try
      {

         const Tokens tokens (line);
         if (tokens.size () == 0) { continue; }

         andrea.parse (tokens);

      }
      catch (const Exception& e)
      {
         cerr << "andrea: " << e << " " << line << endl;
      }
      catch (const std::out_of_range& oor)
      {
         cerr << "andrea: out_of_range " << oor.what () << " " << line << endl;
      }

      add_history (line);
      free (line);

   }

}

