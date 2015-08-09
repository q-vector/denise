#include "andrea.h"

using namespace std;
using namespace denise;
using namespace andrea;

int
main (int argc,
      char** argv)
{

   //const string prompt ("\nandrea> ");
   const string prompt ("");

   Andrea andrea (prompt);
   andrea.loop ();
}

