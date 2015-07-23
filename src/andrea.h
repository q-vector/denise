#include <map>
#include <denise/geodesy.h>

using namespace std;
using namespace denise;

class Journey_Map : public map<string, Journey>
{
};

class Andrea
{

   private:

      Integer
      lat_long_dp;

      Journey_Map
      journey_map;

      void
      assign_journey (const string& variable,
                      const Tokens& arguments);

      void
      journey (const Tokens& arguments);

      void
      distance (const Tokens& arguments);

      void
      azimuth (const Tokens& arguments);

      void
      destination (const Tokens& arguments);

      void
      wind_shear (const Tokens& arguments);

   public:

      Andrea ();

      void
      parse (const Tokens& tokens);

};

