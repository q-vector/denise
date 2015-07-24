#include <map>
#include <denise/geodesy.h>

using namespace std;
using namespace denise;

class Journey_Map : public map<string, Journey>
{
};

class Entity : public string
{

   public:

      Entity (const string& str);

      Real
      value () const;

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
      distance (const Tokens& arguments) const;

      void
      azimuth (const Tokens& arguments) const;

      void
      destination (const Tokens& arguments) const;

      void
      wind_shear (const Tokens& arguments) const;

      void
      print (const Entity& entity) const;

   public:

      Andrea ();

      void
      parse (const Tokens& tokens);

};

