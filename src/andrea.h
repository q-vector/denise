#include <map>
#include <denise/geodesy.h>
#include <denise/thermo.h>

using namespace std;
using namespace denise;

class Journey_Map : public map<string, Journey>
{

   private:

      Integer
      lat_long_dp;

      void
      assign (const string& variable,
              const Tokens& arguments);

      void
      print (const string& variable) const;

      void
      distance (const Tokens& tokens) const;

      void
      azimuth (const Tokens& tokens) const;

      void
      destination (const Tokens& tokens) const;

   public:

      Journey_Map ();

      void
      parse (const Tokens& tokens);

};

class Sounding_Map : public map<string, Sounding>
{

   private:

      Tephigram
      tephigram;

      void
      load (const string& variable,
            const string& file_path);

      void
      print (const string& variable,
             const Tokens& arguments) const;

   public:

      Sounding_Map ();

      void
      parse (const Tokens& tokens);

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

      Journey_Map
      journey_map;

      Sounding_Map
      sounding_map;

      void
      wind_shear (const Tokens& arguments) const;

      void
      print (const Entity& entity) const;

   public:

      Andrea ();

      void
      parse (const Tokens& tokens);

};

