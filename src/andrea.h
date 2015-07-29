#include <map>
#include <denise/geodesy.h>
#include <denise/thermo.h>

using namespace std;
using namespace denise;

class Andrea;

class Journey_Map : public map<string, Journey>
{

   private:

      const Andrea&
      andrea;

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

      Journey_Map (const Andrea& andrea);

      void
      parse (const Tokens& tokens);

};

class Sounding_Map : public map<string, Sounding>
{

   private:

      const Andrea&
      andrea;

      Tephigram
      tephigram;

      void
      load (const string& variable,
            const string& file_path);

      void
      print (const string& variable,
             const Tokens& arguments) const;

   public:

      Sounding_Map (const Andrea& andrea);

      void
      parse (const Tokens& tokens);

};

class Image_Map : public map<string, RefPtr<ImageSurface> >
{

   private:

      const Andrea&
      andrea;

      void
      init (const string& variable,
            const string& geometry);

      void
      save (const string& variable,
            const string& file_path) const;

      void
      title (const string& variable,
             const Tokens& tokens) const;

      void
      tephigram (const string& image_identifier,
                 const string& tephigram_identifier,
                 const Tokens& arguments) const;

      void
      sounding (const Tokens& tokens) const;

      void
      sounding_tephigram (const Tokens& tokens) const;

      void
      sounding_chart (const Tokens& tokens) const;

   public:

      Image_Map (const Andrea& andrea);

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

      Image_Map
      image_map;

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

      const Image_Map&
      get_image_map () const;

      const Journey_Map&
      get_journey_map () const;

      const Sounding_Map&
      get_sounding_map () const;

      void
      parse (const Tokens& tokens);

};

