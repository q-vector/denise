#include <map>
#include <denise/geodesy.h>
#include <denise/thermo.h>

using namespace std;
using namespace denise;

class Andrea;

class Andrea_Package
{

   protected:

      const Andrea&
      andrea;

      Andrea_Package (const Andrea& andrea);

};

class Journey_Package : protected Andrea_Package
{

   protected:

      Integer
      lat_long_dp;

      map<string, Journey>
      journey_map;

      Journey_Package (const Andrea& andrea);

      void
      assign_journey (const string& variable,
                      const Tokens& arguments);

      void
      print_journey (const string& variable) const;

      void
      journey_distance (const Tokens& tokens) const;

      void
      journey_azimuth (const Tokens& tokens) const;

      void
      journey_destination (const Tokens& tokens) const;

      void
      parse_journey (const Tokens& tokens);

   public: 

      const map<string, Journey>&
      get_journey_map () const;

};

class Sounding_Package : public Andrea_Package
{

   protected:

      map<string, Sounding>
      sounding_map;

      Tephigram
      tephigram;

      Sounding_Package (const Andrea& andrea);

      void
      load_sounding (const string& variable,
                     const string& file_path);

      void
      print_sounding (const string& variable,
                      const Tokens& arguments) const;

      void
      parse_sounding (const Tokens& tokens);

   public:

      const map<string, Sounding>&
      get_sounding_map () const;

};

class Image_Package : public Andrea_Package
{

   protected:

      map<string, RefPtr<ImageSurface> >
      image_map;

      Image_Package (const Andrea& andrea);

      void
      init_image (const string& variable,
                  const string& geometry);

      void
      save_image (const string& variable,
                  const string& file_path) const;

      void
      image_title (const string& variable,
                   const Tokens& tokens) const;

      void
      image_tephigram (const string& image_identifier,
                       const string& tephigram_identifier,
                       const Tokens& arguments) const;

      void
      image_sounding (const Tokens& tokens) const;

      void
      image_sounding_tephigram (const Tokens& tokens) const;

      void
      image_sounding_chart (const Tokens& tokens) const;

      void
      parse_image (const Tokens& tokens);

   public:

      const map<string, RefPtr<ImageSurface> >&
      get_image_map () const;

};

class Entity : public string
{

   public:

      Entity (const string& str);

      Real
      value () const;

};

class Andrea : public Journey_Package,
               public Sounding_Package,
               public Image_Package
{

   private:

      void
      wind_shear (const Tokens& arguments) const;

      void
      print (const Entity& entity) const;

   public:

      Andrea ();

      void
      parse (const Tokens& tokens);

};

