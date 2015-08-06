#include <map>
#include <denise/geodesy.h>
#include <denise/gis.h>
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

class Geodesy_Package : protected Andrea_Package
{

   protected:

      Integer
      lat_long_dp;

      map<string, Geodesy>
      geodesy_map;

      Geodesy_Package (const Andrea& andrea);

      void
      geodesy_assign (const string& identifier,
                      const string& str);

      void
      geodesy_print (const string& identifier) const;

      void
      geodesy_distance (const Tokens& tokens) const;

      void
      geodesy_azimuth (const Tokens& tokens) const;

      void
      geodesy_destination (const Tokens& tokens) const;

      void
      geodesy_parse (const Tokens& tokens);

   public: 

      const map<string, Geodesy>&
      get_geodesy_map () const;

      const Geodesy&
      get_geodesy (const string& idenfifier) const;

};

class Journey_Package : protected Andrea_Package
{

   protected:

      map<string, Journey>
      journey_map;

      Journey_Package (const Andrea& andrea);

      void
      journey_assign (const string& identifier,
                      const Tokens& arguments);

      void
      journey_print (const string& identifier) const;

      void
      journey_parse (const Tokens& tokens);

   public: 

      const map<string, Journey>&
      get_journey_map () const;

      const Journey&
      get_journey (const string& identifier) const;

};

class Geodetic_Mesh_Package : public Andrea_Package
{

   protected:

      map<string, Geodetic_Mesh>
      geodetic_mesh_map;

      Geodetic_Mesh_Package (const Andrea& andrea);

      void
      geodetic_mesh_assign (const string& identifier,
                            const Size_2D& size_2d,
                            const Domain_2D& domain_2d);

      void
      geodetic_mesh_add (const string& identifier,
                         const Tokens& arguments);

      void
      geodetic_mesh_print (const string& identifier,
                           const Tokens& arguments) const;

      void
      geodetic_mesh_parse (const Tokens& tokens);

   public:

      const map<string, Geodetic_Mesh>&
      get_geodetic_mesh_map () const;

      const Geodetic_Mesh&
      get_geodetic_mesh (const string& identifier) const;

};

class Geodetic_Transform_Package : public Andrea_Package
{

   protected:

      map<string, string>
      geodetic_transform_str_map;

      Geodetic_Transform_Package (const Andrea& andrea);

      void
      geodetic_transform_assign (const string& identifier,
                                 const string& str);

      void
      geodetic_transform_print (const string& identifier,
                                const Tokens& arguments) const;

      void
      geodetic_transform_parse (const Tokens& tokens);

   public:

      const map<string, string>&
      get_geodetic_transform_str_map () const;

      const Geodetic_Transform*
      get_geodetic_transform_ptr (const string& identifier,
                                  const Point_2D& point) const;

};

class Gshhs_Package : public Andrea_Package
{

   protected:

      map<string, Gshhs*>
      gshhs_ptr_map;

      Gshhs_Package (const Andrea& andrea);

      ~Gshhs_Package ();

      void
      gshhs_load (const string& identifier,
                     const string& file_path);

      void
      gshhs_print (const string& identifier,
                   const Tokens& arguments) const;

      void
      gshhs_parse (const Tokens& tokens);

   public:

      const map<string, Gshhs*>&
      get_gshhs_ptr_map () const;

      const Gshhs*
      get_gshhs_ptr (const string& identifier) const;

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
      sounding_load (const string& identifier,
                     const string& file_path);

      void
      sounding_print (const string& identifier,
                      const Tokens& arguments) const;

      void
      sounding_parse (const Tokens& tokens);

   public:

      const map<string, Sounding>&
      get_sounding_map () const;

      const Sounding&
      get_sounding (const string& identifier) const;

};

class Image_Package : public Andrea_Package
{

   protected:

      map<string, RefPtr<ImageSurface> >
      image_map;

      map<string, RefPtr<Context> >
      cr_map;

      Image_Package (const Andrea& andrea);

      void
      image_init (const string& identifier,
                  const string& geometry);

      void
      image_paint (const string& identifier,
                   const Color& color = Color::white ()) const;

      void
      image_save (const string& identifier,
                  const string& file_path) const;

      void
      image_title (const string& identifier,
                   const Tokens& tokens) const;

      void
      image_sounding (const Tokens& tokens) const;

      void
      image_sounding_tephigram (const Tokens& tokens) const;

      void
      image_sounding_chart (const Tokens& tokens) const;

      void
      image_sounding_chart (const RefPtr<Context>& cr,
                            const Transform_2D& transform,
                            const bool is_p,
                            const Mesh_2D& mesh_2d,
                            const string& fmt_x,
                            const string& fmt_y,
                            const string& unit_x,
                            const string& unit_y,
                            const Sounding& sounding,
                            const Real_Profile& real_profile,
                            const Symbol& symbol,
                            const Color& color) const;

      void
      image_journey (const string& image_identifier,
                     const string& geodetic_transform_identifier,
                     const string& journey_identifier);

      void
      image_geodetic_mesh (const string& image_identifier,
                           const string& geodetic_transform_identifier,
                           const string& geodetic_mesh_identifier);

      void
      image_gshhs (const string& image_identifier,
                   const string& geodetic_transform_identifier,
                   const string& gshhs_identifier,
                   const Tokens& arguments);

      void
      image_parse (const Tokens& tokens);

   public:

      const map<string, RefPtr<ImageSurface> >&
      get_image_map () const;

      const map<string, RefPtr<Context> >&
      get_cr_map () const;

      const RefPtr<ImageSurface>&
      get_image (const string& identifier) const;

      const RefPtr<Context>&
      get_cr (const string& identifier) const;

};

class Entity : public string
{

   public:

      Entity (const string& str);

      Real
      value () const;

};

class Andrea : public Geodesy_Package,
               public Gshhs_Package,
               public Sounding_Package,
               public Image_Package,
               public Journey_Package,
               public Geodetic_Mesh_Package,
               public Geodetic_Transform_Package
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

