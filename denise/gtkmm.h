// gtkmm.h
// 
// Copyright (C) 2005-2013 Simon E. Ching
// 
// This file is part of libdenise.
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

#ifndef DENISE_GTKMM_H
#define DENISE_GTKMM_H

#include <map>
#include <set>
#include <glibmm.h>
#include <gtkmm.h>
#include <gtkmm/drawingarea.h>
#include <denise/basics.h>
#include <denise/geometry.h>
#include <denise/gis.h>
#include <denise/graphics.h>
#include <denise/met.h>
#include <denise/nwp.h>

using namespace std;
using namespace Cairo;
using namespace denise;

namespace denise
{

   class Dwidget;
   class Dcanvas;
   class Dcontainer;
   class Map_Console;
   class Geodetic_Transform;

   enum Orientation
   {
      ORIENTATION_NE,
      ORIENTATION_SE,
      ORIENTATION_NW,
      ORIENTATION_SW,
      NUMBER_OF_ORIENTATIONS
   };

   class Devent
   {

      public:

         class Modifier
         {

            public:

               bool
               shift;

               bool
               control;

               bool
               alt;

               Modifier (const bool shift = false,
                         const bool control = false,
                         const bool alt = false);

               Modifier (const Modifier& modifier);

               bool
               modified () const;

         };

         GdkEventType
         type;

         guint
         state;

         Devent (const bool shift = false,
                 const bool control = false,
                 const bool alt = false);

         Devent (const GdkEventType type,
                 const guint state);

         virtual
         ~Devent () {}

         Modifier
         get_modifier () const;

         bool
         shift () const;

         bool
         control () const;

         bool
         alt () const;

   };

   class Dmouse_Event : public Devent
   {

      public:

         const Point_2D
         point;

         Dmouse_Event (const GdkEventType type,
                       const guint state,
                       const Point_2D& point);

         Dmouse_Event (const bool shift = false,
                       const bool control = false,
                       const bool alt = false,
                       const Point_2D& point = Point_2D (-1, -1));

   };

   class Dmouse_Button_Event : public Dmouse_Event
   {

      public:

         const guint
         button;

         Dmouse_Button_Event (const GdkEventType type,
                              const guint state,
                              const guint button,
                              const Point_2D& point);

         Dmouse_Button_Event (const bool shift = false,
                              const bool control = false,
                              const bool alt = false,
                              const Point_2D& point_2d = Point_2D (-1, -1),
                              const guint button = 0);

   };

   class Dmouse_Scroll_Event : public Dmouse_Event
   {

      public:

         GdkScrollDirection
         direction;

         Dmouse_Scroll_Event (const GdkEventType type,
                              const guint state,
                              const GdkScrollDirection direction,
                              const Point_2D& point);

   };

   class Dmouse_Motion_Event : public Dmouse_Event
   {

      public:

         Dmouse_Motion_Event (const GdkEventType type,
                              const guint state,
                              const Point_2D& point);

   };

   class Dmouse_Crossing_Event : public Dmouse_Event
   {

      public:

         Dmouse_Crossing_Event (const GdkEventType type,
                                const guint state,
                                const Point_2D& point);

   };

   class Dkey_Event : public Devent
   {

      public:

         const guint
         value;

         Dkey_Event (const GdkEventType type,
                     const guint state,
                     const guint value);

   };

   class Image_Buffer : public RefPtr<ImageSurface>
   {

      public:

         RefPtr<ImageSurface>
         image_surface;

         bool
         ready;

         Image_Buffer ();

         void
         initialize (const Dwidget& widget);

         const RefPtr<Context>
         get_cr () const;

         void
         blit (const RefPtr<Context>& cr);

         void
         clear ();

   };

   class Dwidget
   {

      private:

         Point_2D
         prev_point;

         bool
         hidable;

      protected:

         Dcanvas&
         canvas;

         Dcontainer*
         container_ptr;

         Image_Buffer
         background_buffer;

         bool
         activated;

         Point_2D
         anchor;

         Real
         font_size;

         Real
         width;

         Real
         height;

         Real
         preferred_width;

         Real
         preferred_height;

      public:

         Dwidget (const Dcanvas& canvas);

         void
         set_container (Dcontainer& container);

         void
         activate ();

         void
         deactivate ();

         virtual void
         set_hidable (const bool hidable);

         bool
         get_hidable () const;

         const Point_2D&
         get_anchor () const;

         Real
         get_width () const;

         Real
         get_height () const;

         Size_2D
         get_size_2d () const;

         Real
         get_font_size () const;

         void
         set_font_size (const Real font_size);

         virtual void
         being_packed (const Point_2D& anchor,
                       const Real width,
                       const Real height);

         virtual void
         set_preferred_size (const Real preferred_width,
                             const Real preferred_height);

         virtual Size_2D
         get_preferred_size_2d () const;

         virtual Real
         get_preferred_width () const;

         virtual Real
         get_preferred_height () const;

         virtual bool
         out_of_bounds (const Point_2D& point) const;

         virtual bool
         on_key_pressed (const Dkey_Event& event);

         virtual bool
         on_key_released (const Dkey_Event& event);

         virtual bool
         on_mouse_button_pressed (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_button_released (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_scroll (const Dmouse_Scroll_Event& event);

         virtual bool
         on_mouse_motion (const Dmouse_Motion_Event& event);

         virtual bool
         on_mouse_enter (const Dmouse_Crossing_Event& event);

         virtual bool
         on_mouse_exit (const Dmouse_Crossing_Event& event);

         virtual void
         highlight (const RefPtr<Context>& cr,
                    const Color& color,
                    const Real line_width = 5) const;

         void
         set_background_ready (const bool background_ready);

         void
         blit_background_buffer (const RefPtr<Context>& cr);

         virtual void
         render_background_buffer ();

         virtual void
         render_background_buffer (const RefPtr<Context>& cr);

         virtual void
         cairo (const RefPtr<Context>& cr);

   };

   class Dcontainer : public Dwidget
   {

      protected:

         map<Integer, Dwidget*>
         widget_ptr_map;

         bool
         packed;

         bool
         hide_hidable;

      public:

         Dcontainer (const Dcanvas& canvas);

         void
         register_widget (Dwidget& widget);

         void
         clear ();

         virtual void
         pack ();

         virtual void
         set_hidable (const bool hidable);

         virtual void
         set_hide_hidable (const bool hide_hidable);

         bool
         get_hide_hidable () const;

         virtual bool
         on_mouse_button_pressed (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_button_released (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_scroll (const Dmouse_Scroll_Event& event);

         virtual bool
         on_mouse_motion (const Dmouse_Motion_Event& event);

         virtual void
         render_widgets (const RefPtr<Context>& cr);

         virtual void
         cairo (const RefPtr<Context>& cr);

         virtual void
         refresh (const RefPtr<Context>& cr,
                  const Point_2D& point,
                  const Real width,
                  const Real height);

   };

   class Dpack_Box : public Dcontainer
   {

      protected:

         const bool
         horizontal;

         const Real
         margin;

         const Real
         spacing;

         const Real
         separator;

         set<Integer>
         separator_set;

         Real
         get_widget_width_h () const;

         Real
         get_widget_width_v () const;

         Real
         get_widget_height_h () const;

         Real
         get_widget_height_v () const;

         Real
         get_preferred_width_h () const;

         Real
         get_preferred_width_v () const;

         Real
         get_preferred_height_h () const;

         Real
         get_preferred_height_v () const;

         void
         pack_h ();

         void
         pack_v ();

      public:

         Dpack_Box (const Dcanvas& canvas,
                    const Real margin,
                    const Real spacing,
                    const bool horizontal = true);

         void
         pack_front (const Dwidget& widget);

         void
         pack_back (const Dwidget& widget);

         void
         add_separator (const bool front = false);

         virtual Real
         get_widget_width () const;

         virtual Real
         get_widget_height () const;

         virtual Real
         get_preferred_width () const;

         virtual Real
         get_preferred_height () const;

         virtual void
         pack ();

   };

   class Dh_Pack_Box : public Dpack_Box
   {

      public:

         Dh_Pack_Box (const Dcanvas& canvas,
                      const Real margin,
                      const Real spacing);

   };

   class Dv_Pack_Box : public Dpack_Box
   {

      public:

         Dv_Pack_Box (const Dcanvas& canvas,
                      const Real margin,
                      const Real spacing);

   };

   class Dgrid_Box : public Dcontainer
   {

      protected:

         Size_2D
         size_2d;

         const Real
         margin;

         const Real
         spacing;

         map<Index_2D, Integer>
         widget_id_map;

      public:

         Dgrid_Box (const Dcanvas& canvas,
                    const Real margin,
                    const Real spacing,
                    const Size_2D& size_2d = Size_2D (0, 0));

         void
         pack (const Dwidget& widget,
               const Index_2D& index_2d);

         virtual Real
         get_widget_width () const;

         virtual Real
         get_widget_height () const;

         virtual Real
         get_preferred_width () const;

         virtual Real
         get_preferred_height () const;

         virtual void
         pack ();

   };

   class Widget_Panel : public Dh_Pack_Box
   {

      private:

         const Real
         hide_button_size;

      protected:

         virtual Real
         get_widget_width () const;

      public:

         Widget_Panel (const Dcanvas& canvas,
                       const Real margin,
                       const Real spacing,
                       const Real hide_button_size);


         virtual void
         cairo (const RefPtr<Context>& cr);

   };

   class Progress
   {

      public:

         Real
         fraction;

         string
         description;

         Progress (const Real fraction,
                   const string& description);

   };

   class With_Progress
   {

      public:

         virtual void
         update_progress (const Progress& progress) = 0;

   };

   class Dcanvas : public Dcontainer, public Gtk::DrawingArea
   {

      protected:

         Glib::Mutex
         mutex;

         class Title : public denise::Title
         {

            protected:

               Dcanvas&
               dcanvas;

            public:

               Title (Dcanvas& dcanvas,
                      const Color& bg_color = Color (0, 0, 0, 0.5),
                      const Color& fg_color = Color (1, 1, 1),
                      const Color& shadow_color = Color (0, 0, 0, 0.25));

               void
               cairo (const RefPtr<Context>& cr);

         };

         Gtk::Window&
         gtk_window;

         string
         easter_egg;

         Title
         title;

         RefPtr<ImageSurface>
         widget_surface;

         RefPtr<ImageSurface>
         image_surface;

         Image_Buffer
         image_buffer;

         Image_Buffer
         foreground_buffer;

         virtual bool
         on_configure_event (GdkEventConfigure* event);

         virtual bool
         on_draw (const RefPtr<Context>& cr);

         virtual void
         cairo (const RefPtr<Context>& cr);

      public:

         Dcanvas (Gtk::Window& gtk_window);

         const RefPtr<Context>
         get_widget_cr () const;

         const RefPtr<Context>
         get_image_surface_cr () const;

         bool
         is_initialized () const;

         void
         initialize ();

         virtual void
         save_image (const string& file_path);

         virtual bool
         on_key_pressed (const Dkey_Event& event);

         virtual bool
         on_key_press_event (GdkEventKey* event);

         virtual bool
         on_key_release_event (GdkEventKey* event);

         virtual bool
         on_button_press_event (GdkEventButton* event);

         virtual bool
         on_button_release_event (GdkEventButton* event);

         virtual bool
         on_scroll_event (GdkEventScroll* event);

         virtual bool
         on_motion_notify_event (GdkEventMotion* event);

         virtual bool
         on_enter_notify (GdkEventCrossing* event);

         virtual bool
         on_leave_notify (GdkEventCrossing* event);

         virtual bool
         on_delete_event (GdkEventAny* event);

         virtual void
         set_preferred_size (const Real preferred_width,
                             const Real preferred_height);

         virtual void
         refresh (const RefPtr<Context>& cr);

         virtual void
         refresh (const RefPtr<Context>& cr,
                  const Point_2D& point,
                  const Real width,
                  const Real height);

         virtual void
         render ();

         virtual void
         render_queue_draw ();

         void
         set_image_ready (const bool ready);

         void
         blit_image_buffer (const RefPtr<Context>& cr);

         virtual void
         render_image_buffer ();

         virtual void
         render_image_buffer (const RefPtr<Context>& cr);

         void
         set_foreground_ready (const bool ready);

         void
         blit_foreground_buffer (const RefPtr<Context>& cr);

         virtual void
         render_foreground_buffer ();

         virtual void
         render_foreground_buffer (const RefPtr<Context>& cr);

   };

   enum Button_State
   {
      BUTTON_ON,
      BUTTON_PRESSED,
      BUTTON_OFF
   };

   class Dbutton : public Dwidget
   {

      friend class Dtoggle_Button;

      public:

         typedef sigc::signal<void>
         Signal;

         typedef sigc::signal<void, const string&>
         Str_Signal;

         typedef sigc::signal<void, const Dmouse_Button_Event>
         Full_Signal;

         typedef sigc::signal<void, const string&, const Dmouse_Button_Event>
         Full_Str_Signal;

      protected:

         bool
         disabled;

         Button_State
         state;

         const string
         str;

         const string
         separator;

         const Symbol
         symbol;

         Color
         band_color;

         Color
         led_color;

         Signal
         signal;

         Str_Signal
         str_signal;

         Full_Signal
         full_signal;

         Full_Str_Signal
         full_str_signal;

         Color
         get_bg_color () const;

         Color
         get_fg_color () const;

         void
         render_background (const RefPtr<Context>& cr) const;

         void
         render_band (const RefPtr<Context>& cr) const;

         void
         render_context (const RefPtr<Context>& cr) const;

         void
         render_led (const RefPtr<Context>& cr) const;

      public:

         Dbutton (const Dcanvas& canvas,
                  const string& str,
                  const Real font_size,
                  const string& separator = "\t\n");

         Dbutton (const Dcanvas& canvas,
                  const Symbol& symbol);

         virtual string
         get_str () const;

         Signal&
         get_signal ();

         Str_Signal&
         get_str_signal ();

         Full_Signal&
         get_full_signal ();

         Full_Str_Signal&
         get_full_str_signal ();

         virtual bool
         on_mouse_button_pressed (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_button_released (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_enter (const Dmouse_Crossing_Event& event);

         virtual bool
         on_mouse_exit (const Dmouse_Crossing_Event& event);

         virtual void
         clicked (const Dmouse_Button_Event& event);

         virtual void
         cairo (const RefPtr<Context>& cr);

         void
         set_band_color (const Color& band_color);

         void
         set_led_color (const Color& led_color);

         void
         enable ();

         void
         disable ();

   };

   template<class T>
   class Template_Button : public Dbutton
   {

      public:

         typedef sigc::signal<void, const T&>
         T_Signal;

         typedef sigc::signal<void, const T&, const Dmouse_Button_Event&>
         Full_T_Signal;

      protected:

         T
         t;

         T_Signal
         t_signal;

         Full_T_Signal
         full_t_signal;

      public:

         Template_Button (const Dcanvas& canvas,
                          const T t,
                          const string& str,
                          const Real font_size)
            : Dbutton (canvas, str, font_size),
              t (t)
         {
         }

         Template_Button (const Dcanvas& canvas,
                          const T t,
                          const Real font_size)
            : Dbutton (canvas, t.get_string (), font_size),
              t (t)
         {
         }

         const T&
         get_t () const
         {
            return t;
         }

         T_Signal&
         get_t_signal ()
         {
            return t_signal;
         }

         Full_T_Signal&
         get_full_t_signal ()
         {
            return full_t_signal;
         }

         void
         clicked (const Dmouse_Button_Event& event)
         {
            t_signal.emit (t);
            full_t_signal.emit (t, event);
         }

   };

   class Spin_Button : public Dbutton,
                       public Tokens
   {

      public:

         typedef sigc::signal<void>
         Update_Signal;

         typedef sigc::signal<void, const string&>
         Update_Str_Signal;

         typedef sigc::signal<void, const Devent&>
         Full_Update_Signal;

         typedef sigc::signal<void, const string&, const Devent>
         Full_Update_Str_Signal;

      protected:

         Tokens::iterator
         iterator;

         Tokens::iterator
         official_iterator;

         Update_Signal
         update_signal;

         Update_Str_Signal
         update_str_signal;

         Full_Update_Signal
         full_update_signal;

         Full_Update_Str_Signal
         full_update_str_signal;

         bool
         instant_update;

         virtual bool
         on_mouse_scroll (const Dmouse_Scroll_Event& event);

         virtual void
         clicked (const Dmouse_Button_Event& event);

      public:

         Spin_Button (const Dcanvas& canvas,
                      const bool instant_update,
                      const Real font_size,
                      const string& separator = "\t\n");

         string
         get_official_str () const;

         string
         get_str () const;

         Update_Signal&
         get_update_signal ();

         Update_Str_Signal&
         get_update_str_signal ();

         Full_Update_Signal&
         get_full_update_signal ();

         Full_Update_Str_Signal&
         get_full_update_str_signal ();

         void
         set_instant_update (const bool instant_update);

         void
         clear ();

         void
         set (const Tokens& tokens,
              const bool try_preserve = true);

         void
         add_token (const string& token,
                    const bool clear_first = false);

         void
         add_tokens (const Tokens& tokens,
                     const bool clear_first = false);

         void
         increment ();

         void
         decrement ();

         void
         updated (const Devent& event);

   };

   class Dtoggle_Button : public Dbutton
   {

      public:

         typedef sigc::signal<void>
         Signal;

         typedef sigc::signal<void, const string&>
         Str_Signal;

         typedef sigc::signal<void, const Dmouse_Button_Event>
         Full_Signal;

         typedef sigc::signal<void, const string&, const Dmouse_Button_Event>
         Full_Str_Signal;

      protected:

         bool
         switched_on;

         Signal
         signal;

         Str_Signal
         str_signal;

         Full_Signal
         full_signal;

         Full_Str_Signal
         full_str_signal;

      public:

         Dtoggle_Button (const Dcanvas& canvas,
                         const string& str,
                         const Real font_size,
                         const bool switched_on = false);

         Signal&
         get_signal ();

         Str_Signal&
         get_str_signal ();

         Full_Signal&
         get_full_signal ();

         Full_Str_Signal&
         get_full_str_signal ();

         Dbutton::Signal&
         get_click_signal ();

         Dbutton::Str_Signal&
         get_click_str_signal ();

         Dbutton::Full_Signal&
         get_click_full_signal ();

         Dbutton::Full_Str_Signal&
         get_click_full_str_signal ();

         void
         set (const bool switched_on);

         virtual void
         toggle ();

         const bool&
         is_switched_on () const;

         virtual bool
         on_mouse_button_pressed (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_button_released (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_enter (const Dmouse_Crossing_Event& event);

         virtual bool
         on_mouse_exit (const Dmouse_Crossing_Event& event);

         virtual void
         toggled (const Dmouse_Button_Event& event);

   };

   template<class T>
   class Template_Toggle_Button : public Dbutton
   {

      public:

         typedef sigc::signal<void>
         T_Signal;

         typedef sigc::signal<void, const T&, const Dmouse_Button_Event>
         Full_T_Signal;

      protected:

         T
         t;

         T_Signal
         t_signal;

         Full_T_Signal
         full_t_signal;

      public:

         Template_Toggle_Button (const Dcanvas& canvas,
                                 const T t,
                                 const string& str,
                                 const Real font_size,
                                 const bool switched_on)
            : Dtoggle_Button (canvas, str, font_size, switched_on),
              t (t)
         {
         }

         Template_Toggle_Button (const Dcanvas& canvas,
                                 const T t,
                                 const Real font_size,
                                 const bool switched_on)
            : Dtoggle_Button (canvas, t.get_string (), font_size, switched_on),
              t (t)
         {
         }

         const T&
         get_t () const
         {
            return t;
         }

         Signal&
         get_t_signal ()
         {
            return t_signal;
         }

         Full_Signal&
         get_full_t_signal ()
         {
            return full_t_signal;
         }

         void
         toggled (const Dmouse_Button_Event& event)
         {
            t_signal.emit (t);
            full_t_signal.emit (t, event);
         }

   };

   class Radio_Button : public Dtoggle_Button
   {

      public:

         class Group
         {

            private:

               vector<Radio_Button*>
               radio_button_ptr_vector;

            public:

               Group (Radio_Button& radio_button)
               {
                  add (radio_button);
               };

               Integer
               add (Radio_Button& radio_button);

               void
               set_active (const Integer index);

         } group;

         Group*
         group_ptr;

         Integer
         group_index;

         Radio_Button (const Dcanvas& canvas,
                       const string& str,
                       const Real font_size);

         Group&
         get_group ();

         void
         set_group (Group& group);

         virtual void
         set_active (); 

         virtual void
         toggled (const Dmouse_Button_Event& event);

   };

   class Drawer : public Dtoggle_Button
   {

      protected:

         Dcanvas&
         dcanvas;

         Dv_Pack_Box
         widget_panel;

         Dcontainer&
         dcontainer;

         bool
         open_upwards;

         map<Integer, bool>
         ref_only_map;

         map<Integer, Dwidget*>
         widget_ptr_map;

         void
         pack ();

         void
         clear ();

      public:

         Drawer (Dcanvas& dcanvas,
                 Dcontainer& dcontainer,
                 const bool open_upwards,
                 const string& str,
                 const Real font_size);

         ~Drawer ();

         const Dwidget&
         get_widget (const Integer index) const;

         Dwidget&
         get_widget (const Integer index);

         virtual void
         set_hidable (const bool hidable);

         void
         add_widget (Dwidget& widget);

         void
         add_widget_ptr (Dwidget* widget_ptr);

         void
         being_packed (const Point_2D& anchor,
                       const Real width,
                       const Real height);

         void
         toggled (const Dmouse_Button_Event& event);

         void
         toggle ();

         void
         expand ();

         void
         collapse ();

   };

   class Drawer_Panel : public Dh_Pack_Box
   {

      public:

         typedef std::map<string, Drawer*>
         Drawer_Ptr_Map;

         typedef std::map<string, Dbutton*>
         Button_Ptr_Map;

         typedef std::map<string, Dtoggle_Button*>
         Toggle_Button_Ptr_Map;

      protected:

         Dcanvas&
         dcanvas;

         const bool
         open_upwards;

         const Real
         font_size;

         Drawer_Ptr_Map
         drawer_ptr_map;

         Button_Ptr_Map
         button_ptr_map;

         Toggle_Button_Ptr_Map
         toggle_button_ptr_map;

      public:

         Drawer_Panel (Dcanvas& dcanvas,
                       const bool open_upwards,
                       const Real font_size = 12);

         ~Drawer_Panel ();

         Drawer_Ptr_Map::iterator
         add_drawer_ptr (Drawer* drawer_ptr,
                         const bool back = true);

         Drawer_Ptr_Map::iterator
         add_drawer (const string& str,
                     const bool back = true);

         void
         add_widget (const string& drawer_str,
                     Dwidget& widget);

         void
         add_widget_ptr (const string& drawer_str,
                         Dwidget* widget_ptr);

         void
         add_button (const string& str,
                     const bool back = true);

         void
         add_toggle_button (const string& str,
                            const bool switched_on,
                            const bool back = true);

         const Drawer&
         get_drawer (const string& str) const;

         Drawer&
         get_drawer (const string& str);

         const Dbutton&
         get_button (const string& str) const;

         Dbutton&
         get_button (const string& str);

         const Dtoggle_Button&
         get_toggle_button (const string& str) const;

         Dtoggle_Button&
         get_toggle_button (const string& str);

         bool
         is_switched_on (const string& str) const;

         void
         toggle (const string& str);

         void
         set_toggle (const string& str,
                     const bool switched_on);

         bool
         all_collapsed () const;

         void
         expand_all ();

         void
         collapse_all ();

   };

   class Dtitle : public Dwidget
   {

      protected:

         const Real
         margin;

         string
         string_l;

         string
         string_c;

         string
         string_r;

      public:

         Dtitle (const Dcanvas& canvas,
                 const Real font_size);

         Dtitle (const Dcanvas& canvas,
                 const Real font_size,
                 const string& string_l,
                 const string& string_c,
                 const string& string_r);

         void
         set_string_l (const string& string_l);

         void
         set_string_c (const string& string_c);

         void
         set_string_r (const string& string_r);

         virtual Real
         get_preferred_width () const;

         virtual Real
         get_preferred_height () const;

         virtual void
         cairo (const RefPtr<Context>& cr);

   };

   class Popup : public Dwidget
   {

      protected:

         Orientation
         orientation;

         Tokens
         tokens;

         Real
         font_size;

         Real
         margin;

         Integer
         index;

         Integer
         get_index (const Point_2D& point) const;

         void
         cairo (const RefPtr<Context>& cr,
                const Integer index,
                const FontExtents& fe) const;

      public:

         Popup (const Dcanvas& canvas,
                const Real font_size);

         void
         hide ();

         void
         clear ();

         void
         append (const string& str);

         void
         set_tokens (const Tokens& tokens);

         void
         set_shape (const Point_2D& anchor,
                    const Orientation orientation,
                    const Real width,
                    const Real height);

         Real
         get_preferred_width () const;

         Real
         get_preferred_height () const;

         void
         reset_index ();

         void
         cairo (const RefPtr<Context>& cr);

   };

   class Popup_Menu : public Popup
   {

      public:

         typedef sigc::signal<void>
         Signal;

         typedef sigc::signal<void, const string&>
         Str_Signal;

      protected: 

         bool
         fresh;

         map<string, Signal>
         signal_map;

         map<string, Str_Signal>
         str_signal_map;

         virtual void
         emit_signal (const string& str) const;

      public:

         Popup_Menu (const Dcanvas& canvas,
                     const Real font_size);

         const bool
         is_fresh () const;

         void
         set_fresh (const bool fresh);

         void
         clear ();

         virtual void
         append (const string& str);

         Signal&
         get_signal (const string& str);

         Str_Signal&
         get_str_signal (const string& str);

         bool
         is_on () const;

         virtual bool
         on_mouse_motion (const Dmouse_Motion_Event& event);

         virtual bool
         on_mouse_button_released (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_exit (const Dmouse_Crossing_Event& event);

   };

   class Time_Chooser : public Dcontainer
   {

      public:

         class Shape : public set<Dtime>
         {

            public:

               Dtime
               start_time;

               Dtime
               end_time;

               Real
               leap;

               Shape ();

               Shape (const Shape& shape);

               Shape (const set<Dtime>& time_set);

               Real
               get_span_t () const;

               bool
               out_of_bounds (const Dtime& dtime) const;

               Shape::const_iterator
               get_iterator (const Dtime& dtime) const;

               const Dtime&
               get_nearest_time (const Dtime& dtime) const;

               const Dtime&
               get_next_time (const Dtime& dtime,
                              const bool is_leap = false) const;

               const Dtime&
               get_prev_time (const Dtime& dtime,
                              const bool is_leap = false) const;

               bool
               operator== (const Time_Chooser::Shape& shape) const;

         };

         class Data : public Tokens
         {

            private:

               void
               init (const string& str);

               static void
               rectify (string& str);

               static bool
               contains (const string& str,
                         const Dtime& dtime);

            public:

               Dtime
               dtime;

               Dtime
               candidate_time;

               Data ();

               Data (const Data& data);

               Data (const Tokens& tokens);

               Data (const Dtime& dtime);

               Data (const string& str);

               void
               insert (const Dtime& dtime,
                       const bool range);

               void
               clear ();

               void
               set_time (const Dtime& dtime,
                         const bool clear_first,
                         const bool is_range,
                         const bool is_leap);

               void
               set_time (const Dtime& dtime,
                         const Devent& = Devent ());

               Shape
               get_highlighted_shape (const Shape& shape) const;

               bool
               matches (const Dtime& dtime) const;

               bool
               is_candidate (const Dtime& dtime) const;

               bool
               contains (const Dtime& dtime) const;

               void
               conform_to (const Shape& shape);

               void
               dump () const;

         };

         typedef sigc::signal<void>
         Signal;

         typedef sigc::signal<void, const Time_Chooser::Data&>
         Data_Signal;

      protected:

         const Real
         margin;

         bool
         reverse;

         bool
         vertical;

         Box_2D
         button_box;

         Real
         line_interval;

         Shape
         shape;

         Data
         data;

         Integer
         dwell;

         Dtoggle_Button
         anim_button;

         Dbutton
         prev_button;

         Dbutton
         next_button;

         Signal
         signal;

         Data_Signal
         data_signal;

         set<Integer>
         additional_indices;

         Glib::Thread*
         thread_ptr;

         Real
         get_s (const Dtime& dtime) const;

         Dtime
         get_time (const Point_2D& point,
                   const bool snap_to_nearest) const;

         void
         render_background (const RefPtr<Context>& cr);

         void
         render_background_months (const RefPtr<Context>& cr);

         void
         render_month (const RefPtr<Context>& cr,
                       const Integer year,
                       const Integer month) const;

         void
         render_month_lines (const RefPtr<Context>& cr,
                                   const Integer year,
                                   const Integer month,
                                   const string& format = "") const;

         void
         render_background_days (const RefPtr<Context>& cr);

         void
         render_dates (const RefPtr<Context>& cr) const;

         void
         render_lines (const RefPtr<Context>& cr,
                       const Real hours,
                       const string& format = "") const;

         void
         render_nodes (const RefPtr<Context>& cr) const;

         void
         render_now (const RefPtr<Context>& cr) const;

         void
         init (const Real font_size);

         void
         animate ();

      public:

         Time_Chooser (const Dcanvas& canvas,
                       const Real font_size,
                       const bool reverse = false,
                       const bool vertical = true,
                       const Real line_interval = 6);

         ~Time_Chooser ();

         void
         set_vertical (const bool vertical);

         const Time_Chooser::Shape&
         get_shape () const;

         Time_Chooser::Shape&
         get_shape ();

         const Time_Chooser::Data&
         get_data () const;

         Time_Chooser::Data&
         get_data ();

         Signal&
         get_signal ();

         Data_Signal&
         get_data_signal ();

         Shape
         get_highlighted_shape () const;

         void
         set_reverse (const bool reverse);

         void
         set_line_interval (const Real line_interval);

         void
         set_time_chooser (const Time_Chooser& time_chooser);

         bool
         set_data (const Time_Chooser::Data& data);

         void
         set_shape (const Shape& shape);

         const Dtime&
         get_time () const;

         void
         advance (const Devent& event);

         void
         handle_animate (const Devent& event = Devent ());

         void
         handle_dwell (const Devent& event = Devent ());

         virtual void
         step_forward (const Devent& event = Devent ());

         virtual void
         step_backward (const Devent& event = Devent ());

         virtual void
         go_to_first ();

         virtual void
         go_to_last ();

         void
         toggle_anim_button ();

         void
         set_dwell (const Integer dwell);

         void
         set_leap (const Real leap);

         Point_2D
         get_point (const Dtime& dtime) const;

         Rect
         get_rect (const Dtime& dtime) const;

         virtual void
         cairo (const RefPtr<Context>& cr);

         Real
         get_preferred_width () const;

         Real
         get_preferred_height () const;

         virtual void
         pack ();

         bool
         on_mouse_button_pressed (const Dmouse_Button_Event& event);

         bool
         on_mouse_motion (const Dmouse_Motion_Event& event);

         bool
         on_mouse_button_released (const Dmouse_Button_Event& event);

         virtual void
         selected ();

   };

   class Time_Canvas
   {

      protected:

         const Dcanvas&
         canvas;

         Time_Chooser
         time_chooser;

      public:

         Time_Canvas (const Dcanvas& canvas,
                      const Real font_size,
                      const bool reverse = false,
                      const bool vertical = true,
                      const Real line_interval = 6);

         const Time_Chooser&
         get_time_chooser () const;

         Time_Chooser&
         get_time_chooser ();

         bool
         on_key_pressed (const Dkey_Event& event);

   };

   class Level_Panel : public Dv_Pack_Box
   {

      public:

         typedef sigc::signal<void, const Level&>
         Level_Signal;

         typedef sigc::signal<void, const Layer&>
         Layer_Signal;

         typedef sigc::signal<void, const Level&, const Devent&>
         Full_Level_Signal;

         typedef sigc::signal<void, const Layer&, const Devent&>
         Full_Layer_Signal;

      protected:

         Dcanvas&
         dcanvas;

         const Real
         font_size;

         Real
         reference_y;

         Level
         level;

         Level
         candidate_level;

         Real
         z;

         Real
         pressure;

         Real
         theta;

         Real
         sigma;

         bool
         setting_level;

         bool
         setting_level_0;

         bool
         setting_level_1;

         const Real
         start_margin;

         const Real
         end_margin;

         const Level_Tuple
         level_tuple_z;

         const Level_Tuple
         level_tuple_p;

         const Level_Tuple
         level_tuple_theta;

         const Level_Tuple
         level_tuple_sigma;

         vector<Level>
         extra_level_vector;

         Level_Signal
         level_signal;

         Layer_Signal
         layer_signal;

         Full_Level_Signal
         full_level_signal;

         Full_Layer_Signal
         full_layer_signal;

         //Dh_Pack_Box
         //switcher_panel;

         //Dbutton
         //pressure_button;

         //Dbutton
         //theta_button;

         //Dbutton
         //sigma_button;

         const Level_Tuple&
         get_level_tuple () const;

         Integer
         get_extra_level_index () const;

         Level
         get_level (const Real y) const;

         Level
         get_level_extra (const Real y) const;

         Level
         get_level_z (const Real y) const;

         Level
         get_level_p (const Real y) const;

         Level
         get_level_theta (const Real y) const;

         Level
         get_level_sigma (const Real y) const;

         Real
         get_y (const Level& level) const;

         Real
         get_y_z (const Real z) const;

         Real
         get_y_p (const Real p) const;

         Real
         get_y_theta (const Real theta) const;

         Real
         get_y_sigma (const Real sigma) const;

         bool
         is_close_to_level_0 (const Real y) const;

         bool
         is_close_to_level_1 (const Real y) const;

         bool
         is_within_layer (const Real y) const;

         void
         render_background (const RefPtr<Context> cr) const;

         void
         render_level (const RefPtr<Context> cr,
                       const Level& level,
                       const bool with_triangles) const;

         void
         render_layer (const RefPtr<Context> cr,
                       const Level& level) const;

         bool
         on_mouse_button_pressed (const Dmouse_Button_Event& event);

         bool
         on_mouse_motion (const Dmouse_Motion_Event& event);

         bool
         on_mouse_button_released (const Dmouse_Button_Event& event);

      public:

         Level_Panel (Dcanvas& dcanvas,
                      const Real font_size);

         Level_Signal&
         get_level_signal ();

         Full_Level_Signal&
         get_full_level_signal ();

         void
         add_extra_level (const Level& level);

         void
         cairo (const RefPtr<Context>& cr);

         void
         set_level (const Level& level);

         void
         move_level_up (const Devent& event);

         void
         move_level_down (const Devent& event);

         const Level&
         get_level () const;

         Level&
         get_level ();

         void
         clear ();

         void
         set_no_buttons ();

         void
         set_theta_buttons ();

         void
         set_sigma_buttons ();

         void
         set_p_buttons ();

         void
         set_wind_level_buttons ();

         void
         set_temperature_level_buttons ();

   };

   class Level_Canvas
   {

      protected:

         Dcanvas&
         canvas;

         Level_Panel
         level_panel;

      public:

         Level_Canvas (Dcanvas& canvas,
                       const Real font_size);

         const Level_Panel&
         get_level_panel () const;

         Level_Panel&
         get_level_panel ();

         bool
         on_key_pressed (const Dkey_Event& event);

   };

   class Console_2D : public Dcanvas
   {

      protected:

         class Manipulate_Data
         {

            public:

               Point_2D
               point;

               Real
               tilt;

               Manipulate_Data ();

               Manipulate_Data (const Manipulate_Data& manipulate_data);

               Manipulate_Data (const Point_2D& point,
                                const Real tilt);

               void
               reset ();

               void
               set_point (const Point_2D& point);

               void
               set_tilt (const Real tilt);

               bool
               is_translating () const;

               bool
               is_rotating () const;

               bool
               is_manipulating () const;

         };

         class Zoom_Box
         {

            friend class Console_2D;
            friend class Map_Console;

            protected:

               class Scale_Ratio
               {

                  public:

                     const Real
                     fine;

                     const Real
                     coarse;

                     Scale_Ratio (const Real fine,
                                  const Real coarse);

               };

               Console_2D&
               console_2d;

               Color
               color;

               bool
               visible;

               Manipulate_Data
               manipulate_data;

               Point_2D
               centre;

               Real
               zoom_factor_x;

               Real
               zoom_factor_y;

               Real
               tilt;

            public:

               const Scale_Ratio
               scale_ratio;

               Zoom_Box (Console_2D& console_2d);

               ~Zoom_Box ();

               void
               set_color (const Color& color);

               virtual void
               apply_to (Affine_Transform_2D& transform) const;

               virtual void
               pre_manipulate (const Dmouse_Button_Event& event);

               virtual void
               manipulate (const Dmouse_Motion_Event& event);

               virtual void
               post_manipulate ();

               virtual void
               scroll (const Dmouse_Scroll_Event& event);

               virtual void
               reset ();

               virtual void
               toggle_visible ();

               virtual void
               set_visible (const bool visible);

               virtual bool
               is_visible () const;

               virtual bool
               contains (const Point_2D& point) const;

               virtual void
               cairo (const RefPtr<Context>& cr);

         };

      public:

         class Hud
         {

            protected:

               const Integer
               id;

               const Real
               node_size;

               sigc::signal<void>
               moved_signal;

               sigc::signal<void>
               settled_signal;

            public:

               class Popup_Menu : public denise::Popup_Menu
               {

                  public:

                     typedef sigc::signal<void, const Integer>
                     Id_Signal;

                  protected:

                     Integer
                     id;

                     Console_2D&
                     console_2d;

                     sigc::signal<void, const Integer>
                     clear_signal;

                     std::map<string, Id_Signal>
                     id_signal_map;

                     virtual void
                     emit_signal (const string& str) const;

                  public:

                     Popup_Menu (Console_2D& console_ed,
                                 const Real font_size);

                     void
                     clear ();

                     virtual void
                     append (const string& str);

                     Id_Signal&
                     get_id_signal (const string& str);

                     virtual void
                     setup (const Hud& hud,
                            const Point_2D& point);

               };

               Hud (const Integer id,
                    const Real node_size);

               virtual ~Hud ();

               const Integer
               get_id () const;

               virtual bool
               matches (const Transform_2D& transform,
                        const Point_2D& point) const = 0;

               virtual void
               cairo (const RefPtr<Context>& cr,
                      const Console_2D& console) const = 0;

               sigc::signal<void>&
               get_moved_signal ();

               sigc::signal<void>&
               get_settled_signal ();

               virtual void
               double_clicked ();

         };

         class Marker : public Hud,
                        public Point_2D
         {

            protected:

               string
               str;

            public:

               class Popup_Menu : public Hud::Popup_Menu
               {

                  public:

                     Popup_Menu (Console_2D& console_2d,
                                 const Real font_size,
                                 const bool connect_default_signal = true);

               };

               Marker (const Integer id,
                       const Point_2D& point,
                       const Real node_size = 6);

               virtual void
               attract_by (const Attractor& attractor);

               const string&
               get_str () const;

               void
               set_str (const string& str);

               bool
               matches (const Transform_2D& transform,
                        const Point_2D& point) const;

               virtual void
               cairo (const RefPtr<Context>& cr,
                      const Console_2D& console_2d) const;

         };

         class Route : public Hud,
                       public Simple_Polyline
         {

            protected:

               Point_2D
               origin;

            public:

               class Popup_Menu : public Hud::Popup_Menu
               {

                  public:

                     Popup_Menu (Console_2D& console_2d,
                                 const Real font_size,
                                 const bool connect_default_signal = true);

               };

               Route (const Integer id,
                      const Real node_size = 8);

               Route (const Integer id,
                      const Point_2D& point,
                      const Real node_size = 8);

               Route (const Integer id,
                      const Point_2D& point_a,
                      const Point_2D& point_b,
                      const Real node_size = 8);

               virtual bool
               is_too_short () const;

               virtual void
               translate (const Point_2D& from,
                          const Point_2D& to);

               virtual void
               translate (const Point_2D& to);

               void
               set_origin (const Point_2D& point);

               const Point_2D&
               get_origin () const;

               bool
               has_origin () const;

               virtual bool
               matches (const Transform_2D& transform,
                        const Point_2D& point) const;

               virtual void
               cairo (const RefPtr<Context>& cr,
                      const Console_2D& console_2d) const;

               virtual Route::iterator
               add (const Transform_2D& transform,
                    const Point_2D& point);

               virtual Route::const_iterator
               get_node (const Transform_2D& transform,
                         const Point_2D& point) const;

               virtual Route::iterator
               get_node (const Transform_2D& transform,
                         const Point_2D& point);

         };

         class Shape : public Route
         {

            protected:

               Polygon*
               get_polygon_ptr () const;

               bool
               on_node (const Transform_2D& transform,
                        const Point_2D& point) const;

            public:

               class Popup_Menu : public Hud::Popup_Menu
               {

                  public:

                     Popup_Menu (Console_2D& console_2d,
                                 const Real font_size,
                                 const bool connect_default_signal = true);

               };

               Shape (const Integer id,
                      const Real node_size = 5);

               Shape (const Integer id,
                      const Point_2D& point,
                      const Real node_size = 5);

               virtual bool
               is_too_short () const;

               virtual bool
               contains (const Point_2D& point) const;

               virtual bool
               contains (const Transform_2D& transform,
                         const Point_2D& point) const;

               bool
               on_edge (const Transform_2D& transform,
                        const Point_2D& point) const;

               virtual bool
               matches (const Transform_2D& transform,
                        const Point_2D& point) const;

               virtual void
               cairo (const RefPtr<Context>& cr,
                      const Console_2D& console_2d) const;

         };

         class Hud_Store : public std::map<Integer, Hud*>
         {

            protected:

               Integer
               dragging_id;

               Integer
               get_first_available_id () const;

               const Attractor*
               attractor_ptr;

            public:

               Hud_Store ();

               virtual ~Hud_Store ();

               void
               set_attractor (const Attractor& attractor);

               virtual const Attractor&
               get_attractor () const;

               virtual void
               clear ();

               virtual void
               reset ();

               virtual const Hud&
               get_hud (const Integer id) const;

               virtual Hud&
               get_hud (const Integer id);

               Integer
               get_id (const Transform_2D& transform,
                       const Point_2D& point);

               virtual void
               insert (const Integer id,
                       Hud* hud_ptr);

               virtual void
               erase (const Integer id);

               virtual void
               cairo (const RefPtr<Context>& cr,
                      const Console_2D& console_2d) const;

               virtual bool
               show_popup_menu (Console_2D& console_2d,
                                const Point_2D& point,
                                Hud::Popup_Menu& popup_menu);

         };

         class Marker_Store : public Hud_Store
         {

            protected:

               bool
               adding_marker;

            public:

               Marker_Store ();

               virtual ~Marker_Store ();

               void
               reset ();

               virtual const Marker&
               get_marker (const Integer id) const;

               virtual Marker&
               get_marker (const Integer id);

               virtual Integer
               insert (const Point_2D& );

               virtual bool
               button_1_pressed (Console_2D& console_2d,
                                 const Dmouse_Button_Event& event);

               virtual bool
               try_insert (Console_2D& console_2d,
                           const Dmouse_Button_Event& event);

               virtual bool
               button_1_released (Console_2D& console_2d,
                                  const Dmouse_Button_Event& event);

               virtual bool
               button_3_released (Console_2D& console_2d,
                                  const Dmouse_Button_Event& event);

               virtual bool
               try_dragging (Console_2D& console_2d,
                             const Dmouse_Motion_Event& event);

         };

         class Route_Store : public Hud_Store
         {

            protected:

               Route
               dummy;

               Route::iterator
               node;

               virtual Console_2D::Route*
               new_route_ptr (const Integer id,
                              const Point_2D& point,
                              const Real node_size = 8);

               virtual Console_2D::Route*
               new_route_ptr (const Integer id,
                              const Point_2D& point_a,
                              const Point_2D& point_b,
                              const Real node_size = 8);

            public:

               Route_Store ();

               virtual ~Route_Store ();

               void
               reset ();

               virtual const Route&
               get_route (const Integer id) const;

               virtual Route&
               get_route (const Integer id);

               Route::iterator
               get_node (const Transform_2D& transform,
                         const Point_2D& point);

               virtual Integer
               insert (const Point_2D& point);

               virtual Integer
               insert (const Point_2D& point_a,
                       const Point_2D& point_b);

               virtual bool
               button_1_pressed (Console_2D& console_2d,
                                 const Dmouse_Button_Event& event);

               virtual bool
               button_2_pressed (Console_2D& console_2d,
                                 const Dmouse_Button_Event& event);

               virtual bool
               button_1_released (Console_2D& console_2d,
                                  const Dmouse_Button_Event& event);

               virtual bool
               button_2_released (Console_2D& console_2d,
                                  const Dmouse_Button_Event& event);

               virtual bool
               try_dragging (Console_2D& console_2d,
                             const Dmouse_Motion_Event& event);

         };

         class Shape_Store : public Route_Store
         {

            public:

               Shape_Store ();

               virtual const Shape&
               get_shape (const Integer id) const;

               virtual Shape&
               get_shape (const Integer id);

               virtual Integer
               insert (const Point_2D& point);

               virtual bool
               button_1_pressed (Console_2D& console_2d,
                                 const Dmouse_Button_Event& event);

               virtual bool
               button_1_released (Console_2D& console_2d,
                                  const Dmouse_Button_Event& event);

         };

      protected:

         Affine_Transform_2D
         affine_transform;

         Marker::Popup_Menu
         marker_popup_menu;

         Route::Popup_Menu
         route_popup_menu;

         Shape::Popup_Menu
         shape_popup_menu;

         Marker_Store
         marker_store;

         Route_Store
         route_store;

         Shape_Store
         shape_store;

         bool
         ignore_heavy;

         bool
         computer_convention;

         Manipulate_Data
         manipulate_data;

         Zoom_Box
         zoom_box;

         virtual const Zoom_Box&
         get_zoom_box () const;

         virtual Zoom_Box&
         get_zoom_box ();

         virtual void
         reset_transform ();

         virtual void
         pre_manipulate (const Dmouse_Button_Event& event);

         virtual void
         manipulate (const Dmouse_Motion_Event& event);

         virtual void
         post_manipulate ();

         virtual bool
         on_key_pressed (const Dkey_Event& event);

         virtual bool
         on_mouse_button_pressed (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_button_released (const Dmouse_Button_Event& event);

         virtual bool
         on_mouse_motion (const Dmouse_Motion_Event& event);

         virtual bool
         on_mouse_scroll (const Dmouse_Scroll_Event& event);

         virtual Marker::Popup_Menu&
         get_marker_popup_menu ();

         virtual Route::Popup_Menu&
         get_route_popup_menu ();

         virtual Shape::Popup_Menu&
         get_shape_popup_menu ();

         virtual string
         get_string (const Marker& marker) const;

         virtual Tokens
         get_tokens (const Marker& marker) const;

      public:

         Console_2D (Gtk::Window& gtk_window,
                     const Size_2D& size_2d,
                     const bool computer_convention = true);

         ~Console_2D ();

         const Affine_Transform_2D&
         get_affine_transform () const;

         Affine_Transform_2D&
         get_affine_transform ();

         virtual const Transform_2D&
         get_transform () const;

         virtual Transform_2D&
         get_transform ();

         virtual void
         set_hide_hidable (const bool hide_hidable);

         virtual void
         clear_marker (const Integer marker_id);

         virtual void
         clear_route (const Integer route_id);

         virtual void
         clear_shape (const Integer shape_id);

         virtual void
         refresh_all ();

         virtual void
         total_refresh ();

         virtual void
         translate (const Real dx,
                    const Real dy);

         virtual void
         scroll (const Dmouse_Scroll_Event& event);

         bool
         is_translating () const;

         bool
         is_rotating () const;

         bool
         is_manipulating () const;

         virtual void
         apply_zoom_box ();

         virtual const Console_2D::Marker_Store&
         get_marker_store () const;

         virtual Console_2D::Marker_Store&
         get_marker_store ();

         virtual const Console_2D::Route_Store&
         get_route_store () const;

         virtual Console_2D::Route_Store&
         get_route_store ();

         virtual const Console_2D::Shape_Store&
         get_shape_store () const;

         virtual Console_2D::Shape_Store&
         get_shape_store ();

         virtual void
         render_huds_and_widgets (const RefPtr<Context>& cr);

         virtual void
         cairo (const RefPtr<Context>& cr);

         virtual void
         save_image (const string& file_path);

         virtual void
         save_png (const string& file_path);

         virtual void
         save_svg (const string& file_path);

   };

   class Map_Console : public Console_2D
   {

      //protected:
      public:

         class Zoom_Box : public Console_2D::Zoom_Box
         {

            friend class Map_Console;

            protected:

               Map_Console&
               map_console;

            public:

               Zoom_Box (Map_Console& map_console);

               ~Zoom_Box ();

               void
               apply_to (Geodetic_Transform& transform) const;

               void
               pre_manipulate (const Dmouse_Button_Event& event);

               void
               manipulate (const Dmouse_Motion_Event& event);

               void
               post_manipulate ();

               void
               scroll (const Dmouse_Scroll_Event& event);

         };

         class Overlay
         {

            private:

               Geodetic_Cairoable*
               gc_ptr;

               void
               init (Geodetic_Cairoable* gc_ptr, 
                     const bool filled,
                     const Real line_width,
                     const Color& color,
                     const bool heavy);

            public:

               bool
               heavy;

               bool
               filled;

               Real
               line_width;

               Color
               color;

               Overlay (Geodetic_Cairoable* gc_ptr,
                        const bool filled,
                        const Real line_width,
                        const Color& color,
                        const bool heavy = false);

               const void
               cairo (const RefPtr<Context>& cr,
                      const Geodetic_Transform& transform) const;

         };

      public:

         class Overlay_Store : public Tokens
         {

            private:

               class Overlay_Ptr_Map : public std::map<string, Overlay*>
               {
               } overlay_ptr_map;

               class On_Off_Map : public std::map<string, bool>
               {
               } on_off_map;

            public:

               ~Overlay_Store ();

               const Overlay&
               get_overlay (const string& identifier) const;

               void
               add (const string& overlay_string,
                    const bool heavy);

               void
               add (const string& identifier,
                    Geodetic_Cairoable* gc_ptr,
                    const bool filled,
                    const Real line_width,
                    const Color& color,
                    const bool heavy);

               bool
               is_on (const string& identifier) const;

               void
               set_on_off (const string& identifier,
                           const bool on_off);

               void
               toggle_on_off (const string& identifier);

               void
               cairo (const RefPtr<Context>& cr,
                      const Geodetic_Transform& transform,
                      const bool ignore_heavy) const;

         };

         typedef sigc::signal<void, const Geodetic_Transform&>
         Transform_Signal;

      protected:

         class Option_Panel : public Drawer_Panel
         {

            private:

               class Zoom_Drawer : public Drawer
               {

                  public:

                     Zoom_Drawer (Option_Panel& option_panel);

                     void
                     switch_off_zoom_button ();

               };

               class Overlay_Drawer : public Drawer
               {

                  private:

                     std::map<string, Integer>
                     index_map;

                  public:

                     Overlay_Drawer (Option_Panel& option_panel);

                     void
                     add_toggle_button_ptr (const string& str,
                                            Dtoggle_Button* toggle_button_ptr);

                     bool
                     is_on (const string& str) const;

                     void
                     set_on_off (const string& str,
                                 const bool on_off);

               };

               Map_Console&
               map_console;

	       void
               setup_zoom (const Tokens& config_file_content);

            public:   

               Option_Panel (Map_Console& map_console,
                             const Tokens& config_file_content);

               void
               setup_overlay ();

               void
               switch_off_zoom_button ();

               void
               align_overlay_toggles ();

               bool
               overlay_is_on (const string& overlay_str) const;

               const Map_Console::Option_Panel::Zoom_Drawer&
               get_zoom_drawer () const;

               Map_Console::Option_Panel::Zoom_Drawer&
               get_zoom_drawer ();

               const Map_Console::Option_Panel::Overlay_Drawer&
               get_overlay_drawer () const;

               Map_Console::Option_Panel::Overlay_Drawer&
               get_overlay_drawer ();

         };

      public:

         class Route_Store : public Console_2D::Route_Store
         {

            public:

               Console_2D::Route*
               new_route_ptr (const Integer id,
                              const Point_2D& point,
                              const Real node_size = 8);

               Console_2D::Route*
               new_route_ptr (const Integer id,
                              const Point_2D& point_a,
                              const Point_2D& point_b,
                              const Real node_size = 8);

         };

      protected:

         Map_Console::Route_Store
         mc_route_store;

         Geodetic_Transform*
         geodetic_transform_ptr;

         Map_Console::Zoom_Box
         geodetic_zoom_box;

         Transform_Signal
         transform_signal;

         Overlay_Store
         overlay_store;

         Option_Panel
         option_panel;

         Real
         zoom_factor;

         virtual const Zoom_Box&
         get_zoom_box () const;

         virtual Zoom_Box&
         get_zoom_box ();

         virtual bool
         on_key_pressed (const Dkey_Event& event);

         virtual void
         render_overlays (const RefPtr<Context>& cr);

         virtual void
         render_mesh (const RefPtr<Context>& cr);

         const Domain_2D
         get_domain_2d () const;

         virtual Domain_2D
         get_render_domain (const Geodetic_Vector_Data_2D& gvd_2d) const;

      public:

         class Route : public Console_2D::Route
         {

            public:

               Route (const Integer id,
                      const Point_2D& point,
                      const Real node_size = 8);

               Route (const Integer id,
                      const Point_2D& point_a,
                      const Point_2D& point_b,
                      const Real node_size = 8);

               virtual bool
               matches (const Transform_2D& transform,
                        const Point_2D& point) const;

               virtual void
               cairo (const RefPtr<Context>& cr,
                      const Console_2D& console_2d) const;

               virtual Route::iterator
               add (const Transform_2D& transform,
                    const Point_2D& point);

               Real
               get_distance (const Geodesy& geodesy) const;

         };

      protected:

         virtual string
         get_string (const Marker& marker) const;

      public:

         Map_Console (Gtk::Window& gtk_window,
                      const Size_2D& size_2d,
                      const Tokens& config_file_content);

         ~Map_Console ();

         Transform_Signal&
         get_transform_signal ();

         virtual const Map_Console::Overlay_Store&
         get_overlay_store () const;

         virtual Map_Console::Overlay_Store&
         get_overlay_store ();

         const Geodetic_Transform&
         get_geodetic_transform () const;

         Geodetic_Transform&
         get_geodetic_transform ();

         const Transform_2D&
         get_transform () const;

         Transform_2D&
         get_transform ();

         const Geodetic_Transform::Data&
         get_geodetic_transform_data () const;

         virtual void
         set_geodetic_transform_data (const Geodetic_Transform::Data& gtd);

         virtual void
         apply_zoom_box ();

         virtual void
         refresh_all ();

         virtual void
         toggle_overlay (const string& overlay_identifier);

         virtual void
         translate (const Real dx,
                    const Real dy);

         virtual void
         scroll (const Dmouse_Scroll_Event& event);

         virtual void
         align_overlay_toggles ();

         virtual void
         unify_drawers (const Devent& event);

         virtual const Console_2D::Route_Store&
         get_route_store () const;

         virtual Console_2D::Route_Store&
         get_route_store ();

   };

   class Time_Series_Canvas : public Dcanvas
   {

      protected:

         class Transform : public Affine_Transform_2D
         {

            private:

               const Real
               margin_t;

               const Real
               margin_b;

               const Real
               margin_l;

               const Real
               margin_r;

               Domain_2D
               domain_2d;

            public:

               Transform (const Size_2D& size_2d,
                          const Real margin_t,
                          const Real margin_b,
                          const Real margin_l,
                          const Real margin_r);

               void
               update (const Size_2D& size_2d);

               void
               set_domain_2d (const Domain_2D& domain_2d);

               const Domain_1D&
               get_domain_t () const;

               Domain_1D&
               get_domain_t ();

               const Domain_1D&
               get_domain_y () const;

               Domain_1D&
               get_domain_y ();

         };

         Transform
         transform;

         Mesh_2D
         mesh_2d;

         string
         fmt_y;

         void
         update_transform ();

         bool
         on_configure_event (GdkEventConfigure* event);

         virtual void
         render_background_buffer (const RefPtr<Context>& cr);

         virtual void
         render_time_series (const RefPtr<Context>& cr,
                             const Vector_Data_1D& time_series,
                             const Integer vector_element) const;

         virtual void
         render_sun_elevation (const RefPtr<Context>& cr,
                               const Lat_Long& lat_long) const;

      public:

         Time_Series_Canvas (Gtk::Window& gtk_window,
                             const Size_2D& size_2d,
                             const Real margin_t = 0,
                             const Real margin_b = 0,
                             const Real margin_l = 0,
                             const Real margin_r = 0);

         ~Time_Series_Canvas ();

   };

};

#endif /* DENISE_GTKMM_H */

