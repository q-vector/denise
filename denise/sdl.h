//
// sdl.h
// 
// Copyright (C) 2005-2015 Simon E. Ching
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

#ifndef DENISE_SDL_H
#define DENISE_SDL_H

#include <SDL/SDL.h>
#include <denise/cairo.h>

using namespace std;

namespace denise
{

   class Sdl_Cairo : public Raster_Cairo
   {

      protected:

         SDL_Surface*
         sdl_surface_ptr;

         void
         init ();

      public:

         Sdl_Cairo (const Box_2D& box_2d);

         Sdl_Cairo (const Box_2D& box_2d,
                    const uint32_t sdl_flags);
                    

         ~Sdl_Cairo ();

         SDL_Surface*
         get_sdl_surface_ptr () const;

   };

   class Sdl_Widget : public Polygon
   {

      protected:

         bool
         focus;

         Point_2D
         origin;

         list<Sdl_Widget*>
         widget_ptr_list;

         virtual bool
         process_key_down (const SDL_Event& event);

         virtual bool
         process_key_up (const SDL_Event& event);

         virtual bool
         process_mouse_button_down (const SDL_Event& event);

         virtual bool
         process_mouse_button_up (const SDL_Event& event);

         virtual bool
         process_mouse_motion (const SDL_Event& event);

         virtual bool
         process_event (const SDL_Event& event);

         virtual bool
         handle_event (SDL_Event event);

         virtual void
         render (const Sdl_Cairo& sdl_cairo) const = 0;

      public:

         Sdl_Widget (const Point_2D& origin = Point_2D (0, 0));

         Sdl_Widget (const Polygon& polygon,
                     const Point_2D& origin = Point_2D (0, 0));

         virtual void
         paint (const Sdl_Cairo& sdl_cairo) const;

   };

   class Sdl_Label : public Label,
                     public Sdl_Widget
   {

      public:

         Sdl_Label (const string& text,
                    const Point_2D& point_2d,
                    const char justify_h,
                    const char justify_v,
                    const Point_2D& offset = Point_2D (0, 0),
                    const Real margin_x = 3,
                    const Real margin_y = 3,
                    const Point_1D text_angle = 0,
                    const bool outline = false);

         Sdl_Label (const string& text,
                    const Point_2D& point_2d,
                    const char justify_h,
                    const char justify_v,
                    const Real padding,
                    const Real margin_x = 3,
                    const Real margin_y = 3,
                    const Point_1D text_angle = 0,
                    const bool outline = false);

   };

   enum Sdl_Button_State
   {
      SDL_BUTTON_NORMAL,
      SDL_BUTTON_DEPRESSED
   };

   class Sdl_Button : public Sdl_Widget
   {

      protected:

         Sdl_Button_State
         state;

         string
         text;

         Real
         width;

         Real
         height;

         void
         render_normal (const Sdl_Cairo& sdl_cairo) const;

         void
         render_depressed (const Sdl_Cairo& sdl_cairo) const;

         void
         render (const Sdl_Cairo& sdl_cairo) const;

      public:

         Sdl_Button (const string& text,
                     const Real width,
                     const Real height);

   };

   class Sdl_Application : public Sdl_Cairo,
                           public Sdl_Widget
   {

      protected:

         void
         blit (const Sdl_Cairo& cairo) const;

         Point_2D
         get_point (const Integer screen_i,
                    const Integer screen_j) const;

         Point_2D
         get_point (const Real screen_x,
                    const Real screen_y) const;

         Point_2D
         get_point (const Point_2D& screen_point) const;

         virtual void
         main ();

      public:

         Sdl_Application (const Size_2D& size_2d,
                          const uint32_t sdl_flags);

   };

}

/*
   class Sdl_Window : public Sdl_Plotter
   {

      public:

         Sdl_Window (const Size_2D& size_2d,
                     bool full_screen = false,
                     bool resizable = false);

         ~Sdl_Window ();

         void
         blit (const Sdl_Plotter& sdl_plotter,
               const Index_2D& index = Index_2D (0, 0));

         virtual void
         update ();

         virtual void
         main_loop ();

   };

}
*/

#endif /* DENISE_SDL_H */

