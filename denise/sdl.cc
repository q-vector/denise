//
// sdl.cc
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

#include "sdl.h"

using namespace denise;

void
Sdl_Cairo::init ()
{

   const Size_2D& size_2d = get_size_2d ();

   unsigned char* buffer = (unsigned char*)(sdl_surface_ptr->pixels);
   cs_ptr = cairo_image_surface_create_for_data (buffer,
      CAIRO_FORMAT_ARGB32, size_2d.i, size_2d.j, size_2d.i*4);

   cr = cairo_create (cs_ptr);
   cairo_set_fill_rule (cr, CAIRO_FILL_RULE_EVEN_ODD);

}

Sdl_Cairo::Sdl_Cairo (const Box_2D& box_2d)
      : Raster_Cairo (box_2d, false)
{

   const Size_2D& size_2d = get_size_2d ();

   uint32_t sdl_flags = SDL_SWSURFACE;
   sdl_surface_ptr = SDL_CreateRGBSurface (sdl_flags, size_2d.i,
      size_2d.j, 32, 0x00ff0000, 0x0000ff00, 0x000000ff, 0xff000000);

   init ();

}

Sdl_Cairo::Sdl_Cairo (const Box_2D& box_2d,
                      const uint32_t sdl_flags)
      : Raster_Cairo (box_2d, false)
{

   const Size_2D& size_2d = get_size_2d ();

   SDL_Init (SDL_INIT_VIDEO);
   sdl_surface_ptr = SDL_SetVideoMode (size_2d.i, size_2d.j, 32, sdl_flags);

   init ();

}

Sdl_Cairo::~Sdl_Cairo ()
{
   SDL_FreeSurface (sdl_surface_ptr);
}

SDL_Surface*
Sdl_Cairo::get_sdl_surface_ptr () const
{
   return sdl_surface_ptr;
}

bool
Sdl_Widget::process_key_down (const SDL_Event& event)
{
   return false;
}

bool
Sdl_Widget::process_key_up (const SDL_Event& event)
{
   return false;
}

bool
Sdl_Widget::process_mouse_button_down (const SDL_Event& event)
{
   return false;
}

bool
Sdl_Widget::process_mouse_button_up (const SDL_Event& event)
{
   return false;
}

bool
Sdl_Widget::process_mouse_motion (const SDL_Event& event)
{
   return false;
}

bool
Sdl_Widget::process_event (const SDL_Event& event)
{

   switch (event.type)
   {

      default:
         return false;

      case SDL_KEYDOWN:
         return process_key_down (event);

      case SDL_KEYUP:
         return process_key_up (event);

      case SDL_MOUSEBUTTONDOWN:
         return process_mouse_button_down (event);

      case SDL_MOUSEBUTTONUP:
         return process_mouse_button_up (event);

      case SDL_MOUSEMOTION:
         return process_mouse_motion (event);

   }

}

bool
Sdl_Widget::handle_event (SDL_Event event)
{

   bool processed = false;

   for (list<Sdl_Widget*>::iterator iterator = widget_ptr_list.begin ();
        iterator != widget_ptr_list.end (); iterator++)
   {
      Sdl_Widget& widget = **(iterator);
      if (processed = widget.handle_event (event)) { break; }
   }

   if (!processed)
   {

      switch (event.type)
      {

         case SDL_MOUSEMOTION:
         {
            Integer i = event.motion.x - origin.x;
            Integer j = event.motion.y - origin.y;
            if (!contains (i, j)) { return false; }
            break;
         }

         case SDL_MOUSEBUTTONUP:
         case SDL_MOUSEBUTTONDOWN:
         {
            Integer i = event.button.x - origin.x;
            Integer j = event.button.y - origin.y;
            if (!contains (Point_2D (i, j))) { return false; }
            break;
         }

         case SDL_KEYUP:
         case SDL_KEYDOWN:
         {
            if (!focus) { return false; }
            break;
         }

      }

      return process_event (event);

   }

}

Sdl_Widget::Sdl_Widget (const Point_2D& origin)
   : Polygon (),
     origin (origin)
{
}

Sdl_Widget::Sdl_Widget (const Polygon& polygon,
                        const Point_2D& origin)
   : Polygon (polygon),
     origin (origin)
{
}

void
Sdl_Widget::paint (const Sdl_Cairo& sdl_cairo) const
{

   for (list<Sdl_Widget*>::const_iterator iterator = widget_ptr_list.begin ();
        iterator != widget_ptr_list.end (); iterator++)
   {
      const Sdl_Widget& widget = **(iterator);
      sdl_cairo.translate (-widget.origin.x, -widget.origin.y);
      widget.paint (sdl_cairo);
      sdl_cairo.translate (widget.origin.x, widget.origin.y);
   }

   render (sdl_cairo);

}

Sdl_Label::Sdl_Label (const string& text,
                      const Point_2D& point_2d,
                      const char justify_h,
                      const char justify_v,
                      const Point_2D& offset,
                      const Real margin_x,
                      const Real margin_y,
                      const Point_1D text_angle,
                      const bool outline)
   : Label (text, point_2d, justify_h, justify_v, offset,
            margin_x, margin_y, text_angle, outline),
     Sdl_Widget ()
{
   add (get_rectangle ());
}

Sdl_Label::Sdl_Label (const string& text,
                      const Point_2D& point_2d,
                      const char justify_h,
                      const char justify_v,
                      const Real padding,
                      const Real margin_x,
                      const Real margin_y,
                      const Point_1D text_angle,
                      const bool outline)
   : Label (text, point_2d, justify_h, justify_v, padding,
            margin_x, margin_y, text_angle, outline),
     Sdl_Widget ()
{
   add (get_rectangle ());
}

Sdl_Button::Sdl_Button (const string& text,
                        const Real width,
                        const Real height)
   : text (text),
     width (width),
     height (height)
{
}

void
Sdl_Button::render_normal (const Sdl_Cairo& sdl_cairo) const
{

   sdl_cairo.set_color (Color (0, 0, 0, 0.75));
   sdl_cairo.rectangle (Point_2D (4, 0), Point_2D (width - 4, height - 4));

   sdl_cairo.set_color (Color (0.5, 0.5, 0.5));
   sdl_cairo.rectangle (Point_2D (0, 4), Point_2D (width - 4, height - 4));

   sdl_cairo.text (text, Point_2D ((width-4) / 2, (height+4) / 2), 'c', 'c');

}

void
Sdl_Button::render_depressed (const Sdl_Cairo& sdl_cairo) const
{

   sdl_cairo.set_color (Color (0.5, 0.5, 0.5));
   sdl_cairo.rectangle (Point_2D (4, 0), Point_2D (width - 4, height - 4));

   sdl_cairo.text (text, Point_2D ((width+4) / 2, (height-4) / 2), 'c', 'c');

}

void
Sdl_Button::render (const Sdl_Cairo& sdl_cairo) const
{
}

void
Sdl_Application::blit (const Sdl_Cairo& cairo) const
{

   SDL_Rect sdl_rect;

   const Index_2D& index_2d = cairo.get_index_2d ();
   const Index_2D& size_2d = cairo.get_size_2d ();

   sdl_rect.x = index_2d.i;
   sdl_rect.y = sdl_surface_ptr->h - index_2d.j - size_2d.j;

   SDL_BlitSurface (cairo.get_sdl_surface_ptr (), NULL,
      (SDL_Surface*)sdl_surface_ptr, &sdl_rect);

   SDL_UpdateRect ((SDL_Surface*)sdl_surface_ptr, 0, 0, 0, 0);

}

Point_2D
Sdl_Application::get_point (const Integer screen_i,
                            const Integer screen_j) const
{
   return get_point (Point_2D (Real (screen_i), Real (screen_j)));
}

Point_2D
Sdl_Application::get_point (const Real screen_x,
                            const Real screen_y) const
{
   return get_point (Point_2D (screen_x, screen_y));
}

Point_2D
Sdl_Application::get_point (const Point_2D& screen_point) const
{
   return Point_2D (screen_point.x, get_size_2d ().j - screen_point.y - 1);
}

void
Sdl_Application::main ()
{

   SDL_Event event;
   
   while (true)
   {

      SDL_WaitEvent (&event);

      switch (event.type)
      {
         case SDL_MOUSEMOTION:
            event.motion.y = get_size_2d ().j - event.button.y - 1;
            event.motion.yrel *= -1;
            break;
         case SDL_MOUSEBUTTONUP:
         case SDL_MOUSEBUTTONDOWN:
            event.button.y = get_size_2d ().j - event.button.y - 1;
            break;
      }

      handle_event (event);

   }

}

Sdl_Application::Sdl_Application (const Size_2D& size_2d,
                                  const uint32_t sdl_flags)
   : Sdl_Cairo (size_2d, sdl_flags),
     Sdl_Widget (Rectangle (Box_2D (size_2d)))
{
}

/*
void
Sdl_Plotter::lock () const
{
   SDL_LockSurface (sdl_surface_ptr);
}

void
Sdl_Plotter::unlock () const
{
   SDL_UnlockSurface (sdl_surface_ptr);
}

Sdl_Window::Sdl_Window (const Size_2D& size_2d,
                        bool full_screen,
                        bool resizable)
   : Sdl_Plotter (size_2d, false),
     Sdl_Widget (Rectangle (Box_2D (size_2d)))
{

   SDL_Init (SDL_INIT_VIDEO);

   uint32_t sdl_flags = SDL_SWSURFACE | SDL_DOUBLEBUF;
   if (full_screen) { sdl_flags |= SDL_FULLSCREEN; }
   if (resizable) { sdl_flags |= SDL_RESIZABLE; }

   sdl_surface_ptr = SDL_SetVideoMode (size_2d.i, size_2d.j, 32, sdl_flags);

   init ();

}

Sdl_Window::~Sdl_Window ()
{
   SDL_FreeSurface (sdl_surface_ptr);
}

void
Sdl_Window::blit (const Sdl_Plotter& sdl_plotter,
                  const Index_2D& index)
{
   SDL_BlitSurface ( sdl_plotter.sdl_surface_ptr, NULL, sdl_surface_ptr, NULL);
}

void
Sdl_Window::update ()
{
   SDL_Flip (sdl_surface_ptr);
}

void
Sdl_Window::main_loop ()
{

   SDL_Event event;

   while (SDL_WaitEvent (&event))
   {
      switch (event.type)
      {
         case SDL_KEYDOWN:
            if (event.key.keysym.sym == SDLK_q) { SDL_Quit (); exit (0); }
            break;
      }
   }

}

*/

