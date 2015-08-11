//
// gtkmm.cc
// 
// Copyright (C) 2010 Simon E. Ching
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

#include <gtkmm.h>
#include <gtkmm/accelkey.h>
#include <gtkmm/adjustment.h>
#include <gtkmm/calendar.h>
#include <gtkmm/spinbutton.h>
#include <gdkmm/cursor.h>
#include <gtkmm/dialog.h>
#include <gtkmm/main.h>
#include <gtkmm/filechooserdialog.h>
#include <denise/astronomy.h>
#include <denise/geodesy.h>
#include <denise/gtkmm.h>
#include <denise/util.h>

Devent::Modifier::Modifier (const bool shift,
                            const bool control,
                            const bool alt)
   : shift (shift),
     control (control),
     alt (alt)
{
}

Devent::Modifier::Modifier (const Modifier& modifier)
   : shift (modifier.shift),
     control (modifier.control),
     alt (modifier.alt)
{
}

bool
Devent::Modifier::modified () const
{
   return (shift || control || alt);
}

Devent::Devent (const bool shift,
                const bool control,
                const bool alt)
   : type (GDK_NOTHING),
     state (0)
{
   if (shift) { state |= GDK_SHIFT_MASK; }
   if (control) { state |= GDK_CONTROL_MASK; }
   if (alt) { state |= GDK_MOD1_MASK; }
}

Devent::Devent (const GdkEventType type,
                const guint state)
   : type (type),
     state (state)
{
}

Devent::Modifier
Devent::get_modifier () const
{
   return Modifier (shift (), control (), alt ());
}

bool
Devent::shift () const
{
   return (state & GDK_SHIFT_MASK);
}

bool
Devent::control () const
{
   return (state & GDK_CONTROL_MASK);
}

bool
Devent::alt () const
{
   return (state & GDK_MOD1_MASK);
}

Dmouse_Event::Dmouse_Event (const GdkEventType type,
                            const guint state,
                            const Point_2D& point)
   : Devent (type, state),
     point (point)
{
}

Dmouse_Event::Dmouse_Event (const bool shift,
                            const bool control,
                            const bool alt,
                            const Point_2D& point)
   : Devent (shift, control, alt),
     point (point)
{
}

Dmouse_Button_Event::Dmouse_Button_Event (const GdkEventType type,
                                          const guint state,
                                          const guint button,
                                          const Point_2D& point)
   : Dmouse_Event (type, state, point),
     button (button)
{
}

Dmouse_Button_Event::Dmouse_Button_Event (const bool shift,
                                          const bool control,
                                          const bool alt,
                                          const Point_2D& point_2d,
                                          const guint button)
   : Dmouse_Event (shift, control, alt, point),
     button (button)
{
}

Dmouse_Scroll_Event::Dmouse_Scroll_Event (const GdkEventType type,
                                          const guint state,
                                          const GdkScrollDirection direction,
                                          const Point_2D& point)
   : Dmouse_Event (type, state, point),
     direction (direction)
{
}

Dmouse_Motion_Event::Dmouse_Motion_Event (const GdkEventType type,
                                          const guint state,
                                          const Point_2D& point)
   : Dmouse_Event (type, state, point)
{
}

Dmouse_Crossing_Event::Dmouse_Crossing_Event (const GdkEventType type,
                                              const guint state,
                                              const Point_2D& point)
   : Dmouse_Event (type, state, point)
{
}

Dkey_Event::Dkey_Event (const GdkEventType type,
                        const guint state,
                        const guint value)
   : Devent (type, state),
     value (value)
{
}

Image_Buffer::Image_Buffer ()
   : ready (false)
{
}

void
Image_Buffer::initialize (const Dwidget& widget)
{

   if (image_surface == 0)
   {
      const Integer w = Integer (round (widget.get_width ()));
      const Integer h = Integer (round (widget.get_height ()));
      image_surface = denise::get_image_surface (Size_2D (w, h));
   }

   RefPtr<Context> cr = denise::get_cr (image_surface);
   transparent.cairo (cr);
   cr->set_operator (Cairo::OPERATOR_SOURCE);
   cr->paint ();

}

const RefPtr<Context>
Image_Buffer::get_cr () const
{
   return denise::get_cr (image_surface);
}

void
Image_Buffer::blit (const RefPtr<Context>& cr)
{
   if (image_surface != 0)
   {
      cr->set_source (image_surface, 0, 0);
      cr->paint ();
   }
}

void
Image_Buffer::clear ()
{
   image_surface.clear ();
   ready = false;
}

Dwidget::Dwidget (const Dcanvas& canvas)
   : canvas ((Dcanvas&)canvas),
     container_ptr (NULL),
     activated (false),
     hidable (false),
     anchor (GSL_NAN, GSL_NAN),
     width (GSL_NAN),
     height (GSL_NAN),
     preferred_width (GSL_NAN),
     preferred_height (GSL_NAN)
{
}

void
Dwidget::set_container (Dcontainer& container)
{
   this->container_ptr = &container;
}

void
Dwidget::activate ()
{
   activated = true;
}

void
Dwidget::deactivate ()
{
   activated = false;
}

void
Dwidget::set_hidable (const bool hidable)
{
   this->hidable = hidable;
}

bool
Dwidget::get_hidable () const
{
   return hidable;
}

const Point_2D&
Dwidget::get_anchor () const
{
   return anchor;
}

Real
Dwidget::get_width () const
{
   return width;
}

Real
Dwidget::get_height () const
{
   return height;
}

Size_2D
Dwidget::get_size_2d () const
{
   const Integer w = Integer (round (width));
   const Integer h = Integer (round (height));
   return Size_2D (w, h);
}

Real
Dwidget::get_font_size () const
{
   return font_size;
}

void
Dwidget::set_font_size (const Real font_size)
{
   this->font_size = font_size;
}

void
Dwidget::being_packed (const Point_2D& anchor,
                       const Real width,
                       const Real height)
{
   this->anchor = anchor;
   this->width = width;
   this->height = height;
}

void
Dwidget::set_preferred_size (const Real preferred_width,
                             const Real preferred_height)
{
   this->preferred_width = preferred_width;
   this->preferred_height = preferred_height;
}

Size_2D
Dwidget::get_preferred_size_2d () const
{
   const Integer w = Integer (round (preferred_width));
   const Integer h = Integer (round (preferred_height));
   return Size_2D (w, h);
}

Real
Dwidget::get_preferred_width () const
{
   return preferred_width;
}

Real
Dwidget::get_preferred_height () const
{
   return preferred_height;
}

bool
Dwidget::out_of_bounds (const Point_2D& point) const
{
   if (point.is_nap () || !gsl_finite (width) || !gsl_finite (height))
   {
      return true;
   }
   else
   {
      const Real x = point.x;
      const Real y = point.y;
      return (x < 0 || y < 0 || x > width || y > height);
   }
}

bool
Dwidget::on_key_pressed (const Dkey_Event& event)
{
   return false;
}

bool
Dwidget::on_key_released (const Dkey_Event& event)
{
   return false;
}

bool
Dwidget::on_mouse_button_pressed (const Dmouse_Button_Event& event)
{
   return false;
}

bool
Dwidget::on_mouse_button_released (const Dmouse_Button_Event& event)
{
   return false;
}

bool
Dwidget::on_mouse_scroll (const Dmouse_Scroll_Event& event)
{
   return false;
}

bool
Dwidget::on_mouse_motion (const Dmouse_Motion_Event& event)
{

   bool return_value = false;

   const Real this_point_x = event.point.x - anchor.x;
   const Real this_point_y = event.point.y - anchor.y;
   const Point_2D this_point (this_point_x, this_point_y);

   const bool prev_oob = out_of_bounds (prev_point);
   const bool this_oob = out_of_bounds (this_point);

   if (prev_oob && !this_oob)
   {
      const Dmouse_Crossing_Event mce (event.type, event.state, this_point);
      if (on_mouse_enter (mce)) { return_value = true; }
   }
   else
   if (!prev_oob && this_oob)
   {
      const Dmouse_Crossing_Event mce (event.type, event.state, this_point);
      if (on_mouse_exit (mce)) { return_value = true; }
   }

   prev_point = this_point;

   if (return_value)
   {
      canvas.queue_draw ();
   }

   return return_value;

}

bool
Dwidget::on_mouse_enter (const Dmouse_Crossing_Event& event)
{
   return false;
}

bool
Dwidget::on_mouse_exit (const Dmouse_Crossing_Event& event)
{
   return false;
}

void
Dwidget::highlight (const RefPtr<Context>& cr,
                    const Color& color,
                    const Real line_width) const
{

   const Rect rect (Point_2D (0, 0), width, height);

   cr->save ();
   color.cairo (cr);
   cr->set_line_width (line_width);
   cr->translate (anchor.x, anchor.y);

   rect.cairo (cr);

   cr->translate (-anchor.x, -anchor.y);
   cr->restore ();

}

void
Dwidget::set_background_ready (const bool background_ready)
{
   background_buffer.ready = background_ready;
}

void
Dwidget::blit_background_buffer (const RefPtr<Context>& cr)
{
   if (!background_buffer.ready) { render_background_buffer (); }
   background_buffer.blit (cr);
}

void
Dwidget::render_background_buffer ()
{
   background_buffer.initialize (*this);
   const RefPtr<Context> cr = background_buffer.get_cr ();
   render_background_buffer (cr);
   set_background_ready (true);
}

void
Dwidget::render_background_buffer (const RefPtr<Context>& cr)
{
}

void
Dwidget::cairo (const RefPtr<Context>& cr)
{
}

Dcontainer::Dcontainer (const Dcanvas& canvas)
   : Dwidget (canvas),
     packed (false),
     hide_hidable (false)
{
}

void
Dcontainer::register_widget (Dwidget& widget)
{
   const Integer n = widget_ptr_map.size ();
   const Integer id = (n == 0 ?  0 : widget_ptr_map.rbegin ()->first + 1);
   widget_ptr_map.insert (make_pair (id, (Dwidget*)&widget));
   widget.set_container (*this);
}

void
Dcontainer::set_hidable (const bool hidable)
{
   Dwidget::set_hidable (hidable);
   for (map<Integer, Dwidget*>::iterator iterator = widget_ptr_map.begin ();
        iterator != widget_ptr_map.end (); iterator++)
   {
      Dwidget& widget = *(iterator->second);
      widget.set_hidable (hidable);
   }
}

void
Dcontainer::set_hide_hidable (const bool hide_hidable)
{
   this->hide_hidable = hide_hidable;
}

bool
Dcontainer::get_hide_hidable () const
{
   return hide_hidable;
}

bool
Dcontainer::on_mouse_button_pressed (const Dmouse_Button_Event& event)
{

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   const Dmouse_Button_Event e (event.type, event.state, event.button, point);

   for (map<Integer, Dwidget*>::iterator iterator = widget_ptr_map.begin ();
        iterator != widget_ptr_map.end (); iterator++)
   {

      Dwidget& widget = *(iterator->second);
      if (hide_hidable && widget.get_hidable ()) { continue; }

      try
      {
         Dcontainer& container = dynamic_cast<Dcontainer&>(widget);
         if (container.on_mouse_button_pressed (e)) { return true; }
      }
      catch (const std::bad_cast& bc)
      {
         if (widget.on_mouse_button_pressed (e)) { return true; }
      }

   }

   return false;

}

void
Dcontainer::clear ()
{
   widget_ptr_map.clear ();
}

void
Dcontainer::pack ()
{
}

bool
Dcontainer::on_mouse_button_released (const Dmouse_Button_Event& event)
{

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   const Dmouse_Button_Event e (event.type, event.state, event.button, point);

   for (map<Integer, Dwidget*>::iterator iterator = widget_ptr_map.begin ();
        iterator != widget_ptr_map.end (); iterator++)
   {

      Dwidget& widget = *(iterator->second);
      if (hide_hidable && widget.get_hidable ()) { continue; }

      try
      {
         Dcontainer& container = dynamic_cast<Dcontainer&>(widget);
         if (container.on_mouse_button_released (e)) { return true; }
      }
      catch (const std::bad_cast& bc)
      {
         if (widget.on_mouse_button_released (e)) { return true; }
      }

   }

   return false;

}

bool
Dcontainer::on_mouse_scroll (const Dmouse_Scroll_Event& event)
{

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   const Dmouse_Scroll_Event e (event.type, event.state, event.direction, point);

   for (map<Integer, Dwidget*>::iterator iterator = widget_ptr_map.begin ();
        iterator != widget_ptr_map.end (); iterator++)
   {

      Dwidget& widget = *(iterator->second);
      if (hide_hidable && widget.get_hidable ()) { continue; }

      try
      {
         Dcontainer& container = dynamic_cast<Dcontainer&>(widget);
         if (container.on_mouse_scroll (e)) { return true; }
      }
      catch (const std::bad_cast& bc)
      {
         if (widget.on_mouse_scroll (e)) { return true; }
      }

   }

   return false;

}

bool
Dcontainer::on_mouse_motion (const Dmouse_Motion_Event& event)
{

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   const Dmouse_Motion_Event e (event.type, event.state, point);

   bool return_value = false;

   for (map<Integer, Dwidget*>::iterator iterator = widget_ptr_map.begin ();
        iterator != widget_ptr_map.end (); iterator++)
   {

      Dwidget& widget = *(iterator->second);
      if (hide_hidable && widget.get_hidable ()) { continue; }

      try
      {
         Dcontainer& container = dynamic_cast<Dcontainer&>(widget);
         if (container.on_mouse_motion (e))
         {
            return_value = true;
            break;
         }
      }
      catch (const std::bad_cast& bc)
      {
         if (widget.on_mouse_motion (e))
         {
            return_value = true;
            break;
         }
      }

   }

   return return_value;

}

void
Dcontainer::render_widgets (const RefPtr<Context>& cr)
{

   if (anchor.is_nap ()) { return; }

   cr->save ();
   cr->translate (anchor.x, anchor.y);

   cr->save ();
   for (map<Integer, Dwidget*>::const_iterator i = widget_ptr_map.begin ();
        i != widget_ptr_map.end (); i++)
   {
      const Integer id = (i->first);
      Dwidget& widget = *(i->second);
      if (hide_hidable && widget.get_hidable ()) { continue; }
      widget.cairo (cr);
   }

   cr->restore ();

   cr->translate (-anchor.x, -anchor.y);
   cr->restore ();

}

void
Dcontainer::cairo (const RefPtr<Context>& cr)
{
   render_widgets (cr);
}

void
Dcontainer::refresh (const RefPtr<Context>& cr,
                     const Point_2D& point,
                     const Real width,
                     const Real height)
{
   if (container_ptr != NULL)
   {
      Dcontainer& container = *container_ptr;
      container.refresh (cr, anchor + point, width, height);
   }
}

Real
Dpack_Box::get_widget_width_h () const
{

   const Integer n = widget_ptr_map.size ();
   const Integer ns = separator_set.size ();
   const Integer nn = n + ns;

   const Real margins = 2 * margin;
   const Real spacings = (nn - 1) * spacing;
   const Real separators = ns * separator;

   return (width - separators - margins - spacings) / n;

}

Real
Dpack_Box::get_widget_width_v () const
{
   return (width - 2 * margin);
}

Real
Dpack_Box::get_widget_height_h () const
{
   return (height - 2 * margin);
}

Real
Dpack_Box::get_widget_height_v () const
{

   const Integer n = widget_ptr_map.size ();
   const Integer ns = separator_set.size ();
   const Integer nn = n + ns;

   const Real margins = 2 * margin;
   const Real spacings = (nn - 1) * spacing;
   const Real separators = ns * separator;

   return (height - separators - margins - spacings) / n;

}

Real
Dpack_Box::get_preferred_width_h () const
{

   Real preferred_width = 2 * margin;

   for (map<Integer, Dwidget*>::const_iterator i = widget_ptr_map.begin ();
        i != widget_ptr_map.end (); i++)
   {
      const Dwidget& widget = *(i->second);
      const Real w = widget.get_preferred_width ();
      preferred_width += w;
   }

   return preferred_width;

}

Real
Dpack_Box::get_preferred_width_v () const
{

   Real max_width = 0;

   for (map<Integer, Dwidget*>::const_iterator i = widget_ptr_map.begin ();
        i != widget_ptr_map.end (); i++)
   {
      const Dwidget& widget = *(i->second);
      const Real w = widget.get_preferred_width ();
      if (w > max_width) { max_width = w; }
   }

   return max_width + 2 * margin;

}

Real
Dpack_Box::get_preferred_height_h () const
{

   Real max_height = 0;

   for (map<Integer, Dwidget*>::const_iterator i = widget_ptr_map.begin ();
        i != widget_ptr_map.end (); i++)
   {
      const Dwidget& widget = *(i->second);
      const Real h = widget.get_preferred_height ();
      if (h > max_height) { max_height = h; }
   }

   return max_height + 2 * margin;

}

Real
Dpack_Box::get_preferred_height_v () const
{

   Real preferred_height = 2 * margin;

   for (map<Integer, Dwidget*>::const_iterator i = widget_ptr_map.begin ();
        i != widget_ptr_map.end (); i++)
   {
      const Dwidget& widget = *(i->second);
      const Real h = widget.get_preferred_height ();
      preferred_height += h;
   }

   return preferred_height;

}

void
Dpack_Box::pack_h ()
{

   const Integer n = widget_ptr_map.size ();
   const Real widget_w = get_widget_width ();
   const Real widget_h = get_widget_height ();

   Point_2D anchor (margin, margin);

   for (map<Integer, Dwidget*>::iterator iterator = widget_ptr_map.begin ();
        iterator != widget_ptr_map.end (); iterator++)
   {

      const Integer i = distance (widget_ptr_map.begin (), iterator);
      if (separator_set.find (i) != separator_set.end ())
      {
         anchor.x += separator + spacing;
      }

      Dwidget& widget = *(iterator->second);
      widget.being_packed (anchor, widget_w, widget_h);

      try
      {
         Dcontainer& container = dynamic_cast<Dcontainer&>(widget);
         container.pack ();
      }
      catch (const std::bad_cast& bc) { }
      anchor.x += widget_w + spacing;

   }

   packed = true;

}

void
Dpack_Box::pack_v ()
{

   const Integer n = widget_ptr_map.size ();
   const Real widget_w = get_widget_width ();
   const Real widget_h = get_widget_height ();

   Point_2D anchor (margin, margin);

   for (map<Integer, Dwidget*>::iterator iterator = widget_ptr_map.begin ();
        iterator != widget_ptr_map.end (); iterator++)
   {

      const Integer i = distance (widget_ptr_map.begin (), iterator);
      if (separator_set.find (i) != separator_set.end ())
      {
         anchor.y += separator + spacing;
      }

      Dwidget& widget = *(iterator->second);
      widget.being_packed (anchor, widget_w, widget_h);

      try
      {
         Dcontainer& container = dynamic_cast<Dcontainer&>(widget);
         container.pack ();
      }
      catch (const std::bad_cast& bc) { }

      anchor.y += widget_h + spacing;

   }

   packed = true;

}

Dpack_Box::Dpack_Box (const Dcanvas& canvas,
                      const Real margin,
                      const Real spacing,
                      const bool horizontal)
   : Dcontainer (canvas),
     margin (margin),
     spacing (spacing),
     separator (margin),
     horizontal (horizontal)
{
}

void
Dpack_Box::pack_front (const Dwidget& widget)
{
   const Integer n = widget_ptr_map.size ();
   const Integer id = (n == 0 ? 0 : widget_ptr_map.begin ()->first - 1);
   Dwidget* widget_ptr = (Dwidget*)(&widget);
   widget_ptr_map.insert (make_pair (id, widget_ptr));
}

void
Dpack_Box::pack_back (const Dwidget& widget)
{
   const Integer n = widget_ptr_map.size ();
   const Integer id = (n == 0 ? 0 : widget_ptr_map.rbegin ()->first + 1);
   Dwidget* widget_ptr = (Dwidget*)(&widget);
   widget_ptr_map.insert (make_pair (id, widget_ptr));
}

void
Dpack_Box::add_separator (const bool front)
{
   const Integer index = (front ? 0 : widget_ptr_map.size ());
   separator_set.insert (index);
}

Real
Dpack_Box::get_widget_width () const
{
   return (horizontal ? get_widget_width_h () : get_widget_width_v ());
}

Real
Dpack_Box::get_widget_height () const
{
   return (horizontal ? get_widget_height_h () : get_widget_height_v ());
}

Real
Dpack_Box::get_preferred_width () const
{
   return (horizontal ? get_preferred_width_h () : get_preferred_width_v ());
}

Real
Dpack_Box::get_preferred_height () const
{
   return (horizontal ? get_preferred_height_h () : get_preferred_height_v ());
}

void
Dpack_Box::pack ()
{
   if (horizontal) { pack_h (); } else { pack_v (); }
}

Dh_Pack_Box::Dh_Pack_Box (const Dcanvas& canvas,
                          const Real margin,
                          const Real spacing)
   : Dpack_Box (canvas, margin, spacing, true)
{
}

Dv_Pack_Box::Dv_Pack_Box (const Dcanvas& canvas,
                          const Real margin,
                          const Real spacing)
   : Dpack_Box (canvas, margin, spacing, false)
{
}

Dgrid_Box::Dgrid_Box (const Dcanvas& canvas,
                      const Real margin,
                      const Real spacing,
                      const Size_2D& size_2d)
   : Dcontainer (canvas),
     margin (margin),
     spacing (spacing),
     size_2d (size_2d)
{
}

void
Dgrid_Box::pack (const Dwidget& widget,
                 const Index_2D& index_2d)
{

   // Do nothing if the slot is already occupied
   typedef map<Index_2D, Integer>::const_iterator Iterator;
   Iterator iterator = widget_id_map.find (index_2d);
   if (iterator != widget_id_map.end ()) { return; }

   const Integer n = widget_ptr_map.size ();
   const Integer id = (n == 0 ? 0 : widget_ptr_map.begin ()->first - 1);
   Dwidget* widget_ptr = (Dwidget*)(&widget);
   widget_ptr_map.insert (make_pair (id, widget_ptr));

   widget_id_map.insert (make_pair (index_2d, id));

   if (index_2d.i >= size_2d.i) { size_2d.i = index_2d.i + 1; }
   if (index_2d.j >= size_2d.j) { size_2d.j = index_2d.j + 1; }

}

Real
Dgrid_Box::get_widget_width () const
{

   const Integer n = size_2d.i;

   const Real margins = 2 * margin;
   const Real spacings = (n - 1) * spacing;

   return (width - margins - spacings) / n;

}

Real
Dgrid_Box::get_widget_height () const
{

   const Integer n = size_2d.j;

   const Real margins = 2 * margin;
   const Real spacings = (n - 1) * spacing;

   return (height - margins - spacings) / n;

}

Real
Dgrid_Box::get_preferred_width () const
{

   Index_2D index_2d;
   Real* widths = new Real[size_2d.i];
   for (index_2d.i = 0; index_2d.i < size_2d.i; index_2d.i++)
   {
      widths[index_2d.i] = 0;
   }

   typedef map<Index_2D, Integer>::const_iterator Iterator;

   for (index_2d.i = 0; index_2d.i < size_2d.i; index_2d.i++)
   {
      for (index_2d.j = 0; index_2d.j < size_2d.j; index_2d.j++)
      {
         Iterator iterator = widget_id_map.find (index_2d);
         if (iterator == widget_id_map.end ()) { continue; }
         const Integer id = iterator->second;
         const Dwidget& widget = *((widget_ptr_map.find (id))->second);
         const Real w = widget.get_preferred_width ();
         if (w > widths[index_2d.i]) { widths[index_2d.i] = w; }
      }
   }
   
   Real preferred_width = margin + margin;
   for (Integer i = 0; i < size_2d.i; i++)
   {
      preferred_width += widths[i];
      if (i != 0) { preferred_width += spacing; }
   }

   delete[] widths;
   return preferred_width;

}

Real
Dgrid_Box::get_preferred_height () const
{

   Index_2D index_2d;
   Real* heights = new Real[size_2d.j];
   for (index_2d.j = 0; index_2d.j < size_2d.j; index_2d.j++)
   {
      heights[index_2d.j] = 0;
   }

   typedef map<Index_2D, Integer>::const_iterator Iterator;

   for (index_2d.i = 0; index_2d.i < size_2d.i; index_2d.i++)
   {
      for (index_2d.j = 0; index_2d.j < size_2d.j; index_2d.j++)
      {
         Iterator iterator = widget_id_map.find (index_2d);
         if (iterator == widget_id_map.end ()) { continue; }
         const Integer id = iterator->second;
         const Dwidget& widget = *((widget_ptr_map.find (id))->second);
         const Real h = widget.get_preferred_height ();
         if (h > heights[index_2d.j]) { heights[index_2d.j] = h; }
      }
   }
   
   Real preferred_height = margin + margin;
   for (Integer j = 0; j < size_2d.j; j++)
   {
      preferred_height += heights[j];
      if (j != 0) { preferred_height += spacing; }
   }

   delete[] heights;
   return preferred_height;

}

void
Dgrid_Box::pack ()
{

   const Real widget_w = get_widget_width ();
   const Real widget_h = get_widget_height ();
   typedef map<Index_2D, Integer>::const_iterator Iterator;

   for (map<Index_2D, Integer>::iterator iterator = widget_id_map.begin ();
        iterator != widget_id_map.end (); iterator++)
   {

      const Index_2D& index_2d = iterator->first;
      if (widget_id_map.find (index_2d) == widget_id_map.end ()) { continue; }

      const Integer id = iterator->second;
      const Integer i = index_2d.i;
      const Integer j = index_2d.j;
      Dwidget& widget = *(widget_ptr_map[id]);

      const Real x = margin + i * (widget_w + spacing);
      const Real y = margin + j * (widget_h + spacing);
      const Point_2D anchor (x, y);
      widget.being_packed (anchor, widget_w, widget_h);

      // If widget is a container, pack the widget / container
      try
      {
         Dcontainer& container = dynamic_cast<Dcontainer&>(widget);
         container.pack ();
      }
      catch (const std::bad_cast& bc) { }

   }

   packed = true;

}

Real
Widget_Panel::get_widget_width () const
{
   const Real w = width;
   const Real s = hide_button_size;
   const Integer n = widget_ptr_map.size ();
   return (w - (2 * (margin + s)) - (n + 1) * spacing) / n;
}

Widget_Panel::Widget_Panel (const Dcanvas& canvas,
                            const Real margin,
                            const Real spacing,
                            const Real hide_button_size)
   : Dh_Pack_Box (canvas, margin, spacing),
     hide_button_size (hide_button_size)
{
}

void
Widget_Panel::cairo (const RefPtr<Context>& cr)
{
   Dcontainer::cairo (cr);
}

Progress::Progress (const Real fraction,
                    const Dstring& description)
   : fraction (fraction),
     description (description)
{
}

Dcanvas::Title::Title (Dcanvas& dcanvas,
                       const Color& bg_color,
                       const Color& fg_color,
                       const Color& shadow_color)
   : denise::Title (dcanvas.get_size_2d (), bg_color, fg_color, shadow_color),
     dcanvas (dcanvas)
{
}

void
Dcanvas::Title::cairo (const RefPtr<Context>& cr)
{
   set_size_2d (dcanvas.get_size_2d ());
   denise::Title::cairo (cr);
}

bool     
Dcanvas::on_configure_event (GdkEventConfigure* event)
{

   const Size_2D size_2d (event->width, event->height);
   if (get_size_2d () == size_2d) { return false; }

   this->width = size_2d.i;
   this->height = size_2d.j;

   background_buffer.clear ();
   image_buffer.clear ();
   foreground_buffer.clear ();
   initialize ();
   render_queue_draw ();

   return true;

}

bool
Dcanvas::on_draw (const RefPtr<Context>& cr)
{

   //get_window ()->set_cursor (Gdk::Cursor (Gdk::CROSSHAIR));

   if (!is_initialized ())
   {
      initialize ();
      render ();
   }

   //const GdkRectangle& area = event->area;
   //refresh (area.x, area.y, area.width, area.height);
   refresh (cr);
   return true;

}

void
Dcanvas::cairo (const RefPtr<Context>& cr)
{
   blit_background_buffer (cr);
   blit_image_buffer (cr);
   blit_foreground_buffer (cr);
   Dcontainer::cairo (cr);
}

Dcanvas::Dcanvas (Gtk::Window& gtk_window)
   : Dcontainer (*this),
     gtk_window (gtk_window),
     title (*this),
     easter_egg ("")
{
   Glib::Mutex::Lock lock (mutex);
}

const RefPtr<Context>
Dcanvas::get_widget_cr () const
{
   return denise::get_cr (widget_surface);
}

const RefPtr<Context>
Dcanvas::get_image_surface_cr () const
{
   return denise::get_cr (image_surface);
}

bool
Dcanvas::is_initialized () const
{
   return (image_surface && widget_surface);
}

void
Dcanvas::initialize ()
{

   const Integer w = Integer (round (width));
   const Integer h = Integer (round (height));
   const Size_2D size_2d (w, h);

   image_surface.clear ();
   widget_surface.clear ();

   image_surface = denise::get_image_surface (size_2d);
   widget_surface = denise::get_image_surface (size_2d);

}

void
Dcanvas::save_image (const Dstring& file_path)
{

   const Size_2D& size_2d = get_size_2d ();
   RefPtr<Surface> surface = denise::get_surface (size_2d);
   const RefPtr<Context> cr = denise::get_cr (surface);
   const string fp (file_path.begin (), file_path.end ());

   cr->set_source (image_surface, 0, 0);
   cr->paint ();

   cairo (cr);
   surface->write_to_png (fp);

}

bool
Dcanvas::on_key_pressed (const Dkey_Event& event)
{

   switch (event.value)
   {

      case GDK_KEY_space:
      {
         easter_egg += L' ';
         break;
      }

      case GDK_KEY_0:
      case GDK_KEY_1:
      case GDK_KEY_2:
      case GDK_KEY_3:
      case GDK_KEY_4:
      case GDK_KEY_5:
      case GDK_KEY_6:
      case GDK_KEY_7:
      case GDK_KEY_8:
      case GDK_KEY_9:
      {
         easter_egg += L'0' + (event.value - GDK_KEY_0);
         break;
      }

      case GDK_KEY_a:
      case GDK_KEY_b:
      case GDK_KEY_c:
      case GDK_KEY_d:
      case GDK_KEY_e:
      case GDK_KEY_f:
      case GDK_KEY_g:
      case GDK_KEY_h:
      case GDK_KEY_i:
      case GDK_KEY_j:
      case GDK_KEY_k:
      case GDK_KEY_l:
      case GDK_KEY_m:
      case GDK_KEY_n:
      case GDK_KEY_o:
      case GDK_KEY_p:
      case GDK_KEY_q:
      case GDK_KEY_r:
      case GDK_KEY_s:
      case GDK_KEY_t:
      case GDK_KEY_u:
      case GDK_KEY_v:
      case GDK_KEY_w:
      case GDK_KEY_x:
      case GDK_KEY_y:
      case GDK_KEY_z:
      {
         easter_egg += L'a' + (event.value - GDK_KEY_a);
         break;
      }

      case GDK_KEY_A:
      case GDK_KEY_B:
      case GDK_KEY_C:
      case GDK_KEY_D:
      case GDK_KEY_E:
      case GDK_KEY_F:
      case GDK_KEY_G:
      case GDK_KEY_H:
      case GDK_KEY_I:
      case GDK_KEY_J:
      case GDK_KEY_K:
      case GDK_KEY_L:
      case GDK_KEY_M:
      case GDK_KEY_N:
      case GDK_KEY_O:
      case GDK_KEY_P:
      case GDK_KEY_Q:
      case GDK_KEY_R:
      case GDK_KEY_S:
      case GDK_KEY_T:
      case GDK_KEY_U:
      case GDK_KEY_V:
      case GDK_KEY_W:
      case GDK_KEY_X:
      case GDK_KEY_Y:
      case GDK_KEY_Z:
      {
         easter_egg += L'A' + (event.value - GDK_KEY_A);
         break;
      }

      case GDK_KEY_F9:
      {

         Gtk::FileChooserDialog dialog ("Save Image...",
            Gtk::FILE_CHOOSER_ACTION_SAVE);
         dialog.set_transient_for (gtk_window);

         dialog.add_button ("Cancel", Gtk::RESPONSE_CANCEL);
         dialog.add_button ("Save", Gtk::RESPONSE_OK);

         //Glib::RefPtr<Gtk::FileFilter> filter_png = Gtk::FileFilter::create ();
         //filter_png->set_name ("PNG files");
         //filter_png->add_mime_type ("image/png");
         //dialog.add_filter (filter_png);

         if (dialog.run () == Gtk::RESPONSE_OK)
         {
            save_image (dialog.get_filename ());
         }

         return true;
         break;

      }

   }

   return Dwidget::on_key_pressed (event);

}

bool
Dcanvas::on_key_press_event (GdkEventKey* event)
{
   const Dkey_Event e (event->type, event->state, event->keyval);
   return on_key_pressed (e);
}

bool
Dcanvas::on_key_release_event (GdkEventKey* event)
{
   const Dkey_Event e (event->type, event->state, event->keyval);
   return on_key_released (e);
}

bool
Dcanvas::on_button_press_event (GdkEventButton* event)
{
   grab_focus ();
   const GdkEventType& type = event->type;
   const guint& state = event->state;
   const guint& button = event->button;
   const Point_2D point (event->x, event->y);
   const Dmouse_Button_Event e (type, state, button, point);
   return on_mouse_button_pressed (e);
}

bool
Dcanvas::on_button_release_event (GdkEventButton* event)
{
   const GdkEventType& type = event->type;
   const guint& state = event->state;
   const guint& button = event->button;
   const Point_2D point (event->x, event->y);
   const Dmouse_Button_Event e (type, state, button, point);
   return on_mouse_button_released (e);
}

bool
Dcanvas::on_scroll_event (GdkEventScroll* event)
{
   const GdkEventType& type = event->type;
   const guint& state = event->state;
   const GdkScrollDirection& direction = event->direction;
   const Point_2D point (event->x, event->y);
   const Dmouse_Scroll_Event e (type, state, direction, point);
   return on_mouse_scroll (e);
}

bool
Dcanvas::on_motion_notify_event (GdkEventMotion* event)
{
   const Point_2D point (event->x, event->y);
   const Dmouse_Motion_Event e (event->type, event->state, point);
   return on_mouse_motion (e);
}

bool
Dcanvas::on_enter_notify (GdkEventCrossing* event)
{
   const Point_2D point (event->x, event->y);
   const Dmouse_Crossing_Event e (event->type, event->state, point);
   return Dcontainer::on_mouse_enter (e);
}

bool
Dcanvas::on_leave_notify (GdkEventCrossing* event)
{
   const Point_2D point (event->x, event->y);
   const Dmouse_Crossing_Event e (event->type, event->state, point);
   return on_mouse_exit (e);
}

bool
Dcanvas::on_delete_event (GdkEventAny* event)
{
   delete this;
   return true;
}

void
Dcanvas::set_preferred_size (const Real preferred_width,
                             const Real preferred_height)
{
   const Integer pw = Integer (round (preferred_width));
   const Integer ph = Integer (round (preferred_height));
   //Gtk::DrawingArea::set_size_request (pw, ph);
   Dwidget::set_preferred_size (preferred_width, preferred_height);
}

void
Dcanvas::refresh (const RefPtr<Context>& cr)
{

   //gdk_threads_enter ();

   if (image_surface == 0) { return; }
   if (widget_surface == 0) { return; }

   const RefPtr<Context> widget_cr = get_widget_cr ();
   widget_cr->set_source (image_surface, 0, 0);
   widget_cr->paint ();

   cairo (widget_cr);

   cr->set_source (widget_surface, 0, 0);
   cr->paint ();

   if (has_focus ())
   {
      cr->save ();
      Color (0, 1, 0).cairo (cr);
      Rect (Point_2D (0, 0), 20, 2).cairo (cr);
      cr->fill ();
      cr->restore ();
   }

   //gdk_threads_leave ();

}

void
Dcanvas::refresh (const RefPtr<Context>& cr,
                  const Point_2D& point,
                  const Real width,
                  const Real height)
{

   //gdk_threads_enter ();

   if (image_surface == 0) { return; }
   if (widget_surface == 0) { return; }

   const Rect rect (point, width, height);

   const RefPtr<Context> widget_cr = get_widget_cr ();
   widget_cr->set_source (image_surface, 0, 0);
   rect.cairo (widget_cr);
   widget_cr->fill ();

   cairo (widget_cr);

   cr->set_source (widget_surface, 0, 0);
   rect.cairo (cr);
   cr->fill ();

   //gdk_threads_leave ();

}

void
Dcanvas::render ()
{

   if (!packed) { pack (); }

   if (image_surface != 0)
   {
      const RefPtr<Context>& cr = Context::create (image_surface);
      cr->select_font_face ("Verdana", FONT_SLANT_NORMAL, FONT_WEIGHT_NORMAL);
      cr->set_line_cap (LINE_CAP_ROUND);
      cr->set_line_join (LINE_JOIN_ROUND);
      cairo (cr);
   }

}

void
Dcanvas::render_queue_draw ()
{
   set_image_ready (false);
   render ();
   queue_draw ();
}

void
Dcanvas::set_image_ready (const bool ready)
{
   image_buffer.ready = ready;
}

void
Dcanvas::blit_image_buffer (const RefPtr<Context>& cr)
{
   if (!image_buffer.ready) { render_image_buffer (); }
   image_buffer.blit (cr);
}

void
Dcanvas::render_image_buffer ()
{
   image_buffer.initialize (*this);
   const RefPtr<Context> cr = image_buffer.get_cr ();
   render_image_buffer (cr);
   set_image_ready (true);
}

void
Dcanvas::render_image_buffer (const RefPtr<Context>& cr)
{
}

void
Dcanvas::set_foreground_ready (const bool ready)
{
   foreground_buffer.ready = ready;
}

void
Dcanvas::blit_foreground_buffer (const RefPtr<Context>& cr)
{
   if (!foreground_buffer.ready) { render_foreground_buffer (); }
   foreground_buffer.blit (cr);
}

void
Dcanvas::render_foreground_buffer ()
{
   foreground_buffer.initialize (*this);
   const RefPtr<Context> cr = foreground_buffer.get_cr ();
   render_foreground_buffer (cr);
   set_foreground_ready (true);
}

void
Dcanvas::render_foreground_buffer (const RefPtr<Context>& cr)
{
   title.cairo (cr);
}

Color
Dbutton::get_bg_color () const
{

   if (!disabled)
   {
      switch (state)
      {
         case BUTTON_OFF:     return Color (1.0, 1.0, 1.0, 0.9);
         case BUTTON_PRESSED: return Color (0.7, 0.7, 0.7, 0.9);
         case BUTTON_ON:      return Color (0.85, 0.85, 0.85, 0.9);
      }
   }

   if (disabled) { return Color (1.0, 1.0, 1.0, 0.9); }

}

Color
Dbutton::get_fg_color () const
{
   const Real alpha = (disabled ? 0.2 : 0.9);
   return Color (0.0, 0.0, 0.0, alpha);
}

void
Dbutton::render_background (const RefPtr<Context>& cr) const
{
   const Color& bg_color = get_bg_color ();
   const Rect rect (Point_2D (0, 0), width, height);
   bg_color.cairo (cr);
   rect.cairo (cr);
   cr->fill ();
}

void
Dbutton::render_band (const RefPtr<Context>& cr) const
{

   if (band_color.is_nac ()) { return; }

   cr->save ();
   band_color.cairo (cr);
   Rect (Point_2D (0, 0), min (font_size, width), height).cairo (cr);
   cr->fill ();
   cr->restore ();

}

void
Dbutton::render_context (const RefPtr<Context>& cr) const
{

   const Color& fg_color = get_fg_color ();
   const Rect rect (Point_2D (0, 0), width, height);
   const Point_2D center (width/2, height/2);

   fg_color.cairo (cr);
   rect.cairo (cr);
   cr->stroke ();

   const Tokens tokens (get_str (), separator);
   const Real row_height = font_size * 1.15;
   const Integer n = tokens.size ();
   const Integer n_2 = n / 2;
   const bool odd = (n % 2 == 1);

   for (Integer i = 0; i < tokens.size (); i++)
   {
      const Dstring& token = tokens[i];
      const Real dy = (i - n_2 + (odd ? 0 : 0.5)) * row_height;
      const Point_2D p (center.x, center.y + dy);
      Label label (token, p, 'c', 'c');
      label.cairo (cr, true);
   }

   symbol.cairo (cr, center);
   cr->fill ();

}

void
Dbutton::render_led (const RefPtr<Context>& cr) const
{

   if (led_color.is_nac ()) { return; }

   const Real r_0 = height / 10;
   const Real r_1 = height / 6;
   cr->save ();
   led_color.cairo (cr);
   Ring (r_0).cairo (cr, Point_2D (width - r_1, height - r_1));
   cr->fill ();
   cr->restore ();

}

Dbutton::Dbutton (const Dcanvas& canvas,
                  const Dstring& str,
                  const Real font_size,
                  const Dstring& separator)
   : Dwidget (canvas),
     str (str),
     state (BUTTON_OFF),
     disabled (false),
     band_color (Color (GSL_NAN, GSL_NAN, GSL_NAN)),
     led_color (Color (GSL_NAN, GSL_NAN, GSL_NAN)),
     separator (separator)
{
   set_font_size (font_size);
}

Dbutton::Dbutton (const Dcanvas& canvas,
                  const Symbol& symbol)
   : Dwidget (canvas),
     symbol (symbol),
     state (BUTTON_OFF),
     disabled (false),
     band_color (Color (GSL_NAN, GSL_NAN, GSL_NAN)),
     led_color (Color (GSL_NAN, GSL_NAN, GSL_NAN)),
     separator ("\t\n")
{
}

Dstring
Dbutton::get_str () const
{
   return str;
}

Dbutton::Signal&
Dbutton::get_signal ()
{
   return signal;
}

Dbutton::Str_Signal&
Dbutton::get_str_signal ()
{
   return str_signal;
}

Dbutton::Full_Signal&
Dbutton::get_full_signal ()
{
   return full_signal;
}

Dbutton::Full_Str_Signal&
Dbutton::get_full_str_signal ()
{
   return full_str_signal;
}

bool
Dbutton::on_mouse_button_pressed (const Dmouse_Button_Event& event)
{

   if (disabled) { return false; }

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   if (out_of_bounds (point)) { return false; }

   state = BUTTON_PRESSED;
   activate ();
   canvas.queue_draw ();

   return true;

}

bool
Dbutton::on_mouse_button_released (const Dmouse_Button_Event& event)
{

   if (disabled) { return false; }
   if (state != BUTTON_PRESSED) { return false; }

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   if (out_of_bounds (point)) { return false; }

   state = BUTTON_OFF;
   deactivate ();
   canvas.queue_draw ();

   clicked (event);
   return true;

}

bool
Dbutton::on_mouse_enter (const Dmouse_Crossing_Event& event)
{
   if (disabled) { return false; }
   if (activated)
   {
      state = BUTTON_PRESSED;
      canvas.queue_draw ();
      return true;
   }
   return false;
}

bool
Dbutton::on_mouse_exit (const Dmouse_Crossing_Event& event)
{
   if (disabled) { return false; }
   if (state == BUTTON_PRESSED)
   {
      state = BUTTON_OFF;
      canvas.queue_draw ();
      return true;
   }
   return false;
}

void
Dbutton::clicked (const Dmouse_Button_Event& event)
{
   const Dstring& str = get_str ();
   signal.emit ();
   str_signal.emit (str);
   full_signal.emit (event);
   full_str_signal.emit (str, event);
}

void
Dbutton::cairo (const RefPtr<Context>& cr)
{

   cr->save ();
   cr->set_line_width (1);
   cr->set_font_size (font_size);
   cr->translate (anchor.x, anchor.y);

   render_background (cr);
   render_band (cr);
   render_context (cr);
   render_led (cr);

   cr->translate (-anchor.x, -anchor.y);
   cr->restore ();

}

void
Dbutton::set_band_color (const Color& band_color)
{
   this->band_color = band_color;
}

void
Dbutton::set_led_color (const Color& led_color)
{
   this->led_color = led_color;
}

void
Dbutton::enable ()
{
   this->disabled = false;
}

void
Dbutton::disable ()
{
   this->disabled = true;
}

bool
Spin_Button::on_mouse_scroll (const Dmouse_Scroll_Event& event)
{

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   if (out_of_bounds (point)) { return false; }

   bool processed = false;

   switch (event.direction)
   {

      case GDK_SCROLL_UP:
      {
         increment ();
         processed = true;
         canvas.queue_draw ();
         break;
      }

      case GDK_SCROLL_DOWN:
      {
         decrement ();
         processed = true;
         canvas.queue_draw ();
         break;
      }

   }

   if (instant_update)
   {
      if (processed)
      {
         updated (event);
      }
   }

   return true;

}

void
Spin_Button::clicked (const Dmouse_Button_Event& event)
{

   if (!instant_update) { Dbutton::clicked (event); }

   switch (event.button)
   {

      case 1:
      {
         official_iterator = iterator;
         const Color red (1, 0, 0, 1);
         const Color nac (GSL_NAN, GSL_NAN, GSL_NAN, GSL_NAN);
         const bool b = (official_iterator != iterator);
         const Color& led_color = (b ? red : nac);
         set_led_color (led_color);

         canvas.queue_draw ();
         updated (event);
         break;
      }

      case 3:
      {
         iterator = official_iterator;
         const Color red (1, 0, 0, 1);
         const Color nac (GSL_NAN, GSL_NAN, GSL_NAN, GSL_NAN);
         const bool b = (official_iterator != iterator);
         const Color& led_color = (b ? red : nac);
         set_led_color (led_color);

         canvas.queue_draw ();
         break;
      }

   }

}

Spin_Button::Spin_Button (const Dcanvas& canvas,
                          const bool instant_update,
                          const Real font_size,
                          const Dstring& separator)
   : Dbutton (canvas, "", font_size, separator),
     instant_update (instant_update)
{
   official_iterator = end ();
   iterator = official_iterator;
}

Dstring
Spin_Button::get_official_str () const
{
   if (official_iterator == end ()) { return ""; }
   return *(official_iterator);
}

Dstring
Spin_Button::get_str () const
{
   if (iterator == end ()) { return ""; }
   return *(iterator);
}

Spin_Button::Update_Signal&
Spin_Button::get_update_signal ()
{
   return update_signal;
}

Spin_Button::Update_Str_Signal&
Spin_Button::get_update_str_signal ()
{
   return update_str_signal;
}

Spin_Button::Full_Update_Signal&
Spin_Button::get_full_update_signal ()
{
   return full_update_signal;
}

Spin_Button::Full_Update_Str_Signal&
Spin_Button::get_full_update_str_signal ()
{
   return full_update_str_signal;
}

void
Spin_Button::clear ()
{
   Tokens::clear ();
   official_iterator = end ();
   iterator = official_iterator;
}

void
Spin_Button::set (const Tokens& tokens,
                  const bool try_preserve)
{

   const Dstring& str = get_str ();
   const Dstring& official_str = get_official_str ();

   add_tokens (tokens, true);

   if (try_preserve)
   {

      iterator = std::find (begin (), end (), str);
      official_iterator = std::find (begin (), end (), official_str);   

      if (iterator == end ()) { iterator = begin (); }
      if (official_iterator == end ()) { official_iterator = begin (); }

   }

}

void
Spin_Button::add_token (const Dstring& token,
                        const bool clear_first)
{
   const Tokens tokens (token);
   add_tokens (tokens, clear_first);
}

void
Spin_Button::add_tokens (const Tokens& tokens,
                         const bool clear_first)
{

   if (clear_first) { clear (); }

   Dstring default_str = "";

   for (Tokens::const_iterator i = tokens.begin (); i != tokens.end (); i++)
   {

      const Dstring& token = *(i);
      const bool is_default = (token[0] == L'*');

      const Dstring& str = token.substr ((is_default ? 1 : 0));
      Tokens::insert (end (), str);
      if (is_default) { default_str = str; }

   }

   this->official_iterator = begin ();

   if (default_str != "")
   {
      Tokens::iterator i = std::find (begin (), end (), default_str);
      if (i != end ()) { this->official_iterator = i; }
   }

   this->iterator = this->official_iterator;

}

void
Spin_Button::increment ()
{

   if (size () < 2) { return; }

   Tokens::iterator next = iterator;
   next++;
   if (next != end ()) { iterator = next; }

   if (instant_update) { official_iterator = iterator; }

   const Color red (1, 0, 0, 1);
   const Color nac (GSL_NAN, GSL_NAN, GSL_NAN, GSL_NAN);
   const bool b = (official_iterator != iterator);
   const Color& led_color = (b ? red : nac);
   set_led_color (led_color);

}

void
Spin_Button::decrement ()
{

   if (size () < 2) { return; }

   if (iterator != begin ()) { iterator--; }
   if (instant_update) { official_iterator = iterator; }

   const Color red (1, 0, 0, 1);
   const Color nac (GSL_NAN, GSL_NAN, GSL_NAN, GSL_NAN);
   const bool b = (official_iterator != iterator);
   const Color& led_color = (b ? red : nac);
   set_led_color (led_color);

}

void
Spin_Button::updated (const Devent& event)
{
   update_signal.emit ();
   update_str_signal.emit (get_official_str ());
   full_update_signal.emit (event);
   full_update_str_signal.emit (get_official_str (), event);
}

Dtoggle_Button::Dtoggle_Button (const Dcanvas& canvas,
                                const Dstring& str,
                                const Real font_size,
                                const bool switched_on)
   : Dbutton (canvas, str, font_size)
{
   set (switched_on);
}

Dtoggle_Button::Signal&
Dtoggle_Button::get_signal ()
{
   return signal;
}

Dtoggle_Button::Str_Signal&
Dtoggle_Button::get_str_signal ()
{
   return str_signal;
}

Dtoggle_Button::Full_Signal&
Dtoggle_Button::get_full_signal ()
{
   return full_signal;
}

Dtoggle_Button::Full_Str_Signal&
Dtoggle_Button::get_full_str_signal ()
{
   return full_str_signal;
}

Dbutton::Signal&
Dtoggle_Button::get_click_signal ()
{
   Dbutton& dbutton = dynamic_cast<Dbutton&>(*this);
   return dbutton.signal;
}

Dbutton::Str_Signal&
Dtoggle_Button::get_click_str_signal ()
{
   Dbutton& dbutton = dynamic_cast<Dbutton&>(*this);
   return dbutton.str_signal;
}

Dbutton::Full_Signal&
Dtoggle_Button::get_click_full_signal ()
{
   Dbutton& dbutton = dynamic_cast<Dbutton&>(*this);
   return dbutton.full_signal;
}

Dbutton::Full_Str_Signal&
Dtoggle_Button::get_click_full_str_signal ()
{
   Dbutton& dbutton = dynamic_cast<Dbutton&>(*this);
   return dbutton.full_str_signal;
}

void
Dtoggle_Button::set (const bool switched_on)
{
   this->switched_on = switched_on;
   state = (switched_on ? BUTTON_ON : BUTTON_OFF);
}

void
Dtoggle_Button::toggle ()
{
   switched_on = !(switched_on);
   state = (switched_on ? BUTTON_ON : BUTTON_OFF);
}

const bool&
Dtoggle_Button::is_switched_on () const
{
   return switched_on;
}

bool
Dtoggle_Button::on_mouse_button_pressed (const Dmouse_Button_Event& event)
{

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   if (out_of_bounds (point)) { return false; }
   //if (event.button != 1) { return false; }

   state = BUTTON_PRESSED;
   activate ();
   canvas.queue_draw ();
   return true;

}

bool
Dtoggle_Button::on_mouse_button_released (const Dmouse_Button_Event& event)
{

   if (state != BUTTON_PRESSED) { return false; }

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   if (out_of_bounds (point)) { return false; }
   //if (event.button != 1) { return false; }

   if (event.button == 1) { switched_on = !switched_on; }
   state = (switched_on ? BUTTON_ON : BUTTON_OFF);
   deactivate ();

   if (event.button == 1) { toggled (event); }
   else { clicked (event); }

   canvas.queue_draw ();
   return true;

}

bool
Dtoggle_Button::on_mouse_enter (const Dmouse_Crossing_Event& event)
{
   if (activated)
   {
      state = BUTTON_PRESSED;
      canvas.queue_draw ();
      return true;
   }
   return false;
}

bool
Dtoggle_Button::on_mouse_exit (const Dmouse_Crossing_Event& event)
{
   if (state == BUTTON_PRESSED)
   {
      state = (switched_on ? BUTTON_ON : BUTTON_OFF);
      canvas.queue_draw ();
      return true;
   }
   return false;
}

void
Dtoggle_Button::toggled (const Dmouse_Button_Event& event)
{
   signal.emit ();
   str_signal.emit (str);
   full_signal.emit (event);
   full_str_signal.emit (str, event);
}

Integer
Radio_Button::Group::add (Radio_Button& radio_button)
{
   radio_button_ptr_vector.push_back (&radio_button);
   return radio_button_ptr_vector.size () - 1;
}

void
Radio_Button::Group::set_active (const Integer index)
{

   const Integer n = radio_button_ptr_vector.size ();

   for (Integer i = 0; i < n; i++)
   {
      Radio_Button& rb = *(radio_button_ptr_vector[i]);
      rb.set (i == index);
   }

}

Radio_Button::Radio_Button (const Dcanvas& canvas,
                            const Dstring& str,
                            const Real font_size)
   : Dtoggle_Button (canvas, str, font_size, false),
     group (*this),
     group_ptr (&group),
     group_index (0)
{
}

Radio_Button::Group&
Radio_Button::get_group ()
{
   return *group_ptr;
}

void
Radio_Button::set_group (Radio_Button::Group& group)
{
   group_ptr = &group;
   group_index = group.add (*this);
}

void
Radio_Button::set_active ()
{
   group_ptr->set_active (group_index);
}

void
Radio_Button::toggled (const Dmouse_Button_Event& event)
{
   Dtoggle_Button::toggled (event);
   set_active ();
}

void
Drawer::pack ()
{

   widget_panel.clear ();

   if (switched_on)
   {
      for (auto iterator = widget_ptr_map.begin ();
           iterator != widget_ptr_map.end (); iterator++)
      {
         const Dwidget& widget = *(iterator->second);
         widget_panel.pack_back (widget);
      }

   }

   widget_panel.pack ();
   //dcanvas.queue_draw ();

}

void
Drawer::clear ()
{

   for (auto iterator = widget_ptr_map.begin ();
        iterator != widget_ptr_map.end (); iterator++)
   {
      const Integer id = iterator->first;
      if (ref_only_map.at (id)) { continue; }
      Dwidget* widget_ptr = iterator->second;
      delete widget_ptr;
   }

   ref_only_map.clear ();
   widget_ptr_map.clear ();

}

Drawer::Drawer (Dcanvas& dcanvas,
                Dcontainer& dcontainer,
                const bool open_upwards,
                const Dstring& str,
                const Real font_size)
   : Dtoggle_Button (dcanvas, str, font_size),
     dcanvas (dcanvas),
     dcontainer (dcontainer),
     open_upwards (open_upwards),
     widget_panel (dcanvas, 0, 6)
{
   dcanvas.register_widget (widget_panel);
}

Drawer::~Drawer ()
{
   clear ();
}

const Dwidget&
Drawer::get_widget (const Integer index) const
{
   return *(widget_ptr_map.at (index));
}

Dwidget&
Drawer::get_widget (const Integer index)
{
   return *(widget_ptr_map.at (index));
}

void
Drawer::set_hidable (const bool hidable)
{
   Dwidget::set_hidable (hidable);
   widget_panel.set_hidable (hidable);
}

void
Drawer::add_widget (Dwidget& widget)
{
   const Integer id = widget_ptr_map.size ();
   widget_ptr_map.insert (make_pair (id, &widget));
   ref_only_map.insert (make_pair (id, true));
}

void
Drawer::add_widget_ptr (Dwidget* widget_ptr)
{
   const Integer id = widget_ptr_map.size ();
   widget_ptr_map.insert (make_pair (id, widget_ptr));
   ref_only_map.insert (make_pair (id, false));
}

void
Drawer::being_packed (const Point_2D& anchor,
                      const Real width,
                      const Real height)
{

   Dwidget::being_packed (anchor, width, height);

   const Real margin = 6;
   const Integer n = widget_ptr_map.size ();

   const Real wp_width = width;
   const Real wp_height = n * height + (n - 1) * margin;

   Point_2D wp_anchor = dcontainer.get_anchor ();

   wp_anchor.x += anchor.x;
   wp_anchor.y += (open_upwards ? -(margin + wp_height) : (margin + height));

   widget_panel.being_packed (wp_anchor, wp_width, wp_height);
   widget_panel.pack ();

}

void
Drawer::toggled (const Dmouse_Button_Event& event)
{
   pack ();
}

void
Drawer::toggle ()
{
   Dtoggle_Button::toggle ();
   pack ();
}

void
Drawer::expand ()
{
   if (!switched_on)
   {
      Dtoggle_Button::toggle ();
      pack ();
   }
}

void
Drawer::collapse ()
{
   if (switched_on)
   {
      Dtoggle_Button::toggle ();
      pack ();
   }
}

Drawer_Panel::Drawer_Panel (Dcanvas& dcanvas,
                            const bool open_upwards,
                            const Real font_size)
   : Dh_Pack_Box (dcanvas, 0, font_size / 2),
     dcanvas (dcanvas),
     open_upwards (open_upwards),
     font_size (font_size)
{
}

Drawer_Panel::~Drawer_Panel ()
{

   typedef Button_Ptr_Map::iterator B_Iterator;
   typedef Toggle_Button_Ptr_Map::iterator Tb_Iterator;

   // for drawers as well?

   for (B_Iterator iterator = button_ptr_map.begin ();
        iterator != button_ptr_map.end (); iterator++)
   {
      delete iterator->second;
   }

   for (Tb_Iterator iterator = toggle_button_ptr_map.begin ();
        iterator != toggle_button_ptr_map.end (); iterator++)
   {
      delete iterator->second;
   }

}

Drawer_Panel::Drawer_Ptr_Map::iterator
Drawer_Panel::add_drawer_ptr (Drawer* drawer_ptr,
                              const bool back)
{
   const Dstring& drawer_str = drawer_ptr->get_str ();
   if (back) { pack_back (*drawer_ptr); } else { pack_front (*drawer_ptr); }
   return (drawer_ptr_map.insert (make_pair (drawer_str, drawer_ptr))).first;
}

Drawer_Panel::Drawer_Ptr_Map::iterator
Drawer_Panel::add_drawer (const Dstring& str,
                          const bool back)
{

   Drawer_Ptr_Map::iterator iterator = drawer_ptr_map.find (str);
   if (iterator != drawer_ptr_map.end ()) { return iterator; }

   Drawer* drawer_ptr = new Drawer (
      dcanvas, *this, open_upwards, str, font_size);
   return add_drawer_ptr (drawer_ptr, back);

}

void
Drawer_Panel::add_widget (const Dstring& drawer_str,
                          Dwidget& widget)
{

   Drawer_Ptr_Map::iterator iterator = drawer_ptr_map.find (drawer_str);
   const bool cannot_find = (iterator == drawer_ptr_map.end ());
   if (cannot_find) { iterator = add_drawer (drawer_str, true); }

   Drawer& drawer = *(iterator->second);
   drawer.add_widget (widget);

}

void
Drawer_Panel::add_widget_ptr (const Dstring& drawer_str,
                              Dwidget* widget_ptr)
{

   Drawer_Ptr_Map::iterator iterator = drawer_ptr_map.find (drawer_str);
   const bool cannot_find = (iterator == drawer_ptr_map.end ());
   if (cannot_find) { iterator = add_drawer (drawer_str, true); }

   Drawer& drawer = *(iterator->second);
   drawer.add_widget_ptr (widget_ptr);

}

void
Drawer_Panel::add_button (const Dstring& str,
                          const bool back)
{
   Dbutton* button_ptr = new Dbutton (dcanvas, str, font_size);
   if (back) { pack_back (*button_ptr); } else { pack_front (*button_ptr); }
   button_ptr_map.insert (make_pair (str, button_ptr));
}

void
Drawer_Panel::add_toggle_button (const Dstring& str,
                                 const bool switched_on,
                                 const bool back)
{
   Dtoggle_Button* button_ptr = new Dtoggle_Button (
      dcanvas, str, font_size, switched_on);
   if (back) { pack_back (*button_ptr); } else { pack_front (*button_ptr); }
   toggle_button_ptr_map.insert (make_pair (str, button_ptr));
}

const Drawer&
Drawer_Panel::get_drawer (const Dstring& str) const
{
   Drawer_Ptr_Map::const_iterator iterator = drawer_ptr_map.find (str);
   const bool cannot_find = (iterator == drawer_ptr_map.end ());
   if (cannot_find) { throw Exception ("No Such Drawer"); }
   return *(iterator->second);
}

Drawer&
Drawer_Panel::get_drawer (const Dstring& str)
{
   Drawer_Ptr_Map::iterator iterator = drawer_ptr_map.find (str);
   const bool cannot_find = (iterator == drawer_ptr_map.end ());
   if (cannot_find) { throw Exception ("No Such Drawer"); }
   return *(iterator->second);
}

const Dbutton&
Drawer_Panel::get_button (const Dstring& str) const
{
   Button_Ptr_Map::const_iterator iterator = button_ptr_map.find (str);
   const bool cannot_find = (iterator == button_ptr_map.end ());
   if (cannot_find) { throw Exception ("No Such Button"); }
   return *(iterator->second);
}

Dbutton&
Drawer_Panel::get_button (const Dstring& str)
{
   Button_Ptr_Map::iterator iterator = button_ptr_map.find (str);
   const bool cannot_find = (iterator == button_ptr_map.end ());
   if (cannot_find) { throw Exception ("No Such Button"); }
   return *(iterator->second);
}

const Dtoggle_Button&
Drawer_Panel::get_toggle_button (const Dstring& str) const
{
   Toggle_Button_Ptr_Map::const_iterator iterator = toggle_button_ptr_map.find (str);
   const bool cannot_find = (iterator == toggle_button_ptr_map.end ());
   if (cannot_find) { throw Exception ("No Such Toggle_Button"); }
   return *(iterator->second);
}

Dtoggle_Button&
Drawer_Panel::get_toggle_button (const Dstring& str)
{
   Toggle_Button_Ptr_Map::iterator iterator = toggle_button_ptr_map.find (str);
   const bool cannot_find = (iterator == toggle_button_ptr_map.end ());
   if (cannot_find) { throw Exception ("No Such Toggle_Button"); }
   return *(iterator->second);
}

bool
Drawer_Panel::is_switched_on (const Dstring& str) const
{
   const Dtoggle_Button& toggle_button = get_toggle_button (str);
   return toggle_button.is_switched_on ();
}

void
Drawer_Panel::toggle (const Dstring& str)
{
   Dtoggle_Button& toggle_button = get_toggle_button (str);
   toggle_button.toggle ();
}

void
Drawer_Panel::set_toggle (const Dstring& str,
                          const bool switched_on)
{
   Dtoggle_Button& toggle_button = get_toggle_button (str);
   toggle_button.set (switched_on);
}

bool
Drawer_Panel::all_collapsed () const
{
   for (Drawer_Ptr_Map::const_iterator iterator = drawer_ptr_map.begin ();
        iterator != drawer_ptr_map.end (); iterator++)
   {
      const Drawer& drawer = *(iterator->second);
      if (drawer.is_switched_on ()) { return false; }
   }
   return true;
}

void
Drawer_Panel::expand_all ()
{
   for (Drawer_Ptr_Map::iterator iterator = drawer_ptr_map.begin ();
        iterator != drawer_ptr_map.end (); iterator++)
   {
      Drawer& drawer = *(iterator->second);
      drawer.expand ();
   }
}

void
Drawer_Panel::collapse_all ()
{
   for (Drawer_Ptr_Map::iterator iterator = drawer_ptr_map.begin ();
        iterator != drawer_ptr_map.end (); iterator++)
   {
      Drawer& drawer = *(iterator->second);
      drawer.collapse ();
   }
}

Dtitle::Dtitle (const Dcanvas& canvas,
                const Real font_size)
   : Dwidget (canvas),
     margin (font_size / 2)
{
   set_font_size (font_size);
}

Dtitle::Dtitle (const Dcanvas& canvas,
                const Real font_size,
                const Dstring& str_l,
                const Dstring& str_c,
                const Dstring& str_r)
   : Dwidget (canvas),
     margin (font_size / 2),
     string_l (string_l),
     string_c (string_c),
     string_r (string_r)
{
   set_font_size (font_size);
}

void
Dtitle::set_string_l (const Dstring& string_l)
{
   this->string_l = string_l;
}

void
Dtitle::set_string_c (const Dstring& string_c)
{
   this->string_c = string_c;
}

void
Dtitle::set_string_r (const Dstring& string_r)
{
   this->string_r = string_r;
}

Real
Dtitle::get_preferred_width () const
{
   return 5;
}

Real
Dtitle::get_preferred_height () const
{

   FontExtents fe;
   const RefPtr<Context> cr = canvas.get_widget_cr ();

   cr->save ();
   cr->set_font_size (font_size);
   cr->get_font_extents (fe);
   cr->restore ();

   return fe.height + 2 * margin;

}

void
Dtitle::cairo (const RefPtr<Context>& cr)
{

   const Point_2D point_l (margin, height/2);
   const Point_2D point_c (width/2, height/2);
   const Point_2D point_r (width - margin, height/2);

   const Rect rect (Point_2D (0, 0), width, height);

   Label label_l (string_l, point_l, 'l', 'c');
   Label label_c (string_c, point_c, 'c', 'c');
   Label label_r (string_r, point_r, 'r', 'c');

   cr->save ();
   cr->set_line_width (2);
   cr->set_font_size (font_size);
   cr->translate (anchor.x, anchor.y);

   Color (1, 1, 1, 0.8).cairo (cr);
   rect.cairo (cr);
   cr->fill_preserve ();
   Color (0, 0, 0, 0.8).cairo (cr);
   cr->stroke ();

   label_l.cairo (cr, true);
   label_c.cairo (cr, true);
   label_r.cairo (cr, true);

   cr->translate (-anchor.x, -anchor.y);
   cr->restore ();

}

Integer
Popup::get_index (const Point_2D& point) const
{

   FontExtents fe;
   const RefPtr<Context> cr = canvas.get_widget_cr ();
   
   cr->save ();
   cr->set_font_size (font_size);
   cr->get_font_extents (fe);
   cr->restore ();

   const Real& y = point.y;
   
   const bool ne = (orientation == ORIENTATION_NE);
   const bool nw = (orientation == ORIENTATION_NW);
   const bool north = (ne || nw);
   const Real dy = (north ? y - margin + height : y - margin);
   const Integer i = Integer (floor (dy / fe.height)); 

   return i;

}

void
Popup::cairo (const RefPtr<Context>& cr,
              const Integer index,
              const FontExtents& fe) const
{

   const Real fe_h = fe.height;
   const Dstring& str = tokens[index];
   const Real dy = (index + 0.5) * fe_h;

   switch (orientation)
   {

      case ORIENTATION_NE:
      {
         const Point_2D point (fe_h / 2, dy + margin - height);
         Label label (str, point, 'l', 'c');
         label.cairo (cr, true);
         break;
      }

      case ORIENTATION_SE:
      {
         const Point_2D point (fe_h / 2, dy + margin);
         Label label (str, point, 'l', 'c');
         label.cairo (cr, true);
         break;
      }

      case ORIENTATION_NW:
      {
         const Point_2D point (-fe_h / 2, dy + margin - height);
         Label label (str, point, 'r', 'c');
         label.cairo (cr, true);
         break;
      }

      case ORIENTATION_SW:
      {
         const Point_2D point (-fe_h / 2, dy + margin);
         Label label (str, point, 'r', 'c');
         label.cairo (cr, true);
         break;
      }

   }

}


Popup::Popup (const Dcanvas& canvas,
              const Real font_size)
   : Dwidget (canvas),
     font_size (font_size),
     margin (font_size / 2),
     orientation (ORIENTATION_SE),
     index (-1)
{
}

void
Popup::hide ()
{
   this->anchor.x = GSL_NAN;
   this->anchor.y = GSL_NAN;
   canvas.queue_draw ();
}

void
Popup::clear ()
{
   tokens.clear ();
   hide ();
}

void
Popup::append (const Dstring& str)
{
   tokens.push_back (str);
}

void
Popup::set_tokens (const Tokens& tokens)
{
   this->tokens.clear ();
   for (Tokens::const_iterator iterator = tokens.begin ();
        iterator != tokens.end (); iterator++)
   {
      const Dstring& str = *(iterator);
      this->tokens.push_back (str);
   }
}

void
Popup::set_shape (const Point_2D& anchor,
                  const Orientation orientation,
                  const Real width,
                  const Real height)
{
   this->orientation = orientation;
   being_packed (anchor, width, height);
}

Real
Popup::get_preferred_width () const
{
   return 14 * font_size;
}

Real
Popup::get_preferred_height () const
{

   FontExtents fe;
   const RefPtr<Context> cr = canvas.get_widget_cr ();
   
   cr->save ();
   cr->set_font_size (font_size);
   cr->get_font_extents (fe);
   cr->restore ();
   
   return tokens.size () * fe.height + 2 * margin;

}

void
Popup::reset_index ()
{
   this->index = -1;
}

void
Popup::cairo (const RefPtr<Context>& cr)
{

   if (anchor.is_nap ()) { return; }

   cr->save ();
   cr->translate (anchor.x, anchor.y);

   switch (orientation)
   {

      case ORIENTATION_NE:
      {
         Rect (Point_2D (0, -height), width, height).cairo (cr);
         break;
      }

      case ORIENTATION_SE:
      {
         Rect (Point_2D (0, 0), width, height).cairo (cr);
         break;
      }

      case ORIENTATION_NW:
      {
         Rect (Point_2D (-width, -height), width, height).cairo (cr);
         break;
      }

      case ORIENTATION_SW:
      {
         Rect (Point_2D (-width, 0), width, height).cairo (cr);
         break;
      }

   }

   Color (1, 1, 1, 0.7).cairo (cr);
   cr->fill_preserve ();
   Color (0, 0, 0).cairo (cr);
   cr->set_line_width (font_size / 6);
   cr->stroke ();

   FontExtents fe;
   cr->set_font_size (font_size);
   cr->get_font_extents (fe);
   const Real fe_h = fe.height;

   for (Integer i = 0; i < tokens.size (); i++)
   {
      if (i == index) { Color (1, 0, 0, 1).cairo (cr); }
      else            { Color (0, 0, 0, 1).cairo (cr); }
      cairo (cr, i, fe);
   }

   cr->restore ();

}

void
Popup_Menu::emit_signal (const Dstring& str) const
{
   signal_map.at (str).emit ();
   str_signal_map.at (str).emit (tokens[index]);
}

Popup_Menu::Popup_Menu (const Dcanvas& canvas,
                        const Real font_size)
   : Popup (canvas, font_size)
{
}

const bool
Popup_Menu::is_fresh () const
{
   return fresh;
}

void
Popup_Menu::set_fresh (const bool fresh)
{
   this->fresh = fresh;
}

void
Popup_Menu::clear ()
{
   Popup::clear ();
   signal_map.clear ();
   str_signal_map.clear ();
}

void
Popup_Menu::append (const Dstring& str)
{
   Popup::append (str);
   Popup_Menu::Signal signal;
   Popup_Menu::Str_Signal str_signal;
   signal_map.insert (make_pair (str, signal));
   str_signal_map.insert (make_pair (str, str_signal));
}

Popup_Menu::Signal&
Popup_Menu::get_signal (const Dstring& str)
{
   return signal_map.at (str);
}

Popup_Menu::Str_Signal&
Popup_Menu::get_str_signal (const Dstring& str)
{
   return str_signal_map.at (str);
}

bool
Popup_Menu::is_on () const
{
   return !(gsl_isnan (anchor.x) || gsl_isnan (anchor.y));
}

bool
Popup_Menu::on_mouse_motion (const Dmouse_Motion_Event& event)
{

   //if (Gtk::Main::events_pending ()) { return true; }
   if (Dwidget::on_mouse_motion (event)) { return true; }

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   if (out_of_bounds (point)) { return false; }

   const Integer i = get_index (point);
   if (i != index)
   {
      index = i;
      canvas.queue_draw ();
   }

   return true;
}

bool
Popup_Menu::on_mouse_button_released (const Dmouse_Button_Event& event)
{

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   if (out_of_bounds (point)) { return false; }

   const Integer i = get_index (point);
   if (i >= 0 && i < tokens.size ()) { emit_signal (tokens[i]); }

   hide ();
   return true;

}

bool
Popup_Menu::on_mouse_exit (const Dmouse_Crossing_Event& event)
{
   reset_index ();
   canvas.queue_draw ();
   return true;
}

Time_Chooser::Shape::Shape ()
   : start_time (GSL_NAN),
     end_time (GSL_NAN),
     leap (24)
{
}

Time_Chooser::Shape::Shape (const Shape& shape)
   : set<Dtime> (shape),
     start_time (shape.start_time),
     end_time (shape.end_time),
     leap (shape.leap)
{

   if (shape.size () == 0) { return; }

   start_time = *(shape.begin ());
   end_time = *(shape.rbegin ());
   if (start_time == end_time) { end_time.t = start_time.t + 3; }

   //const bool short_span = (get_span_t () < 36);
   //leap = (short_span ? 3 : 24);

}

Time_Chooser::Shape::Shape (const set<Dtime>& time_set)
   : set<Dtime> (time_set),
     start_time (GSL_NAN),
     end_time (GSL_NAN),
     leap (24)
{

   if (time_set.size () == 0) { return; }

   start_time = *(time_set.begin ());
   end_time = *(time_set.rbegin ());
   if (start_time == end_time) { end_time.t = start_time.t + 3; }

   //const bool short_span = (get_span_t () < 36);
   //leap = (short_span ? 3 : 24);

}

Real
Time_Chooser::Shape::get_span_t () const
{
   return end_time.t - start_time.t;
}

bool
Time_Chooser::Shape::out_of_bounds (const Dtime& dtime) const
{
   return (dtime < start_time || dtime > end_time);
}

Time_Chooser::Shape::const_iterator
Time_Chooser::Shape::get_iterator (const Dtime& dtime) const
{

   if (size () == 0) { throw Exception ("Time_Chooser::Data empty"); }

   Real smallest_fabs_dt = GSL_POSINF;
   Shape::const_iterator iterator = end ();

   for (Shape::const_iterator i = begin (); i != end (); i++)
   {
      const Real fabs_dt = fabs (i->t - dtime.t);
      if (fabs_dt < smallest_fabs_dt)
      {
         iterator = i;
         smallest_fabs_dt = fabs_dt;
      }
   }

   return iterator;

}

const Dtime&
Time_Chooser::Shape::get_nearest_time (const Dtime& dtime) const
{
   return *(get_iterator (dtime));
}

const Dtime&
Time_Chooser::Shape::get_next_time (const Dtime& dtime,
                                    const bool is_leap) const
{

   if (is_leap) { return get_nearest_time (dtime.t + leap); }

   Shape::const_iterator iterator = get_iterator (dtime);
   iterator++;
   if (iterator == end ()) { iterator = begin (); }
   return *(iterator);

}

const Dtime&
Time_Chooser::Shape::get_prev_time (const Dtime& dtime,
                                    const bool is_leap) const
{

   if (is_leap) { return get_nearest_time (dtime.t - leap); }

   Shape::const_iterator iterator = get_iterator (dtime);
   if (iterator == begin ()) { return *(rbegin ()); }
   else { return *(--iterator); }

}

bool
Time_Chooser::Shape::operator== (const Time_Chooser::Shape& shape) const
{

   if (start_time != shape.start_time) { return false; }
   if (end_time != shape.end_time) { return false; }
   if (leap != shape.leap) { return false; }

   typedef set<Dtime> St;
   const St& t = dynamic_cast<const St&>(*this);
   const St& s = dynamic_cast<const St&>(shape);
   if (t != s) { return false; }

   return true;

}

void
Time_Chooser::Data::init (const Dstring& str)
{
   const Tokens& tokens (str);
   if (tokens.size () > 0) { dtime = Dtime (tokens[0]); }
}

void
Time_Chooser::Data::rectify (Dstring& str)
{
   const Tokens tokens (str, ":");
   const bool b = (tokens.size () > 1 && tokens[0] > tokens[1]);
   if (b) { str = tokens[1] + ":" + tokens[0]; } }

bool
Time_Chooser::Data::contains (const Dstring& str,
                              const Dtime& dtime)
{

   const Tokens tokens (str, ":");
   if (tokens.size () == 0) { return false; }
   if (tokens.size () == 1) { return dtime == Dtime (tokens[0]); }

   const Dtime start_time (tokens[0]);
   const Dtime end_time (tokens[1]);
   const Real dt_s = dtime.t - start_time.t;
   const Real dt_e = dtime.t - end_time.t;
   return (dt_s * dt_e) <= 0;

}

Time_Chooser::Data::Data ()
   : dtime (GSL_NAN),
     candidate_time (GSL_NAN)
{
}

Time_Chooser::Data::Data (const Data& data)
   : Tokens (data),
     dtime (data.dtime),
     candidate_time (data.candidate_time)
{
}

Time_Chooser::Data::Data (const Tokens& tokens)
   : Tokens (tokens),
     dtime (GSL_NAN),
     candidate_time (GSL_NAN)
{
   if (tokens.size () > 0) { init (tokens[0]); }
}

Time_Chooser::Data::Data (const Dtime& dtime)
   : Tokens (dtime.get_string ("%Y%m%d%H%M")),
     dtime (GSL_NAN),
     candidate_time (GSL_NAN)
{
   init (dtime.get_string ("%Y%m%d%H%M"));
}

Time_Chooser::Data::Data (const Dstring& str)
   : Tokens (str),
     dtime (GSL_NAN),
     candidate_time (GSL_NAN)
{
   init (str);
}

void
Time_Chooser::Data::clear ()
{
   dtime = GSL_NAN;
   candidate_time = GSL_NAN;
   Tokens::clear ();
}

void
Time_Chooser::Data::set_time (const Dtime& dtime,
                              const bool clear_first,
                              const bool is_range,
                              const bool is_leap)
{

   const Dstring fmt ("%Y%m%d%H%M");

   const Dstring& old_time_str = this->dtime.get_string (fmt);
   const Dstring& new_time_str = dtime.get_string (fmt);

   if (clear_first) { clear (); }

   if (!is_range)
   {
      const bool no_need = contains (dtime);
      const bool needy = (!clear_first || size () == 0);
      if (!no_need && needy) { push_back (new_time_str); }
   }
   else
   {
      for (Data::iterator iterator = begin ();
           iterator != end (); iterator++)
      {
         Dstring& str = *(iterator);
         const Tokens tokens (str, ":");
         const bool singleton = (tokens.size () == 1);
         if (singleton)
         {
            if (str == old_time_str)
            {
               str = old_time_str + ":" + new_time_str;
               rectify (str);
               break;
            }
         }
         else
         {
            if (!contains (str, this->dtime)) { continue; }
            const Dstring& str_a = std::min (tokens[0], new_time_str);
            const Dstring& str_b = std::max (tokens[1], new_time_str);
            str = str_a + ":" + str_b;
            rectify (str);
            break;
         }
      }
   }

   this->dtime = dtime;

}

void
Time_Chooser::Data::set_time (const Dtime& dtime,
                              const Devent& event)
{

   const bool is_range = event.shift ();
   const bool clear_first = !(event.shift () || event.control ());

   set_time (dtime, clear_first, is_range, false);

}

Time_Chooser::Shape
Time_Chooser::Data::get_highlighted_shape (const Shape& shape) const
{

   set<Dtime> time_set;

   for (Shape::const_iterator iterator = shape.begin ();
        iterator != shape.end (); iterator++)
   {
      const Dtime& dtime = *(iterator);
      if (contains (dtime)) { time_set.insert (dtime); }
   }

   Time_Chooser::Shape highlighted_shape (time_set);
   highlighted_shape.leap = shape.leap;
   return highlighted_shape;

}

bool
Time_Chooser::Data::matches (const Dtime& dtime) const
{
   return this->dtime == dtime;
}

bool
Time_Chooser::Data::is_candidate (const Dtime& dtime) const
{
   return this->candidate_time == dtime;
}

bool
Time_Chooser::Data::contains (const Dtime& dtime) const
{
   for (Data::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Dstring& str = *(iterator);
      if (contains (str, dtime)) { return true; }
   }
   return false;
}

void
Time_Chooser::Data::conform_to (const Shape& shape)
{

   Tokens new_tokens;
   const Dstring fmt ("%Y%m%d%H%M");

   for (Data::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Tokens tokens (*(iterator), ":");
      const bool singleton = (tokens.size () == 1);

      if (singleton)
      {
         const Dtime& new_time = shape.get_nearest_time (tokens[0]);
         new_tokens.push_back (new_time.get_string (fmt));
      }
      else
      {
         const Dtime& s (std::min (tokens[0], tokens[1]));
         const Dtime& e (std::max (tokens[0], tokens[1]));
         const Dtime& st = shape.start_time;
         const Dtime& et = shape.end_time;
         if (st > e || et < s) { return; }
         const Dstring& str_s = std::max (s, st).get_string (fmt);
         const Dstring& str_e = std::min (e, et).get_string (fmt);
         new_tokens.push_back (str_s + ":" + str_e);
      }

   }

   if (new_tokens.size () == 0) { set_time (*(shape.begin ())); }
   else
   {
      clear ();
      for (Tokens::const_iterator i = new_tokens.begin ();
           i != new_tokens.end (); i++)
      {
         push_back (*(i));
      }
      init (new_tokens[0]);
      candidate_time = GSL_NAN;
   }

}

void
Time_Chooser::Data::dump () const
{
   cout << "===" << endl;
   for (Data::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      cout << *(iterator) << endl;
   }
   cout << "===" << endl;
}

Real
Time_Chooser::get_s (const Dtime& dtime) const
{

   const Real start_x = button_box.index_2d.i + button_box.size_2d.i + margin;
   const Real start_y = button_box.index_2d.j + button_box.size_2d.j + margin;
   const Real span_x = width - start_x - margin;
   const Real span_y = height - start_y - margin;
   const Real span = (vertical ? span_y : span_x);
   const Real span_t = shape.get_span_t ();
   const Real start = (vertical ? start_y : start_x);
   const Real ratio = span / span_t;

   const Real sign = boolean_sign (!reverse);
   const Real dt = sign * (dtime.t - shape.start_time.t);
   return start + dt * ratio;

}

Dtime
Time_Chooser::get_time (const Point_2D& point,
                        const bool snap_to_nearest) const
{

   const Real w = 4;
   const Real centre = (vertical ? width - margin : margin);
   const Real x = (vertical ? point.x : point.y);
   if (fabs (x - centre) > w) { return Dtime (GSL_NAN); }

   const Real start_x = button_box.index_2d.i + button_box.size_2d.i + margin;
   const Real start_y = button_box.index_2d.j + button_box.size_2d.j + margin;
   const Real span_x = width - start_x - margin;
   const Real span_y = height - start_y - margin;
   const Real span = (vertical ? span_y : span_x);
   const Real span_t = shape.get_span_t ();
   const Real start = (vertical ? start_y : start_x);
   const Real ratio = span / span_t;

   const Real s = (vertical ? point.y : point.x);
   const Real dt = (s - start) / ratio;
   const Dtime dtime (shape.start_time.t + dt);

   return (snap_to_nearest ? shape.get_nearest_time (dtime) : dtime);

}

void
Time_Chooser::render_background (const RefPtr<Context>& cr)
{
   const Real span_t = shape.get_span_t ();
   const Real number_of_days = span_t / 24;
   const bool in_days = (number_of_days < 21);
   const bool in_months = (number_of_days > 50);
   if (in_months) { render_background_months (cr); }
   else { render_background_days (cr); }
}

void
Time_Chooser::render_background_months (const RefPtr<Context>& cr)
{

   cr->save ();

   const Color bg_color (1.0, 1.0, 1.0, 0.7);
   const Color fg_color (0.0, 0.0, 0.0, 1.0);
   const Color color (0.5, 0.5, 0.5, 1.0);
   const Rect rect (Point_2D (0, 0), width, height);

   bg_color.cairo (cr);
   rect.cairo (cr);
   cr->fill_preserve ();
   fg_color.cairo (cr);
   cr->stroke ();

   cr->set_font_size (font_size);

   cr->set_line_width (1);

   const Integer start_year = shape.start_time.get_year ();
   const Integer end_year = shape.end_time.get_year ();

   for (Integer year = start_year; year <= end_year; year++)
   {
      for (Integer month = 1; month <= 12; month++)
      {
         color.cairo (cr);
         render_month_lines (cr, year, month, "%d");
         fg_color.cairo (cr);
         render_month (cr, year, month);
      }
   }

   cr->restore ();

}

void
Time_Chooser::render_month (const RefPtr<Context>& cr,
                            const Integer year,
                            const Integer month) const
{

   typedef Point_2D P;
   const Dtime dtime (year, month, 1);
   if (shape.out_of_bounds (dtime)) { return; }

   cr->save ();

   const Point_2D& p = get_point (dtime);
   const Real x = p.x;
   const Real y = p.y;
   const Real w = width;
   const Real h = height;
   const Real m = margin;
   const bool v = vertical;

   const Point_2D& p_a = (v ? P (w - 6.5 * m, y) : P (x, h - 1.5 * m));
   const Point_2D& p_b = (v ? P (w - 5.0 * m, y) : P (x, h - 3.0 * m));
   const Point_2D& p_c = (v ? P (w - 4.5 * m, y) : P (x, h - 3.5 * m));
   const Point_2D& p_d = (v ? P (w - 7.0 * m, y) : P (x, h - 1.0 * m));
   const Dstring& str_a = dtime.get_string (v ? "%b" : "%B");
   const Dstring& str_b = dtime.get_string ("%Y");

   cr->set_line_width (1);
   cr->move_to (p_c.x, p_c.y);
   cr->line_to (p_d.x, p_d.y);
   cr->stroke ();

   const bool align_left = (vertical == reverse);
   const char align_char = (align_left ? 'l' : 'r');

   Label label_a (str_a, p_a, align_char, 'c', margin);
   Label label_b (str_b, p_b, align_char, 'c', margin);

   if (vertical)
   {
      label_a.set_text_angle (M_PI_2 * 3);
      label_b.set_text_angle (M_PI_2 * 3);
   }

   label_a.cairo (cr);
   label_b.cairo (cr);

   cr->restore ();

}

void
Time_Chooser::render_month_lines (const RefPtr<Context>& cr,
                                  const Integer year,
                                  const Integer month,
                                  const Dstring& format) const
{

   if (format == "") { return; }

   for (Integer day = 1; day <= 21; day += 10)
   {

      const Dtime dtime (year, month, day);
      if (shape.out_of_bounds (dtime)) { return; }

      Point_2D point = get_point (dtime);
      if (vertical) { point.x = width - 2.0 * margin; }
      else          { point.y = 1.5 * margin; }

      const Dstring str = dtime.get_string (format);
      Label label (str, point, 'r', 'c');
      if (!vertical) { label.set_text_angle (M_PI_2 * 3); }
      label.cairo (cr);

   }

}

void
Time_Chooser::render_background_days (const RefPtr<Context>& cr)
{

   cr->save ();

   const Rect rect (Point_2D (0, 0), width, height);

   Color (1, 1, 1, 0.7).cairo (cr);
   rect.cairo (cr);
   cr->fill_preserve ();
   Color (0, 0, 0).cairo (cr);
   cr->stroke ();

   cr->set_font_size (font_size);

   cr->set_line_width (1);
   Color (0.5, 0.5, 0.5).cairo (cr);
   render_lines (cr, line_interval, "%HZ");

   cr->set_line_width (1);
   Color (0, 0, 0).cairo (cr);
   render_dates (cr);

   cr->restore ();

}

void
Time_Chooser::render_dates (const RefPtr<Context>& cr) const
{

   const Integer hours = 24;
   const Dtime st (shape.start_time.t - fmod (shape.start_time.t, hours));
   const Dstring time_format_str ("%b %d");

   typedef Point_2D P;
   const Real w = width;
   const Real h = height;
   const Real m = margin;
   const bool v = vertical;

   cr->save ();

   for (Dtime dtime = st; dtime <= shape.end_time; dtime.t += hours)
   {

      if (dtime < shape.start_time) { continue; }

      const Point_2D& point = get_point (dtime);
      const Real x = point.x;
      const Real y = point.y;

      const Point_2D p_a = (v ? P (w - 6.5 * margin, y) : P (x, h - 1.5 * m));
      const Point_2D p_b = (v ? P (w - 5.0 * margin, y) : P (x, h - 3.0 * m));
      const Point_2D p_c = (v ? P (w - 4.5 * margin, y) : P (x, h - 3.5 * m));
      const Point_2D p_d = (v ? P (w - 7.0 * margin, y) : P (x, h - 1.0 * m));

      const Dstring str_a = dtime.get_string ("%b %d");
      const Dstring str_b = dtime.get_string ("(%a)");

      cr->set_line_width (1);
      cr->move_to (p_c.x, p_c.y);
      cr->line_to (p_d.x, p_d.y);
      cr->stroke ();

      const bool align_left = (vertical == reverse);
      const char align_char = (align_left ? 'l' : 'r');

      Label label_a (str_a, p_a, align_char, 'c', margin);
      Label label_b (str_b, p_b, align_char, 'c', margin);

      if (vertical)
      {
         label_a.set_text_angle (M_PI_2 * 3);
         label_b.set_text_angle (M_PI_2 * 3);
      }

      label_a.cairo (cr);
      label_b.cairo (cr);

   }

   cr->restore ();

}

void
Time_Chooser::render_lines (const RefPtr<Context>& cr,
                            const Real hours,
                            const Dstring& format) const
{

   if (format == "") { return; }
   const Dtime& start_time = shape.start_time;
   const Dtime& end_time = shape.end_time;
   const Dtime st (start_time.t - fmod (start_time.t, double (hours)));

   for (Dtime dtime = st; dtime <= end_time; dtime.t += hours)
   {

      if (dtime < start_time) { continue; }

      Point_2D point = get_point (dtime);

      if (vertical) { point.x = width - 2.0 * margin; }
      else          { point.y = 1.5 * margin; }

      if (vertical)
      {
         Ring (1).cairo (cr, Point_2D (point.x + margin/3, point.y));
         cr->stroke ();
      }

      const Dstring& str = dtime.get_string (format);
      Label label (str, point, 'r', 'c');
      if (!vertical) { label.set_text_angle (M_PI_2 * 3); }
      label.cairo (cr);

   }

}

void
Time_Chooser::render_nodes (const RefPtr<Context>& cr) const
{

   typedef Point_2D P;
   const bool v = vertical;
   const Real split = 0.25;
   const Real node_size = margin / 3;

   cr->set_line_width (1);

   Dtime this_time, last_time;

   for (Shape::const_iterator iterator = shape.begin ();
        iterator != shape.end (); iterator++)
   {

      const Dtime& dtime = *(iterator);
      if (shape.out_of_bounds (dtime)) { return; }

      Rect rect = get_rect (dtime);
      rect.grow (vertical ? 0 : -split, vertical ? -split : 0);

      if (dtime.get_minute () % 1 == 0)
      {

         rect.cairo (cr);

         if (activated && data.is_candidate (dtime))
         {
            Color (0, 0.8, 0, 0.5).cairo (cr);
         }
         else
         if (data.matches (dtime))
         {
            Color (1, 0, 0, 0.5).cairo (cr);
         }
         else
         if (data.contains (dtime))
         {
            Color (0.7, 0.4, 0.4, 0.5).cairo (cr);
         }
         else
         {
            Color (0, 0, 0, 0.2).cairo (cr);
         }

         cr->fill ();

         if (data.matches (dtime))
         {
            rect.grow (1.5, 1.5);
            rect.cairo (cr);
            cr->stroke ();
         }

//         cr->fill_preserve ();
//         Color (0, 0, 0.4).cairo (cr);
//         cr->stroke ();

      }

/*
      if (data.matches (dtime))
      {

         const Real d = 2*node_size;
         const Triangle triangle_a (node_size, vertical ? 0 : M_PI_2);
         const Triangle triangle_b (node_size, vertical ? M_PI : M_PI_2 * 3);

         triangle_a.cairo (cr, point + P (v ? -d : 0, v ? 0 : -d));
         triangle_b.cairo (cr, point + P (v ? d : 0, v ? 0 : d));

         Color (1, 0, 0, 0.5).cairo (cr);
         cr->fill ();

      }
*/

      last_time = this_time;

   }

}

void
Time_Chooser::render_now (const RefPtr<Context>& cr) const
{
   const Dtime now;
   const Real s = get_s (now);

   cr->save ();
   cr->set_line_width (2.5);
   cr->set_font_size (12);
   Color (0.9, 0.2, 0.2, 0.5).cairo (cr);

   if (vertical)
   {
      const Point_2D point (width - 3.0 * margin, s);
      const Arrow arrow (0, 2*margin);
      arrow.cairo (cr, point);
      cr->stroke ();
      Label ("NOW", point + Point_2D (-1.5 * margin, 0), 'r', 'c').cairo (cr);
   }
   else
   {
      const Point_2D point (s, 4.0 * margin);
      const Arrow arrow (-M_PI/2, 2*margin);
      arrow.cairo (cr, point);
      cr->stroke ();
      Label ("NOW", point + Point_2D (0, 1.5 * margin), 'c', 't').cairo (cr);
   }

   cr->restore ();

}

void
Time_Chooser::init (const Real font_size)
{

   set_font_size (font_size);

   anim_button.get_full_signal ().connect (
      sigc::mem_fun (*this, &Time_Chooser::handle_animate));
   anim_button.get_click_full_signal ().connect (
      sigc::mem_fun (*this, &Time_Chooser::handle_dwell));
   prev_button.get_full_signal ().connect (
      sigc::mem_fun (*this, &Time_Chooser::step_backward));
   next_button.get_full_signal ().connect (
      sigc::mem_fun (*this, &Time_Chooser::step_forward));

   register_widget (anim_button);
   register_widget (prev_button);
   register_widget (next_button);

}

void
Time_Chooser::animate ()
{
   while (anim_button.is_switched_on ())
   {
      Glib::usleep (dwell);
      step_forward ();
   }
}

Time_Chooser::Time_Chooser (const Dcanvas& canvas,
                            const Real font_size,
                            const bool reverse,
                            const bool vertical,
                            const Real line_interval)
   : Dcontainer (canvas),
     anim_button (canvas, "ANIM", font_size, false),
     prev_button (canvas, "PREV", font_size),
     next_button (canvas, "NEXT", font_size),
     margin (font_size),
     thread_ptr (NULL),
     dwell (400000),
     reverse (reverse),
     vertical (vertical),
     line_interval (line_interval)
{
   data.clear ();
   init (font_size);
}

Time_Chooser::~Time_Chooser ()
{
   if (thread_ptr != NULL)
   {
      anim_button.set (false);
      thread_ptr->join ();
      thread_ptr = NULL;
   }
}

void
Time_Chooser::set_vertical (const bool vertical)
{
   this->vertical = vertical;
}

const Time_Chooser::Shape&
Time_Chooser::get_shape () const
{
   return shape;
}

Time_Chooser::Shape&
Time_Chooser::get_shape ()
{
   return shape;
}

const Time_Chooser::Data&
Time_Chooser::get_data () const
{
   return data;
}

Time_Chooser::Data&
Time_Chooser::get_data ()
{
   return data;
}

Time_Chooser::Signal&
Time_Chooser::get_signal ()
{
   return signal;
}

Time_Chooser::Data_Signal&
Time_Chooser::get_data_signal ()
{
   return data_signal;
}

Time_Chooser::Shape
Time_Chooser::get_highlighted_shape () const
{
   return data.get_highlighted_shape (shape);
}

void
Time_Chooser::set_reverse (const bool reverse)
{
   this->reverse = reverse;
}

void
Time_Chooser::set_line_interval (const Real line_interval)
{
   this->line_interval = line_interval;
}

void
Time_Chooser::set_time_chooser (const Time_Chooser& time_chooser)
{
   shape = time_chooser.shape;
   data = time_chooser.data;
}

bool
Time_Chooser::set_data (const Time_Chooser::Data& data)
{
   this->data = data;
   return true;
}

void
Time_Chooser::set_shape (const Shape& shape)
{
   this->shape = shape;
   const bool short_span = (shape.get_span_t () < 36);
   set_line_interval (short_span ? 1 : 6);

   data.conform_to (shape);
}

const Dtime&
Time_Chooser::get_time () const
{
   return data.dtime;
}

void
Time_Chooser::handle_animate (const Devent& event)
{
   if (anim_button.is_switched_on ())
   {
      thread_ptr = Glib::Thread::create (
         sigc::mem_fun (*this, &Time_Chooser::animate), true);
   }
   else
   {
      thread_ptr->join ();
      thread_ptr = NULL;
   }
}

void
Time_Chooser::handle_dwell (const Devent& event)
{

   Gtk::Dialog dialog ("Change Dwell...", true);

   Glib::RefPtr<Gtk::Adjustment> a (Gtk::Adjustment::create (
      dwell / 1000.0, 200, 4000, 200, 500));
   Gtk::SpinButton dwell_spin_button (a);

   dialog.get_vbox ()->add (dwell_spin_button);
   dialog.add_button ("Ok", Gtk::RESPONSE_OK);
   dialog.add_button ("Cancel", Gtk::RESPONSE_CANCEL);

   dialog.show_all_children ();
   guint result = dialog.run ();

   switch (result)
   {
      case Gtk::RESPONSE_OK:
      {
         //guint yyyy_f, mm_f, dd_f;
         //guint yyyy_t, mm_t, dd_t;
         //calendar_from.get_date (yyyy_f, mm_f, dd_f);
         //calendar_to.get_date (yyyy_t, mm_t, dd_t);

         //shape.start_time = Dtime (yyyy_f, mm_f + 1, dd_f);
         //shape.end_time = Dtime (yyyy_t, mm_t + 1, dd_t);
         dwell = Integer (round (dwell_spin_button.get_value ())) * 1000;
         cout << "dwell is now " << dwell << endl;
         canvas.queue_draw ();
         break;
      }
   }

}

void
Time_Chooser::step_forward (const Devent& event)
{

   const bool is_range = event.shift ();
   const bool is_leap = event.control ();

   try
   {

      if (is_range)
      {
         const Dtime& new_time = shape.get_next_time (data.dtime, is_leap);
         data.set_time (new_time, false, is_range, is_leap);
      }
      else
      {
         const Shape& highlighted_shape = get_highlighted_shape ();
         const bool singleton = (highlighted_shape.size () == 1);
         const Shape& s = (singleton ? shape : highlighted_shape);
         const Dtime& new_time = s.get_next_time (data.dtime, is_leap);
         data.set_time (new_time, singleton, is_range, is_leap);
      }

      signal.emit ();
      data_signal.emit (data);

   }
   catch (const Exception& e)
   {
   }

}

void
Time_Chooser::step_backward (const Devent& event)
{

   const bool is_range = event.shift ();
   const bool is_leap = event.control ();

   try
   {

      if (is_range)
      {
         const Dtime& new_time = shape.get_prev_time (data.dtime, is_leap);
         data.set_time (new_time, false, is_range, is_leap);
      }
      else
      {
         const Shape& highlighted_shape = get_highlighted_shape ();
         const bool singleton = (highlighted_shape.size () == 1);
         const Shape& s = (singleton ? shape : highlighted_shape);
         const Dtime& new_time = s.get_prev_time (data.dtime, is_leap);
         data.set_time (new_time, singleton, is_range, is_leap);
      }

      signal.emit ();
      data_signal.emit (data);

   }
   catch (const Exception& e)
   {
   }

}

void
Time_Chooser::go_to_first ()
{
   data.set_time (*(shape.begin ()));
   signal.emit ();
   data_signal.emit (data);
}

void
Time_Chooser::go_to_last ()
{
   data.set_time (*(shape.rbegin ()));
   signal.emit ();
   data_signal.emit (data);
}

void
Time_Chooser::toggle_anim_button ()
{
   anim_button.toggle ();
}

void
Time_Chooser::set_dwell (const Integer dwell)
{
   this->dwell = dwell;
}

void
Time_Chooser::set_leap (const Real leap)
{
   this->shape.leap = leap;
}

Point_2D
Time_Chooser::get_point (const Dtime& dtime) const
{
   const Real s = get_s (dtime);
   return (vertical ? Point_2D (width - margin, s) : Point_2D (s, margin));
}

Rect
Time_Chooser::get_rect (const Dtime& dtime) const
{

   const Real w = 4;
   const Real centre = (vertical ? width - margin : margin);

   Shape::const_iterator iterator = shape.get_iterator (dtime);
   Shape::const_iterator next_iterator = iterator; next_iterator++;

   const bool first = (iterator == shape.begin ());
   const bool last = (next_iterator == shape.end ());

   Dtime prev_time, next_time;

   if (first)
   {
      prev_time = *(iterator);
      next_time.t = (next_iterator->t + iterator->t) * 0.5;
   }
   else
   if (last)
   {
      Shape::const_iterator prev_iterator = iterator; prev_iterator--;
      prev_time.t = (prev_iterator->t + iterator->t) * 0.5;
      next_time = *(iterator);
   }
   else
   {
      Shape::const_iterator prev_iterator = iterator; prev_iterator--;
      prev_time.t = (prev_iterator->t + iterator->t) * 0.5;
      next_time.t = (next_iterator->t + iterator->t) * 0.5;
   }

   const Real prev_s = get_s (prev_time);
   const Real next_s = get_s (next_time);

   return vertical ?
      Rect (Point_2D (centre - w, prev_s), Point_2D (centre + w, next_s)) :
      Rect (Point_2D (prev_s, centre - w), Point_2D (next_s, centre + w));

}

void
Time_Chooser::cairo (const RefPtr<Context>& cr)
{

   cr->save ();
   cr->translate (anchor.x, anchor.y);

   render_background (cr);
   render_nodes (cr);
   render_now (cr);

   cr->translate (-anchor.x, -anchor.y);
   cr->restore ();

   Dcontainer::cairo (cr);

}

bool
Time_Chooser::on_mouse_button_pressed (const Dmouse_Button_Event& event)
{

   if (Dcontainer::on_mouse_button_pressed (event)) { return true; }

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   if (out_of_bounds (point)) { return false; }

   if (event.button == 1)
   {

      if (event.type == GDK_3BUTTON_PRESS)
      {

         Gtk::Dialog dialog ("Change Time...", true);
         Gtk::Calendar calendar_from;
         Gtk::Calendar calendar_to;

         const Integer start_year = shape.start_time.get_year ();
         const Integer end_year = shape.end_time.get_year ();
         const Integer start_month = shape.start_time.get_month ();
         const Integer end_month = shape.end_time.get_month ();
         const Integer start_day = shape.start_time.get_day ();
         const Integer end_day = shape.end_time.get_day ();

         calendar_from.select_month (start_month - 1, start_year);
         calendar_from.select_day (start_day);
         calendar_to.select_month (end_month - 1, end_year);
         calendar_to.select_day (end_day);

         dialog.get_vbox ()->add (calendar_from);
         dialog.get_vbox ()->add (calendar_to);
         dialog.add_button ("Ok", Gtk::RESPONSE_OK);
         dialog.add_button ("Cancel", Gtk::RESPONSE_CANCEL);

         dialog.show_all_children ();
         guint result = dialog.run ();

         switch (result)
         {
            case Gtk::RESPONSE_OK:
            {
               guint yyyy_f, mm_f, dd_f;
               guint yyyy_t, mm_t, dd_t;
               calendar_from.get_date (yyyy_f, mm_f, dd_f);
               calendar_to.get_date (yyyy_t, mm_t, dd_t);

               shape.start_time = Dtime (yyyy_f, mm_f + 1, dd_f);
               shape.end_time = Dtime (yyyy_t, mm_t + 1, dd_t);
               canvas.queue_draw ();
               break;
            }
         }

         return true;
      }

      const Dtime& dtime = get_time (point, true);
      if (!dtime.is_nat ())
      {
         activate ();
         data.candidate_time = dtime;
         const Integer ii = Integer (round (anchor.x));
         const Integer jj = Integer (round (anchor.y));
         const Integer w = Integer (round (width));
         const Integer h = Integer (round (height));
         canvas.queue_draw_area (ii, jj, w, h);
      }

      return true;

   }

   return false;

}

bool
Time_Chooser::on_mouse_motion (const Dmouse_Motion_Event& event)
{

   if (Dcontainer::on_mouse_motion (event)) { return true; }
   if (!activated) { return false; }

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   //if (out_of_bounds (point)) { return false; }

   if (get_time (point, true) != data.candidate_time)
   {
      deactivate ();
   }
   else
   {
      activate ();
   }

   const Integer i = Integer (round (anchor.x));
   const Integer j = Integer (round (anchor.y));
   const Integer w = Integer (round (width));
   const Integer h = Integer (round (height));
   canvas.queue_draw_area (i, j, w, h);
   return true;

}

Real
Time_Chooser::get_preferred_width () const
{
   if (vertical) { return 8 * font_size; }
   else
   {
      const Real button_width =  3 * font_size + margin;
      return button_width + 2 * margin;
   }
}

Real
Time_Chooser::get_preferred_height () const
{
   if (!vertical) { return 2 * font_size + 5 * margin; }
   else
   {
      const Real button_height = font_size + margin;
      return button_height + 2 * margin;
   }
}

void
Time_Chooser::pack ()
{

   if (vertical)
   {
      // vertical

      const Real button_w = (width - 1.5 * margin) / 2;
      const Real button_h = font_size * 1.5;
      const Real anchor_x = 0.5 * margin;
      const Real anchor_y = 0.5 * margin;

      button_box.index_2d.i = Integer (round (anchor_x));
      button_box.index_2d.j = Integer (round (anchor_y));
      button_box.size_2d.i = Integer (round (width - margin));


      const Real anim_button_w = Real (button_box.size_2d.i);

      const Point_2D anim_anchor (anchor_x, anchor_y);
      const Point_2D prev_anchor (anchor_x, margin + button_h);
      const Point_2D next_anchor (margin + button_w, margin + button_h);

      anim_button.being_packed (anim_anchor, anim_button_w, button_h);
      prev_button.being_packed (prev_anchor, button_w, button_h);
      next_button.being_packed (next_anchor, button_w, button_h);

      button_box.size_2d.j = Integer (round (2 * button_h + 0.5 * margin));

   }
   else
   {
      // horizontal

      const Real button_w = 3 * font_size + margin;
      const Real anchor_x = 0.75 * margin;
      const Real anchor_y = 0.5 * margin;

      button_box.index_2d.i = Integer (round (anchor_x));
      button_box.index_2d.j = Integer (round (anchor_y));
      button_box.size_2d.i = Integer (round (button_w));
      button_box.size_2d.j = Integer (round (height - 2 * margin));

      const Real button_h = Real (button_box.size_2d.j) / 3;

      const Point_2D anim_anchor (anchor_x, 0.5 * margin);
      const Point_2D prev_anchor (anchor_x, margin + button_h);
      const Point_2D next_anchor (anchor_x, 1.5*margin + 2*button_h);

      anim_button.being_packed (anim_anchor, button_w, button_h);
      prev_button.being_packed (prev_anchor, button_w, button_h);
      next_button.being_packed (next_anchor, button_w, button_h);

   }

}

bool
Time_Chooser::on_mouse_button_released (const Dmouse_Button_Event& event)
{

   if (Dcontainer::on_mouse_button_released (event)) { return true; }
   if (!activated) { return false; }

   const Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   if (out_of_bounds (point)) { return false; }

   switch (event.button)
   {

      case 1:
      {

         const Dtime& dtime = get_time (point, true);

         if (dtime.is_nat ())
         {
            data.candidate_time = GSL_NAN;
            deactivate ();
            canvas.queue_draw ();
            return true;
         }

         if (dtime == data.candidate_time)
         {
            selected ();
            data.set_time (dtime, event);
            deactivate ();
            canvas.queue_draw ();
            signal.emit ();
            data_signal.emit (data);
            return true;
         }

      }

   }

   return false;

}

void
Time_Chooser::selected ()
{
}

Time_Canvas::Time_Canvas (const Dcanvas& canvas,
                          const Real font_size,
                          const bool reverse,
                          const bool vertical,
                          const Real line_interval)
   : canvas (canvas),
     time_chooser (canvas, font_size, reverse, vertical, line_interval)
{
}

const Time_Chooser&
Time_Canvas::get_time_chooser () const
{
   return time_chooser;
}

Time_Chooser&
Time_Canvas::get_time_chooser ()
{
   return time_chooser;
}

bool
Time_Canvas::on_key_pressed (const Dkey_Event& event)
{

   switch (event.value)
   {

      case GDK_KEY_Home:
      {
         time_chooser.go_to_first ();
         return true;
         break;
      };

      case GDK_KEY_End:
      {
         time_chooser.go_to_last ();
         return true;
         break;
      };

      case GDK_KEY_Left:
      {
         time_chooser.step_backward (event);
         return true;
         break;
      }

      case GDK_KEY_Right:
      {
         time_chooser.step_forward (event);
         return true;
         break;
      }

      case GDK_KEY_a:
      case GDK_KEY_A:
      {
         time_chooser.toggle_anim_button ();
         time_chooser.handle_animate (event);
         return true;
         break;
      }

   }

   return false;

}

const Level_Tuple&
Level_Panel::get_level_tuple () const
{

   switch (level.type)
   {
      default:
      case Level::HEIGHT:   return level_tuple_z;
      case Level::PRESSURE: return level_tuple_p;
      case Level::SIGMA:    return level_tuple_sigma;
      case Level::THETA:    return level_tuple_theta;
   }

   throw Nwp_Exception ("Level_Tuple::get_level_tuple confused");

}

Integer
Level_Panel::get_extra_level_index () const
{

   for (Integer i = 0; i < extra_level_vector.size (); i++)
   {
      const Level& l = extra_level_vector[i];
      if (level.get_string () == l.get_string ()) { return i; }
   }

   return -1;

}

Level
Level_Panel::get_level (const Real y) const
{

   const Real fs2 = font_size / 2;

   if (y < start_margin - fs2 || y > height - end_margin + fs2)
   {
      for (Integer i = 0; i < extra_level_vector.size (); i++)
      {
         const Real y_i = height - (2 * i + 1.5) * font_size;
         if (fabs (y_i - y) < font_size / 2)
         {
            const Level& level = extra_level_vector[i];
            return level;
         }
      }
      Level level;
      level.type = Level::NAL;
      return level;
   }
   else
   {
      switch (level.type)
      {
         default:
            return Level (Level::NAL, GSL_NAN);
         case Level::SCREEN:
         case Level::TEN_METRE:
         case Level::FIFTY_METRE:
         case Level::HEIGHT:
            return get_level_z (y);
         case Level::PRESSURE:
            return get_level_p (y);
         case Level::THETA:
            return get_level_theta (y);
         case Level::SIGMA:
            return get_level_sigma (y);
      }
   }

}

Level
Level_Panel::get_level_z (const Real y) const
{

   const Real z0 = level_tuple_z.front ();
   const Real z1 = level_tuple_z.back ();

   const Real c = 400;

   const Real span_y = height - end_margin - start_margin;
   const Real x0 = log (fabs (c - z0));
   const Real x1 = log (fabs (c - z1));
   const Real m = span_y / (x1 - x0);

   const Real z = (exp ((y - start_margin) / m + x0) - c);

   Real nearest_y, nearest_z;
   Real min_delta_y = GSL_POSINF;
   for (Level_Tuple::const_iterator iterator = level_tuple_z.begin ();
        iterator != level_tuple_z.end (); iterator++)
   {
      const Real z_i = *(iterator);
      const Real y_i = get_y_z (z_i);
      const Real delta_y = fabs (y - y_i);
      if (min_delta_y > delta_y)
      {
         nearest_y = y_i;
         nearest_z = z_i;
         min_delta_y = delta_y;
      }
   }

   const bool close = (fabs (nearest_y - y) < font_size / 2);
   return Level::z_level (close ? nearest_z : z);

}

Level
Level_Panel::get_level_p (const Real y) const
{

   const Real p0 = level_tuple_p.front ();
   const Real p1 = level_tuple_p.back ();

   const Real c = 1400e2;

   const Real span_y = height - end_margin - start_margin;
   const Real x0 = log (c - p0);
   const Real x1 = log (c - p1);
   const Real m = span_y / (x1 - x0);

   const Real p = -(exp ((y - start_margin) / m + x0) - c);

   Real nearest_y, nearest_p;
   Real min_delta_y = GSL_POSINF;
   for (Level_Tuple::const_iterator iterator = level_tuple_p.begin ();
        iterator != level_tuple_p.end (); iterator++)
   {
      const Real p_i = *(iterator);
      const Real y_i = get_y_p (p_i);
      const Real delta_y = fabs (y - y_i);
      if (min_delta_y > delta_y)
      {
         nearest_y = y_i;
         nearest_p = p_i;
         min_delta_y = delta_y;
      }
   }

   const bool close = (fabs (nearest_y - y) < font_size / 2);
   return Level::pressure_level (close ? nearest_p : p);

}

Level
Level_Panel::get_level_theta (const Real y) const
{

   const Real theta0 = level_tuple_theta.front ();
   const Real theta1 = level_tuple_theta.back ();

   const Real span_y = height - end_margin - start_margin;
   const Real x0 = theta0;
   const Real x1 = theta1;
   const Real m = span_y / (x1 - x0);

   const Real theta = (y - start_margin) / m + x0;

   Real nearest_y, nearest_theta;
   Real min_delta_y = GSL_POSINF;
   for (Level_Tuple::const_iterator iterator = level_tuple_theta.begin ();
        iterator != level_tuple_theta.end (); iterator++)
   {
      const Real theta_i = *(iterator);
      const Real y_i = get_y_theta (theta_i);
      const Real delta_y = fabs (y - y_i);
      if (min_delta_y > delta_y)
      {
         nearest_y = y_i;
         nearest_theta = theta_i;
         min_delta_y = delta_y;
      }
   }

   const Real font_size = 12;
   const bool close = (fabs (nearest_y - y) < font_size / 2);
   return Level::theta_level (close ? nearest_theta : theta);

}

Level
Level_Panel::get_level_sigma (const Real y) const
{

   const Real sigma0 = level_tuple_sigma.front ();
   const Real sigma1 = level_tuple_sigma.back ();

   const Real c = 1.02;

   const Real span_y = height - end_margin - start_margin;
   const Real x0 = log (c - sigma0);
   const Real x1 = log (c - sigma1);
   const Real m = span_y / (x1 - x0);

   const Real sigma = -(exp ((y - start_margin) / m + x0) - c);

   Real nearest_y, nearest_sigma;
   Real min_delta_y = GSL_POSINF;
   for (Level_Tuple::const_iterator iterator = level_tuple_sigma.begin ();
        iterator != level_tuple_sigma.end (); iterator++)
   {
      const Real sigma_i = *(iterator);
      const Real y_i = get_y_sigma (sigma_i);
      const Real delta_y = fabs (y - y_i);
      if (min_delta_y > delta_y)
      {
         nearest_y = y_i;
         nearest_sigma = sigma_i;
         min_delta_y = delta_y;
      }
   }

   const bool close = (fabs (nearest_y - y) < font_size / 2);
   return Level::sigma_level (close ? nearest_sigma : sigma);

}

Real
Level_Panel::get_y (const Level& level) const
{

   for (Integer i = 0; i < extra_level_vector.size (); i++)
   {
      const Level& l = extra_level_vector[i];
      const Dstring& str = l.get_string ();

      if (str == level.get_string ())
      {
         const Real y = height - (2 * i + 1.5) * font_size;
         return y;
      }
   }

   switch (level.type)
   {
      default:
      case Level::HEIGHT:
         return get_y_z (level.value);
      case Level::PRESSURE:
         return get_y_p (level.value);
      case Level::THETA:
         return get_y_theta (level.value);
      case Level::SIGMA:
         return get_y_sigma (level.value);
   }
}

Real
Level_Panel::get_y_z (const Real z) const
{

   const Real z0 = level_tuple_z.front ();
   const Real z1 = level_tuple_z.back ();

   const Real c = 400;

   const Real span_y = height - end_margin - start_margin;
   const Real x0 = log (fabs (c - z0));
   const Real x1 = log (fabs (c - z1));
   const Real m = span_y / (x1 - x0);

   return (log (c + z) - x0) * m + start_margin; 

}

Real
Level_Panel::get_y_p (const Real p) const
{

   const Real p0 = level_tuple_p.front ();
   const Real p1 = level_tuple_p.back ();

   const Real c = 1400e2;

   const Real span_y = height - end_margin - start_margin;
   const Real x0 = log (c - p0);
   const Real x1 = log (c - p1);
   const Real m = span_y / (x1 - x0);

   return (log (c - p) - x0) * m + start_margin; 

}

Real
Level_Panel::get_y_theta (const Real theta) const
{

   const Real theta0 = level_tuple_theta.front ();
   const Real theta1 = level_tuple_theta.back ();

   const Real span_y = height - end_margin - start_margin;
   const Real x0 = theta0;
   const Real x1 = theta1;
   const Real m = span_y / (x1 - x0);

   return (theta - x0) * m + start_margin; 

}

Real
Level_Panel::get_y_sigma (const Real sigma) const
{

   const Real sigma0 = level_tuple_sigma.front ();
   const Real sigma1 = level_tuple_sigma.back ();

   const Real c = 1.02;

   const Real span_y = height - end_margin - start_margin;
   const Real x0 = log (c - sigma0);
   const Real x1 = log (c - sigma1);
   const Real m = span_y / (x1 - x0);

   return (log (c - sigma) - x0) * m + start_margin; 

}

bool
Level_Panel::is_close_to_level_0 (const Real y) const
{
   const Level& level_0 = level.get_level_0 ();
   const Real y_0 = get_y (level_0);
   return (fabs (y - y_0) < font_size / 2);
}

bool
Level_Panel::is_close_to_level_1 (const Real y) const
{
   const Level& level_1 = level.get_level_1 ();
   const Real y_1 = get_y (level_1);
   return (fabs (y - y_1) < font_size / 2);
}

bool
Level_Panel::is_within_layer (const Real y) const
{
   const Level& level_0 = level.get_level_0 ();
   const Level& level_1 = level.get_level_1 ();
   const Real y_0 = get_y (level_0);
   const Real y_1 = get_y (level_1);
   return ((y - y_0) * (y - y_1) <= 0);
}

void
Level_Panel::render_background (const RefPtr<Context> cr) const
{

   cr->set_font_size (font_size);
   Color (0, 0, 0, 0.5).cairo (cr);

//   if (level.type == Level::PRESSURE ||
//       level.type == Level::THETA ||
//       level.type == Level::SIGMA)
   {


      Level::Type level_type = level.type;

      switch (level.type)
      {
         case Level::SCREEN:   
         case Level::FIFTY_METRE:
         case Level::TEN_METRE:
            level_type = Level::PRESSURE;
            break;
      }

      const Level_Tuple& level_tuple = get_level_tuple ();

      for (Level_Tuple::const_iterator iterator = level_tuple.begin ();
           iterator != level_tuple.end (); iterator++)
      {
         const Level l (level_type, *(iterator), GSL_NAN);
         const Real y = get_y (l);
         const Dstring& str = l.get_string ();
         Label label (str, Point_2D (width / 2, y), 'c', 'c');
         label.cairo (cr);
      }

   }

   for (Integer i = 0; i < extra_level_vector.size (); i++)
   {
      const Level& l = extra_level_vector[i];
      const Dstring& str = l.get_string ();
      const Real y = height - (2 * i + 1.5) * font_size;
      Label label (str, Point_2D (width / 2, y), 'c', 'c');
      label.cairo (cr);
   }

}

void
Level_Panel::render_level (const RefPtr<Context> cr,
                           const Level& level,
                           const bool with_triangles) const
{

   if (level.is_nal ()) { return; }

   const Real ns = 7;
   const Real y = get_y (level);
   const Dstring& str = level.get_string ();

   cr->set_font_size (font_size);
   Color (0, 0, 0, 0.8).cairo (cr);
   Label label (str, Point_2D (width / 2, y), 'c', 'c');
   label.cairo (cr);

   if (with_triangles)
   {
      Triangle (ns).cairo (cr, Point_2D (ns, y));
      Triangle (ns, M_PI).cairo (cr, Point_2D (width-ns, y));
      Color (1, 0, 0, 0.5).cairo (cr);
      cr->fill ();
   }

}

void
Level_Panel::render_layer (const RefPtr<Context> cr,
                           const Level& level) const
{

   if (!level.is_layer ()) { return; }

   const Real ns = font_size / 2;

   const Level& level_0 = level.get_level_0 ();
   const Level& level_1 = level.get_level_1 ();
   const Real y_0 = get_y (level_0);
   const Real y_1 = get_y (level_1);
   const Dstring& str_0 = level_0.get_string ();
   const Dstring& str_1 = level_1.get_string ();

   const Color color_0 (0.6, 0.6, 0, 0.2);
   const Color color_1 (0.0, 0.0, 0, 0.2);
   Stripped (color_0, color_1, 20).cairo (cr);

   Rect (Point_2D (5, y_0), width - 10, y_1 - y_0).cairo (cr);
   cr->fill ();

   cr->set_font_size (font_size);
   Color (0, 0, 0, 0.8).cairo (cr);
   Label label_0 (str_0, Point_2D (width / 2, y_0), 'c', 'c');
   Label label_1 (str_1, Point_2D (width / 2, y_1), 'c', 'c');
   label_0.cairo (cr);
   label_1.cairo (cr);

   if (true)
   {
      Triangle (ns).cairo (cr, Point_2D (ns, y_0));
      Triangle (ns).cairo (cr, Point_2D (ns, y_1));
      Triangle (ns, M_PI).cairo (cr, Point_2D (width-ns, y_0));
      Triangle (ns, M_PI).cairo (cr, Point_2D (width-ns, y_1));
      Color (1, 0, 0, 0.5).cairo (cr);
      cr->fill ();
   }

}

bool
Level_Panel::on_mouse_button_pressed (const Dmouse_Button_Event& event)
{

   if (level.type == Level::NAL && !level.is_layer ())
   {
      return false;
   }

   Point_2D point (event.point.x - anchor.x, event.point.y - anchor.y);
   if (out_of_bounds (point)) { return false; }

   if (level.is_layer ())
   {

      // It is a layer
      setting_level = true;

      if (is_within_layer (point.y))
      {
         reference_y = point.y;
         return true;
      }

      if (is_close_to_level_0 (point.y))
      {
         setting_level_0 = true;
      }
      else
      if (is_close_to_level_1 (point.y))
      {
         setting_level_1 = true;
      }
      else
      {
         const Level& l = get_level (point.y);
         level.value = l.get_value ();
         level.value_ = l.get_value ();
      }
      canvas.queue_draw ();
      return true;
   }
   else
   {
      // Not a layer, just a simple level
      setting_level = true;
      candidate_level = get_level (point.y);
      canvas.queue_draw ();
      return true;
   }

}

bool
Level_Panel::on_mouse_motion (const Dmouse_Motion_Event& event)
{

   const Real fs = font_size;
   Real y = event.point.y - anchor.y;

   if (level.is_layer ())
   {

      if (!setting_level) { return false; }

      const Level& l = get_level (y);

      if (!gsl_isnan (reference_y))
      {
         const Level& reference_level = get_level (reference_y);
         const Level& level = get_level (y);
         const Real d_value = level.value - reference_level.value;
         this->level.value += d_value;
         this->level.value_ += d_value;
         reference_y = y;
      }
      else
      if (setting_level_0)
      {
         level.value = l.get_value ();
      }
      else
      {
         level.value_ = l.get_value ();
      }

      level_signal.emit (get_level ());
      full_level_signal.emit (get_level (), event);

      canvas.queue_draw ();
      return true;

   }
   else
   {
      if (setting_level)
      {
         candidate_level = get_level (y);
         if (candidate_level.type != Level::NAL)
         {
            level = candidate_level;
            level_signal.emit (get_level ());
            full_level_signal.emit (get_level (), event);
            canvas.queue_draw ();
         }
         return true;
      }
   }

   return false;

}

bool
Level_Panel::on_mouse_button_released (const Dmouse_Button_Event& event)
{

   const Real fs = font_size;
   Real y = event.point.y - anchor.y;

   if (level.is_layer ())
   {

      if (!setting_level) { return false; }

      if (!gsl_isnan (reference_y))
      {
         reference_y = GSL_NAN;
      }
      else
      {
         const Level& l = get_level (y);
         if (setting_level_0) { this->level.value = l.get_value (); }
         else                 { this->level.value_ = l.get_value (); }

         setting_level_0 = false;
         setting_level_1 = false;
         level.order ();
      }

      setting_level = false;
      level_signal.emit (get_level ());
      full_level_signal.emit (get_level (), event);
      canvas.queue_draw ();
      return true;

   }
   else
   {
      if (setting_level)
      {

         setting_level = false;
         candidate_level.type = Level::NAL;
         const Real y = event.point.y - anchor.y;
         const Level& l = get_level (y);

         if (l.type != Level::NAL)
         {
            level = l;
            level_signal.emit (get_level ());
            full_level_signal.emit (get_level (), event);
            canvas.queue_draw ();
            candidate_level.type = Level::NAL;
            return true;
         }

      }
   }

   return false;

}

Level_Panel::Level_Panel (Dcanvas& dcanvas,
                          const Real font_size)
   : Dv_Pack_Box (dcanvas, 0, 6/*margin*/),
     dcanvas (dcanvas),
     level_tuple_z (Level::HEIGHT, "40000:36000:32000:28000:24000:20000:18000:16000:14000:12000:10000:9000:8000:7000:6000:5000:4500:4000:3500:3000:2500:2000:1800:1600:1400:1200:1000:900:800:700:600:500:400:350:300:250:200:150:100:75:45:20:5:0"),
     level_tuple_sigma (Level::SIGMA, "0.2:0.25:0.3:0.35:0.4:0.45:0.5:0.55:0.6:0.65:0.7:0.75:0.8:0.85:0.875:0.9:0.925:0.95:0.975:0.9943:0.9975:0.9988"),
     level_tuple_theta (Level::THETA, "355:350:345:340:335:330:325:320:315:310:305:300:295:290:285:280:275:270"),
     level_tuple_p (Level::PRESSURE, "200e2:250e2:300e2:350e2:400e2:450e2:500e2:550e2:600e2:650e2:700e2:750e2:800e2:850e2:900e2:925e2:950e2:975e2:995e2:1000e2"),
     start_margin (50),
     end_margin (90),

     font_size (font_size),

     setting_level_0 (false),
     setting_level_1 (false),

     z (500),
     pressure (850e2),
     theta (315),
     sigma (0.5),
     setting_level (false),
     level (Level::screen_level ()),

     reference_y (GSL_NAN)

     //switcher_panel (dcanvas, 0, 6),
     //pressure_button (dcanvas, "P", font_size),
     //theta_button (dcanvas, "TH", font_size),
     //sigma_button (dcanvas, "S", font_size)
{

   level.type = Level::NAL;
   candidate_level.type = Level::NAL;

   //pressure_button.get_signal ().connect (
   //   sigc::mem_fun (*this, &Level_Panel::set_p_buttons));
   //theta_button.get_signal ().connect (
   //   sigc::mem_fun (*this, &Level_Panel::set_theta_buttons));
   //sigma_button.get_signal ().connect (
   //   sigc::mem_fun (*this, &Level_Panel::set_sigma_buttons));

   //switcher_panel.pack_back (pressure_button);
   //switcher_panel.pack_back (theta_button);
   //switcher_panel.pack_back (sigma_button);

//   level_set_signal.connect (sigc::mem_fun (
//      *this, &Level_Panel::set_level));
//   layer_set_signal.connect (sigc::mem_fun (
//      *this, &Level_Panel::set_layer));

}

Level_Panel::Level_Signal&
Level_Panel::get_level_signal ()
{
   return level_signal;
}

Level_Panel::Full_Level_Signal&
Level_Panel::get_full_level_signal ()
{
   return full_level_signal;
}

void
Level_Panel::add_extra_level (const Level& level)
{
   extra_level_vector.push_back (level);
}

void
Level_Panel::cairo (const RefPtr<Context>& cr)
{
//   if (level.type == Level::NAL && !level.is_layer ())
   if (level.type == Level::NAL)
   {
      return;
   }

   Dv_Pack_Box::cairo (cr);

   cr->save ();
   cr->translate (anchor.x, anchor.y);

   const Rect rect (Point_2D (0, 0), width, height);
   cr->set_line_width (2);
   Color (1, 1, 1, 0.7).cairo (cr);
   rect.cairo (cr);
   cr->fill_preserve ();
   Color (0, 0, 0, 0.8).cairo (cr);
   cr->stroke ();

   render_background (cr);

   if (level.is_layer ())
   {
      render_layer (cr, level);
   }
   else
   {
      render_level (cr, candidate_level, false);
      render_level (cr, level, true);
   }

   cr->translate (-anchor.x, -anchor.y);
   cr->restore ();

}

void
Level_Panel::set_level (const Level& level)
{
   this->level = level;
   switch (level.type)
   {
      case Level::PRESSURE: pressure = level.value; break;
      case Level::THETA: theta = level.value; break;
      case Level::SIGMA: sigma = level.value; break;
   }
}

const Level&
Level_Panel::get_level () const
{
   return this->level;
}

Level&
Level_Panel::get_level ()
{
   return this->level;
}

void
Level_Panel::clear ()
{
   Dv_Pack_Box::clear ();
   pack ();
}

void
Level_Panel::set_no_buttons ()
{

   level.type = Level::NAL;
   level.value_ = GSL_NAN;
   extra_level_vector.clear ();

//   clear ();
//   pack_back (switcher_panel);

//   pack ();

}

void
Level_Panel::set_theta_buttons ()
{

   if (level.type != Level::THETA)
   {
      level.value = theta;
   }

   level.type = Level::THETA;
   level.value_ = GSL_NAN;
   extra_level_vector.clear ();

//   clear ();
//   pack_back (switcher_panel);

//   pack ();

}

void
Level_Panel::set_sigma_buttons ()
{

   if (level.type != Level::SIGMA)
   {
      level.value = sigma;
   }

   level.type = Level::SIGMA;
   level.value_ = GSL_NAN;
   extra_level_vector.clear ();

//   clear ();
//   pack_back (switcher_panel);

//   pack ();

}

void
Level_Panel::set_p_buttons ()
{

   if (level.type == Level::HEIGHT ||
       level.type == Level::THETA ||
       level.type == Level::SIGMA ||
       level.type == Level::NAL)
   {
      level.value = pressure;
      level.type = Level::PRESSURE;
      level.value_ = GSL_NAN;
   }

   extra_level_vector.clear ();
//   clear ();
//   pack_back (switcher_panel);

//   pack ();

}

void
Level_Panel::set_wind_level_buttons ()
{

   set_p_buttons ();
   extra_level_vector.clear ();
   extra_level_vector.push_back (Level ("10m"));
   extra_level_vector.push_back (Level ("50m"));

   if (level.type == Level::SCREEN)
   {
      level.type = Level::TEN_METRE;
   }


   //extra_level_vector.push_back (Level ("0.995"));

//   clear ();
//   pack_back (switcher_panel);
//   typedef vector<Level_Button*>::const_iterator Iterator;
//
//   for (Iterator iterator = pressure_level_button_ptr_vector.begin ();
//        iterator != pressure_level_button_ptr_vector.end (); iterator++)
//   {
//      Level_Button* button_ptr = *(iterator);
//      pack_back (*button_ptr);
//   }
//
//   pack_back (sigma_995_button);
//   pack_back (fifty_metre_button);
//   pack_back (ten_metre_button);
//   pack ();

}

void
Level_Panel::set_temperature_level_buttons ()
{

   set_p_buttons ();
   extra_level_vector.clear ();
   extra_level_vector.push_back (Level ("Screen"));

   if (level.type == Level::TEN_METRE ||
       level.type == Level::FIFTY_METRE)
   {
      level.type = Level::SCREEN;
   }

//   clear ();
//   pack_back (switcher_panel);
//   typedef vector<Level_Button*>::const_iterator Iterator;
//
//   for (Iterator iterator = pressure_level_button_ptr_vector.begin ();
//        iterator != pressure_level_button_ptr_vector.end (); iterator++)
//   {
//      Level_Button* button_ptr = *(iterator);
//      pack_back (*button_ptr);
//   }
//
//   pack_back (screen_button);
//   pack ();

}

void
Level_Panel::move_level_up (const Devent& event)
{

   const Integer n = extra_level_vector.size ();
   const Integer i = get_extra_level_index ();

   if (i < 0 || i >= n)
   {
      try
      {
         const Level_Tuple& level_tuple = get_level_tuple ();
         level.value = level_tuple.get_next_up (level.value);
         switch (level.type)
         {
            case Level::PRESSURE: pressure = level.value; break;
            case Level::THETA: theta = level.value; break;
            case Level::SIGMA: sigma = level.value; break;
         }
      }
      catch (const Nwp_Exception& ne)
      {
      }
   }
   else
   {
      if (i == n - 1)
      {
         const Level_Tuple& level_tuple = get_level_tuple ();
         level = Level (level_tuple.type, level_tuple.back ());
      }
      else
      {
         level = extra_level_vector.at (i + 1);
      }
   }

   level_signal.emit (get_level ());
   full_level_signal.emit (get_level (), event);

}

void
Level_Panel::move_level_down (const Devent& event)
{

   const Integer n = extra_level_vector.size ();
   const Integer i = get_extra_level_index ();

   if (i < 0 || i >= n)
   {
      const Level_Tuple& level_tuple = get_level_tuple ();
      if (extra_level_vector.size () > 0 && level.value >= level_tuple.back ())
      {
         level = extra_level_vector.back ();
      }
      else
      {
         try
         {
            level.value = level_tuple.get_next_down (level.value);
            switch (level.type)
            {
               case Level::PRESSURE: pressure = level.value; break;
               case Level::THETA: theta = level.value; break;
               case Level::SIGMA: sigma = level.value; break;
            }
         }
         catch (const Nwp_Exception& ne)
         {
         }
      }
   }
   else
   {
      level = extra_level_vector.at (std::max (0, i - 1));
   }

   level_signal.emit (get_level ());
   full_level_signal.emit (get_level (), event);

}

Level_Canvas::Level_Canvas (Dcanvas& canvas,
                            const Real font_size)
   : canvas (canvas),
     level_panel (canvas, font_size)
{
}

const Level_Panel&
Level_Canvas::get_level_panel () const
{
   return level_panel;
}

Level_Panel&
Level_Canvas::get_level_panel ()
{
   return level_panel;
}

bool
Level_Canvas::on_key_pressed (const Dkey_Event& event)
{

   switch (event.value)
   {

      case GDK_KEY_Up:
      {

         level_panel.move_level_up (event);
         //const Level& level = level_panel.get_level ();
         //ink.set_level (level);
         return true;
         break;
      }

      case GDK_KEY_Down:
      {
         level_panel.move_level_down (event);
         //const Level& level = level_panel.get_level ();
         //ink.set_level (level);
         break;
      }

   }

   return false;

}

Console_2D::Manipulate_Data::Manipulate_Data ()
   : point (GSL_NAN, GSL_NAN),
     tilt (GSL_NAN)
{
}

Console_2D::Manipulate_Data::Manipulate_Data (const Point_2D& point,
                                              const Real tilt)
   : point (point),
     tilt (tilt)
{
}

Console_2D::Manipulate_Data::Manipulate_Data (const Manipulate_Data& md)
   : point (md.point),
     tilt (md.tilt)
{
}

void
Console_2D::Manipulate_Data::set_point (const Point_2D& point)
{
   this->point.x = point.x;
   this->point.y = point.y;
}

void
Console_2D::Manipulate_Data::reset ()
{
   tilt = GSL_NAN;
   point.x = GSL_NAN;
   point.y = GSL_NAN;
}

void
Console_2D::Manipulate_Data::set_tilt (const Real tilt)
{
   this->tilt = tilt;
}

bool
Console_2D::Manipulate_Data::is_rotating () const
{
   return !gsl_isnan (tilt);
}

bool
Console_2D::Manipulate_Data::is_translating () const
{
   return !point.is_nap ();
}

bool
Console_2D::Manipulate_Data::is_manipulating () const
{
   return (is_rotating () || is_translating ());
}

Console_2D::Zoom_Box::Scale_Ratio::Scale_Ratio (const Real fine,
                                                const Real coarse)
   : fine (fine),
     coarse (coarse)
{
}

Console_2D::Zoom_Box::Zoom_Box (Console_2D& console_2d)
   : console_2d (console_2d),
     color (1, 0, 0, 0.5),
     scale_ratio (1.01, 1.1),
     visible (false),
     centre (0, 0),
     zoom_factor_x (1.5),
     zoom_factor_y (1.5),
     tilt (0)
{
//   reset (true);
}

Console_2D::Zoom_Box::~Zoom_Box ()
{
}

void
Console_2D::Zoom_Box::set_color (const Color& color)
{
   this->color = color;
}

void
Console_2D::Zoom_Box::apply_to (Affine_Transform_2D& affine_transform) const
{

   // centre is in screen coordinates
   const Size_2D& size_2d = console_2d.get_size_2d ();
   const Real w2 = Real (size_2d.i) / 2;
   const Real h2 = Real (size_2d.j) / 2;
   const Real theta = atan2 (h2, w2);
   const Real w2zf = w2 / zoom_factor_x;
   const Real h2zf = h2 / zoom_factor_x;
   const Real r = sqrt (w2zf*w2zf + h2zf*h2zf);
   const Real phi = tilt + theta;

   const Real dx = (r * cos (phi) - centre.x);
   const Real dy = (r * sin (phi) - centre.y);

   affine_transform.translate (dx, dy);
   affine_transform.rotate (-tilt);
   affine_transform.scale (zoom_factor_x, zoom_factor_y);
}

void
Console_2D::Zoom_Box::pre_manipulate (const Dmouse_Button_Event& event)
{

   const Point_2D& point = event.point;

   if (!event.shift ()) { manipulate_data.set_point (point); } 
   else
   {
      const Real dx = point.x - centre.x;
      const Real dy = point.y - centre.y;
      manipulate_data.set_tilt (atan2 (dy, dx));
   }

}

void
Console_2D::Zoom_Box::manipulate (const Dmouse_Motion_Event& event)
{

   const Point_2D& point = event.point;

   if (manipulate_data.is_rotating ())
   {
      const Real dx = point.x - centre.x;
      const Real dy = point.y - centre.y;
      const Real theta = atan2 (dy, dx);
      const Real d_theta = theta - manipulate_data.tilt;
      tilt += d_theta;
      manipulate_data.tilt += d_theta;
   }

   if (manipulate_data.is_translating ())
   {
      centre.x += point.x - manipulate_data.point.x;
      centre.y += point.y - manipulate_data.point.y;
      manipulate_data.set_point (point);
   }

   console_2d.queue_draw ();

}

void
Console_2D::Zoom_Box::post_manipulate ()
{
   manipulate_data.reset ();
   console_2d.queue_draw ();
}

void
Console_2D::Zoom_Box::scroll (const Dmouse_Scroll_Event& event)
{

   const Real f = (event.control () ? scale_ratio.fine : scale_ratio.coarse);
   const bool up = (event.direction == GDK_SCROLL_UP);
   const Real zf = (up ? 1.0 / f : f);

   if (event.alt ())
   {
      this->zoom_factor_x /= zf;
      this->zoom_factor_y *= zf;
   }
   else
   {
      this->zoom_factor_x *= zf;
      this->zoom_factor_y *= zf;
   }

}

void
Console_2D::Zoom_Box::reset ()
{
   const Size_2D& size_2d = console_2d.get_size_2d ();
   this->centre.x = size_2d.i / 2;
   this->centre.y = size_2d.j / 2;
   this->zoom_factor_x = 1.5;
   this->zoom_factor_y = 1.5;
   this->tilt = 0;
}

void
Console_2D::Zoom_Box::toggle_visible ()
{
   visible = !visible;
}

void
Console_2D::Zoom_Box::set_visible (const bool visible)
{
   this->visible = visible;
}

bool
Console_2D::Zoom_Box::is_visible () const
{
   return this->visible;
}

bool
Console_2D::Zoom_Box::contains (const Point_2D& point) const
{

   const Size_2D& size_2d = console_2d.get_size_2d ();
   const Real width = Real (size_2d.i / zoom_factor_x);
   const Real height = Real (size_2d.j / zoom_factor_y);
   const Real w = width / 2;
   const Real h = height / 2;
   const Real phi = atan2 (height, width);
   const Real theta_a = tilt + phi;
   const Real theta_b = M_PI - phi + tilt;
   const Real theta_c = M_PI + phi + tilt;
   const Real theta_d = -phi + tilt;

   const Real r = sqrt (w*w + h*h);
   const Real x = centre.x;
   const Real y = centre.y;
   const Real xa = x + r * cos (theta_a), ya = y + r * sin (theta_a);
   const Real xb = x + r * cos (theta_b), yb = y + r * sin (theta_b);
   const Real xc = x + r * cos (theta_c), yc = y + r * sin (theta_c);
   const Real xd = x + r * cos (theta_d), yd = y + r * sin (theta_d);

   Polygon polygon;
   polygon.add (Point_2D (xa, ya));
   polygon.add (Point_2D (xb, yb));
   polygon.add (Point_2D (xc, yc));
   polygon.add (Point_2D (xd, yd));

   return polygon.contains (point);

}

void
Console_2D::Zoom_Box::cairo (const RefPtr<Context>& cr)
{

   if (!visible) { return; }

   const Real font_size = 12;
   const Size_2D& size_2d = console_2d.get_size_2d ();
   const Real width = Real (size_2d.i / zoom_factor_x);
   const Real height = Real (size_2d.j / zoom_factor_y);
   const Real w = width / 2;
   const Real h = height / 2;
   const Real phi = atan2 (height, width);
   const Real theta_a = phi + tilt;
   const Real theta_b = M_PI - phi + tilt;
   const Real theta_c = M_PI + phi + tilt;
   const Real theta_d = -phi + tilt;

   const Real r = sqrt (w*w + h*h);
   const Real wa = r * cos (theta_a), ha = r * sin (theta_a);
   const Real wb = r * cos (theta_b), hb = r * sin (theta_b);
   const Real wc = r * cos (theta_c), hc = r * sin (theta_c);
   const Real wd = r * cos (theta_d), hd = r * sin (theta_d);

   const Point_2D point_2a (centre.x + wa, centre.y + ha);
   const Point_2D point_2b (centre.x + wb, centre.y + hb);
   const Point_2D point_2c (centre.x + wc, centre.y + hc);
   const Point_2D point_2d (centre.x + wd, centre.y + hd);

   cr->save ();
   color.cairo (cr);

   Dashes ("4:2").cairo (cr);

   cr->set_line_width (1);
   cr->move_to (centre.x + wa, centre.y + ha);
   cr->line_to (centre.x + wc, centre.y + hc);
   cr->stroke ();
   cr->move_to (centre.x + wb, centre.y + hb);
   cr->line_to (centre.x + wd, centre.y + hd);
   cr->stroke ();

   cr->move_to (centre.x + wa / 2, centre.y + ha / 2);
   cr->line_to (centre.x + wb / 2, centre.y + hb / 2);
   cr->line_to (centre.x + wc / 2, centre.y + hc / 2);
   cr->line_to (centre.x + wd / 2, centre.y + hd / 2);
   cr->close_path ();
   cr->stroke ();

   cr->unset_dash ();

   cr->set_line_width (4);
   cr->set_line_cap (LINE_CAP_ROUND);
   cr->move_to (point_2a.x, point_2a.y);
   cr->line_to (point_2b.x, point_2b.y);
   cr->line_to (point_2c.x, point_2c.y);
   cr->line_to (point_2d.x, point_2d.y);
   cr->close_path ();
   cr->stroke ();

   cr->set_line_width (2);
   Ring (6).cairo (cr, point_2c);
   cr->fill ();
   cr->set_font_size (font_size);
   Label label ("TOP LEFT", point_2c, 'l', 't', font_size / 2);
   label.set_text_angle (tilt);
   label.cairo (cr);

   cr->restore ();

}

void
Console_2D::Hud::Popup_Menu::emit_signal (const Dstring& str) const
{
   if (index == 0) { clear_signal.emit (id); }
   str_signal_map.at (str).emit (str);
   id_signal_map.at (str).emit (id);
}

Console_2D::Hud::Popup_Menu::Popup_Menu (Console_2D& console_2d,
                                         const Real font_size)
   : denise::Popup_Menu (console_2d, font_size),
     console_2d (console_2d)
{
   clear ();
   append ("Clear");
}

void
Console_2D::Hud::Popup_Menu::clear ()
{
   denise::Popup_Menu::clear ();
   id_signal_map.clear ();
}

void
Console_2D::Hud::Popup_Menu::append (const Dstring& str)
{
   denise::Popup_Menu::append (str);
   Console_2D::Hud::Popup_Menu::Id_Signal id_signal;
   id_signal_map.insert (make_pair (str, id_signal));
}

Console_2D::Hud::Popup_Menu::Id_Signal&
Console_2D::Hud::Popup_Menu::get_id_signal (const Dstring& str)
{
   return id_signal_map.at (str);
}

void
Console_2D::Hud::Popup_Menu::setup (const Hud& hud,
                                    const Point_2D& point)
{

   this->id = hud.get_id ();

   const Real width = get_preferred_width ();
   const Real height = get_preferred_height ();
   const Point_2D anchor = point + Point_2D (6, 6);
   Popup_Menu::set_shape (anchor, orientation, width, height);

}

Console_2D::Hud::Hud (const Integer id,
                      const Real node_size)
   : id (id),
     node_size (node_size)
{
}

Console_2D::Hud::~Hud ()
{
}

const Integer
Console_2D::Hud::get_id () const
{
   return id;
}

sigc::signal<void>&
Console_2D::Hud::get_moved_signal ()
{
   return moved_signal;
}

sigc::signal<void>&
Console_2D::Hud::get_settled_signal ()
{
   return settled_signal;
}

void
Console_2D::Hud::double_clicked ()
{
}

Console_2D::Marker::Popup_Menu::Popup_Menu (Console_2D& console_2d,
                                            const Real font_size,
                                            const bool connect_default_signal)
   : Console_2D::Hud::Popup_Menu (console_2d, font_size)
{
   if (connect_default_signal)
   {
      Console_2D& c = console_2d;
      clear_signal.connect (sigc::mem_fun (c, &Console_2D::clear_marker));
   }
}

Console_2D::Marker::Marker (const Integer id,
                            const Point_2D& point,
                            const Real node_size)
   : Hud (id, node_size),
     Point_2D (point),
     str ("")
{
}

void
Console_2D::Marker::attract_by (const Attractor& attractor)
{
   attractor.attract (*this);
}

const Dstring&
Console_2D::Marker::get_str () const
{
   return str;
}

void
Console_2D::Marker::set_str (const Dstring& str)
{
   this->str = str;
}

bool
Console_2D::Marker::matches (const Transform_2D& transform,
                             const Point_2D& point) const
{
   const Point_2D& p = transform.transform (*this);
   const Real dx = p.x - point.x;
   const Real dy = p.y - point.y;
   return (dx * dx + dy * dy < node_size * node_size);
}

void
Console_2D::Marker::cairo (const RefPtr<Context>& cr,
                           const Console_2D& console_2d) const
{

   const Transform_2D& transform = console_2d.get_transform ();
   const Point_2D& point = transform.transform (*this);

   cr->save ();

   if (console_2d.get_tokens (*this).size () > 0)
   {

      Popup popup (console_2d, 12);
      popup.set_tokens (console_2d.get_tokens (*this));
      const Real popup_w = popup.get_preferred_width ();
      const Real popup_h = popup.get_preferred_height ();
      const Size_2D& size_2d = console_2d.get_size_2d ();

      Point_2D off (6, 6);
      Orientation o = ORIENTATION_SE;

      const Real x = point.x;
      const Real y = point.y;
      const Real w = size_2d.i;
      const Real h = size_2d.j;

      const bool over_x = (x + popup_w + 2 * off.x) > w;
      const bool over_y = (y + popup_h + 2 * off.y) > h;

      if (over_x && over_y) { o = ORIENTATION_NW; off.x *= -1; off.y *= -1; }
      else if (over_x) { o = ORIENTATION_SW; off.x *= -1; }
      else if (over_y) { o = ORIENTATION_NE; off.y *= -1; }

      const Point_2D anchor (x + off.x, y + off.y);
      popup.set_shape (anchor, o, popup_w, popup_h);
      popup.cairo (cr);

   }

   const Real hue = (str == "" ? 0.667 : 0.000);
   const Color& bg_color = Color::hsb (hue, 0.5, 1, 0.7);
   const Color fg_color (0, 0, 0);
   Ring (node_size).cairo (cr, point);
   bg_color.cairo (cr);
   cr->fill_preserve ();
   fg_color.cairo (cr);
   cr->stroke ();

   cr->restore ();

}

Console_2D::Route::Popup_Menu::Popup_Menu (Console_2D& console_2d,
                                           const Real font_size,
                                           const bool connect_default_signal)
   : Console_2D::Hud::Popup_Menu (console_2d, font_size)
{
   if (connect_default_signal)
   {
      Console_2D& c = console_2d;
      clear_signal.connect (sigc::mem_fun (c, &Console_2D::clear_route));
   }
}

Console_2D::Route::Route (const Integer id,
                          const Real node_size)
   : Hud (id, node_size),
     origin (GSL_NAN, GSL_NAN)
{
}

Console_2D::Route::Route (const Integer id,
                          const Point_2D& point,
                          const Real node_size)
   : Hud (id, node_size),
     origin (GSL_NAN, GSL_NAN)
{
   push_back (point);
   push_back (point);
}

Console_2D::Route::Route (const Integer id,
                          const Point_2D& point_a,
                          const Point_2D& point_b,
                          const Real node_size)
   : Hud (id, node_size),
     origin (GSL_NAN, GSL_NAN)
{
   push_back (point_a);
   push_back (point_b);
}

bool
Console_2D::Route::is_too_short () const
{
   return (size () < 2);
}

void
Console_2D::Route::translate (const Point_2D& from,
                              const Point_2D& to)
{
   const Real dx = to.x - from.x;
   const Real dy = to.y - from.y;
   Simple_Polyline::translate (dx, dy);
}

void
Console_2D::Route::translate (const Point_2D& to)
{
   const Real dx = to.x - origin.x;
   const Real dy = to.y - origin.y;
   Simple_Polyline::translate (dx, dy);
}

void
Console_2D::Route::set_origin (const Point_2D& point)
{
   this->origin.x = point.x;
   this->origin.y = point.y;
}

const Point_2D&
Console_2D::Route::get_origin () const
{
   return origin;
}

bool
Console_2D::Route::has_origin () const
{
   return (!origin.is_nap ());
}

bool
Console_2D::Route::matches (const Transform_2D& transform,
                            const Point_2D& point) const
{
   return (get_iterator (transform, point, node_size) != end ());

/*
   Simple_Polyline sp;
 
   for (Simple_Polyline::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Point_2D& p = transform.transform (*(iterator));
      sp.push_back (p);
   }

   return (sp.get_iterator (point, node_size) != sp.end ());
*/

}

void
Console_2D::Route::cairo (const RefPtr<Context>& cr,
                          const Console_2D& console_2d) const
{

   cr->save ();
   cr->set_line_cap (LINE_CAP_ROUND);
   cr->set_line_join (LINE_JOIN_ROUND);

   const Real node_size = 16;
   const Color& bg_color = Color::hsb (0.167, 0.2, 0.5, 0.7);
   const Color& fg_color = Color::hsb (0.167, 0.2, 1.0, 1.0);

   const Transform_2D& transform = console_2d.get_transform ();

   cr->set_line_width (node_size);
   bg_color.cairo (cr);
   Simple_Polyline::cairo (cr, transform);
   cr->stroke ();

   fg_color.cairo (cr);
   cr->set_line_width (2);
   for (Simple_Polyline::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Point_2D& p = transform.transform (*(iterator));
      Ring (node_size/2).cairo (cr, p);
      cr->stroke ();
   }

   cr->restore ();

}

Console_2D::Route::iterator
Console_2D::Route::add (const Transform_2D& transform,
                        const Point_2D& point)
{
   return Simple_Polyline::implant (transform, point, node_size);
}

Console_2D::Route::const_iterator
Console_2D::Route::get_node (const Transform_2D& transform,
                             const Point_2D& point) const
{

   const Real threshold = node_size * node_size;

   for (Console_2D::Route::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Point_2D& p = transform.transform (*(iterator));
      const Real dx = p.x - point.x;
      const Real dy = p.y - point.y;
      if (dx * dx + dy * dy < threshold) { return iterator; }
   }

   return end ();

}

Console_2D::Route::iterator
Console_2D::Route::get_node (const Transform_2D& transform,
                             const Point_2D& point)
{

   const Real threshold = node_size * node_size;

   for (Console_2D::Route::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Point_2D& p = transform.transform (*(iterator));
      const Real dx = p.x - point.x;
      const Real dy = p.y - point.y;
      if (dx * dx + dy * dy < threshold) { return iterator; }
   }

   return end ();

}

Polygon*
Console_2D::Shape::get_polygon_ptr () const
{

   Polygon* polygon_ptr = new Polygon ();

   for (Simple_Polyline::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Point_2D& p = *(iterator);
      polygon_ptr->add (p);
   }

   return polygon_ptr;

}

bool
Console_2D::Shape::contains (const Point_2D& point) const
{
   Polygon* polygon_ptr = get_polygon_ptr ();
   const bool b = polygon_ptr->contains (point, true);
   delete polygon_ptr;
   return b;
}

bool
Console_2D::Shape::contains (const Transform_2D& transform,
                             const Point_2D& point) const
{
   Polygon* polygon_ptr = get_polygon_ptr ();
   const bool b = polygon_ptr->contains (transform.reverse (point), true);
   delete polygon_ptr;
   return b;
}

bool
Console_2D::Shape::on_node (const Transform_2D& transform,
                            const Point_2D& point) const
{
   return (Console_2D::Route::get_node (transform, point) != end ());
}

Console_2D::Shape::Popup_Menu::Popup_Menu (Console_2D& console_2d,
                                           const Real font_size,
                                           const bool connect_default_signal)
   : Console_2D::Hud::Popup_Menu (console_2d, font_size)
{
   if (connect_default_signal)
   {
      Console_2D& c = console_2d;
      clear_signal.connect (sigc::mem_fun (c, &Console_2D::clear_shape));
   }
}

Console_2D::Shape::Shape (const Integer id,
                          const Real node_size)
   : Route (id, node_size)
{
   set_closed (true);
}

Console_2D::Shape::Shape (const Integer id,
                          const Point_2D& point,
                          const Real node_size)
   : Route (id, node_size)
{

   const Real x = point.x;
   const Real y = point.y;
   push_back (Point_2D (x + 4, y + 4));
   push_back (Point_2D (x + 4, y - 4));
   push_back (Point_2D (x - 4, y - 4));
   push_back (Point_2D (x - 4, y + 4));

   set_closed (true);

}

bool
Console_2D::Shape::is_too_short () const
{
   return (size () < 3);
}

bool
Console_2D::Shape::on_edge (const Transform_2D& transform,
                            const Point_2D& point) const
{
   return Console_2D::Route::matches (transform, point);
}

bool
Console_2D::Shape::matches (const Transform_2D& transform,
                            const Point_2D& point) const
{
   const bool m = on_edge (transform, point);
   const bool ii = contains (transform, point);
   const bool on = on_node (transform, point);
   return (m || ii || on);
}

void
Console_2D::Shape::cairo (const RefPtr<Context>& cr,
                          const Console_2D& console_2d) const
{

   const Real dash_width = 10;
   const Dashes dashes (Tuple (2, dash_width));
   const Transform_2D& transform = console_2d.get_transform ();

   const Color& bg_a = Color::hsb (0.667, 0.2, 0.8, 0.4);
   const Color& bg_b = Color::hsb (0.900, 0.2, 0.8, 0.4);
   const Color& fg_a = Color::hsb (0.667, 0.4, 0.4, 0.8);
   const Color& fg_b = Color::hsb (0.900, 0.4, 0.4, 0.8);

   cr->save ();
   cr->set_line_cap (LINE_CAP_BUTT);
   cr->set_line_join (LINE_JOIN_ROUND);
   Stripped (bg_a, bg_b, dash_width).cairo (cr);

   Polygon* polygon_ptr = get_polygon_ptr ();
   polygon_ptr->cairo (cr, transform);
   cr->fill_preserve ();
   delete polygon_ptr;

   cr->set_line_width (4);
   dashes.cairo (cr);
   fg_a.cairo (cr);
   cr->stroke_preserve ();
   dashes.cairo (cr);
   fg_b.cairo (cr);
   cr->stroke ();

   Stripped (fg_a, fg_b, 2).cairo (cr);

   for (Console_2D::Shape::const_iterator i = begin (); i != end (); i++)
   {
      const Point_2D& p = transform.transform (*i);
      Ring (node_size).cairo (cr, p);
      cr->fill ();
   }

   cr->restore ();

}

Integer
Console_2D::Hud_Store::get_first_available_id () const
{
   for (Integer candidate_id = 0; true; candidate_id++)
   {
      if (find (candidate_id) == end ()) { return candidate_id; }
   }
}

Console_2D::Hud_Store::Hud_Store ()
   : dragging_id (-1),
     attractor_ptr (NULL)
{
}

Console_2D::Hud_Store::~Hud_Store ()
{
   for (Hud_Store::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Hud* hud_ptr = iterator->second;
      delete hud_ptr;
   }
   std::map<Integer, Console_2D::Hud*>::clear ();
}

void
Console_2D::Hud_Store::set_attractor (const Attractor& attractor)
{
   this->attractor_ptr = &attractor;
}

const Attractor&
Console_2D::Hud_Store::get_attractor () const
{
   if (attractor_ptr != NULL) { return *attractor_ptr; }
   else throw Exception ("attractor_ptr is NULL.");
}

void
Console_2D::Hud_Store::clear ()
{
   for (Hud_Store::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      Hud* hud_ptr = iterator->second;
      delete hud_ptr;
   }
   std::map<Integer, Hud*>::clear ();
}

void
Console_2D::Hud_Store::reset ()
{
   dragging_id = -1;
}

const Console_2D::Hud&
Console_2D::Hud_Store::get_hud (const Integer id) const
{
   return *(find (id)->second);
}

Console_2D::Hud&
Console_2D::Hud_Store::get_hud (const Integer id)
{
   return *(find (id)->second);
}

Integer
Console_2D::Hud_Store::get_id (const Transform_2D& transform,
                               const Point_2D& point)
{

   for (Marker_Store::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Integer& id = iterator->first;
      const Hud& hud = (*(iterator->second));
      if (hud.matches (transform, point)) { return id; }
   }

   return -1;

}

void
Console_2D::Hud_Store::insert (const Integer id,
                               Hud* hud_ptr)
{
   std::map<Integer, Hud*>::insert (make_pair (id, hud_ptr));
}

void
Console_2D::Hud_Store::erase (const Integer id)
{
   delete (find (id)->second);
   std::map<Integer, Console_2D::Hud*>::erase (id);
}

void
Console_2D::Hud_Store::cairo (const RefPtr<Context>& cr,
                              const Console_2D& console_2d) const
{

   cr->save ();
   cr->unset_dash ();
   cr->set_fill_rule (FILL_RULE_EVEN_ODD);

   for (Hud_Store::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Hud& hud = (*(iterator->second));
      hud.cairo (cr, console_2d);
   }

   cr->restore ();

}

bool
Console_2D::Hud_Store::show_popup_menu (Console_2D& console_2d,
                                        const Point_2D& point,
                                        Hud::Popup_Menu& popup_menu)
{

   const Transform_2D& transform = console_2d.get_transform ();
   const Integer id = get_id (transform, point);

   if (id == -1) { return false; }

   // show pop up for the hud
   Hud& hud = get_hud (id);
   popup_menu.setup (hud, point);
   console_2d.queue_draw ();
   return true;

}

Console_2D::Marker_Store::Marker_Store ()
   : adding_marker (false)
{
}

Console_2D::Marker_Store::~Marker_Store ()
{
}

void
Console_2D::Marker_Store::reset ()
{
   Console_2D::Hud_Store::reset ();
   adding_marker = false;
}

const Console_2D::Marker&
Console_2D::Marker_Store::get_marker (const Integer id) const
{
   return dynamic_cast<const Marker&>(*(find (id)->second));
}

Console_2D::Marker&
Console_2D::Marker_Store::get_marker (const Integer id)
{
   return dynamic_cast<Marker&>(*(find (id)->second));
}

Integer
Console_2D::Marker_Store::insert (const Point_2D& point)
{
   const Integer& id = get_first_available_id ();
   Marker* marker_ptr = new Marker (id, point);
   Hud_Store::insert (id, marker_ptr);
   return id;
}

bool                                                               
Console_2D::Marker_Store::button_1_pressed (Console_2D& console_2d,
                                            const Dmouse_Button_Event& event)
{

   const Transform_2D& transform = console_2d.get_transform ();
   const Point_2D& point = event.point;
   const Integer id = get_id (transform, point);

   // not clicked on the route, pass
   if (id == -1) { return false; }

   this->dragging_id = id;
   console_2d.queue_draw ();
   return true;

}          
                                                                                          
bool                                                               
Console_2D::Marker_Store::try_insert (Console_2D& console_2d,
                                      const Dmouse_Button_Event& event)
{

   const Transform_2D& transform = console_2d.get_transform ();
   const Point_2D& point = event.point;
   const Integer id = get_id (transform, point);

   // clicked on the marker, pass
   if (id != -1) { return false; }

   // add a new marker
   const Point_2D& p = transform.reverse (point);
   dragging_id = insert (p);
   adding_marker = true;
   console_2d.refresh_all ();
   return true;

}          
               
bool                                                               
Console_2D::Marker_Store::button_1_released (Console_2D& console_2d,
                                             const Dmouse_Button_Event& event)
{
   if (dragging_id == -1) { return false; }
   get_marker (dragging_id).get_settled_signal ().emit ();
   reset ();
   return true;
}          
                                                                                          
bool                                                               
Console_2D::Marker_Store::button_3_released (Console_2D& console_2d,
                                             const Dmouse_Button_Event& event)
{
   if (!adding_marker) { return false; }
   get_marker (dragging_id).get_settled_signal ().emit ();
   reset ();
   return true;
}          

bool
Console_2D::Marker_Store::try_dragging (Console_2D& console_2d,
                                        const Dmouse_Motion_Event& event)
{

   if (dragging_id == -1) { return false; }

   Marker& marker = get_marker (dragging_id);
   const Transform_2D& transform = console_2d.get_transform ();
   const Point_2D& point = event.point;

   transform.reverse (marker.x, marker.y, point.x, point.y);
   marker.set_str ("");

   if (event.control ())
   {
      try
      {
         const Attractor& attractor = get_attractor ();
         marker.attract_by (attractor);
      }
      catch (const Exception& e) { }
   }

   marker.get_moved_signal ().emit ();
   console_2d.refresh_all ();
   return true;

}                    

Console_2D::Route*
Console_2D::Route_Store::new_route_ptr (const Integer id,
                                        const Point_2D& point,
                                        const Real node_size)
{
   return new Console_2D::Route (id, point, node_size);
}

Console_2D::Route*
Console_2D::Route_Store::new_route_ptr (const Integer id,
                                        const Point_2D& point_a,
                                        const Point_2D& point_b,
                                        const Real node_size)
{
   return new Console_2D::Route (id, point_a, point_b, node_size);
}

Console_2D::Route_Store::Route_Store ()
   : dummy (-2),
     node (dummy.end ())
{
}

Console_2D::Route_Store::~Route_Store ()
{
}

void
Console_2D::Route_Store::reset ()
{
   Console_2D::Hud_Store::reset ();
   node = dummy.end ();
}

const Console_2D::Route&
Console_2D::Route_Store::get_route (const Integer id) const
{
   return dynamic_cast<const Route&>(*(find (id)->second));
}

Console_2D::Route&
Console_2D::Route_Store::get_route (const Integer id)
{
   return dynamic_cast<Route&>(*(find (id)->second));
}

Console_2D::Route::iterator
Console_2D::Route_Store::get_node (const Transform_2D& transform,
                                   const Point_2D& point)
{

   for (Route_Store::iterator i = begin (); i != end (); i++)
   {
      Route& route = dynamic_cast<Route&>(*(i->second));
      Console_2D::Route::iterator node = route.get_node (transform, point);
      if (node != route.end ()) { return node; }
   }

   return dummy.end ();

}

Integer
Console_2D::Route_Store::insert (const Point_2D& point)
{

   const Integer& id = get_first_available_id ();
   Route* route_ptr = this->new_route_ptr (id, point);
   Hud_Store::insert (id, route_ptr);

   // node should still be modified because node is not truely defined yet
   node = (--route_ptr->end ());

   return id;

}

Integer
Console_2D::Route_Store::insert (const Point_2D& point_a,
                                 const Point_2D& point_b)
{

   const Integer& id = get_first_available_id ();
   Route* route_ptr = this->new_route_ptr (id, point_a, point_b);
   Hud_Store::insert (id, route_ptr);

   // node should still be dummy.end () becasue the node is inserted already

   return id;

}

bool
Console_2D::Route_Store::button_1_pressed (Console_2D& console_2d,
                                           const Dmouse_Button_Event& event)
{

   const Transform_2D& transform = console_2d.get_transform ();
   const Point_2D& point = event.point;
   const Integer id = get_id (transform, point);

   // not clicked on the route, pass
   if (id == -1) { return false; }

   this->dragging_id = id;
   Route& route = get_route (id);
   node = route.get_node (transform, point);

   const bool on_node = (node != route.end ());

   if (on_node)
   {
      // clicked on a route node
      transform.reverse (node->x, node->y, point.x, point.y);
      console_2d.queue_draw ();
   }
   else
   {
      // clicked on a route but not a node
      route.set_origin (transform.reverse (point));
   }

   return true;

}

bool
Console_2D::Route_Store::button_2_pressed (Console_2D& console_2d,
                                           const Dmouse_Button_Event& event)
{

   const Transform_2D& transform = console_2d.get_transform ();
   const Point_2D& point = event.point;

   const Integer id = get_id (transform, point);

   // add a new route
   if (id == -1)
   {
      const Point_2D& p = transform.reverse (point);
      dragging_id = insert (p);
      console_2d.refresh_all ();
      return true;
   }

   dragging_id = id;
   Route& route = get_route (id);
   node = route.get_node (transform, point);

   // an existing route node, erase this node
   if (node != route.end ())
   {
      route.erase (node);
      if (route.is_too_short ()) { erase (id); reset (); }
      console_2d.refresh_all ();
      return true;
   }

   // a new node is added to the route
   node = route.add (transform, point);
   console_2d.refresh_all ();
   return true;

}

bool
Console_2D::Route_Store::button_1_released (Console_2D& console_2d,
                                            const Dmouse_Button_Event& event)
{

   if (dragging_id != -1)
   {
      Route& route = get_route (dragging_id);
      if (route.has_origin ())
      {
         route.set_origin (Point_2D (GSL_NAN, GSL_NAN));
         route.get_settled_signal ().emit ();
         reset ();
         return true;
      }
   }

   if (node != dummy.end ())
   {
      get_route (dragging_id).get_settled_signal ().emit ();
      reset ();
      return true;
   }

   return false;

}

bool
Console_2D::Route_Store::button_2_released (Console_2D& console_2d,
                                            const Dmouse_Button_Event& event)
{

   if (node != dummy.end ())
   {
      Route& route = get_route (dragging_id);
      route.get_settled_signal ().emit ();
      reset ();
      return true;
   }

   return false;

}

bool
Console_2D::Route_Store::try_dragging (Console_2D& console_2d,
                                       const Dmouse_Motion_Event& event)
{

   if (dragging_id == -1) { return false; }

   Route& route = get_route (dragging_id);
   const Transform_2D& transform = console_2d.get_transform ();
   const Point_2D& point = event.point;

   // moving a route node
   if (node != route.end ())
   {

      Point_2D& p = *(node);
      transform.reverse (p.x, p.y, point.x, point.y);

      if (event.control ())
      {
         try
         {
            const Attractor& attractor = get_attractor ();
            attractor.attract (p);
         }
         catch (const Exception& e) { }
      }

      route.get_moved_signal ().emit ();
      console_2d.refresh_all ();
      return true;

   }

   // route is grabbed and is moved
   if (route.has_origin ())
   {

      const Point_2D& p = transform.reverse (point);
      route.translate (p);

      if (event.control ())
      {
         try
         {
            const Attractor& attractor = get_attractor ();
            route.attract_by (attractor);
         }
         catch (const Exception& e) { }
      }

      route.set_origin (p);
      route.get_moved_signal ().emit ();
      console_2d.refresh_all ();
      return true;

   }

   return false;

}

Console_2D::Shape_Store::Shape_Store ()
   : Console_2D::Route_Store ()
{
}

const Console_2D::Shape&
Console_2D::Shape_Store::get_shape (const Integer id) const
{
   return dynamic_cast<const Shape&>(*(find (id)->second));
}

Console_2D::Shape&
Console_2D::Shape_Store::get_shape (const Integer id)
{
   return dynamic_cast<Shape&>(*(find (id)->second));
}

Integer
Console_2D::Shape_Store::insert (const Point_2D& point)
{
   const Integer& id = get_first_available_id ();
   Shape* shape_ptr = new Shape (id, point);
   Hud_Store::insert (id, shape_ptr);
   return id;
}

bool
Console_2D::Shape_Store::button_1_pressed (Console_2D& console_2d,
                                           const Dmouse_Button_Event& event)
{

   const Transform_2D& transform = console_2d.get_transform ();
   const Point_2D& point = event.point;
   const Integer id = get_id (transform, point);
   const bool double_clicked = (event.type == GDK_2BUTTON_PRESS);

   // not clicked on the shape
   if (id == -1)
   {

      // not double clicked either, pass
      if (!double_clicked) { return false; }

      // double clicked on thin air, initiate a new shape
      transform.reverse (node->x, node->y, point.x, point.y);
      insert (*node);
      console_2d.refresh_all ();
      return true;

   }

   this->dragging_id = id;
   Shape& shape = get_shape (id);
   node = shape.get_node (transform, point);

   const bool on_node = (node != shape.end ());
   const bool on_edge = shape.on_edge (transform, point);

   if (on_node)
   {

      // clicked on a node
      if (double_clicked)
      {
         // an existing route node, erase this node
         shape.erase (node);
         if (shape.is_too_short ()) { erase (id); reset (); }
         console_2d.refresh_all ();
         return true;
      }
      else
      {
         // Grab the node
         transform.reverse (node->x, node->y, point.x, point.y);
         console_2d.queue_draw ();
         return true;
      }

   }
   else
   {

      if (on_edge)
      {
         node = shape.add (transform, point);
         // a new node is added to the route
         if (node != shape.end ())
         {
            console_2d.refresh_all ();
            return true;
         }
      }

      if (double_clicked) { /* shape.double_clicked (); */  }
      else { shape.set_origin (transform.reverse (point)); }
      return true;

   }

}

bool
Console_2D::Shape_Store::button_1_released (Console_2D& console_2d,
                                            const Dmouse_Button_Event& event)
{

   if (dragging_id != -1)
   {
      Shape& shape = get_shape (dragging_id);
      if (shape.has_origin ())
      {
         shape.set_origin (Point_2D (GSL_NAN, GSL_NAN));
         shape.get_settled_signal ().emit ();
         reset ();
         return true;
      }
   }

   if (node != dummy.end ())
   {
      get_shape (dragging_id).get_settled_signal ().emit ();
      reset ();
      return true;
   }

   return false;

}

const Console_2D::Zoom_Box&
Console_2D::get_zoom_box () const
{
   return zoom_box;
}

Console_2D::Zoom_Box&
Console_2D::get_zoom_box ()
{
   return zoom_box;
}

void
Console_2D::reset_transform ()
{

   affine_transform.set_identity ();

   if (!computer_convention)
   {
      const Size_2D& size_2d = get_size_2d ();
      affine_transform.scale (1, -1);
      affine_transform.translate (0, size_2d.j);
   }

   set_background_ready (false);
   set_foreground_ready (false);

}

void
Console_2D::pre_manipulate (const Dmouse_Button_Event& event)
{

   ignore_heavy = true;
   const Point_2D& point = event.point;

   const Size_2D& size_2d = get_size_2d ();
   const Real w2 = Real (size_2d.i) / 2;
   const Real h2 = Real (size_2d.j) / 2;

   if (!event.shift ()) { manipulate_data.set_point (point); } 
   else
   {
      const Real dx = point.x - w2;
      const Real dy = point.y - h2;
      manipulate_data.set_tilt (atan2 (dy, dx));
   }
}

void
Console_2D::manipulate (const Dmouse_Motion_Event& event)
{

   ignore_heavy = true;
   const Point_2D& point = event.point;

   const Size_2D& size_2d = get_size_2d ();
   const Real w2 = Real (size_2d.i) / 2;
   const Real h2 = Real (size_2d.j) / 2;

   if (manipulate_data.is_rotating ())
   {
      const Real dx = point.x - w2;
      const Real dy = point.y - h2;
      const Real theta = atan2 (dy, dx);
      const Real d_theta = theta - manipulate_data.tilt;
      affine_transform.translate (-w2, -h2);
      affine_transform.rotate (d_theta);
      affine_transform.translate (w2, h2);
      manipulate_data.set_tilt (theta);
   }

   if (manipulate_data.is_translating ())
   {
      const Real dx = point.x - manipulate_data.point.x;
      const Real dy = point.y - manipulate_data.point.y;
      translate (dx, dy);
      manipulate_data.set_point (point);
   }

   set_background_ready (false);
   set_foreground_ready (false);
   render_queue_draw ();

}


void
Console_2D::post_manipulate ()
{
   ignore_heavy = false;
   manipulate_data.reset ();
   set_background_ready (false);
   set_foreground_ready (false);
   render_queue_draw ();
}

bool
Console_2D::on_key_pressed (const Dkey_Event& event)
{

   Console_2D::Zoom_Box& zoom_box = get_zoom_box ();

   switch (event.value)
   {

      case GDK_KEY_h:
      case GDK_KEY_H:
      {
         const bool b = !get_hide_hidable ();
         set_hide_hidable (b);
         return true;
         break;
      }

      case GDK_KEY_x:
      case GDK_KEY_X:
      {
         get_marker_store ().clear ();
         get_route_store ().clear ();
         get_shape_store ().clear ();
         refresh_all ();
         return true;
         break;
      }

      case GDK_KEY_z:
      case GDK_KEY_Z:
      {
         zoom_box.toggle_visible ();
         queue_draw ();
         return true;
         break;
      }

      case GDK_KEY_BackSpace:
      {

         if (event.control ())
         {
            reset_transform ();
            render_queue_draw ();
            return true;
         }

         zoom_box.reset ();
         queue_draw ();
         return true;
         break;
      }

   }

   return Dcanvas::on_key_pressed (event);

}

bool
Console_2D::on_mouse_button_pressed (const Dmouse_Button_Event& event)
{

   if (Dcontainer::on_mouse_button_pressed (event)) { return true; }

   const Point_2D& point = event.point;
   const bool double_clicked = (event.type == GDK_2BUTTON_PRESS);

   Marker_Store& marker_store = get_marker_store ();
   Route_Store& route_store = get_route_store ();
   Shape_Store& shape_store = get_shape_store ();

   const Transform_2D& transform = get_transform ();
   Console_2D::Zoom_Box& zoom_box = get_zoom_box ();

   switch (event.button)
   {

      case 1:
      {

         if (zoom_box.is_visible () && zoom_box.contains (point))
         {
            if (double_clicked)
            {
               apply_zoom_box ();
               return true;
            }

            zoom_box.pre_manipulate (event);
            return true;

         }

         if (marker_store.button_1_pressed (*this, event)) { return true; }
         if (route_store.button_1_pressed (*this, event)) { return true; }
         if (shape_store.button_1_pressed (*this, event)) { return true; }

         if (event.control ())
         {
            pre_manipulate (event);
            return true;
         }

         break;

      }

      case 2:
      {
         if (route_store.button_2_pressed (*this, event)) { return true; }
         break;
      }

      case 3:
      {

         Marker::Popup_Menu& mpm = get_marker_popup_menu ();
         if (marker_store.show_popup_menu (*this, point, mpm)) { return true; }

         Route::Popup_Menu& rpm = get_route_popup_menu ();
         if (route_store.show_popup_menu (*this, point, rpm)) { return true; }

         Shape::Popup_Menu& spm = get_shape_popup_menu ();
         if (shape_store.show_popup_menu (*this, point, spm)) { return true; }

         if (marker_store.try_insert (*this, event)) { return true; }
         break;

      }

   }

   return false;

}

bool
Console_2D::on_mouse_button_released (const Dmouse_Button_Event& event)
{

   if (Dcontainer::on_mouse_button_released (event)) { return true; }

   const Point_2D& point = event.point;
   Marker_Store& marker_store = get_marker_store ();
   Route_Store& route_store = get_route_store ();
   Shape_Store& shape_store = get_shape_store ();

   Console_2D::Zoom_Box& zoom_box = get_zoom_box ();

   switch (event.button)
   {

      case 1:
      {

         if (zoom_box.manipulate_data.is_manipulating ())
         {
            zoom_box.post_manipulate ();
            return true;
         }

         if (manipulate_data.is_manipulating ())
         {
            post_manipulate ();
            return true;
         } 

         if (marker_store.button_1_released (*this, event)) { return true; }
         if (route_store.button_1_released (*this, event)) { return true; }
         if (shape_store.button_1_released (*this, event)) { return true; }
         break;
 
      }

      case 2:
      {
         if (route_store.button_2_released (*this, event)) { return true; }
         break;
      }

      case 3:
      {
         if (marker_store.button_3_released (*this, event)) { return true; }
         break;
      }

   }

   Marker::Popup_Menu& marker_popup_menu = get_marker_popup_menu ();
   Route::Popup_Menu& route_popup_menu = get_route_popup_menu ();
   Shape::Popup_Menu& shape_popup_menu = get_shape_popup_menu ();

   if (marker_popup_menu.is_on ())
   {
      marker_popup_menu.hide ();
      marker_popup_menu.reset_index ();
      queue_draw ();
   }

   if (route_popup_menu.is_on ())
   {
      route_popup_menu.hide (); 
      route_popup_menu.reset_index ();
      queue_draw ();
   }

   if (shape_popup_menu.is_on ())
   {
      shape_popup_menu.hide ();
      shape_popup_menu.reset_index ();
      queue_draw ();
   }

   return false;

}

bool
Console_2D::on_mouse_motion (const Dmouse_Motion_Event& event)
{

   //if (Gtk::Main::events_pending ()) { return true; }
   if (Dcontainer::on_mouse_motion (event)) { return true; }

   const Point_2D& point = event.point;
   const Transform_2D& transform = get_transform ();
   Console_2D::Zoom_Box& zoom_box = get_zoom_box ();

   // moving zoom box around
   if (zoom_box.is_visible () && zoom_box.manipulate_data.is_manipulating ())
   {
      zoom_box.manipulate (event);
      return true;
   }

   // shift drag: translating
   if (manipulate_data.is_manipulating ())
   {
      manipulate (event);
      return true;
   }

   Marker_Store& marker_store = get_marker_store ();
   Route_Store& route_store = get_route_store ();
   Shape_Store& shape_store = get_shape_store ();

   // Moving a marker, route or shape
   if (marker_store.try_dragging (*this, event)) { return true; }
   if (route_store.try_dragging (*this, event)) { return true; }
   if (shape_store.try_dragging (*this, event)) { return true; }

   return false;

}

bool
Console_2D::on_mouse_scroll (const Dmouse_Scroll_Event& event)
{

   if (Dcontainer::on_mouse_scroll (event)) { return true; }

   Console_2D::Zoom_Box& zoom_box = get_zoom_box ();
   if (zoom_box.is_visible ())
   {
      zoom_box.scroll (event);
      queue_draw ();
      return true;
   }

   if (event.control ())
   {
      scroll (event);
      return true;
   }

   return false;

}

Console_2D::Marker::Popup_Menu&
Console_2D::get_marker_popup_menu ()
{
   return marker_popup_menu;
}

Console_2D::Route::Popup_Menu&
Console_2D::get_route_popup_menu ()
{
   return route_popup_menu;
}

Console_2D::Shape::Popup_Menu&
Console_2D::get_shape_popup_menu ()
{
   return shape_popup_menu;
}

Dstring
Console_2D::get_string (const Marker& marker) const
{
   return Dstring::render ("%f %f", marker.x, marker.y);
}

Tokens
Console_2D::get_tokens (const Marker& marker) const
{
   Tokens tokens;
   tokens.push_back (get_string (marker));
   return tokens;
}

Console_2D::Console_2D (Gtk::Window& gtk_window,
                        const Size_2D& size_2d,
                        const bool computer_convention)
   : Dcanvas (gtk_window),
     ignore_heavy (false),
     computer_convention (computer_convention),
     zoom_box (*this),
     marker_popup_menu (*this, 10),
     route_popup_menu (*this, 10),
     shape_popup_menu (*this, 10),
     marker_store (),
     route_store (),
     shape_store ()
{

   Glib::Mutex::Lock lock (mutex);

   Gdk::EventMask event_mask = (Gdk::SCROLL_MASK);
   event_mask |= (Gdk::POINTER_MOTION_MASK);
   event_mask |= (Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK);
   event_mask |= (Gdk::KEY_PRESS_MASK | Gdk::KEY_RELEASE_MASK);
   set_events (event_mask);

   // This basically sets the minimum size
   set_size_request (100, 100);
   set_can_focus ();

   set_hexpand (true);
   set_vexpand (true);

   const Real width = size_2d.i;
   const Real height = size_2d.j;
   being_packed (Point_2D (0, 0), width, height);

   marker_popup_menu.set_hidable (false);
   route_popup_menu.set_hidable (false);
   shape_popup_menu.set_hidable (false);

   register_widget (marker_popup_menu);
   register_widget (route_popup_menu);
   register_widget (shape_popup_menu);

   reset_transform ();

   Console_2D::Zoom_Box& zoom_box = get_zoom_box ();
   zoom_box.reset ();

}

Console_2D::~Console_2D ()
{
}

const Affine_Transform_2D&
Console_2D::get_affine_transform () const
{
   return affine_transform;
}

Affine_Transform_2D&
Console_2D::get_affine_transform ()
{
   return affine_transform;
}

const Transform_2D&
Console_2D::get_transform () const
{
   return affine_transform;
}

Transform_2D&
Console_2D::get_transform ()
{
   return affine_transform;
}

void
Console_2D::set_hide_hidable (const bool hide_hidable)
{
   Dcontainer::set_hide_hidable (hide_hidable);
   queue_draw ();
}

void
Console_2D::clear_marker (const Integer marker_id)
{
   Marker_Store& marker_store = get_marker_store ();
   marker_store.erase (marker_id);
   refresh_all ();
}

void
Console_2D::clear_route (const Integer route_id)
{
   Route_Store& route_store = get_route_store ();
   route_store.erase (route_id);
   route_store.reset ();
   refresh_all ();
}

void
Console_2D::clear_shape (const Integer shape_id)
{
   Shape_Store& shape_store = get_shape_store ();
   shape_store.erase (shape_id);
   shape_store.reset ();
   refresh_all ();
}

void
Console_2D::refresh_all ()
{
   queue_draw ();
}

void
Console_2D::total_refresh ()
{
   set_background_ready (false);
   set_foreground_ready (false);
   render_queue_draw ();
}

void
Console_2D::translate (const Real dx,
                       const Real dy)
{
   affine_transform.translate (dx, dy);
}

void
Console_2D::scroll (const Dmouse_Scroll_Event& event)
{

   const bool up = (event.direction == GDK_SCROLL_UP);
   const Console_2D::Zoom_Box& zoom_box = get_zoom_box ();
   const Real f = zoom_box.scale_ratio.coarse;
   const Real zf = (up ? f : 1.0 / f);

   const Size_2D& size_2d = get_size_2d ();
   const Real w2 = Real (size_2d.i) / 2;
   const Real h2 = Real (size_2d.j) / 2;

   const Point_2D& point = event.point;
   const Point_2D& nw = affine_transform.reverse (point.x, point.y);
   //const Point_2D& nw = affine_transform.reverse (event.point.x - w2,
   //                                               event.point.y - h2);

   affine_transform.translate (-point.x, -point.y);
   affine_transform.scale ((event.alt () ? 1.0 / zf : zf), zf);
   affine_transform.translate (point.x, point.y);

   set_background_ready (false);
   set_foreground_ready (false);
   render_queue_draw ();

}

void
Console_2D::apply_zoom_box ()
{

   Console_2D::Zoom_Box& zoom_box = get_zoom_box ();
   Affine_Transform_2D& affine_transform = get_affine_transform ();

   zoom_box.apply_to (affine_transform);

   zoom_box.set_visible (false);
   zoom_box.reset ();
   zoom_box.post_manipulate ();

   set_background_ready (false);
   set_foreground_ready (false);
   render_queue_draw ();

}

const Console_2D::Marker_Store&
Console_2D::get_marker_store () const
{
   return marker_store;
}

Console_2D::Marker_Store&
Console_2D::get_marker_store ()
{
   return marker_store;
}

const Console_2D::Route_Store&
Console_2D::get_route_store () const
{
   return route_store;
}

Console_2D::Route_Store&
Console_2D::get_route_store ()
{
   return route_store;
}

const Console_2D::Shape_Store&
Console_2D::get_shape_store () const
{
   return shape_store;
}

Console_2D::Shape_Store&
Console_2D::get_shape_store ()
{
   return shape_store;
}

void
Console_2D::render_huds_and_widgets (const RefPtr<Context>& cr)
{

   Shape_Store& shape_store = get_shape_store ();
   shape_store.cairo (cr, *this);

   Route_Store& route_store = get_route_store ();
   route_store.cairo (cr, *this);

   Marker_Store& marker_store = get_marker_store ();
   marker_store.cairo (cr, *this);

   Console_2D::Zoom_Box& zoom_box = get_zoom_box ();
   zoom_box.cairo (cr);

   Dcontainer::cairo (cr);

}

void
Console_2D::cairo (const RefPtr<Context>& cr)
{
   Glib::Mutex::Lock lock (mutex);
   blit_background_buffer (cr);
   blit_image_buffer (cr);
   blit_foreground_buffer (cr);
   render_huds_and_widgets (cr);
}

void
Console_2D::save_image (const Dstring& file_path)
{

   const Tokens tokens (file_path, ".");
   const Dstring& extension = tokens.back ().get_lower_case ();

   if (extension == "png") { save_png (file_path); return; }
   if (extension == "svg") { save_svg (file_path); return; }

   cerr << "Unrecognized file format" << endl;

}

void
Console_2D::save_png (const Dstring& file_path)
{
   Dcanvas::save_image (file_path);
}

void
Console_2D::save_svg (const Dstring& file_path)
{

   const Size_2D& size_2d = get_size_2d ();
   const string fp (file_path.begin (), file_path.end ());
   Cairo::RefPtr<Cairo::SvgSurface> surface =
       Cairo::SvgSurface::create (fp, size_2d.i, size_2d.j);
   Cairo::RefPtr<Cairo::Context> cr = denise::get_cr (surface);

   render_background_buffer (cr);
   render_image_buffer (cr);
   render_foreground_buffer (cr);
   render_huds_and_widgets (cr);
   cr->show_page ();

}

Map_Console::Zoom_Box::Zoom_Box (Map_Console& map_console)
   : Console_2D::Zoom_Box (map_console),
     map_console (map_console)
{
}

Map_Console::Zoom_Box::~Zoom_Box ()
{
}

void
Map_Console::Zoom_Box::apply_to (Geodetic_Transform& transform) const
{
}

void
Map_Console::Zoom_Box::pre_manipulate (const Dmouse_Button_Event& event)
{
   manipulate_data.set_point (event.point);
}

void
Map_Console::Zoom_Box::manipulate (const Dmouse_Motion_Event& event)
{

   if (manipulate_data.is_translating ())
   {

      const Point_2D& point = event.point;
      const Geodetic_Transform& gt = map_console.get_geodetic_transform ();

      centre.x += point.x - manipulate_data.point.x;
      centre.y += point.y - manipulate_data.point.y;
      manipulate_data.set_point (point);

      Real latitude, longitude;
      gt.reverse (latitude, longitude, centre.x, centre.y);
      tilt = M_PI/2 + gt.get_theta (0, 1, latitude, longitude);

   }

   map_console.queue_draw ();

}

void
Map_Console::Zoom_Box::post_manipulate ()
{
   manipulate_data.reset ();
   map_console.queue_draw ();
}

void
Map_Console::Zoom_Box::scroll (const Dmouse_Scroll_Event& event)
{

   const Real f = (event.control () ? scale_ratio.fine : scale_ratio.coarse);
   const bool up = (event.direction == GDK_SCROLL_UP);
   const Real zf = (up ? 1.0 / f : f);

   this->zoom_factor_x *= zf;
   this->zoom_factor_y *= zf;

}

Domain_2D
Map_Console::get_render_domain (const Geodetic_Vector_Data_2D& gvd_2d) const
{

   const Domain_2D& model_domain = gvd_2d.get_domain_2d ();

   const Geodetic_Transform& transform = get_geodetic_transform ();
   const Geodetic_Transform::Genre gtg = transform.data.genre;
   if (gtg == Geodetic_Transform::POLAR_STEREOGRAPHIC_NORTH ||
       gtg == Geodetic_Transform::POLAR_STEREOGRAPHIC_SOUTH)
   {
      return model_domain;
   }

   const Real& model_start_latitude = model_domain.domain_x.start;
   const Real& model_end_latitude = model_domain.domain_x.end;
   const Real& model_start_longitude = model_domain.domain_y.start;
   const Real& model_end_longitude = model_domain.domain_y.end;

   Domain_2D domain_2d = get_domain_2d ();
   Real& start_latitude = domain_2d.domain_x.start;
   Real& end_latitude   = domain_2d.domain_x.end;
   Real& start_longitude = domain_2d.domain_y.start;
   Real& end_longitude   = domain_2d.domain_y.end;

   start_latitude = max (start_latitude, model_start_latitude);
   end_latitude   = min (end_latitude, model_end_latitude);
   start_longitude = max (start_longitude, model_start_longitude);
   end_longitude   = min (end_longitude, model_end_longitude);

   return domain_2d;

}

const Map_Console::Zoom_Box&
Map_Console::get_zoom_box () const
{
   return geodetic_zoom_box;
}

Map_Console::Zoom_Box&
Map_Console::get_zoom_box ()
{
   return geodetic_zoom_box;
}

bool
Map_Console::on_key_pressed (const Dkey_Event& event)
{

   switch (event.value)
   {

      case GDK_KEY_Escape:
      {
         unify_drawers (event);
         return true;
         break;
      }

   }

   return Console_2D::on_key_pressed (event);

}

void
Map_Console::render_overlays (const RefPtr<Context>& cr)
{
   cr->save ();
   const Geodetic_Transform& transform = get_geodetic_transform ();
   const Map_Console::Overlay_Store& os = get_overlay_store ();
   os.cairo (cr, transform, ignore_heavy);
   cr->restore ();
}

void
Map_Console::render_mesh (const RefPtr<Context>& cr)
{

   const Size_2D size_2d (100, 100);
   const Domain_2D domain_2d = get_domain_2d ();
   const Geodetic_Transform& transform = get_geodetic_transform ();

   //const Simple_Mesh_2D sm0_2 (Color::black(0.05), 0.2, 0.2);
   const Geodetic_Mesh mesh_small (1.0, 1.0, Color::black (0.10),
      10., 10., Color::black (0.40), size_2d, domain_2d);
   const Geodetic_Mesh mesh_large (5.0, 5.0, Color::black (0.10),
      30., 30., Color::black (0.40), size_2d, domain_2d);

   const Real latitude_span = domain_2d.domain_x.get_span ();
   const Real longitude_span = domain_2d.domain_y.get_span ();
   const Real span = std::min (latitude_span, longitude_span);
   const Geodetic_Mesh& mesh = (span > 90 ? mesh_large : mesh_small);

   cr->save ();
   mesh.cairo (cr, transform);
   cr->restore ();

}

const Domain_2D
Map_Console::get_domain_2d () const
{
   const Geodetic_Transform& transform = get_geodetic_transform ();
   const Size_2D& size_2d = get_size_2d ();
   return transform.get_domain_2d (size_2d);
}

void
Map_Console::Overlay::init (Geodetic_Cairoable* gc_ptr,
                            const bool filled,
                            const Real line_width,
                            const Color& color,
                            const bool heavy)
{
   this->gc_ptr = gc_ptr;
   this->filled = filled;
   this->line_width = line_width;
   this->color = color;
   this->heavy = heavy;
}

Map_Console::Overlay::Overlay (Geodetic_Cairoable* gc_ptr,
                               const bool filled,
                               const Real line_width,
                               const Color& color,
                               const bool heavy)
{
   init (gc_ptr, filled, line_width, color, heavy);
}

const void
Map_Console::Overlay::cairo (const RefPtr<Context>& cr,
                             const Geodetic_Transform& transform) const
{
   cr->set_line_width (line_width);
   color.cairo (cr);
   gc_ptr->cairo (cr, transform);
   if (filled) { cr->fill (); } else { cr->stroke (); }
}

Map_Console::Overlay_Store::~Overlay_Store ()
{
   for (Overlay_Ptr_Map::iterator iterator = overlay_ptr_map.begin ();
        iterator != overlay_ptr_map.end (); iterator++)
   {
      Overlay* overlay_ptr = iterator->second;
      delete overlay_ptr;
   }
}

const Map_Console::Overlay&
Map_Console::Overlay_Store::get_overlay (const Dstring& identifier) const
{
   Overlay_Ptr_Map::const_iterator iterator = overlay_ptr_map.find (identifier);
   const Overlay* overlay_ptr = iterator->second;
   return *overlay_ptr;
}

void
Map_Console::Overlay_Store::add (const Dstring& overlay_string,
                                 const bool heavy)
{

   const Tokens tokens (overlay_string, ":");
   if (tokens.size () < 7) { return; }

   const Dstring& genre = tokens[0];
   const Dstring& identifier = tokens[1];
   const Real& line_width = stof (tokens[2]);
   const Real& r = stof (tokens[3]);
   const Real& g = stof (tokens[4]);
   const Real& b = stof (tokens[5]);
   const Real& a = stof (tokens[6]);
   const Color color (r, g, b, a);

   if (genre == "GSHHS")
   {
      const bool filled = (tokens[7] == "filled");
      const Dstring& file_path = tokens[8];
      Gshhs* gshhs_ptr = new Gshhs (file_path);
      add (identifier, gshhs_ptr, filled, line_width, color, heavy);
      return;
   }

   if (genre == "OUT")
   {
      const Dstring& file_path = tokens[7];
      Outl* outl_ptr = new Outl (file_path);
      add (identifier, outl_ptr, false, line_width, color, heavy);
      return;
   }

   if (genre == "KIDNEY")
   {
      const Dstring& file_path = tokens[7];
      Kidney* kidney_ptr = new Kidney (file_path);
      add (identifier, kidney_ptr, false, line_width, color, heavy);
      return;
   }

}

void
Map_Console::Overlay_Store::add (const Dstring& identifier,
                                 Geodetic_Cairoable* gc_ptr,
                                 const bool filled, 
                                 const Real line_width,
                                 const Color& color,
                                 const bool heavy)
{
   push_back (identifier);
   on_off_map.insert (make_pair (identifier, !heavy));
   Overlay* o_ptr = new Overlay (gc_ptr, filled, line_width, color, heavy);
   overlay_ptr_map.insert (make_pair (identifier, o_ptr));
}

bool
Map_Console::Overlay_Store::is_on (const Dstring& identifier) const
{
   On_Off_Map::const_iterator iterator = on_off_map.find (identifier);
   if (iterator == on_off_map.end ())
   {
      throw Exception ("Invalid overlay identifier");
   }
   return iterator->second;
}

void
Map_Console::Overlay_Store::set_on_off (const Dstring& identifier,
                                        const bool on_off)
{
   on_off_map[identifier] = on_off;
}

void
Map_Console::Overlay_Store::toggle_on_off (const Dstring& identifier)
{
   const bool b = on_off_map[identifier];
   on_off_map[identifier] = !b;
}

void
Map_Console::Overlay_Store::cairo (const RefPtr<Context>& cr,
                                   const Geodetic_Transform& transform,
                                   const bool ignore_heavy) const
{

   cr->save ();
   cr->set_line_join (LINE_JOIN_ROUND);

   for (Overlay_Store::const_iterator iterator = begin ();
        iterator != end (); iterator++)
   {

      const Dstring& identifier = *(iterator);
      On_Off_Map::const_iterator i = on_off_map.find (identifier);
      if (i == on_off_map.end ()) { cout << "cont a " << endl; continue; }

      const bool on = i->second;
      if (!on) { cout << "cont b " << endl; continue; }

      const Overlay& overlay = get_overlay (identifier);
      if (ignore_heavy && overlay.heavy) { cout << "cont c " << endl; continue; }
      overlay.cairo (cr, transform);

   }

   cr->restore ();

}

Map_Console::Option_Panel::Zoom_Drawer::Zoom_Drawer (Option_Panel& option_panel)
   : Drawer (option_panel.map_console, option_panel, false, "Zoom", 12) 
{
}

void
Map_Console::Option_Panel::Zoom_Drawer::switch_off_zoom_button ()
{

   for (auto iterator = widget_ptr_map.begin ();
        iterator != widget_ptr_map.begin (); iterator++)
   {

      Dwidget& widget = *(iterator->second);

      try
      {
         Dtoggle_Button& tb = dynamic_cast<Dtoggle_Button&>(widget);
         if (tb.get_str () != "Custom") { continue; }
         tb.set (false);
         break;
      }
      catch (const std::bad_cast& bc)
      {
      }

   }

}

Map_Console::Option_Panel::Overlay_Drawer::Overlay_Drawer (Option_Panel& option_panel)
   : Drawer (option_panel.map_console, option_panel, false, "Overlay", 12) 
{
}

void
Map_Console::Option_Panel::Overlay_Drawer::add_toggle_button_ptr (const Dstring& str,
                                                                  Dtoggle_Button* toggle_button_ptr)
{
   add_widget_ptr (toggle_button_ptr);
   index_map.insert (make_pair (str, index_map.size ()));
}

bool
Map_Console::Option_Panel::Overlay_Drawer::is_on (const Dstring& str) const
{

   for (auto iterator = widget_ptr_map.begin ();
        iterator != widget_ptr_map.end (); iterator++)
   {

      const Dwidget& widget = *(iterator->second);

      try
      {
         const Dtoggle_Button& tb = dynamic_cast<const Dtoggle_Button&>(widget);
         if (tb.get_str () != str) { continue; }
         return tb.is_switched_on ();
      }
      catch (const std::bad_cast& bc)
      {
      }

   }

   return false;

}

void
Map_Console::Option_Panel::Overlay_Drawer::set_on_off (const Dstring& str,
                                                       const bool on_off)
{
   
   std::map<Dstring, Integer>::iterator iterator = index_map.find (str);
   if (iterator == index_map.end ()) { return ; }
   
   const Integer index = iterator->second;
   Dtoggle_Button& tb = (Dtoggle_Button&)(*(widget_ptr_map.at (index)));
   tb.set (on_off);

}

void
Map_Console::Option_Panel::setup_zoom (const Tokens& config_file_content)
{

   const Real font_size = 12;
   typedef Geodetic_Transform::Data Gtd;
   typedef Template_Button<Gtd> Gtdb;
   Map_Console& mc = map_console;

   // Make Drawer
   Zoom_Drawer* drawer_ptr = new Zoom_Drawer (*this);
   add_drawer_ptr (drawer_ptr, true);

   // Make Custom Zoom Button
   Dbutton* button_ptr = new Dbutton (map_console, "Custom", font_size);
   Dbutton::Signal& signal = button_ptr->get_signal ();
   Console_2D::Zoom_Box& zoom_box = map_console.get_zoom_box ();
   signal.connect (sigc::mem_fun (zoom_box, &Console_2D::Zoom_Box::toggle_visible));
   drawer_ptr->add_widget_ptr (button_ptr);

   // Make Zoom buttons for default Zooms
   for (auto iterator = config_file_content.begin ();
        iterator != config_file_content.end (); iterator++)
   {

      const Tokens tokens (*(iterator));
      if (tokens.size () != 3) { continue; }
      if (tokens[0] != "geodetic_transform") { continue; }

      const Dstring& identifier = tokens[1];
      const Geodetic_Transform::Data gtd (tokens[2]);

      Gtdb* gtdb_ptr = new Gtdb (map_console, gtd, identifier, font_size);
      Gtdb::T_Signal& signal = gtdb_ptr->get_t_signal ();
      signal.connect (sigc::mem_fun (
         mc, &Map_Console::set_geodetic_transform_data));

      drawer_ptr->add_widget_ptr (gtdb_ptr);

   }

}

Map_Console::Option_Panel::Option_Panel (Map_Console& map_console,
                                         const Tokens& config_file_content)
   : Drawer_Panel (map_console, false, 12),
     map_console (map_console)
{
   Overlay_Drawer* drawer_ptr = new Overlay_Drawer (*this);
   add_drawer_ptr (drawer_ptr, true);
   setup_zoom (config_file_content);
}

void
Map_Console::Option_Panel::setup_overlay ()
{

   const Overlay_Store& overlay_store = map_console.get_overlay_store ();
   Overlay_Drawer& overlay_drawer = get_overlay_drawer ();

   // Make Drawer
   const Real font_size = 12;
   Map_Console& mc = map_console;
   typedef Dtoggle_Button Tb;

   // Make Overlay Buttons
   for (Overlay_Store::const_iterator iterator = overlay_store.begin ();
        iterator != overlay_store.end (); iterator++)
   {

      const Dstring& identifier = *(iterator);
      const bool on = overlay_store.is_on (identifier);

      Tb* tb_ptr = new Tb (map_console, identifier, font_size, on);
      Tb::Str_Signal& signal = tb_ptr->get_str_signal ();
      signal.connect (sigc::mem_fun (mc, &Map_Console::toggle_overlay));

      overlay_drawer.add_toggle_button_ptr (identifier, tb_ptr);

   }

}

void
Map_Console::Option_Panel::switch_off_zoom_button ()
{
   Zoom_Drawer& zoom_drawer = get_zoom_drawer ();
   zoom_drawer.switch_off_zoom_button ();
}

void
Map_Console::Option_Panel::align_overlay_toggles ()
{

   const Overlay_Store& overlay_store = map_console.get_overlay_store ();
   Overlay_Drawer& overlay_drawer = get_overlay_drawer ();

   // Make Drawer
   Map_Console& mc = map_console;
   typedef Dtoggle_Button Tb;

   // Make Overlay Buttons
   for (Overlay_Store::const_iterator iterator = overlay_store.begin ();
        iterator != overlay_store.end (); iterator++)
   {
      const Dstring& identifier = *(iterator);
      const bool on_off = overlay_store.is_on (identifier);
      overlay_drawer.set_on_off (identifier, on_off);
   }

}


bool
Map_Console::Option_Panel::overlay_is_on (const Dstring& overlay_str) const
{
   const Overlay_Drawer& overlay_drawer = get_overlay_drawer ();
   return overlay_drawer.is_on (overlay_str);
}

const Map_Console::Option_Panel::Zoom_Drawer&
Map_Console::Option_Panel::get_zoom_drawer () const
{
   const Drawer& drawer = *(drawer_ptr_map.find ("Zoom")->second);
   return dynamic_cast<const Zoom_Drawer&>(drawer);
}

Map_Console::Option_Panel::Zoom_Drawer&
Map_Console::Option_Panel::get_zoom_drawer ()
{
   Drawer& drawer = *(drawer_ptr_map.find ("Zoom")->second);
   return dynamic_cast<Zoom_Drawer&>(drawer);
}

const Map_Console::Option_Panel::Overlay_Drawer&
Map_Console::Option_Panel::get_overlay_drawer () const
{
   const Drawer& drawer = *(drawer_ptr_map.find ("Overlay")->second);
   return dynamic_cast<const Overlay_Drawer&>(drawer);
}

Map_Console::Option_Panel::Overlay_Drawer&
Map_Console::Option_Panel::get_overlay_drawer ()
{
   Drawer& drawer = *(drawer_ptr_map.find ("Overlay")->second);
   return dynamic_cast<Overlay_Drawer&>(drawer);
}

Console_2D::Route*
Map_Console::Route_Store::new_route_ptr (const Integer id,
                                         const Point_2D& point,
                                         const Real node_size)
{
   return new Map_Console::Route (id, point, node_size);
}

Console_2D::Route*
Map_Console::Route_Store::new_route_ptr (const Integer id,
                                         const Point_2D& point_a,
                                         const Point_2D& point_b,
                                         const Real node_size)
{
   return new Map_Console::Route (id, point_a, point_b, node_size);
}

Map_Console::Route::Route (const Integer id,
                           const Point_2D& point,
                           const Real node_size)
   : Console_2D::Route (id, point, node_size)
{
}

Map_Console::Route::Route (const Integer id,
                           const Point_2D& point_a,
                           const Point_2D& point_b,
                           const Real node_size)
   : Console_2D::Route (id, point_a, point_b, node_size)
{
}

bool
Map_Console::Route::matches (const Transform_2D& transform,
                             const Point_2D& point) const
{
   const Geodesy geodesy;
   const Journey j (*this);
   return (j.get_iterator (transform, point, geodesy, node_size) != j.end ());
}

void
Map_Console::Route::cairo (const RefPtr<Context>& cr,
                           const Console_2D& console_2d) const
{

   const Point_2D& first_p = *(begin ());
   if (first_p.is_nap ()) { return; }

   const Transform_2D& transform = console_2d.get_transform ();

   const Journey journey (*this);
   journey.cairo (cr, transform);

}

Console_2D::Route::iterator
Map_Console::Route::add (const Transform_2D& transform,
                         const Point_2D& point)
{
   const Geodesy geodesy;
   Journey journey (*this);
   Console_2D::Route::iterator i = journey.implant (
      transform, point, geodesy, node_size);
   set (journey);

   for (Console_2D::Route::iterator iterator = begin ();
        iterator != end (); iterator++)
   {
      const Point_2D& p = transform.transform (*(iterator));
      const Real dx = p.x - point.x;
      const Real dy = p.y - point.y;
      if (dx*dx + dy*dy < node_size*node_size/4) { return iterator; }
   }

   return end ();

}

Real
Map_Console::Route::get_distance (const Geodesy& geodesy) const
{
   const Journey journey (*this);
   return journey.get_distance (geodesy);
}

Dstring
Map_Console::get_string (const Marker& marker) const
{
   const Lat_Long lat_long (marker);
   return lat_long.get_string (false, Dstring ("%.3f\u00b0"));
}

Map_Console::Map_Console (Gtk::Window& gtk_window,
                          const Size_2D& size_2d,
                          const Tokens& config_file_content)
   : Console_2D (gtk_window, size_2d),
     geodetic_zoom_box (*this),
     option_panel (*this, config_file_content)
{

   Glib::Mutex::Lock lock (mutex);

   typedef Geodetic_Transform Gt;
   typedef Geodetic_Transform::Data Gtd;
   typedef Template_Button<Gtd> Gtdb;
   const Point_2D centre (size_2d.i/2, size_2d.j/2);

//   auto zoom_drawer = option_panel.get_zoom_drawer ();
//   const Gtdb& zoom_widget = (const Gtdb&)(zoom_drawer.get_widget (0));
//   const Gt::Data& gtd = zoom_widget.get_t ();
//   geodetic_transform_ptr = Gt::get_transform_ptr (gtd, centre);

   for (auto iterator = config_file_content.begin ();
        iterator != config_file_content.end (); iterator++)
   {

      const Tokens tokens (*(iterator));
      if (tokens.size () != 3) { continue; }
      if (tokens[0] != "geodetic_transform") { continue; }

      const Dstring& identifier = tokens[1];
      const Geodetic_Transform::Data gtd (tokens[2]);
      geodetic_transform_ptr = Gt::get_transform_ptr (gtd, centre);

      break;

   }



   Map_Console::Zoom_Box& zoom_box = get_zoom_box ();
   zoom_box.reset ();

}

Map_Console::~Map_Console ()
{
   delete geodetic_transform_ptr;
}

Map_Console::Transform_Signal&
Map_Console::get_transform_signal ()
{
   return transform_signal;
}

const Map_Console::Overlay_Store&
Map_Console::get_overlay_store () const
{
   return overlay_store;
}

Map_Console::Overlay_Store&
Map_Console::get_overlay_store ()
{
   return overlay_store;
}

const Geodetic_Transform&
Map_Console::get_geodetic_transform () const
{
   return *geodetic_transform_ptr;
}

Geodetic_Transform&
Map_Console::get_geodetic_transform ()
{
   return *geodetic_transform_ptr;
}

const Transform_2D&
Map_Console::get_transform () const
{
   return *geodetic_transform_ptr;
}

Transform_2D&
Map_Console::get_transform ()
{
   return *geodetic_transform_ptr;
}

void
Map_Console::set_geodetic_transform_data (const Geodetic_Transform::Data& gtd)
{

   typedef Geodetic_Transform Gt;
   const Size_2D& size_2d = get_size_2d ();
   const Point_2D middle (size_2d.i / 2, size_2d.j / 2);

   delete geodetic_transform_ptr;
   geodetic_transform_ptr =
      Gt::get_transform_ptr (gtd.genre, gtd.scale, gtd.lat_long, middle);
   cout << geodetic_transform_ptr->data.get_string () << endl;

   zoom_box.reset ();
   set_background_ready (false);
   set_foreground_ready (false);
   render_queue_draw ();

}

void
Map_Console::apply_zoom_box ()
{

   Console_2D::Zoom_Box& zoom_box = get_zoom_box ();

   const Geodetic_Transform& transform = get_geodetic_transform ();
   const Geodetic_Transform::Genre genre = transform.data.genre;
   const Real scale = transform.get_scale () / zoom_box.zoom_factor_x;
   const Lat_Long lat_long = transform.get_lat_long (zoom_box.centre);

   zoom_box.set_visible (false);
   zoom_box.reset ();
   zoom_box.post_manipulate ();

   const Geodetic_Transform::Data gtd (genre, scale, lat_long);
   set_geodetic_transform_data (gtd);

}

void
Map_Console::refresh_all ()
{
   queue_draw ();
}

void
Map_Console::toggle_overlay (const Dstring& overlay_identifier)
{
   Map_Console::Overlay_Store& os = get_overlay_store ();
   os.toggle_on_off (overlay_identifier);
   total_refresh ();
}

void
Map_Console::translate (const Real dx,
                        const Real dy)
{

   typedef Geodetic_Transform Gt;
   const Point_2D centre (width / 2, height / 2);

   const Geodetic_Transform& transform = get_geodetic_transform ();
   const Geodetic_Transform::Genre genre = transform.data.genre;
   const Real scale = transform.get_scale ();
   Lat_Long lat_long = transform.get_lat_long ();

   const Geodetic_Transform::Data& data = transform.data;
   Gt* t_ptr = Gt::get_transform_ptr (data, centre + Point_2D (dx, dy));
   t_ptr->reverse (lat_long, centre);
   delete t_ptr;

   const Geodetic_Transform::Data gtd (genre, scale, lat_long);
   set_geodetic_transform_data (gtd);

}

void
Map_Console::scroll (const Dmouse_Scroll_Event& event)
{

   const Console_2D::Zoom_Box& zoom_box = get_zoom_box ();
   const bool up = (event.direction == GDK_SCROLL_UP);
   const Real f = zoom_box.scale_ratio.coarse;
   const Real zf = (up ? 1.0/f : f);

   typedef Geodetic_Transform Gt;
   const Point_2D centre (width / 2, height / 2);

   const Geodetic_Transform& transform = get_geodetic_transform ();
   const Geodetic_Transform::Genre genre = transform.data.genre;
   const Real scale = transform.get_scale () * zf;
   Lat_Long lat_long = transform.get_lat_long ();

   const Geodetic_Transform::Data& data = transform.data;
   Gt* t_ptr = Gt::get_transform_ptr (data, centre);
   t_ptr->reverse (lat_long, centre);
   delete t_ptr;

   const Geodetic_Transform::Data gtd (genre, scale, lat_long);
   set_geodetic_transform_data (gtd);

}

void
Map_Console::align_overlay_toggles ()
{
   option_panel.align_overlay_toggles ();
}

void
Map_Console::unify_drawers (const Devent& event)
{

   const bool o = option_panel.all_collapsed ();
   const bool expand = (o);

   if (expand) { option_panel.collapse_all (); }
   else { option_panel.expand_all (); }

   queue_draw ();

}

const Console_2D::Route_Store&
Map_Console::get_route_store () const
{
   return mc_route_store;
}

Console_2D::Route_Store&
Map_Console::get_route_store ()
{
   return mc_route_store;
}

/*
Time_Series_Canvas::Transform::Transform (const Size_2D& size_2d,
                                          const Real margin_t,
                                          const Real margin_b,
                                          const Real margin_l,
                                          const Real margin_r)
   : margin_t (margin_t),
     margin_b (margin_b),
     margin_l (margin_l),
     margin_r (margin_r)
{
   update (size_2d);
}

void
Time_Series_Canvas::Transform::update (const Size_2D& size_2d)
{

   const Real w = size_2d.i - margin_l - margin_r;
   const Real h = size_2d.j - margin_t - margin_b;

   const Domain_1D& domain_x = domain_2d.domain_x;
   const Domain_1D& domain_y = domain_2d.domain_y;

   set_identity ();
   translate (-domain_x.start, -domain_y.start);
   scale (w / domain_x.get_span (), h / domain_y.get_span ());
   translate (margin_l, margin_t);

}

void
Time_Series_Canvas::Transform::set_domain_2d (const Domain_2D& domain_2d)
{
   this->domain_2d = domain_2d;
}

const Domain_1D&
Time_Series_Canvas::Transform::get_domain_t () const
{
   return domain_2d.domain_x;
}

Domain_1D&
Time_Series_Canvas::Transform::get_domain_t ()
{
   return domain_2d.domain_x;
}

const Domain_1D&
Time_Series_Canvas::Transform::get_domain_y () const
{
   return domain_2d.domain_y;
}

Domain_1D&
Time_Series_Canvas::Transform::get_domain_y ()
{
   return domain_2d.domain_y;
}

bool
Time_Series_Canvas::on_configure_event (GdkEventConfigure* event)
{

   const Size_2D size_2d (event->width, event->height);
   if (get_size_2d () == size_2d) { return false; }

   this->width = size_2d.i;
   this->height = size_2d.j;

   transform.update (size_2d);

   background_buffer.clear ();
   image_buffer.clear ();
   foreground_buffer.clear ();
   initialize ();
   render_queue_draw ();

   return true;

}

void
Time_Series_Canvas::render_background_buffer (const RefPtr<Context>& cr)
{

   cr->set_line_width (2);
   Color (1, 1, 1).cairo (cr);
   cr->paint ();

   const Dstring fmt_t ("%HZ");

   mesh_2d.render (cr, transform);

   Color (0, 0, 0).cairo (cr);

   const Domain_2D& domain_2d = mesh_2d.get_domain_2d ();
   const Real x = domain_2d.domain_x.start;
   const Real y = domain_2d.domain_y.end;

   mesh_2d.render_label_x (cr, transform, 1, y,
      fmt_t, NUMBER_TIME, 'l', 't', 5);
   mesh_2d.render_label_y (cr, transform, 1, x,
      fmt_y, NUMBER_REAL, 'r', 'c', 5);


   const Real margin = 6;
   const Integer hours = 24;
   const Dtime start_time (domain_2d.domain_x.start);
   const Dtime end_time (domain_2d.domain_x.end);
   const Dtime st (start_time.t - fmod (start_time.t, hours));

   cr->save ();

   for (Dtime dtime = st; dtime <= end_time; dtime.t += hours)
   {

      if (dtime < start_time) { continue; }

      const Point_2D p (dtime.t, domain_2d.domain_y.end);
      const Point_2D& point = transform.transform (p);
      const Point_2D point_a (point.x, point.y + 15);
      const Point_2D point_b (point.x, point.y + 25);
      const Dstring str_a = dtime.get_string ("%b %d");
      const Dstring str_b = dtime.get_string ("(%a)");

      const Point_2D point_c (point.x, point.y + 5);
      const Point_2D point_d (point.x, point.y + 40);
      cr->set_line_width (1);
      cr->move_to (point_c.x, point_c.y);
      cr->line_to (point_d.x, point_d.y);
      cr->stroke ();

      Label label_a (str_a, point_a, 'l', 't', 5);
      Label label_b (str_b, point_b, 'l', 't', 5);
      label_a.cairo (cr);
      label_b.cairo (cr);

   }

   cr->restore ();

}

void
Time_Series_Canvas::render_time_series (const RefPtr<Context>& cr,
                                        const Vector_Data_1D& time_series,
                                        const Integer vector_element) const
{

   const Real dt = 0.1;
   const Domain_1D& domain_t = transform.get_domain_t ();
   const Integer n = Integer (round (domain_t.get_span () / dt)) + 1;

   const Real node_size = 2;

   Point_2D point;
   Simple_Polyline simple_polyline;

   for (Integer i = 0; i < n; i++)
   {

      const Real t = domain_t.start + i * dt;
      if (time_series.out_of_bounds (0, t)) { continue; }
      const Real datum = time_series.evaluate (vector_element, t);

      transform.transform (point.x, point.y, t, datum);
      simple_polyline.add (point);

   }

   simple_polyline.cairo (cr);
   cr->stroke ();

   for (Integer i = 0; i < time_series.size (); i++)
   {

      const Real t = time_series.get_coordinate (0, i);
      if (time_series.node_out_of_bounds (0, i)) { continue; }
      const Real datum = time_series.get_datum (vector_element, i);

      transform.transform (point.x, point.y, t, datum);
      Ring (node_size).cairo (cr, point);
      cr->fill ();

   }

}

void
Time_Series_Canvas::render_sun_elevation (const RefPtr<Context>& cr,
                                          const Lat_Long& lat_long) const
{

   const Real hours = 0.5;
   const Domain_1D& domain_t = transform.get_domain_t ();
   const Dtime st (domain_t.start - fmod (domain_t.start, hours));

   const Real& latitude = lat_long.latitude;
   const Real& longitude = lat_long.longitude;

   Color (0.8, 0.7, 0.6, 0.2).cairo (cr);

//   const Real start_y = margin_t;
//   const Real end_y = height - margin_b;
//   const Real chart_height = height - margin_t - margin_b;
//   const Real horizon_y = start_y + (1 - 0) * chart_height / 1.34;
//
//   Point_2D p;
//
//   for (Dtime dtime = st; dtime <= domain_t.end; dtime.t += hours)
//   {
//
//      if (dtime < domain_t.start) { continue; }
//
//      const Zenith_Field zenith_field (Sun (), dtime);
//      const Real zenith = zenith_field.evaluate (0, latitude, longitude);
//      const Real elevation = 90 - zenith;
//
//      if (elevation < -20) { continue; }
//
//      const Real cos_z = cos (zenith * DEGREE_TO_RADIAN);
//      transform.transform (p.x, p.y, dtime.t, 500e2);
//      p.y = start_y + (1 - cos_z) * chart_height / 1.34;
//
//      Ring (6).cairo (cr, p);
//      cr->fill ();
//
//   }
//
//   const Real start_x = margin_l;
//   const Real end_x = width - margin_r;
//
//   Polygon polygon;
//   polygon.add (Point_2D (start_x, end_y));
//   polygon.add (Point_2D (start_x, horizon_y));
//   polygon.add (Point_2D (end_x, horizon_y));
//   polygon.add (Point_2D (end_x, end_y));
//
//   Color (0, 0, 0, 0.08).cairo (cr);
//   polygon.cairo (cr);
//   cr->fill ();

}

Time_Series_Canvas::Time_Series_Canvas (Gtk::Window& gtk_window,
                                        const Size_2D& size_2d,
                                        const Real margin_t,
                                        const Real margin_b,
                                        const Real margin_l,
                                        const Real margin_r)
   : Dcanvas (gtk_window),
     transform (size_2d, margin_t, margin_b, margin_l, margin_r)
{
}

Time_Series_Canvas::~Time_Series_Canvas ()
{
}

*/
