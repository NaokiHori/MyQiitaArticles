#if defined(ENABLE_XWINDOW)

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <X11/Xlib.h>
#include "memory.h"
#include "visualiser.h"

// store XWindow-related stuffs
typedef struct {
  Display * display;
  Window window;
  GC gc;
  XSegment * lines;
  XArc * arcs;
} xobjs_t;

// "global" variables which are only accessible inside this source

// schedulers
static double g_rate = 0.;
static double g_next = 0.;

// XWindow objects
static xobjs_t g_xobjs;

// update display
static int execute(
    const double time,
    const p_object_t * p_object
){
  const size_t nitems = p_solver.get_nitems(p_object);
  Display * display = g_xobjs.display;
  Window window     = g_xobjs.window;
  GC gc             = g_xobjs.gc;
  XSegment * lines  = g_xobjs.lines;
  XArc * arcs       = g_xobjs.arcs;
  // configure positions and sizes nicely
  //   so that the whole system is in the display
  XWindowAttributes attr = {0};
  XGetWindowAttributes(display, window, &attr);
  const unsigned short width  = attr.width;
  const unsigned short height = attr.height;
  const unsigned short pend_originx = attr.x + width  / 2;
  const unsigned short pend_originy = attr.y + height / 2;
  // compute amplification factor
  const double sum_length = 1. * (nitems + 1);
  const unsigned short ref_length
    = width < height
    ? width  / 2
    : height / 2;
  const double factor = ref_length / sum_length;
  const unsigned short radius = (unsigned short)(0.25 * factor);
  // initialise objects to be drawn
  double x1 = 0.;
  double y1 = 0.;
  for (size_t n = 0; n < nitems; n++) {
    // lines connecting gravity centers
    const double p = p_solver.get_pos(p_object, n);
    const double x2 = x1 + cos(p);
    const double y2 = y1 + sin(p);
    lines[n].x1 = (unsigned short)(pend_originx + factor * x1);
    lines[n].y1 = (unsigned short)(pend_originy + factor * y1);
    lines[n].x2 = (unsigned short)(pend_originx + factor * x2);
    lines[n].y2 = (unsigned short)(pend_originy + factor * y2);
    // gravity center to left-top corner
    arcs[n].x      = lines[n].x2 - radius;
    arcs[n].y      = lines[n].y2 - radius;
    arcs[n].width  = 2 * radius;
    arcs[n].height = 2 * radius;
    // angle, from 0 to 360 (circle)
    // multiply 64 to convert degree
    arcs[n].angle1 =   0 * 64;
    arcs[n].angle2 = 360 * 64;
    x1 = x2;
    y1 = y2;
  }
  // clear display before drawing
  XClearWindow(display, window);
  // draw lines
  XSetForeground(display, gc, 0xFFFFFF);
  XDrawSegments(display, window, gc, lines, nitems);
  // draw particles
  XSetForeground(display, gc, 0xFF0000);
  XFillArcs(display, window, gc, arcs, nitems);
  // write time
  const char prefix[] = {"time: "};
  const int ndigits = 5;
  const size_t nchars = strlen(prefix) + (size_t)ndigits + 1;
  char * string = memory_calloc(nchars, sizeof(char));
  snprintf(string, nchars, "%s% *.1f", prefix, ndigits, time);
  // position
  const int time_x = 10;
  const int time_y = 10;
  XSetForeground(display, gc, 0xFFFFFF);
  XDrawString(display, window, gc, time_x, time_y, string, nchars - 1);
  memory_free(string);
  // update display
  XFlush(display);
  g_next += g_rate;
  return 0;
}

// constructor
static int init(
    const double rate,
    const p_object_t * p_object
){
  // set up X connection for in-situ visualisation
  Display ** display = &g_xobjs.display;
  Window * window    = &g_xobjs.window;
  GC * gc            = &g_xobjs.gc;
  XSegment ** lines  = &g_xobjs.lines;
  XArc ** arcs       = &g_xobjs.arcs;
  // open display connection
  *display = XOpenDisplay(NULL);
  if (NULL == *display) {
    printf("XOpenDisplay error\n");
    return 1;
  }
  const int screen = DefaultScreen(*display);
  // create window
  *window = XCreateSimpleWindow(
      /* display   */ *display,
      /* window    */ RootWindow(*display, screen),
      /* x         */ 100,
      /* y         */ 100,
      /* width     */ 800,
      /* height    */ 800,
      /* border_w  */ 1,
      /* border    */ BlackPixel(*display, screen),
      /* backgroud */ BlackPixel(*display, screen)
  );
  if (BadAlloc == *window || BadMatch == *window || BadValue == *window || BadWindow == *window) {
    printf("XCreateSimpleWindow error\n");
    return 1;
  }
  if (BadWindow == XMapWindow(*display, *window)) {
    printf("XMapWindow error\n");
    return 1;
  }
  // graphic context
  *gc = XCreateGC(
      *display,
      *window,
      0,
      NULL
  );
  XStoreName(*display, *window, "Pendulum");
  // allocate memory for drawn objects
  const size_t nitems = p_solver.get_nitems(p_object);
  *lines = memory_calloc(nitems, sizeof(XSegment));
  *arcs  = memory_calloc(nitems, sizeof(    XArc));
  // assume current time is 0,
  //   then the next event will happen at "rate"
  g_rate = rate;
  g_next = rate;
  // render initial field
  execute(0., p_object);
  return 0;
}

// destructor
static int finalise(
    void
){
  memory_free(g_xobjs.lines);
  memory_free(g_xobjs.arcs);
  XFreeGC(g_xobjs.display, g_xobjs.gc);
  XDestroyWindow(g_xobjs.display, g_xobjs.window);
  XCloseDisplay(g_xobjs.display);
  return 0;
}

// getter
static double get_next_time(
    void
){
  return g_next;
}

const visualiser_t visualiser = {
  .init          = init,
  .finalise      = finalise,
  .get_next_time = get_next_time,
  .execute       = execute,
};

#else

extern char dummy;

#endif // ENABLE_XWINDOW
