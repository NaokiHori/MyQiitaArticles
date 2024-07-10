#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <math.h>
#include "pendulum.h"
#include "logger.h"

#define FNAMEMAX 128

// "global" variables which are only accessible inside this source

// schedulers
static double g_rate = 0.;
static double g_next = 0.;

// wrapper of fopen
static FILE * my_fopen(
    const char * restrict path,
    const char * restrict mode
){
  errno = 0;
  FILE * stream = fopen(path, mode);
  if (NULL == stream) {
    perror(path);
    return NULL;
  }
  return stream;
}

// wrapper of fclose
static int my_fclose(
    FILE * stream
){
  fclose(stream);
  return 0;
}

static int output_energy(
    const double time,
    const p_object_t * p_object
){
  const size_t nitems = p_solver.get_nitems(p_object);
  // add gauge value
  struct {
    double t;
    double u;
  } energies = {
    .t = 0.,
    .u = 0.5 * nitems * (nitems + 1),
  };
  for (size_t i = 0; i < nitems; i++) {
    const double pi = p_solver.get_pos(p_object, i);
    const double vi = p_solver.get_vel(p_object, i);
    for (size_t j = 0; j < nitems; j++) {
      const double pj = p_solver.get_pos(p_object, j);
      const double vj = p_solver.get_vel(p_object, j);
      const double m = nitems - fmax(i, j);
      energies.t += 0.5 * m * vi * vj * cos(pi - pj);
    }
    energies.u -= (nitems - i) * sin(pi);
  }
  // normalise by the total energy
  energies.t /= nitems * (nitems + 1);
  energies.u /= nitems * (nitems + 1);
  //
  const char filename[] = {"output/energy.dat"};
  static bool is_called = false;
  FILE * fp = my_fopen(filename, is_called ? "a" : "w");
  if (NULL == fp) {
    return 1;
  }
  fprintf(
      fp,
      "% .15e % .15e % .15e % .15e\n",
      time,
      energies.t,
      energies.u,
      energies.t + energies.u
  );
  my_fclose(fp);
  is_called = true;
  return 0;
}

static int output_objects(
    const double time,
    const p_object_t * p_object
){
  // output xy positions
  const size_t nitems = p_solver.get_nitems(p_object);
  const char prefix[] = {"output/mass"};
  const char suffix[] = {".dat"};
  const int ndigits = 3;
  static bool is_called = false;
  double x = 0.;
  double y = 0.;
  for (size_t n = 0; n < nitems; n++) {
    char filename[FNAMEMAX] = {'\0'};
    snprintf(filename, FNAMEMAX, "%s%0*zu%s", prefix, ndigits, n, suffix);
    FILE * fp = my_fopen(filename, is_called ? "a" : "w");
    if (NULL == fp) {
      return 1;
    }
    const double p = p_solver.get_pos(p_object, n);
    x += cos(p);
    y += sin(p);
    fprintf(fp, "% .15e % .15e % .15e\n", time, x, y);
    my_fclose(fp);
  }
  is_called = true;
  return 0;
}

// constructor
static int init(
    const double rate
){
  // schedule logging events
  // assume current time is 0,
  //   then the next event will happen at "rate"
  g_rate = rate;
  g_next = rate;
  return 0;
}

// destructor
static int finalise(
    void
){
  return 0;
}

// getter
static double get_next_time(
    void
){
  return g_next;
}

static int execute(
    const double time,
    const double dt,
    const p_object_t * p_object
){
  printf("time: % .2e dt: % .2e\n", time, dt);
  g_next += g_rate;
  //
  int retval = 0;
  retval += output_energy(time, p_object);
  retval += output_objects(time, p_object);
  return retval;
}

const logger_t logger = {
  .init          = init,
  .finalise      = finalise,
  .get_next_time = get_next_time,
  .execute       = execute,
};

