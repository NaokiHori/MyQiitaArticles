// MIT License
//
// Copyright (c) 2023 NaokiHori
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

// USAGE
// cc -std=c99 -o a.out mandelbrot.c && ./a.out
// or
// cc -std=c99 -fopenmp -o a.out mandelbrot.c && ./a.out

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

typedef uint_fast16_t iter_t;
static const iter_t iter_t_max = UINT_FAST16_MAX;

typedef double real;

static int kernel(
    const iter_t maxiter,
    const real x0,
    const real y0,
    iter_t * iter
){
  // [xy]0: given   complex number (c)
  // [xy]1: current complex number (z^n)
  // [xy]2: next    complex number (z^{n+1})
  // z^0 = (0, 0)
  real x1 = 0.;
  real y1 = 0.;
  for(*iter = 1; *iter < maxiter; (*iter)++){
    const real x2 = x0 + x1 * x1 - y1 * y1;
    const real y2 = y0 + 2. * x1 * y1;
    const real z2 = x2 * x2 + y2 * y2;
    if(z2 > 4.){
      // diverged
      return 0;
    }
    x1 = x2;
    y1 = y2;
  }
  // not diverged (at least for now)
  return 0;
}

static int solve(
    const iter_t maxiter,
    const real delta,
    const size_t resols[2],
    const real center[2],
    iter_t * const iters
){
  const size_t nitems = resols[0] * resols[1];
#pragma omp parallel for
  for(size_t n = 0; n < nitems; n++){
    const size_t i = n % resols[0];
    const size_t j = n / resols[0];
    const real x0 = center[0] - 0.5 * delta * resols[0] + 1. * i * delta;
    const real y0 = center[1] - 0.5 * delta * resols[1] + 1. * j * delta;
    kernel(maxiter, x0, y0, iters + n);
  }
  return 0;
}

static int output(
    const char fname[],
    const size_t resols[2],
    const iter_t * iters
){
  uint_fast8_t * pixels = malloc(3 * resols[0] * resols[1] * sizeof(uint_fast8_t));
  // convert number of iterations to the corresponding color
  // check extrema and decide pixel colour
  iter_t miniter = iter_t_max;
  iter_t maxiter = 0;
  for(size_t n = 0; n < resols[0] * resols[1]; n++){
    const iter_t iter = iters[n];
    miniter = miniter > iter ? iter : miniter;
    maxiter = maxiter < iter ? iter : maxiter;
  }
#pragma omp parallel for
  for(size_t n = 0; n < resols[0] * resols[1]; n++){
    const iter_t iter = iters[n];
    real val = (1. * (iter - miniter)) / (1. * (maxiter - miniter));
    val = val < 0. ? 0. : val;
    val = val > 1. ? 1. : val;
    pixels[3 * n + 0] = (uint_fast8_t)(255. * val);
    pixels[3 * n + 1] = (uint_fast8_t)(255. * val);
    pixels[3 * n + 2] = (uint_fast8_t)(255. * val);
  }
  // write to file
  FILE * fp = fopen(fname, "w");
  if(NULL == fp){
    fprintf(stderr, "file open error: %s\n", fname);
    goto abort;
  }
  fprintf(fp, "P6\n%zu %zu\n255\n", resols[0], resols[1]);
  fwrite(pixels, sizeof(uint_fast8_t), 3 * resols[0] * resols[1], fp);
  fclose(fp);
abort:
  free(pixels);
  return 0;
}

int main(
    void
){
  // image resolution
  const size_t resols[2] = {2560, 1600};
  // maximum number of iteration
  // when reached, I consider the recurrence relation is converged
  // NOTE: should be smaller than 1 << 16 - 1
  const iter_t maxiter = 1024;
  // grid size
  const real delta = 1.e-8;
  // center of the image
  const real center[2] = {-7.431689909813917e-1, -1.263133777274951e-1};
  // solve the recurrence relation
  iter_t * iters = malloc(resols[0] * resols[1] * sizeof(iter_t));
  if(0 != solve(maxiter, delta, resols, center, iters)){
    goto abort;
  }
  // name of the resulting image, which should ends with ".ppm"
  const char fname[] = {"image.ppm"};
  if(0 != output(fname, resols, iters)){
    goto abort;
  }
abort:
  // clean-up
  free(iters);
  return 0;
}

