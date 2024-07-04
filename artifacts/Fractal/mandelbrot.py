# MIT License
#
# Copyright (c) 2023 NaokiHori
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# USAGE
# python3 mandelbrot.py

import numpy as np
from matplotlib import pyplot

def is_inside_mandelbrot(x, y):
    c = x + y * 1j
    cnt_max = 32
    z = 0. + 0.j
    for cnt in range(cnt_max):
        z = np.power(z, 2.) + c
        if np.abs(z) > 2.:
            return False
    return True

def mandelbrot():
    corner_x = -2.5
    corner_y = -1.5
    Nx = 800
    Ny = 600
    grid_size = 5.e-3
    array = np.full((Ny, Nx), False)
    for m in range(Ny):
        y = corner_y + 1. * m * grid_size
        for n in range(Nx):
            x = corner_x + 1. * n * grid_size
            array[m, n] = is_inside_mandelbrot(x, y)
    pyplot.contourf(array)
    pyplot.show()
    pyplot.close()

mandelbrot()
