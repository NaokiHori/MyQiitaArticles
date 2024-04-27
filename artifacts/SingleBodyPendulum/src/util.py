import numpy

def compute_energy(m, l, g, vel, pos):
    e_t = 0.5 * m * l**2. * vel**2.
    e_u = 1. - m * g * l * numpy.sin(pos)
    return e_t + e_u

def compute_phase_space(m, l, vel, pos):
    p = m * l**2. * vel
    q = pos
    return p, q

