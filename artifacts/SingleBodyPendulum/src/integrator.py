import numpy

niter = 20

def kernel0(m, l, g, dt, vel, pos):
    # vel: Euler-explicit
    # pos: Euler-explicit
    dvel = dt * g / l * numpy.cos(pos)
    vel_new = vel + dvel
    dpos = dt * vel
    pos_new = pos + dpos
    return vel_new, pos_new

def kernel1(m, l, g, dt, vel, pos):
    # vel: Euler-explicit
    # pos: Crank-Nicolson
    dvel = dt * g / l * numpy.cos(pos)
    vel_new = vel + dvel
    dpos = dt * (0.5 * vel + 0.5 * vel_new)
    pos_new = pos + dpos
    return vel_new, pos_new

def kernel2(m, l, g, dt, vel, pos):
    # vel: Euler-explicit
    # pos: Euler-implicit
    dvel = dt * g / l * numpy.cos(pos)
    vel_new = vel + dvel
    dpos = dt * vel_new
    pos_new = pos + dpos
    return vel_new, pos_new

def kernel3(m, l, g, dt, vel, pos):
    # vel: Crank-Nicolson
    # pos: Crank-Nicolson
    dpos = 0.
    for _ in range(niter):
        pos_new = pos + dpos
        dvel = dt * g / l * (0.5 * numpy.cos(pos) + 0.5 * numpy.cos(pos_new))
        vel_new = vel + dvel
        dpos = dt * (0.5 * vel + 0.5 * vel_new)
    return vel_new, pos_new

def kernel4(m, l, g, dt, vel, pos):
    # vel: Crank-Nicolson
    # pos: Euler-implicit
    dpos = 0.
    for _ in range(niter):
        pos_new = pos + dpos
        dvel = dt * g / l * (0.5 * numpy.cos(pos) + 0.5 * numpy.cos(pos_new))
        vel_new = vel + dvel
        dpos = dt * vel_new
    return vel_new, pos_new

def kernel5(m, l, g, dt, vel, pos):
    # vel: Euler-implicit
    # pos: Euler-implicit
    dpos = 0.
    for _ in range(niter):
        pos_new = pos + dpos
        dvel = dt * g / l * numpy.cos(pos_new)
        vel_new = vel + dvel
        dpos = dt * vel_new
    return vel_new, pos_new

def kernel6(m, l, g, dt, vel, pos):
    def mysinc(x):
        # sin(x) / x
        # NOTE: original numpy.sinc is
        #   normalised differently
        return numpy.sinc(x / numpy.pi)
    # vel: Energy-conserving
    # pos: Energy-conserving
    dpos = 0.
    for _ in range(niter):
        pos_new = pos + dpos
        dvel = dt * g / l \
            * numpy.cos(0.5 * pos + 0.5 * pos_new) \
            * mysinc(0.5 * pos_new - 0.5 * pos)
        vel_new = vel + dvel
        dpos = dt * (0.5 * vel + 0.5 * vel_new)
    return vel_new, pos_new

