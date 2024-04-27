import numpy
import integrator

def integrate(m, l, g, dt, tmax, kernel):
    vel = 0.
    pos = 0.
    t = 0.
    while t < tmax:
        vel, pos = kernel(m, l, g, dt, vel, pos)
        t += dt
    return pos

# Richardson extrapolation
# ref: https://na-inet.jp/nasoft/chap17.pdf
def extrapolate(xs, ys):
    n = len(xs)
    r = numpy.zeros((n, n))
    for i in range(n):
        r[i, 0] = ys[i]
    for i in range(1, n):
        for j in range(1, i + 1):
            r[i, j] = r[i, j - 1] + (r[i, j - 1] - r[i - 1, j - 1]) / (xs[i - j] / xs[i] - 1.)
    return r[n - 1, n - 1]

def try_different_dt(m, l, g, dts, tmax, fname, kernel):
    last_vels = list()
    for dt in dts:
        last_vel = integrate(m, l, g, dt, tmax, kernel)
        last_vels.append(last_vel)
    theory = extrapolate(dts, last_vels)
    last_vels = numpy.abs(last_vels - theory)
    numpy.savetxt(fname, numpy.transpose([dts, last_vels]), fmt="% .15e", delimiter=" ")

def main():
    m = +1.
    l = +1.
    g = -1.
    tmax = 1.e+1
    dts = 1. / numpy.array([16., 32., 64., 128., 256., 512., 1024., 2048.])
    try_different_dt(m, l, g, dts, tmax, "data/convergence_euler_explicit_euler_implicit.dat", integrator.kernel2)
    try_different_dt(m, l, g, dts, tmax, "data/convergence_crank_nicolson_crank_nicolson.dat", integrator.kernel3)
    try_different_dt(m, l, g, dts, tmax, "data/convergence_energy_conserving.dat",             integrator.kernel6)

if __name__ == "__main__":
    main()
