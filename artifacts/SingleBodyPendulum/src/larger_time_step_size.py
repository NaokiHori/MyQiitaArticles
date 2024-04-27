import integrator
import util

def integrate(m, l, g, dt, tmax, logname, kernel):
    vel = 0.
    pos = 0.
    # save initial kinetc / potential energies
    e0 = util.compute_energy(m, l, g, vel, pos)
    with open(logname, "w") as f:
        t = 0.
        while t < tmax:
            # itegrate one time step
            vel, pos = kernel(m, l, g, dt, vel, pos)
            t += dt
            # check current energy
            e1 = util.compute_energy(m, l, g, vel, pos)
            # normalised energy
            ne = e1 / e0
            p, q = util.compute_phase_space(m, l, vel, pos)
            f.write(f"{t: .15e} {ne: .15e} {p: .15e} {q: .15e}\n")
    return

def main():
    m = +1.
    l = +1.
    g = -1.
    dt = 1.e-1
    tmax = 1.e+2
    integrate(m, l, g, dt, tmax, "data/larger_time_step_size_crank_nicolson_crank_nicolson.dat", integrator.kernel3)
    integrate(m, l, g, dt, tmax, "data/larger_time_step_size_energy_conserving.dat",             integrator.kernel6)

if __name__ == "__main__":
    main()
