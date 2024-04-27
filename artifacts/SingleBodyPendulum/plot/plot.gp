{
  reset
  set terminal epslatex standalone color size 5.0,3.5 font ',12'
  set output 'explicit_velocity.tex'
  set xlabel 'Time'
  set ylabel 'Normalised energy'
  set xrange [0:100]
  set yrange [0.95:1.75]
  set xtics 0., 50., 100.
  set ytics 0.2
  set format x '$% .0f$'
  set format y '$% .1f$'
  set style line 1 lc rgb '#ff0000' lw 5
  set style line 2 lc rgb '#0000ff' lw 5
  set style line 3 lc rgb '#33aa00' lw 5
  set key left top
  plot \
    '../data/euler_explicit_euler_explicit.dat' u 1:2 every 10 t '$\omega$: E.-Exp., $\theta$: E.-Exp.' ls 1 w l, \
    '../data/euler_explicit_crank_nicolson.dat' u 1:2 every 10 t '$\omega$: E.-Exp., $\theta$: C.-N.'   ls 2 w l, \
    '../data/euler_explicit_euler_implicit.dat' u 1:2 every 10 t '$\omega$: E.-Exp., $\theta$: E.-Imp.' ls 3 w l
}

{
  reset
  set terminal epslatex standalone color size 5.0,3.5 font ',12'
  set output 'non_explicit_velocity.tex'
  set xlabel 'Time'
  set ylabel 'Normalised energy'
  set xrange [0:100]
  set yrange [0.35:1.05]
  set xtics 0., 50., 100.
  set ytics 0.2
  set format x '$% .0f$'
  set format y '$% .1f$'
  set style line 1 lc rgb '#ff0000' lw 5
  set style line 2 lc rgb '#0000ff' lw 5
  set style line 3 lc rgb '#33aa00' lw 5
  set key left bottom
  plot \
    '../data/crank_nicolson_crank_nicolson.dat' u 1:2 every 10 t '$\omega$: C.-N., $\theta$: C.-N.'     ls 1 w l, \
    '../data/crank_nicolson_euler_implicit.dat' u 1:2 every 10 t '$\omega$: C.-N., $\theta$: E.-Imp.'   ls 2 w l, \
    '../data/euler_implicit_euler_implicit.dat' u 1:2 every 10 t '$\omega$: E.-Imp., $\theta$: E.-Imp.' ls 3 w l
}

{
  reset
  set terminal epslatex standalone color size 5.0,3.5 font ',12'
  set output 'crank_nicolson_expanded.tex'
  set xlabel 'Time'
  set ylabel 'Normalised energy'
  yr = 5e-5
  set xrange [0:100]
  set yrange [1.-yr:1.+yr]
  set xtics 0., 50., 100.
  set ytics 1e-5
  set format x '$% .0f$'
  set style line 1 lc rgb '#ff0000' lw 5
  set key left top
  plot \
    '../data/crank_nicolson_crank_nicolson.dat' u 1:2 every 10 t '$\omega$: C.-N., $\theta$: C.-N.' ls 1 w l
}

{
  reset
  set terminal epslatex standalone color size 5.0,3.5 font ',12'
  set output 'energy_conserving_expanded.tex'
  set xlabel 'Time'
  set ylabel 'Normalised energy'
  yr = 1e-14
  set xrange [0:100]
  set yrange [-yr:+yr]
  set xtics 0., 50., 100.
  set ytics 1e-14
  set format x '$% .0f$'
  set format y '$% .1e$'
  set style line 1 lc rgb '#ff0000' lw 5
  set key left top
  plot \
    '../data/energy_conserving.dat' u 1:($2-1.) every 10 t 'Energy-conserving' ls 1 w l
}

{
  reset
  set terminal epslatex standalone color size 5.0,3.5 font ',12'
  set output 'compare_crank_nicolson_energy_conserving.tex'
  set xlabel 'Time'
  set ylabel 'Normalised energy'
  set xrange [0:100]
  set yrange [1e-16:1e0]
  set xtics 0., 50., 100.
  set ytics 1e4
  set logscale y
  set format x '$% .0f$'
  set format y '$10^{%L}$'
  set style line 1 lc rgb '#ff0000' lw 5
  set style line 2 lc rgb '#0000ff' lw 5
  set key left top
  plot \
    '../data/larger_time_step_size_crank_nicolson_crank_nicolson.dat' u 1:(abs(($2-1.))) t '$\omega$: C.-N., $\theta$: C.-N.' ls 1 w l, \
    '../data/larger_time_step_size_energy_conserving.dat' u 1:(abs(($2-1.))) t 'Energy-conserving' ls 2 w l
}

{
  reset
  set terminal epslatex standalone color size 5.0,3.5 font ',12'
  set output 'convergence.tex'
  set xlabel 'Time-step size'
  set ylabel 'Error'
  set xrange [1. / 2048. / sqrt(2.) : 1. / 16. * sqrt(2.)]
  set yrange [:]
  set xtics 1e1
  set ytics 1e2
  set logscale x
  set logscale y
  set format x '$10^{%L}$'
  set format y '$10^{%L}$'
  set style line 1 lc rgb '#ff0000' lw 5
  set style line 2 lc rgb '#0000ff' lw 5
  set style line 3 lc rgb '#33aa00' lw 5
  set style line 4 lc rgb '#000000' lw 5 dt 2
  set key right bottom
  plot \
    '../data/convergence_euler_explicit_euler_implicit.dat' u 1:2 t '$\omega$: E.-Exp., $\theta$: E.-Exp.' ls 1 pt 7 ps 1.5 w lp, \
    '../data/convergence_crank_nicolson_crank_nicolson.dat' u 1:2 t '$\omega$: C.-N., $\theta$: C.-N.' ls 2 pt 7 ps 1.5 w lp, \
    '../data/convergence_energy_conserving.dat' u 1:2 t 'Energy-conserving' ls 3 pt 7 ps 1.5 w lp, \
    x**1 notitle ls 4 w l, \
    x**2 notitle ls 4 w l
}
