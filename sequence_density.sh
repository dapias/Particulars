#! /bin/bash

echo ""
echo "Beginning on " ; date

list=(0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5)

for i in ${list[*]} ; do

    echo ""
    echo "$i density"
    echo ""

    sed -i '/rho=/c\    REAL, PARAMETER :: rho='$i'          ! Particle density' parameters.f90

    #ifort -O3 -xHost -real-size 64 -qopenmp -assume bscc -assume realloc_lhs parameters.f90 random.f90 lattice.f90 thermostat.f90 force.f90 plot.f90 integration.f90 diagnostics.f90 main_cell_2D.f90
    gfortran -O2 -march=native -fdefault-real-8 -fopenmp -fbackslash parameters.f90 random.f90 lattice.f90 thermostat.f90 force.f90 plot.f90 integration.f90 diagnostics.f90 main_cell_2D.f90
    ./a.out

done

echo ""
echo "Ending on " ; date
echo ""

#notify-send "The program has finished running."
