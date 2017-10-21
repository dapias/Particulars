#! /bin/bash

echo ""
echo "Beginning on " ; date

list=(1000 5000 10000)

for i in ${list[*]} ; do

    echo ""
    echo "$i particles"
    echo ""

    sed -i '/npart=/c\    INTEGER, PARAMETER :: npart='$i'    ! Number of particles' parameters.f90

    #ifort -O3 -xHost -real-size 64 -qopenmp -assume bscc -assume realloc_lhs parameters.f90 random.f90 lattice.f90 thermostat.f90 force.f90 plot.f90 integration.f90 diagnostics.f90 main_cell_2D.f90
    gfortran -O2 -march=native -fdefault-real-8 -fopenmp -fbackslash parameters.f90 random.f90 lattice.f90 thermostat.f90 force.f90 plot.f90 integration.f90 diagnostics.f90 main_cell_2D.f90

    ./a.out

done

echo ""
echo "Ending on " ; date
echo ""

#notify-send "The program has finished running."
