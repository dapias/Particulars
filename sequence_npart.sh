#! /bin/bash
#
#   Copyright (C) 2015-2018  Vishnu V. Krishnan : vishnugb@gmail.com
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <https://www.gnu.org/licenses/>.

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
