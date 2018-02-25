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

#list=(0.05 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0 1.1 1.2 1.3 1.4 1.5)
list=(1000 5000 10000)

for i in ${list[*]} ; do

    echo ""
    echo "$i particles"
    echo ""

    #sed -i '/rho=/c\    REAL, PARAMETER :: rho='$i'          ! Particle density' parameters.f90
    sed -i '/npart=/c\    INTEGER, PARAMETER :: npart='$i'    ! Number of particles' parameters.f90

    ifort -O3 -xHost -real-size 64 -assume bscc -assume realloc_lhs parameters.f90 random.f90 lattice.f90 thermostat.f90 force.f90 plot.f90 integration.f90 diagnostics.f90 main_cell_2D.f90

    ./a.out

done

echo ""
echo "Ending on " ; date
echo ""

#notify-send "The program has finished running."
