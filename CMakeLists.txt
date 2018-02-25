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

# THis cmake file is based on:
# http://fortranwiki.org/fortran/show/Build+tools
# https://cmake.org/Wiki/CMakeForFortranExample

cmake_minimum_required(VERSION 3.0)

project(Particulars)
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(checks "-check")
    set(warnings "-warn -diag-enable=all")
    set(debug "-O0 -g")
    set(optimize "-xHost -O3")
    set(openmp "-qopenmp")
    #set(lapack_include "-I./lapack95_modules_ifort")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    #set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(dialect "-ffree-form -fimplicit-none")
    set(checks "-fcheck=all")
    set(warnings "-Wall -Wextra -pedantic")
    set(debug "-O0 -g")
    set(optimize "-march=native -O3")
    set(openmp "-fopenmp")
    #set(lapack_include "-I./lapack95_modules_gfort")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(checks "-C")
endif()

set(CMAKE_Fortran_FLAGS "${dialect}")
set(CMAKE_Fortran_FLAGS_DEBUG "${debug} ${checks} ${warnings}")
set(CMAKE_Fortran_FLAGS_RELEASE "${optimize}")

#
# Compile.
#

set(sources parameters.f90 random_init.F90 lattice.f90 thermostat.f90 force.f90 plot.f90 integration.f90 diagnostics.f90 main_cell_pack_3D_KobAnd.f90)
add_executable(prog ${sources})

#find_package(LAPACK REQUIRED)
#if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
#    target_link_libraries(prog -lblas -llapack ./lapack95_modules_ifort/lapack95.a)
#elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
#    target_link_libraries(prog -lblas -llapack ./lapack95_modules_gfort/lapack95.a)
#endif()
