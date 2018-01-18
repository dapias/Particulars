!-----------------------------------------------------------------------------------------------!
!                                              README                                           !
!-----------------------------------------------------------------------------------------------!

 * All module file names indicate their function.
 * Procedures names begin with the name of the module thay belong to.
 * Use one 'main_' file. It calls procedures from the modules.
 * Unless specified, all procedures work in 3 dimensions.
 * In the parameters file, you may switch from density control to box side length control

!-----------------------------------------------------------------------------------------------!

 COMPILE TIME OPTIONS:

 Add module files before the main Program in order of increasing dependencies (i.e. module with no deps first ...)

 Use the following compiler option to change to double precision:
 intel :   -real-size 64
 gnu   :   -fdefault-real-8

 Use the following compiler option to change to prevent integer overflow:
 intel :   -integer-size 64
 gnu   :   -fdefault-integer-8

 Use the following compiler option if you want the 'percentage complete' indication. Else, comment out relevant code:
 intel :   -assume bscc
 gnu   :   -fbackslash

 Add the following OpenMP compiler options only if the supported OMP is .GE. 3.1 . This due to array REDUCTION support
 intel :   -qopenmp
 gnu   :   -fopenmp

 Use the following in ifort, to enable allocation upon assignment (no need to allocate/deallocate). Not necessary in gfortran
 intel:    -assume realloc_lhs

 To use the apptopriate memory model use the following flags:
 -mcmodel=<x>
   where <x> : small(default), medium, large

 To optimize for speed, use the following:
 intel :   -O3 -xHost
    very agressive: -fast
 gnu   :   -O3 -march=native

 Set OpenMP environment variables in .bashrc as follows:
 export OMP_PLACES=<place>
   where <place> : threads, cores, sockets
 export OMP_NUM_THREADS=<number>
   where <number> is less than or equal to number of cores

 Use the following compiler option to use the pre-processor:
 intel :   -fpp
 gnu   :   -cpp
 NOT NECESSARY if the relevant file has the '.F90' extension

 Use the following during initial compilations, to see all generated warnings:
 intel	: -check -warn -diag-enable=all
 gnu	: -Wall -Wextra

!-----------------------------------------------------------------------------------------------!

 CURRENTLY USED OPTIONS:

 ifort -I./lapack95_modules_ifort -O3 -xHost -real-size 64 -qopenmp -assume bscc -assume realloc_lhs parameters.f90 random_init.F90 lattice.f90 thermostat.f90 force.f90 plot.f90 integration.f90 diagnostics.f90 main.f90 -lblas -llapack ./lapack95_modules_ifort/lapack95.a

 gfortran -I./lapack95_modules_gfort -O3 -march=native -fdefault-real-8 -fopenmp -fbackslash parameters.f90 random_init.F90 lattice.f90 thermostat.f90 force.f90 plot.f90 integration.f90 diagnostics.f90 main.f90 -lblas -llapack ./lapack95_modules_gfort/lapack95.a

!-----------------------------------------------------------------------------------------------!

 MISCELLANEOUS:

 Use the following to prevent hang up signals (due to logouts, terminal closures etc.) from killing the program:
 nohup <compile command> &

 You may make videos using:
 ffmpeg -framerate 7 -pattern_type glob -i '*.png' -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4
UNTESTED:
 ffmpeg -framerate 7 -pattern_type glob -i '*.png' -c:v libvpx-vp9 -b:v 1M -threads 8 -speed 1 -tile-columns 6 -frame-parallel 1 -f webm out.webm

 Create and update archived backups with the following:
   c:create, x:extract, u:update, v:verbose, p:permissions, f:file_name
 tar -cvpf code.tar *.f90 *.sh *.txt
 tar -xvpf code.tar

!-----------------------------------------------------------------------------------------------!

   TO DO:

 - convert velocity arrays in MODULE thermostat to INTENT(OUT)
 - examine using 2D arrays instead of separate arrays for each dimension
 - convert from 3 - 1D arrays to assumed shape arrays with the CONTIGUOUS attribute
 - look for possiblities of DO CONCURRENT
 - use 'set term tikz latex' in plot for LaTeX


!-----------------------------------------------------------------------------------------------!
