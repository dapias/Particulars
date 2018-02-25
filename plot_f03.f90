!----------------------------------------------------------------------
!                        DATA & PLOTTING                              !
!----------------------------------------------------------------------
!
!   Copyright (C) 2015-2018  Vishnu V. Krishnan : vishnugb@gmail.com
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <https://www.gnu.org/licenses/>.

MODULE plot

    USE parameters, ONLY : dt, npart, box, num_b, size_b, norm_b, box2, dimnsn

    IMPLICIT NONE

CONTAINS



    SUBROUTINE plot_energy_init()

        OPEN(UNIT=10,FILE="data/energy.dat",ACTION="write",STATUS="replace")
        WRITE(10,'(4(A10))') "time", "tot", "pot", "kin"
        CLOSE(10)

    END SUBROUTINE plot_energy_init


    SUBROUTINE plot_energy_data(k, en, pot_en, kin_en)

        INTEGER, INTENT(IN) :: k
        REAL, INTENT(IN) :: en, pot_en, kin_en

        OPEN(UNIT=10,FILE="data/energy.dat",ACTION="write",ACCESS="append")
        WRITE(10,'(4(F10.5))') k*dt, en, pot_en, kin_en
        CLOSE(10)

    END SUBROUTINE plot_energy_data


    SUBROUTINE plot_energy()

        OPEN(UNIT=11,FILE="data/gnuplot.style",ACTION="write",STATUS="replace") ! energy plot
        WRITE(11,*)"set terminal pngcairo"
        WRITE(11,*)"set output 'energy.png'"
        !WRITE(11,*)"set terminal epslatex color"
        !WRITE(11,*)"set format '$%g$'"
        !WRITE(11,*)"set output 'energy.tex'"
        WRITE(11,*)"set key outside"
        WRITE(11,*)"set xlabel 'time'"
        WRITE(11,*)"set ylabel 'energy'"
        WRITE(11,*)"plot for [col=2:4] 'energy.dat' using 1:col with lines title columnheader"
        CLOSE(11)
        CALL EXECUTE_COMMAND_LINE("cd data/ && gnuplot gnuplot.style")

    END SUBROUTINE plot_energy



    SUBROUTINE plot_positions(k, temp, en, Rx, Ry, Rz, Vsq)

        ! Particle positions in a box

        INTEGER, INTENT(IN) :: k
        REAL, INTENT(IN) :: temp, en
        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Rz, Vsq
        INTEGER :: i

        OPEN(UNIT=1,FILE="data/Position/pos.dat",ACTION="write",STATUS="replace")! Write particle positions
        DO i = 1, npart
            WRITE(1,'(F10.5,3(",",F10.5))') Rx(i), Ry(i), Rz(i), Vsq(i)
        END DO
        CLOSE(1)

        OPEN(UNIT=2,FILE="data/Position/gnuplot.style",ACTION="write",STATUS="replace")  ! Write the gnuplot style file
        WRITE(2,*)"set terminal pngcairo"
        WRITE(2,'(A12,I0.5,A5)')"set output '",k,".png'"
        !WRITE(2,*)"set terminal epslatex color"
        !WRITE(2,*)"set format '$%g$'"
        !WRITE(2,'(A25)')"set output 'position.tex'"
        WRITE(2,*)"set key off"
        WRITE(2,*)"set border 4095"
        WRITE(2,'(A17,F0.3,A9,F0.0)')"set label 'Temp =",temp,"' at 0,0,",box+(0.5*box)
        WRITE(2,'(A19,F0.4,A9,F0.0)')"set label 'Energy =",en,"' at 0,0,",box+(0.7*box)
        WRITE(2,'(A14,F0.0,A1)')"set xrange [0:",box,"]"
        WRITE(2,'(A14,F0.0,A1)')"set yrange [0:",box,"]"
        WRITE(2,'(A14,F0.0,A1)')"set zrange [0:",box,"]"
        !WRITE(2,*)"set cbrange [0:1]"
        WRITE(2,'(A75)')"splot 'pos.dat' using 1:2:3:4 with points palette pointsize 1.2 pointtype 7"
        CLOSE(2)

        CALL EXECUTE_COMMAND_LINE("cd data/Position/ && gnuplot gnuplot.style")          ! Plotting and Saving the image


    END SUBROUTINE plot_positions


    SUBROUTINE plot_positions_KobAnd(k, temp, en, Rx, Ry, Rz, species)

        ! Particle positions in a box

        INTEGER, INTENT(IN) :: k
        REAL, INTENT(IN) :: temp, en

        LOGICAL, DIMENSION(:) , INTENT(IN) :: species
        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Rz

        INTEGER :: i
        INTEGER, DIMENSION(npart) :: colour

        WHERE (species)
            colour = 2
        ELSEWHERE
            colour = 5
        ENDWHERE

        OPEN(UNIT=1,FILE="data/Position/pos.dat",ACTION="write",STATUS="replace")! Write particle positions
        DO i = 1, npart
            WRITE(1,'(F10.5,2(",",F10.5),(",",I6.1))') Rx(i), Ry(i), Rz(i), colour(i)
        END DO
        CLOSE(1)

        OPEN(UNIT=2,FILE="data/Position/gnuplot.style",ACTION="write",STATUS="replace")  ! Write the gnuplot style file
        WRITE(2,*)"set terminal pngcairo"
        WRITE(2,'(A12,I0.5,A5)')"set output '",k,".png'"
        !WRITE(2,*)"set terminal epslatex color"
        !WRITE(2,*)"set format '$%g$'"
        !WRITE(2,'(A25)')"set output 'position.tex'"
        WRITE(2,*)"set key off"
        WRITE(2,*)"set border 4095"
        WRITE(2,'(A17,F0.3,A9,F0.0)')"set label 'Temp =",temp,"' at 0,0,",box+(0.5*box)
        WRITE(2,'(A19,F0.4,A9,F0.0)')"set label 'Energy =",en,"' at 0,0,",box+(0.7*box)
        WRITE(2,'(A14,F0.0,A1)')"set xrange [0:",box,"]"
        WRITE(2,'(A14,F0.0,A1)')"set yrange [0:",box,"]"
        WRITE(2,'(A14,F0.0,A1)')"set zrange [0:",box,"]"
        WRITE(2,*)"set cbrange [0:9]"
        WRITE(2,*)"unset colorbox"
        WRITE(2,'(A75)')"splot 'pos.dat' using 1:2:3:4 with points palette pointsize 1.2 pointtype 7"
        CLOSE(2)

        CALL EXECUTE_COMMAND_LINE("cd data/Position/ && gnuplot gnuplot.style")          ! Plotting and Saving the image


    END SUBROUTINE plot_positions_KobAnd



    SUBROUTINE plot_positions_2D(k, temp, en, Rx, Ry, Vsq)

        ! Particle positions in a 2D box

        INTEGER, INTENT(IN) :: k
        REAL, INTENT(IN) :: temp, en
        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Vsq
        INTEGER :: i

        OPEN(UNIT=1,FILE="data/Position/pos.dat",ACTION="write",STATUS="replace")! Write particle positions
        DO i = 1, npart
            WRITE(1,'(F10.5,2(",",F10.5))') Rx(i), Ry(i), Vsq(i)
        END DO
        CLOSE(1)

        OPEN(UNIT=2,FILE="data/Position/gnuplot.style",ACTION="write",STATUS="replace")  ! Write the gnuplot style file
        WRITE(2,*)"set terminal pngcairo"
        WRITE(2,'(A12,I0.5,A5)')"set output '",k,".png'"
        !WRITE(2,*)" set terminal epslatex color"
        !WRITE(2,*)" set format '$%g$'"
        !WRITE(2,'(A29)')" set output 'position_2D.tex'"
        WRITE(2,'(A17,F0.3,A7,F0.0)')"set label 'Temp =",temp,"' at 0,",box+(0.05*box)
        WRITE(2,'(A19,F0.4,A7,F0.0)')"set label 'Energy =",en,"' at 0,",box+(0.1*box)
        WRITE(2,'(A14,F0.0,A1)')"set xrange [0:",box,"]"
        WRITE(2,'(A14,F0.0,A1)')"set yrange [0:",box,"]"
        !WRITE(2,*)"set cbrange [0:1]"
        WRITE(2,'(A75)')"plot 'pos.dat' using 1:2:3 with points palette pointsize 0.85 pointtype 7"
        CLOSE(2)

        CALL EXECUTE_COMMAND_LINE("cd data/Position/ && gnuplot gnuplot.style")          ! Plotting and Saving the image


    END SUBROUTINE plot_positions_2D




    SUBROUTINE plot_radial_dist(k, g)

        ! The radial distribution function

        INTEGER, INTENT(IN) :: k
        REAL, DIMENSION(:), INTENT(IN) :: g
        INTEGER :: i
        REAL :: r_b, vol_b, gr

        OPEN(UNIT=4,FILE="data/RadDistFunc/RadDistFunc.dat",ACTION="write",STATUS="replace")
        DO CONCURRENT (i = 1 : num_b)                                                   ! normalizing g(r)
            r_b = (i-0.5)*size_b
            vol_b = norm_b * ( (i**dimnsn) - ((i-1)**dimnsn) )
            gr = g(i) / vol_b
            WRITE(4,'(F7.3,",",F15.5)') r_b, gr
        END DO
        CLOSE(4)

        OPEN(UNIT=9,FILE="data/RadDistFunc/gnuplot.style",ACTION="write",STATUS="replace")  ! Write the gnuplot style file
        WRITE(9,*)"set terminal pngcairo"
        WRITE(9,'(A12,I0.5,A5)')"set output '",k,".png'"
        !WRITE(9,*)"set terminal epslatex color"
        !WRITE(9,*)"set format '$%g$'"
        !WRITE(9,'(A25)')"set output 'radist_2D.tex'"
        WRITE(9,*)"set key off"
        WRITE(9,*)"set xlabel 'r' "
        WRITE(9,*)"set ylabel 'g(r)' "
        WRITE(9,*)"set style fill solid 0.5 "
        WRITE(9,'(A15,F0.0,A1)')" set xrange [0:",box2,"]"
        WRITE(9,*)"plot 'RadDistFunc.dat' with boxes "
        CLOSE(9)

        CALL EXECUTE_COMMAND_LINE("cd data/RadDistFunc/ && gnuplot gnuplot.style")          ! Plotting and Saving the image

    END SUBROUTINE plot_radial_dist



END MODULE plot
