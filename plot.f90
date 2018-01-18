!----------------------------------------------------------------------
!                        DATA & PLOTTING                              !
!----------------------------------------------------------------------

MODULE plot

    USE parameters, ONLY : dt, npart, box, num_b, size_b, norm_b, box2, dimnsn

    IMPLICIT NONE

CONTAINS



    SUBROUTINE plot_energy_init()

        INTEGER :: u

        OPEN(FILE="data/energy.dat",NEWUNIT=u,ACTION="write",STATUS="replace")
        WRITE(u,'(4(A10))') "time", "tot", "pot", "kin"
        CLOSE(u)

    END SUBROUTINE plot_energy_init


    SUBROUTINE plot_energy_data(k, en, pot_en, kin_en)

        INTEGER, INTENT(IN) :: k
        REAL, INTENT(IN) :: en, pot_en, kin_en

        INTEGER :: u

        OPEN(FILE="data/energy.dat",NEWUNIT=u,ACTION="write",POSITION="append")
        WRITE(u,'(4(F10.5))') k*dt, en, pot_en, kin_en
        CLOSE(u)

    END SUBROUTINE plot_energy_data


    SUBROUTINE plot_energy()

        INTEGER :: u

        OPEN(FILE="data/gnuplot.style",NEWUNIT=u,ACTION="write",STATUS="replace") ! energy plot
        WRITE(u,*)"set terminal pngcairo"
        WRITE(u,*)"set output 'energy.png'"
        !WRITE(u,*)"set terminal epslatex color"
        !WRITE(u,*)"set terminal cairolatex pdf color"
        !WRITE(u,*)"set format '$%g$'"
        !WRITE(u,*)"set output 'energy.tex'"
        WRITE(u,*)"set key outside"
        WRITE(u,*)"set xlabel 'time'"
        WRITE(u,*)"set ylabel 'energy'"
        WRITE(u,*)"plot for [col=2:4] 'energy.dat' using 1:col with lines title columnheader"
        CLOSE(u)
        CALL EXECUTE_COMMAND_LINE("cd data/ && gnuplot gnuplot.style")

    END SUBROUTINE plot_energy



    SUBROUTINE plot_positions(k, temp, en, Rx, Ry, Rz, Vsq)

        ! Particle positions in a box

        INTEGER, INTENT(IN) :: k
        REAL, INTENT(IN) :: temp, en
        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Rz, Vsq

        INTEGER :: i, u

        OPEN(FILE="data/Position/pos.dat",NEWUNIT=u,ACTION="write",STATUS="replace")! Write particle positions
        DO i = 1, npart
            WRITE(u,'(F10.5,3(",",F10.5))') Rx(i), Ry(i), Rz(i), Vsq(i)
        END DO
        CLOSE(u)

        OPEN(FILE="data/Position/gnuplot.style",NEWUNIT=u,ACTION="write",STATUS="replace")  ! Write the gnuplot style file
        WRITE(u,*)"set terminal pngcairo"
        WRITE(u,'(A12,I0.5,A5)')"set output '",k,".png'"
        !WRITE(u,*)"set terminal epslatex color"
        !WRITE(u,*)"set terminal cairolatex pdf color"
        !WRITE(u,*)"set format '$%g$'"
        !WRITE(u,'(A25)')"set output 'position.tex'"
        WRITE(u,*)"set key off"
        WRITE(u,*)"set border 4095"
        WRITE(u,'(A17,F0.3,A9,F0.0)')"set label 'Temp =",temp,"' at 0,0,",box+(0.5*box)
        WRITE(u,'(A19,F0.4,A9,F0.0)')"set label 'Energy =",en,"' at 0,0,",box+(0.7*box)
        WRITE(u,'(A14,F0.0,A1)')"set xrange [0:",box,"]"
        WRITE(u,'(A14,F0.0,A1)')"set yrange [0:",box,"]"
        WRITE(u,'(A14,F0.0,A1)')"set zrange [0:",box,"]"
        !WRITE(u,*)"set cbrange [0:1]"
        WRITE(u,'(A75)')"splot 'pos.dat' using 1:2:3:4 with points palette pointsize 1.2 pointtype 7"
        CLOSE(u)

        CALL EXECUTE_COMMAND_LINE("cd data/Position/ && gnuplot gnuplot.style")          ! Plotting and Saving the image


    END SUBROUTINE plot_positions


    SUBROUTINE plot_positions_KobAnd(k, temp, en, Rx, Ry, Rz, species)

        ! Particle positions in a box

        INTEGER, INTENT(IN) :: k
        REAL, INTENT(IN) :: temp, en

        LOGICAL, DIMENSION(:) , INTENT(IN) :: species
        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Rz

        INTEGER :: i, u, v
        INTEGER, DIMENSION(npart) :: colour

        WHERE (species)
            colour = 2
        ELSEWHERE
            colour = 5
        ENDWHERE

        OPEN(FILE="data/Position/pos.dat",NEWUNIT=u,ACTION="write",STATUS="replace")! Write particle positions
        DO i = 1, npart
            WRITE(u,'(F10.5,2(",",F10.5),(",",I6.1))') Rx(i), Ry(i), Rz(i), colour(i)
        END DO
        CLOSE(u)

        OPEN(FILE="data/Position/gnuplot.style",NEWUNIT=v,ACTION="write",STATUS="replace")  ! Write the gnuplot style file
        WRITE(v,*)"set terminal pngcairo"
        WRITE(v,'(A12,I0.5,A5)')"set output '",k,".png'"
        !WRITE(v,*)"set terminal epslatex color"
        !WRITE(u,*)"set terminal cairolatex pdf color"
        !WRITE(v,*)"set format '$%g$'"
        !WRITE(v,'(A25)')"set output 'position.tex'"
        WRITE(v,*)"set key off"
        WRITE(v,*)"set border 4095"
        WRITE(v,'(A17,F0.3,A9,F0.0)')"set label 'Temp =",temp,"' at 0,0,",box+(0.5*box)
        WRITE(v,'(A19,F0.4,A9,F0.0)')"set label 'Energy =",en,"' at 0,0,",box+(0.7*box)
        WRITE(v,'(A14,F0.0,A1)')"set xrange [0:",box,"]"
        WRITE(v,'(A14,F0.0,A1)')"set yrange [0:",box,"]"
        WRITE(v,'(A14,F0.0,A1)')"set zrange [0:",box,"]"
        WRITE(v,*)"set cbrange [0:9]"
        WRITE(v,*)"unset colorbox"
        WRITE(v,'(A75)')"splot 'pos.dat' using 1:2:3:4 with points palette pointsize 1.2 pointtype 7"
        CLOSE(v)

        CALL EXECUTE_COMMAND_LINE("cd data/Position/ && gnuplot gnuplot.style")          ! Plotting and Saving the image


    END SUBROUTINE plot_positions_KobAnd



    SUBROUTINE plot_positions_2D(k, temp, en, Rx, Ry, Vsq)

        ! Particle positions in a 2D box

        INTEGER, INTENT(IN) :: k
        REAL, INTENT(IN) :: temp, en
        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Vsq

        INTEGER :: i, u, v

        OPEN(FILE="data/Position/pos.dat",NEWUNIT=u,ACTION="write",STATUS="replace")! Write particle positions
        DO i = 1, npart
            WRITE(u,'(F10.5,2(",",F10.5))') Rx(i), Ry(i), Vsq(i)
        END DO
        CLOSE(u)

        OPEN(FILE="data/Position/gnuplot.style",NEWUNIT=v,ACTION="write",STATUS="replace")  ! Write the gnuplot style file
        WRITE(v,*)"set terminal pngcairo"
        WRITE(v,'(A12,I0.5,A5)')"set output '",k,".png'"
        !WRITE(v,*)" set terminal epslatex color"
        !WRITE(u,*)"set terminal cairolatex pdf color"
        !WRITE(v,*)" set format '$%g$'"
        !WRITE(v,'(A29)')" set output 'position_2D.tex'"
        WRITE(v,'(A17,F0.3,A7,F0.0)')"set label 'Temp =",temp,"' at 0,",box+(0.05*box)
        WRITE(v,'(A19,F0.4,A7,F0.0)')"set label 'Energy =",en,"' at 0,",box+(0.1*box)
        WRITE(v,'(A14,F0.0,A1)')"set xrange [0:",box,"]"
        WRITE(v,'(A14,F0.0,A1)')"set yrange [0:",box,"]"
        !WRITE(v,*)"set cbrange [0:1]"
        WRITE(v,'(A75)')"plot 'pos.dat' using 1:2:3 with points palette pointsize 0.85 pointtype 7"
        CLOSE(v)

        CALL EXECUTE_COMMAND_LINE("cd data/Position/ && gnuplot gnuplot.style")          ! Plotting and Saving the image


    END SUBROUTINE plot_positions_2D




    SUBROUTINE plot_radial_dist(k, g)

        ! The radial distribution function

        INTEGER, INTENT(IN) :: k
        REAL, DIMENSION(:), INTENT(IN) :: g

        INTEGER :: i, u, v
        REAL :: r_b, vol_b, gr

        OPEN(FILE="data/RadDistFunc/RadDistFunc.dat",NEWUNIT=u,ACTION="write",STATUS="replace")
        DO CONCURRENT (i = 1 : num_b)                                                   ! normalizing g(r)
            r_b = (i-0.5)*size_b
            vol_b = norm_b * ( (i**dimnsn) - ((i-1)**dimnsn) )
            gr = g(i) / vol_b
            WRITE(u,'(F7.3,",",F15.5)') r_b, gr
        END DO
        CLOSE(u)

        OPEN(FILE="data/RadDistFunc/gnuplot.style",NEWUNIT=v,ACTION="write",STATUS="replace")  ! Write the gnuplot style file
        WRITE(v,*)"set terminal pngcairo"
        WRITE(v,'(A12,I0.5,A5)')"set output '",k,".png'"
        !WRITE(v,*)"set terminal epslatex color"
        !WRITE(u,*)"set terminal cairolatex pdf color"
        !WRITE(v,*)"set format '$%g$'"
        !WRITE(v,'(A25)')"set output 'radist_2D.tex'"
        WRITE(v,*)"set key off"
        WRITE(v,*)"set xlabel 'r' "
        WRITE(v,*)"set ylabel 'g(r)' "
        WRITE(v,*)"set style fill solid 0.5 "
        WRITE(v,'(A15,F0.0,A1)')" set xrange [0:",box2,"]"
        WRITE(v,*)"plot 'RadDistFunc.dat' with boxes "
        CLOSE(v)

        CALL EXECUTE_COMMAND_LINE("cd data/RadDistFunc/ && gnuplot gnuplot.style")          ! Plotting and Saving the image

    END SUBROUTINE plot_radial_dist



END MODULE plot
