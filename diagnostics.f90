!----------------------------------------------------------------------
!                            DIAGNOSTICS                              !
!----------------------------------------------------------------------

MODULE diagnostics

    USE parameters, ONLY : npart, num_cell_edge

CONTAINS



    SUBROUTINE diagnostics_pos_vel(Vsq, Vx, Vy, Vz, Rx, Ry, Rz, Fx, Fy, Fz, ptemp, temp)

        ! Writes out a diagnostics file incase of an error.

        REAL, INTENT(IN) :: ptemp, temp
        REAL, DIMENSION(:), INTENT(IN) :: Vsq, Vx, Vy, Vz, Rx, Ry, Rz, Fx, Fy, Fz

        INTEGER :: i
        REAL :: sumVx, sumVy, sumVz

        sumVx = 0 ; sumVy = 0 ; sumVz = 0

        sumVx = SUM(Vx)/REAL(npart)
        sumVy = SUM(Vy)/REAL(npart)
        sumVz = SUM(Vz)/REAL(npart)

        OPEN(UNIT=3,FILE="data/error.dat",ACTION="write",STATUS="replace")
        WRITE(3,'(A15,9(",",A15))') "Vsq", "Rx", "Vx", "Fx", "Ry", "Vy", "Fy", "Rz", "Vz", "Fz"
        DO i = 1, npart
            WRITE(3,'(EN15.5,9(",",F15.5))') Vsq(i), Rx(i), Vx(i), Fx(i), Ry(i), Vy(i), Fy(i), Rz(i), Vz(i), Fz(i)
        END DO
        WRITE(3,*) ""
        WRITE(3,*) "previous temp. =", ptemp, "  and   current temp. = ", temp
        WRITE(3,*) ""
        WRITE(3,*) "sumVx:", sumVx, "sumVy:", sumVy, "sumVz:", sumVz
        CLOSE(3)

        ERROR STOP "*** Annomalous fluctuation detected ; Error file written ***"

    END SUBROUTINE diagnostics_pos_vel


    SUBROUTINE diagnostics_pos_vel_2D(Vsq, Vx, Vy, Rx, Ry, Fx, Fy, ptemp, temp)

        ! Writes out a diagnostics file incase of an error.

        REAL, INTENT(IN) :: ptemp, temp
        REAL, DIMENSION(:), INTENT(IN) :: Vsq, Vx, Vy, Rx, Ry, Fx, Fy

        INTEGER :: i
        REAL :: sumVx, sumVy

        sumVx = 0 ; sumVy = 0

        sumVx = SUM(Vx)/REAL(npart)
        sumVy = SUM(Vy)/REAL(npart)

        OPEN(UNIT=3,FILE="data/error.dat",ACTION="write",STATUS="replace")
        WRITE(3,'(A15,9(",",A15))') "Vsq", "Rx", "Vx", "Fx", "Ry", "Vy", "Fy"
        DO i = 1, npart
            WRITE(3,'(EN15.5,6(",",F15.5))') Vsq(i), Rx(i), Vx(i), Fx(i), Ry(i), Vy(i), Fy(i)
        END DO
        WRITE(3,*) ""
        WRITE(3,*) "previous temp. =", ptemp, "  and   current temp. = ", temp
        WRITE(3,*) ""
        WRITE(3,*) "sumVx:", sumVx, "sumVy:", sumVy
        CLOSE(3)

        ERROR STOP "*** Annomalous fluctuation detected ; Error file written ***"

    END SUBROUTINE diagnostics_pos_vel_2D


    SUBROUTINE diagnostics_NCE_check()

        ! A check for cell-list-based algorithms to ensure that the number of cells along
        ! an edge is always 3 or more, because otherwise, the force calculation will see the
        ! same cell as being adjacent in both directions.
        !   Also, at such low numbers, cell-lists are an unnecessary overhead

        IF (num_cell_edge .LT. 3) THEN
            ERROR STOP "*** Number of CELLs too low; Do NOT use cell-lists ***"
        END IF

    END SUBROUTINE diagnostics_NCE_check


END MODULE diagnostics
