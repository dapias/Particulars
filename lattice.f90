!-------------------------------------------------------------------------------!
!                                     LATTICES                                  !
!-------------------------------------------------------------------------------!

MODULE lattice

    USE parameters, ONLY : npart, box, npart_edge, npart_face, space, num_cell,&
        dimnsn, npart_cell, npart_d_cell
    USE random_init

    IMPLICIT NONE


CONTAINS



    SUBROUTINE lattice_random(Rx, Ry, Rz)

        ! A box of randomly placed partices

        REAL, DIMENSION(:), INTENT(OUT) :: Rx, Ry, Rz

        CALL init_random_seed()

        CALL RANDOM_NUMBER(Rx)
        Rx = box*Rx
        CALL RANDOM_NUMBER(Ry)
        Ry = box*Ry
        CALL RANDOM_NUMBER(Rz)
        Rz = box*Rz


    END SUBROUTINE lattice_random


    SUBROUTINE lattice_simple_cubic(Rx, Ry, Rz)

        ! A box of particles placed on a simple cubic lattice

        REAL, DIMENSION(:), INTENT(OUT) :: Rx, Ry, Rz

        INTEGER :: i, x, y, z

        x=0 ; y=0 ; z=0

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y,z)
        !$OMP DO
        DO i = 1, npart
            x = MODULO(i-1, npart_edge)
            Rx(i) = space*x
            y = MODULO( ( (i-1) / npart_edge ), npart_edge)
            Ry(i) = space*y
            z = (i-1) / npart_face
            Rz(i) = space*z
        END DO
        !$OMP END DO
        !$OMP END PARALLEL


    END SUBROUTINE lattice_simple_cubic



    SUBROUTINE lattice_cell_simple_cubic(R, location)

        ! particles on a lattice, cell-wise

        REAL, DIMENSION(:), INTENT(OUT) :: R
        INTEGER, DIMENSION(:), INTENT(OUT) :: location

        INTEGER :: i

        location(1) = 1

        DO CONCURRENT (i = 1 : (num_cell-1))
            location(i) = 1 + (i * npart_d_cell)
        END DO

        R = 0

    END SUBROUTINE



    SUBROUTINE lattice_square(Rx, Ry)

        ! A box of particles placed on a square lattice

        REAL, DIMENSION(:), INTENT(OUT) :: Rx, Ry

        INTEGER :: i, x, y

        x=0 ; y=0

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(x,y)
        !$OMP DO
        DO i = 1, npart
            x = MODULO(i-1, npart_edge)
            Rx(i) = space*x
            y = MODULO( ( (i-1) / npart_edge ), npart_edge)
            Ry(i) = space*y
        END DO
        !$OMP END DO
        !$OMP END PARALLEL


    END SUBROUTINE lattice_square



    SUBROUTINE lattice_speciate(species)

        ! particle species assignment

        LOGICAL, DIMENSION(:), INTENT(OUT) :: species

        REAL, DIMENSION(npart) :: r

        CALL init_random_seed()

        CALL RANDOM_NUMBER(r)

        WHERE (r .LE. 0.2)
            species = .TRUE.
        ELSEWHERE
            species = .FALSE.
        ENDWHERE

    END SUBROUTINE lattice_speciate



    SUBROUTINE lattice_velocities_random(Vx, Vy, Vz)

        REAL, DIMENSION(:), INTENT(OUT) :: Vx, Vy, Vz

        REAL :: sumVx, sumVy, sumVz

        CALL init_random_seed()

        CALL RANDOM_NUMBER(Vx)
        Vx = Vx - 0.5
        sumVx = SUM(Vx)/REAL(npart)
        CALL RANDOM_NUMBER(Vy)
        Vy = Vy - 0.5
        sumVy = SUM(Vy)/REAL(npart)
        CALL RANDOM_NUMBER(Vz)
        Vz = Vz - 0.5
        sumVz = SUM(Vz)/REAL(npart)

        Vx = Vx - sumVx                                                         ! setting Vcm to zero
        Vy = Vy - sumVy
        Vz = Vz - sumVz

    END SUBROUTINE lattice_velocities_random


    SUBROUTINE lattice_velocities_random_2D(Vx, Vy)

        REAL, DIMENSION(:), INTENT(OUT) :: Vx, Vy

        REAL :: sumVx, sumVy

        CALL init_random_seed()

        CALL RANDOM_NUMBER(Vx)
        Vx = Vx - 0.5
        sumVx = SUM(Vx)/REAL(npart)
        CALL RANDOM_NUMBER(Vy)
        Vy = Vy - 0.5
        sumVy = SUM(Vy)/REAL(npart)

        Vx = Vx - sumVx                                                         ! setting Vcm to zero
        Vy = Vy - sumVy

    END SUBROUTINE lattice_velocities_random_2D



END MODULE lattice
