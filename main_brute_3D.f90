!-------------------------------------------------------------------------------------!
!                                 MOLECULAR DYNAMICS                                  !
!-------------------------------------------------------------------------------------!
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

PROGRAM main_brute

    USE parameters, ONLY : n, npart, num_b
    USE lattice
    USE thermostat
    USE force
    USE plot
    USE integration
    USE diagnostics

    INTEGER :: k, AllocateStatus, initiate, terminate, clock_rate

    REAL :: temp, en, time_cpu, pot_en, ptemp, start, finish, kin_en, time_real

    REAL, DIMENSION(:), ALLOCATABLE :: Fx, Fy, Fz, Rx, Ry, Rz, Vx, Vy, Vz, Vsq
    REAL, DIMENSION(:), ALLOCATABLE :: g

    CALL CPU_TIME(start)
    CALL SYSTEM_CLOCK(initiate, clock_rate)

    CALL EXECUTE_COMMAND_LINE("rm data/*.png")
    CALL EXECUTE_COMMAND_LINE("rm data/Position/*.png")
    CALL EXECUTE_COMMAND_LINE("rm data/RadDistFunc/*.png")

    ALLOCATE(Rx(npart), Ry(npart), Rz(npart), Fx(npart), &
        Fy(npart), Fz(npart), Vx(npart), Vy(npart), Vz(npart), Vsq(npart), &
        g(num_b), STAT = AllocateStatus)

    IF (AllocateStatus /= 0) THEN
        STOP "*** Not enough memory ***"
    END IF


    CALL lattice_simple_cubic(Rx, Ry, Rz)

    CALL lattice_velocities_random(Vx, Vy, Vz)

    CALL thermostat_full(Vx, Vy, Vz)                                        ! rescale all vels using set_temp

    Vsq = ( (Vx*Vx) + (Vy*Vy) + (Vz*Vz) )/2.0                               ! normalized particle kinetic energies
    kin_en = SUM(Vsq)/REAL(npart)
    temp = (2*kin_en)/3.0
    !PRINT*, temp
    ptemp = temp


    CALL force_brute(Rx, Ry, Rz, Fx, Fy, Fz, pot_en, g)

    en = (pot_en/REAL(npart)) + kin_en                                      ! energy per particle

    k=0

    CALL plot_positions(k, temp, en, Rx, Ry, Rz, Vsq)

    CALL plot_radial_dist(k, g)

    CALL plot_energy_init()                                                 ! write titles of energy data file

    CALL thermostat_pid_init()                                              ! initialize PID variables


    DO k = 1, n                                                             ! multiple runs

        CALL integration_VelVerlet_pos(Rx, Ry, Rz, Vx, Vy, Vz, Fx, Fy, Fz)  ! positions from Velocity Verlet & PBC

        CALL integration_VelVerlet_vel(Vx, Vy, Vz, Fx, Fy, Fz)              ! first half-step of velocities from Velocity Verlet

        CALL force_brute(Rx, Ry, Rz, Fx, Fy, Fz, pot_en, g)                 ! force & g(r) calculation

        CALL integration_VelVerlet_vel(Vx, Vy, Vz, Fx, Fy, Fz)              ! second half-step of velocities from Velocity Verlet


        Vsq = (Vx*Vx + Vy*Vy + Vz*Vz)/2.0                                   ! normalized particle kinetic energies
        kin_en = SUM(Vsq)/REAL(npart)
        temp = (2.0*kin_en)/3.0                                             ! instantaneous temperature
        pot_en = pot_en/REAL(npart)
        en = pot_en + kin_en                                                ! energy per particle

        CALL plot_energy_data(k, en, pot_en, kin_en)

        IF (temp .gt. (100*ptemp)) THEN                                     ! excessive temp fluctuation check
            CALL plot_energy()
            PRINT*, ""
            CALL EXECUTE_COMMAND_LINE("fortune -o")
            PRINT*, ""
            CALL diagnostics_pos_vel(Vsq, Vx, Vy, Vz, Rx, Ry, Rz, Fx, Fy, Fz, ptemp, temp)
        ELSE
            ptemp = temp
        END IF

        IF ( k .GT. n+1 ) THEN                                              ! switch on the thermostat
            CALL thermostat_full_pid(temp, Vx, Vy, Vz)
        END IF

        IF ( MODULO(k, n+1) .EQ. 0 ) THEN                                   ! sample rate for images
            CALL plot_positions(k, temp, en, Rx, Ry, Rz, Vsq)
        END IF

        IF ( MODULO(k, n+1) .EQ. 0 ) THEN                                    ! sample rate for rdf images
            CALL plot_radial_dist(k, g)
        END IF

        !WRITE(6,'(4(a))',ADVANCE="NO") "\b","\b","\b","\b"                  ! backslash control character: backspace
        !CALL FLUSH(6)
        !WRITE(6,'(I3,"%")',ADVANCE="NO") INT( (k*100)/REAL(n) )             ! percentage complete indicator to stdout
        !CALL FLUSH(6)

    END DO

    CALL plot_energy()

    DEALLOCATE(Fx, Fy, Fz, Rx, Ry, Rz, Vx, Vy, Vz, Vsq, g)

    PRINT*, ""

    CALL SYSTEM_CLOCK(terminate)
    CALL CPU_TIME(finish)
    time_cpu = finish - start
    time_real = (terminate - initiate)/REAL(clock_rate)
    PRINT*, "CPU time =", time_cpu
    PRINT*, "Real time =", time_real
    PRINT*, ""

    !CALL EXECUTE_COMMAND_LINE("fortune")
    !PRINT*, ""

END PROGRAM main_brute

