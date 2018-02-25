!----------------------------------------------------------------------
!                           MOLECULAR DYNAMICS                        !
!----------------------------------------------------------------------

PROGRAM CrysMelt

    USE parameters, ONLY : n, npart_d, num_cell, num_b, set_temp
    USE lattice
    USE thermostat
    USE force
    USE plot
    USE integration
    USE diagnostics

    INTEGER :: k, AllocateStatus, initiate, terminate, clock_rate

    INTEGER, DIMENSION(:), ALLOCATABLE :: location

    REAL :: temp, en, time_cpu, pot_en, ptemp, start, finish, kin_en, time_real

    REAL, DIMENSION(:), ALLOCATABLE :: F, R, V
    REAL, DIMENSION(:), ALLOCATABLE :: Vsq
    REAL, DIMENSION(:), ALLOCATABLE :: g

    CALL CPU_TIME(start)
    CALL SYSTEM_CLOCK(initiate, clock_rate)

    CALL EXECUTE_COMMAND_LINE("rm data/*.png")
    CALL EXECUTE_COMMAND_LINE("rm data/Position/*.png")
    CALL EXECUTE_COMMAND_LINE("rm data/RadDistFunc/*.png")

    ALLOCATE(R(npart_d), F(npart_d), V(npart_d), Vsq(npart), &
        g(num_b), location(num_cell), STAT = AllocateStatus)

    IF (AllocateStatus /= 0) THEN
        STOP "*** Not enough memory ***"
    END IF


    CALL lattice_cell_simple_cubic(R, location)

    CALL lattice_velocities_random(V)

    CALL thermostat_full(V)                                                 ! rescale all vels using set_temp

    Vsq = ( (Vx*Vx) + (Vy*Vy) + (Vz*Vz) )/2.0                               ! normalized particle kinetic energies
    kin_en = SUM(Vsq)/REAL(npart)
    temp = (2*kin_en)/3.0
    ptemp = temp


    CALL force_create_cell_list(Rx, Ry, Rz)

    CALL force_cell_list(Rx, Ry, Rz, Fx, Fy, Fz, pot_en)

    CALL force_rdf(Rx, Ry, Rz, g)

    en = (pot_en/REAL(npart)) + kin_en                                      ! energy per particle

    k=0

    CALL plot_positions(k, temp, en, Rx, Ry, Rz, Vsq)

    CALL plot_radial_dist(k, g)

    CALL plot_energy_init()                                                 ! write titles of energy data file

    CALL thermostat_pid_init()                                              ! initialize PID variables

    !CALL diagnostics_20_43_init()

    DO k = 1, n                                                             ! multiple runs

        CALL integration_VelVerlet_pos(Rx, Ry, Rz, Vx, Vy, Vz, Fx, Fy, Fz)  ! positions from Velocity Verlet & PBC

        CALL integration_VelVerlet_vel(Vx, Vy, Vz, Fx, Fy, Fz)              ! first half-step of velocities from Velocity Verlet

        !IF ( (MODULO(k, 10) .EQ. 0) .OR. (k<10) ) THEN                      ! sample rate for images
        CALL force_create_cell_list(Rx, Ry, Rz)
        !END IF

        CALL force_cell_list(Rx, Ry, Rz, Fx, Fy, Fz, pot_en)

        CALL integration_VelVerlet_vel(Vx, Vy, Vz, Fx, Fy, Fz)              ! second half-step of velocities from Velocity Verlet

        Vsq = (Vx*Vx + Vy*Vy + Vz*Vz)/2.0                                   ! normalized particle kinetic energies
        kin_en = SUM(Vsq)/REAL(npart)
        temp = (2.0*kin_en)/3.0                                             ! instantaneous temperature
        pot_en = pot_en/REAL(npart)
        en = pot_en + kin_en                                                ! energy per particle

        CALL plot_energy_data(k, en, pot_en, kin_en)


        IF (temp .GT. (100*ptemp)) THEN                                     ! excessive temp fluctuation check
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

        !IF ( MODULO(k, n+1) .EQ. 0 ) THEN                                    ! sample rate for rdf images
        !    CALL force_rdf(Rx, Ry, Rz, g)
        !    CALL plot_radial_dist(k, g)
        !END IF

        WRITE(6,'(4(a))',ADVANCE="no") "\b","\b","\b","\b"                  ! backslash control character: backspace
        CALL FLUSH(6)
        WRITE(6,'(I3,"%")',ADVANCE="no") INT( (k*100)/REAL(n) )             ! percentage complete indicator to stdout
        CALL FLUSH(6)

    END DO

    CALL plot_energy()

    DEALLOCATE(Fx, Fy, Fz, Rx, Ry, Rz, Vx, Vy, Vz, Vsq)

    PRINT*, ""

    CALL SYSTEM_CLOCK(terminate)
    CALL CPU_TIME(finish)
    time_cpu = finish - start
    time_real = (terminate - initiate)/REAL(clock_rate)
    PRINT*, "CPU time =", time_cpu
    PRINT*, "Real time =", time_real
    PRINT*, ""

    !CALL EXECUTE_COMMAND_LINE("fortune")
    PRINT*, ""

END PROGRAM CrysMelt

