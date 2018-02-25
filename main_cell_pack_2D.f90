!----------------------------------------------------------------------
!                           MOLECULAR DYNAMICS                        !
!----------------------------------------------------------------------

PROGRAM CrysMelt

    USE parameters, ONLY : n, npart, num_b, set_temp
    USE lattice
    USE thermostat
    USE force
    USE plot
    USE integration
    USE diagnostics

    INTEGER :: k, AllocateStatus, initiate, terminate, clock_rate

    REAL :: temp, en, time_cpu, pot_en, ptemp, start, finish, kin_en, time_real

    REAL, DIMENSION(:), ALLOCATABLE :: Fx, Fy, Rx, Ry, Vx, Vy, Vsq
    REAL, DIMENSION(:), ALLOCATABLE :: g

    CALL diagnostics_NCE_check()

    CALL CPU_TIME(start)
    CALL SYSTEM_CLOCK(initiate, clock_rate)

    CALL EXECUTE_COMMAND_LINE("rm data/*.png")
    CALL EXECUTE_COMMAND_LINE("rm data/Position/*.png")
    CALL EXECUTE_COMMAND_LINE("rm data/RadDistFunc/*.png")

    ALLOCATE(Rx(npart), Ry(npart), Fx(npart), Fy(npart), &
    Vx(npart), Vy(npart), Vsq(npart), g(num_b), STAT = AllocateStatus)

    IF (AllocateStatus /= 0) THEN
        STOP "*** Not enough memory ***"
    END IF


    CALL lattice_square(Rx, Ry)

    CALL lattice_velocities_random_2D(Vx, Vy)

    CALL thermostat_full_2D(Vx, Vy)                               ! rescale all vels using set_temp

    Vsq = ( (Vx*Vx) + (Vy*Vy) )/2.0                               ! normalized particle kinetic energies
    kin_en = SUM(Vsq)/REAL(npart)
    temp = kin_en
    ptemp = temp


    CALL force_create_cell_list_2D(Rx, Ry)

    CALL force_cell_list_2D(Rx, Ry, Fx, Fy, pot_en)

    CALL force_rdf_2D(Rx, Ry, g)

    en = (pot_en/REAL(npart)) + kin_en                                      ! energy per particle

    k=0

    CALL plot_positions_2D(k, temp, en, Rx, Ry, Vsq)

    CALL plot_radial_dist(k, g)

    CALL plot_energy_init()                                                 ! write titles of energy data file

    CALL thermostat_pid_init()                                              ! initialize PID variables


    DO k = 1, n                                                             ! multiple runs

        CALL integration_VelVerlet_pos_2D(Rx, Ry, Vx, Vy, Fx, Fy)           ! positions from Velocity Verlet & PBC

        CALL integration_VelVerlet_vel_2D(Vx, Vy, Fx, Fy)                   ! first half-step of velocities from Velocity Verlet

        !IF ( (MODULO(k, 10) .EQ. 0) .OR. (k<10) ) THEN                      ! sample rate for images
        CALL force_create_cell_list_2D(Rx, Ry)
        !END IF

        CALL force_cell_list_2D(Rx, Ry, Fx, Fy, pot_en)

        CALL integration_VelVerlet_vel_2D(Vx, Vy, Fx, Fy)                      ! second half-step of velocities from Velocity Verlet

        Vsq = (Vx*Vx + Vy*Vy)/2.0                                           ! normalized particle kinetic energies
        kin_en = SUM(Vsq)/REAL(npart)
        temp = kin_en                                                       ! instantaneous temperature
        pot_en = pot_en/REAL(npart)
        en = pot_en + kin_en                                                ! energy per particle

        CALL plot_energy_data(k, en, pot_en, kin_en)


        IF (temp .GT. (100*ptemp)) THEN                                     ! excessive temp fluctuation check
            CALL plot_energy()
            PRINT*, ""
            CALL EXECUTE_COMMAND_LINE("fortune -o")
            PRINT*, ""
            CALL diagnostics_pos_vel_2D(Vsq, Vx, Vy, Rx, Ry, Fx, Fy, ptemp, temp)
        ELSE
            ptemp = temp
        END IF

        IF ( k .GT. n+1 ) THEN                                              ! switch on the thermostat
            CALL thermostat_full_pid_2D(temp, Vx, Vy)
        END IF

        IF ( MODULO(k, n+1) .EQ. 0 ) THEN                                   ! sample rate for images
            CALL plot_positions_2D(k, temp, en, Rx, Ry, Vsq)
        END IF

        !IF ( MODULO(k, n+1) .EQ. 0 ) THEN                                    ! sample rate for rdf images
        !    CALL force_rdf_2D(Rx, Ry, g)
        !    CALL plot_radial_dist(k, g)
        !END IF

        WRITE(6,'(4(a))',ADVANCE="no") "\b","\b","\b","\b"                  ! backslash control character: backspace
        CALL FLUSH(6)
        WRITE(6,'(I3,"%")',ADVANCE="no") INT( (k*100)/REAL(n) )             ! percentage complete indicator to stdout
        CALL FLUSH(6)

    END DO

    CALL plot_energy()

    DEALLOCATE(Fx, Fy, Rx, Ry, Vx, Vy, Vsq)

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

