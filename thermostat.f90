!----------------------------------------------------------------------
!                              THERMOSTATS                            !
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

MODULE thermostat

    USE parameters, ONLY : npart, set_temp, dt, Kp, Ki, Kd

    REAL, PRIVATE :: proportional, integral, differential

CONTAINS



    PURE SUBROUTINE thermostat_full(Vx, Vy, Vz)

        ! A thermostat simulating an ambient bath/fluid that acts on the
        ! particles as and when called. It rescales all velocities equally.

        REAL, DIMENSION(:), INTENT(INOUT) :: Vx, Vy, Vz

        REAL :: sumVx2, sumVy2, sumVz2, sfx, sfy, sfz

        sumVx2 = 0 ; sumVy2 = 0 ; sumVz2 = 0
        sfx = 0 ; sfy = 0 ; sfz = 0

        sumVx2 = SUM(Vx*Vx)/REAL(npart)                             ! mean squared velocities
        sumVy2 = SUM(Vy*Vy)/REAL(npart)
        sumVz2 = SUM(Vz*Vz)/REAL(npart)

        sfx = set_temp/sumVx2                                       ! temperature scale factors
        sfy = set_temp/sumVy2
        sfz = set_temp/sumVz2

        Vx = Vx * SQRT(sfx)                                         ! scaling velocities
        Vy = Vy * SQRT(sfy)
        Vz = Vz * SQRT(sfz)

    END SUBROUTINE thermostat_full


    PURE SUBROUTINE thermostat_full_2D(Vx, Vy)

        ! A thermostat simulating an ambient bath/fluid that acts on the
        ! particles as and when called. It rescales all velocities equally.

        REAL, DIMENSION(:), INTENT(INOUT) :: Vx, Vy

        REAL :: sumVx2, sumVy2, sfx, sfy

        sumVx2 = 0 ; sumVy2 = 0
        sfx = 0 ; sfy = 0

        sumVx2 = SUM(Vx*Vx)/REAL(npart)                             ! mean squared velocities
        sumVy2 = SUM(Vy*Vy)/REAL(npart)

        sfx = set_temp/sumVx2                                       ! temperature scale factors
        sfy = set_temp/sumVy2

        Vx = Vx * SQRT(sfx)                                         ! scaling velocities
        Vy = Vy * SQRT(sfy)

    END SUBROUTINE thermostat_full_2D




    SUBROUTINE thermostat_pid_init()

        proportional=0 ; integral=0 ; differential=0

    END SUBROUTINE thermostat_pid_init


    SUBROUTINE thermostat_full_pid(current_temp, Vx, Vy, Vz)

        ! A thermostat_full with a PID controlling approach to the set temperature.
        ! Set gains in MODULE parameters.

        REAL, INTENT(IN) :: current_temp
        REAL, DIMENSION(:), INTENT(INOUT) :: Vx, Vy, Vz

        REAL :: temp, error, sumVx2, sumVy2, sumVz2, sfx, sfy, sfz

        error = set_temp - current_temp
        integral = integral + (error*dt)
        differential = (error - proportional)/dt

        proportional = error

        temp = current_temp + (Kp*proportional) + (Ki*integral) + (Kd*differential)


        sumVx2 = 0 ; sumVy2 = 0 ; sumVz2 = 0
        sfx = 0 ; sfy = 0 ; sfz = 0

        sumVx2 = SUM(Vx*Vx)/REAL(npart)                             ! mean squared velocities
        sumVy2 = SUM(Vy*Vy)/REAL(npart)
        sumVz2 = SUM(Vz*Vz)/REAL(npart)

        sfx = temp/sumVx2                                           ! temperature scale factors
        sfy = temp/sumVy2
        sfz = temp/sumVz2

        Vx = Vx * SQRT(sfx)                                         ! scaling velocities
        Vy = Vy * SQRT(sfy)
        Vz = Vz * SQRT(sfz)

    END SUBROUTINE thermostat_full_pid


    SUBROUTINE thermostat_full_pid_2D(current_temp, Vx, Vy)

        ! A thermostat_full with a PID controlling approach to the set temperature.
        ! Set gains in MODULE parameters.

        REAL, INTENT(IN) :: current_temp
        REAL, DIMENSION(:), INTENT(INOUT) :: Vx, Vy

        REAL :: temp, error, sumVx2, sumVy2, sfx, sfy

        error = set_temp - current_temp
        integral = integral + (error*dt)
        differential = (error - proportional)/dt

        proportional = error

        temp = current_temp + (Kp*proportional) + (Ki*integral) + (Kd*differential)


        sumVx2 = 0 ; sumVy2 = 0
        sfx = 0 ; sfy = 0

        sumVx2 = SUM(Vx*Vx)/REAL(npart)                             ! mean squared velocities
        sumVy2 = SUM(Vy*Vy)/REAL(npart)

        sfx = temp/sumVx2                                           ! temperature scale factors
        sfy = temp/sumVy2

        Vx = Vx * SQRT(sfx)                                         ! scaling velocities
        Vy = Vy * SQRT(sfy)

    END SUBROUTINE thermostat_full_pid_2D




    SUBROUTINE thermostat_wall()

        ! A particle passing through a wall results in a rescaling of its
        ! velocity commensurate with a set temperature.


    END SUBROUTINE thermostat_wall



END MODULE thermostat
