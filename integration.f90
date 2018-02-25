!----------------------------------------------------------------------
!                   INTEGRATING EQUATIONS OF MOTION                   !
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

MODULE integration

    USE parameters, ONLY : npart, box, dt

    REAL, DIMENSION(npart), PRIVATE, PARAMETER :: box_array = box

CONTAINS



    SUBROUTINE integration_VelVerlet_pos(Rx, Ry, Rz, Vx, Vy, Vz, Fx, Fy, Fz)

        !  Positions from the Velocity Verlet algorithm with PBCs

        REAL, DIMENSION(:), INTENT(IN) :: Vx, Vy, Vz, Fx, Fy, Fz
        REAL, DIMENSION(:), INTENT(INOUT) :: Rx, Ry, Rz

        Rx = MODULO(Rx + (Vx*dt) + (0.5*Fx*(dt**2)), box_array)
        Ry = MODULO(Ry + (Vy*dt) + (0.5*Fy*(dt**2)), box_array)
        Rz = MODULO(Rz + (Vz*dt) + (0.5*Fz*(dt**2)), box_array)

    END SUBROUTINE integration_VelVerlet_pos

    SUBROUTINE integration_VelVerlet_pos_2D(Rx, Ry, Vx, Vy, Fx, Fy)

        !  Positions from the Velocity Verlet algorithm with PBCs

        REAL, DIMENSION(:), INTENT(IN) :: Vx, Vy, Fx, Fy
        REAL, DIMENSION(:), INTENT(INOUT) :: Rx, Ry

        Rx = MODULO(Rx + (Vx*dt) + (0.5*Fx*(dt**2)), box_array)
        Ry = MODULO(Ry + (Vy*dt) + (0.5*Fy*(dt**2)), box_array)

    END SUBROUTINE integration_VelVerlet_pos_2D



    SUBROUTINE integration_VelVerlet_vel(Vx, Vy, Vz, Fx, Fy, Fz)

        !  A velocity half-step from the Velocity Verlet algorithm

        REAL, DIMENSION(:), INTENT(IN) :: Fx, Fy, Fz
        REAL, DIMENSION(:), INTENT(INOUT) :: Vx, Vy, Vz

        Vx = Vx + (0.5*Fx*dt)
        Vy = Vy + (0.5*Fy*dt)
        Vz = Vz + (0.5*Fz*dt)

    END SUBROUTINE integration_VelVerlet_vel


    SUBROUTINE integration_VelVerlet_vel_2D(Vx, Vy, Fx, Fy)

        !  A velocity half-step from the Velocity Verlet algorithm

        REAL, DIMENSION(:), INTENT(IN) :: Fx, Fy
        REAL, DIMENSION(:), INTENT(INOUT) :: Vx, Vy

        Vx = Vx + (0.5*Fx*dt)
        Vy = Vy + (0.5*Fy*dt)

    END SUBROUTINE integration_VelVerlet_vel_2D



END MODULE integration
