!----------------------------------------------------------------------
!                      RANDOM NUMBER INITIATORS                       !
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

MODULE random_init

    !!DIR$ IF DEFINED (__INTEL_COMPILER)
#ifdef __INTEL_COMPILER
    USE IFPORT
#endif
    !!DIR$ END IF

CONTAINS

    SUBROUTINE init_random_seed()

        USE iso_fortran_env, ONLY: int64

        INTEGER, DIMENSION(:), ALLOCATABLE :: seed
        INTEGER :: i, n, un, istat, dt(8), pid
        INTEGER(int64) :: t

        CALL random_seed(size = n)
        ALLOCATE(seed(n))
        ! First try if the OS provides a random number generator

        OPEN(NEWUNIT=un, FILE="/dev/urandom", ACCESS="stream", &
             FORM="unformatted", ACTION="read", STATUS="old", IOSTAT=istat)

        IF (istat == 0) THEN

            READ(un) seed
            CLOSE(un)

        ELSE

            ! Fallback to XOR:ing the current time and pid. The PID is
            ! useful in case one launches multiple instances of the same
            ! program in parallel.
            CALL SYSTEM_CLOCK(t)
            IF (t == 0) THEN
                CALL DATE_AND_TIME(values=dt)
                t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
                   + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
                   + dt(3) * 24_int64 * 60 * 60 * 1000 &
                   + dt(5) * 60 * 60 * 1000 &
                   + dt(6) * 60 * 1000 + dt(7) * 1000 &
                   + dt(8)
            END IF
            pid = GETPID()
            t = IEOR(t, INT(pid, KIND(t)))
            DO i = 1, n
                seed(i) = lcg(t)
            END DO

        END IF
        CALL random_seed(put=seed)

        CONTAINS
          ! This simple PRNG might not be good enough for real work, but is
          ! sufficient for seeding a better PRNG.
            FUNCTION lcg(s)
              INTEGER :: lcg
              INTEGER(int64) :: s
              IF (s == 0) THEN
                 s = 104729
              ELSE
                 s = MOD(s, 4294967296_int64)
              END IF
              s = MOD(s * 279470273_int64, 4294967291_int64)
              lcg = INT(MOD(s, INT(HUGE(0), int64)), KIND(0))
            end function lcg

    END SUBROUTINE init_random_seed


    SUBROUTINE init_random_seed_primitive()

        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed = clock + 37 * [(i - 1, i = 1, n)]
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
    END SUBROUTINE


    SUBROUTINE init_random_seed_fixed()

        INTEGER :: i, n
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        seed = 37 * [(i - 1, i = 1, n)]
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
    END SUBROUTINE

END MODULE random_init
