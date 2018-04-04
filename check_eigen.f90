PROGRAM check_eigen
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

    !USE la_precision, ONLY: WP => DP
    USE f95_lapack, ONLY: LA_SYEVR

    INTEGER :: IU
    REAL, DIMENSION(5,5) :: A
    REAL, DIMENSION(5) :: W

    A(:,1) = [1.0, 0.0, 0.0, 0.0, 0.0]
    A(:,2) = [2.0, 9.0, 0.0, 0.0, 0.0]
    A(:,3) = [-3.0, 0.0, -2.3, 0.0, 0.0]
    A(:,4) = [75.0, 0.0, 0.0, 91.0, 0.0]
    A(:,5) = [-0.01, 1.0, 21.0, -107.0, 0.05]

    !CALL LA_SYEVX(A,W,IU=1)
    CALL LA_SYEVR(A,W,IU=1)

    PRINT*, W


END PROGRAM check_eigen
