PROGRAM check_eigen

    !USE la_precision, ONLY: WP => DP
    USE f95_lapack, ONLY: LA_SYEVX

    IMPLICIT NONE

    INTEGER :: IU
    REAL, DIMENSION(5,5) :: A
    REAL, DIMENSION(5) :: W

    A(1,:) = [1.0, 0.0, 0.0, 0.0, 0.0]
    A(2,:) = [0.0, 9.0, 0.0, 0.0, 0.0]
    A(3,:) = [0.0, 0.0, -2.3, 0.0, 0.0]
    A(4,:) = [0.0, 0.0, 0.0, 91.0, 0.0]
    A(5,:) = [0.0, 0.0, 0.0, 0.0, 0.05]

    CALL LA_SYEVX(A,W,IU=1)

    PRINT*, W


END PROGRAM check_eigen
