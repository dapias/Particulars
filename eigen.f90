!----------------------------------------------------------------------
!                            EIGENVALUES                              !
!----------------------------------------------------------------------

MODULE eigen

    !USE la_precision, ONLY: WP => DP
    USE f95_lapack, ONLY: LA_SYEV, LA_SYEVD, LA_SYEVX, LA_SYEVR

    IMPLICIT NONE

    ! The matrix 'A' is to be upper triangular; not the full symmetric matrix

CONTAINS


    SUBROUTINE eigen_syev(A, W, INFO)

        ! COMPUTE ALL EIGENVALUES AND, OPTIONALLY, ALL EIGENVECTORS OF A REAL SYMMETRIC MATRIX

        REAL, ALLOCATABLE, INTENT(IN) :: A(:,:)

        INTEGER, INTENT(OUT) :: INFO
        REAL, ALLOCATABLE, INTENT(OUT) :: W(:)

        CALL LA_SYEV(A, W, INFO)

    END SUBROUTINE eigen_syev


    SUBROUTINE eigen_syevd(A, W, INFO)

        ! USE A DIVIDE AND CONQUER ALGORITHM. IF EIGENVECTORS ARE DESIRED, CAN BE MUCH FASTER THAN syev FOR LARGE MATRICES, BUT USE MORE WORKSPACE.

        REAL, ALLOCATABLE, INTENT(IN) :: A(:,:)

        INTEGER, INTENT(OUT) :: INFO
        REAL, ALLOCATABLE, INTENT(OUT) :: W(:)

        CALL LA_SYEVD(A, W, INFO)

    END SUBROUTINE eigen_syevd


    SUBROUTINE eigen_syevx(A, W, INFO)

        ! COMPUTE SELECTED EIGENVALUES AND, OPTIONALLY, THE CORRESPONDING EIGENVECTORS OF A REAL SYMMETRIC HERMITIAN MATRIX

        REAL, ALLOCATABLE, INTENT(IN) :: A(:,:)

        INTEGER, INTENT(OUT) :: INFO
        REAL, ALLOCATABLE, INTENT(OUT) :: W(:)

        CALL LA_SYEVX(A, W, INFO)
        !CALL LA_SYEVX( A, W, IU=1)

    END SUBROUTINE eigen_syevx


    SUBROUTINE eigen_syevr(A, W, INFO)

        ! COMPUTE SELECTED EIGENVALUES AND, OPTIONALLY, THE CORRESPONDING EIGENVECTORS OF A REAL SYMMETRIC HERMITIAN MATRIX, USING A RELATIVELY ROBUST REPRESENTATION (RRR) ALGORITHM. IT IS USUALLY THE FASTEST ALGORITHM OF ALL AND USES THE LEAST WORKSPACE.

        REAL, ALLOCATABLE, INTENT(IN) :: A(:,:)

        INTEGER, INTENT(OUT) :: INFO
        REAL, ALLOCATABLE, INTENT(OUT) :: W(:)

        !CALL LA_SYEVR(A, W, INFO)
        CALL LA_SYEVR( A, W, IU=1, INFO)
        ! ^ Finds the minimum eigen-value

    END SUBROUTINE eigen_syevr


END MODULE eigen
