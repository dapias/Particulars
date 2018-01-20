!-------------------------------------------------------------------------------!
!                                    PARAMETERS                                 !
!-------------------------------------------------------------------------------!

MODULE parameters

    IMPLICIT NONE

    INTEGER, PARAMETER :: n=1000        ! Number of iterations
    INTEGER, PARAMETER :: npart=1000    ! Number of particles
    INTEGER, PARAMETER :: num_b=50      ! Number of bins for g(r)
    INTEGER, PARAMETER :: dimnsn=3      ! Dimensions of the box

    REAL, PARAMETER :: real_zero=REAL(0)! A real zero
    REAL, PARAMETER :: Eps=1            ! Potential parameter
    REAL, PARAMETER :: Eps00=1          ! Potential parameter for Kob-Andersen
    REAL, PARAMETER :: Eps01=1.5        ! Potential parameter for Kob-Andersen
    REAL, PARAMETER :: Eps11=0.5        ! Potential parameter for Kob-Andersen
    REAL, PARAMETER :: Sig=1            ! Potential parameter
    REAL, PARAMETER :: Sig00=1          ! Potential parameter for Kob-Andersen
    REAL, PARAMETER :: Sig01=0.8        ! Potential parameter for Kob-Andersen
    REAL, PARAMETER :: Sig11=0.88       ! Potential parameter for Kob-Andersen
    REAL, PARAMETER :: set_temp=0.55    ! Set temperature
    REAL, PARAMETER :: Kp=0.1           ! Proportion gain for the thermostat PID
    REAL, PARAMETER :: Ki=0.05          ! Integration gain for the thermostat PID
    REAL, PARAMETER :: Kd=0.0005        ! Difference gain for the thermostat PID
    REAL, PARAMETER :: dt=0.0005         ! Time step
    REAL, PARAMETER :: Rc=2.8           ! Cut-off radius for force
    REAL, PARAMETER, PRIVATE :: pi=4.0*ATAN(1.0) ! Pi

    REAL, PARAMETER :: rho=0.5                                      ! Particle density
    REAL, PARAMETER :: box=(npart/rho)**(1/REAL(dimnsn))            ! Box side length
    !REAL, PARAMETER :: box=12                                       ! Box side length
    !REAL, PARAMETER :: rho = npart/(box**dimnsn)                    ! particle density

    INTEGER, PARAMETER :: npart_d = dimnsn * npart                  ! required array length
    INTEGER, PARAMETER :: npart_edge = CEILING(npart**(1/REAL(dimnsn))) ! particle number along edge
    INTEGER, PARAMETER :: npart_face = npart_edge**2                ! particle number on a face
    INTEGER, PARAMETER :: num_cell_edge = INT(box/Rc)               ! number of cells along edge
    INTEGER, PARAMETER :: num_cell_face = num_cell_edge**2          ! number of cells along face
    INTEGER, PARAMETER :: num_cell = num_cell_edge ** dimnsn        ! total number of cells
    !INTEGER, PARAMETER :: npart_cell_edge = npart_edge/ ! particles per cell, along edge
    INTEGER, PARAMETER :: npart_cell = CEILING(npart/REAL(num_cell))! particles per cell
    INTEGER, PARAMETER :: npart_d_cell = dimnsn * npart_cell        ! space on array per cell

    REAL, PARAMETER :: SigSq = Sig**2
    REAL, PARAMETER :: SigSq00 = Sig00**2
    REAL, PARAMETER :: SigSq01 = Sig01**2
    REAL, PARAMETER :: SigSq11 = Sig11**2
    REAL, PARAMETER :: box2 = box/2.0
    REAL, PARAMETER :: space = box / REAL(npart_edge)               ! lattice spacing
    REAL, PARAMETER :: RcSq = Rc**2
    REAL, PARAMETER, PRIVATE :: ratio = (SigSq/Rcsq)**3
    REAL, PARAMETER :: pot_en_cut = 4*Eps*ratio*(ratio-1)           ! excess pot en due to Rc
    REAL, PARAMETER :: size_b = box/(2.0*num_b)                     ! bin size for g(r)
    REAL, PARAMETER :: norm_b = ((dimnsn+1)/3.0) * pi * rho * (size_b**dimnsn)
    ! normlization factor for g(r) with a coefficient; works for 2D & 3D
    REAL, PARAMETER :: cell = box / REAL(num_cell_edge)             ! edge length of a cell box

END MODULE parameters
