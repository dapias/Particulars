!-------------------------------------------------------------------------------!
!                                FORCE CALCULATION                              !
!-------------------------------------------------------------------------------!

MODULE force

    USE parameters, ONLY : npart, box, Rc, RcSq, SigSq, SigSq00, SigSq01, SigSq11, Eps, Eps00, &
        Eps01, Eps11, Sig, Sig00, Sig01, Sig11, pot_en_cut, box2, size_b, &
        num_cell_edge, num_cell_face, num_cell, cell, real_zero

    IMPLICIT NONE

    LOGICAL, DIMENSION(npart, num_cell), PRIVATE :: mask
    REAL, PARAMETER, PRIVATE :: cell2 = 2*cell
    REAL, PARAMETER, PRIVATE :: cell2m = (-1)*cell2

CONTAINS



    SUBROUTINE force_brute(Rx, Ry, Rz, Fx, Fy, Fz, pot_en, g)

        ! Regular force calculation, checking all partices, and calculating forces if within a cut-off radius

        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Rz
        INTEGER :: i, j, ri
        REAL :: Rxi, Ryi, Rzi, dRx, dRy, dRz, dR, dRsq, RF2, RF6, F
        REAL, INTENT(OUT) :: pot_en
        REAL, DIMENSION(:), INTENT(OUT) :: Fx, Fy, Fz, g

        F = 0
        Fx = 0 ; Fy = 0 ; Fz = 0
        pot_en = 0

        !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Rxi,Ryi,Rzi,dRx,dRy,dRz,dRsq,dR,RF2,RF6,F,ri)
        !$OMP DO REDUCTION(+:Fx,Fy,Fz,pot_en,g)
        DO i = 1, npart-1                                                   ! pair-wise force calculation

            Rxi = Rx(i) ; Ryi = Ry(i) ; Rzi = Rz(i)

            DO j = i+1, npart

                dRx = Rxi - Rx(j)                                           ! x-distance between two particles
                dRx = dRx - (box*NINT(dRx/box))                             ! choosing nearest distance, due to periodic boudaries
                dRy = Ryi - Ry(j)
                dRy = dRy - (box*NINT(dRy/box))
                dRz = Rzi - Rz(j)
                dRz = dRz - (box*NINT(dRz/box))

                dRsq = (dRx)**2 + (dRy)**2 + (dRz)**2                       ! distance between particles i and j
                dR = SQRT(dRsq)

                IF (dR .lt. Rc) THEN
                    RF2 = SigSq/dRsq
                    RF6 = RF2 ** 3
                    F = 48 * (Eps/Sig) * RF2 * RF6 * (RF6 - 0.5)            ! force

                    Fx(i) = Fx(i) + (F*dRx)                                 ! pairwise force addition
                    Fy(i) = Fy(i) + (F*dRy)
                    Fz(i) = Fz(i) + (F*dRz)
                    Fx(j) = Fx(j) - (F*dRx)
                    Fy(j) = Fy(j) - (F*dRy)
                    Fz(j) = Fz(j) - (F*dRz)

                    pot_en = pot_en + (4 * RF6 * (RF6 - 1)) - pot_en_cut    ! potential energy
                END IF

                IF (dR .lt. box2) THEN
                    ri = CEILING(dR/size_b)                                 ! populating g(r)
                    g(ri) = g(ri) + 2
                END IF

            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE force_brute



    SUBROUTINE force_neighbour_list()

        ! Implementing Bekker's form of the Verlet list algorithm.

    END SUBROUTINE force_neighbour_list



    SUBROUTINE force_create_cell_list(Rx, Ry, Rz)

        ! Making a full cell list

        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Rz

        INTEGER :: i
        INTEGER, DIMENSION(npart) :: c

        c=0

        mask = .FALSE.

        c = cell_id(  INT(Rz/cell), INT(Ry/cell), INT(Rx/cell)  )

        DO i = 1, npart
            mask( i, c(i) ) = .TRUE.
        END DO

    END SUBROUTINE force_create_cell_list


    SUBROUTINE force_create_cell_list_2D(Rx, Ry)

        ! Making a full cell list

        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry

        INTEGER :: i
        INTEGER, DIMENSION(npart) :: c

        c=0

        mask = .FALSE.

        c = cell_id_2D(  INT(Ry/cell), INT(Rx/cell)  )

        DO CONCURRENT (i = 1 : npart)
            mask( i, c(i) ) = .TRUE.
        END DO

    END SUBROUTINE force_create_cell_list_2D


    SUBROUTINE force_cell_list(Rx, Ry, Rz, Fx, Fy, Fz, pot_en)

        ! Cell list algorithm

        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Rz

        INTEGER :: i, j, m, p, q, s, l, o, r
        REAL :: Rxp, Ryp, Rzp, dRx, dRy, dRz, dRsq, RF2, RF6, F
        LOGICAL :: border_i_0, border_i_n, border_j_0, border_j_n, border_m
        LOGICAL, DIMENSION(npart) :: mask_l_force, mask_s_force
        REAL, DIMENSION(:), ALLOCATABLE :: cell_Rx, cell_Ry, cell_Rz, cell_Fx, cell_Fy, cell_Fz

        REAL, INTENT(OUT) :: pot_en
        REAL, DIMENSION(:), INTENT(OUT) :: Fx, Fy, Fz

        mask_l_force = .FALSE.
        mask_s_force = .FALSE.
        Fx = 0 ; Fy = 0 ; Fz = 0
        pot_en = 0


        !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(mask_force,s,cell_Rx,cell_Ry,cell_Rz,dRx,dRy,dRz,dRsq,dR,RF2,RF6,F,ri)
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(mask,Rx,Ry,Rz,Fx,Fy,Fz,pot_en)
        !$OMP DO REDUCTION(+:Fx,Fy,Fz,pot_en)
        DO m = 0, num_cell_edge-1
            border_m = ( m == (num_cell_edge - 1) )                                     ! a check for when near boundary
            DO j = 0, num_cell_edge-1
                border_j_0 =  (j == 0)
                border_j_n = ( j == (num_cell_edge - 1) )
                DO i = 0, num_cell_edge-1
                    border_i_0 = (i == 0)
                    border_i_n = ( i == (num_cell_edge - 1) )

                    F = 0

                    mask_l_force = mask(:,  cell_id(   m,   j,   i )  )                 ! mask of the cell

                    l = COUNT(mask_l_force)                                             ! no. of particles in the cell


                    mask_s_force = mask(:,  cell_id(   m,   j, i+1 )  ) &               ! mask of half the surrounding cells
                              .OR. mask(:,  cell_id(   m, j+1, i-1 )  ) &
                              .OR. mask(:,  cell_id(   m, j+1,   i )  ) &
                              .OR. mask(:,  cell_id(   m, j+1, i+1 )  ) &
                              .OR. mask(:,  cell_id( m+1, j-1, i-1 )  ) &
                              .OR. mask(:,  cell_id( m+1, j-1,   i )  ) &
                              .OR. mask(:,  cell_id( m+1, j-1, i+1 )  ) &
                              .OR. mask(:,  cell_id( m+1,   j, i-1 )  ) &
                              .OR. mask(:,  cell_id( m+1,   j,   i )  ) &
                              .OR. mask(:,  cell_id( m+1,   j, i+1 )  ) &
                              .OR. mask(:,  cell_id( m+1, j+1, i-1 )  ) &
                              .OR. mask(:,  cell_id( m+1, j+1,   i )  ) &
                              .OR. mask(:,  cell_id( m+1, j+1, i+1 )  )

                    s = COUNT(mask_s_force)                                             ! no. of part. in half the surr. cells

                    r = l+s

                    ALLOCATE( cell_Rx(r), cell_Ry(r), cell_Rz(r) )

                    cell_Rx(1:l) = PACK(Rx, mask_l_force)                               ! packing the cell
                    cell_Ry(1:l) = PACK(Ry, mask_l_force)
                    cell_Rz(1:l) = PACK(Rz, mask_l_force)

                    cell_Rx(l+1:r) = PACK(Rx, mask_s_force)                             ! packing half the surrounding cells
                    cell_Ry(l+1:r) = PACK(Ry, mask_s_force)
                    cell_Rz(l+1:r) = PACK(Rz, mask_s_force)

                    ALLOCATE( cell_Fx(r), cell_Fy(r), cell_Fz(r) )

                    cell_Fx = [(real_zero, o=1, r)]
                    cell_Fy = [(real_zero, o=1, r)]
                    cell_Fz = [(real_zero, o=1, r)]


                    DO p = 1, l                                                         ! iterating over particles in the cell

                        Rxp = cell_Rx(p) ; Ryp = cell_Ry(p) ; Rzp = cell_Rz(p)

                        DO q = p+1, r                                                   ! iterate over all particles considered

                            dRx = Rxp - cell_Rx(q)                                      ! x-distance between two particles
                            dRy = Ryp - cell_Ry(q)
                            dRz = Rzp - cell_Rz(q)

                            IF ( border_i_n ) THEN
                                IF ( dRx .GT. cell2 ) THEN
                                    dRx = dRx - box                                     ! choosing nearest distance, due to PBCs
                                END IF
                            ELSEIF ( border_i_0 ) THEN
                                IF ( dRx .LT. cell2m ) THEN
                                    dRx = dRx + box
                                END IF
                            END IF

                            IF ( border_j_n ) THEN
                                IF ( dRy .GT. cell2 ) THEN
                                    dRy = dRy - box
                                END IF
                            ELSEIF ( border_j_0 ) THEN
                                IF ( dRy .LT. cell2m ) THEN
                                    dRy = dRy + box
                                END IF
                            END IF

                            IF ( border_m ) THEN
                                IF ( dRz .GT. cell2 ) THEN
                                    dRz = dRz - box
                                END IF
                            END IF

                            dRsq = (dRx)**2 + (dRy)**2 + (dRz)**2                       ! distance between particles p and q
                            IF (dRsq .lt. RcSq) THEN

                                RF2 = SigSq/dRsq
                                RF6 = RF2 ** 3
                                F = 48 * (Eps/Sig) * RF2 * RF6 * (RF6 - 0.5)            ! force

                                cell_Fx(p) = cell_Fx(p) + (F*dRx)                       ! pairwise force addition
                                cell_Fy(p) = cell_Fy(p) + (F*dRy)
                                cell_Fz(p) = cell_Fz(p) + (F*dRz)
                                cell_Fx(q) = cell_Fx(q) - (F*dRx)
                                cell_Fy(q) = cell_Fy(q) - (F*dRy)
                                cell_Fz(q) = cell_Fz(q) - (F*dRz)

                                pot_en = pot_en + (4 * RF6 * (RF6 - 1)) - pot_en_cut    ! potential energy
                            END IF

                        END DO

                    END DO

                    Fx = Fx + UNPACK(cell_Fx(1:l), mask_l_force, real_zero)                 ! unpack into the force array
                    Fy = Fy + UNPACK(cell_Fy(1:l), mask_l_force, real_zero)
                    Fz = Fz + UNPACK(cell_Fz(1:l), mask_l_force, real_zero)

                    Fx = Fx + UNPACK(cell_Fx(l+1:r), mask_s_force, real_zero)
                    Fy = Fy + UNPACK(cell_Fy(l+1:r), mask_s_force, real_zero)
                    Fz = Fz + UNPACK(cell_Fz(l+1:r), mask_s_force, real_zero)

                    DEALLOCATE( cell_Rx, cell_Ry, cell_Rz, cell_Fx, cell_Fy, cell_Fz )

                END DO
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL


    END SUBROUTINE force_cell_list



    SUBROUTINE force_cell_list_KobAnd(Rx, Ry, Rz, species, Fx, Fy, Fz, pot_en)

        ! Cell list algorithm for a Kob-Andersen binary mixture

        LOGICAL, DIMENSION(:), INTENT(IN) :: species
        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Rz

        INTEGER :: i, j, m, p, q, s, l, o, r
        REAL :: Rxp, Ryp, Rzp, dRx, dRy, dRz, dRsq, RF2, RF6, F, Eps_b, SigSq_b, Sig_b
        LOGICAL :: border_i_0, border_i_n, border_j_0, border_j_n, border_m
        LOGICAL :: species_p, species_q
        LOGICAL, DIMENSION(npart) :: mask_l_force, mask_s_force
        LOGICAL, DIMENSION(:), ALLOCATABLE :: cell_species, cell_speciesl, cell_speciess
        REAL, DIMENSION(:), ALLOCATABLE :: cell_Rx, cell_Ry, cell_Rz, cell_Fx, cell_Fy, cell_Fz
        REAL, DIMENSION(:), ALLOCATABLE :: cell_Rxl, cell_Ryl, cell_Rzl, cell_Fxl, cell_Fyl, cell_Fzl
        REAL, DIMENSION(:), ALLOCATABLE :: cell_Rxs, cell_Rys, cell_Rzs, cell_Fxs, cell_Fys, cell_Fzs

        REAL, INTENT(OUT) :: pot_en
        REAL, DIMENSION(:), INTENT(OUT) :: Fx, Fy, Fz

        mask_l_force = .FALSE.
        mask_s_force = .FALSE.
        Fx = 0 ; Fy = 0 ; Fz = 0
        pot_en = 0


        !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(mask_force,s,cell_Rx,cell_Ry,cell_Rz,dRx,dRy,dRz,dRsq,dR,RF2,RF6,F,ri)
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(mask,Rx,Ry,Rz,species,Fx,Fy,Fz,pot_en)
        !$OMP DO REDUCTION(+:Fx,Fy,Fz,pot_en)
        DO m = 0, num_cell_edge-1
            border_m = ( m == (num_cell_edge - 1) )                                     ! a check for when near boundary
            DO j = 0, num_cell_edge-1
                border_j_0 =  (j == 0)
                border_j_n = (j == (num_cell_edge - 1))
                DO i = 0, num_cell_edge-1
                    border_i_0 = (i == 0)
                    border_i_n = (i == (num_cell_edge - 1))

                    F = 0

                    mask_l_force = mask(:,  cell_id(   m,   j,   i )  )                 ! mask of the cell

                    !ALLOCATE( cell_Rxl(l), cell_Ryl(l), cell_Rzl(l), cell_Fxl(l), cell_Fyl(l), cell_Fzl(l) )

                    cell_Rxl = PACK(Rx, mask_l_force)                                   ! packing the cell
                    cell_Ryl = PACK(Ry, mask_l_force)
                    cell_Rzl = PACK(Rz, mask_l_force)
                    cell_speciesl = PACK(species, mask_l_force)

                    l = SIZE(cell_Rxl)                                                  ! no. of particles in the cell

                    cell_Fxl = [(real_zero, o=1, l)]
                    cell_Fyl = [(real_zero, o=1, l)]
                    cell_Fzl = [(real_zero, o=1, l)]


                    mask_s_force = mask(:,  cell_id(   m,   j, i+1 )  ) &               ! mask of half the surrounding cells
                              .OR. mask(:,  cell_id(   m, j+1, i-1 )  ) &
                              .OR. mask(:,  cell_id(   m, j+1,   i )  ) &
                              .OR. mask(:,  cell_id(   m, j+1, i+1 )  ) &
                              .OR. mask(:,  cell_id( m+1, j-1, i-1 )  ) &
                              .OR. mask(:,  cell_id( m+1, j-1,   i )  ) &
                              .OR. mask(:,  cell_id( m+1, j-1, i+1 )  ) &
                              .OR. mask(:,  cell_id( m+1,   j, i-1 )  ) &
                              .OR. mask(:,  cell_id( m+1,   j,   i )  ) &
                              .OR. mask(:,  cell_id( m+1,   j, i+1 )  ) &
                              .OR. mask(:,  cell_id( m+1, j+1, i-1 )  ) &
                              .OR. mask(:,  cell_id( m+1, j+1,   i )  ) &
                              .OR. mask(:,  cell_id( m+1, j+1, i+1 )  )


                    !ALLOCATE( cell_Rx(s), cell_Ry(s), cell_Rz(s), cell_Fx(s), cell_Fy(s), cell_Fz(s) )

                    cell_Rxs = PACK(Rx, mask_s_force)                                   ! packing half the surrounding cells
                    cell_Rys = PACK(Ry, mask_s_force)
                    cell_Rzs = PACK(Rz, mask_s_force)
                    cell_speciess = PACK(species, mask_s_force)

                    s = SIZE(cell_Rxs)                                                  ! no. of part. in half the surr. cells

                    cell_Fxs = [(real_zero, o=1, s)]
                    cell_Fys = [(real_zero, o=1, s)]
                    cell_Fzs = [(real_zero, o=1, s)]

                    r = l+s

                    !ALLOCATE( cell_Rx(r), cell_Ry(r), cell_Rz(r), cell_Fx(r), cell_Fy(r), cell_Fz(r) )

                    cell_Rx = [cell_Rxl, cell_Rxs]                                      ! concatenate l & s arrays
                    cell_Ry = [cell_Ryl, cell_Rys]
                    cell_Rz = [cell_Rzl, cell_Rzs]
                    cell_species = [cell_speciesl, cell_speciess]

                    cell_Fx = [cell_Fxl, cell_Fxs]
                    cell_Fy = [cell_Fyl, cell_Fys]
                    cell_Fz = [cell_Fzl, cell_Fzs]


                    DO p = 1, l                                                         ! iterating over particles in the cell

                        Rxp = cell_Rx(p) ; Ryp = cell_Ry(p) ; Rzp = cell_Rz(p)
                        species_p = cell_species(p)

                        DO q = p+1, r                                                   ! iterate over all particles considered

                            species_q = cell_species(q)

                            dRx = Rxp - cell_Rx(q)                                      ! x-distance between two particles
                            dRy = Ryp - cell_Ry(q)
                            dRz = Rzp - cell_Rz(q)

                            IF ( border_i_n ) THEN
                                IF ( dRx .GT. cell2 ) THEN
                                    dRx = dRx - box                                     ! choosing nearest distance, due to PBCs
                                END IF
                            ELSEIF ( border_i_0 ) THEN
                                IF ( dRx .LT. cell2m ) THEN
                                    dRx = dRx + box
                                END IF
                            END IF

                            IF ( border_j_n ) THEN
                                IF ( dRy .GT. cell2 ) THEN
                                    dRy = dRy - box
                                END IF
                            ELSEIF ( border_j_0 ) THEN
                                IF ( dRy .LT. cell2m ) THEN
                                    dRy = dRy + box
                                END IF
                            END IF

                            IF ( border_m ) THEN
                                IF ( dRz .GT. cell2 ) THEN
                                    dRz = dRz - box
                                END IF
                            END IF

                            dRsq = (dRx)**2 + (dRy)**2 + (dRz)**2                       ! distance between particles p and q
                            IF (dRsq .lt. RcSq) THEN

                                IF ( .NOT. (species_p .OR. species_q) ) THEN
                                    Eps_b = Eps00
                                    SigSq_b = SigSq00
                                    Sig_b = Sig00
                                ELSEIF ( species_p .NEQV. species_q ) THEN
                                    Eps_b = Eps01
                                    SigSq_b = SigSq01
                                    Sig_b = Sig01
                                ELSE
                                    Eps_b = Eps11
                                    SigSq_b = SigSq11
                                    Sig_b = Sig11
                                ENDIF


                                RF2 = SigSq_b/dRsq
                                RF6 = RF2 ** 3
                                F = 48 * (Eps_b/Sig_b) * RF2 * RF6 * (RF6 - 0.5)            ! force

                                cell_Fx(p) = cell_Fx(p) + (F*dRx)                       ! pairwise force addition
                                cell_Fy(p) = cell_Fy(p) + (F*dRy)
                                cell_Fz(p) = cell_Fz(p) + (F*dRz)
                                cell_Fx(q) = cell_Fx(q) - (F*dRx)
                                cell_Fy(q) = cell_Fy(q) - (F*dRy)
                                cell_Fz(q) = cell_Fz(q) - (F*dRz)

                                pot_en = pot_en + (4 * RF6 * (RF6 - 1)) - pot_en_cut    ! potential energy
                            END IF

                        END DO

                    END DO

                    cell_Fxl = cell_Fx(1:l)                                             ! de-concatenate into l & s
                    cell_Fxs = cell_Fx(l+1:r)
                    cell_Fyl = cell_Fy(1:l)
                    cell_Fys = cell_Fy(l+1:r)
                    cell_Fzl = cell_Fz(1:l)
                    cell_Fzs = cell_Fz(l+1:r)

                    Fx = Fx + UNPACK(cell_Fxl, mask_l_force, real_zero)                 ! unpack into the force array
                    Fy = Fy + UNPACK(cell_Fyl, mask_l_force, real_zero)
                    Fz = Fz + UNPACK(cell_Fzl, mask_l_force, real_zero)

                    Fx = Fx + UNPACK(cell_Fxs, mask_s_force, real_zero)
                    Fy = Fy + UNPACK(cell_Fys, mask_s_force, real_zero)
                    Fz = Fz + UNPACK(cell_Fzs, mask_s_force, real_zero)

                    !DEALLOCATE( cell_Rx, cell_Ry, cell_Rz, cell_Fx, cell_Fy, cell_Fz )

                END DO
            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL


    END SUBROUTINE force_cell_list_KobAnd



    SUBROUTINE force_cell_list_2D(Rx, Ry, Fx, Fy, pot_en)

        ! Cell list algorithm

        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry

        INTEGER :: i, j, p, q, s, l, o, r
        REAL :: Rxp, Ryp, dRx, dRy, dRsq, RF2, RF6, F
        LOGICAL :: border_i_0, border_i_n, border_j
        LOGICAL, DIMENSION(npart) :: mask_l_force, mask_s_force
        REAL, DIMENSION(:), ALLOCATABLE :: cell_Rx, cell_Ry, cell_Fx, cell_Fy

        REAL, INTENT(OUT) :: pot_en
        REAL, DIMENSION(:), INTENT(OUT) :: Fx, Fy

        mask_l_force = .FALSE.
        mask_s_force = .FALSE.
        Fx = 0 ; Fy = 0
        pot_en = 0


        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(mask,Rx,Ry,Fx,Fy,pot_en)
        !$OMP DO REDUCTION(+:Fx,Fy,pot_en)
        DO j = 0, num_cell_edge-1
            border_j = (j == (num_cell_edge - 1))
            DO i = 0, num_cell_edge-1
                border_i_0 = (i == 0)
                border_i_n = (i == (num_cell_edge - 1))

                F = 0

                mask_l_force = mask(:,  cell_id_2D( j,   i )  )                     ! mask of the cell

                l = COUNT(mask_l_force)                                             ! no. of particles in the cell

                mask_s_force = mask(:,  cell_id_2D(    j, i+1 )  ) &               ! mask of half the surrounding cells
                          .OR. mask(:,  cell_id_2D(  j+1, i-1 )  ) &
                          .OR. mask(:,  cell_id_2D(  j+1,   i )  ) &
                          .OR. mask(:,  cell_id_2D(  j+1, i+1 )  )


                s = COUNT(mask_s_force)                                             ! no. of part. in half the surr. cells

                r = l+s

                ALLOCATE( cell_Rx(r), cell_Ry(r) )

                cell_Rx(1:l) = PACK(Rx, mask_l_force)                               ! packing the cell
                cell_Ry(1:l) = PACK(Ry, mask_l_force)

                cell_Rx(l+1:r) = PACK(Rx, mask_s_force)                             ! packing half the surrounding cells
                cell_Ry(l+1:r) = PACK(Ry, mask_s_force)

                ALLOCATE( cell_Fx(r), cell_Fy(r) )

                cell_Fx = [(real_zero, o=1, r)]
                cell_Fy = [(real_zero, o=1, r)]


                DO p = 1, l                                                         ! iterating over particles in the cell

                    Rxp = cell_Rx(p) ; Ryp = cell_Ry(p)

                    DO q = p+1, r                                                   ! iterate over all particles considered

                        dRx = Rxp - cell_Rx(q)                                      ! x-distance between two particles
                        dRy = Ryp - cell_Ry(q)

                        IF ( border_i_n ) THEN
                            IF ( dRx .GT. cell2 ) THEN
                                dRx = dRx - box                                     ! choosing nearest distance, due to PBCs
                            END IF
                        ELSEIF ( border_i_0 ) THEN
                            IF ( dRx .LT. cell2m ) THEN
                                dRx = dRx + box
                            END IF
                        END IF

                        IF ( border_j ) THEN
                            IF ( dRy .GT. cell2 ) THEN
                                dRy = dRy - box
                            END IF
                        END IF


                        dRsq = (dRx)**2 + (dRy)**2                                  ! distance between particles p and q
                        IF (dRsq .lt. RcSq) THEN
                            RF2 = SigSq/dRsq
                            RF6 = RF2 ** 3
                            F = 48 * (Eps/Sig) * RF2 * RF6 * (RF6 - 0.5)            ! force

                            cell_Fx(p) = cell_Fx(p) + (F*dRx)                       ! pairwise force addition
                            cell_Fy(p) = cell_Fy(p) + (F*dRy)
                            cell_Fx(q) = cell_Fx(q) - (F*dRx)
                            cell_Fy(q) = cell_Fy(q) - (F*dRy)

                            pot_en = pot_en + (4 * RF6 * (RF6 - 1)) - pot_en_cut    ! potential energy
                        END IF

                    END DO

                END DO

                Fx = Fx + UNPACK(cell_Fx(1:l), mask_l_force, real_zero)                 ! unpack into the force array
                Fy = Fy + UNPACK(cell_Fy(1:l), mask_l_force, real_zero)

                Fx = Fx + UNPACK(cell_Fx(l+1:r), mask_s_force, real_zero)
                Fy = Fy + UNPACK(cell_Fy(l+1:r), mask_s_force, real_zero)

                DEALLOCATE( cell_Rx, cell_Ry, cell_Fx, cell_Fy )


            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL


    END SUBROUTINE force_cell_list_2D



    INTEGER ELEMENTAL FUNCTION cell_id(c, b, a)

        ! Detemining the cell ID of a particle from its position

        INTEGER, INTENT(IN) :: a, b, c

        cell_id = (MODULO(a,num_cell_edge)+1) + (MODULO(b,num_cell_edge)*num_cell_edge) + (MODULO(c,num_cell_edge)*num_cell_face)

    END FUNCTION


    INTEGER ELEMENTAL FUNCTION cell_id_2D(b, a)

        ! Detemining the cell ID of a particle from its position

        INTEGER, INTENT(IN) :: a, b

        cell_id_2D = (MODULO(a,num_cell_edge)+1) + (MODULO(b,num_cell_edge)*num_cell_edge)

    END FUNCTION



    SUBROUTINE force_rdf(Rx, Ry, Rz, g)

        ! The radial distribution function

        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry, Rz
        INTEGER :: i, j, ri
        REAL :: Rxi, Ryi, Rzi, dRx, dRy, dRz, dR, dRsq
        REAL, DIMENSION(:), INTENT(OUT) :: g

        g = real_zero

        !!$OMP PARALLEL DEFAULT(SHARED) PRIVATE(Rxi,Ryi,Rzi,dRx,dRy,dRz,dRsq,dR,RF2,RF6,F,ri)
        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Rx,Ry,Rz,g)
        !$OMP DO REDUCTION(+:g)
        DO i = 1, npart-1                                                   ! pair-wise force calculation

            Rxi = Rx(i) ; Ryi = Ry(i) ; Rzi = Rz(i)

            DO j = i+1, npart

                dRx = Rxi - Rx(j)                                           ! x-distance between two particles
                dRx = dRx - (box*NINT(dRx/box))                             ! choosing nearest distance, due to periodic boudaries
                dRy = Ryi - Ry(j)
                dRy = dRy - (box*NINT(dRy/box))
                dRz = Rzi - Rz(j)
                dRz = dRz - (box*NINT(dRz/box))

                dRsq = (dRx)**2 + (dRy)**2 + (dRz)**2                       ! distance between particles i and j
                dR = SQRT(dRsq)

                IF (dR .lt. box2) THEN
                    ri = CEILING(dR/size_b)                                 ! populating g(r)
                    g(ri) = g(ri) + 2
                END IF

            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE force_rdf


    SUBROUTINE force_rdf_2D(Rx, Ry, g)

        ! The radial distribution function

        REAL, DIMENSION(:), INTENT(IN) :: Rx, Ry
        INTEGER :: i, j, ri
        REAL :: Rxi, Ryi, dRx, dRy, dR, dRsq
        REAL, DIMENSION(:), INTENT(OUT) :: g

        g = real_zero

        !$OMP PARALLEL DEFAULT(PRIVATE) SHARED(Rx,Ry,g)
        !$OMP DO REDUCTION(+:g)
        DO i = 1, npart-1                                                   ! pair-wise force calculation

            Rxi = Rx(i) ; Ryi = Ry(i)

            DO j = i+1, npart

                dRx = Rxi - Rx(j)                                           ! x-distance between two particles
                dRx = dRx - (box*NINT(dRx/box))                             ! choosing nearest distance, due to periodic boudaries
                dRy = Ryi - Ry(j)
                dRy = dRy - (box*NINT(dRy/box))

                dRsq = (dRx)**2 + (dRy)**2                                  ! distance between particles i and j
                dR = SQRT(dRsq)

                IF (dR .lt. box2) THEN
                    ri = CEILING(dR/size_b)                                 ! populating g(r)
                    g(ri) = g(ri) + 2
                END IF

            END DO
        END DO
        !$OMP END DO
        !$OMP END PARALLEL

    END SUBROUTINE force_rdf_2D


END MODULE force
