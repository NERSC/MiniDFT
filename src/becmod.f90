!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!
MODULE becmod
  !
  ! ... *bec* contain <beta|psi> - used in h_psi, s_psi, many other places
  ! ... calbec( npw, beta, psi, betapsi [, nbnd ] ) is an interface calculating
  ! ...    betapsi(i,j)  = <beta(i)|psi(j)>   (the sum is over npw components)
  ! ... or betapsi(i,s,j)= <beta(i)|psi(s,j)> (s=polarization index)
  !
  USE kinds,            ONLY : DP
  USE gvect,            ONLY : gstart
  !
  SAVE
  !
#ifdef __STD_F95
  TYPE bec_type
     COMPLEX(DP),POINTER :: k(:)    ! appropriate for generic k
     INTEGER :: comm
     INTEGER :: nproc
     INTEGER :: mype
     INTEGER :: ibnd_begin
  END TYPE bec_type
#else
  TYPE bec_type
     COMPLEX(DP),ALLOCATABLE :: k(:)    ! appropriate for generic k
     INTEGER :: comm
     INTEGER :: nproc
     INTEGER :: mype
     INTEGER :: ibnd_begin
  END TYPE bec_type
#endif
  !
  TYPE (bec_type) :: becp  ! <beta|psi>

  PRIVATE

  COMPLEX(DP), ALLOCATABLE ::  &
       becp_k (:)    !  as above for complex wavefunctions
  !
  INTERFACE calbec
     !
     MODULE PROCEDURE calbec_k, calbec_bec_type
     !
  END INTERFACE

  !
  PUBLIC :: bec_type, becp, allocate_bec_type, deallocate_bec_type, calbec, &
            beccopy, becscal
  !
CONTAINS
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_bec_type ( npw, beta, psi, betapsi, ibnd )
    !-----------------------------------------------------------------------
    !_
    USE mp_global, ONLY: intra_bgrp_comm
    USE mp, ONLY: mp_size, mp_rank, mp_get_comm_null
    !
    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    TYPE (bec_type), INTENT (inout) :: betapsi ! NB: must be INOUT otherwise
                                               !  the allocatd array is lost
    INTEGER, INTENT (in) :: npw
    INTEGER, INTENT (in) :: ibnd
    !
    INTEGER, EXTERNAL :: ldim_block, lind_block, gind_block
    INTEGER :: nproc, mype, m_loc, m_begin, m_max, ip
    REAL(DP), ALLOCATABLE :: dtmp(:,:)
    !

       !
       CALL  calbec_k ( npw, beta, psi, betapsi%k, ibnd )
       !
    !
    RETURN
    !
  END SUBROUTINE calbec_bec_type
  !
  !-----------------------------------------------------------------------
  SUBROUTINE calbec_k ( npw, beta, psi, betapsi, ibnd )
    !-----------------------------------------------------------------------
    !
    ! ... matrix times matrix with summation index (k=1,npw) running on
    ! ... G-vectors or PWs : betapsi(i,j) = \sum_k beta^*(i,k) psi(k,j)
    !
    USE mp_global, ONLY : intra_bgrp_comm
    USE mp,        ONLY : mp_sum

    IMPLICIT NONE
    COMPLEX (DP), INTENT (in) :: beta(:,:), psi(:,:)
    COMPLEX (DP), INTENT (out) :: betapsi(:)
    INTEGER, INTENT (in) :: npw
    INTEGER, INTENT (in) :: ibnd
    !
    INTEGER :: nkb, npwx, m
    !
    nkb = size (beta, 2)
    IF ( nkb == 0 ) RETURN
    !
    CALL start_clock( 'calbec' )
    npwx= size (beta, 1)
    IF ( npwx /= size (psi, 1) ) CALL errore ('calbec', 'size mismatch', 1)
    IF ( npwx < npw ) CALL errore ('calbec', 'size mismatch', 2)
    IF ( nkb /= size (betapsi,1) ) &
      CALL errore ('calbec', 'size mismatch', 3)

       !
       CALL ZGEMV( 'C', npw, nkb, (1.0_DP,0.0_DP), beta, npwx, psi(:,ibnd), 1, &
                   (0.0_DP, 0.0_DP), betapsi, 1 )

    CALL mp_sum( betapsi( : ), intra_bgrp_comm )
    !
    CALL stop_clock( 'calbec' )
    !
    RETURN
    !
  END SUBROUTINE calbec_k
  !
  !-----------------------------------------------------------------------
  SUBROUTINE allocate_bec_type ( nkb, bec, comm )
    !-----------------------------------------------------------------------
    USE mp, ONLY: mp_size, mp_rank, mp_get_comm_null
    IMPLICIT NONE
    TYPE (bec_type) :: bec
    INTEGER, INTENT (in) :: nkb
    INTEGER, INTENT (in), OPTIONAL :: comm
    INTEGER :: ierr
    INTEGER, EXTERNAL :: ldim_block, lind_block, gind_block
    !
#ifdef __STD_F95
    NULLIFY(bec%k)
#endif
    !
    bec%comm = mp_get_comm_null()
    bec%mype = 0
    bec%nproc = 1
    bec%ibnd_begin = 1
    !
       !
       ALLOCATE( bec%k( nkb ), STAT=ierr )
       IF( ierr /= 0 ) &
          CALL errore( ' allocate_bec_type ', ' cannot allocate bec%k ', ABS(ierr) )
       !
       bec%k(:)=(0.0D0,0.0D0)
       !
    !
    RETURN
    !
  END SUBROUTINE allocate_bec_type
  !
  !-----------------------------------------------------------------------
  SUBROUTINE deallocate_bec_type (bec)
    !-----------------------------------------------------------------------
    !
    USE mp, ONLY: mp_get_comm_null
    IMPLICIT NONE
    TYPE (bec_type) :: bec
    !
    bec%comm = mp_get_comm_null()
    !
#ifdef __STD_F95
    IF (associated(bec%k))  DEALLOCATE(bec%k)
#else
    IF (allocated(bec%k))  DEALLOCATE(bec%k)
#endif
    !
    RETURN
    !
  END SUBROUTINE deallocate_bec_type

  SUBROUTINE beccopy(bec, bec1, nkb)
    IMPLICIT NONE
    TYPE(bec_type), INTENT(in) :: bec
    TYPE(bec_type)  :: bec1
    INTEGER, INTENT(in) :: nkb

       CALL zcopy(nkb, bec%k, 1, bec1%k, 1)

    RETURN
  END SUBROUTINE beccopy

  SUBROUTINE becscal_nck(alpha, bec, nkb)
    IMPLICIT NONE
    TYPE(bec_type), INTENT(INOUT) :: bec
    COMPLEX(DP), INTENT(IN) :: alpha
    INTEGER, INTENT(IN) :: nkb

       CALL zscal(nkb, alpha, bec%k, 1)

    RETURN
  END SUBROUTINE becscal_nck


END MODULE becmod
