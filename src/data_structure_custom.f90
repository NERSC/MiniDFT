
!
! Copyright (C) 2001-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE data_structure_custom(fc, gamma_only)
  !-----------------------------------------------------------------------
  ! this routine sets the data structure for the custom fft array
  ! In the parallel case, it distributes columns to processes, too
  !
  USE kinds,      ONLY : DP
  USE cell_base,  ONLY : bg, tpiba, tpiba2
  USE klist,      ONLY : xk, nks
  USE mp,         ONLY : mp_sum, mp_max,mp_barrier
  USE mp_global,  ONLY : mpime, me_bgrp, nproc_bgrp, inter_bgrp_comm,&
                         & intra_bgrp_comm, root_bgrp, inter_pool_comm
  USE mp_global,  ONLY : get_ntask_groups 
  USE stick_set,  ONLY : pstickset_custom
  USE fft_custom, ONLY : fft_cus, gvec_init
  !
  !
  IMPLICIT NONE
  
  TYPE(fft_cus) :: fc
  LOGICAL :: gamma_only
  REAL (DP) :: gkcut
  INTEGER :: ik, ngm_, ngs_, ngw_ , nogrp
  INTEGER :: me, nproc, inter_comm, intra_comm, root

  INTEGER :: kpoint
  ! sticks coordinates
  
  !
  !  Subroutine body
  !

  !
  ! compute gkcut calling an internal procedure
  !

  me = me_bgrp
  nproc = nproc_bgrp
  inter_comm = inter_pool_comm
  intra_comm = intra_bgrp_comm
  root = root_bgrp

  nogrp = get_ntask_groups()

    IF (nks == 0) THEN
       !
       ! if k-points are automatically generated (which happens later)
       ! use max(bg)/2 as an estimate of the largest k-point
       !
       gkcut = 0.5d0 * MAX ( &
            &SQRT (SUM(bg (1:3, 1)**2) ), &
            &SQRT (SUM(bg (1:3, 2)**2) ), &
            &SQRT (SUM(bg (1:3, 3)**2) ) )
    ELSE
       gkcut = 0.0d0
       DO kpoint = 1, nks
          gkcut = MAX (gkcut, SQRT ( SUM(xk (1:3, kpoint)**2) ) )
       ENDDO
    ENDIF
    gkcut = (SQRT (fc%ecutt) / tpiba + gkcut)**2

  !
  ! ... find maximum value among all the processors
  !
  CALL mp_max (gkcut, inter_comm )
  !
  ! ... set up fft descriptors, including parallel stuff: sticks, planes, etc.
  !
  nogrp = get_ntask_groups()
  !
  CALL pstickset_custom( gamma_only, bg, fc%gcutmt, gkcut, &
                  fc%dfftt, ngw_ , ngm_ , me, root, nproc, intra_comm,   &
                  nogrp )
  !
  !     on output, ngm_ and ngs_ contain the local number of G-vectors
  !     for the two grids. Initialize local and global number of G-vectors
  !
  CALL gvec_init (fc, ngm_ , intra_comm )

  
END SUBROUTINE data_structure_custom
