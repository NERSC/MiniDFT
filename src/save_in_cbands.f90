!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!-----------------------------------------------------------------------
subroutine save_in_cbands (iter, ik_, dr2)
  !-----------------------------------------------------------------------
  USE kinds,         ONLY: DP
  USE io_files,      ONLY: iunres, prefix, seqopn
  USE klist,         ONLY: nks
  USE control_flags, ONLY: io_level, tr2, ethr
  USE wvfct,         ONLY: nbnd, et
  implicit none
  character :: where * 20
  ! are we in the right place?
  integer :: ik, ibnd, ik_, iter
  ! counters
  ! last completed kpoint
  ! last completed iteration
  logical :: exst

  real(DP) :: dr2
  !
  ! open recover file
  !
  call seqopn (iunres, 'restart', 'unformatted', exst)
  !
  ! save restart information
  !
  where = 'ELECTRONS'
  write (iunres) where
  write (iunres) ( (et (ibnd, ik), ibnd = 1, nbnd), ik = 1, nks)

  write (iunres) iter, ik_, dr2, tr2, ethr

  close (unit = iunres, status = 'keep')
  !
  return
end subroutine save_in_cbands
