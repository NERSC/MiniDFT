!
! Copyright (C) 2005-2012 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------
MODULE exx
  !--------------------------------------
  !
  USE kinds,                ONLY : DP
  USE coulomb_vcut_module,  ONLY : vcut_init, vcut_type, vcut_info, &
                                   vcut_get,  vcut_spheric_get
  USE io_global,            ONLY : ionode, stdout
  USE fft_custom,           ONLY : fft_cus
  !
  IMPLICIT NONE
  SAVE

  !
  ! general purpose vars
  !
  REAL(DP):: exxalfa=0.d0                ! 1 if exx, 0 elsewhere
  INTEGER :: exx_nwordwfc, ji

  !
  ! variables defining the auxiliary k-point grid 
  ! used in X BZ integration
  !
  INTEGER :: nq1=1, nq2=1, nq3=1         ! integers defining the X integration mesh
  INTEGER :: nqs=1                       ! number of points in the q-gridd
  INTEGER :: nkqs                        ! total number of different k+q
  !
  REAL(DP),    ALLOCATABLE :: xkq(:,:)   ! xkq(3,nkqs) the auxiliary k+q set
  REAL(DP),    ALLOCATABLE :: x_occupation(:,:)           
                                         ! x_occupation(nbnd,nks) the weight of 
                                         ! auxiliary functions in the density matrix
  COMPLEX(DP), ALLOCATABLE :: exxbuff(:,:,:)
                                         ! temporay buffer to store wfc 
  !
  ! let xk(:,ik) + xq(:,iq) = xkq(:,ikq) = S(isym)*xk(ik') + G
  ! 
  !     index_xkq(ik,iq) = ikq
  !     index_xk(ikq)    = ik'
  !     index_sym(ikq)   = isym
  !
  INTEGER, ALLOCATABLE :: index_xkq(:,:) ! index_xkq(nks,nqs) 
  INTEGER, ALLOCATABLE :: index_xk(:)    ! index_xk(nkqs)  
  INTEGER, ALLOCATABLE :: index_sym(:)   ! index_sym(nkqs)
!
!  Used for k points pool parallelization. All pools need these quantities.
!  They are allocated only if needed.
!
  REAL(DP),    ALLOCATABLE :: xk_collect(:,:)
  REAL(DP),    ALLOCATABLE :: wk_collect(:)
  REAL(DP),    ALLOCATABLE :: wg_collect(:,:)
  LOGICAL :: pool_para=.FALSE.
  !
  ! variables to deal with Coulomb divergence
  ! and related issues
  !
  REAL (DP)         :: eps =1.d-6
  REAL (DP)         :: exxdiv = 0.d0
  CHARACTER(80)     :: exxdiv_treatment 
  !
  ! x_gamma_extrapolation
  LOGICAL           :: x_gamma_extrapolation =.TRUE.
  LOGICAl           :: on_double_grid =.FALSE.
  REAL (DP)         :: grid_factor = 8.d0/7.d0 
  !
  ! Gygi-Baldereschi 
  LOGICAL           :: use_regularization = .TRUE.
  !
  ! yukawa method
  REAL (DP)         :: yukawa = 0.d0
  !
  ! erfc screening
  REAL (DP)         :: erfc_scrlen = 0.d0
  !
  ! erf screening
  REAL (DP)         :: erf_scrlen = 0.d0
  ! cutoff techniques
  LOGICAL           :: use_coulomb_vcut_ws = .FALSE.
  LOGICAL           :: use_coulomb_vcut_spheric = .FALSE.
  REAL (DP)         :: ecutvcut
  TYPE(vcut_type)   :: vcut

  !
  ! energy related variables
  !
  REAL(DP) :: fock0 = 0.0_DP, & !   sum <phi|Vx(phi)|phi>
              fock1 = 0.0_DP, & !   sum <psi|vx(phi)|psi>
              fock2 = 0.0_DP, & !   sum <psi|vx(psi)|psi>
              dexx  = 0.0_DP    !   fock1  - 0.5*(fock2+fock0)

  !
  ! custom fft grids
  !
  TYPE(fft_cus) exx_fft_g2r     ! Grid for wfcs -> real space
  TYPE(fft_cus) exx_fft_r2g     ! Grid for real space -> restricted G space
  REAL(DP)  :: ecutfock         ! energy cutoff for custom grid
  REAL(DP)  :: exx_dual = 4.0_DP! dual for the custom grid
CONTAINS
  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_convert( psi, npw, fft, psi_t, sign, igkt )
  !------------------------------------------------------------------------

    USE mp_global,  ONLY : me_bgrp, nproc_bgrp, intra_bgrp_comm, root_bgrp
    USE mp_wave,    ONLY : mergewf, splitwf
    USE gvect, ONLY : ig_l2g

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: npw
    COMPLEX(kind=DP), INTENT(IN) :: psi(npw)
    COMPLEX(kind=DP), INTENT(INOUT) :: psi_t(:)
    INTEGER, OPTIONAL, INTENT(INOUT) :: igkt(:)
    INTEGER, INTENT(IN) :: sign
    TYPE(fft_cus), INTENT(IN) :: fft
    
    COMPLEX(kind=DP), ALLOCATABLE :: evc_g(:)
    INTEGER :: ig
    
    ALLOCATE( evc_g( fft%ngmt_g ) )
    
    IF(sign > 0 .AND. PRESENT(igkt) ) THEN
       DO ig=1, fft%ngmt
          igkt(ig)=ig
       ENDDO
    ENDIF
    
    IF( fft%dual_t==4.d0) THEN
       psi_t(1:fft%npwt)=psi(1:fft%npwt)
    ELSE
       IF (sign > 0 ) THEN
          CALL mergewf(psi, evc_g, npw, ig_l2g, me_bgrp, nproc_bgrp,&
               & root_bgrp, intra_bgrp_comm)  
          CALL splitwf(psi_t(:), evc_g, fft%npwt, fft%ig_l2gt, me_bgrp,&
               & nproc_bgrp, root_bgrp, intra_bgrp_comm)  
       ELSE
          CALL mergewf(psi, evc_g, fft%npwt, fft%ig_l2gt, me_bgrp,&
               & nproc_bgrp, root_bgrp, intra_bgrp_comm)  
          CALL splitwf(psi_t, evc_g, npw, ig_l2g, me_bgrp, nproc_bgrp,&
               & root_bgrp, intra_bgrp_comm)  
       ENDIF
    ENDIF

    DEALLOCATE( evc_g ) 

  END SUBROUTINE exx_grid_convert
  !------------------------------------------------------------------------
  SUBROUTINE exx_fft_create ()
  !------------------------------------------------------------------------
    
    USE wvfct,        ONLY : ecutwfc
    USE gvect,        ONLY : ecutrho

    IMPLICIT NONE

    IF(ecutfock <= 0.0_DP) ecutfock = ecutrho

    IF(ecutfock < ecutwfc) CALL errore('exx_fft_create', 'ecutfock can&
         &not be smaller than ecutwfc!', 1) 

    ! Initalise the g2r grid that allows us to put the wavefunction
    ! onto the new (smaller) grid for rho.
    exx_fft_g2r%ecutt=ecutwfc
    exx_fft_g2r%dual_t=ecutfock/ecutwfc
    CALL allocate_fft_custom(exx_fft_g2r)

    ! Initalise the r2g grid that we then use when applying the Fock
    ! operator in our new restricted space.
    exx_fft_r2g%ecutt=ecutfock/exx_dual
    exx_fft_r2g%dual_t=exx_dual
    CALL allocate_fft_custom(exx_fft_r2g)

  END SUBROUTINE exx_fft_create
  !------------------------------------------------------------------------
  SUBROUTINE exx_fft_destroy ()
  !------------------------------------------------------------------------
    USE fft_custom,  ONLY : deallocate_fft_custom

    IMPLICIT NONE

    CALL deallocate_fft_custom(exx_fft_g2r)
    CALL deallocate_fft_custom(exx_fft_r2g)

  END SUBROUTINE exx_fft_destroy
  !------------------------------------------------------------------------
  SUBROUTINE deallocate_exx ()
  !------------------------------------------------------------------------
  !
  IF ( ALLOCATED (index_xkq) ) DEALLOCATE (index_xkq)
  IF ( ALLOCATED (index_xk ) ) DEALLOCATE (index_xk )
  IF ( ALLOCATED (index_sym) ) DEALLOCATE (index_sym)
  IF ( ALLOCATED (x_occupation) ) DEALLOCATE (x_occupation)
  IF ( ALLOCATED (xkq) ) DEALLOCATE (xkq)
  IF ( ALLOCATED (exxbuff) ) DEALLOCATE (exxbuff)
  !
  CALL exx_fft_destroy()
  !
  !  Pool variables deallocation
  !
  IF ( ALLOCATED (xk_collect) )  DEALLOCATE ( xk_collect )
  IF ( ALLOCATED (wk_collect) )  DEALLOCATE ( wk_collect )
  IF ( ALLOCATED (wg_collect) )  DEALLOCATE ( wg_collect )
  !
  !
  END SUBROUTINE deallocate_exx
  !------------------------------------------------------------------------
  subroutine exx_grid_init()
  !------------------------------------------------------------------------
  !
  USE symm_base,  ONLY : nsym, s
  USE cell_base,  ONLY : bg, at, alat
  USE lsda_mod,   ONLY : nspin
  USE spin_orb,   ONLY : domag
  USE klist,      ONLY : xk, wk, nkstot, nks
  USE wvfct,      ONLY : nbnd
  USE io_global,  ONLY : stdout
  !
  USE mp_global,  ONLY : nproc, npool, nimage
  !
  IMPLICIT NONE
  !
  CHARACTER(13) :: sub_name='exx_grid_init'
  integer       :: iq1, iq2, iq3, isym, ik, ikq, iq, max_nk, temp_nkqs
  integer,   allocatable :: temp_index_xk(:), temp_index_sym(:)
  integer,   allocatable :: temp_index_ikq(:), new_ikq(:)
  real (DP), allocatable :: temp_xkq(:,:)
  logical       :: xk_not_found
  real (DP)     :: sxk(3), dxk(3), xk_cryst(3)
  real (DP)     :: dq1, dq2, dq3
  logical       :: no_pool_para
  integer       :: find_current_k

  CALL start_clock ('exx_grid')

  !
  ! definitions and checks
  !
#ifdef EXXDEBUG
  IF (ionode) WRITE(stdout,'(/,2x,a,3i4)') "EXX : q-grid dimensions are ", nq1,nq2,nq3
#endif
  !
  grid_factor = 1.d0
  !
  IF (x_gamma_extrapolation) THEN
      !
#ifdef EXXDEBUG
      IF (ionode) WRITE (stdout,'(2x,a)') "EXX : q->0 dealt with 8/7 -1/7 trick"
#endif
      grid_factor = 8.d0/7.d0
      !
  ENDIF
  

  !
  nqs = nq1 * nq2 * nq3
  !
  ! all processors need to have access to all k+q points
  !
  pool_para =  npool>1
  if (pool_para ) then
     IF ( .NOT.ALLOCATED (xk_collect) )  ALLOCATE (xk_collect(3,nkstot))
     IF ( .NOT.ALLOCATED (wk_collect) )  ALLOCATE (wk_collect(nkstot))
     CALL xk_wk_collect(xk_collect, wk_collect, xk, wk, nkstot, nks)
  end if
  !
  ! set a safe limit as the maximum number of auxiliary points we may need
  ! and allocate auxiliary arrays
  !
  max_nk = nkstot * min(48, 2 * nsym)
  allocate ( temp_index_xk(max_nk), temp_index_sym(max_nk) )
  allocate ( temp_index_ikq(max_nk), new_ikq(max_nk) )
  allocate ( temp_xkq(3,max_nk) )
  !
  ! find all k-points equivalent by symmetry to the points in the k-list
  !
  temp_nkqs = 0
  do isym=1,nsym
     do ik =1, nkstot
        IF (pool_para) THEN
           xk_cryst(:) = at(1,:)*xk_collect(1,ik) + at(2,:)*xk_collect(2,ik)&
                       + at(3,:)*xk_collect(3,ik)
        ELSE
           xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) &
                       + at(3,:)*xk(3,ik)
        ENDIF
        sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                 s(:,2,isym)*xk_cryst(2) + &
                 s(:,3,isym)*xk_cryst(3)
        ! add sxk to the auxiliary list if it is not already present
        xk_not_found = .true.
        do ikq=1, temp_nkqs
           if (xk_not_found ) then
              dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
              if ( abs(dxk(1)).le.eps .and. &
                   abs(dxk(2)).le.eps .and. &
                   abs(dxk(3)).le.eps ) xk_not_found = .false.
           end if
        end do
        if (xk_not_found) then
           temp_nkqs                 = temp_nkqs + 1
           temp_xkq(:,temp_nkqs)     = sxk(:)
           temp_index_xk(temp_nkqs)  = ik
           temp_index_sym(temp_nkqs) = isym 
        end if

        sxk(:) = - sxk(:)
        xk_not_found = .true.
        do ikq=1, temp_nkqs
           if (xk_not_found ) then
              dxk(:) = sxk(:) - temp_xkq(:,ikq) - nint(sxk(:) - temp_xkq(:,ikq))
              if ( abs(dxk(1)).le.eps .and. &
                   abs(dxk(2)).le.eps .and. &
                   abs(dxk(3)).le.eps ) xk_not_found = .false.
           end if
        end do
        if (xk_not_found .and. .not. (domag) ) then
           temp_nkqs                 = temp_nkqs + 1
           temp_xkq(:,temp_nkqs)     = sxk(:)
           temp_index_xk(temp_nkqs)  = ik
           temp_index_sym(temp_nkqs) =-isym 
        end if

     end do
  end do

  !
  ! define the q-mesh step-sizes
  !
  dq1= 1.d0/DBLE(nq1)
  dq2= 1.d0/DBLE(nq2)
  dq3= 1.d0/DBLE(nq3)
  !
  ! allocate and fill the array index_xkq(nkstot,nqs)
  !
  if(.not.ALLOCATED(index_xkq)) allocate ( index_xkq(nkstot,nqs) )
  if(.not.ALLOCATED(x_occupation)) allocate ( x_occupation(nbnd,nkstot) )
  nkqs = 0
  new_ikq(:) = 0
  do ik=1,nkstot 
     IF (pool_para) THEN
        xk_cryst(:) = at(1,:)*xk_collect(1,ik) + at(2,:)*xk_collect(2,ik)&
                    + at(3,:)*xk_collect(3,ik)
     ELSE
        xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)
     ENDIF

     iq = 0
     do iq1=1, nq1
       sxk(1) = xk_cryst(1) + (iq1-1) * dq1
       do iq2 =1, nq2
         sxk(2) = xk_cryst(2) + (iq2-1) * dq2
         do iq3 =1, nq3
            sxk(3) = xk_cryst(3) + (iq3-1) * dq3
            iq = iq + 1

            xk_not_found = .true.
            do ikq=1, temp_nkqs
               if ( xk_not_found ) then
                  dxk(:) = sxk(:)-temp_xkq(:,ikq) - nint(sxk(:)-temp_xkq(:,ikq))
                  if ( abs(dxk(1)).le.eps .and. &
                       abs(dxk(2)).le.eps .and. &
                       abs(dxk(3)).le.eps ) then

                       xk_not_found = .false.

                       if( new_ikq(ikq) == 0) then
                           nkqs = nkqs + 1
                           temp_index_ikq(nkqs) = ikq
                           new_ikq(ikq) = nkqs
                       end if
                       index_xkq(ik,iq) = new_ikq(ikq)

                  end if
               end if
            end do
            if (xk_not_found) then
               write (*,*) ik, iq, temp_nkqs
               write (*,*) sxk(:)
               call errore('exx_grid_init', &
                           ' k + q is not an S*k ', (ik-1) * nqs + iq )
            end if

         end do
       end do
     end do

  end do
  !
  ! allocate and fill the arrays xkq(3,nkqs), index_xk(nkqs) and index_sym(nkqs)
  !
  allocate ( xkq(3,nkqs), index_xk(nkqs),  &
             index_sym(nkqs) )

  do ik =1, nkqs
     ikq = temp_index_ikq(ik)
     xkq(:,ik) = bg(:,1)*temp_xkq(1,ikq) + &
                 bg(:,2)*temp_xkq(2,ikq) + &
                 bg(:,3)*temp_xkq(3,ikq)
     index_xk(ik)  = temp_index_xk(ikq)
     index_sym(ik) = temp_index_sym(ikq)
  end do

  IF (nspin == 2) THEN
     DO ik = 1, nkstot/2
        DO iq =1, nqs
           index_xkq(nkstot/2+ik,iq) = index_xkq(ik,iq) + nkqs
        END DO
     ENDDO
     do ikq=1,nkqs
        xkq(:,ikq + nkqs)     = xkq(:,ikq)
        index_xk(ikq + nkqs)  = index_xk(ikq) + nkstot/2
        index_sym(ikq + nkqs) = index_sym(ikq)
     end do
     nkqs = 2 * nkqs
  ENDIF
  !
  ! clean up
  !
  deallocate (temp_index_xk, temp_index_sym, temp_index_ikq, new_ikq, temp_xkq)
  !
  ! check that everything is what it should be
  !
  call exx_grid_check () 

  CALL stop_clock ('exx_grid')
  !
  RETURN
  END SUBROUTINE exx_grid_init

  !------------------------------------------------------------------------
  subroutine exx_div_check()
  !------------------------------------------------------------------------
  !
  USE cell_base,  ONLY : bg, at, alat
  USE io_global,  ONLY : stdout
  USE funct,      ONLY : get_screening_parameter
  !
  IMPLICIT NONE
  !
  REAL (DP)     :: atws(3,3)
  CHARACTER(13) :: sub_name='exx_div_check'

  !
  ! EXX singularity treatment
  !
  SELECT CASE ( TRIM(exxdiv_treatment) ) 
  CASE ( "gygi-baldereschi", "gygi-bald", "g-b" )
     !
     use_regularization = .TRUE.
     !
     !
  CASE ( "vcut_ws" )
     !
     use_coulomb_vcut_ws = .TRUE.
     IF ( x_gamma_extrapolation ) &
          CALL errore(sub_name,'cannot use x_gamm_extrap and vcut_ws', 1)
     !
  CASE ( "vcut_spherical" ) 
     !
     use_coulomb_vcut_spheric = .TRUE.
     IF ( x_gamma_extrapolation ) &
          CALL errore(sub_name,'cannot use x_gamm_extrap and vcut_spherical', 1)
     !
  CASE ( "none" )
     use_regularization = .FALSE.
     !
  CASE DEFAULT
     CALL errore(sub_name,'invalid exxdiv_treatment: '//TRIM(exxdiv_treatment), 1)
  END SELECT
  !
#ifdef EXXDEBUG
  IF ( ionode ) WRITE (stdout,'(2x,"EXX : q->0 dealt with ",a, " trick" )') &
                TRIM(exxdiv_treatment)
#endif
  !
  ! <AF>
  ! Set variables for Coulomb vcut
  ! NOTE: some memory is allocated inside this routine (in the var vcut)
  !       and should be deallocated somewehre, at the end of the run
  !
  IF ( use_coulomb_vcut_ws .OR. use_coulomb_vcut_spheric ) THEN
      !
      ! build the superperiodicity direct lattice
      !
      atws = alat * at
      !
      atws(:,1) = atws(:,1) * nq1
      atws(:,2) = atws(:,2) * nq2
      atws(:,3) = atws(:,3) * nq3
      !
      !CALL start_clock ('exx_vcut_init')
      CALL vcut_init( vcut, atws, ecutvcut )
      !CALL stop_clock ('exx_vcut_init')
      !
      IF ( ionode ) CALL vcut_info( stdout, vcut )
      !          
  ENDIF
#ifdef EXXDEBUG
  write (stdout,"(2x,'EXX : exx div treatment check successful')")
#endif 
  RETURN
  END SUBROUTINE exx_div_check 


  !------------------------------------------------------------------------
  SUBROUTINE exx_grid_check ( )
  !------------------------------------------------------------------------
  USE symm_base, ONLY : s
  USE cell_base, ONLY : bg, at
  USE lsda_mod,  ONLY : nspin
  USE io_global, ONLY : stdout
  USE klist,     ONLY : nkstot, xk
  implicit none
  real (DP) :: sxk(3), dxk(3), xk_cryst(3), xkk_cryst(3)
  integer :: iq1, iq2, iq3, isym, ik, ikk, ikq, iq
  real (DP) :: eps, dq1, dq2, dq3
  eps = 1.d-6
  dq1= 1.d0/DBLE(nq1)
  dq2= 1.d0/DBLE(nq2)
  dq3= 1.d0/DBLE(nq3)

  do ik =1, nkstot
     IF (pool_para) THEN
        xk_cryst(:) = at(1,:)*xk_collect(1,ik) + at(2,:)*xk_collect(2,ik) + &
                      at(3,:)*xk_collect(3,ik)
     ELSE
        xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + at(3,:)*xk(3,ik)
     ENDIF

     iq = 0
     do iq1=1, nq1
       sxk(1) = xk_cryst(1) + (iq1-1) * dq1
       do iq2 =1, nq2
         sxk(2) = xk_cryst(2) + (iq2-1) * dq2
         do iq3 =1, nq3
            sxk(3) = xk_cryst(3) + (iq3-1) * dq3
            iq = iq + 1
            
            ikq  = index_xkq(ik,iq) 
            ikk  = index_xk(ikq)
            isym = index_sym(ikq)

            IF (pool_para) THEN
               xkk_cryst(:) = at(1,:)*xk_collect(1,ikk)+at(2,:)* &
                                xk_collect(2,ikk)+at(3,:)*xk_collect(3,ikk)
            ELSE
               xkk_cryst(:) = at(1,:)*xk(1,ikk)+at(2,:)*xk(2,ikk) &
                            + at(3,:)*xk(3,ikk)
            ENDIF

            if (isym < 0 ) xkk_cryst(:) = - xkk_cryst(:)
            isym = abs (isym)
            dxk(:) = s(:,1,isym)*xkk_cryst(1) + &
                     s(:,2,isym)*xkk_cryst(2) + &
                     s(:,3,isym)*xkk_cryst(3) - sxk(:)
            dxk(:) = dxk(:) - nint(dxk(:))
            if ( .not. ( abs(dxk(1)).le.eps .and. &
                         abs(dxk(2)).le.eps .and. &
                         abs(dxk(3)).le.eps )   ) then
                 write(*,*) ik,iq
                 write(*,*) ikq,ikk,isym
                 write(*,*) dxk(:)
                 call errore('exx_grid_check', &
                             'something wrong', 1 )
            end if

         end do
       end do
     end do
  end do
#ifdef EXXDEBUG
  write (stdout,"(2x,'EXX : grid check successful')")
#endif
  return

  end subroutine exx_grid_check

  !------------------------------------------------------------------------
  subroutine exx_restart(l_exx_was_active)
  !------------------------------------------------------------------------
    !This subroutine is called when restarting an exx calculation
    use funct,                ONLY : get_exx_fraction, start_exx, exx_is_active, &
                                     get_screening_parameter
    USE fft_base,             ONLY : dffts
    USE io_global,            ONLY : stdout

    implicit none
    logical, intent(in) :: l_exx_was_active
    logical :: exst

    if (.not. l_exx_was_active ) return ! nothing had happpened yet
    !!
    exx_nwordwfc=2*dffts%nnr
    !iunexx = find_free_unit()
    !call diropn(iunexx,'exx', exx_nwordwfc, exst) 
    erfc_scrlen = get_screening_parameter()
    exxdiv = exx_divergence() 
    exxalfa = get_exx_fraction()
#ifdef EXXDEBUG
    write (stdout,*) " ! EXXALFA SET TO ", exxalfa
#endif
    call start_exx
    call weights()
    call exxinit()
    fock0 = exxenergy2()
 
    return
  end subroutine exx_restart

!------------------------------------------------------------------------
subroutine exxinit()
!------------------------------------------------------------------------

    !This subroutine is run before the first H_psi() of each iteration.
    !It saves the wavefunctions for the right density matrix. in real space
    !It saves all the wavefunctions in a single file called prefix.exx
    !
    USE wavefunctions_module, ONLY : evc  
    USE io_files,             ONLY : nwordwfc, iunwfc, iunigk, &
                                     tmp_dir, prefix
    USE io_global,            ONLY : stdout
    USE buffers,              ONLY : get_buffer
    USE gvecs,              ONLY : nls, nlsm, doublegrid
    USE wvfct,                ONLY : nbnd, npwx, npw, igk, wg, et
    USE control_flags,        ONLY : gamma_only
    USE klist,                ONLY : wk, ngk, nks, nkstot
    USE symm_base,            ONLY : nsym, s, sr, ftau

    use mp_global,            ONLY : nproc_pool, me_pool, nproc_bgrp, me_bgrp, &
                                     init_index_over_band, inter_bgrp_comm, &
                                     mpime, inter_pool_comm, my_bgrp_id, &
                                     ibnd_start, ibnd_end, world_comm, nbgrp
    use mp,                   ONLY : mp_sum, mp_barrier, mp_bcast
    use funct,                ONLY : get_exx_fraction, start_exx, exx_is_active, &
                                     get_screening_parameter 
    USE fft_base,             ONLY : cgather_smooth, cscatter_smooth,&
         & dffts, cgather_custom, cscatter_custom
    use fft_interfaces,       ONLY : invfft

    implicit none
    integer :: ik,ibnd, i, j, k, ir, ri, rj, rk, isym, ikq
    integer :: h_ibnd, half_nbnd
    integer :: ipol, jpol
    COMPLEX(DP),allocatable :: temppsic(:), psic(:), tempevc(:,:)
    INTEGER :: nxxs, nrxxs, nr1x,nr2x,nr3x,nr1,nr2,nr3
#ifdef __MPI
    COMPLEX(DP),allocatable :: temppsic_all(:), psic_all(:)
#endif
    COMPLEX(DP) :: d_spin(2,2,48)
    INTEGER :: current_ik
    logical, allocatable :: present(:)
    logical :: exst
    INTEGER, ALLOCATABLE :: rir(:,:)
    
    COMPLEX(kind=DP), ALLOCATABLE :: state_fc_t(:,:),evc_g(:)
    integer       :: find_current_k


    call start_clock ('exxinit')

    !
    !  prepare the symmetry matrices for the spin part
    !
    ! Beware: not the same as nrxxs in parallel case
    nxxs = dffts%nr1x * dffts%nr2x * dffts%nr3x
    nrxxs= dffts%nnr
    nr1 = dffts%nr1
    nr2 = dffts%nr2
    nr3 = dffts%nr3
    nr1x = dffts%nr1x
    nr2x = dffts%nr2x
    nr3x = dffts%nr3x

#ifdef __MPI
       ALLOCATE(psic_all(nxxs), temppsic_all(nxxs) )
#endif

! JRD this creates ibnd_start and ibnd_end
    CALL init_index_over_band(inter_bgrp_comm,nbnd)

    ALLOCATE(temppsic(nrxxs), psic(nrxxs))

    if( .not. allocated( exxbuff ) ) allocate( exxbuff( nrxxs, nkqs, nbnd ) )

    allocate(present(nsym),rir(nxxs,nsym))
    allocate(tempevc( npwx, nbnd ))

    exx_nwordwfc=2*nrxxs
    if (.not.exx_is_active()) then 
       erfc_scrlen = get_screening_parameter()
       exxdiv = exx_divergence() 
       exxalfa = get_exx_fraction()
#ifdef EXXDEBUG
       write (stdout,*) " ! EXXALFA SET TO ", exxalfa
write(stdout,*) "exxinit, erfc_scrlen set to: ", erfc_scrlen
write(stdout,*) "exxinit, yukawa set to: ", yukawa
#endif
       call start_exx
    endif

#ifdef __MPI
    IF (pool_para) THEN
       IF ( .NOT.ALLOCATED (wg_collect) ) ALLOCATE(wg_collect(nbnd,nkstot))
       CALL wg_all(wg_collect, wg, nkstot, nks)
    ENDIF
#endif

    IF ( nks > 1 ) REWIND( iunigk )

    present(1:nsym) = .false.
    do ikq =1,nkqs
       isym = abs(index_sym(ikq))
       if (.not. present(isym) ) then
          present(isym) = .true.
          if ( mod (s (2, 1, isym) * nr1, nr2) .ne.0 .or. &
               mod (s (3, 1, isym) * nr1, nr3) .ne.0 .or. &
               mod (s (1, 2, isym) * nr2, nr1) .ne.0 .or. &
               mod (s (3, 2, isym) * nr2, nr3) .ne.0 .or. &
               mod (s (1, 3, isym) * nr3, nr1) .ne.0 .or. &
               mod (s (2, 3, isym) * nr3, nr2) .ne.0 ) then
             call errore ('exxinit',' EXX + smooth grid is not working',isym)
          end if
          do ir=1, nxxs
             rir(ir,isym) = ir
          end do
          do k = 1, nr3
             do j = 1, nr2
                do i = 1, nr1
                   call ruotaijk (s(1,1,isym), ftau(1,isym), i, j, k, &
                                  nr1,nr2,nr3, ri, rj , rk )
                   ir =   i + ( j-1)*nr1x + ( k-1)*nr1x*nr2x
                   rir(ir,isym) = ri + (rj-1)*nr1x + (rk-1)*nr1x*nr2x
                end do
             end do
          end do

       end if
    end do

    exxbuff=(0.0_DP,0.0_DP)
    ! set appropriately the x_occupation
    do ik =1,nkstot
       IF (pool_para) THEN
          x_occupation(1:nbnd,ik) = wg_collect (1:nbnd, ik) / wk_collect(ik)
       ELSE
          x_occupation(1:nbnd,ik) = wg(1:nbnd, ik) / wk(ik)
       ENDIF
    end do

!   This is parallelized over pool. Each pool computes only its k-points
    DO ik = 1, nks
       npw = ngk (ik)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          CALL get_buffer (tempevc, nwordwfc, iunwfc, ik)
       ELSE
          tempevc(1:npwx,1:nbnd) = evc(1:npwx,1:nbnd)
       ENDIF
       IF (pool_para) THEN
          current_ik=find_current_k(ik, nkstot, nks)
       ELSE
          current_ik=ik
       ENDIF

! Potentially too much work and memory. 
! May want to thread this loop
       do ibnd =1, nbnd     
             temppsic(:) = ( 0.D0, 0.D0 )
             temppsic(nls(igk(1:npw))) = tempevc(1:npw,ibnd)
             CALL invfft ('Wave', temppsic, dffts)

             do ikq=1,nkqs
                if (index_xk(ikq) .ne. current_ik) cycle

                isym = abs(index_sym(ikq) )
#ifdef __MPI
                call cgather_smooth(temppsic,temppsic_all)
                IF ( me_bgrp == 0 ) &
                    psic_all(1:nxxs) = temppsic_all(rir(1:nxxs,isym))
                call cscatter_smooth(psic_all,psic)
#else
                psic(1:nrxxs) = temppsic(rir(1:nrxxs,isym))
#endif
                if (index_sym(ikq) < 0 ) psic(1:nrxxs) = CONJG(psic(1:nrxxs))
                do ji=1, nrxxs
! JRD Reorder Exxbuf for below?
                   exxbuff(ji,ikq,ibnd)=psic(ji)
                enddo
             end do
       end do

    end do

! JRD Below required for reproducibility
    if(my_bgrp_id>0) then
       exxbuff=(0.0_DP,0.0_DP)
    endif
    if (nbgrp>1) then
       CALL mp_bcast(exxbuff,0,inter_bgrp_comm)
    endif

!   All pools have the complete set of wavefunctions
    IF (pool_para) THEN
          CALL mp_sum(exxbuff, inter_pool_comm)
    END IF

    deallocate(tempevc)
    deallocate(present,rir)
       deallocate(temppsic, psic)
#ifdef __MPI
       deallocate(temppsic_all, psic_all)
#endif 

    call stop_clock ('exxinit')  

end subroutine exxinit
  
!-----------------------------------------------------------------------
SUBROUTINE vexx(lda, n, m, psi, hpsi)
!-----------------------------------------------------------------------

    ! This routine calculates V_xx \Psi
    
    ! ... This routine computes the product of the Hamiltonian
    ! ... matrix with m wavefunctions contained in psi
    !
    ! ... input:
    ! ...    lda   leading dimension of arrays psi, spsi, hpsi
    ! ...    n     true dimension of psi, spsi, hpsi
    ! ...    m     number of states psi
    ! ...    psi
    !
    ! ... output:
    ! ...    hpsi  Vexx*psi
    !
    USE constants, ONLY : fpi, e2, pi
    USE cell_base, ONLY : alat, omega, bg, at, tpiba
    USE symm_base, ONLY : nsym, s
    USE gvect,     ONLY : ngm
    USE gvecs,   ONLY : nls, nlsm, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, current_k
    USE control_flags, ONLY : gamma_only
    USE klist,     ONLY : xk, nks, nkstot
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE gvect,     ONLY : g, nl
    use fft_base,  ONLY : dffts
    use fft_interfaces, ONLY : fwfft, invfft

    USE parallel_include  
    USE mp_global, ONLY : nproc, inter_image_comm, me_image, my_image_id,&
         & nimage, nproc_image, ibnd_start, ibnd_end, mpime, inter_bgrp_comm, intra_bgrp_comm,&
         & my_bgrp_id, nbgrp, me_bgrp, world_comm
    USE mp,        ONLY : mp_sum, mp_barrier, mp_bcast
    USE gvect,        ONLY : ecutrho
    USE wavefunctions_module, ONLY : psic

    IMPLICIT NONE

    INTEGER                     :: lda, n, m
    COMPLEX(DP)                 :: psi(lda,m) 
    COMPLEX(DP)                 :: hpsi(lda,m)

    INTEGER          :: nqi, myrank, mysize

    ! local variables
    COMPLEX(DP), allocatable :: tempphic(:), temppsic(:), result(:)
    COMPLEX(DP),ALLOCATABLE :: gresult(:)

    COMPLEX(DP), allocatable :: rhoc(:), vc(:)
    REAL (DP),   ALLOCATABLE :: fac(:), facb(:)
    INTEGER          :: ibnd, ik, im , ig, ikq, iq, isym, iqi, ipol
    INTEGER          :: h_ibnd, half_nbnd, ierr, nrxxs
    INTEGER          :: current_ik
    REAL(DP) :: x1, x2, dtmp
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), x, q(3)
    ! <LMS> temp array for vcut_spheric
    REAL(DP) :: atws(3,3)
    integer       :: find_current_k, ijk

    COMPLEX(kind=DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER, ALLOCATABLE :: igkt(:)

    CALL start_clock ('vexx')

    ALLOCATE( fac(ngm) )
    nrxxs = dffts%nnr
    ALLOCATE( facb(nrxxs) )

    ALLOCATE ( tempphic(nrxxs), temppsic(nrxxs), result(nrxxs) )
    ALLOCATE( gresult(npw) )

    ALLOCATE (rhoc(nrxxs), vc(nrxxs))

    ! write (*,*) exx_nwordwfc,lda,n,m, lda*n

    nqi=nqs

    IF (pool_para) THEN
       current_ik=find_current_k(current_k,nkstot,nks)
    ELSE
       current_ik = current_k
    ENDIF

    if(my_bgrp_id>0) then
       hpsi=(0.0_DP,0.0_DP)
       psi=(0.0_DP,0.0_DP)
    endif
    if (nbgrp>1) then
       CALL mp_bcast(hpsi,0,inter_bgrp_comm)
       CALL mp_bcast(psi,0,inter_bgrp_comm)
    endif

    DO im=1,m !for each band of psi (the k cycle is outside band)
       temppsic(:) = ( 0.D0, 0.D0 )

       temppsic(nls(igk(1:npw))) = psi(1:npw,im)
       CALL invfft ('Wave', temppsic, dffts)

       result(:)   = (0.d0,0.d0)

       DO iqi=1,nqi
          iq=iqi

          ikq  = index_xkq(current_ik,iq)
          ik   = index_xk(ikq)
          isym = ABS(index_sym(ikq))

          IF (pool_para) THEN
             xk_cryst(:) = at(1,:)*xk_collect(1,ik) + at(2,:)*xk_collect(2,ik)&
                         + at(3,:)*xk_collect(3,ik)
          ELSE
             xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) &
                         + at(3,:)*xk(3,ik)
          ENDIF
          IF (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
          sxk(:) = s(:,1,isym)*xk_cryst(1) + &
               s(:,2,isym)*xk_cryst(2) + &
               s(:,3,isym)*xk_cryst(3) 
          xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

          ! calculate the 1/|r-r'| factor and place it in fac
          CALL g2_convolution(ngm, g, xk(:,current_k), xkq, fac)

          facb = 0D0
          do ijk = 1, ngm
            !if (nls(ijk) .le. 0 .or. nls(ijk) .gt. nrxxs) write(*,*) "Ooops, outside bounds", nls(ijk), nrxxs
            facb(nls(ijk)) = fac(ijk)
          enddo

          DO ibnd=ibnd_start,ibnd_end !for each band of psi
              !IF ( ABS(x_occupation(ibnd,ik) - 1D0) > 1.d-6) write(*,*) 'WARNING',  my_bgrp_id, me_bgrp, ibnd, ik, x_occupation(ibnd,ik)
              IF ( ABS(x_occupation(ibnd,ik)) < 1.d-6) CYCLE

              !loads the phi from file

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ijk) SCHEDULE(STATIC)
              do ijk = 1, nrxxs
                tempphic(ijk)=exxbuff(ijk,ikq,ibnd)
                rhoc(ijk)=CONJG(tempphic(ijk))*temppsic(ijk) / omega
              enddo
!$OMP END PARALLEL DO
                
              !brings it to G-space
              CALL fwfft ('Smooth', rhoc, dffts)

              dtmp = x_occupation(ibnd,ik) / nqs
   
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ijk) SCHEDULE(STATIC)
              do ijk = 1, nrxxs
                !vc(nls(1:ngm)) = fac(1:ngm) * rhoc(nls(1:ngm)) * dtmp
                vc(ijk) = 0D0
                vc(ijk) = facb(ijk) * rhoc(ijk) * dtmp
              enddo
!$OMP END PARALLEL DO

              !brings back v in real space
              CALL invfft ('Smooth', vc, dffts) 
                
              !accumulates over inner bands and k points
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ijk) SCHEDULE(STATIC)
              do ijk = 1, nrxxs
                result(ijk)=result(ijk)+vc(ijk)*tempphic(ijk)
              enddo
!$OMP END PARALLEL DO
          END DO

       END DO

       !brings back result in G-space
       CALL fwfft ('Wave', result, dffts)

       gresult = 0D0
       DO ig = 1, npw
         gresult(ig) = result(nls(igk(ig)))
       ENDDO

       CALL mp_sum( gresult(1:npw), inter_bgrp_comm)

       DO ig = 1, npw
         hpsi(ig,im)=hpsi(ig,im) - exxalfa*gresult(ig)
       ENDDO

    END DO
    
    DEALLOCATE (tempphic, temppsic, result, gresult) 

    DEALLOCATE (rhoc, vc, fac, facb )

    CALL stop_clock ('vexx')

END SUBROUTINE vexx

!-----------------------------------------------------------------------
SUBROUTINE g2_convolution(ngm, g, xk, xkq, fac)
!-----------------------------------------------------------------------
  ! This routine calculates the 1/|r-r'| part of the exact exchange 
  ! expression in reciprocal space (the G^-2 factor).
  USE kinds,     ONLY : DP
  USE cell_base, ONLY : tpiba, at
  USE constants, ONLY : fpi, e2, pi
  
  IMPLICIT NONE
  
  INTEGER,  INTENT(IN)    :: ngm   ! Number of G vectors
  REAL(DP), INTENT(IN)    :: g(3,ngm) ! Cartesian components of G vectors
  REAL(DP), INTENT(IN)    :: xk(3) ! current k vector
  REAL(DP), INTENT(IN)    :: xkq(3) ! current q vector
  
  REAL(DP), INTENT(INOUT) :: fac(ngm) ! Calculated convolution
  
  
  !Local variables
  INTEGER :: ig !Counters 
  REAL(DP) :: q(3), qq, x
  
  !CALL start_clock ('vexx_ngmloop')
  DO ig=1,ngm
     !
     q(:)= xk(:) - xkq(:) + g(:,ig)
     !
     q = q * tpiba
     !
     qq = SUM(q(:)**2.0_DP) 
     !
     IF (x_gamma_extrapolation) THEN
        on_double_grid = .TRUE.
        x= 0.5d0/tpiba*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
        on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
        x= 0.5d0/tpiba*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
        on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
        x= 0.5d0/tpiba*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
        on_double_grid = on_double_grid .AND. (ABS(x-NINT(x))<eps)
     ENDIF
     
     IF ( use_coulomb_vcut_ws ) THEN
        !
        fac(ig) = vcut_get(vcut,q)
        !
     ELSE IF ( use_coulomb_vcut_spheric ) THEN
        !
        fac(ig) = vcut_spheric_get(vcut,q)
        !
     ELSE IF (qq.GT.1.d-8) THEN
        !
        IF ( erfc_scrlen > 0  ) THEN
           fac(ig)=e2*fpi/qq*(1.d0-EXP(-qq/4.d0/erfc_scrlen**2)) * grid_factor
        ELSEIF( erf_scrlen > 0 ) THEN
           fac(ig)=e2*fpi/qq*(EXP(-qq/4.d0/erf_scrlen**2)) * grid_factor
        ELSE
           fac(ig)=e2*fpi/( qq + yukawa ) * grid_factor
        END IF
        IF (on_double_grid) fac(ig) = 0.d0
        !
     ELSE
        !
        fac(ig)= - exxdiv ! or rather something else (see F.Gygi)
        !
        IF ( yukawa > 0.d0.AND. .NOT. x_gamma_extrapolation ) &
             fac(ig) = fac(ig) + e2*fpi/( qq + yukawa )
        IF( erfc_scrlen > 0.d0.AND. .NOT. x_gamma_extrapolation ) fac(ig) = fac(ig) + e2*pi/(erfc_scrlen**2)
        !
     ENDIF
     !
  ENDDO
  !CALL stop_clock ('vexx_ngmloop')
END SUBROUTINE g2_convolution

!-----------------------------------------------------------------------
function exxenergy2()
!-----------------------------------------------------------------------

    USE constants, ONLY : fpi, e2, pi
    USE io_files,  ONLY : iunigk,iunwfc, nwordwfc
    USE buffers,   ONLY : get_buffer
    USE cell_base, ONLY : alat, omega, bg, at, tpiba
    USE symm_base,ONLY : nsym, s
    USE gvect,     ONLY : ngm, gstart
    USE gvecs,   ONLY : nls, nlsm, doublegrid
    USE wvfct,     ONLY : nbnd, npwx, npw, igk, wg, current_k
    USE control_flags, ONLY : gamma_only
    USE wavefunctions_module, ONLY : evc
    USE klist,     ONLY : xk, ngk, nks, nkstot
    USE lsda_mod,  ONLY : lsda, current_spin, isk
    USE gvect,     ONLY : g, nl
    USE mp_global, ONLY : inter_pool_comm, inter_image_comm, inter_bgrp_comm, intra_bgrp_comm, nbgrp
    USE mp_global, ONLY : my_image_id, nimage, ibnd_start, ibnd_end
    USE mp,        ONLY : mp_sum
    use fft_base,  ONLY : dffts
    use fft_interfaces, ONLY : fwfft, invfft
    USE gvect,     ONLY : ecutrho

    IMPLICIT NONE
    REAL (DP)   :: exxenergy2,  energy

    ! local variables
    COMPLEX(DP), allocatable :: tempphic(:), temppsic(:)
    COMPLEX(DP), ALLOCATABLE :: rhoc(:)
    REAL (DP),   ALLOCATABLE :: fac(:)
    integer          :: jbnd, ibnd, ik, ikk, ig, ikq, iq, isym
    integer          :: half_nbnd, h_ibnd, nqi, iqi, nrxxs, ijk
    real(DP)    :: x1, x2
    REAL(DP) :: qq, xk_cryst(3), sxk(3), xkq(3), vc, x, q(3)
    ! temp array for vcut_spheric
    real(DP) :: atws(3,3) 
    integer       :: find_current_k

    COMPLEX(kind=DP), ALLOCATABLE :: psi_t(:), prod_tot(:)
    INTEGER, ALLOCATABLE :: igkt(:)
    INTEGER :: ngm_fft

    call start_clock ('exxen2')

    nrxxs = dffts%nnr
    ALLOCATE( fac(ngm) )

    ALLOCATE (tempphic(nrxxs), temppsic(nrxxs)) 
    ALLOCATE ( rhoc(nrxxs) )

    energy=0.d0

    nqi=nqs

    IF ( nks > 1 ) REWIND( iunigk )
    do ikk=1,nks
       IF (pool_para) THEN
          current_k=find_current_k(ikk,nkstot,nks)
       ELSE
          current_k = ikk
       ENDIF
       IF ( lsda ) current_spin = isk(ikk)
       npw = ngk (ikk)
       IF ( nks > 1 ) THEN
          READ( iunigk ) igk
          call get_buffer (evc, nwordwfc, iunwfc, ikk)
       END IF
       do jbnd=ibnd_start, ibnd_end !for each band of psi (the k cycle is outside band)
          temppsic(:) = ( 0.D0, 0.D0 )
          temppsic(nls(igk(1:npw))) = evc(1:npw,jbnd)
          CALL invfft ('Wave', temppsic, dffts)
       
          do iqi=1,nqi
             iq=iqi

             ikq  = index_xkq(current_k,iq)
             ik   = index_xk(ikq)
             isym = abs(index_sym(ikq))

             IF (pool_para) THEN
                xk_cryst(:) = at(1,:)*xk_collect(1,ik) + &
                              at(2,:)*xk_collect(2,ik) + &
                              at(3,:)*xk_collect(3,ik)
             ELSE
                xk_cryst(:) = at(1,:)*xk(1,ik) + at(2,:)*xk(2,ik) + &
                              at(3,:)*xk(3,ik)
             ENDIF

             if (index_sym(ikq) < 0 ) xk_cryst = - xk_cryst
             sxk(:) = s(:,1,isym)*xk_cryst(1) + &
                      s(:,2,isym)*xk_cryst(2) + &
                      s(:,3,isym)*xk_cryst(3) 
             xkq(:) = bg(:,1)*sxk(1) + bg(:,2)*sxk(2) + bg(:,3)*sxk(3)

             CALL g2_convolution(ngm, g, xk(:,current_k), xkq, fac)

                do ibnd=1,nbnd !for each band of psi
                   if ( abs(x_occupation(ibnd,ik)) < 1.d-6) cycle
                   !
                   !loads the phi from file
                   !

!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ijk) SCHEDULE(STATIC)
                   do ijk = 1, nrxxs
                        tempphic(ijk)=exxbuff(ijk,ikq,ibnd)

                        !calculate rho in real space
                        rhoc(ijk)=CONJG(tempphic(ijk))*temppsic(ijk) / omega
                   enddo
!$OMP END PARALLEL DO

                   !brings it to G-space
                   CALL fwfft ('Smooth', rhoc, dffts)

                   vc = 0.D0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ig) SCHEDULE(STATIC) REDUCTION(+:vc)
                   do ig=1,ngm
                      vc = vc + fac(ig) * rhoc(nls(ig)) * CONJG(rhoc(nls(ig)))
                   end do
!$OMP END PARALLEL DO
   
                   vc = vc * omega * x_occupation(ibnd,ik) / nqs

                   energy = energy - exxalfa * vc * wg(jbnd,ikk)
                end do

          end do
       end do
    end do

    DEALLOCATE (tempphic, temppsic) 

    DEALLOCATE (rhoc, fac )

    call mp_sum( energy, inter_bgrp_comm )
    call mp_sum( energy, intra_bgrp_comm )
    call mp_sum( energy, inter_pool_comm )

    exxenergy2 = energy

    call stop_clock ('exxen2')

end function  exxenergy2

!-----------------------------------------
function exx_divergence ()
!-----------------------------------------

     USE constants, ONLY : fpi, e2, pi
     USE cell_base, ONLY : bg, at, alat, omega
     USE gvect,     ONLY : ngm, g
     USE wvfct,     ONLY : ecutwfc
     USE io_global, ONLY : stdout
     USE control_flags, ONLY : gamma_only
     USE mp_global, ONLY : intra_bgrp_comm
     USE mp,        ONLY : mp_sum

     implicit none
     real(DP) :: exx_divergence

     ! local variables
     integer :: iq1,iq2,iq3, ig
     real(DP) :: div, dq1, dq2, dq3, xq(3), q_, qq, tpiba2, alpha, x, q(3)

     integer :: nqq, iq
     real(DP) :: aa, dq

     call start_clock ('exx_div')

     tpiba2 = (fpi / 2.d0 / alat) **2

     alpha  = 10.d0 * tpiba2 / ecutwfc

     IF ( .NOT. use_regularization ) THEN
        exx_divergence = 0.d0
        return
     END IF

     dq1= 1.d0/DBLE(nq1)
     dq2= 1.d0/DBLE(nq2) 
     dq3= 1.d0/DBLE(nq3) 

     div = 0.d0
     do iq1=1,nq1
        do iq2=1,nq2
           do iq3=1,nq3
              xq(:) = bg(:,1) * (iq1-1) * dq1 + &
                      bg(:,2) * (iq2-1) * dq2 + &
                      bg(:,3) * (iq3-1) * dq3 
              do ig=1,ngm
                 q(1)= xq(1) + g(1,ig)
                 q(2)= xq(2) + g(2,ig)
                 q(3)= xq(3) + g(3,ig)
                 qq = ( q(1)*q(1) + q(2)*q(2) + q(3)*q(3) ) 
                 if (x_gamma_extrapolation) then
                    on_double_grid = .true.
                    x= 0.5d0*(q(1)*at(1,1)+q(2)*at(2,1)+q(3)*at(3,1))*nq1
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0*(q(1)*at(1,2)+q(2)*at(2,2)+q(3)*at(3,2))*nq2
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                    x= 0.5d0*(q(1)*at(1,3)+q(2)*at(2,3)+q(3)*at(3,3))*nq3
                    on_double_grid = on_double_grid .and. (abs(x-nint(x))<eps)
                 end if
                 if (.not.on_double_grid) then
                    if ( qq > 1.d-8 ) then
                       if ( erfc_scrlen > 0 ) then
                          div = div + exp( -alpha * qq) / qq * &
                                (1.d0-exp(-qq*tpiba2/4.d0/erfc_scrlen**2)) * grid_factor
                       elseif ( erf_scrlen >0 ) then
                          div = div + exp( -alpha * qq) / qq * &
                                (exp(-qq*tpiba2/4.d0/erf_scrlen**2)) * grid_factor
                       else

                          div = div + exp( -alpha * qq) / (qq + yukawa/tpiba2) &
                                                     * grid_factor
                       endif
                    end if
                 end if
              end do
           end do
        end do
     end do
     call mp_sum(  div, intra_bgrp_comm )

     if ( .not. x_gamma_extrapolation ) then
        if ( yukawa > 0.d0) then
           div = div + tpiba2/yukawa
        elseif( erfc_scrlen > 0.d0 ) then
           div = div + tpiba2/4.d0/erfc_scrlen**2
        else
           div = div - alpha
        end if
     end if

     div = div * e2 * fpi / tpiba2 / nqs

     alpha = alpha / tpiba2

     nqq = 100000
     dq = 5.0d0 / sqrt(alpha) /nqq
     aa = 0.d0
     do iq=0,  nqq
        q_ = dq * (iq+0.5d0)
        qq = q_ * q_
        if ( erfc_scrlen > 0 ) then
           aa = aa  -exp( -alpha * qq) * exp(-qq/4.d0/erfc_scrlen**2) * dq
        elseif ( erf_scrlen > 0 ) then
           aa = 0.d0
        else
           aa = aa - exp( -alpha * qq) * yukawa / (qq + yukawa) * dq
        end if
     end do
     aa = aa * 8.d0 /fpi
     aa = aa + 1.d0/sqrt(alpha*0.25d0*fpi) 
     if( erf_scrlen > 0) aa = 1.d0/sqrt((alpha+1.d0/4.d0/erf_scrlen**2)*0.25d0*fpi)
#ifdef EXXDEBUG
     write (stdout,*) aa, 1.d0/sqrt(alpha*0.25d0*fpi)
#endif   
     div = div - e2*omega * aa

     exx_divergence = div * nqs
#ifdef EXXDEBUG
     write (stdout,'(a,i4,a,3f12.4)') 'EXX divergence (',nq1,')= ', &
                                  div, alpha
#endif
     call stop_clock ('exx_div')

     return
  end function exx_divergence 
  
end module exx
