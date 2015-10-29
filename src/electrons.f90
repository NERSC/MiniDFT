!
! Copyright (C) 2001-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
SUBROUTINE electrons()
  !----------------------------------------------------------------------------
  !
  ! ... This routine is a driver of the self-consistent cycle.
  ! ... It uses the routine c_bands for computing the bands at fixed
  ! ... Hamiltonian, the routine sum_band to compute the charge density,
  ! ... the routine v_of_rho to compute the new potential and the routine
  ! ... mix_rho to mix input and output charge densities.
  ! ... It prints on output the total energy and its decomposition in
  ! ... the separate contributions.
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : eps8, pi
  USE io_global,            ONLY : stdout, ionode
  USE cell_base,            ONLY : at, bg, alat, omega, tpiba2
  USE ions_base,            ONLY : zv, nat, nsp, ityp, tau, compute_eextfor
  USE basis,                ONLY : starting_pot, starting_wfc
  USE fft_base,             ONLY : dfftp
  USE gvect,                ONLY : ngm, gstart, nl, nlm, g, gg, gcutm
  USE gvecs,                ONLY : doublegrid, ngms
  USE klist,                ONLY : xk, wk, nelec, ngk, nks, nkstot, lgauss
  USE lsda_mod,             ONLY : lsda, nspin, magtot, absmag, isk
  USE vlocal,               ONLY : strf
  USE wvfct,                ONLY : nbnd, et, npwx, ecutwfc
  USE ener,                 ONLY : etot, hwf_energy, eband, deband, ehart, &
                                   vtxc, etxc, etxcc, ewld, demet, epaw
  USE scf,                  ONLY : scf_type, scf_type_COPY, &
                                   create_scf_type, destroy_scf_type, &
                                   rho, rho_core, rhog_core, &
                                   v, vltot, vrs, kedtau, vnew
  USE control_flags,        ONLY : mixing_beta, tr2, ethr, niter, nmix, &
                                   iprint, istep, lscf, lmd, conv_elec, &
                                   restart, io_level, do_makov_payne,  &
                                   iverbosity, textfor 
  USE io_files,             ONLY : iunwfc, iunocc, nwordwfc, output_drho, &
                                   iunefield, iunpaw
  USE buffers,              ONLY : save_buffer
  USE extfield,             ONLY : tefield, etotefield
  USE exx,                  ONLY : exxinit, exxenergy2, &
                                   fock0, fock1, fock2, dexx, exx_restart
  USE funct,                ONLY : dft_is_hybrid, exx_is_active
  USE wavefunctions_module, ONLY : evc, psic
  USE spin_orb,             ONLY : domag
  USE control_flags,        ONLY : adapt_thr, tr2_init, tr2_multi
  USE mp_global,            ONLY : intra_bgrp_comm, nbgrp, mpime, &
                                   inter_bgrp_comm, my_bgrp_id
  USE mp,                   ONLY : mp_sum
  !
  !
  USE uspp_param,           ONLY : nh, nhm ! used for PAW
  USE dfunct,                 only : newd
  !
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !
  real(dp) :: v_kin_r !subsitute dummy variable for v%kin_r
  real(dp) :: eth     !usbsitutue dummy variable for ldaU%eth
  REAL(DP) :: &
      dr2,          &! the norm of the diffence between potential
      charge,       &! the total charge
      deband_hwf,   &! deband for the Harris-Weinert-Foulkes functional
      mag           ! local magnetization
  INTEGER :: &
      i,            &! counter on polarization
      idum,         &! dummy counter on iterations
      iter,         &! counter on iterations
      ik_,          &! used to read ik from restart file
      kilobytes
  REAL(DP) :: &
      tr2_min,     &! estimated error on energy coming from diagonalization
      descf,       &! correction for variational energy
      en_el=0.0_DP,&! electric field contribution to the total energy
      eext=0.0_DP   ! external forces contribution to the total energy
  LOGICAL :: &
      first
  !
  ! ... auxiliary variables for calculating and storing temporary copies of
  ! ... the charge density and of the HXC-potential
  !
  type (scf_type), save :: rhoin ! used to store rho_in of current/next iteration
  !
  ! ... external functions
  !
  REAL(DP), EXTERNAL :: ewald, get_clock
  REAL(DP) :: tr2_final ! final threshold for exx minimization 
                        ! when using adaptive thresholds.
  iter = 0
  ik_  = 0
  dr2  = 0.0_dp
  tr2_final = tr2
  IF (dft_is_hybrid() .AND. adapt_thr ) THEN
     tr2= tr2_init
  ENDIF
  !
  WRITE( stdout, 9000 ) get_clock( 'MiniDFT' )
  !
  CALL memstat( kilobytes )
  !
  IF ( kilobytes > 0 ) WRITE( stdout, 9001 ) kilobytes/1000.0
  !
  CALL flush_unit( stdout )
  !
  !
  CALL start_clock( 'electrons' )
  !
  !  
  CALL flush_unit( stdout )
  !
  ! ... calculates the ewald contribution to total energy
  !
  ewld = ewald( alat, nat, nsp, ityp, zv, at, bg, tau, &
                omega, g, gg, ngm, gcutm, gstart, .false., strf )
  !
  call create_scf_type ( rhoin )
  !
10 CONTINUE

  !write(*,*) "About to start loop in electrons", dft_is_hybrid(), exx_is_active()

  ! ... Convergence threshold for iterative diagonalization

  ! ... for the first scf iteration of each ionic step (after the first),
  ! ... the threshold is fixed to a default value of 1.D-6

  IF ( istep > 0 ) ethr = 1.D-6

  WRITE( stdout, 9002 )

  CALL flush_unit( stdout )

  !%%%%%%%%%%%%%%%%%%%%          iterate !          %%%%%%%%%%%%%%%%%%%%%
  DO idum = 1, niter

     IF ( check_stop_now() ) RETURN

     iter = iter + 1

     WRITE( stdout, 9010 ) iter, ecutwfc, mixing_beta

     CALL flush_unit( stdout )

     ! ... Convergence threshold for iterative diagonalization is
     ! ... automatically updated during self consistency
     IF ( iter > 1 .AND. ik_ == 0 ) THEN
        IF ( iter == 2 ) ethr = 1.D-2
        ethr = MIN( ethr, 0.1D0*dr2 / MAX( 1.D0, nelec ) )
        ! ... do not allow convergence threshold to become too small:
        ! ... iterative diagonalization may become unstable
        ethr = MAX( ethr, 1.D-13 )
     END IF

     first = ( iter == 1 )

     ! ... deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v> is calculated a
     ! ... first time here using the input density and potential ( to be
     ! ... used to calculate the Harris-Weinert-Foulkes energy )
     deband_hwf = delta_e()

     ! save input current density in rhoin
     call scf_type_COPY( rho, rhoin )

     scf_step: DO

        ! ... tr2_min is set to an estimate of the error on the energy
        ! ... due to diagonalization - used only for the first scf iteration
        tr2_min = 0.D0

        IF ( first ) tr2_min = ethr*MAX( 1.D0, nelec ) 

        ! ... diagonalization of the KS hamiltonian
        CALL c_bands( iter, ik_, dr2 )

        IF ( check_stop_now() ) RETURN

        ! ... xk, wk, isk, et, wg are distributed across pools;
        ! ... the first node has a complete copy of xk, wk, isk,
        ! ... while eigenvalues et and weights wg must be
        ! ... explicitely collected to the first node
        ! ... this is done here for et, in sum_band for wg
        CALL poolrecover( et, nbnd, nkstot, nks )

        ! ... the new density is computed here
        CALL sum_band()

        ! ... the Harris-Weinert-Foulkes energy is computed here using only
        ! ... quantities obtained from the input density
        hwf_energy = eband + deband_hwf + (etxc - etxcc) + ewld + ehart + demet

        ! ... eband  = \sum_v \epsilon_v    is calculated by sum_band
        ! ... deband = - \sum_v <\psi_v | V_h + V_xc |\psi_v>
        ! ... eband + deband = \sum_v <\psi_v | T + Vion |\psi_v>
        deband = delta_e()

        ! ... mix_rho mixes several quantities: rho in g-space, tauk (for meta-gga)
        ! ... ns (for lda+u) and becsum (for paw)
        CALL mix_rho( rho, rhoin, mixing_beta, dr2, tr2_min, iter, nmix, conv_elec )

        ! ... if convergence is achieved or if the self-consistency error
        ! ... (dr2) is smaller than the estimated error due to diagonalization
        ! ... (tr2_min), rhoin and rho are unchanged: rhoin contains the input
        ! ...  density and rho contains the output density
        ! ... In the other cases rhoin contains the mixed charge density 
        ! ... (the new input density) while rho is unchanged
        IF ( first .and. nat > 0) THEN

           ! ... first scf iteration: check if the threshold on diagonalization
           ! ... (ethr) was small enough wrt the error in self-consistency (dr2)
           ! ... if not, perform a new diagonalization with reduced threshold
           first = .FALSE.

           IF ( dr2 < tr2_min ) THEN

              WRITE( stdout, '(/,5X,"Threshold (ethr) on eigenvalues was ", &
                               &    "too large:",/,5X,                      &
                               & "Diagonalizing with lowered threshold",/)' )

              ethr = 0.1D0*dr2 / MAX( 1.D0, nelec )

              CYCLE scf_step

           END IF
        END IF

        not_converged_electrons : &
        IF ( .NOT. conv_elec ) THEN

           ! ... no convergence yet: calculate new potential from mixed
           ! ... charge density (i.e. the new estimate)
           CALL v_of_rho( rhoin, rho_core, rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v)

           ! ... estimate correction needed to have variational energy:
           ! ... T + E_ion (eband + deband) are calculated in sum_band
           ! ... and delta_e using the output charge density rho;
           ! ... E_H (ehart) and E_xc (etxc) are calculated in v_of_rho
           ! ... above, using the mixed charge density rhoin%of_r.
           ! ... delta_escf corrects for this difference at first order
           descf = delta_escf()

           ! ... now copy the mixed charge density in R- and G-space in rho
           CALL scf_type_COPY( rhoin, rho )

        ELSE not_converged_electrons

           ! ... convergence reached:
           ! ... 1) the output HXC-potential is saved in vr
           ! ... 2) vnew contains V(out)-V(in) ( used to correct the forces ).
           vnew%of_r(:,:) = v%of_r(:,:)

           CALL v_of_rho( rho,rho_core,rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v)

           vnew%of_r(:,:) = v%of_r(:,:) - vnew%of_r(:,:)

           ! ... note that rho is here the output, not mixed, charge density
           ! ... so correction for variational energy is no longer needed
           descf = 0._dp

        END IF not_converged_electrons

        IF ( exx_is_active() ) THEN
           fock1 = exxenergy2()
           fock2 = fock0
        ELSE
           fock0 = 0.D0
           fock1 = 0.D0
           fock2 = 0.D0
        END IF

        ! ... if we didn't cycle before we can exit the do-loop
        EXIT scf_step

     END DO scf_step

     ! ... define the total local potential (external + scf)
     CALL set_vrs( vrs, vltot, v%of_r, kedtau, v_kin_r, dfftp%nnr, nspin, doublegrid )

     ! ... in the US case we have to recompute the self-consistent
     ! ... term in the nonlocal potential
     ! ... PAW: newd contains PAW updates of NL coefficients
     !write(*,*) "electrons.f90:480 this call to newd() might be moot if .not.okpaw"
     CALL newd()

     en_el=0.d0 

     WRITE( stdout, 9000 ) get_clock( 'MiniDFT' )

     IF ( conv_elec ) WRITE( stdout, 9101 )

     IF ( conv_elec .OR. MOD( iter, iprint ) == 0 ) THEN
        CALL print_ks_energies()
     END IF

     etot = eband + ( etxc - etxcc ) + ewld + ehart + deband + demet + descf +en_el
     IF( textfor ) THEN
        eext =  compute_eextfor()
        etot = etot + eext
     END IF

     etot = etot - 0.5D0*fock0
     hwf_energy = hwf_energy -0.5D0*fock0

     IF ( dft_is_hybrid() .AND. conv_elec ) THEN

        first = .NOT. exx_is_active()

        CALL exxinit()

        IF ( first ) THEN
           fock0 = exxenergy2()

           CALL v_of_rho( rho, rho_core,rhog_core, &
                          ehart, etxc, vtxc, eth, etotefield, charge, v)

           CALL set_vrs( vrs, vltot, v%of_r, kedtau, v_kin_r, dfftp%nnr, nspin, doublegrid )

           conv_elec = .false.
           iter = 0
           WRITE( stdout,'(5x,"EXX: now go back to refine exchange calculation")')
           WRITE( stdout, * ) fock0

           GO TO 10

        END IF

        fock2 = exxenergy2()

        dexx = fock1 - 0.5D0*( fock0 + fock2 )

        etot = etot  - dexx
        hwf_energy = hwf_energy - dexx

        WRITE( stdout, * ) fock0, fock1, fock2
        WRITE( stdout, 9066 ) dexx

        fock0 = fock2

     END IF

     IF ( tefield ) THEN
        etot = etot + etotefield
        hwf_energy = hwf_energy + etotefield
     END IF

     IF ( ( conv_elec .OR. MOD( iter, iprint ) == 0 ) .AND. .NOT. lmd ) THEN

        IF ( dr2 > eps8 ) THEN
           WRITE( stdout, 9081 ) etot, hwf_energy, dr2
        ELSE
           WRITE( stdout, 9083 ) etot, hwf_energy, dr2
        END IF
        IF ( dft_is_hybrid()) THEN
           WRITE( stdout, 9062 ) - fock1
           WRITE( stdout, 9064 ) 0.5D0*fock2
        ENDIF

        WRITE( stdout, 9060 ) &
            ( eband + deband ), ehart, ( etxc - etxcc ), ewld

        IF ( textfor)             WRITE( stdout, &
            '(/5x,"Energy of the external Forces = ", F18.8)' ) eext
        IF ( tefield )            WRITE( stdout, 9061 ) etotefield
        IF ( ABS (descf) > eps8 ) WRITE( stdout, 9069 ) descf

        ! ... With Fermi-Dirac population factor, etot is the electronic
        ! ... free energy F = E - TS , demet is the -TS contribution

        IF ( lgauss ) WRITE( stdout, 9070 ) demet

     ELSE IF ( conv_elec .AND. lmd ) THEN
        IF ( dr2 > eps8 ) THEN
           WRITE( stdout, 9081 ) etot, hwf_energy, dr2
        ELSE
           WRITE( stdout, 9083 ) etot, hwf_energy, dr2
        END IF
     ELSE
        IF ( dr2 > eps8 ) THEN
           WRITE( stdout, 9080 ) etot, hwf_energy, dr2
        ELSE
           WRITE( stdout, 9082 ) etot, hwf_energy, dr2
        END IF
     END IF

     IF ( lsda ) WRITE( stdout, 9017 ) magtot, absmag

     CALL flush_unit( stdout )

     IF ( conv_elec ) THEN

        IF ( dft_is_hybrid() .AND. dexx > tr2_final ) THEN

           conv_elec = .false.
           iter = 0

           WRITE (stdout,*) " NOW GO BACK TO REFINE HYBRID CALCULATION"

           IF ( adapt_thr ) THEN
              tr2 = MAX(tr2_multi * dexx, tr2_final)
              WRITE( stdout, 9121 ) tr2
           ENDIF

           GO TO 10
        END IF

        WRITE( stdout, 9110 ) iter

        ! ... jump to the end
        IF ( output_drho /= ' ' ) CALL remove_atomic_rho()

        CALL stop_clock( 'electrons' )

        call destroy_scf_type ( rhoin )

        RETURN

     END IF

  END DO

  WRITE( stdout, 9101 )
  WRITE( stdout, 9120 ) iter

  CALL flush_unit( stdout )

  IF ( output_drho /= ' ' ) CALL remove_atomic_rho()

  CALL stop_clock( 'electrons' )

  RETURN

9000 FORMAT(/'     total cpu time spent up to now is ',F10.1,' secs' )
9001 FORMAT(/'     per-process dynamical memory: ',f7.1,' Mb' )
9002 FORMAT(/'     Self-consistent Calculation' )
9010 FORMAT(/'     iteration #',I3,'     ecut=', F9.2,' Ry',5X,'beta=',F4.2 )
9017 FORMAT(/'     total magnetization       =', F9.2,' Bohr mag/cell', &
            /'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9018 FORMAT(/'     total magnetization       =',3F9.2,' Bohr mag/cell' &
       &   ,/'     absolute magnetization    =', F9.2,' Bohr mag/cell' )
9050 FORMAT(/'     WARNING: integrated charge=',F15.8,', expected=',F15.8 )
9060 FORMAT(/'     The total energy is the sum of the following terms:',/,&
            /'     one-electron contribution =',F17.8,' Ry' &
            /'     hartree contribution      =',F17.8,' Ry' &
            /'     xc contribution           =',F17.8,' Ry' &
            /'     ewald contribution        =',F17.8,' Ry' )
9061 FORMAT( '     electric field correction =',F17.8,' Ry' )
9062 FORMAT( '     - averaged Fock potential =',F17.8,' Ry' )
9064 FORMAT( '     + Fock energy             =',F17.8,' Ry' )
9065 FORMAT( '     Hubbard energy            =',F17.8,' Ry' )
9066 FORMAT( '     est. exchange err (dexx)  =',F17.8,' Ry' )
9069 FORMAT( '     scf correction            =',F17.8,' Ry' )
9070 FORMAT( '     smearing contrib. (-TS)   =',F17.8,' Ry' )
9071 FORMAT( '     Magnetic field            =',3F12.7,' Ry' )
9072 FORMAT( '     Magnetic field            =',F12.7, ' Ry' )
9073 FORMAT( '     lambda                    =',F11.2,' Ry' )
9074 FORMAT( '     Dispersion Correction     =',F17.8,' Ry' )
9080 FORMAT(/'     total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',0PF17.8,' Ry' )
9081 FORMAT(/'!    total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',0PF17.8,' Ry' )
9082 FORMAT(/'     total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',1PE17.1,' Ry' )
9083 FORMAT(/'!    total energy              =',0PF17.8,' Ry' &
            /'     Harris-Foulkes estimate   =',0PF17.8,' Ry' &
            /'     estimated scf accuracy    <',1PE17.1,' Ry' )
9085 FORMAT(/'     total all-electron energy =',0PF17.6,' Ry' )
9101 FORMAT(/'     End of self-consistent calculation' )
9110 FORMAT(/'     convergence has been achieved in ',i3,' iterations' )
9120 FORMAT(/'     convergence NOT achieved after ',i3,' iterations: stopping' )
9121 FORMAT(/'     scf convergence threshold =',1PE17.1,' Ry' )

  CONTAINS
     !
     !-----------------------------------------------------------------------
     !
     !-----------------------------------------------------------------------
     FUNCTION check_stop_now()
       !-----------------------------------------------------------------------
       !
       USE check_stop,    ONLY : global_check_stop_now => check_stop_now
       !
       IMPLICIT NONE
       !
       LOGICAL :: check_stop_now
       INTEGER :: unit
       !
       unit = stdout
       !
       check_stop_now = global_check_stop_now( unit )
       !
       IF ( check_stop_now ) conv_elec = .FALSE.
       !
       RETURN
       !
     END FUNCTION check_stop_now
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_e()
       !-----------------------------------------------------------------------
       ! ... delta_e = - \int rho%of_r(r)  v%of_r(r)
       !               - \int rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !               - \sum rho%ns       v%ns       [for LDA+U]
       !               - \sum becsum       D1_Hxc     [for PAW]
       IMPLICIT NONE
       REAL(DP) :: delta_e, delta_e_hub
       !
       delta_e = - SUM( rho%of_r(:,:)*v%of_r(:,:) )
       !
       !
       delta_e = omega * delta_e / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
       CALL mp_sum( delta_e, intra_bgrp_comm )
       !
       !
       !
       RETURN
       !
     END FUNCTION delta_e
     !
     !-----------------------------------------------------------------------
     FUNCTION delta_escf()
       !-----------------------------------------------------------------------
       !
       ! ... delta_escf = - \int \delta rho%of_r(r)  v%of_r(r)
       !                  - \int \delta rho%kin_r(r) v%kin_r(r) [for Meta-GGA]
       !                  - \sum \delta rho%ns       v%ns       [for LDA+U]
       !                  - \sum \delta becsum       D1         [for PAW]
       ! ... calculates the difference between the Hartree and XC energy
       ! ... at first order in the charge density difference \delta rho(r)
       IMPLICIT NONE
       !
       REAL(DP) :: delta_escf, delta_escf_hub
       !
       delta_escf = - SUM( ( rhoin%of_r(:,:)-rho%of_r(:,:) )*v%of_r(:,:) )
       !
       !
       delta_escf = omega * delta_escf / ( dfftp%nr1*dfftp%nr2*dfftp%nr3 )
       !
       CALL mp_sum( delta_escf, intra_bgrp_comm )
       !


       RETURN
       !
     END FUNCTION delta_escf
     !
     !-----------------------------------------------------------------------
     !
END SUBROUTINE electrons
