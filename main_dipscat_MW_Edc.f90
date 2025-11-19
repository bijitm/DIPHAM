      program dipolar_main
        use dipole_module, only: npair, ipair, nmax, nmonbasis, iqn, &
                           nqn, frot, ffld, lrstrct, evalref, itrns, &
                           etrns, fdlt, mxlam, nvlblk, nphmn, nphmx, &
                           lambda, jlevel, theta, phi, xi, npol, &
                           nconst, zero, npower
        use constants, only: MHz_to_invcm, fstark, pi
        implicit none
        integer :: ia, n, mn, ii, k, ia1, ia2, nph, np, mnp, iref
        integer :: itime_start, itime_end, istep
        real*8  :: brot, dipole, omega, delta, eflddc, efldac, &
                vcoup, vconst(nconst), a(100), p(100), rmin, rmax, &
                dr, r, veff, eref, omega_pi, delta_pi
        real*8  :: count_rate
        real*8, allocatable :: elevel(:), thrshvals(:), wmat(:,:), &
                eval(:), work(:)
        logical :: leffMW, leffdual
        integer :: idate_time(8)
        character*10 :: b(3)
        integer :: info

        ! functions
        real*8 :: c1func

        namelist/params/nmax, npair, ipair, nphmn, nphmx, npol, brot, &
                        dipole, efldac, delta, omega, xi, theta, phi, &
                        eflddc, mxlam, lambda, npower, a, rmin, rmax, &
                        dr, iref, leffMW, leffdual, itrns, omega_pi, &
                        delta_pi

        ! Get date and time
        call date_and_time (b(1),b(2),b(3),idate_time)
        write(*,'(i2,":",i2,":",i2,"  ",i2,"/",i2,"/",i4)') &
                idate_time(5:7),idate_time(3),idate_time(2),&
                idate_time(1)

        write(*,1)
    1   format(/' Program to calculate adiabats as a function of R',&
     &         ' and theta for dipolar scattering'/&
     &         ' including static electric and MW fields')

        write(*,2)
    2   format(/" Basis set is symmetrized product of dc field",&
     &          " dressed states of diatoms and photon states,"/& 
     &          " which in general reads |til{n}1 mn1>|til{n}2 mn2>|N>")

        ! Define default values of input parameters
        nmax     = 1      ! Molecule nmax
        npair    = 0      ! No. of restricted basis functions. 0 for full.
        if (.not. allocated(ipair)) allocate(ipair(1000))
        ipair    = -999   ! List of pair states: 3 integers for each pair state:
                          ! a1, a2, nph, separated by commas
        nphmn    = 0      ! Min photon quantum number for the Fock state
        nphmx    = 0      ! Max photon quantum number for the Fock state
        brot     = -999d0 ! Rotational constant in MHz
        dipole   = -999d0 ! Dipole moment in debye
        eflddc   = 0.d0   ! DC electric field in kV/cm
        efldac   = 0.d0   ! AC electric field in kV/cm for the MW field
        npol     = 0      ! Polarization: 1:RHC_x-y, 0:Li_z
        delta    = 0.d0   ! Detuning frequncy (MHz) of the MW field
        omega    = 0.d0   ! Rabi frequncy (MHz) of the MW field
        delta_pi = 0.d0   ! Detuning frequncy (MHz) of the 2nd linear z MW field
        omega_pi = 0.d0   ! Rabi frequncy (MHz) of the 2nd linear z MW field
        xi       = 0.d0   ! Ellipticity (degrees) of the MW field polarization
        theta    = 0.d0   ! Interaction anisotropy (degrees)
        phi      = 0.d0   ! Interaction anisotropy (degrees)
        itrns    = -999   ! ITRNS(2) array takes two integer elements that define
                          ! the two monomer dressed states for the MW transition
        mxlam    = 1      ! Number of R-dependent potential terms
        if (.not. allocated(lambda)) allocate(lambda(1000))
        lambda   = -999   ! Groups of three indices to describe diatom-diatom
                          ! potential terms
        npower   = 0      ! Array of mxlam elements: R^{npower} describes the 
                          ! potential scaling 
        a        = 0d0    ! Array of mxlam elements: Coefficient of the potential
                          ! term in cm-1
        rmin     = 50d0   ! all r in bohr
        rmax     = 1000d0
        dr       = 5d0    
        iref     = 0      ! Reference pair state which will be 0 in
                          ! energy, if zero then no reference level
        leffMW   = .false.! Logical flag for writing Deng et al potential for
                          ! single MW
        leffdual = .false.! Logical flag for writing Deng et al potential for
                          ! double MW

        read(5,params)
        nmonbasis = (nmax+1)**2

        ! Sanity checks on a few input params
        if (brot<0d0) stop "Rotational constant not supplied"
        if (dipole<0d0) stop "Dipole moment not supplied"
        if (.not.(npol==0.or.npol==1)) stop &
          " STOP: Incorrect npol param, only allowed values"// &
          " are 0 (Linear z) or 1 (RHC x-y plane)"
        if (xi<0d0.or.xi>45d0) stop "xi must be between 0 and 45."
        if (abs(xi)>zero.and.npol==0) stop &
        "xi is non-zero: cannot have z-polarization (npol=0)"
        if (npair>0.and.any(ipair(:npair*nqn)==-999)) &
     &    stop "incompatible inputs for npair and ipair"
        if (nphmx<nphmn) stop "nphmx < nphmn"

        ! Write out input params provided

        write(6,102) nmax
 102    format(/' Basis set uses monomer diatom nmax =',i3)

        write(6,103) eflddc
 103    format(/' Dressed basis set will be constructed at',&
     &          ' reference dc electric field =', f8.4, ' kV/cm')

        write(*,105) nmonbasis
 105    format(/' Number of monomer basis functions =',i4)

        write(*,109) theta, phi
 109    format(/' Theta and phi in degrees are ',f5.1," and ",f5.1)

        ! Important prefactors (in cm-1)
        frot = brot*MHz_to_invcm
        ffld = -eflddc*dipole*fstark  ! (kV/cm)*debye*(1/kV/debye)
        fdlt = delta*MHz_to_invcm

        ! Convert angles from degrees to radians
        theta = theta*pi/180d0
        phi   = phi*pi/180d0
        xi    = xi*pi/180d0

        ! Flag for restricted basis set
        lrstrct = .false.
        if (npair>0) lrstrct = .true.
 
        ! Defines the basis functions for monomers 
        ! Creates an array that stores individual monomer functions
        call monomer_mapping

        write(*,107)
 107    format(/" Assigning a non-negative integer quantum number",&
     &       " 'a' to represent each monomer state")

        write(*,106)
 106    format(/" Mapping of monomer quantum numbers:",&
     &         /" ----------------------------- ",&
     &         /"       a       n       mn      ",&
     &         /" ----------------------------- ")

        do ia = 1, nmonbasis
           n  = iqn(ia,1)
           mn = iqn(ia,2)
           write(*,'(5i8)') ia-1,n,mn
        enddo

        ! Creates the pair basis set
        ! Calculates the dc-field dressed energies (eigenvalues) and eigenvectors
        call pairbasis_builder

        ! If MW field is on
        if (efldac>zero .or. omega>zero) then
          ! First calculate angular part
          write(6,48) itrns(1), itrns(2)
  48      format(" MW transition from states a_i =",i4," to a_f =",i4)
          if(npol==1) write(6,'(/" Right-hand circularly polarized")')
          if(npol==0) write(6,'(/" Linear polarization along z" )')
          write(6,49) xi*180d0/pi
  49      format(/" Ellipticity xi in polarization (degrees) =",f6.1)
          n   = iqn(itrns(1)+1,1)
          mn  = iqn(itrns(1)+1,2)
          np  = iqn(itrns(2)+1,1)
          mnp = iqn(itrns(2)+1,2)
          if (npol==0) call vmfdr(n,mn,0,np,mnp,-1,0,vcoup)
          if (npol/=0) call vmfdr(n,mn,0,np,mnp,(mn-mnp)/npol,npol,vcoup)
          if (abs(vcoup)<zero) stop "No Rabi coupling. Check ITRNS array and/or NPOL"
          write(6,53) vcoup
  53      format(/" MW transition dipole moment (in units of mu):", es10.2)

          if (efldac>zero) then
            write(6,*)
            write(6,*)"Since EFLDAC is provided, OMEGA will be calculated based on it."
            omega = abs(efldac*dipole*fstark*vcoup)/MHz_to_invcm ! in MHz
          endif
          if (efldac<zero .and. omega>zero) &
                  efldac = abs(omega*MHz_to_invcm/dipole/fstark/vcoup)  ! in kV/cm
          write(6,50) omega
  50      format(/" Rabi frequency for the transition:", es10.2, " MHz")
          write(6,555) delta
 555      format(/" Detuning frequency for the transition:", es10.2, " MHz")
          write(6,540) efldac*1d3
 540      format(/' AC electric field for the MW field:', es10.2, ' V/cm')
          write(6,505) (etrns+fdlt)/brot/MHz_to_invcm ! in units of brot
 505      format(/" MW transition frequency (in units of brot):", es12.4)
        endif

        write(*,*)"-------------------------------"
        write(*,*)"Efield dressed monomer energies"
        write(*,*)"-------------------------------"
        write(*,*)"     a       Energy (brot)     "
        write(*,*)"-------------------------------"
        do ia = 1, nmonbasis
          write(*,'(i7, es18.4)')ia-1, evalref(ia)/brot/MHz_to_invcm
        enddo

        ! NPAIR and IPAIR get updated from BASIS_BUILDER
        if (.not. allocated(jlevel)) allocate(jlevel(npair*nqn))
        if (.not. allocated(elevel)) allocate(elevel(npair))
        jlevel = -999
        elevel = 0.d0

        write(*,*)
        write(*,*)"-------------------------------------------------------------"
        write(*,*)"       Pair state quantum numbers and diagonal energies      "
        write(*,*)"-------------------------------------------------------------"
        write(*,*)"    State      a1       a2      nph      Energy (MHz)        "
        write(*,*)"-------------------------------------------------------------"
        k = 0
        do ii = 1, npair
          ia1 = ipair(k+1)
          ia2 = ipair(k+2)
          nph = ipair(k+3)
          jlevel(nqn*(ii-1)+1) = ia1
          jlevel(nqn*(ii-1)+2) = ia2
          jlevel(nqn*(ii-1)+3) = nph
          k = k+nqn
          ! energies in cm-1
          elevel(ii) = evalref(ia1+1)+evalref(ia2+1) &
                  + (etrns+fdlt)*nph
          write(*,'(4i9,es18.6)')ii,ia1,ia2,nph, elevel(ii)/MHz_to_invcm
        enddo

        nvlblk = mxlam+nconst ! Total number of Hamiltonian blocks

        write(*,345) mxlam, nconst
 345    format(/" Potential terms will be constructed from ",&
     &       i2," R-dependent term(s) and ",i2," asymptotic term(s).")

        ! define vconst array in cm-1
        vconst(1)   = 1d0
        vconst(2:3) = -efldac*dipole*fstark/2d0  ! (kV/cm)*debye*(1/kV/debye)

        write(6,'(/" VCONST (cm-1) for MW field is ",es16.8)') &
                vconst(2)

        ! Create coupling matrix elements
        call coupling

        ! Create potential coeff array
        p = 0d0
        p(:mxlam) = a(:mxlam)
        p(mxlam+1:nvlblk) = vconst

        ! Calculate thresholds
        if (.not. allocated(thrshvals)) allocate(thrshvals(npair))
        call thresholds(nvlblk,p,npair,elevel,thrshvals)
        write (6,'(" The thresholds (MHz) are"/,50es14.4)') &
                thrshvals/MHz_to_invcm

        if (iref<0.or.iref>npair) stop "invalid IREF"
        ! Define reference energy
        eref = 0d0
        if (iref/=0) eref = thrshvals(iref)

        ! Allocate a few arrays before final diagonalization
        if (.not. allocated(wmat)) allocate(wmat(npair,npair))
        if (.not. allocated(eval)) allocate(eval(npair))
        if (.not. allocated(work)) allocate(work(3*npair-1))
        wmat = 0d0
        eval = -999d0
        work = 0d0

        call system_clock(itime_start,count_rate)

        open (10,file="calculated_adiabats.dat",status="unknown")
        write(10,('("# Adiabats:",/"# R (bohr)  theta (rad) Energies (MHz)",/ &
     &              "#-------------------------")'))
        
        if (leffMW) then
        ! Writes Deng et al effective potential
          open (20,file="Deng_et_al_veff.dat",status="unknown")
          write(20,('("# V_effective:",/"# R (bohr)  theta (rad) Energy (MHz)",/ &
     &                "#-------------------------")'))
        endif

        if (leffdual) then
        ! Writes Deng et al effective potential
          open (21,file="Deng_et_al_veff_dualMW.dat",status="unknown")
          write(21,('("# V_effective:",/"# R (bohr)  theta (rad) Energy (MHz)",/ &
     &                "#-------------------------")'))
        endif

        r = rmin
        istep = 0
        do while (r<=rmax)
          ! Get Hamiltonian at a given value of r-grid
          call hammat(r,nvlblk,p,npair,elevel,wmat) 
          call dsyev('N','U',npair,wmat,npair,eval,work(:3*npair-1),&
                  3*npair-1,info)
          write(10,'(1x,f10.2, f18.10,50es14.4)') r, theta, &
                  (eval-eref)/MHz_to_invcm
          if (leffMW) call microwave_eff_pot(r,theta,phi,xi,omega,delta,&
                  dipole,veff)
          if (leffMW) write(20,'(1x,f10.2, f14.6, es14.4)') r,theta,veff/MHz_to_invcm
          if (leffdual) call dual_MW_eff_pot(r,theta,omega,delta,omega_pi,delta_pi, &
                  dipole,veff)
          if (leffdual) write(21,'(1x,f10.2, f14.6, es14.4)') r,theta,veff/MHz_to_invcm
          r = r+dr
          istep = istep+1
        enddo

        write(6,"(' No. of R steps:',i10)") istep

        call system_clock(itime_end)
        write(6,'(/" Time (s) in R loop:",es12.2,/)') &
                dble(itime_end-itime_start)/count_rate

        write(6,*)"Adiabats written to file 'calculated_adiabats.dat'"
        if (leffMW) write(6,*)"Deng et al. effective potential written", &
     &          " in file 'Deng_et_al_veff.dat'"
        if (leffdual) write(6,*)"Deng et al. effective potential for dual MW written", &
     &          " in file 'Deng_et_al_veff_dualMW.dat'"

        write(6,'(/" Program run successfully")')

        ! Now deallocate arrays
        if (allocated(ipair))  deallocate(ipair)
        if (allocated(lambda)) deallocate(lambda)
        if (allocated(wmat))   deallocate(wmat)
        if (allocated(eval))   deallocate(eval)
        if (allocated(work))   deallocate(work)
        if (allocated(jlevel)) deallocate(jlevel)
        if (allocated(elevel)) deallocate(elevel)

      end program dipolar_main
