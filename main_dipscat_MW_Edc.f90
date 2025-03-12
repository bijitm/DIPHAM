      program dipolar_main
        use dipole_module, only: npair, ipair, nmax, nmonbasis, iqn, &
                           nqn, frot, ffld, lrstrct, nphmx, evalref, &
                           itrns, etrns, fdlt, fomg, mxlam, nvlblk, &
                           lambda, jlevel, theta, phi, xi, npol, &
                           nconst, zero, npower
        use constants, only: MHz_to_invcm, fstark, pi
        implicit none
        integer :: ia, n, mn, ii, k, ia1, ia2, nph, np, mnp, iref
        integer :: itime_start, itime_end, istep
        real*8  :: brot, dipole, omega, delta, eflddc, efldac, &
                vcoup, vconst(nconst), a(100), p(100), rmin, rmax, &
                dr, r, veff, eref
        real*8  :: count_rate
        real*8, allocatable :: elevel(:), thrshvals(:), wmat(:,:), &
                eval(:), work(:)
        logical :: ldeng
        integer :: idate_time(8)
        character*10 :: b(3)
        integer :: info

        namelist/params/nmax, npair, ipair, nphmx, npol, brot, &
                        dipole, omega, delta, xi, theta, phi, eflddc, &
                        mxlam, lambda, npower, a, rmin, rmax, dr, &
                        iref, ldeng, itrns

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
        nmax   = 1      ! Molecule nmax
        npair  = 0      ! No. of restricted basis functions. 0 for full.
        if (.not. allocated(ipair)) allocate(ipair(1000))
        ipair  = -999   ! List of pair states
        nphmx  = 0      ! Max photon quantum number
        brot   = -999d0 ! Rotational constant in MHz
        dipole = -999d0 ! Dipole moment in debye
        npol   = 0      ! Polarization: 1:RHC_x-y,0:Li_z
        omega  = 0.d0   ! Rabi frequency for n=0-1 transition (2pixMHz)
        delta  = 0.d0   ! MW detuning (2pixMHz)
        xi     = 0.d0   ! in units of degree
        theta  = 0.d0   ! in units of degree
        phi    = 0.d0   ! in units of degree
        eflddc = 0.d0   ! DC electric field in kV/cm
        itrns  = (/0,1/)! ITRNS array defines the two monomer dressed 
                        ! states for the MW transition
        mxlam  = 1      ! Number of R-dependent potential terms
        if (.not. allocated(lambda)) allocate(lambda(1000))
        lambda = -999   ! Groups of three indices to describe diatom-diatom
                        ! potential terms
        npower = 0      ! Array of mxlam elements: R^{npower} describes the 
                        ! potential scaling 
        a      = 0d0    ! Array of mxlam elements: Coefficient of the potential
                        ! term in cm-1
        rmin   = 50d0   ! all r in bohr
        rmax   = 1000d0
        dr     = 5d0    
        iref   = 0      ! Reference pair state which will be 0 in
                        ! energy, if zero then no reference level
        ldeng  = .false.! Logical flag for writing Deng et al potential

        read(5,params)
        nmonbasis = (nmax+1)**2

        ! Some checks on input params
        if (brot<0d0) stop "Rotational constant not supplied"
        if (dipole<0d0) stop "Dipole moment not supplied"
        if (.not.(npol==0.or.npol==1)) stop &
          " STOP: Incorrect npol param, only allowed values"// &
          " are 0 (Linear z) or 1 (RHC x-y plane)"

        if (xi<0d0.or.xi>45d0) stop "xi must be between 0 and 45."

        if (abs(xi)>zero.and.npol==0) stop &
        "xi is non-zero: cannot have z-polarization (npol=0)"

        ! Important prefactors (in cm-1)
        frot = brot*MHz_to_invcm
        ffld = -eflddc*dipole*fstark  ! (kV/cm)*debye*(1/kV/debye)
        fomg = omega*MHz_to_invcm
        fdlt = delta*MHz_to_invcm

        write(6,102) nmax
 102    format(/' Basis set uses monomer diatom nmax =',i3)

        write(6,103) eflddc
 103    format(/' Dressed basis set will be constructed at',&
     &          ' reference electric field =', f8.4, ' kV/cm')

        write(*,104) nmonbasis
 104    format(/' Number of monomer basis functions =',i4)

        write(*,109) theta, phi
 109    format(/' Theta and phi in degrees are ',f5.1," and ",f5.1)

        ! Convert angles from degrees to radians
        theta = theta*pi/180d0
        phi   = phi*pi/180d0
        xi    = xi*pi/180d0

        ! Flag for restricted basis set
        lrstrct = .false.
        if (npair>0) lrstrct = .true.
 
        ! Define basis functions for 1 molecule
        ! Create an array that stores individual monomer functions
        call monomer_mapping

        write(*,105)
 105    format(/" Assigning a non-negative integer quantum number",&
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

        ! Sanity check on IPAIR array
        if (npair>0.and.any(ipair(:npair*nqn)==-999)) &
     &    stop "incompatible inputs for npair and ipair"

        ! Create pair basis set
        call pairbasis_builder

        write(*,*)"-------------------------------"
        write(*,*)"Efield dressed monomer energies"
        write(*,*)"-------------------------------"
        write(*,*)"     a       Energy (brot)     "
        write(*,*)"-------------------------------"
        do ia = 1, nmonbasis
          write(*,'(i7, es18.4)')ia-1, evalref(ia)/brot/MHz_to_invcm
        enddo

        ! Check IREF input parameter
        if (iref<0.or.iref>npair) stop "invalid IREF"

        ! NPAIR and IPAIR get updated from BASIS_BUILDER
        if (.not. allocated(jlevel)) allocate(jlevel(npair*nqn))
        if (.not. allocated(elevel)) allocate(elevel(npair))
        jlevel = -999
        elevel = 0.d0

        write(*,*)
        write(*,*)"-----------------------------------------"
        write(*,*)"       Pair state quantum numbers        "
        write(*,*)"-----------------------------------------"
        write(*,*)"    State      a1       a2      nph      "
        write(*,*)"-----------------------------------------"
        k = 0
        do ii = 1, npair
          ia1 = ipair(k+1)
          ia2 = ipair(k+2)
          nph = ipair(k+3)
          jlevel(nqn*(ii-1)+1) = ia1
          jlevel(nqn*(ii-1)+2) = ia2
          jlevel(nqn*(ii-1)+3) = nph
          k = k+nqn
          write(*,'(4i9)')ii,ia1,ia2,nph
          ! energies in cm-1
          elevel(ii) = evalref(ia1+1)+evalref(ia2+1) &
                  + (etrns+fdlt)*nph
        enddo

        nvlblk = mxlam+nconst ! Total number of Hamiltonian blocks

        write(*,107) mxlam, nconst
 107    format(/" Potential terms will be constructed from ",&
     &       i2," R-dependent term(s) and ",i2," asymptotic term(s).")

        ! Calculate ac electric field if MW field is on
        efldac = 0d0
        if (abs(omega)>zero) then
          ! First calculate angular part
          write(6,48) itrns(1), itrns(2)
  48      format(/" MW transition from states a_i =",i4," to a_f =",i4)
          write(6,49) xi*180d0/pi
  49      format(/" Ellipticity xi in degrees =",f6.1)
          n   = iqn(itrns(1)+1,1)
          mn  = iqn(itrns(1)+1,2)
          np  = iqn(itrns(2)+1,1)
          mnp = iqn(itrns(2)+1,2)
          call vmfdr(n,mn,0,np,mnp,-1,npol,vcoup)
          if (abs(vcoup)<zero) stop "Check ITRNS array and/or NPOL"
          efldac = abs(fomg/dipole/fstark/vcoup) ! in kV/cm
          write(6,50) efldac
  50      format(/" AC electric field (kV/cm) for Rabi coupling is:",&
                  es10.2)
        endif
        
        ! define vconst array in cm-1
        vconst(1)   = 1d0
        vconst(2:3) = -efldac*dipole*fstark/2d0  ! (kV/cm)*debye*(1/kV/debye)

        write(6,'(/" VCONST (cm-1) for MW coupling is ",es16.8)') &
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
        write(10,('("# Adiabats:",/"# R (bohr)  Energies (MHz)",/ &
     &              "#-------------------------")'))
        
        if (ldeng) then
        ! Writes Deng et al effective potential
          open (20,file="Deng_et_al_veff.dat",status="unknown")
          write(20,('("# V_effective:",/"# R (bohr)  Energy (MHz)",/ &
     &                "#-------------------------")'))
        endif

        r = rmin
        istep = 0
        do while (r<=rmax)
          ! Get Hamiltonian at a given value of r-grid
          call hammat(r,nvlblk,p,npair,elevel,wmat) 
          call dsyev('N','U',npair,wmat,npair,eval,work(:3*npair-1),&
                  3*npair-1,info)
          write(10,'(1x,f10.2,50es14.4)') r, (eval-eref)/MHz_to_invcm
          if (ldeng) call effective_pot(r,theta,phi,xi,omega,delta,&
                  dipole,veff)
          if (ldeng) write(20,'(1x,f10.2,es14.4)') r,veff/MHz_to_invcm
          r = r+dr
          istep = istep+1
        enddo

        write(6,"(' No. of R steps:',i10)") istep

        call system_clock(itime_end)
        write(6,'(/" Time (s) in R loop:",es12.2,/)') &
                dble(itime_end-itime_start)/count_rate

        write(6,*)"Adiabats written to file 'calculated_adiabats.dat'"
        if (ldeng) write(6,*)"Deng et al. effective potential written", &
     &          " in file 'Deng_et_al_veff.dat'"

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
