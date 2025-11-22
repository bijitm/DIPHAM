      subroutine coupling
        ! Calculates the angular coupling matrix elements of the 
        ! Hamiltonian for different symmetry blocks 
        use dipole_module, only: npair, mxlam, nvlblk, lambda, vl, &
                jlevel, nqn, iqn, lrstrct, etrns, fdlt, mpair, &
                evalref, imlst, zero, xi, npol
        implicit none
        integer :: iham, ii, lm1, lm2, lm, icol, irow, &
                ia1c, ia2c, nphc, ia1r, ia2r, nphr, &
                n1c, mn1c, n2c, mn2c, n1r, mn1r, n2r, mn2r, &
                iint, ia1k, ia2k, nphk, n1k, mn1k, n2k, mn2k
        integer :: nblck, nremn, iblck, istrt
        integer :: itime_start, itime_end
        real*8  :: p1, anormc, anormr, enrc, enrr, enrk, anormk, &
                edifc, edifr, hvvck, hvvkr, term
        real*8  :: count_rate
        logical :: lcoldif, lrowdif, lintdif

        call system_clock(itime_start,count_rate)

        if (.not.allocated(vl)) allocate(vl(nvlblk,npair,npair))
        vl = 0.d0

        do iham = 1, mxlam
          
          ii = (iham-1)*3
          lm1 = lambda(ii+1)
          lm2 = lambda(ii+2)
          lm  = lambda(ii+3)

          ! Electronic C_6 isotropic term
          if (lm1==0.and.lm2==0.and.lm==0) then
            do icol = 1, npair
              vl(iham,icol,icol) = 1.d0
            enddo
    
          ! Dipole-dipole term
          elseif (lm1==1.and.lm2==1.and.lm==2) then
            
            p1 = -sqrt(30d0)

            do icol = 1, npair
              
              ia1c = jlevel(nqn*(icol-1)+1)+1
              ia2c = jlevel(nqn*(icol-1)+2)+1
              nphc = jlevel(nqn*(icol-1)+3)
              lcoldif = .true.
              anormc = sqrt(0.5d0)
              if (ia1c==ia2c) lcoldif = .false.
              if (ia1c==ia2c) anormc  = 1d0

              ! Get actual quantum numbers from IQN array
              n1c  = iqn(ia1c,1)
              mn1c = iqn(ia1c,2)
              n2c  = iqn(ia2c,1)
              mn2c = iqn(ia2c,2)

              do irow = 1, icol

                ia1r = jlevel(nqn*(irow-1)+1)+1
                ia2r = jlevel(nqn*(irow-1)+2)+1
                nphr = jlevel(nqn*(irow-1)+3)
                lrowdif = .true.
                anormr = sqrt(0.5d0)
                if (ia1r==ia2r) lrowdif = .false.
                if (ia1r==ia2r) anormr  = 1d0

                n1r  = iqn(ia1r,1)
                mn1r = iqn(ia1r,2)
                n2r  = iqn(ia2r,1)
                mn2r = iqn(ia2r,2)

                ! Dipole-dipole matrix elements are diagonal in photon states
                if (nphc/=nphr) cycle

                call vdd_dressed(n1c,mn1c,n2c,mn2c,n1r,mn1r,n2r,mn2r,&
                        lcoldif,anormc,lrowdif,anormr,vl(iham,irow,icol))

                vl(iham,irow,icol) = p1*vl(iham,irow,icol)

              enddo ! irow loop
            enddo ! icol loop

          ! -----------------------------------------------
          ! Potential 1st order term evaluation ends here
          ! -----------------------------------------------

          elseif (lm1==1.and.lm2==1.and.lm<0.and.lrstrct) then
          ! LM = -1: This is Van-Vleck in dipole-dipole
          ! LM = -2: 2nd order Van-Vleck in mixed operators:
          !          1st combination: dipole-dipole then MW
          !          2nd combination: MW then dipole-dipole          

            p1 = -sqrt(30d0)

            ! Van Vleck: loop over Class-I functions: explicit basis set
            do icol = 1, npair
              
              ia1c = jlevel(nqn*(icol-1)+1)+1
              ia2c = jlevel(nqn*(icol-1)+2)+1
              nphc = jlevel(nqn*(icol-1)+3)
              lcoldif = .true.
              anormc = sqrt(0.5d0)
              if (ia1c==ia2c) lcoldif = .false.
              if (ia1c==ia2c) anormc  = 1d0

              ! Get actual quantum numbers from IQN array
              n1c  = iqn(ia1c,1)
              mn1c = iqn(ia1c,2)
              n2c  = iqn(ia2c,1)
              mn2c = iqn(ia2c,2)

              ! Get the diagonal energy
              enrc = evalref(ia1c)+evalref(ia2c)+(etrns+fdlt)*nphc

              do irow = 1, icol

                ia1r = jlevel(nqn*(irow-1)+1)+1
                ia2r = jlevel(nqn*(irow-1)+2)+1
                nphr = jlevel(nqn*(irow-1)+3)
                lrowdif = .true.
                anormr = sqrt(0.5d0)
                if (ia1r==ia2r) lrowdif = .false.
                if (ia1r==ia2r) anormr  = 1d0

                n1r  = iqn(ia1r,1)
                mn1r = iqn(ia1r,2)
                n2r  = iqn(ia2r,1)
                mn2r = iqn(ia2r,2)

                enrr = evalref(ia1r)+evalref(ia2r)+(etrns+fdlt)*nphr

                ! Van-Vleck dipole-dipole matrix elements are diagonal in photon states
                if (lm==-1.and.nphc/=nphr) cycle

                ! Van Vleck: loop over Class-II functions
                ! This loop will calculate terms for lm=-1 and first
                ! combination of lm=-2
                do iint = 1, mpair

                  ia1k = imlst(iint,1)
                  ia2k = imlst(iint,2)
                  nphk = imlst(iint,3)

                  n1k  = iqn(ia1k,1)
                  mn1k = iqn(ia1k,2)
                  n2k  = iqn(ia2k,1)
                  mn2k = iqn(ia2k,2)

                  ! Exploit some of 3j-symbol restrictions to skip
                  ! intermediate pair functions
                  if(mn1k==mn2k.and.(abs(mn1k-mn1c)>1.or. &
                          abs(mn1k-mn2c)>1)) cycle
                  if(mn1c==mn2c.and.(abs(mn1c-mn1k)>1.or. &
                          abs(mn1c-mn2k)>1)) cycle
                  if(abs(mn1k-mn1c)>1.and.abs(mn1k-mn2c)>1) cycle
                  if(abs(mn2k-mn1c)>1.and.abs(mn2k-mn2c)>1) cycle

                  ! The first operator of all Van-Vleck term in this
                  ! loop is dipole-dipole operator which cannot change 
                  ! the photon state between column pair function and
                  ! the intermediate pair function

                  ! Van vleck dipole-dipole diagonal in photon states
                  if (nphc/=nphk) cycle

                  ! Photon states must differ by 1 quantum for VV-MW
                  if (lm==-2.and.abs(nphk-nphr)/=1) cycle

                  ! Atleast one monomer function should match for 
                  ! Van-Vleck MW operator
                  if (lm==-2.and.(.not.(ia1k==ia1r.or.ia1k==ia2r&
                          .or.ia2k==ia1r.or.ia2k==ia2r))) cycle

                  ! Van-Vleck MW operator will conserve total projection
                  ! for the intermediate state with row basis function                  
                  if (lm==-2.and.abs(xi)<zero.and.(mn1k+mn2k+npol*nphk &
                      /=mn1r+mn2r+npol*nphr)) cycle

                  if (lm==-2.and.abs(xi)>zero.and.(.not.(&
                       (mn1k+mn2k+nphk==mn1r+mn2r+nphr).or.&
                       (mn1k+mn2k-nphk==mn1r+mn2r-nphr)))) cycle

                  enrk = evalref(ia1k)+evalref(ia2k)+(etrns+fdlt)*nphk
                  edifc = enrc-enrk
                  edifr = enrr-enrk

                  if (abs(edifc)<zero.or.abs(edifr)<zero) then
                    write(*,100)ia1k-1,ia2k-1,nphk
  100               format(/" Near-degenerate state (a1, a2, nph):",&
                           3i7,"  encountered. Not included in", &
                           " the Van-Vleck treatment")
                    cycle
                  endif

                  lintdif = .true.
                  anormk  = sqrt(0.5d0)
                  if (ia1k==ia2k) lintdif = .false.
                  if (ia1k==ia2k) anormk  = 1d0

                  ! Calculate dipole-dipole matrix element between col
                  ! and intermediate basis functions                  
                  call vdd_dressed(n1c,mn1c,n2c,mn2c,n1k,mn1k,n2k,mn2k,&
                       lcoldif,anormc,lintdif,anormk,hvvck)

                  if (abs(hvvck)<zero) cycle
                  hvvck = p1*hvvck
                  hvvkr = 0d0

                  if (lm==-1) then
                    ! Van-Vleck term for 2nd order dipole-dipole

                    call vdd_dressed(n1k,mn1k,n2k,mn2k,n1r,mn1r,n2r,&
                         mn2r,lintdif,anormk,lrowdif,anormr,hvvkr)

                    hvvkr = p1*hvvkr

                  elseif (lm==-2) then
                    ! Van-Vleck term for 2nd order dipole-dipole + MW

                    ! Hm: Monomer Hamiltonian
                    ! <1K|Hm|1R>*krondel(2K,2R)
                    if (ia2k==ia2r) then
                      term = 0d0
                      if (mn1k-mn1r+npol*(nphk-nphr)==0) then
                        call vmfdr(n1k,mn1k,nphk,n1r,mn1r,nphr,npol,term)
                        term = term*cos(xi)
                      endif
                      if (abs(xi)>zero.and.abs(mn1k-mn1r+npol*(nphk-nphr))==2) then
                        call vmfdr(n1k,mn1k,nphk,n1r,mn1r,nphr,-npol,term)
                        term = -term*sin(xi)
                      endif
                      hvvkr = hvvkr+term*anormk*anormr*2d0
                    endif
                    
                    if (lintdif) then
                      ! <2K|Hm|1R>*krondel(1K,2R)
                      if (ia1k==ia2r) then
                        term = 0d0
                        if (mn2k-mn1r+npol*(nphk-nphr)==0) then
                          call vmfdr(n2k,mn2k,nphk,n1r,mn1r,nphr,npol,term)
                          term = term*cos(xi)
                        endif
                        if (abs(xi)>zero.and.abs(mn2k-mn1r+npol*(nphk-nphr))==2) then
                          call vmfdr(n2k,mn2k,nphk,n1r,mn1r,nphr,-npol,term)
                          term = -term*sin(xi)
                        endif
                        hvvkr = hvvkr+term*anormk*anormr*2d0
                      endif
                    endif
                    
                    if (lrowdif) then
                      ! <1K|Hm|2R>*krondel(2K,1R)
                      if (ia2k==ia1r) then
                        term  = 0d0
                        if (mn1k-mn2r+npol*(nphk-nphr)==0) then
                          call vmfdr(n1k,mn1k,nphk,n2r,mn2r,nphr,npol,term)
                          term = term*cos(xi)
                        endif
                        if (abs(xi)>zero.and.abs(mn1k-mn2r+npol*(nphk-nphr))==2) then
                          call vmfdr(n1k,mn1k,nphk,n2r,mn2r,nphr,-npol,term)
                          term = -term*sin(xi)
                        endif
                        hvvkr = hvvkr+term*anormk*anormr*2d0
                      endif
                    endif
                    
                    if (lintdif.and.lrowdif) then
                      ! <2K|Hm|2R>*krondel(1K,1R)
                      if (ia1k==ia1r) then
                        term  = 0d0
                        if (mn2k-mn2r+npol*(nphk-nphr)==0) then
                          call vmfdr(n2k,mn2k,nphk,n2r,mn2r,nphr,npol,term)
                          term = term*cos(xi)
                        endif
                        if (abs(xi)>zero.and.abs(mn2k-mn2r+npol*(nphk-nphr))==2) then
                          call vmfdr(n2k,mn2k,nphk,n2r,mn2r,nphr,-npol,term)
                          term = -term*sin(xi)
                        endif
                        hvvkr = hvvkr+term*anormk*anormr*2d0
                      endif
                    endif ! lintdif and lrowdif if block

                  endif ! lm conditions if block

                  term = 0.5d0*hvvck*hvvkr*((1d0/edifc)+(1d0/edifr))

                  vl(iham,irow,icol) = vl(iham,irow,icol)+term

                enddo ! iint loop

                if (lm==-1) cycle
                ! Now consider the other combination of the mixed VV-MW
                ! LM = -2: 2nd combination: MW then dipole-dipole

                ! Now start the sum over intermediate left over basis
                ! functions again                
                do iint = 1, mpair

                  ia1k = imlst(iint,1)
                  ia2k = imlst(iint,2)
                  nphk = imlst(iint,3)

                  n1k  = iqn(ia1k,1)
                  mn1k = iqn(ia1k,2)
                  n2k  = iqn(ia2k,1)
                  mn2k = iqn(ia2k,2)

                  ! Exploit some of 3j-symbol restrictions to skip
                  ! intermediate pair functions
                  if(mn1k==mn2k.and.(abs(mn1k-mn1c)>1.or. &
                          abs(mn1k-mn2c)>1)) cycle
                  if(mn1c==mn2c.and.(abs(mn1c-mn1k)>1.or. &
                          abs(mn1c-mn2k)>1)) cycle
                  if(abs(mn1k-mn1c)>1.and.abs(mn1k-mn2c)>1) cycle
                  if(abs(mn2k-mn1c)>1.and.abs(mn2k-mn2c)>1) cycle
                  if(abs(mn1k-mn1r)>1.and.abs(mn1k-mn2r)>1) cycle
                  if(abs(mn2k-mn1r)>1.and.abs(mn2k-mn2r)>1) cycle

                  ! Van vleck dipole-dipole diagonal in photon states
                  if (nphk/=nphr) cycle

                  ! Photon states must differ by 1 quantum for VV-MW
                  if (abs(nphc-nphk)/=1) cycle

                  ! Atleast one monomer function should match for 
                  ! Van-Vleck MW operator
                  if (.not.(ia1c==ia1k.or.ia1c==ia2k &
                        .or.ia2c==ia1k.or.ia2c==ia2k)) cycle

                  ! Van-Vleck MW operator will conserve total projection
                  ! for the intermediate state with col basis function                  
                  if (abs(xi)<zero.and.(mn1c+mn2c+npol*nphc &
                      /=mn1k+mn2k+npol*nphk)) cycle

                  if (abs(xi)>zero.and.(.not.(&
                       (mn1c+mn2c+nphc==mn1k+mn2k+nphk).or.&
                       (mn1c+mn2c-nphc==mn1k+mn2k-nphk)))) cycle

                  enrk = evalref(ia1k)+evalref(ia2k)+(etrns+fdlt)*nphk
                  edifc = enrc-enrk
                  edifr = enrr-enrk

                  if (abs(edifc)<zero.or.abs(edifr)<zero) cycle

                  lintdif = .true.
                  anormk  = sqrt(0.5d0)
                  if (ia1k==ia2k) lintdif = .false.
                  if (ia1k==ia2k) anormk  = 1d0

                  ! Calculate dipole-dipole matrix element between row
                  ! and intermediate basis functions                  
                  call vdd_dressed(n1k,mn1k,n2k,mn2k,n1r,mn1r,n2r,mn2r,&
                       lintdif,anormk,lrowdif,anormr,hvvkr)

                  if (abs(hvvkr)<zero) cycle
                  hvvkr = p1*hvvkr
                  hvvck = 0d0

                  ! Hm: Monomer Hamiltonian
                  ! <1C|Hm|1K>*krondel(2C,2K)
                  if (ia2c==ia2k) then
                    term = 0d0
                    if (mn1c-mn1k+npol*(nphc-nphk)==0) then
                      call vmfdr(n1c,mn1c,nphc,n1k,mn1k,nphk,npol,term)
                      term = term*cos(xi)
                    endif
                    if (abs(xi)>zero.and.abs(mn1c-mn1k+npol*(nphc-nphk))==2) then
                      call vmfdr(n1c,mn1c,nphc,n1k,mn1k,nphk,-npol,term)
                      term = -term*sin(xi)
                    endif
                    hvvck = hvvck+term*anormc*anormk*2d0
                  endif
                  
                  if (lcoldif) then
                    ! <2C|Hm|1K>*krondel(1C,2K)
                    if (ia1c==ia2k) then
                      term = 0d0
                      if (mn2c-mn1k+npol*(nphc-nphk)==0) then
                        call vmfdr(n2c,mn2c,nphc,n1k,mn1k,nphk,npol,term)
                        term = term*cos(xi)
                      endif
                      if (abs(xi)>zero.and.abs(mn2c-mn1k+npol*(nphc-nphk))==2) then
                        call vmfdr(n2c,mn2c,nphc,n1k,mn1k,nphk,-npol,term)
                        term = -term*sin(xi)
                      endif
                      hvvck = hvvck+term*anormc*anormk*2d0
                    endif
                  endif
                  
                  if (lintdif) then
                    ! <1C|Hm|2K>*krondel(2C,1K)
                    if (ia2c==ia1k) then
                      term = 0d0
                      if (mn1c-mn2k+npol*(nphc-nphk)==0) then
                        call vmfdr(n1c,mn1c,nphc,n2k,mn2k,nphk,npol,term)
                        term = term*cos(xi)
                      endif
                      if (abs(xi)>zero.and.abs(mn1c-mn2k+npol*(nphc-nphk))==2) then
                        call vmfdr(n1c,mn1c,nphc,n2k,mn2k,nphk,-npol,term)
                        term = -term*sin(xi)
                      endif
                      hvvck = hvvck+term*anormc*anormk*2d0
                    endif
                  endif
                  
                  if (lcoldif.and.lintdif) then
                    ! <2C|Hm|2K>*krondel(1C,1K)
                    if (ia1c==ia1k) then
                      term = 0d0
                      if (mn2c-mn2k+npol*(nphc-nphk)==0) then
                        call vmfdr(n2c,mn2c,nphc,n2k,mn2k,nphk,npol,term)
                        term = term*cos(xi)
                      endif
                      if (abs(xi)>zero.and.abs(mn2c-mn2k+npol*(nphc-nphk))==2) then
                        call vmfdr(n2c,mn2c,nphc,n2k,mn2k,nphk,-npol,term)
                        term = -term*sin(xi)
                      endif
                      hvvck = hvvck+term*anormc*anormk*2d0
                    endif
                  endif ! lintdif and lrowdif if block

                  term = 0.5d0*hvvck*hvvkr*((1d0/edifc)+(1d0/edifr))

                  vl(iham,irow,icol) = vl(iham,irow,icol)+term

                enddo ! iint loop

              enddo ! irow loop
            enddo ! icol loop

          else
            write(6,'(/" The set of lambda terms:",3i4,&
     &              " will be ignored")')lm1,lm2,lm
          endif

        enddo ! iham loop

        !--------------------------------------------------------------
        ! Now build VL matrix elements for asymptotic Hamiltonian terms
        !--------------------------------------------------------------        

        do icol = 1, npair

          ia1c = jlevel(nqn*(icol-1)+1)+1
          ia2c = jlevel(nqn*(icol-1)+2)+1
          nphc = jlevel(nqn*(icol-1)+3)
          lcoldif = .true.
          anormc = sqrt(0.5d0)
          if (ia1c==ia2c) lcoldif = .false.
          if (ia1c==ia2c) anormc  = 1d0

          ! Get actual quantum numbers from IQN array
          n1c  = iqn(ia1c,1)
          mn1c = iqn(ia1c,2)
          n2c  = iqn(ia2c,1)
          mn2c = iqn(ia2c,2)

          do irow = 1, icol

            ia1r = jlevel(nqn*(irow-1)+1)+1
            ia2r = jlevel(nqn*(irow-1)+2)+1
            nphr = jlevel(nqn*(irow-1)+3)
            lrowdif = .true.
            anormr = sqrt(0.5d0)
            if (ia1r==ia2r) lrowdif = .false.
            if (ia1r==ia2r) anormr  = 1d0

            n1r  = iqn(ia1r,1)
            mn1r = iqn(ia1r,2)
            n2r  = iqn(ia2r,1)
            mn2r = iqn(ia2r,2)

            ! If no ellipticity then total projection of pair state 
            ! should be conserved
            if (abs(xi)<zero .and. (mn1c+mn2c+npol*nphc /= &
                    mn1r+mn2r+npol*nphr)) cycle

            ! If ellipticity, then check total projection with both
            ! + and - circular polarization
            if (abs(xi)>zero .and. (.not.((mn1c+mn2c+nphc == &
                    mn1r+mn2r+nphr).or.(mn1c+mn2c-nphc == &
                    mn1r+mn2r-nphr)))) cycle

            do iham = mxlam+1, nvlblk
              if (iham==mxlam+1.and.icol==irow) then
                ! dc Stark+rotation+MW detuning
                vl(iham,irow,icol) = 1d0
              endif

              if (iham==mxlam+1) cycle

              ! Photon states must differ by 1 quantum
              if (abs(nphc-nphr)/=1) cycle

              ! Hm: Monomer Hamiltonian
              ! <1C|Hm|1R>*krondel(2C,2R)
              if (ia2c==ia2r) then
                term = 0d0
                if (mn1c-mn1r+npol*(nphc-nphr)==0) then
                  call vmfdr(n1c,mn1c,nphc,n1r,mn1r,nphr,npol,term)
                  term = term*cos(xi)
                endif
                if (abs(xi)>zero.and.abs(mn1c-mn1r+npol*(nphc-nphr))==2) then
                  call vmfdr(n1c,mn1c,nphc,n1r,mn1r,nphr,-npol,term)
                  term = -term*sin(xi)
                endif
                vl(iham,irow,icol) = vl(iham,irow,icol) &
                        +term*anormc*anormr*2d0
              endif

              if (lcoldif) then
                ! <2C|Hm|1R>*krondel(1C,2R)
                if (ia1c==ia2r) then
                  term = 0d0
                  if (mn2c-mn1r+npol*(nphc-nphr)==0) then
                    call vmfdr(n2c,mn2c,nphc,n1r,mn1r,nphr,npol,term)
                    term = term*cos(xi)
                  endif
                  if (abs(xi)>zero.and.abs(mn2c-mn1r+npol*(nphc-nphr))==2) then
                    call vmfdr(n2c,mn2c,nphc,n1r,mn1r,nphr,-npol,term)
                    term = -term*sin(xi)
                  endif
                  vl(iham,irow,icol) = vl(iham,irow,icol) &
                        +term*anormc*anormr*2d0
                endif
              endif

              if (lrowdif) then
                ! <1C|Hm|2R>*krondel(2C,1R)
                if (ia2c==ia1r) then
                  term  = 0d0
                  if (mn1c-mn2r+npol*(nphc-nphr)==0) then
                    call vmfdr(n1c,mn1c,nphc,n2r,mn2r,nphr,npol,term)
                    term = term*cos(xi)
                  endif
                  if (abs(xi)>zero.and.abs(mn1c-mn2r+npol*(nphc-nphr))==2) then
                    call vmfdr(n1c,mn1c,nphc,n2r,mn2r,nphr,-npol,term)
                    term = -term*sin(xi)
                  endif
                  vl(iham,irow,icol) = vl(iham,irow,icol) &
                        +term*anormc*anormr*2d0
                endif
              endif

              if (lcoldif.and.lrowdif) then
                ! <2C|Hm|2R>*krondel(1C,1R)
                if (ia1c==ia1r) then
                  term  = 0d0
                  if (mn2c-mn2r+npol*(nphc-nphr)==0) then
                    call vmfdr(n2c,mn2c,nphc,n2r,mn2r,nphr,npol,term)
                    term = term*cos(xi)
                  endif
                  if (abs(xi)>zero.and.abs(mn2c-mn2r+npol*(nphc-nphr))==2) then
                    call vmfdr(n2c,mn2c,nphc,n2r,mn2r,nphr,-npol,term)
                    term = -term*sin(xi)
                  endif
                  vl(iham,irow,icol) = vl(iham,irow,icol) &
                        +term*anormc*anormr*2d0
                endif
              endif

            enddo ! iham loop

          enddo ! irow loop
        enddo ! icol loop

        ! Symmetrize VL matrix
        do icol = 1, npair
          do irow = icol+1, npair
            vl(:,irow,icol) = vl(:,icol,irow)
          enddo
        enddo

        ! Print vl matrix elements
        do iham = 1, nvlblk
          write(6,"(/' Potential term',i3/)") iham
          if (.not. any(abs(vl(iham,:,:))>zero)) then
            write(6,*)"All elements for this potential term are zero"
            cycle
          endif

          nblck = npair/10
          nremn = mod(npair,10)

          do iblck = 1, nblck
            write(6,*)
            istrt=(iblck-1)*10
            write(6,"(2x,10i12)") (icol,icol=istrt+1,istrt+10)
            do icol = istrt+1, npair
              write(6,"(i5,10es12.2)") icol, (vl(iham,irow,icol),&
                      irow=istrt+1,min(istrt+10,icol))
            enddo
          enddo

          if (nremn == 0) cycle
          write(6,"(2x,9i12)") (icol,icol=10*nblck+1,npair)
          do icol = 10*nblck+1, npair
            write(6,"(i5,9es12.2)") icol, (vl(iham,irow,icol),&
                    irow=10*nblck+1,icol)
          enddo

        enddo

        call system_clock(itime_end)
        write(*,'(/" Time (s) in coupling:",es12.2,/)') &
                dble(itime_end-itime_start)/count_rate

        return
      end subroutine coupling

      subroutine vdd_dressed(n1c,mn1c,n2c,mn2c,n1r,mn1r,n2r,mn2r,&
                      lcoldif,anormc,lrowdif,anormr,vcoup)
        ! Calculates the dipole-dipole matrix element in the dressed
        ! rotor basis set
        use dipole_module, only: mnidx, ndim, nnidx, jqn, nmax, &
                evecref, zero, theta, phi
        implicit none
        integer, intent(in) :: n1c, mn1c, n2c, mn2c, n1r, mn1r, n2r, mn2r
        real*8, intent(in)  :: anormc, anormr
        logical, intent(in) :: lcoldif, lrowdif
        real*8, intent(out) :: vcoup

        integer :: ib1c, ib2c, ib1r, ib2r, nd1c, nd2c, &
                in1c, in2c, in1r, in2r, i1c, i2c, i1r, i2r, &
                n1, mn1, n2, mn2, n1p, mn1p, n2p, mn2p, &
                mn1dif, mn2dif
        real*8  :: p1, p2, p3, p4, p5, elem, term

        ! Functions
        real*8  :: threej, thrj, c2func

        vcoup = 0d0

        ! Pointer to the block for this MN1C
        ib1c = mnidx(mn1c)
        ! Size of the block for this MN1C
        nd1c = ndim(ib1c)
        ! Pointer to the index for this N1C within block IB1C
        in1c = nnidx(n1c,ib1c)

        ib2c = mnidx(mn2c)
        nd2c = ndim(ib2c)
        in2c = nnidx(n2c,ib2c)

        ib1r = mnidx(mn1r)
        in1r = nnidx(n1r,ib1r)

        ib2r = mnidx(mn2r)
        in2r = nnidx(n2r,ib2r)        

        ! Now loops over primitive basis set
        ! Just loop over the relevant block for MN1C        
        do i1c = 1, nd1c
         ! Get N1 from JQN array
         n1  = jqn(i1c,ib1c)
         mn1 = mn1c
         if (abs(evecref(i1c,in1c,ib1c))<zero) cycle
         p1  = sqrt(dble(2*n1+1))

         do i2c = 1, nd2c
           n2  = jqn(i2c,ib2c)
           mn2 = mn2c
           if (abs(evecref(i2c,in2c,ib2c))<zero) cycle
           p2  = p1*sqrt(dble(2*n2+1))
           if (mod(mn1+mn2,2)/=0) p2 = -p2

           ! Loop over primed quantum numbers according to 3-j symbol restrictions
           do n1p = abs(n1-1), n1+1

             ! Considering <AC|BC|V|AR|BR> term
             if (n1p>nmax) cycle
             mn1p = mn1r
             if (n1p<abs(mn1p)) cycle
             ! Pointer to the index for this N1P within block IB1R
             i1r = nnidx(n1p,ib1r)
             if (abs(evecref(i1r,in1r,ib1r))<zero) cycle
             ! 3j restriction
             if (n1==n1p) cycle
             p3 = p2*threej(n1,1,n1p)*sqrt(dble(2*n1p+1))
             if (abs(p3)<zero) cycle
             mn1dif = mn1-mn1p
             p4 = p3*thrj( dble(n1), 1d0, dble(n1p), &
                     -dble(mn1), dble(mn1dif), dble(mn1p))
             if(abs(p4)<zero) cycle

             do n2p = abs(n2-1), n2+1

               if (n2p>nmax) cycle
               mn2p = mn2r
               if (n2p<abs(mn2p)) cycle
               i2r = nnidx(n2p,ib2r)
               if (abs(evecref(i2r,in2r,ib2r))<zero) cycle
               ! 3j restriction
               if (n2==n2p) cycle
               p5 = p4*threej(n2,1,n2p)*sqrt(dble(2*n2p+1))
               if (abs(p5)<zero) cycle
               mn2dif = mn2-mn2p

               elem = p5*thrj( dble(n2), 1d0, dble(n2p), &
                        -dble(mn2), dble(mn2dif), dble(mn2p)) &
                     *thrj(1d0, 1d0, 2d0, &
                     dble(mn1dif), dble(mn2dif), -dble(mn1dif+mn2dif))

               if(abs(elem)<zero) cycle

               term = evecref(i1c,in1c,ib1c)*evecref(i2c,in2c,ib2c)* &
                      evecref(i1r,in1r,ib1r)*evecref(i2r,in2r,ib2r)* &
                      elem*anormc*anormr* &
                      c2func(-(mn1dif+mn2dif),theta,phi)

               vcoup = vcoup+term

             enddo ! n2p loop
           enddo ! n1p loop

           ! Next to consider exchange symmetry if both column and
           ! row pair states are non-identical
           ! If one of column or row pair state is identical then
           ! skip the rest
           if ((.not.lcoldif).or.(.not.lrowdif)) cycle

           ! Recalculate the primed primitive pair states by
           ! exchanging them against the original explicit row basis
           ! set           

           do n1p = abs(n1-1), n1+1

             ! Considering <AC|BC|V|BR|AR> term
             if (n1p>nmax) cycle
             mn1p = mn2r
             if (n1p<abs(mn1p)) cycle
             i1r = nnidx(n1p,ib2r)
             if (abs(evecref(i1r,in2r,ib2r))<zero) cycle
             ! 3j restriction
             if (n1==n1p) cycle
             p3 = p2*threej(n1,1,n1p)*sqrt(dble(2*n1p+1))
             if (abs(p3)<zero) cycle
             mn1dif = mn1-mn1p
             p4 = p3*thrj( dble(n1), 1d0, dble(n1p), &
                     -dble(mn1), dble(mn1dif), dble(mn1p))
             if(abs(p4)<zero) cycle

             do n2p = abs(n2-1), n2+1

               if (n2p>nmax) cycle
               mn2p = mn1r
               if (n2p<abs(mn2p)) cycle
               i2r = nnidx(n2p,ib1r)
               if (abs(evecref(i2r,in1r,ib1r))<zero) cycle
               ! 3j restriction
               if (n2==n2p) cycle
               p5 = p4*threej(n2,1,n2p)*sqrt(dble(2*n2p+1))
               if (abs(p5)<zero) cycle
               mn2dif = mn2-mn2p

               elem = p5*thrj( dble(n2), 1d0, dble(n2p), &
                        -dble(mn2), dble(mn2dif), dble(mn2p)) &
                     *thrj(1d0, 1d0, 2d0, &
                     dble(mn1dif), dble(mn2dif), -dble(mn1dif+mn2dif))

               if(abs(elem)<zero) cycle

               term = evecref(i1c,in1c,ib1c)*evecref(i2c,in2c,ib2c)* &
                      evecref(i1r,in2r,ib2r)*evecref(i2r,in1r,ib1r)* &
                      elem*anormc*anormr* &
                      c2func(-(mn1dif+mn2dif),theta,phi)

               vcoup = vcoup+term

             enddo ! n2p loop
           enddo ! n1p loop

         enddo ! i2c loop
        enddo ! i1c loop

        ! Due to the symmetry of the dipole-dipole
        ! potential form, exchanging the column or row or
        ! both pair states will have same value
        if (lcoldif.or.lrowdif) vcoup = 2d0*vcoup

        return
      end subroutine vdd_dressed

      subroutine vmfdr(ndc,mnc,nphc,ndr,mnr,nphr,npol,vcoup)
        ! Calculates the angular part of MW coupling using dressed
        ! rotor basis sets
        use dipole_module, only: mnidx, ndim, nnidx, jqn, evecref, &
                zero
        implicit none
        integer, intent(in) :: ndc, mnc, nphc, ndr, mnr, nphr, npol
        real*8, intent(out) :: vcoup

        integer :: ibc, ibr, ndimc, ndimr, inc, inr, ic, ir, nc, nr
        real*8  :: term

        ! Function
        real*8  :: c1func

        vcoup = 0d0
        ! Photon quantum numbers must differ by 1
        if (abs(nphc-nphr)/=1) return

        ! Pointers to the block for MNC and MNR
        ibc = mnidx(mnc)
        ibr = mnidx(mnr)
        ! Sizes of the blocks for MNC and MNR
        ndimc = ndim(ibc)
        ndimr = ndim(ibr)
        ! Pointers to the indices for NDC and NDR
        inc = nnidx(ndc,ibc)
        inr = nnidx(ndr,ibr)

        do ic = 1, ndimc
          if (abs(evecref(ic,inc,ibc))<zero) cycle
          ! now get the primitive quantum number
          nc = jqn(ic,ibc)
          do ir = 1, ndimr
            if (abs(evecref(ir,inr,ibr))<zero) cycle
            nr = jqn(ir,ibr)
            ! n should differ by 1
            if (abs(nc-nr)/=1) cycle

            if (nphc==nphr-1) then
               term = c1func(nc,mnc,nr,mnr,npol)
!              if (npol==1) term = -term
            elseif (nphc==nphr+1) then
               term = c1func(nc,mnc,nr,mnr,-npol)
!              if (-npol==1) term = -term
               if (abs(npol)==1) term = -term
            endif
            vcoup = vcoup+evecref(ic,inc,ibc)*evecref(ir,inr,ibr)*term
          enddo ! ir loop
        enddo ! ic loop

        return
      end subroutine vmfdr

