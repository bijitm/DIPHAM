      subroutine pairbasis_builder
        ! Creates dressed pair basis set, calculates eigenvalues and
        ! eigenvectors of the dressed pair basis set.
        ! Create Class-II function list to be used for Van Vleck
        ! transformations
        use dipole_module, only: nmax, nblmn, ndim, evecref, evalref, &
                nmonbasis, lrstrct, nqn, npair, ipair, nphmn, nphmx, &
                mpair, imlst, itrns, etrns
        implicit none
        integer :: ifunc, ibl, nn, ii, jj, nsym, nph, ia1, ia2, kk, &
                im1, im2, idx, itmp
        integer :: itime_start, itime_end
        real*8, allocatable :: hmon(:,:), evec(:,:), eval(:), work(:)
        real*8  :: count_rate
        integer, allocatable :: ipair_idx(:)
        integer :: info

        ! Function
        integer :: func_cantor

        call system_clock(itime_start,count_rate)

        if (.not. allocated(hmon)) allocate(hmon(nmax+1,nmax+1))
        if (.not. allocated(evec)) allocate(evec(nmax+1,nmax+1))
        if (.not. allocated(eval)) allocate(eval(nmax+1))
        if (.not. allocated(evalref)) allocate(evalref(nmonbasis))
        if (.not. allocated(evecref)) &
                allocate(evecref(nmax+1,nmax+1,nblmn))
        if (.not. allocated(work)) allocate(work(3*(nmax+1)-1))

        hmon    = 0.d0
        evec    = 0.d0
        eval    = 0.d0
        evalref = 0.d0
        evecref = 0.d0
        work    = 0.d0

        ifunc = 1
        do ibl = 1, nblmn
          nn = ndim(ibl)
          ii = ifunc
          jj = ifunc+nn-1
          if (ibl==1 .or. mod(ibl,2)==0) then
            call hamstark(ibl,nn,hmon(:nn,:nn))
            call dsyev('V','U',nn,hmon(:nn,:nn),nn,eval(:nn),&
                    work(:3*nn-1),3*nn-1,info)
            evec(:nn,:nn) = hmon(:nn,:nn)
          endif
          evecref(:nn,:nn,ibl) = evec(:nn,:nn)
          evalref(ii:jj)       = eval(:nn)
          ifunc                = ifunc+nn
        enddo

        etrns = 0d0
        if (.not. any(itrns==-999)) etrns = evalref(itrns(2)+1)-evalref(itrns(1)+1)

        if (allocated(work)) deallocate(work)
        if (allocated(hmon)) deallocate(hmon)
        if (allocated(evec)) deallocate(evec)
        if (allocated(eval)) deallocate(eval)


        ! Calculation of monomer eigenvalues and eigenvectors finishes here
        ! ----------------------------------------------------------------

        ! ----------------------------------------------------------------
        ! Assigns monomer quantum numbers and energies to pair states
        ! ----------------------------------------------------------------
        ! Pair states selection can be done by two mechanisms:
        ! 1. All pair states are included (usually a large basis set will
        !    be created) - LRSTRCT must be false.
        ! 2. Subset of pair states are to be selected using NPAIR
        !    mechanism. NPAIR should be a positive integer specifying
        !    total number of required pair states. IPAIR then should
        !    contain sets of a1, a2 quantum numbers for those pair states.

        ! Total number of symmetrized pairs allowed
        nsym = (nphmx-nphmn+1)*nmonbasis*(nmonbasis+1)/2

        if (.not. lrstrct) then

          write(6,109)
 109      format(/" Pair Hamiltonian basis set will be constructed ",&
     &           "using full direct product of monomer basis ",&
     &           "functions.")

          if (allocated(ipair)) deallocate(ipair)
          allocate(ipair(nsym*nqn+1))
          ipair = -999

          ii = 0
          jj = 0
          do ia1 = 1, nmonbasis
            do ia2 = ia1, nmonbasis
              do nph = nphmn, nphmx
                ii          = ii+1
                ipair(jj+1) = ia1-1
                ipair(jj+2) = ia2-1
                ipair(jj+3) = nph
                jj          = jj+nqn
              enddo
            enddo
          enddo
          npair = ii

        else

          write(6,112)
 112      format(/" Pair Hamiltonian basis set will be constructed ",&
     &           "using user input npair and ipair.")

          if (.not. allocated(ipair_idx)) allocate(ipair_idx(npair))
          ipair_idx = -999

          kk = 0
          do ii = 1, npair
            ia1 = min(ipair(kk+1),ipair(kk+2))
            ia2 = max(ipair(kk+1),ipair(kk+2))
            nph = ipair(kk+3)
            if (nph<nphmn .or. nph>nphmx) then
               stop "In pair_basis_builder: nph in ipair is out of range"
            endif
            ! Assign a tag to the triples using Cantor's pairing function
            itmp = func_cantor(ia1,ia2) 
            ipair_idx(ii) = func_cantor(itmp,nph)
            kk = kk+nqn
          enddo

          ! Number of functions excluded from the actual basis functions
          ! functions
          mpair = nsym-npair
          if (.not. allocated(imlst)) allocate(imlst(mpair,nqn))
          imlst = -999

          ! Now store the information of the functions excluded from
          ! the restricted basis set
          kk = 0
          do im1 = 1, nmonbasis
            do im2 = im1, nmonbasis
              do nph = nphmn, nphmx
                itmp = func_cantor(im1-1,im2-1)
                idx = func_cantor(itmp,nph)
                ! Inclusion of this function if tag is not present in
                ! the ipair_idx array
                if (any(ipair_idx-idx==0)) cycle
                kk = kk+1
                imlst(kk,1) = im1
                imlst(kk,2) = im2
                imlst(kk,3) = nph
              enddo
            enddo
          enddo          

        endif ! IF (.NOT. LRSTRCT) block

        call system_clock(itime_end)
        write(*,'(/" Time (s) in pairbasis_builder:",es12.2,/)') &
                dble(itime_end-itime_start)/count_rate

        ! Now deallocate array
        if (allocated(ipair_idx)) deallocate(ipair_idx)

        return
      end subroutine pairbasis_builder

      subroutine hamstark(ib,nham,ham)
        use dipole_module, only: jqn, nmax, frot, ffld
        implicit none
        integer, intent(in) :: ib, nham
        real*8, intent(out) :: ham(nham,nham)
        real*8, dimension(nham,nham) :: hrot, hstark
        integer :: kc, kr, icol, irow, nc, nr, mnc, mnr

        ! Functions
        real*8 :: c1func

        hrot = 0.d0
        hstark = 0.d0

        kc = 0
        do icol = 1, nham
          kc  = kc+1
          nc  = jqn(kc,ib)
          mnc = jqn(kc+nmax+1,ib)
          kr  = 0
          do irow = 1, icol
            kr  = kr+1
            nr  = jqn(kr,ib)
            mnr = jqn(kr+nmax+1,ib)
            if (irow==icol) hrot(irow,icol) = dble((nc+1)*nc)
            if (mnc==mnr .and. abs(nc-nr)==1) &
              hstark(irow,icol) = c1func(nc,mnc,nr,mnc,0)
          enddo
        enddo

        do icol = 1, nham
          do irow = icol+1, nham
            hstark(irow,icol) = hstark(icol,irow)
          enddo
        enddo

        ! In cm-1 now
        hrot   = frot*hrot
        hstark = ffld*hstark
        ham    = hrot+hstark

        return
      end subroutine hamstark

      integer function func_cantor(k1,k2)
        implicit none
        integer, intent(in) :: k1, k2
        func_cantor = k2+(k1+k2)*(k1+k2+1)/2
        return
      end function func_cantor
