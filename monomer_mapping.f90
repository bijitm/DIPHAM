      subroutine monomer_mapping
         ! Maps monomer rotor quantum numbers with indices later used
         ! as pointers to arrays
         use dipole_module, only: nmax, nmonbasis, nblmn, iqn, jqn, &
                                  ndim, nnidx, mnidx, nqn
         implicit none
         integer :: ifunc, kk, nn, ii, ntot
         integer, allocatable :: mnn(:)

         ! Total number of monomer blocks of given mn.
         nblmn = 2*nmax+1
         if(.not. allocated(mnn)) allocate(mnn(nblmn))
         mnn = -999

         ! Total number of monomer n quantum nos.
         ntot = nmax+1
         mnn(1) = 0
         kk = 1
         do ii = 1, nmax
           mnn(kk+1) = -ii
           mnn(kk+2) =  ii
           kk        =  kk+2
         enddo

         if(.not. allocated(iqn)) allocate(iqn(nmonbasis,nqn-1))
         iqn = -999
         if(.not. allocated(jqn)) allocate(jqn(2*ntot,nblmn))
         jqn = -999

         if(.not. allocated(ndim)) allocate(ndim(nblmn))
         ndim = -999
         if(.not. allocated(nnidx)) allocate(nnidx(0:nmax,nblmn))
         nnidx = -999
         if(.not. allocated(mnidx)) allocate(mnidx(-nmax:nmax))
         mnidx = -999

        ! Store the quantum no. n,mn in arrays IQN and JQN to be used
        ! later.
        ! Store dimension of each block in array NDIM.
        ifunc = 0
        do ii = 1, nblmn
           kk = 0
           do nn = 0, nmax
              if (nn<abs(mnn(ii))) cycle
              ifunc = ifunc+1
              iqn(ifunc,1)    = nn
              iqn(ifunc,2)    = mnn(ii)
              kk              = kk+1
              jqn(kk,ii)      = nn
              jqn(kk+ntot,ii) = mnn(ii)
              nnidx(nn,ii)    = kk
           enddo
           ndim(ii) = kk
           mnidx(mnn(ii)) = ii
        enddo

        if(allocated(mnn)) deallocate(mnn)

        return
      end subroutine monomer_mapping

