      subroutine hammat(r,nvlblk,p,npair,elevel,wmat)
        ! Creates the total R-dependent Hamiltonian
        use dipole_module, only: mxlam, npower, vl
        implicit none
        integer, intent(in) :: nvlblk, npair
        real*8, intent(in)  :: r, p(nvlblk), elevel(npair)
        real*8, intent(out) :: wmat(npair,npair)
        integer :: iblk, icol, irow

        wmat = 0.d0

        do icol = 1, npair
          do irow = 1, icol

            do iblk = 1, mxlam
              wmat(irow,icol) = wmat(irow,icol) + &
                      p(iblk)*vl(iblk,irow,icol)*r**(npower(iblk))
            enddo ! iblk loop

            if (icol==irow) wmat(irow,icol) = wmat(irow,icol) + &
                      p(mxlam+1)*elevel(icol)*vl(mxlam+1,irow,icol)

            do iblk = mxlam+2, nvlblk
              wmat(irow,icol) = wmat(irow,icol) + &
                      p(iblk)*vl(iblk,irow,icol)
            enddo ! iblk loop

          enddo ! irow loop
        enddo ! icol loop

        ! Symmetrize wmat
        do icol = 1, npair
          do irow = icol+1, npair
            wmat(irow,icol) = wmat(icol,irow)
          enddo
        enddo

        return
      end subroutine hammat

