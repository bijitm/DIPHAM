      subroutine thresholds(nvlblk,p,npair,elevel,thrshvals)
        use dipole_module, only: mxlam, vl
        implicit none
        integer, intent(in) :: nvlblk, npair
        real*8, intent(in)  :: p(nvlblk), elevel(npair)
        real*8, intent(out) :: thrshvals(npair)
        integer :: icol, irow, iblk
        real*8  :: wmat(npair,npair)
        real*8  :: work(3*npair-1)
        integer :: info

        wmat = 0d0
        do icol = 1, npair
          do irow = 1, icol

            if (icol==irow) wmat(irow,icol) = wmat(irow,icol) + &
                      p(mxlam+1)*elevel(icol)*vl(mxlam+1,irow,icol)

            do iblk = mxlam+2,nvlblk
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

        call dsyev('N','U',npair,wmat,npair,thrshvals,work(:3*npair-1),&
                3*npair-1,info)

        return
      end subroutine thresholds



