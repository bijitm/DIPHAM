      real*8 function c1func(n,mn,np,mnp,mq)
        implicit none
        integer, intent(in) :: n, mn, np, mnp, mq
        real*8 :: thrj, threej
        c1func = 0d0
        if (mn-mnp/=mq) return
        c1func = sqrt(dble((2*n+1)*(2*np+1))) &
              * threej(n,1,np) &
              * thrj(dble(n),1.d0,dble(np), &
                    -dble(mn),dble(mq),dble(mnp))
        if (mod(mn,2)/=0) c1func = -c1func
        return
      end function c1func

