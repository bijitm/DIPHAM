      real*8 function c2func(iq,theta)
        ! This is reduced spherical harmonic C_{2q} function
        ! which depends both on theta and phi. For non-zero phi
        ! the functions are complex, so we assume phi = 0 here.
        implicit none
        integer, intent(in) :: iq
        real*8, intent(in) :: theta ! in radians

        c2func = 0d0
        if (iq==0) then
           c2func = 0.5d0*(3d0*cos(theta)**2-1d0)
        elseif (abs(iq)==1) then
           c2func = -(-1d0)**iq*sqrt(3d0/2d0)*cos(theta)*sin(theta)
        elseif (abs(iq)==2) then
           c2func = sqrt(3d0/8d0)*sin(theta)**2
        else
           stop "In c2func: invalid argument"
        endif

        return
      end function c2func
