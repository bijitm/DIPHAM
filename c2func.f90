      real*8 function c2func(iq,theta,phi)
        ! This is reduced spherical harmonic C_{2q} function
        ! which depends both on theta and phi. For non-zero phi
        ! the functions are complex. Here we only consider the 
        ! real part of the functions
        use dipole_module, only: zero
        implicit none
        integer, intent(in) :: iq
        real*8, intent(in) :: theta, phi ! in radians

        c2func = 0d0
        if (iq==0) then
           c2func = 0.5d0*(3d0*cos(theta)**2-1d0)
        elseif (abs(iq)==1) then
           c2func = -(-1d0)**iq*sqrt(3d0/2d0)*cos(theta)*sin(theta)&
                   *cos(phi)    ! i*sin(phi) term ignored
        elseif (abs(iq)==2) then
           c2func = sqrt(3d0/8d0)*sin(theta)**2*cos(2d0*phi)
           ! i*sin(2*phi) term ignored
        else
           stop "In c2func: invalid argument"
        endif

        if (abs(c2func)<zero) c2func = 0d0

        return
      end function c2func
