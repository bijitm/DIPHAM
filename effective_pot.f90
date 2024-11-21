      subroutine effective_pot (r,theta,phi,xi,omega,delta,dipole,veff)
        ! This is effective potential proposed by Deng et al in
        ! PRL 130, 183001 (2023)
        use constants, only: Debye_in_au, hartree_to_invcm, &
                MHz_to_invcm
        implicit none
        real*8, intent(in)  :: r              ! in bohr
        real*8, intent(in)  :: theta, phi, xi ! in radians
        real*8, intent(in)  :: omega, delta   ! in MHz
        real*8, intent(in)  :: dipole         ! in debye
        real*8, intent(out) :: veff           ! in cm-1
        real*8 :: del_r, c3, c6, fcurly, c3ang, c6ang

        del_r = abs(delta/omega)
        c3 = (dipole*Debye_in_au)**2*hartree_to_invcm/ &
                (12d0*(1d0+del_r**2))
        c6 = ((dipole*Debye_in_au)**2*hartree_to_invcm)**2/ &
                (8d0*omega*MHz_to_invcm*sqrt((1d0+del_r**2)**3))
        fcurly = sin(2d0*xi)*cos(2d0*phi)
        c3ang = 3d0*cos(theta)**2-1d0+3d0*fcurly*sin(theta)**2
        c6ang = sin(theta)**2*(1d0-fcurly**2+(1d0-fcurly)**2* &
                cos(theta)**2)
        veff = c3*c3ang/r**3 + c6*c6ang/r**6
        return
      end subroutine effective_pot

