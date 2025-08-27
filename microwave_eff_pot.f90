      subroutine microwave_eff_pot (r,theta,phi,xi,omega,delta,dipole,veff)
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
        real*8 :: del_r, sinsq, cossq
        real*8 :: dipfc, denom
        real*8 :: c3, c6, felli, c3ang, c6ang
        real*8 :: r3inv, r6inv

        if (omega<1d-16) stop "microwave_eff_pot: omega is zero"

        del_r = abs(delta/omega)
        sinsq = sin(theta)**2
        cossq = 1d0-sinsq
        felli = sin(2d0*xi)*cos(2d0*phi)

        c3ang = 3d0*cossq-1d0+3d0*felli*sinsq
        c6ang = 1d0-felli**2+(1d0-felli)**2*cossq
        c6ang = c6ang*sinsq

        dipfc = (dipole*Debye_in_au)**2*hartree_to_invcm ! cm-1 bohr^3
        denom = 1d0+del_r**2

        c3    = dipfc/12d0/denom                   ! cm-1 bohr^3
        c6    = dipfc**2/8d0/omega/sqrt(denom**3) &
              / MHz_to_invcm                       ! cm-1 bohr^6

        r3inv = 1d0/r**3
        r6inv = r3inv**2

        veff = c3*c3ang*r3inv + c6*c6ang*r6inv
        return
      end subroutine microwave_eff_pot

