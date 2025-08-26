      subroutine dual_MW_eff_pot (r,theta,omg_sg,del_sg,omg_pi,del_pi,dipole,veff)
        ! This is effective potential for dual MW (RHC+Linear z) proposed by Deng 
        ! et al in arXiv:2501.0521 (2025)
        use constants, only: Debye_in_au, hartree_to_invcm, &
                MHz_to_invcm, pi
        implicit none
        integer, parameter  :: nst = 4
        real*8, intent(in)  :: r              ! in bohr
        real*8, intent(in)  :: theta          ! in radians
        real*8, intent(in)  :: omg_sg, del_sg ! in MHz
        real*8, intent(in)  :: omg_pi, del_pi ! in MHz
        real*8, intent(in)  :: dipole         ! in debye
        real*8, intent(out) :: veff           ! in cm-1

        ! dsyev params
        real*8, dimension(nst,nst) :: hmat, evec
        real*8, dimension(nst)     :: eval
        real*8, dimension(3*nst-1) :: work
        integer                    :: info

        ! potential params
        real*8 :: epls, emis, emi1, ezro, omga
        real*8 :: alph, beta, gamm
        real*8 :: sinalph, cosalph, cos2alph, cos3alph
        real*8 :: sinalpsq, cosalpsq, sinalph4, cosalph4
        real*8 :: sinbeta, cosbeta
        real*8 :: sinbetsq, cosbetsq, sin2beta, cos2beta
        real*8 :: singamm, cosgamm, tangamm
        real*8 :: singamsq, cosgamsq, sin2gamm, cos2gamm, cosgamm4
        real*8 :: term1, term2, term3, term4, term5, term6, term7
        real*8 :: term8, term9, term10, term11, term12, term13
        real*8 :: pisq
        real*8 :: w0, w1, w2
        real*8 :: eta
        real*8 :: C3, C6, C3ang, C6ang
        real*8 :: sintheta, costheta, sinthesq, costhesq, sinthet4
        real*8 :: r3, rinv3, rinv6

        ! Eq A1
        hmat = 0d0
        hmat(1,2) = omg_sg/2d0
        hmat(2,1) = hmat(1,2)
        hmat(1,3) = omg_pi/2d0
        hmat(3,1) = hmat(1,3)
        hmat(2,2) = del_sg
        hmat(3,3) = del_pi
        hmat(4,4) = del_sg

        evec = hmat
        call dsyev('V','U',nst,evec,nst,eval,work,3*nst-1,info)

        ! Eq A2 but columns reordered in ascending order of energy
        ! Columns 1,2,3,4 (in paper) -> 4,3,2,1 (dsyev output)
        alph = acos(evec(1,4))
        if (abs(alph)<1d-16) stop 'dual_MW_eff_pot: alpha is zero'
        sinalph = sin(alph)
        beta = acos(evec(3,4)/sinalph)
        gamm = acos(-evec(1,2)/sinalph)

        epls = eval(4)
        emi1 = eval(3)
        emis = eval(2)
        ezro = eval(1)
        omga = del_pi-del_sg

        cosalph  = cos(alph)
        sinalpsq = sinalph**2
        cosalpsq = 1d0-sinalpsq
        cos2alph = 2d0*cosalpsq-1d0
        cos3alph = 4d0*cosalph**3-3d0*cosalph
        sinalph4 = sinalpsq**2
        cosalph4 = cosalpsq**2

        sinbeta  = sin(beta)
        cosbeta  = cos(beta)
        sinbetsq = sinbeta**2
        cosbetsq = 1d0-sinbetsq
        sin2beta = 2d0*sinbeta*cosbeta
        cos2beta = 2d0*cosbetsq-1d0

        singamm  = sin(gamm)
        cosgamm  = cos(gamm)
        singamsq = singamm**2
        cosgamsq = 1d0-singamsq
        sin2gamm = 2d0*singamm*cosgamm
        cos2gamm = 2d0*cosgamm**2-1d0
        cosgamm4 = cosgamsq**2
        tangamm  = singamm/cosgamm

        pisq = pi**2

        ! Calculation of w0 (Eq C1)
        term1 = 3d0*cosalpsq*sin2beta*cosgamm &
              + 0.5d0*(cosalph+cos3alph)*(3d0*cos2beta-1d0)*singamm
        term1 = sinalpsq*term1**2/(ezro-epls)

        term2 = cos2alph*(1d0-3d0*cos2beta)*cosgamm &
              + 3d0*cosalph*sin2beta*singamm
        term2 = cosalpsq*sinalpsq*term2**2/(emis-epls)

        term3 = 3d0*sin2beta*cosgamm &
              + cosalph*(3d0*cos2beta-1d0)*singamm
        term3 = sinalph4*cosalpsq*singamsq*term3**2/(ezro-epls)

        term4 = 3d0*sin2beta*cos2gamm &
              + cosalph*(3d0*cos2beta-1d0)*sin2gamm
        term4 = sinalph4*cosalpsq*term4**2/(ezro+emis-2d0*epls)

        term5 = 1d0-3d0*cos2beta+3d0*sin2beta*tangamm/cosalph
        term5 = sinalph4*cosalph4*cosgamm4*term5**2/(emis-epls)

        w0 = term1+term2+term3+term4+term5
        w0 = -w0/288d0/pisq  ! MHz^-1

        ! Calculation of w1 (Eq C2)
        term1  = sinalpsq*cosalph4*sinbetsq/(epls-emi1-omga)
        term2  = sinalph4*cosalpsq*sinbetsq*cosgamsq & 
               / (2d0*epls-emis-emi1-omga)
        term3  = sinalph4*cosalpsq*sinbetsq*singamsq & 
               / (2d0*epls-ezro-emi1-omga)
               
        term4  = cosalph*sinbeta*cosgamm+cosbeta*singamm
        term4  = 2d0*sinalph4*cosalpsq*cosbetsq*cosgamsq*term4**2 & 
               / (2d0*epls-2d0*emis+omga)
               
        term5  = sinbeta*cosgamm+cosalph*cosbeta*singamm
        term5  = 2d0*sinalph4*cosalpsq*sinbetsq*singamsq*term5**2 & 
               / (2d0*epls-2d0*ezro-omga)
               
        term6  = cosalph*cosbeta*cosgamm-sinbeta*singamm
        term6  = 2d0*sinalph4*cosalpsq*sinbetsq*cosgamsq*term6**2 & 
               / (2d0*epls-2d0*emis-omga)
               
        term7  = cosbeta*cosgamm-cosalph*sinbeta*singamm
        term7  = 2d0*sinalph4*cosalpsq*cosbetsq*singamsq*term7**2 & 
               / (2d0*epls-2d0*ezro+omga)
               
        term8  = cos2alph*sinbeta*cosgamm+cosalph*cosbeta*singamm
        term8  = sinalpsq*cosalpsq*cosbetsq*term8**2 & 
               / (epls-emis+omga)
               
        term9  = cos2alph*cosbeta*singamm+cosalph*sinbeta*cosgamm
        term9  = sinalpsq*cosalpsq*sinbetsq*term9**2 & 
               / (epls-ezro-omga)

        term10 = cos2alph*cosbeta*cosgamm-cosalph*sinbeta*singamm
        term10 = sinalpsq*cosalpsq*sinbetsq*term10**2 & 
               / (epls-emis-omga)

        term11 = sinbeta*cos2gamm+cosalph*cosbeta*sin2gamm
        term11 = sinalph4*cosalpsq*sinbetsq*term11**2 & 
               / (2d0*epls-ezro-emis-omga)
               
        term12 = cos2alph*sinbeta*singamm-cosalph*cosbeta*singamm
        term12 = sinalpsq*cosalpsq*cosbetsq*term12**2 & 
               / (epls-ezro+omga)

        term13 = cosalph*sinbeta*sin2gamm-cos2gamm*cosbeta
        term13 = sinalph4*cosalpsq*cosbetsq*term13**2 & 
               / (2d0*epls-ezro-emis+omga)

        w1 = term1+term2+term3+term4+term5+term6+term7+term8+term9 &
           + term10+term11+term12+term13
        w1 = w1/4d0/pisq  ! MHz-1
               
        ! Calculation of w2 (Eq C3)
        term1 = cosalph4*cosbetsq*sinalpsq/(epls-emi1)
        term2 = sinalph4*cosalpsq*cosbetsq*cosgamsq/(2d0*epls-emis-emi1)
        term3 = sinalph4*cosalpsq*cosbetsq*singamsq/(2d0*epls-ezro-emi1)
        
        w2 = (term1+term2+term3)/8d0/pisq  ! MHz^-1

        ! Below Eq (2)
        eta = sqrt(8d0*pi/15d0)*(dipole*Debye_in_au)**2*(4d0*pi)  ! hartree bohr^3
        eta = eta*hartree_to_invcm/MHz_to_invcm                   ! MHz bohr^3

        ! Eq (21)
        C3 = (3d0*cos2beta-1d0)*cosalpsq*sinalpsq
        C3 = sqrt(15d0/2d0/pi)*eta/48d0/pi*C3                     ! MHz bohr^3

        ! Below Eq (23)
        C6 = 15d0*eta**2*w2/32d0/pi                               ! MHz bohr^6

        sintheta = sin(theta)
        costheta = cos(theta)
        sinthesq = sintheta**2
        costhesq = 1d0-sinthesq
        sinthet4 = sinthesq**2

        ! Eq (23)
        C3ang = 3d0*costhesq-1d0
        C6ang = sinthet4+sinthesq*costhesq*(w1/w2) & 
              + (3d0*costhesq-1d0)**2*(w0/w2)

        r3    = r**3
        rinv3 = 1d0/r3
        rinv6 = rinv3**2

        veff = C3*C3ang*rinv3+C6*C6ang*rinv6

        veff = veff*MHz_to_invcm

        return
      end subroutine dual_MW_eff_pot

