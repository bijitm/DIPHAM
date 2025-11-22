      module dipole_module
         implicit none
         integer, parameter :: nqn = 3, nconst = 2
         ! NCONST: 1. Diagonal Stark+Rot+MW detuning;
         ! NCONST: 2. MW Rabi coupling
         real*8,  parameter :: zero = 1d-12
         integer :: nmax, npair, nmonbasis, nblmn, nphmn, nphmx, &
                    mpair, itrns(2), mxlam, nvlblk, npol, npower(100)
         real*8  :: frot, ffld, fdlt, etrns, theta, phi, xi
         integer, allocatable, dimension(:) :: ipair, ndim, mnidx, &
                 lambda, jlevel
         integer, allocatable, dimension(:,:) :: iqn, jqn, nnidx, &
                 imlst
         real*8, allocatable :: evalref(:), evecref(:,:,:), vl(:,:,:)
         logical :: lrstrct
      end module dipole_module

