module constants
    implicit none
    real(kind=8), parameter :: hcross = 0.06350779929588898d0, & ! in eps*tau;  1 eps = 1 amu*ang^2/tau^2; 1 tau = 1.e-14 sec
&                              bohr_to_ang = 0.529177210903d0, &
&                              hartree_to_eV = 27.211386245988034d0, &
&                              hartree_to_inv_cm = 2.1947463136320d5, &
&                              eps_to_eV = 1.0364269652680504d0, &
&                              C6au_to_C6invang = hartree_to_inv_cm*bohr_to_ang**6,& ! = 4819.37950023638
&                              Hz_to_inv_cm = 3.335640951d-11, &
&                              MHz_to_inv_cm = Hz_to_inv_cm * 1.d6, & ! = 3.335640951d-5
&                              K_to_inv_cm = 0.6950348004d0, &
&                              K_to_MHz = K_to_inv_cm/MHz_to_inv_cm, & ! = 20836.619126876765 
&                              eV_to_inv_cm = 8065.543937349212d0, &
&                              eps_to_inv_cm = eps_to_eV * eV_to_inv_cm, &  ! = 8359.347226222966
&                              eps_to_MHz = eps_to_inv_cm / MHz_to_inv_cm, &  ! = 250606925.29622823
&                              eps_to_K = eps_to_inv_cm / K_to_inv_cm, &    ! = 12027.235501606641
&                              eV_to_MHz = eV_to_inv_cm / MHz_to_inv_cm, &  ! = 241798924.27964178
&                              hartree_to_eps = hartree_to_eV / eps_to_eV, & ! = 26.25499640387142
&                              bfct = 16.85762919164018d0, & ! in cm^-1 amu ang^2. It is basically hcross^2*eps_to_inv_cm/2.
&                              atomic_mass_constant = 1.66053906660d-27, & ! in Kg.
&                              electron_mass = 9.1093837015d-31, & ! in Kg.
&                              amu_to_au = atomic_mass_constant/electron_mass, & ! = 1822.8884862173131
&                              speed_of_light_in_SI = 2.99792458d8, &  ! ms^-1
&                              Planck_constant_in_SI = 6.62607015d-34, &  ! Js
&                              Boltzmann_constant_in_SI = 1.380649d-23, & ! JK^-1
&                              Joule_to_inv_cm = 1d-2/Planck_constant_in_SI/speed_of_light_in_SI, &
&                              Debye_in_SI = Hz_to_inv_cm/1d19, &   ! Cm;  1 D = 3.335640951e-30 Cm
&                              elementary_charge_in_SI = 1.6021766340d-19, &  ! C
&                              bohr_in_SI = 0.529177210903d-10, & ! m
&                              Debye_in_au = Debye_in_SI/elementary_charge_in_SI/bohr_in_SI, &  ! 1 D = 0.3934302694041317 a.u. 
&                              fstark = 1.679200537480329d-2, & ! 1/(kV*debye); multiply electric field (kV/cm) and dipole moment
                                                                ! (debye) with fstark to get Stark energy in cm^-1
&                              bohr_magneton = 4.6686447783d-5, & ! cm^-1/G
&                              nuclear_magneton = 2.54262341353d-8, & ! cm^-1/G
&                              pi = atan(1.d0) * 4.d0
end module constants
