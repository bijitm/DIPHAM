module constants
    implicit none
    real(kind=8), parameter ::  &
    pi = atan(1d0)*4d0, &

    ! NIST CODATA 2022 values
    Planck_constant_in_SI = 6.62607015d-34, &  ! Js
    speed_of_light_in_SI = 2.99792458d8, &  ! ms-1
    electron_mass = 9.1093837139d-31, & ! Kg
    atomic_mass_constant = 1.66053906892d-27, & ! Kg
    Boltzmann_constant_in_SI = 1.380649d-23, & ! JK^-1
    fine_structure_constant = 7.2973525643d-3, &  ! dimensionless
    inverse_fine_structure_const = 137.035999177d0, &  ! dimensionless
    bohr_magneton_in_SI = 9.2740100657d-24, &  ! JT-1
    nuclear_magneton_in_SI = 5.0507837393d-27, &  ! JT-1
    elementary_charge_in_SI = 1.6021766340d-19, &  ! C
    bohr_in_SI = 5.29177210544d-11, & ! m


    ! Derived units
    epsilon_zero_in_SI = elementary_charge_in_SI**2/(2d0*fine_structure_constant  &
                       * Planck_constant_in_SI*speed_of_light_in_SI), & ! 8.854187818840059e-12 C^2 Kg^-1 m^-3 s^2
    amu_to_au = atomic_mass_constant/electron_mass, & ! = 1822.8884862827601 au_mass/amu
    hbar_in_SI = Planck_constant_in_SI/2d0/pi, &  ! Js
    bohr_to_ang = bohr_in_SI*1d10, &  ! 0.529177210544 ang/bohr
    hartree_in_SI = hbar_in_SI**2/electron_mass/bohr_in_SI**2, & ! 4.359744716288645e-18 J/hartree
    autime_in_SI = hbar_in_SI/hartree_in_SI, & ! 2.4188843265857192e-17 s/au_time
    autime_to_fs = autime_in_SI*1d14, & ! 2.4188843265857192e-3 fs/au_time
    eps_to_hartree = amu_to_au*(autime_to_fs/bohr_to_ang)**2, &  ! 0.03808798851141269 hartree/eps
    hbar = autime_to_fs/eps_to_hartree, & ! in eps*fs ! 0.06350779920715995 eps fs
    hartree_to_eV = hartree_in_SI/elementary_charge_in_SI, &  ! 27.21138624596899 eV/hartree
    hartree_to_invcm = 1d-2*hartree_in_SI/Planck_constant_in_SI/speed_of_light_in_SI, &  ! 219474.6313630429 cm-1/hartree
    eps_to_eV = eps_to_hartree*hartree_to_eV, & ! 1.0364269667168056 eV/eps
    Hz_to_invcm = 1d-2/speed_of_light_in_SI, & ! 3.335640951981521e-11 cm-1/Hz
    MHz_to_invcm = Hz_to_invcm*1d6, & ! 3.335640951981521e-5 cm-1/MHz
    K_to_invcm = Boltzmann_constant_in_SI*hartree_to_invcm/hartree_in_SI, & ! 0.695034800486625 cm-1/K
    K_to_MHz = K_to_invcm/MHz_to_invcm, & ! 20836.61912332757 MHz/K
    eV_to_invcm = hartree_to_invcm/hartree_to_eV, & ! 8065.543937349212 cm-1/eV
    eV_to_MHz = eV_to_invcm/MHz_to_invcm, &  ! 241798924.20849177 MHz/eV
    eps_to_invcm = eps_to_eV*eV_to_invcm, &  ! 8359.347237902113 cm-1/eps
    eps_to_MHz = eps_to_invcm/MHz_to_invcm, &  ! 250606925.57261848 MHz/eps
    eps_to_K = eps_to_invcm/K_to_invcm, &    ! 12027.235516919935 K/eps
    hartree_to_eps = 1d0/eps_to_hartree, & ! 26.254996367171234 eps/hartree
    bfct = hbar**2*eps_to_invcm/2d0, & ! 16.85762919164018 cm^-1 amu ang^2
    Joule_to_invcm = 1d-2/Planck_constant_in_SI/speed_of_light_in_SI, &
    Debye_in_SI = Hz_to_invcm/1d19, &   ! 3.335640951981521e-30 Cm/D
    Debye_in_au = Debye_in_SI/elementary_charge_in_SI/bohr_in_SI, &  ! 0.3934302697868072 au/D
    bohr_magneton = bohr_magneton_in_SI*Joule_to_invcm*1d-4, & ! 4.668644771929823e-05 cm-1/G
    nuclear_magneton = nuclear_magneton_in_SI*Joule_to_invcm*1d-4, & ! 2.54262341353e-8 cm-1/G
    Efld_au_to_SI = elementary_charge_in_SI/bohr_in_SI**2/(4d0*pi*epsilon_zero_in_SI), &
                    ! 5.142206751090901e11 Vm-1/au
    Efld_au_to_kVinvcm = 1d-5*Efld_au_to_SI, &  ! 5.142206751090901e6 kVcm-1/au
    fstark = Debye_in_au*hartree_to_invcm/Efld_au_to_kVinvcm ! 0.01679200537983106 kV-1 D-1
             ! multiply electric field (kV/cm) and dipole moment (D) with fstark to get energy in cm-1
    
end module constants
