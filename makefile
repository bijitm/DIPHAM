# A simple makefile to compile the code
main = main_dipscat_MW_Edc.f90
modules = dipole_module.f90 constants_module.f90
subroutines = monomer_mapping.f90 pairbasis_builder.f90 coupling.f90 \
	      hammat.f90 thresholds.f90 microwave_eff_pot.f90 dual_MW_eff_pot.f90 \
	      threej.f thrj.f
functions = c1func.f90 c2func.f90
all_files = $(modules) $(subroutines) $(functions) $(main)
#flags = -Wall -fbacktrace
flags = -fbacktrace
lapack_lib = -L~/lapack-3.11/ -llapack -lblas

dipham: $(all_files) 
	gfortran $(all_files) $(flags) $(lapack_lib) -o dipham

