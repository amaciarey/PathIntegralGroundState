#F90=gfortran
#FFLAGS=-Ofast -fno-range-check -flto -funroll-loops -fwhole-program -fopenmp

F90=ifort
#FFLAGS=-O3 -ipo -no-prec-div -p -g 
FFLAGS=-fast -qopenmp

vpi: \
 	random_mod.o\
 	global_mod.o\
 	system_mod.o\
 	sample_mod.o\
 	vpi_mod.o\
 	vpi.o\
 	r8_gamma.o\
 	interpolate.o
	$(F90) $(FFLAGS) -o vpi *.o

vpi.o: \
	vpi.f90\
 	random_mod.mod\
 	global_mod.mod\
 	system_mod.mod\
 	vpi_mod.mod\
 	pbc_mod.mod
	$(F90) $(FFLAGS) -c vpi.f90

sample_mod.o sample_mod.mod: \
	sample_mod.f90\
 	global_mod.mod\
 	system_mod.mod\
	pbc_mod.mod\
 	vpi_mod.mod
	$(F90) $(FFLAGS) -c sample_mod.f90

vpi_mod.o vpi_mod.mod: \
	vpi_mod.f90\
 	random_mod.mod\
 	global_mod.mod\
 	system_mod.mod\
 	pbc_mod.mod\
 	interpolate.f90 
	$(F90) $(FFLAGS) -c vpi_mod.f90

system_mod.o system_mod.mod: \
	system_mod.f90\
 	global_mod.mod\
 	bessel_mod.mod
	$(F90) $(FFLAGS) -c system_mod.f90

pbc_mod.o pbc_mod.mod: \
	pbc_mod.f90\
 	global_mod.mod
	$(F90) $(FFLAGS) -c pbc_mod.f90

global_mod.o global_mod.mod: \
	global_mod.f90
	$(F90) $(FFLAGS) -c global_mod.f90

bessel_mod.o bessel_mod.mod: \
	bessel_mod.f90
	$(F90) $(FFLAGS) -c bessel_mod.f90

random_mod.o random_mod.mod: \
	random_mod.f90
	$(F90) $(FFLAGS) -c random_mod.f90

r8_gamma.o: \
	r8_gamma.f90
	$(F90) $(FFLAGS) -c r8_gamma.f90

interpolate.o: \
	interpolate.f90
	$(F90) $(FFLAGS) -c interpolate.f90

clean:
	rm -f *.o *.mod
