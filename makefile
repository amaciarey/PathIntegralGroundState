#f90=ifort -fast 
f90=gfortran -Ofast -Wall -fno-range-check 

vpi: random_mod.mod\
 	global_mod.mod\
 	system_mod.mod\
 	sample_mod.mod\
 	vpi_mod.mod\
 	pbc_mod.mod\
 	random_mod.o\
 	global_mod.o\
 	system_mod.o\
 	sample_mod.o\
 	vpi_mod.o\
 	vpi.o\
 	r8_gamma.o\
 	interpolate.o
	$(f90) -o vpi *.o

vpi.o: vpi.f90\
 	random_mod.mod\
 	global_mod.mod\
 	system_mod.mod\
 	vpi_mod.mod\
 	pbc_mod.mod
	$(f90) -c vpi.f90

sample_mod.o sample_mod.mod: sample_mod.f90\
 	global_mod.mod\
 	pbc_mod.mod\
 	vpi_mod.mod
	$(f90) -c sample_mod.f90

vpi_mod.o vpi_mod.mod: vpi_mod.f90\
 	random_mod.mod\
 	global_mod.mod\
 	system_mod.mod\
 	pbc_mod.mod\
 	interpolate.f90 
	$(f90) -c vpi_mod.f90

system_mod.o system_mod.mod: system_mod.f90\
 	global_mod.mod\
 	bessel_mod.mod
	$(f90) -c system_mod.f90

pbc_mod.o pbc_mod.mod: pbc_mod.f90\
 	global_mod.mod
	$(f90) -c pbc_mod.f90

global_mod.o global_mod.mod: global_mod.f90
	$(f90) -c global_mod.f90

bessel_mod.o bessel_mod.mod: bessel_mod.f90
	$(f90) -c bessel_mod.f90

random_mod.o random_mod.mod: random_mod.f90
	$(f90) -c random_mod.f90

r8_gamma.o: r8_gamma.f90
	$(f90) -c r8_gamma.f90

interpolate.o: interpolate.f90
	$(f90) -c interpolate.f90

clean:
	rm -f *.o *.mod 
