# Makefile for amr2cube

OPTIONS = -cpp -DDISPLAY_PROGRESS_PERCENT

F90     = gfortran -O3 -ffree-line-length-none -fopenmp -fbacktrace ###-ffree-form 
DEBUG   = -g -O2 -ffpe-trap=invalid -fbacktrace -fbounds-check -frange-check -Wall
FFLAGS += $(OPTIONS)



all: ramses_info.o hydro_reader.o amr2cube.o
	$(F90) ramses_info.o hydro_reader.o amr2cube.o -o amr2cube

%.o: %.f90
	$(F90C) $(FFLAGS) $(HOPTIONS) -c $*.f90

%.s: %.f90
	$(F90C) $(FFLAGS) $(GOPTIONS) -c $*.f90

#############################################################################
%.o:%.f90
	$(F90)  $(FFLAGS) -c $^ -o $@
#############################################################################
clean :
	rm *.o *.mod 
#############################################################################
